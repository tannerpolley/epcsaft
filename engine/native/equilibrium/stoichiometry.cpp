#include "equilibrium/stoichiometry.hpp"

#include "core/declared.hpp"
#include "core/physical_constants.hpp"
#include "core/trace.hpp"
#include "properties/thermochemistry.hpp"

#include <Highs.h>

#include <Eigen/Dense>

#include <algorithm>
#include <array>
#include <cmath>
#include <iterator>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>

namespace epcsaft::equilibrium {
namespace {

void require_nonnegative(const std::vector<double>& values, const char* name) {
    if (values.empty() || !std::all_of(values.begin(), values.end(), [](double value) {
            return std::isfinite(value) && value >= 0.0;
        })) {
        throw std::invalid_argument(std::string("invalid ") + name);
    }
}

void require_positive(double value, const char* name) {
    if (!std::isfinite(value) || value <= 0.0) {
        throw std::invalid_argument(std::string(name) + " must be positive");
    }
}

// A fraction feed is admitted within the composition policy and never renormalized.
void require_unit_sum(const std::vector<double>& values, const char* name) {
    if (!(std::abs(sum(values) - 1.0) <= policy::COMPOSITION_SUM_TOLERANCE.value)) {
        throw std::invalid_argument(std::string(name) + " must sum to one");
    }
}

// Callers validate the set first (validate_reaction_set).
Eigen::MatrixXd reaction_matrix(
    std::size_t components, const std::vector<ReactionSpec>& reactions
) {
    Eigen::MatrixXd matrix(
        static_cast<Eigen::Index>(reactions.size()), static_cast<Eigen::Index>(components)
    );
    for (std::size_t row = 0; row < reactions.size(); ++row) {
        const auto& values = reactions[row].stoichiometry;
        matrix.row(static_cast<Eigen::Index>(row)) = Eigen::Map<const Eigen::VectorXd>(
            values.data(), static_cast<Eigen::Index>(components)).transpose();
    }
    return matrix;
}

Eigen::Index numerical_rank(const Eigen::MatrixXd& matrix) {
    if (matrix.rows() == 0 || matrix.cols() == 0) return 0;
    const auto svd = matrix.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV);
    const double largest = svd.singularValues().size() == 0
        ? 0.0 : svd.singularValues()(0);
    const double scale = std::numeric_limits<double>::epsilon()
        * static_cast<double>(std::max(matrix.rows(), matrix.cols()))
        * std::max(1.0, largest);
    Eigen::Index rank = 0;
    for (const double singular : svd.singularValues()) {
        if (singular > scale) ++rank;
    }
    return rank;
}

}  // namespace

bool exact_zero_dot(const std::vector<double>& left, const std::vector<double>& right) {
    if (left.size() != right.size()) throw std::invalid_argument("exact dot width");
    std::vector<double> expansion;
    const auto add = [&](double value) {
        std::vector<double> next;
        for (double term : expansion) {
            const double sum = value + term, recovered = sum - value;
            const double error = (value - (sum - recovered)) + (term - recovered);
            if (error != 0.0) next.push_back(error);
            value = sum;
        }
        if (value != 0.0) next.push_back(value);
        expansion = std::move(next);
    };
    for (std::size_t i = 0; i < left.size(); ++i) {
        if (left[i] == 0.0 || right[i] == 0.0) continue;
        const double product = left[i] * right[i];
        // Decline certification if an error-free product residual could underflow.
        if (!std::isfinite(product) || std::abs(product) < std::ldexp(std::numeric_limits<double>::min(), 53)) return false;
        add(std::fma(left[i], right[i], -product));
        add(product);
    }
    return expansion.empty();
}

std::vector<double> reconstruct_conserved_functional(std::vector<double> candidate,
    const std::vector<ReactionSpec>& reactions) {
    const auto conserved = [&](const auto& row) {
        return std::all_of(reactions.begin(), reactions.end(),
            [&](const auto& reaction) { return exact_zero_dot(row, reaction.stoichiometry); });
    };
    if (conserved(candidate)) return candidate;
    double largest = 0.0, scale = std::numeric_limits<double>::infinity();
    for (double value : candidate) largest = std::max(largest, std::abs(value));
    // Rounding proposes a rational functional; exact original stoichiometry alone admits it.
    for (double value : candidate) if (std::abs(value) > std::sqrt(std::numeric_limits<double>::epsilon()) * largest)
        scale = std::min(scale, std::abs(value));
    if (!std::isfinite(scale)) return {};
    for (double& value : candidate) value = std::round(value / scale);
    return conserved(candidate) ? candidate : std::vector<double>{};
}

void restrict_inventory_balances(CompiledSystem& compiled, const std::vector<bool>& active) {
    std::vector<std::size_t> columns;
    for (std::size_t i = 0; i < active.size(); ++i) if (active[i]) columns.push_back(i);
    Eigen::MatrixXd restricted(compiled.B.size(), columns.size());
    for (std::size_t k = 0; k < compiled.B.size(); ++k)
        for (std::size_t j = 0; j < columns.size(); ++j) restricted(k, j) = compiled.B[k][columns[j]];
    const auto rank = static_cast<std::size_t>(numerical_rank(restricted));
    auto basis = compiled.B;
    // The shared SVD owner decides rank; full-column pivot operations only realise its L B.
    for (std::size_t k = 0; k < rank; ++k) {
        std::size_t pivot = k, column = columns.front();
        double largest = 0.0;
        for (std::size_t r = k; r < basis.size(); ++r) for (auto i : columns)
            if (std::abs(basis[r][i]) > largest) { largest = std::abs(basis[r][i]); pivot = r; column = i; }
        if (!(largest > 0.0)) throw std::invalid_argument("rank-deficient inventory face");
        std::swap(basis[pivot], basis[k]);
        const double divisor = basis[k][column];
        for (double& value : basis[k]) value /= divisor;
        for (std::size_t r = 0; r < basis.size(); ++r) if (r != k) {
            const double factor = basis[r][column];
            for (std::size_t i = 0; i < active.size(); ++i) basis[r][i] -= factor * basis[k][i];
            basis[r][column] = 0.0;
        }
        columns.erase(std::find(columns.begin(), columns.end(), column));
    }
    for (std::size_t k = rank; k < basis.size(); ++k) {
        auto restriction = reconstruct_conserved_functional(basis[k], compiled.problem.reactions);
        bool state_zero = !restriction.empty();
        if (state_zero) for (std::size_t i = 0; i < active.size(); ++i)
            state_zero = state_zero && (!active[i] || restriction[i] == 0.0);
        if (restriction.empty() || !exact_zero_dot(restriction, compiled.n_feed)
            || !state_zero) throw std::invalid_argument("infeasible_inventory: discarded state balance");
        compiled.inventory_restrictions.push_back(std::move(restriction));
    }
    basis.resize(rank);
    compiled.B = basis;
    compiled.rank = rank;
    compiled.balance = matvec(compiled.B, compiled.n_feed);
    for (auto& row : basis) {
        std::vector<double> state;
        for (std::size_t i = 0; i < active.size(); ++i) if (active[i]) state.push_back(row[i]);
        row = std::move(state);
    }
    compiled.conserves_total_moles = conserves_total_moles(basis);
}

std::vector<InventoryBound> inventory_bounds(const CompiledSystem& c,
    const std::vector<std::size_t>& indices, bool fraction, double floor) {
    core::Scope scope("inventory_bounds");
    const std::size_t n = c.n_feed.size(), width = n + fraction;
    std::vector<double> lower(width, 0.0), upper(width, 0.0);
    for (const auto& phase : c.problem.phases)
        for (const auto i : phase.support) upper[i] = kHighsInf;
    const auto count = std::count(upper.begin(), upper.end(), kHighsInf);
    if (!std::isfinite(floor) || floor < 0.0 || count * floor > 1.0 || (!fraction && floor != 0.0))
        throw std::invalid_argument("invalid conserved-inventory fraction floor");
    for (const auto i : indices) if (i >= n) throw std::invalid_argument("inventory species out of range");
    for (std::size_t i = 0; i < n; ++i) if (upper[i] != 0.0 && fraction) lower[i] = floor;
    if (fraction) upper[n] = kHighsInf;
    Highs lp;
    lp.setOptionValue("output_flag", false);
    lp.setOptionValue("threads", 1);
    bool assembled = lp.addVars(static_cast<HighsInt>(width), lower.data(), upper.data()) != HighsStatus::kError;
    // Solve in q=tau_scale*tau: a tiny nonzero inventory must not be dropped by HiGHS.
    double tau_scale = 1.0;
    if (fraction) for (const double value : c.balance)
        if (value != 0.0) tau_scale = std::min(tau_scale, std::abs(value / c.feed_total));
    std::vector<HighsInt> columns(width);
    std::iota(columns.begin(), columns.end(), 0);
    for (std::size_t k = 0; k < c.rank; ++k) {
        auto row = c.B[k];
        const double rhs = c.balance[k] / c.feed_total;
        if (fraction) row.push_back(-rhs / tau_scale);
        assembled = lp.addRow(fraction ? 0.0 : rhs, fraction ? 0.0 : rhs,
            static_cast<HighsInt>(width), columns.data(), row.data()) != HighsStatus::kError && assembled;
    }
    if (fraction) {
        std::vector<double> row(width, 1.0);
        row[n] = 0.0;
        assembled = lp.addRow(1.0, 1.0, static_cast<HighsInt>(width), columns.data(), row.data())
            != HighsStatus::kError && assembled;
    }
    lp.changeObjectiveSense(ObjSense::kMaximize);
    const auto run = [&]() {
        if (!assembled) return std::string("numerical_failure");
        const auto outcome = lp.run();
        const auto status = lp.getModelStatus();
        if (outcome == HighsStatus::kError) return std::string("numerical_failure");
        if (status == HighsModelStatus::kOptimal) return std::string("optimal");
        if (status == HighsModelStatus::kInfeasible) return std::string("infeasible_inventory");
        if (status == HighsModelStatus::kUnbounded) return std::string("unbounded_amount");
        if (status == HighsModelStatus::kTimeLimit || status == HighsModelStatus::kIterationLimit)
            return std::string("exhausted_work");
        return std::string("numerical_failure");
    };
    const auto diagnose = [&](InventoryBound& r, const std::vector<double>& x) {
        r.witness = x;
        if (x.size() != width || !std::all_of(x.begin(), x.end(), [](double v) { return std::isfinite(v); }))
            return false;
        r.balance_residual = r.normalization_residual = r.lower_bound_residual = 0.0;
        for (std::size_t k = 0; k < c.rank; ++k) {
            const double value = std::inner_product(c.B[k].begin(), c.B[k].end(), x.begin(), 0.0);
            r.balance_residual = std::max(r.balance_residual,
                std::abs(value - c.balance[k] / c.feed_total * (fraction ? x[n] : 1.0)));
        }
        if (fraction) r.normalization_residual = std::abs(std::accumulate(x.begin(), x.begin() + n, 0.0) - 1.0);
        for (std::size_t j = 0; j < width; ++j)
            r.lower_bound_residual = std::max(r.lower_bound_residual,
                std::max(lower[j] - x[j], upper[j] == 0.0 ? std::abs(x[j]) : 0.0));
        return std::max({r.balance_residual, r.normalization_residual, r.lower_bound_residual})
            <= policy::EQUILIBRIUM_ROW_RESIDUAL.value;
    };
    InventoryBound finite;
    std::string feasibility = "optimal", feasibility_cause;
    if (fraction) {
        lp.changeColCost(static_cast<HighsInt>(n), 1.0);
        feasibility = run();
        finite.lp_error = lp.getInfo().max_primal_infeasibility;
        feasibility_cause = lp.modelStatusToString(lp.getModelStatus());
        if (feasibility == "optimal" || feasibility == "unbounded_amount") {
            auto x = lp.getSolution().col_value;
            if (feasibility == "unbounded_amount") {
                bool has_ray = false;
                std::vector<double> ray(width);
                lp.getPrimalRay(has_ray, ray.data());
                if (has_ray && ray[n] > 0.0) for (std::size_t j = 0; j < width; ++j) x[j] += ray[j];
            }
            x[n] /= tau_scale;
            const bool valid = diagnose(finite, x);
            feasibility = !valid ? "numerical_failure" : x[n] > 0.0 ? "optimal"
                : feasibility == "optimal" ? "infeasible_inventory" : "numerical_failure";
        }
        lp.changeColCost(static_cast<HighsInt>(n), 0.0);
    }
    std::vector<InventoryBound> results;
    for (const auto i : indices) {
        InventoryBound r = finite;
        r.species = i;
        r.supported = upper[i] != 0.0;
        r.floor = floor;
        r.finite_witness = feasibility == "optimal" ? finite.witness : std::vector<double>{};
        r.feasibility_termination = feasibility_cause;
        r.finite_residual = std::max({finite.balance_residual, finite.normalization_residual, finite.lower_bound_residual});
        r.status = feasibility;
        r.termination = feasibility_cause;
        if (feasibility == "optimal") {
            lp.changeColCost(static_cast<HighsInt>(i), 1.0);
            r.status = run();
            r.termination = lp.modelStatusToString(lp.getModelStatus());
            r.lp_error = lp.getInfo().max_primal_infeasibility;
            if (r.status == "optimal") {
                auto x = lp.getSolution().col_value;
                if (fraction) x[n] /= tau_scale;
                if (!fraction && !diagnose(r, x)) {
                    r.witness.clear();
                    r.balance_residual = r.normalization_residual = r.lower_bound_residual = NAN;
                    r.status = lp.clearSolver() == HighsStatus::kError ? "numerical_failure" : run();
                    r.termination = lp.modelStatusToString(lp.getModelStatus());
                    r.lp_error = lp.getInfo().max_primal_infeasibility;
                    if (r.status == "optimal") x = lp.getSolution().col_value;
                }
                if (r.status != "optimal" || !diagnose(r, x)) r.status = "numerical_failure";
                else {
                    r.maximum = x[i] * (fraction ? 1.0 : c.feed_total);
                    if (fraction && x[i] > 0.0) r.log_fraction = std::log(x[i]);
                    if (!fraction) r.feed_normalized_amount = x[i];
                    if (!fraction && x[i] == 0.0) {
                        r.conserved_functional.assign(n, 0.0);
                        const auto& dual = lp.getSolution().row_dual;
                        for (std::size_t k = 0; k < c.rank; ++k)
                            for (std::size_t j = 0; j < n; ++j)
                                r.conserved_functional[j] += dual[k] * c.B[k][j];
                    }
                    r.limiting = fraction && x[n] == 0.0;
                    if (!r.limiting) for (std::size_t j = 0; j < n; ++j)
                        r.amounts_mol.push_back(x[j] * c.feed_total / (fraction ? x[n] : 1.0));
                    if (upper[i] == 0.0) r.status = "structural_zero";
                }
            }
            lp.changeColCost(static_cast<HighsInt>(i), 0.0);
        }
        results.push_back(std::move(r));
    }
    return results;
}

std::vector<double> resolve_feed(const HelmholtzModel& model, const FeedSpec& feed) {
    const std::size_t components = model.parameters().component_count();
    if (feed.values.size() != components) {
        throw std::invalid_argument("feed size does not match components");
    }
    require_nonnegative(feed.values, "feed");
    std::vector<double> result(components, 0.0);
    if (feed.kind == FeedKind::Amounts) {
        result = feed.values;
    } else if (feed.kind == FeedKind::MoleFractions) {
        require_unit_sum(feed.values, "mole fractions");
        result = feed.values;
    } else if (feed.kind == FeedKind::MassFractions) {
        const auto& masses = model.parameters().molar_mass_kg_per_mol;
        require_unit_sum(feed.values, "mass fractions");
        for (std::size_t i = 0; i < components; ++i) {
            require_positive(masses[i], "molar mass");
            result[i] = feed.values[i] / masses[i];
        }
    } else {
        if (feed.solvent.empty()) throw std::invalid_argument("molality requires solvent");
        require_positive(feed.solvent_mass_kg, "solvent mass");
        std::vector<bool> is_solvent(components, false);
        double solvent_mass = 0.0;
        std::size_t supplied_solvents = 0;
        for (const std::size_t index : feed.solvent) {
            if (index >= components || is_solvent[index]) {
                throw std::invalid_argument("invalid solvent set");
            }
            is_solvent[index] = true;
            const double amount = feed.values[index];
            if (amount > 0.0) {
                require_positive(model.parameters().molar_mass_kg_per_mol[index], "molar mass");
                solvent_mass += amount * model.parameters().molar_mass_kg_per_mol[index];
                ++supplied_solvents;
            }
        }
        if (supplied_solvents == 0 && feed.solvent.size() > 1) {
            throw std::invalid_argument("molality feed omits the relative solvent amounts");
        }
        if (supplied_solvents == 0) solvent_mass = feed.solvent_mass_kg;
        if (!(solvent_mass > 0.0)) throw std::invalid_argument("solvent has zero mass");
        const double mass_scale = feed.solvent_mass_kg / solvent_mass;
        for (std::size_t i = 0; i < components; ++i) {
            if (is_solvent[i]) {
                result[i] = supplied_solvents == 0
                    ? feed.solvent_mass_kg / static_cast<double>(feed.solvent.size())
                        / model.parameters().molar_mass_kg_per_mol[i]
                        : feed.values[i] * mass_scale;
            } else {
                result[i] = feed.values[i] * feed.solvent_mass_kg;
            }
        }
    }
    const double total = std::accumulate(result.begin(), result.end(), 0.0);
    if (!(std::isfinite(total) && total > 0.0)) {
        throw std::invalid_argument("feed has zero total");
    }
    return result;
}

std::vector<double> resolve_feed(const IdealMoleFraction& model, const FeedSpec& feed) {
    if (model.activity_basis == "molarity" && feed.kind != FeedKind::Amounts)
        throw std::invalid_argument("molarity feed requires absolute mol amounts");
    if (feed.kind != FeedKind::Amounts && feed.kind != FeedKind::MoleFractions)
        throw std::invalid_argument("ideal mole-fraction feed requires mol or mole fractions");
    require_nonnegative(feed.values, "feed");
    if (feed.values.size() != model.component_ids.size()) throw std::invalid_argument("feed size does not match species");
    if (feed.kind == FeedKind::MoleFractions) require_unit_sum(feed.values, "mole fractions");
    require_positive(sum(feed.values), "feed total");
    return feed.values;
}

std::vector<std::vector<double>> balance_basis(
    std::size_t components, const std::vector<ReactionSpec>& reactions
) {
    if (components == 0) throw std::invalid_argument("no components");
    validate_reaction_set(components, reactions);
    const auto matrix = reaction_matrix(components, reactions);
    const Eigen::Index rank = numerical_rank(matrix);
    if (rank != matrix.rows()) {
        for (std::size_t row = 0; row < reactions.size(); ++row) {
            Eigen::MatrixXd prefix = matrix.topRows(static_cast<Eigen::Index>(row + 1));
            if (numerical_rank(prefix) <= static_cast<Eigen::Index>(row)) {
                const std::string label = reaction_name(reactions[row], row);
                throw std::invalid_argument("dependent reaction " + label);
            }
        }
        throw std::invalid_argument("dependent reaction set");
    }
    if (rank == 0) {
        std::vector<std::vector<double>> result(components, std::vector<double>(components, 0.0));
        for (std::size_t i = 0; i < components; ++i) result[i][i] = 1.0;
        return result;
    }
    // Reduce the SVD null space to echelon form: each row is one conserved moiety with a unit
    // pivot, so a trace moiety (Li+ + LiA) is not mixed into rows of major species, whose
    // rounding would swamp its balance and actions (#172). Pivots are the species taking part
    // in the most reactions (shared H+ or carrier; ties by index), leaving each reaction's
    // peripheral species off the pivots; this is structural, not feed dependent. A candidate
    // is a pivot when the remaining species still carry the full reaction rank, tested on the
    // exact stoichiometry: the SVD columns of dependent species differ only by rounding.
    // Coefficients below sqrt(eps) are rounding of the rational stoichiometry and set to zero.
    const Eigen::Index rows = static_cast<Eigen::Index>(components) - rank;
    const Eigen::MatrixXd space = matrix.jacobiSvd(Eigen::ComputeFullV).matrixV().rightCols(rows).transpose();
    std::vector<Eigen::Index> order(components), pivots, rest(components);
    std::iota(order.begin(), order.end(), 0);
    std::iota(rest.begin(), rest.end(), 0);
    std::stable_sort(order.begin(), order.end(), [&](Eigen::Index a, Eigen::Index b) {
        return (matrix.col(a).array() != 0.0).count() > (matrix.col(b).array() != 0.0).count();
    });
    for (const Eigen::Index species : order) {
        rest.erase(std::find(rest.begin(), rest.end(), species));
        if (numerical_rank(matrix(Eigen::all, rest)) == rank) pivots.push_back(species);
        else rest.push_back(species);
    }
    const Eigen::MatrixXd echelon = space(Eigen::all, pivots).partialPivLu().solve(space);
    std::vector<std::vector<double>> result(
        static_cast<std::size_t>(rows), std::vector<double>(components, 0.0)
    );
    for (Eigen::Index row = 0; row < rows; ++row) {
        for (std::size_t column = 0; column < components; ++column) {
            const double value = echelon(row, static_cast<Eigen::Index>(column));
            result[static_cast<std::size_t>(row)][column] =
                std::abs(value) > std::sqrt(std::numeric_limits<double>::epsilon()) ? value : 0.0;
        }
    }
    return result;
}

std::vector<std::vector<double>> reaction_target_series(
    double temperature, const std::vector<ReactionSpec>& reactions, std::size_t order,
    const ThermochemistryRecord* thermochemistry, const std::string& standard_state_id
) {
    require_positive(temperature, "temperature");
    if (order > 3) throw std::invalid_argument("reaction target order");
    const double R = GAS_CONSTANT_J_PER_MOL_K;
    std::vector<std::vector<double>> result(order + 1, std::vector<double>(reactions.size(), 0.0));
    for (std::size_t row = 0; row < reactions.size(); ++row) {
        const auto& reaction = reactions[row];
        if (reaction.uses_species_thermochemistry()) {
            const auto label = "reaction " + reaction_name(reaction, row) + " at T=" + std::to_string(temperature) + " K: ";
            if (!thermochemistry) throw std::invalid_argument(label + "MissingRecord");
            const auto source = reaction_thermochemistry_density(*thermochemistry, reaction.stoichiometry, temperature);
            if (source.status != ThermochemistryStatus::Available) throw std::invalid_argument(label + source.cause);
            for (std::size_t k = 0; k <= order; ++k) result[k][row] = source.gibbs_derivatives[k];
            continue;
        }
        // l[k] = d^k ln K / dT^k; t = -R T ln K, so t^(k) = -R (k l[k-1] + T l[k]).
        std::array<double, 4> l{};
        if (reaction.affinity_j_per_mol.has_value()) {
            result[0][row] = *reaction.affinity_j_per_mol;
            continue;
        }
        if (reaction.correlation.has_value()) {
            if (reaction.correlation->standard_state_id != standard_state_id) {
                throw std::invalid_argument("reference-basis-mismatch in reaction "
                                            + reaction_name(reaction, row));
            }
            const auto properties = evaluate_reaction_properties(*reaction.correlation, temperature);
            l = {properties.ln_k, properties.derivative, properties.second_derivative,
                 properties.third_derivative};
        } else {
            l[0] = std::log(*reaction.equilibrium_constant);
        }
        for (std::size_t k = 0; k <= order; ++k) {
            result[k][row] = -R * ((k == 0 ? 0.0 : static_cast<double>(k) * l[k - 1])
                                   + temperature * l[k]);
        }
    }
    return result;
}

std::optional<std::size_t> thermal_reference_missing(
    double temperature, const std::vector<ReactionSpec>& reactions, const std::vector<double>& supplied, const ThermochemistryRecord* thermochemistry
) {
    if (supplied.empty() || reactions.empty()) return std::nullopt;
    const auto targets = reaction_target_series(temperature, reactions, 3, thermochemistry);
    for (std::size_t row = 0; row < reactions.size(); ++row) {
        const double scale = std::sqrt(std::numeric_limits<double>::epsilon()) * std::max(1.0, std::abs(targets[0][row]));
        double power = 1.0;
        for (std::size_t k = 1; k <= 3; ++k) {
            power *= temperature;
            if (std::abs(targets[k][row]) * power > scale) return row;
        }
    }
    return std::nullopt;
}

std::vector<std::vector<double>> standard_state_series(
    double temperature, std::size_t components,
    const std::vector<ReactionSpec>& reactions,
    const std::vector<double>& supplied, std::size_t order, const ThermochemistryRecord* thermochemistry
) {
    require_positive(temperature, "temperature");
    validate_reaction_set(components, reactions);
    std::vector<std::vector<double>> result(order + 1, std::vector<double>(components, 0.0));
    if (reactions.empty() && supplied.empty()) return result;
    const auto targets = reaction_target_series(temperature, reactions, order, thermochemistry);
    if (reactions.empty()) {
        if (supplied.size() != components) throw std::invalid_argument("standard state size");
        if (!all_numeric(supplied)) throw std::invalid_argument("invalid standard state");
        result[0] = supplied;
        return result;
    }
    const auto matrix = reaction_matrix(components, reactions);
    Eigen::VectorXd target(static_cast<Eigen::Index>(reactions.size()));
    for (std::size_t row = 0; row < reactions.size(); ++row) {
        target(static_cast<Eigen::Index>(row)) = targets[0][row];
    }
    const double criterion = std::sqrt(std::numeric_limits<double>::epsilon())
        * std::max(1.0, target.norm());
    if (!supplied.empty()) {
        if (supplied.size() != components) throw std::invalid_argument("standard state size");
        if (!all_numeric(supplied)) throw std::invalid_argument("invalid standard state");
        const Eigen::Map<const Eigen::VectorXd> gibbs(supplied.data(),
                                                      static_cast<Eigen::Index>(components));
        if ((matrix * gibbs - target).norm() > criterion) {
            throw std::invalid_argument("inconsistent standard state");
        }
        result[0] = supplied;
        return result;
    }
    const auto svd = matrix.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV);
    for (std::size_t k = 0; k <= order; ++k) {
        Eigen::VectorXd right(static_cast<Eigen::Index>(reactions.size()));
        for (std::size_t row = 0; row < reactions.size(); ++row) {
            right(static_cast<Eigen::Index>(row)) = targets[k][row];
        }
        const Eigen::VectorXd solution = svd.solve(right);
        if (k == 0 && (matrix * solution - target).norm() > criterion) {
            throw std::invalid_argument("inconsistent standard state");
        }
        result[k].assign(solution.data(), solution.data() + solution.size());
    }
    return result;
}

std::vector<double> standard_state(
    double temperature, std::size_t components,
    const std::vector<ReactionSpec>& reactions,
    const std::vector<double>& supplied, const ThermochemistryRecord* thermochemistry
) {
    return standard_state_series(temperature, components, reactions, supplied, 0, thermochemistry).front();
}

bool conserves_total_moles(const std::vector<std::vector<double>>& basis) {
    if (basis.empty()) throw std::invalid_argument("empty balance basis");
    const std::size_t components = basis.front().size();
    Eigen::MatrixXd extended(static_cast<Eigen::Index>(basis.size() + 1),
                             static_cast<Eigen::Index>(components));
    for (std::size_t row = 0; row < basis.size(); ++row) {
        extended.row(static_cast<Eigen::Index>(row)) = Eigen::Map<const Eigen::VectorXd>(
            basis[row].data(), static_cast<Eigen::Index>(components)).transpose();
    }
    extended.row(static_cast<Eigen::Index>(basis.size())).setConstant(
        1.0 / std::sqrt(static_cast<double>(components)));
    return numerical_rank(extended) == numerical_rank(extended.topRows(
        static_cast<Eigen::Index>(basis.size())));
}

std::vector<bool> reaction_participants(std::size_t components, const std::vector<ReactionSpec>& reactions) {
    std::vector<bool> result(components, false);
    for (const auto& reaction : reactions) {
        if (reaction.stoichiometry.size() != components) throw std::invalid_argument("stoichiometry size");
        for (std::size_t i = 0; i < components; ++i) result[i] = result[i] || reaction.stoichiometry[i] != 0.0;
    }
    return result;
}

std::vector<std::size_t> resolve_support(const PhaseSpec& phase, const std::vector<std::size_t>& candidates,
    const std::vector<int>& z, const Parameters* parameters) {
    const bool vapor = phase.kind == PhaseKind::Vapor;
    auto result = phase.support;
    if (result.empty() && !phase.support_declared)
        std::copy_if(candidates.begin(), candidates.end(), std::back_inserter(result),
            [&](std::size_t i) { return !vapor || z[i] == 0; });
    if (result.empty()) throw core::DeclarationError("empty support in phase " + phase.name);
    std::sort(result.begin(), result.end());
    if (std::adjacent_find(result.begin(), result.end()) != result.end())
        throw core::DeclarationError("duplicate support index");
    if (result.back() >= z.size()) throw core::DeclarationError("support index out of range");
    for (const std::size_t i : result)
        if (vapor && z[i] != 0 && !(parameters && parameters->permittivity_form.is(PermittivityForm::DensityTemperature)))
            throw core::DeclarationError("vapor phase " + phase.name + " support contains charged component index "
                + std::to_string(i) + " (z=" + std::to_string(z[i])
                + (parameters && parameters->permittivity_form.is(PermittivityForm::DipolarPerturbation)
                    ? "); explicit charged vapor under DipolarPerturbation is not numerically qualified"
                    : "); charged vapor outside the condensed domain of the selected permittivity"));
    return result;
}

std::vector<double> matvec(
    const std::vector<std::vector<double>>& matrix, const std::vector<double>& vector
) {
    std::vector<double> result(matrix.size(), 0.0);
    for (std::size_t row = 0; row < matrix.size(); ++row) {
        if (matrix[row].size() != vector.size()) throw std::invalid_argument("matrix width");
        result[row] = std::inner_product(matrix[row].begin(), matrix[row].end(), vector.begin(), 0.0);
    }
    return result;
}

double sum(const std::vector<double>& values) {
    return std::accumulate(values.begin(), values.end(), 0.0);
}

std::string reaction_name(const ReactionSpec& reaction, std::size_t index) {
    return reaction.name.empty() ? std::to_string(index) : reaction.name;
}

void validate_reaction_set(
    std::size_t components, const std::vector<ReactionSpec>& reactions
) {
    if (components == 0) throw std::invalid_argument("no components");
    for (std::size_t row = 0; row < reactions.size(); ++row) {
        const auto& reaction = reactions[row];
        const std::string label = reaction_name(reaction, row);
        if (reaction.stoichiometry.size() != components) {
            throw std::invalid_argument("stoichiometry size does not match components");
        }
        if (!all_numeric(reaction.stoichiometry)) {
            throw std::invalid_argument("invalid stoichiometry in reaction " + label);
        }
        const int targets = reaction.equilibrium_constant.has_value()
            + reaction.affinity_j_per_mol.has_value() + reaction.correlation.has_value();
        if (targets > 1 || (targets == 0 && reaction.reference)) {
            throw std::invalid_argument("inconsistent standard state: reaction " + label
                                        + " needs at most one target; a source reference requires a target");
        }
        if (reaction.affinity_j_per_mol.has_value()
            && !std::isfinite(*reaction.affinity_j_per_mol)) {
            throw std::invalid_argument("invalid affinity in reaction " + label);
        }
        if (reaction.equilibrium_constant.has_value()
            && (!std::isfinite(*reaction.equilibrium_constant)
                || *reaction.equilibrium_constant <= 0.0)) {
            throw std::invalid_argument("equilibrium constant is not positive in reaction "
                                        + label);
        }
    }
}

std::size_t reaction_rank(
    std::size_t components, const std::vector<ReactionSpec>& reactions
) {
    validate_reaction_set(components, reactions);
    return static_cast<std::size_t>(numerical_rank(reaction_matrix(components, reactions)));
}

bool all_numeric(const std::vector<double>& values) {
    return std::all_of(values.begin(), values.end(),
                       [](double value) { return std::isfinite(value); });
}

}  // namespace epcsaft::equilibrium
