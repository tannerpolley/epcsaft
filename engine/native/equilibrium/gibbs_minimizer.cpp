#include "equilibrium/gibbs_minimizer.hpp"
#include "equilibrium/stoichiometry.hpp"
#include "equilibrium/ipopt.hpp"
#include "contributions/helmholtz.hpp"

#include "core/trace.hpp"

#include "core/numerical_policy.hpp"
#include "core/physical_constants.hpp"

#include <Highs.h>
#include <IpIpoptApplication.hpp>

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <exception>
#include <limits>
#include <memory>
#include <numeric>
#include <string>

namespace epcsaft::equilibrium {
namespace {

double maximum_abs(const std::vector<double>& values) {
    double result = 0.0;
    for (const double value : values) {
        const double magnitude = std::abs(value);
        if (!std::isfinite(magnitude)) return std::numeric_limits<double>::infinity();
        result = std::max(result, magnitude);
    }
    return result;
}

double specified_amount(const PhaseSpec& phase, const std::vector<double>& feed) {
    if (phase.amount.kind == AmountKind::Free) return 1.0;
    double total = 0.0;
    for (const double value : feed) total += value;
    return phase.amount.fraction * total;
}

// EqID reduced_hessian_stability on the balance nullspace. `minimum` is the smallest eigenvalue in the solved
// chart. The decision reads the congruent equilibration D = diag(S)^-1/2, S = |objective| + |weighted constraint|
// curvature, whose inertia is the same, against the rounding band eps (n ||DSD|| + 2 kappa(AD) ||DHD||)
// (ALG-EQUILIBRIUM, #159): trace log-chart directions have curvature far below eps ||H||.
struct Curvature final { double minimum = 0.0, equilibrated = 1.0, band = 0.0; };

Curvature reduced_curvature_at(
    const CompiledSystem& compiled, PhasePoint prepared,
    const std::vector<double>& multipliers
) {
    if (compiled.pinned_count != 0 || prepared.point.empty()) return {};
    prepared.point.resize(compiled.primal_size);
    using Dense = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    const auto n = static_cast<Eigen::Index>(prepared.point.size());
    const auto jacobian = compiled.jacobian_at(prepared);
    Eigen::MatrixXd constraints(0, n);
    for (std::size_t row = 0; row < compiled.rows.size(); ++row) {
        if (compiled.rows[row].kind == RowKind::Chemical || compiled.rows[row].kind == RowKind::Pressure) continue;
        constraints.conservativeResize(constraints.rows() + 1, n);
        constraints.bottomRows(1) = Eigen::Map<const Dense>(jacobian.data() + row * prepared.point.size(), 1, n);
    }
    const Eigen::MatrixXd objective = Eigen::Map<const Dense>(compiled.hessian_at(prepared, {}).data(), n, n);
    const Eigen::MatrixXd weighted = Eigen::Map<const Dense>(
        compiled.constraint_hessian_at(prepared, multipliers).data(), n, n);
    const Eigen::MatrixXd hessian = objective + weighted;
    const auto reduce = [n](const Eigen::MatrixXd& h, const Eigen::MatrixXd& a) -> std::pair<double, double> {
        const auto svd = a.jacobiSvd(Eigen::ComputeFullV);
        const Eigen::Index rank = svd.rank();
        if (rank >= n) return {1.0, 1.0};
        const Eigen::MatrixXd z = svd.matrixV().rightCols(n - rank);
        const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigen(z.transpose() * h * z, Eigen::EigenvaluesOnly);
        return {eigen.info() == Eigen::Success ? eigen.eigenvalues()(0) : std::numeric_limits<double>::quiet_NaN(),
                rank > 0 ? svd.singularValues()(0) / svd.singularValues()(rank - 1) : 1.0};
    };
    const Eigen::MatrixXd magnitude = objective.cwiseAbs() + weighted.cwiseAbs();
    const Eigen::VectorXd d = magnitude.diagonal().unaryExpr([](double s) { return s > 0.0 ? 1.0 / std::sqrt(s) : 1.0; });
    const Eigen::MatrixXd scaled = d.asDiagonal() * hessian * d.asDiagonal();
    const auto [equilibrated, kappa] = reduce(scaled, constraints * d.asDiagonal());
    const double rounding = static_cast<double>(n) * (d.asDiagonal() * magnitude * d.asDiagonal()).norm();
    return {reduce(hessian, constraints).first, equilibrated,
            std::numeric_limits<double>::epsilon() * (rounding + 2.0 * kappa * scaled.norm())};
}

double largest_difference(
    const std::vector<double>& left, const std::vector<double>& right
) {
    if (left.size() != right.size()) return std::numeric_limits<double>::infinity();
    double result = 0.0;
    for (std::size_t index = 0; index < left.size(); ++index) {
        result = std::max(result, std::abs(left[index] - right[index]));
    }
    return result;
}

double michelsen_deficit(
    const CompiledSystem& compiled,
    const std::vector<double>& composition,
    std::size_t phase,
    const std::vector<double>& residuals
) {
    double sum = 0.0;
    for (std::size_t row = 0; row < compiled.rows.size(); ++row) {
        const auto& item = compiled.rows[row];
        if (item.kind != RowKind::Chemical || item.phase != phase) continue;
        if (item.index >= composition.size() || row >= residuals.size()) continue;
        sum += composition[item.index] * std::exp(-residuals[row]);
    }
    return 1.0 - sum;
}

void clear_state(EquilibriumResult& result) {
    result.temperature.clear();
    result.pressure.clear();
    result.amounts.clear();
    result.mole_fractions.clear();
    result.molar_densities.clear();
    result.packing_fractions.clear();
    result.multipliers.clear();
    result.constraint_multipliers.clear();
    result.coordinates.clear();
}

void require_solution_shape(
    const EquilibriumResult& result, std::size_t phases, std::size_t components
) {
    if (!all_numeric(result.residuals) || result.rows.empty()) {
        throw std::runtime_error("equilibrium result has invalid residuals");
    }
    if (result.amounts.size() != phases || result.molar_densities.size() != phases
        || result.mole_fractions.size() != phases * components) {
        throw std::runtime_error("equilibrium result has invalid topology");
    }
    for (const double amount : result.amounts) {
        if (!std::isfinite(amount) || amount < 0.0) {
            throw std::runtime_error("equilibrium result has invalid amount");
        }
    }
    for (const double density : result.molar_densities) {
        if (!std::isfinite(density) || density <= 0.0) {
            throw std::runtime_error("equilibrium result has invalid density");
        }
    }
    for (std::size_t phase = 0; phase < phases; ++phase) {
        double total = 0.0;
        for (std::size_t component = 0; component < components; ++component) {
            const double value = result.mole_fractions[phase * components + component];
            if (!std::isfinite(value) || value < 0.0) {
                throw std::runtime_error("equilibrium result has invalid composition");
            }
            total += value;
        }
        if (std::abs(total - 1.0) > policy::COMPOSITION_SUM_TOLERANCE.value) {
            throw std::runtime_error("equilibrium result is not normalized");
        }
    }
}

std::vector<double> add_multipliers(
    const CompiledSystem& compiled, const std::vector<double>& point
) {
    if (point.size() != compiled.primal_size) return point;
    const std::size_t width = compiled.rank + compiled.charge_row_count;
    std::vector<double> result = point;
    result.resize(point.size() + width, 0.0);
    if (width == 0) return result;
    // Recover the balance and charge multipliers by the KKT-weighted least squares of the
    // chemical rows: the free-amount minimizer's stationarity conditions carry the weights
    // N_p x_i / N_feed, so the weighted fit is the one whose classified residual is minimal.
    const auto prepared = compiled.prepare(result, true);
    const auto values = compiled.evaluate_at(prepared).values;
    std::vector<std::size_t> chemical;
    std::vector<std::size_t> phi_phase;
    for (std::size_t row = 0; row < compiled.rows.size(); ++row) {
        if (compiled.rows[row].kind == RowKind::Chemical) chemical.push_back(row);
        if (compiled.rows[row].kind == RowKind::Electro) phi_phase.push_back(compiled.rows[row].phase);
    }
    Eigen::MatrixXd matrix(
        static_cast<Eigen::Index>(chemical.size()), static_cast<Eigen::Index>(width)
    );
    Eigen::VectorXd rhs(static_cast<Eigen::Index>(chemical.size()));
    const auto& z = compiled.model->parameters().charge_number;
    for (std::size_t row = 0; row < chemical.size(); ++row) {
        const auto& item = compiled.rows[chemical[row]];
        const double weight = prepared.amount[item.phase] * prepared.composition[item.phase][item.index]
            / compiled.feed_total;
        rhs(static_cast<Eigen::Index>(row)) = weight * values[chemical[row]];
        for (std::size_t column = 0; column < compiled.rank; ++column) {
            matrix(static_cast<Eigen::Index>(row), static_cast<Eigen::Index>(column)) =
                weight * compiled.B[column][item.index];
        }
        for (std::size_t column = 0; column < phi_phase.size(); ++column) {
            matrix(static_cast<Eigen::Index>(row), static_cast<Eigen::Index>(compiled.rank + column)) =
                weight * z[item.index] * (item.phase == phi_phase[column] ? 1.0 : 0.0);
        }
    }
    const Eigen::VectorXd solved = matrix.colPivHouseholderQr().solve(rhs);
    for (std::size_t column = 0; column < width; ++column) {
        result[point.size() + column] = solved(static_cast<Eigen::Index>(column));
    }
    return result;
}

bool valid_previous_phase(
    const EquilibriumResult& start, std::size_t phase, std::size_t components
) {
    if (!start.success || phase >= start.phase_count
        || start.mole_fractions.size() != start.phase_count * components
        || start.molar_densities.size() != start.phase_count) return false;
    if (!(std::isfinite(start.molar_densities[phase])
          && start.molar_densities[phase] > 0.0)) return false;
    for (std::size_t species = 0; species < components; ++species) {
        const double value = start.mole_fractions[phase * components + species];
        if (!std::isfinite(value) || value < 0.0) return false;
    }
    return true;
}

void copy_previous_phase(
    const CompiledSystem& compiled, const EquilibriumResult& start,
    std::size_t phase, std::vector<double>& point,
    const std::vector<std::size_t>& offsets
) {
    const auto& spec = compiled.problem.phases[phase];
    for (std::size_t i = 0; i < spec.support.size(); ++i) {
        point[offsets[phase] + i] = std::log(start.mole_fractions[
            phase * compiled.model->parameters().component_count() + spec.support[i]]);
    }
    std::size_t cursor = offsets[phase] + spec.support.size();
    if (spec.amount.kind == AmountKind::Free && start.amounts.size() == start.phase_count) {
        if (std::isfinite(start.amounts[phase]) && start.amounts[phase] > 0.0) {
            point[cursor] = start.amounts[phase];
        }
        ++cursor;
    }
    std::vector<double> composition(compiled.model->parameters().component_count(), 0.0);
    for (std::size_t species = 0; species < composition.size(); ++species)
        composition[species] = start.mole_fractions[phase * composition.size() + species];
    double density = start.molar_densities[phase];
    if (spec.amount.kind == AmountKind::Pinned && !start.temperature.empty()
        && !start.pressure.empty()) {
        const auto branch = spec.kind == PhaseKind::Vapor
            ? DensityBranch::Vapor : DensityBranch::Liquid;
        density = solve_density(*compiled.model, start.temperature[0], start.pressure[0],
                                composition, branch, std::nullopt).molar_density;
    }
    point[cursor] = std::log(compiled.model->packing_fraction_coefficient(
        start.temperature.empty() ? compiled.problem.temperature.value : start.temperature[0],
        composition) * density);
}

std::vector<double> start_point(
    const CompiledSystem& compiled, const EquilibriumResult& start,
    std::vector<double> point
) {
    const std::size_t components = compiled.model->parameters().component_count();
    if (start.phase_count != compiled.problem.phases.size()) return point;
    const auto offsets = compiled.phase_offsets();
    for (std::size_t phase = 0; phase < start.phase_count; ++phase) {
        if (valid_previous_phase(start, phase, components)) {
            copy_previous_phase(compiled, start, phase, point, offsets);
        }
    }
    std::size_t intensive = compiled.primal_size - compiled.free_intensive_count;
    if (compiled.problem.temperature.free && !start.temperature.empty()) point[intensive++] = start.temperature[0];
    if (compiled.problem.pressure.free && !start.pressure.empty()) point[intensive++] = start.pressure[0];
    if (compiled.pinned_count != 0 && start.multipliers.size() == compiled.rank + compiled.charge_row_count) {
        std::copy(start.multipliers.begin(), start.multipliers.end(), point.begin() + compiled.primal_size);
    }
    return point;
}

bool usable_start(const CompiledSystem& compiled, const EquilibriumResult& start) {
    if (!supports_split_seed(compiled) || !start.success
        || start.phase_count != compiled.problem.phases.size()
        || start.amounts.size() != start.phase_count) return false;
    const std::size_t components = compiled.model->parameters().component_count();
    for (std::size_t phase = 0; phase < start.phase_count; ++phase) {
        if (!valid_previous_phase(start, phase, components)
            || !(std::isfinite(start.amounts[phase]) && start.amounts[phase] > 0.0)) return false;
        std::vector<double> composition(components);
        for (std::size_t i = 0; i < components; ++i) {
            composition[i] = start.mole_fractions[phase * components + i];
            if (!(composition[i] > 0.0)) return false;
        }
        const double packing = compiled.model->packing_fraction_coefficient(
            start.temperature.empty() ? compiled.problem.temperature.value : start.temperature[0],
            composition) * start.molar_densities[phase];
        if (!(std::isfinite(packing) && packing >= policy::PACKING_FRACTION_FLOOR.value
              && packing <= policy::PACKING_FRACTION_CEILING.value)) return false;
    }
    return true;
}

double active_constraint_residual(
    const CompiledSystem& compiled, const std::vector<double>& values
) {
    double result = 0.0;
    for (std::size_t row = 0; row < values.size(); ++row) {
        if (compiled.pinned_count == 0 && compiled.rows[row].kind == RowKind::Chemical) continue;
        result = std::max(result, std::abs(values[row]));
    }
    return result;
}

struct ResidualSummary final {
    double maximum = 0.0;
    double intensive = 0.0;
    double balance = 0.0;
    double normalization = 0.0;
    double charge = 0.0;
};

ResidualSummary summarize_residuals(
    const std::vector<Row>& rows, const std::vector<double>& values
) {
    ResidualSummary result;
    const auto keep = [](double& slot, double value) {  // a nonfinite row stays infinite
        slot = std::isfinite(value) ? std::max(slot, std::abs(value)) : std::numeric_limits<double>::infinity();
    };
    for (std::size_t index = 0; index < std::min(rows.size(), values.size()); ++index) {
        const RowKind kind = rows[index].kind;
        keep(result.maximum, values[index]);
        if (kind == RowKind::Chemical || kind == RowKind::Pressure) keep(result.intensive, values[index]);
        else if (kind == RowKind::Moiety) keep(result.balance, values[index]);
        else if (kind == RowKind::Normalization) keep(result.normalization, values[index]);
        else if (kind == RowKind::Electro) keep(result.charge, values[index]);
    }
    return result;
}

}  // namespace

GibbsMinimizer::GibbsMinimizer(
    const CompiledSystem& compiled, std::vector<double> start, core::Scope* scope
)
    : compiled_(compiled), start_(std::move(start)), scope_(scope) {}

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
                if (!diagnose(r, x)) r.status = "numerical_failure";
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

bool GibbsMinimizer::get_nlp_info(
    Ipopt::Index& n, Ipopt::Index& m, Ipopt::Index& nnz_jac_g,
    Ipopt::Index& nnz_h_lag, IndexStyleEnum& index_style
) {
    n = static_cast<Ipopt::Index>(compiled_.pinned_count == 0
        ? compiled_.primal_size : compiled_.primal_size + compiled_.rank + compiled_.charge_row_count);
    m = static_cast<Ipopt::Index>(constraint_rows().size());
    nnz_jac_g = n * m;
    nnz_h_lag = n * (n + 1) / 2;
    index_style = C_STYLE;
    return true;
}

bool GibbsMinimizer::get_bounds_info(
    Ipopt::Index n, Ipopt::Number* x_l, Ipopt::Number* x_u,
    Ipopt::Index m, Ipopt::Number* g_l, Ipopt::Number* g_u
) {
    const double infinity = 1.0e19;
    std::fill(x_l, x_l + n, -infinity);
    std::fill(x_u, x_u + n, infinity);
    // Pinned square cells impose chemical equality and need the full physical support. Other
    // cells cap the box at a fraction of each species' attainable inventory u_i/N_feed, the LP
    // max y_i s.t. B y = b/N_feed, y >= 0 on the supports (B has mixed signs, so no row ratio
    // bounds y_i): the phases are never forced to hold a trace species that does not exist (#172).
    const std::size_t species = compiled_.n_feed.size();
    std::vector<double> box(species, policy::STABILITY_COMPOSITION_FLOOR.value);
    if (compiled_.pinned_count == 0 || !compiled_.is_square()) {
        for (std::size_t i = 0; i < species; ++i)
            box[i] = std::max(box[i], std::min(policy::EQUILIBRIUM_COMPOSITION_BOX.value,
                policy::EQUILIBRIUM_BOX_INVENTORY_FRACTION.value * compiled_.inventory_upper_bounds[i]));
    }
    std::size_t offset = 0;
    for (const auto& phase : compiled_.problem.phases) {
        for (const std::size_t i : phase.support) {
            x_l[offset] = std::log(box[i]);
            x_u[offset++] = 0.0;
        }
        if (phase.amount.kind == AmountKind::Free) x_l[offset++] = 0.0;
        x_l[offset] = std::log(policy::PACKING_FRACTION_FLOOR.value);
        x_u[offset++] = std::log(policy::PACKING_FRACTION_CEILING.value);
    }
    if (compiled_.problem.temperature.free) x_l[offset++] = std::numeric_limits<double>::epsilon();
    if (compiled_.problem.pressure.free) x_l[offset++] = std::numeric_limits<double>::epsilon();
    std::fill(g_l, g_l + m, 0.0);
    std::fill(g_u, g_u + m, 0.0);
    return true;
}

bool GibbsMinimizer::get_starting_point(
    Ipopt::Index n, bool init_x, Ipopt::Number* x,
    bool init_z, Ipopt::Number* z_L, Ipopt::Number* z_U,
    Ipopt::Index m, bool init_lambda, Ipopt::Number* lambda
) {
    if (!init_x) return false;
    const std::vector<double> point = compiled_.pinned_count == 0
        ? std::vector<double>(start_.begin(), start_.begin() + compiled_.primal_size)
        : start_;
    for (Ipopt::Index i = 0; i < n; ++i) x[i] = i < static_cast<Ipopt::Index>(point.size()) ? point[i] : 0.0;
    if (init_z) {
        std::fill(z_L, z_L + n, 0.0);
        std::fill(z_U, z_U + n, 0.0);
    }
    if (init_lambda) std::fill(lambda, lambda + m, 0.0);
    return true;
}

std::vector<std::size_t> GibbsMinimizer::constraint_rows() const {
    std::vector<std::size_t> result;
    for (std::size_t row = 0; row < compiled_.rows.size(); ++row) {
        if (compiled_.pinned_count != 0
            || (compiled_.rows[row].kind != RowKind::Chemical
                && compiled_.rows[row].kind != RowKind::Pressure)) {
            result.push_back(row);
        }
    }
    return result;
}

const PhasePoint& GibbsMinimizer::prepared(Ipopt::Index n, const Ipopt::Number* x) {
    std::vector<double> candidate(x, x + n);
    if (!prepared_.has_value() || prepared_->point != candidate) {
        try {
            prepared_ = compiled_.prepare(candidate, true);
        } catch (const std::exception& error) {
            prepared_.reset();
            rows_.reset();
            jacobian_.reset();
            gradient_.reset();
            objective_hessian_.reset();
            constraint_hessian_.reset();
            preparation_error_ = error.what();
            throw;
        }
        rows_.reset();
        jacobian_.reset();
        gradient_.reset();
        objective_hessian_.reset();
        constraint_hessian_.reset();
        ++work_.prepared_points;
    }
    return *prepared_;
}

const PointEvaluation& GibbsMinimizer::rows_at(Ipopt::Index n, const Ipopt::Number* x) {
    const auto& point = prepared(n, x);
    if (!rows_.has_value()) {
        rows_ = compiled_.evaluate_at(point);
        ++work_.row_evaluations;
    }
    return *rows_;
}

bool GibbsMinimizer::eval_f(
    Ipopt::Index n, const Ipopt::Number* x, bool /*new_x*/, Ipopt::Number& value
) {
    ++evaluations_;
    try {
        value = compiled_.pinned_count == 0 ? rows_at(n, x).objective : 0.0;
        return std::isfinite(value);
    } catch (const std::exception&) {
        return false;
    }
}

bool GibbsMinimizer::eval_grad_f(
    Ipopt::Index n, const Ipopt::Number* x, bool /*new_x*/, Ipopt::Number* grad_f
) {
    try {
        const auto& point = prepared(n, x);
        if (!gradient_.has_value()) {
            gradient_ = compiled_.objective_gradient_at(point);
            ++work_.gradients;
        }
        std::copy(gradient_->begin(), gradient_->end(), grad_f);
        return true;
    } catch (const std::exception&) {
        return false;
    }
}

bool GibbsMinimizer::eval_g(
    Ipopt::Index n, const Ipopt::Number* x, bool /*new_x*/,
    Ipopt::Index m, Ipopt::Number* g
) {
    try {
        const auto& values = rows_at(n, x).values;
        const auto rows = constraint_rows();
        for (Ipopt::Index row = 0; row < m; ++row) g[row] = values[rows[row]];
        return true;
    } catch (const std::exception&) {
        return false;
    }
}

bool GibbsMinimizer::eval_jac_g(
    Ipopt::Index n, const Ipopt::Number* x, bool /*new_x*/,
    Ipopt::Index m, Ipopt::Index nele_jac, Ipopt::Index* iRow,
    Ipopt::Index* jCol, Ipopt::Number* values
) {
    const auto rows = constraint_rows();
    if (values == nullptr) {
        Ipopt::Index entry = 0;
        for (Ipopt::Index row = 0; row < m; ++row) {
            for (Ipopt::Index column = 0; column < n; ++column) {
                iRow[entry] = row;
                jCol[entry++] = column;
            }
        }
        return entry == nele_jac;
    }
    try {
        const auto& point = prepared(n, x);
        if (!jacobian_.has_value()) {
            jacobian_ = compiled_.jacobian_at(point);
            ++work_.jacobians;
        }
        Ipopt::Index entry = 0;
        for (const std::size_t row : rows) {
            for (Ipopt::Index column = 0; column < n; ++column) {
                values[entry++] = (*jacobian_)[row * static_cast<std::size_t>(n) + static_cast<std::size_t>(column)];
            }
        }
        return entry == nele_jac;
    } catch (const std::exception&) {
        return false;
    }
}

bool GibbsMinimizer::eval_h(
    Ipopt::Index n, const Ipopt::Number* x, bool /*new_x*/,
    Ipopt::Number obj_factor, Ipopt::Index m, const Ipopt::Number* lambda,
    bool /*new_lambda*/, Ipopt::Index nele_hess, Ipopt::Index* iRow,
    Ipopt::Index* jCol, Ipopt::Number* values
) {
    ++hessian_calls_;
    if (values == nullptr) {
        Ipopt::Index entry = 0;
        for (Ipopt::Index row = 0; row < n; ++row) {
            for (Ipopt::Index column = 0; column <= row; ++column) {
                iRow[entry] = row;
                jCol[entry++] = column;
            }
        }
        return entry == nele_hess;
    }
    std::fill(values, values + nele_hess, 0.0);
    try {
        const auto& point = prepared(n, x);
        // Objective curvature is multiplier independent; the weighted constraint curvature is
        // recombined whenever the optimizer multipliers change at the same point.
        if (!objective_hessian_.has_value()) {
            objective_hessian_ = compiled_.hessian_at(point, {});
            ++work_.objective_hessians;
        }
        const std::vector<double> multipliers(lambda, lambda + m);
        if (!constraint_hessian_.has_value() || constraint_hessian_multipliers_ != multipliers) {
            constraint_hessian_ = compiled_.constraint_hessian_at(point, multipliers);
            constraint_hessian_multipliers_ = multipliers;
            ++work_.constraint_hessians;
        }
        Ipopt::Index entry = 0;
        for (Ipopt::Index row = 0; row < n; ++row) {
            for (Ipopt::Index column = 0; column <= row; ++column) {
                const std::size_t index = static_cast<std::size_t>(row) * n
                    + static_cast<std::size_t>(column);
                values[entry++] = obj_factor * (*objective_hessian_)[index]
                    + (*constraint_hessian_)[index];
            }
        }
        return entry == nele_hess;
    } catch (const std::exception&) {
        return false;
    }
}

bool GibbsMinimizer::intermediate_callback(
    Ipopt::AlgorithmMode /*mode*/, Ipopt::Index iter, Ipopt::Number /*obj_value*/,
    Ipopt::Number inf_pr, Ipopt::Number inf_du, Ipopt::Number mu,
    Ipopt::Number d_norm, Ipopt::Number regularization_size,
    Ipopt::Number alpha_du, Ipopt::Number alpha_pr, Ipopt::Index ls_trials,
    const Ipopt::IpoptData* /*ip_data*/, Ipopt::IpoptCalculatedQuantities* /*ip_cq*/
) {
    iterations_ = static_cast<std::size_t>(iter);
    regularization_size_ = regularization_size;
    if (scope_ != nullptr) scope_->probe("ipopt_iteration",
        std::vector<double>{inf_pr, inf_du, mu, d_norm, regularization_size,
                            alpha_pr, alpha_du, static_cast<double>(ls_trials)},
        static_cast<std::uint32_t>(iter));
    return true;
}

void GibbsMinimizer::finalize_solution(
    Ipopt::SolverReturn status, Ipopt::Index n, const Ipopt::Number* x,
    const Ipopt::Number* /*z_L*/, const Ipopt::Number* /*z_U*/, Ipopt::Index m,
    const Ipopt::Number* /*g*/, const Ipopt::Number* lambda,
    Ipopt::Number /*objective*/,
    const Ipopt::IpoptData* /*ip_data*/, Ipopt::IpoptCalculatedQuantities* /*ip_cq*/
) {
    status_ = status;
    solution_.assign(x, x + n);
    constraint_multipliers_.assign(lambda, lambda + m);
}

namespace {

EquilibriumResult structural_result(const CompiledSystem& compiled, core::Scope& scope) {
    EquilibriumResult result;
    result.model_identity = compiled.model;
    result.declaration = std::make_shared<const EquilibriumProblem>(compiled.problem);
    result.rows = compiled.rows;
    result.reference_diagnostics = compiled.reference_diagnostics;
    result.phase_count = compiled.problem.phases.size();
    result.unknown_count = compiled.unknown_count();
    result.jacobian_rank = compiled.jacobian_rank;
    result.generic_nullity = compiled.generic_nullity;
    if (!compiled.feasible) {
        // Structural rank loss at the generic point: no start, no optimizer, no candidate.
        result.solver_status = "NOT_RUN";
        result.message = "rank-deficient declared equations; nonlinear solve not run";
        result.convergence.residual = std::numeric_limits<double>::quiet_NaN();
        result.convergence.stop = core::Stop::infeasible;
        scope.probe("generic_nullity", static_cast<double>(compiled.generic_nullity));
        scope.finish(result.convergence);
    }
    return result;
}

// The one Ipopt lifecycle shared by the public seed route and the private coordinate route.
EquilibriumResult run_driver(
    const HelmholtzModel& model, const CompiledSystem& compiled,
    std::vector<double> initial, EquilibriumResult result, core::Scope& scope,
    std::vector<double>* attempt = nullptr, EquilibriumResult* audit = nullptr
) {
    const double row_tolerance = compiled.problem.row_tolerance;
    // The optimizer works one policy decade below the caller criterion, which alone classifies.
    // A pinned square system has a constant objective and is Newton-polished below, so it stops there.
    const double inner_tolerance = std::min(row_tolerance, policy::EQUILIBRIUM_ROW_RESIDUAL.value)
        * (compiled.pinned_count != 0 && compiled.is_square() ? 1.0 : policy::EQUILIBRIUM_INNER_TOLERANCE_FACTOR.value);
    Ipopt::SmartPtr<GibbsMinimizer> nlp = new GibbsMinimizer(compiled, std::move(initial), &scope);
    const auto [initialized, returned] = run_ipopt(nlp, inner_tolerance,
        static_cast<int>(policy::EQUILIBRIUM_MAX_ITERATIONS.value), 0.0, true, 1.0 / compiled.feed_total);
    const Ipopt::SolverReturn status = nlp->status();
    // ALG-TERMINATION: a failed initialization (status stays UNASSIGNED) or fatal return leaves no candidate.
    const bool fatal = initialized != Ipopt::Solve_Succeeded
        || (returned <= Ipopt::Invalid_Problem_Definition && returned != Ipopt::Invalid_Number_Detected);
    // The one accepted-outcome rule: only a non-error stop polishes, classifies its rows and publishes.
    const bool non_error = !fatal && (status == Ipopt::SUCCESS
        || status == Ipopt::STOP_AT_ACCEPTABLE_POINT || status == Ipopt::STOP_AT_TINY_STEP
        || status == Ipopt::FEASIBLE_POINT_FOUND || status == Ipopt::RESTORATION_FAILURE
        || status == Ipopt::MAXITER_EXCEEDED || status == Ipopt::CPUTIME_EXCEEDED);
    result.solver_status = solver_status_name(status);  // both exact names tell a solver failure from infeasibility
    if (fatal) result.solver_status += (initialized != Ipopt::Solve_Succeeded ? "; Initialize returned "
        : "; OptimizeTNLP returned ") + application_return_name(returned);
    const std::string ipopt_failure = fatal ? "; Ipopt " + result.solver_status : "";
    if (audit) {
        audit->solver_status = result.solver_status;
        audit->iterations = nlp->iterations(); audit->evaluations = nlp->evaluations();
        audit->hessian_calls = nlp->hessian_calls();
    }
    std::vector<double> point = nlp->solution();
    if (point.empty()) point = boundary_start(compiled);
    std::vector<double> full = compiled.pinned_count == 0 ? add_multipliers(compiled, point) : point;
    // The optimizer's amount-weighted rows underresolve intensive trace defects.
    // Exact Newton steps contract raw rows inside the physical composition cutoff.
    // Curvature may raise classified rows temporarily; return the best classified
    // point visited, including the optimizer start, rather than the last step.
    std::size_t polish_evaluations = 0;
    std::string rejected_preparation = nlp->preparation_error().value_or("");
    // A budget stop polishes too: an unattainable request can exhaust Ipopt after passing a
    // better point, and the polish never returns a worse classified point than its start.
    if (non_error && compiled.is_square()) {
        auto at = compiled.prepare(full, true);
        auto rows_before = compiled.evaluate_at(at).values;
        std::vector<double> best = full;
        double best_classified = maximum_abs(compiled.classified_residuals(at, rows_before));
        ++polish_evaluations;
        // Past the step budget only raw rows above the caller criterion continue, under a separate polish cap.
        for (std::size_t step = 0; step < static_cast<std::size_t>(policy::EQUILIBRIUM_MAX_ITERATIONS.value); ++step) {
            const double defect = maximum_abs(rows_before);
            if (step >= static_cast<std::size_t>(policy::EQUILIBRIUM_POLISH_STEPS.value) && defect <= row_tolerance) break;
            const double classified_before = maximum_abs(compiled.classified_residuals(at, rows_before));
            if (!(defect > 0.0) || !std::isfinite(defect)) break;
            const auto jacobian = compiled.jacobian_at(at);
            const auto n = static_cast<Eigen::Index>(full.size());
            if (jacobian.size() != full.size() * full.size() || rows_before.size() != full.size()) break;
            const Eigen::VectorXd delta = Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic,
                Eigen::RowMajor>>(jacobian.data(), n, n).fullPivLu().solve(
                    -Eigen::Map<const Eigen::VectorXd>(rows_before.data(), n));
            bool accepted = false;
            // Bounded damping retains raw contraction and the physical cutoff.
            // Best-point selection below protects the returned classified residual.
            for (const double alpha : {1.0, 0.5, 0.25}) {
                std::vector<double> candidate = full;
                for (Eigen::Index i = 0; i < n; ++i)
                    candidate[static_cast<std::size_t>(i)] += alpha * delta(i);
                PhasePoint at_candidate;
                std::vector<double> rows_after;
                try {
                    at_candidate = compiled.prepare(candidate, true);
                    rows_after = compiled.evaluate_at(at_candidate).values;
                    ++polish_evaluations;
                } catch (const std::exception& error) {
                    rejected_preparation = error.what();
                    continue;
                }
                bool inside = true;
                for (const auto& phase : at_candidate.composition)
                    for (const double x : phase)
                        inside = inside && (x == 0.0 || x >= policy::STABILITY_COMPOSITION_FLOOR.value);
                const double after = maximum_abs(rows_after);
                const double classified_after = maximum_abs(compiled.classified_residuals(at_candidate, rows_after));
                scope.probe("polish_step", {defect, after, classified_before, classified_after, alpha},
                    static_cast<std::uint32_t>(step));
                if (!inside || !(after < defect)) continue;
                full = std::move(candidate);
                at = std::move(at_candidate);
                rows_before = std::move(rows_after);
                if (classified_after <= best_classified) {
                    best = full;
                    best_classified = classified_after;
                }
                accepted = true;
                break;
            }
            if (!accepted) break;
        }
        full = std::move(best);
    }
    if (audit) audit->evaluations += polish_evaluations;
    PhasePoint prepared;
    try {
        prepared = compiled.prepare(full);
    } catch (const std::invalid_argument& error) {
        if (std::string(error.what()).rfind("reference_unavailable:", 0) != 0) throw;
        result.success = false;
        result.message = error.what() + ipopt_failure;
        result.convergence.stop = core::Stop::left_domain;
        result.convergence.residual = std::numeric_limits<double>::quiet_NaN();
        result.iterations = result.convergence.iterations = nlp->iterations();
        result.evaluations = result.convergence.evaluations = nlp->evaluations() + polish_evaluations;
        clear_state(result);
        scope.finish(result.convergence);
        return result;
    }
    if (attempt && !fatal && !nlp->solution().empty()) *attempt = full;
    const auto evaluation = compiled.evaluate_at(prepared);
    result.residuals = evaluation.values;
    result.classified_residuals = compiled.classified_residuals(prepared, evaluation.values);
    result.objective = evaluation.objective;
    const auto residuals = summarize_residuals(result.rows, result.residuals);
    scope.probe("maximum_residual", residuals.maximum);
    scope.probe("intensive_residual", residuals.intensive);
    scope.probe("balance_residual", residuals.balance);
    scope.probe("normalization_residual", residuals.normalization);
    scope.probe("charge_residual", residuals.charge);
    scope.probe("row_count", static_cast<double>(result.rows.size()));
    scope.probe("primal_size", static_cast<double>(compiled.primal_size));
    scope.probe("balance_rank", static_cast<double>(compiled.rank));
    scope.probe("pinned_count", static_cast<double>(compiled.pinned_count));
    scope.probe("active_residual", active_constraint_residual(compiled, evaluation.values));
    scope.probe("objective", result.objective);
    scope.probe("hessian_calls", static_cast<double>(nlp->hessian_calls()));
    scope.probe("prepared_points", static_cast<double>(nlp->work().prepared_points));
    core::Convergence& convergence = result.convergence;
    convergence = {};
    convergence.residual = maximum_abs(result.classified_residuals);
    if (row_tolerance == policy::EQUILIBRIUM_ROW_RESIDUAL.value) {
        convergence.criterion = &policy::EQUILIBRIUM_ROW_RESIDUAL;
    } else {
        // The caller-valued criterion lives with the record so copies outlive the request.
        convergence.criterion_owner = std::make_shared<const policy::Tolerance>(policy::Tolerance{
            "EQUILIBRIUM_REQUESTED_ROW_RESIDUAL", row_tolerance,
            policy::Category::AlgorithmTolerance, "caller-supplied equilibrium row tolerance"});
        convergence.criterion = convergence.criterion_owner.get();
    }
    if (evaluation.floor > 0.0) convergence.floor = evaluation.floor;
    convergence.amplification = nlp->regularization_size();
    result.iterations = convergence.iterations = nlp->iterations();
    result.evaluations = convergence.evaluations = nlp->evaluations() + polish_evaluations;
    result.hessian_calls = nlp->hessian_calls();
    switch (status) {
    case Ipopt::INVALID_NUMBER_DETECTED: convergence.stop = core::Stop::nonfinite; break;
    case Ipopt::DIVERGING_ITERATES: convergence.stop = core::Stop::left_domain; break;
    case Ipopt::ERROR_IN_STEP_COMPUTATION: convergence.stop = core::Stop::singular; break;
    case Ipopt::LOCAL_INFEASIBILITY:
    case Ipopt::TOO_FEW_DEGREES_OF_FREEDOM: convergence.stop = core::Stop::infeasible; break;
    default:
        // No candidate after a fatal return: infeasible like an unrun structural solve, whatever the start's rows.
        convergence.stop = !non_error ? core::Stop::infeasible : core::classify(convergence.residual,
            *convergence.criterion, convergence.floor, convergence.contraction,
            status == Ipopt::MAXITER_EXCEEDED || status == Ipopt::CPUTIME_EXCEEDED, convergence.amplification);
    }
    if (std::isfinite(convergence.residual)) {
        result.requested_tolerance_met = convergence.residual <= row_tolerance;
    }
    result.message = core::describe(convergence);
    scope.finish(convergence);
    const auto offsets = compiled.phase_offsets();
    std::vector<std::vector<double>> compositions;
    std::vector<double> amounts;
    std::vector<double> volumes;
    std::size_t intensive_offset = compiled.primal_size - compiled.free_intensive_count;
    const double candidate_temperature = compiled.problem.temperature.free
        ? full[intensive_offset] : compiled.problem.temperature.value;
    bool topology_ok = true;
    std::string topology_failure;
    const auto topology = check_topology(compiled, full);
    result.topology_event = topology.event;
    result.topology_phase = topology.phase;
    result.topology_pair = topology.pair;
    for (std::size_t phase = 0; phase < compiled.problem.phases.size(); ++phase) {
        const auto& spec = compiled.problem.phases[phase];
        std::vector<double> composition(model.parameters().component_count(), 0.0);
        double total = 0.0;
        for (std::size_t i = 0; i < spec.support.size(); ++i) {
            const double value = std::exp(full[offsets[phase] + i]);
            composition[spec.support[i]] = value;
            total += value;
        }
        if (!(total > 0.0)) topology_ok = false;
        for (double& value : composition) value /= std::max(total, std::numeric_limits<double>::min());
        compositions.push_back(composition);
        std::size_t cursor = offsets[phase] + spec.support.size();
        const double amount = spec.amount.kind == AmountKind::Free
            ? full[cursor++] : specified_amount(spec, compiled.n_feed);
        amounts.push_back(amount);
        volumes.push_back(model.packing_fraction_coefficient(
            candidate_temperature, composition) / std::exp(full[cursor]));
        const double deficit = michelsen_deficit(
            compiled, compositions[phase], phase, result.residuals);
        scope.probe("michelsen_deficit", deficit, static_cast<std::uint32_t>(phase));
        if (topology.event == TopologyEvent::PhaseDisappeared && topology.phase == phase) {
            topology_ok = false;
            topology_failure = "phase " + spec.name + " vanished; Michelsen deficit "
                + std::to_string(deficit);
        }
    }
    for (const double distance : topology.distances) scope.probe("phase_distance", distance);
    if (topology.event == TopologyEvent::PhasesCoalesced) {
        topology_ok = false;
        topology_failure = "phases_coincide: phase distance between "
            + compiled.problem.phases[topology.pair->first].name + " and "
            + compiled.problem.phases[topology.pair->second].name + " is too small: "
            + std::to_string(topology.distance);
    }
    // The recomputed rows classify: a criterion-meeting point is accepted after any non-error
    // stop (error stops never classify as converged); a floor-limited point needs success.
    const bool accepted = non_error && topology_ok && (result.convergence.stop == core::Stop::converged
        || (result.convergence.stop == core::Stop::converged_to_floor && status == Ipopt::SUCCESS));
    if (!accepted) {
        result.success = false;
        if (!topology_ok) result.message = topology_failure;
        result.message += ipopt_failure;
        if (!rejected_preparation.empty()) result.message += "; trial preparation failed: " + rejected_preparation;
        clear_state(result);
        return result;
    }
    if (result.requested_tolerance_met.has_value() && !*result.requested_tolerance_met) {
        result.message += "; requested row tolerance " + std::to_string(row_tolerance)
            + " not met by an available floor-limited state";
    }
    const double temperature = compiled.problem.temperature.free
        ? full[intensive_offset++] : compiled.problem.temperature.value;
    const double pressure = compiled.problem.pressure.free
        ? full[intensive_offset++] : compiled.problem.pressure.value;
    result.success = true;
    result.temperature = {temperature};
    result.pressure = {pressure};
    if (!prepared.standard.reference_diagnostics.empty()) {
        result.reference_diagnostics = prepared.standard.reference_diagnostics;
    }
    result.amounts = amounts;
    result.molar_densities.reserve(volumes.size());
    result.packing_fractions.reserve(volumes.size());
    for (std::size_t phase = 0; phase < volumes.size(); ++phase) {
        result.molar_densities.push_back(1.0 / volumes[phase]);
        result.packing_fractions.push_back(
            model.packing_fraction_coefficient(temperature, compositions[phase]) / volumes[phase]);
        result.mole_fractions.insert(result.mole_fractions.end(),
                                     compositions[phase].begin(), compositions[phase].end());
    }
    result.multipliers.assign(full.begin() + static_cast<std::ptrdiff_t>(compiled.primal_size), full.end());
    // Extensive-objective constraint multipliers at the returned (polished) point, from the row
    // multipliers: a material row (B n - b)/N_feed carries -N_feed lambda_j, a charge row sum z x
    // carries -n_p phi_p, and a free-amount phase's normalization row is zero because the
    // objective and material rows obey joint extensive homogeneity in amount, composition
    // and volume, together with stationarity in the free amount. The raw
    // optimizer values stay a trace probe.
    scope.probe("optimizer_constraint_multipliers", nlp->constraint_multipliers());
    if (compiled.pinned_count == 0) {
        result.constraint_multipliers.clear();
        std::size_t electro = 0;
        for (const auto& row : compiled.rows) {
            if (row.kind == RowKind::Moiety) {
                result.constraint_multipliers.push_back(-compiled.feed_total * full[compiled.primal_size + row.index]);
            } else if (row.kind == RowKind::Electro) {
                result.constraint_multipliers.push_back(
                    -prepared.amount[row.phase] * full[compiled.primal_size + compiled.rank + electro++]);
            } else if (row.kind != RowKind::Chemical && row.kind != RowKind::Pressure) {
                result.constraint_multipliers.push_back(0.0);
            }
        }
    } else {
        result.constraint_multipliers = nlp->constraint_multipliers();
    }
    result.coordinates = full;
    const auto curvature = reduced_curvature_at(compiled, std::move(prepared), result.constraint_multipliers);
    result.reduced_hessian_minimum = curvature.minimum;
    scope.probe("reduced_hessian_equilibrated", curvature.equilibrated);
    scope.probe("reduced_hessian_band", curvature.band);
    if (compiled.pinned_count == 0 && !(curvature.equilibrated > curvature.band)) {
        // Inside the band the inertia is unresolved at working precision: a typed singular outcome, not a noise sign.
        const bool singular = curvature.equilibrated >= -curvature.band;
        if (singular || std::isnan(curvature.equilibrated))
            result.convergence.stop = singular ? core::Stop::singular : core::Stop::nonfinite;
        result.success = false;
        result.message = std::string("reduced_hessian_stability ") + (singular ? "is indeterminate" : "is not positive definite")
            + ": equilibrated minimum " + core::optional_number(curvature.equilibrated)
            + ", rounding band " + core::optional_number(curvature.band);
        if (!rejected_preparation.empty()) {
            result.message += "; trial preparation failed: " + rejected_preparation;
        }
        clear_state(result);
        return result;
    }
    result.rows = compiled.rows;
    require_solution_shape(result, compiled.problem.phases.size(),
                           model.parameters().component_count());
    scope.probe("reduced_hessian_minimum", result.reduced_hessian_minimum);
    scope.probe("phase_count", static_cast<double>(result.phase_count));
    scope.probe("residual", result.convergence.residual);
    scope.probe("iterations", static_cast<double>(result.iterations));
    scope.probe("evaluations", static_cast<double>(result.evaluations));
    scope.probe("rows", static_cast<double>(result.rows.size()));
    scope.probe("multipliers", static_cast<double>(result.multipliers.size()));
    scope.probe("amount_total", std::accumulate(result.amounts.begin(), result.amounts.end(), 0.0));
    scope.probe("topology_ok", topology_ok ? 1.0 : 0.0);
    scope.probe("feasible", compiled.feasible ? 1.0 : 0.0);
    scope.probe("square", compiled.is_square() ? 1.0 : 0.0);
    return result;
}

bool zero_amount(const PhaseSpec& phase) {
    return phase.amount.kind == AmountKind::Pinned && phase.amount.fraction == 0.0;
}

// A cold reactive boundary at fixed T and free P: the phases present (nonzero amount) alone, in reaction
// equilibrium at the pressure seed, are the solution's present phases up to their pressure dependence. Their
// compositions replace the caller's guesses before the ordinary boundary start, and an unguessed zero-amount
// phase starts from the first present phase on its support (#152); none if any of it fails.
std::optional<std::vector<double>> speciated_start(const CompiledSystem& compiled) {
    const auto& problem = compiled.problem;
    EquilibriumProblem present = problem;
    present.phases.erase(std::remove_if(present.phases.begin(), present.phases.end(), zero_amount), present.phases.end());
    if (!compiled.has_reaction || !problem.pressure.free || problem.temperature.free || present.phases.empty()
        || present.phases.size() == problem.phases.size()) return std::nullopt;
    present.pressure = FreeValue::fixed(problem.pressure.seed_or(1.0e5));
    try {
        const auto solved = solve_equilibrium(*compiled.model, present);
        if (!solved.success) return std::nullopt;
        const std::size_t components = compiled.model->parameters().component_count();
        CompiledSystem seeded = compiled;
        std::size_t index = 0;
        for (auto& phase : seeded.problem.phases) {
            const bool guessed = phase.composition_guess.size() == phase.support.size();
            if (zero_amount(phase) && guessed) continue;
            const std::size_t source = zero_amount(phase) ? 0 : index++;
            phase.composition_guess.clear();
            for (const std::size_t i : phase.support) phase.composition_guess.push_back(solved.mole_fractions[source * components + i]);
        }
        return boundary_start(seeded);
    } catch (const std::exception&) {
        return std::nullopt;
    }
}

}  // namespace

TopologyCheck check_topology(const CompiledSystem& compiled, const std::vector<double>& full) {
    TopologyCheck result;
    if (full.size() < compiled.primal_size) return result;
    const std::size_t components = compiled.model->parameters().component_count();
    const auto offsets = compiled.phase_offsets();
    std::vector<std::vector<double>> compositions;
    std::vector<std::size_t> eta_index;
    for (std::size_t phase = 0; phase < compiled.problem.phases.size(); ++phase) {
        const auto& spec = compiled.problem.phases[phase];
        std::vector<double> composition(components, 0.0);
        double total = 0.0;
        for (std::size_t i = 0; i < spec.support.size(); ++i) {
            composition[spec.support[i]] = std::exp(full[offsets[phase] + i]);
            total += composition[spec.support[i]];
        }
        for (double& value : composition) value /= std::max(total, std::numeric_limits<double>::min());
        std::size_t cursor = offsets[phase] + spec.support.size();
        const double amount = spec.amount.kind == AmountKind::Free
            ? full[cursor++] : specified_amount(spec, compiled.n_feed);
        eta_index.push_back(cursor);
        compositions.push_back(std::move(composition));
        if (spec.amount.kind == AmountKind::Free
            && amount / compiled.feed_total <= policy::STABILITY_PHASE_AMOUNT_FLOOR.value
            && result.event == TopologyEvent::None) {
            result.event = TopologyEvent::PhaseDisappeared;
            result.phase = phase;
        }
    }
    for (std::size_t left = 0; left + 1 < compositions.size(); ++left) {
        for (std::size_t right = left + 1; right < compositions.size(); ++right) {
            const double distance = std::max(
                std::abs(full[eta_index[left]] - full[eta_index[right]]),
                largest_difference(compositions[left], compositions[right]));
            result.distances.push_back(distance);
            if (distance <= policy::PHASE_COINCIDENCE_DISTANCE.value
                && result.event == TopologyEvent::None) {
                result.event = TopologyEvent::PhasesCoalesced;
                result.pair = std::make_pair(left, right);
                result.distance = distance;
            }
        }
    }
    return result;
}

EquilibriumResult solve_equilibrium(const HelmholtzModel& model, const EquilibriumProblem& request,
                                    const std::optional<EquilibriumResult>& start) {
    core::Scope scope("solve_equilibrium");  // compilation belongs to the solve's span
    if (request.discovery) return discover_phases(model, request);
    return solve_compiled(compile_problem(model, request), start, &scope);
}

EquilibriumResult solve_compiled(const CompiledSystem& compiled, const std::optional<EquilibriumResult>& start,
                                 core::Scope* outer) {
    std::optional<core::Scope> own;
    auto& scope = outer != nullptr ? *outer : own.emplace("solve_equilibrium");
    const auto& model = *compiled.model;
    const EquilibriumResult structural = structural_result(compiled, scope);
    if (!compiled.feasible) return structural;
    const bool use_start = start.has_value() && usable_start(compiled, *start);
    // A previous result that start_point copies replaces the start's composition; skip its speciation.
    const bool copied = start.has_value() && start->phase_count == compiled.problem.phases.size()
        && valid_previous_phase(*start, 0, model.parameters().component_count());
    auto seeded = start.has_value() ? std::nullopt : speciated_start(compiled);
    EquilibriumResult result;
    // The speciated seed when there is one; the caller's start alone or after the seed's failure.
    for (int attempt = seeded ? 0 : 1; attempt < 2 && !result.success; ++attempt) {
        std::vector<double> initial;
        if (attempt == 0) {
            initial = std::move(*seeded);
        } else {
            initial = use_start ? compiled.generic_point() : boundary_start(compiled, !copied);
            if (start.has_value()) initial = start_point(compiled, *start, std::move(initial));
        }
        try {
            result = run_driver(model, compiled, std::move(initial), structural, scope);
        } catch (const std::exception&) {
            if (attempt == 0) continue;  // the seeded attempt's exception falls back to the declared start
            throw;
        }
        if (!result.success) continue;
        BranchCheck branches;
        try {
            branches = check_density_branches(model, compiled, result, {}, true, scope);
        } catch (const std::exception& error) {
            branches.message = std::string("density branch admission failed: ") + error.what();
        }
        if (!branches.accepted) {
            result.success = false;
            result.message = branches.message;
            clear_state(result);
        }
    }
    return result;
}

EquilibriumResult solve_equilibrium_from_coordinates(
    const HelmholtzModel& model,
    const EquilibriumProblem& request,
    const std::vector<double>& coordinates, std::vector<double>* attempt
) {
    core::Scope scope("solve_equilibrium");
    CompiledSystem compiled = compile_problem(model, request);
    EquilibriumResult result = structural_result(compiled, scope);
    if (!compiled.feasible) return result;
    if (coordinates.size() != compiled.unknown_count() || !all_numeric(coordinates)) {
        throw std::invalid_argument("continuation coordinates do not match the compiled layout");
    }
    return run_driver(model, compiled, coordinates, std::move(result), scope, attempt);
}

double present_phase_residual_distance(double chemical_contraction, double pressure_row) {
    return chemical_contraction - pressure_row;
}
bool descent_resolved(double tpd, double slope, double error_estimate) {
    return std::isfinite(tpd) && std::isfinite(slope) && std::isfinite(error_estimate)
        && error_estimate >= 0.0 && tpd < 0.0 && slope < -error_estimate;
}
namespace {
using Quantities = std::map<std::string, double>;

// Explicit searches alone capture the shared owners' work, including initialization,
// reference roots, failed roots, branch admission and compiler evaluations.
void count_stability_work(Quantities& work, const core::Trace& trace, std::size_t begin = 0) {
    std::size_t roots = 0, density_evaluations = 0, energy_evaluations = 0;
    for (std::size_t i = begin; i < trace.spans.size(); ++i) {
        const auto& span = trace.spans[i];
        if (std::string(span.name) == "solve_density") ++roots;
        if (std::string(span.name) != "evaluate") continue;
        ++energy_evaluations;
        auto parent = span.parent;
        while (parent) {
            if (std::string(trace.spans[*parent].name) == "solve_density") { ++density_evaluations; break; }
            parent = trace.spans[*parent].parent;
        }
    }
    work["root_calls"] = roots; work["root_evaluations"] = density_evaluations;
    work["energy_evaluations"] = energy_evaluations;
}

void count_driver(Quantities& work, const EquilibriumResult& result, const char* role) {
    ++work["driver_invocations"]; ++work[role];
    work["iterations"] += result.iterations;
    work["evaluations"] += result.evaluations;
    work["hessian_calls"] += result.hessian_calls;
}

DensityRoot stability_root(const HelmholtzModel& model, double T, double P,
    const std::vector<double>& w, PhaseKind kind, Quantities& work) {
    ++work["root_calls"];
    const auto root = solve_density(model, T, P, w, kind == PhaseKind::Liquid
        ? DensityBranch::Liquid : DensityBranch::Vapor, std::nullopt);
    work["root_evaluations"] += root.convergence.evaluations;
    return root;
}

void packing_admission(double eta) {
    if (!(std::isfinite(eta) && eta >= policy::PACKING_FRACTION_FLOOR.value
          && eta <= policy::PACKING_FRACTION_CEILING.value))
        throw std::invalid_argument("stability point outside packing domain");
}

// Reduced energy and Euler/association estimates use unit-mol amount columns and the SI volume column.
Quantities phase_energy(const HelmholtzModel& model, double T, double P, double v,
    const std::vector<double>& w, const std::vector<double>& standard, Quantities& work,
    HelmholtzDerivatives* retained = nullptr) {
    packing_admission(model.packing_fraction_coefficient(T, w) / v);
    ++work["energy_evaluations"];
    auto d = model.evaluate(T, v, w, DerivativeOrder::IsothermalPhaseHessian);
    const double rt = GAS_CONSTANT_J_PER_MOL_K * T;
    const double ideal = ideal_mechanical_helmholtz(w, v, 1.0);
    double euler = v * d.gradient[1], magnitude = std::abs(d.value);
    for (std::size_t i = 0; i < w.size(); ++i) if (w[i] > 0.0) {
        euler += w[i] * d.gradient[2 + i];
        magnitude += std::abs(w[i] * (std::log(w[i] / v) - 1.0 + standard[i]));
    }
    const double pv = P * v / rt;
    const double value = d.value + ideal + pv
        + std::inner_product(w.begin(), w.end(), standard.begin(), 0.0);
    // gamma covers final summation only, not the EOS. Euler defect and inner error are separate.
    const double m = 2.0 * w.size() + 4.0;
    const double gamma = m * std::numeric_limits<double>::epsilon()
        / (1.0 - m * std::numeric_limits<double>::epsilon());
    Quantities q{{"g_over_rt", value}, {"volume_m3_per_mol", v},
        {"pressure_row", v * (1.0 / v - d.gradient[1] - P / rt)},
        {"gradient_floor_max_si_columns", d.gradient_floor},
        {"association_energy_estimate_over_rt", (1.0 + std::abs(v)) * d.gradient_floor},
        {"euler_defect_over_rt", std::abs(d.value - euler)},
        {"summation_estimate_over_rt", gamma * (magnitude + std::abs(pv))}};
    if (!model.parameters().site_component.empty()) {
        ++work["association_diagnostic_evaluations"];
        const auto diagnostics = model.density_actions(T, 1.0 / v, w,
            {PropertyObservable::ResidualEntropy}, {DensityActionDirection{0.0, 0.0,
                std::vector<double>(w.size(), 0.0), {}}}).diagnostics;
        q["association_residual"] = diagnostics.association_residual;
        q["association_site_pivot_ratio"] = diagnostics.site_pivot_ratio;
        q["association_site_rank"] = static_cast<double>(diagnostics.site_rank);
    }
    if (retained) *retained = std::move(d);
    return q;
}

double energy_estimate(const Quantities& q) {
    return q.at("association_energy_estimate_over_rt") + q.at("euler_defect_over_rt")
        + q.at("summation_estimate_over_rt");
}

std::vector<PhaseSpec> stability_supports(const CompiledSystem& c, const StabilityOptions& options) {
    const auto& z = c.model->parameters().charge_number;
    std::vector<std::size_t> united, zero_charge;
    for (std::size_t i = 0; i < c.n_feed.size(); ++i) {
        bool participant = c.n_feed[i] > 0.0;
        for (const auto& reaction : c.problem.reactions) participant |= reaction.stoichiometry[i] != 0.0;
        if (participant) { united.push_back(i); if (z[i] == 0.0) zero_charge.push_back(i); }
    }
    std::vector<PhaseSpec> phases;
    if (options.trials) phases = *options.trials;
    else {
        PhaseSpec liquid; liquid.name = "trial-liquid"; liquid.support = united;
        phases.push_back(liquid);
        if (!zero_charge.empty()) {
            liquid.support = zero_charge; phases.push_back(liquid);
            liquid.name = "trial-vapor"; liquid.kind = PhaseKind::Vapor; phases.push_back(liquid);
        }
    }
    if (phases.empty()) throw std::invalid_argument("empty stability trial family");
    std::vector<PhaseSpec> admitted;
    for (auto phase : phases) {
        if (phase.support.empty()) throw std::invalid_argument("empty stability support");
        std::sort(phase.support.begin(), phase.support.end());
        if (std::adjacent_find(phase.support.begin(), phase.support.end()) != phase.support.end())
            throw std::invalid_argument("duplicate stability support species");
        bool positive = false, negative = false;
        for (auto i : phase.support) {
            if (std::find(united.begin(), united.end(), i) == united.end())
                throw std::invalid_argument("stability support outside fed/reaction union");
            positive |= z[i] > 0; negative |= z[i] < 0;
        }
        if (positive != negative || (positive && phase.kind == PhaseKind::Vapor))
            throw std::invalid_argument("inadmissible ionic stability support/label");
        const auto& form = c.model->parameters().permittivity_form;
        if (positive && !form.is(PermittivityForm::DensityTemperature)
            && !form.is(PermittivityForm::DipolarPerturbation))
            throw std::invalid_argument("permittivity model has no low-density limit; declare phases");
        phase.support_declared = true; phase.amount = AmountSpec::free();
        if (std::none_of(admitted.begin(), admitted.end(), [&](const auto& p) {
                return p.kind == phase.kind && p.support == phase.support; })) admitted.push_back(std::move(phase));
    }
    return admitted;
}

struct Candidate final {
    PhasePoint physical;
    std::vector<Quantities> quantities;
    double correction_norm = 0.0, feasibility_estimate = 0.0;
};

// Equality repair acts on supported amounts. Redundant charge/balance rows are solved together;
// the achieved defect, positivity and structural support decide admission, not a row-count guess.
Candidate stability_candidate(const CompiledSystem& c, const EquilibriumResult& result,
    StabilityReport& report) {
    if (!result.has_state() || result.model_identity != c.model || !result.declaration)
        throw std::invalid_argument("stability test requires this model's accepted declared result");
    const auto& original = *result.declaration;
    if (original.phases.size() != c.problem.phases.size() || original.feed.kind != c.problem.feed.kind
        || original.feed.values != c.problem.feed.values || original.reactions.size() != c.problem.reactions.size()
        || original.standard_gibbs_j_per_mol != c.problem.standard_gibbs_j_per_mol
        || original.temperature.free != c.problem.temperature.free || original.pressure.free != c.problem.pressure.free
        || (!original.temperature.free && original.temperature.value != c.problem.temperature.value)
        || (!original.pressure.free && original.pressure.value != c.problem.pressure.value))
        throw std::invalid_argument("stability result/declaration mismatch");
    for (std::size_t p = 0; p < original.phases.size(); ++p) {
        const auto& a = original.phases[p]; const auto& b = c.problem.phases[p];
        if (a.support != b.support || a.kind != b.kind || a.amount.kind != b.amount.kind
            || a.amount.fraction != b.amount.fraction) throw std::invalid_argument("stability phase declaration mismatch");
    }
    for (std::size_t r = 0; r < original.reactions.size(); ++r)
        if (!same_reaction(original.reactions[r], c.problem.reactions[r]))
            throw std::invalid_argument("stability reaction declaration mismatch");
    Candidate candidate;
    auto& at = candidate.physical;
    at = c.prepare(result.coordinates);
    report.work["energy_evaluations"] += at.composition.size();
    const auto raw = c.evaluate_at(at);
    if (raw.values != result.residuals) throw std::invalid_argument("stability result/input identity mismatch");
    report.standard_over_rt = at.standard.derivative[0];
    report.plane_over_rt.assign(c.n_feed.size(), 0.0);
    for (std::size_t k = 0; k < c.rank; ++k)
        for (std::size_t i = 0; i < c.n_feed.size(); ++i)
            report.plane_over_rt[i] += c.B[k][i] * result.multipliers[k];
    std::vector<std::pair<std::size_t, std::size_t>> cells;
    for (std::size_t p = 0; p < at.amount.size(); ++p) if (at.amount[p] > 0.0)
        for (auto i : c.problem.phases[p].support) cells.emplace_back(p, i);
    Eigen::MatrixXd J = Eigen::MatrixXd::Zero(c.rank + at.amount.size(), cells.size());
    Eigen::VectorXd rhs = Eigen::VectorXd::Zero(J.rows());
    for (std::size_t k = 0; k < c.rank; ++k) rhs(k) = c.balance[k];
    const auto& z = c.model->parameters().charge_number;
    for (std::size_t j = 0; j < cells.size(); ++j) {
        const auto [p, i] = cells[j];
        for (std::size_t k = 0; k < c.rank; ++k) J(k, j) = c.B[k][i];
        J(c.rank + p, j) = z[i];
        rhs -= J.col(j) * at.amount[p] * at.composition[p][i];
    }
    const Eigen::VectorXd delta = J.completeOrthogonalDecomposition().solve(rhs);
    if (!delta.allFinite() || (J * delta - rhs).lpNorm<Eigen::Infinity>() > policy::ADDITION_MATERIAL_BALANCE.value * std::max(1.0, c.feed_total))
        throw std::invalid_argument("stability amount correction rank/consistency failure");
    std::vector<std::vector<double>> inventory(at.amount.size(), std::vector<double>(c.n_feed.size(), 0.0));
    for (std::size_t j = 0; j < cells.size(); ++j) {
        const auto [p, i] = cells[j];
        inventory[p][i] = at.amount[p] * at.composition[p][i] + delta(j);
        if (!(inventory[p][i] > 0.0)) throw std::invalid_argument("stability correction destroys positive support");
        const double mu = std::log(at.composition[p][i] / at.volume[p])
            + at.derivative[p].gradient[2 + i] + report.standard_over_rt[i];
        candidate.feasibility_estimate += std::abs(delta(j) * mu);
        candidate.correction_norm += std::abs(delta(j));
    }
    report.candidate_coordinates.clear();
    for (std::size_t p = 0; p < at.amount.size(); ++p) {
        if (at.amount[p] > 0.0) {
            at.amount[p] = std::accumulate(inventory[p].begin(), inventory[p].end(), 0.0);
            for (auto i : c.problem.phases[p].support) at.composition[p][i] = inventory[p][i] / at.amount[p];
        }
        auto q = phase_energy(*c.model, at.temperature, at.pressure, at.volume[p], at.composition[p],
            report.standard_over_rt, report.work, &at.derivative[p]);
        double normalization = -1.0, charge = 0.0, potential = 0.0, residual = 0.0;
        for (auto i : c.problem.phases[p].support) {
            normalization += at.composition[p][i]; charge += z[i] * at.composition[p][i];
            potential += at.composition[p][i] * report.plane_over_rt[i];
            residual += at.composition[p][i] * (std::log(at.composition[p][i] / at.volume[p])
                + at.derivative[p].gradient[2 + i] + report.standard_over_rt[i] - report.plane_over_rt[i]);
            report.candidate_coordinates.push_back(std::log(at.composition[p][i]));
        }
        if (std::abs(normalization) > policy::COMPOSITION_SUM_TOLERANCE.value || std::abs(charge) > policy::COMPOSITION_SUM_TOLERANCE.value)
            throw std::invalid_argument("stability candidate arithmetic feasibility failure");
        report.candidate_coordinates.push_back(at.amount[p]);
        report.candidate_coordinates.push_back(std::log(c.model->packing_fraction_coefficient(at.temperature,
            at.composition[p]) / at.volume[p]));
        q["tpd_over_rt"] = q.at("g_over_rt") - potential;
        q["chemical_plus_electric_contraction"] = residual;
        q["present_identity_defect_over_rt"] = q["tpd_over_rt"] - present_phase_residual_distance(residual, q.at("pressure_row"));
        q["normalization_residual"] = normalization; q["charge_residual"] = charge;
        q["correction_l1_mol"] = candidate.correction_norm;
        q["feasibility_estimate_over_rt"] = candidate.feasibility_estimate;
        candidate.quantities.push_back(std::move(q));
    }
    report.present_phases = candidate.quantities;
    return candidate;
}

std::vector<double> trial_coordinates(const CompiledSystem& c, const std::vector<double>& w, double eta) {
    auto point = c.generic_point();
    const auto offsets = c.phase_offsets();
    for (std::size_t i = 0; i < c.problem.phases[0].support.size(); ++i)
        point[offsets[0] + i] = std::log(w[c.problem.phases[0].support[i]]);
    point[offsets[0] + c.problem.phases[0].support.size()] = 1.0;
    point[offsets[0] + c.problem.phases[0].support.size() + 1] = std::log(eta);
    return point;
}
} // namespace

// Charge rescaling followed by a two-pivot arithmetic correction (largest entry, ties by index).
std::vector<double> reconstruct_trial_composition(std::vector<double> w,
    const std::vector<std::size_t>& support, const std::vector<int>& z) {
    if (support.empty() || z.size() != w.size()) throw std::invalid_argument("invalid reconstruction support");
    double positive = 0.0, negative = 0.0;
    std::optional<std::size_t> plus, minus;
    for (auto i : support) {
        if (i >= w.size() || !(std::isfinite(w[i]) && w[i] > 0.0))
            throw std::invalid_argument("nonpositive/nonfinite trial composition");
        if (z[i] > 0.0) { positive += z[i] * w[i]; if (!plus || w[i] > w[*plus]) plus = i; }
        if (z[i] < 0.0) { negative -= z[i] * w[i]; if (!minus || w[i] > w[*minus]) minus = i; }
    }
    if (plus.has_value() != minus.has_value()) throw std::invalid_argument("one-sign trial support");
    if (plus) {
        const double ratio = std::sqrt(negative / positive);
        for (auto i : support) { if (z[i] > 0) w[i] *= ratio; if (z[i] < 0) w[i] /= ratio; }
    }
    double total = 0.0;
    for (auto i : support) total += w[i];
    if (!(std::isfinite(total) && total > 0.0)) throw std::invalid_argument("trial normalization failed");
    for (auto i : support) w[i] /= total;
    if (plus) {
        double sum = 0.0, charge = 0.0;
        for (auto i : support) { sum += w[i]; charge += z[i] * w[i]; }
        const double dp = (-charge + z[*minus] * (sum - 1.0)) / (z[*plus] - z[*minus]);
        w[*plus] += dp; w[*minus] += 1.0 - sum - dp;
    } else {
        const auto pivot = *std::max_element(support.begin(), support.end(), [&](auto a, auto b) { return w[a] < w[b]; });
        double sum = 0.0; for (auto i : support) sum += w[i]; w[pivot] += 1.0 - sum;
    }
    double sum = 0.0, charge = 0.0;
    for (auto i : support) {
        if (!(std::isfinite(w[i]) && w[i] > 0.0)) throw std::invalid_argument("reconstruction destroys support");
        sum += w[i]; charge += z[i] * w[i];
    }
    if (std::abs(sum - 1.0) > policy::COMPOSITION_SUM_TOLERANCE.value || std::abs(charge) > policy::COMPOSITION_SUM_TOLERANCE.value)
        throw std::invalid_argument("trial reconstruction arithmetic residual");
    return w;
}
namespace {
void donor_descent(const CompiledSystem& c, const Candidate& candidate, StabilityTrial& trial,
    StabilityReport& report, const StabilityOptions& options) {
    const auto& at = candidate.physical;
    const auto& w = trial.composition;
    auto& q = trial.quantities;
    double capacity = 0.0;
    for (std::size_t p = 0; p < at.amount.size(); ++p) {
        if (!(at.amount[p] > 0.0)) continue;
        double beta = std::numeric_limits<double>::infinity();
        for (auto i : trial.phase.support) {
            if (at.composition[p][i] <= 0.0) { beta = 0.0; break; }
            beta = std::min(beta, at.amount[p] * at.composition[p][i] / w[i]);
        }
        if (std::isfinite(beta) && beta > capacity) { capacity = beta; trial.donor = p; }
    }
    if (!trial.donor) { trial.evidence = "negative_plane_without_positive_donor"; return; }
    const auto p = *trial.donor;
    const auto& donor = candidate.quantities[p];
    double chemical = 0.0, magnitude = std::abs(q.at("g_over_rt"));
    double standard_error = 0.0;
    for (auto i : trial.phase.support) {
        const double term = w[i] * (std::log(at.composition[p][i] / at.volume[p])
            + at.derivative[p].gradient[2 + i] + report.standard_over_rt[i]);
        chemical += term; magnitude += std::abs(term);
        if (!options.standard_uncertainty_over_rt.empty())
            standard_error += 2.0 * w[i] * options.standard_uncertainty_over_rt[i];
    }
    const double slope = q.at("g_over_rt") - chemical + donor.at("pressure_row");
    const double gamma = (w.size() + 3.0) * std::numeric_limits<double>::epsilon();
    const double arithmetic = gamma * (magnitude + std::abs(donor.at("pressure_row"))) / (1.0 - gamma);
    const double uncertainty = energy_estimate(q) + (1.0 + std::abs(at.volume[p]))
        * at.derivative[p].gradient_floor + candidate.feasibility_estimate + arithmetic + standard_error;
    q["donor_capacity_mol"] = capacity;
    q["donor_chemical_contraction_over_rt"] = chemical;
    q["donor_pressure_row"] = donor.at("pressure_row");
    q["donor_correction_over_rt"] = slope - q.at("tpd_over_rt");
    q["donor_slope_over_rt"] = slope;
    q["descent_error_estimate_over_rt"] = uncertainty;
    q["standard_error_estimate_over_rt"] = standard_error;
    q["feasibility_error_estimate_over_rt"] = candidate.feasibility_estimate;
    q["arithmetic_error_estimate_over_rt"] = arithmetic;
    trial.evidence = "indeterminate";
    if (!descent_resolved(q.at("tpd_over_rt"), slope, uncertainty)) return;
    const auto& z = c.model->parameters().charge_number;
    double beta = capacity / 2.0;
    for (std::size_t check = 0; check <= 8; ++check, beta /= 2.0) {
        Quantities split{{"beta_mol", beta}, {"admitted", 0.0}};
        try {
            const double remaining = at.amount[p] - beta;
            std::vector<double> x(w.size(), 0.0), error(w.size(), 0.0);
            double sum = 0.0, charge = 0.0, inventory_error = 0.0;
            for (auto i : c.problem.phases[p].support) {
                const double n = at.amount[p] * at.composition[p][i] - beta * w[i];
                if (!(n > 0.0 && remaining > 0.0)) throw std::invalid_argument("split destroys donor support");
                x[i] = n / remaining; sum += x[i]; charge += z[i] * x[i];
                error[i] = remaining * x[i] + beta * w[i] - at.amount[p] * at.composition[p][i];
                inventory_error = std::max(inventory_error, std::abs(error[i]));
            }
            const auto moiety_error = matvec(c.B, error);
            split["component_error_mol"] = inventory_error;
            split["moiety_error_mol"] = maximum_abs(moiety_error);
            split["normalization_residual"] = sum - 1.0; split["charge_residual"] = charge;
            if (std::abs(sum - 1.0) > policy::COMPOSITION_SUM_TOLERANCE.value || std::abs(charge) > policy::COMPOSITION_SUM_TOLERANCE.value
                || inventory_error > policy::ADDITION_MATERIAL_BALANCE.value * std::max(1.0, c.feed_total)
                || maximum_abs(moiety_error) > policy::ADDITION_MATERIAL_BALANCE.value * std::max(1.0, c.feed_total))
                throw std::invalid_argument("finite split arithmetic conservation failed");
            ++report.work["finite_split_evaluations"];
            const auto changed = phase_energy(*c.model, at.temperature, at.pressure, at.volume[p], x,
                report.standard_over_rt, report.work);
            const double delta = remaining * changed.at("g_over_rt") + beta * q.at("g_over_rt")
                - at.amount[p] * donor.at("g_over_rt");
            const double error_estimate = remaining * energy_estimate(changed) + beta * energy_estimate(q)
                + at.amount[p] * energy_estimate(donor) + candidate.feasibility_estimate
                + gamma * (std::abs(remaining * changed.at("g_over_rt")) + std::abs(beta * q.at("g_over_rt"))
                    + std::abs(at.amount[p] * donor.at("g_over_rt"))) / (1.0 - gamma) + beta * standard_error;
            split["delta_g_over_rt"] = delta; split["error_estimate_over_rt"] = error_estimate;
            split["admitted"] = 1.0;
            if (descent_resolved(q.at("tpd_over_rt"), delta, error_estimate)) {
                trial.resolved_descent = true;
                trial.evidence = "numerically_resolved_descending_split_with_estimated_uncertainty";
                q["split_beta_mol"] = beta; q["split_delta_g_over_rt"] = delta;
                q["split_error_estimate_over_rt"] = error_estimate;
                for (std::size_t phase = 0; phase < at.amount.size(); ++phase) {
                    const auto& composition = phase == p ? x : at.composition[phase];
                    for (auto i : c.problem.phases[phase].support) trial.split_coordinates.push_back(std::log(composition[i]));
                    trial.split_coordinates.push_back(phase == p ? remaining : at.amount[phase]);
                    trial.split_coordinates.push_back(std::log(c.model->packing_fraction_coefficient(at.temperature,
                        composition) / at.volume[phase]));
                }
                for (auto i : trial.phase.support) trial.split_coordinates.push_back(std::log(w[i]));
                trial.split_coordinates.push_back(beta);
                trial.split_coordinates.push_back(std::log(c.model->packing_fraction_coefficient(at.temperature, w)
                    / q.at("volume_m3_per_mol")));
                trial.split_checks.push_back(std::move(split));
                return;
            }
        } catch (const std::exception& error) { trial.message += "; finite split: " + std::string(error.what()); }
        trial.split_checks.push_back(std::move(split));
    }
}

StabilityReport test_snapshot(const CompiledSystem& c, const EquilibriumResult& result,
    const StabilityOptions& options, const std::vector<std::pair<PhaseSpec, std::vector<double>>>& removed,
    const core::Trace& trace) {
    StabilityReport report;
    const auto begin = trace.spans.size();
    const auto supports = stability_supports(c, options);
    const auto& z = c.model->parameters().charge_number;
    if (!options.standard_uncertainty_over_rt.empty()
        && (options.standard_uncertainty_over_rt.size() != z.size()
            || !std::all_of(options.standard_uncertainty_over_rt.begin(), options.standard_uncertainty_over_rt.end(),
                [](double x) { return std::isfinite(x) && x >= 0.0; })))
        throw std::invalid_argument("invalid standard-term uncertainty estimates");
    const auto candidate = stability_candidate(c, result, report);
    const auto& at = candidate.physical;
    std::vector<double> shift = report.standard_over_rt;
    for (std::size_t i = 0; i < shift.size(); ++i) shift[i] -= report.plane_over_rt[i];
    for (const auto& phase : supports) {
        std::vector<double> base(z.size(), 0.0), remainder;
        for (auto i : phase.support) base[i] = 1.0;
        base = reconstruct_trial_composition(base, phase.support, z);
        remainder = base;
        for (std::size_t p = 0; p < at.amount.size(); ++p) {
            if (std::all_of(phase.support.begin(), phase.support.end(), [&](auto i) { return at.composition[p][i] > 0.0; })) {
                std::vector<double> restriction(z.size(), 0.0);
                for (auto i : phase.support) restriction[i] = at.composition[p][i];
                remainder = reconstruct_trial_composition(std::move(restriction), phase.support, z); break;
            }
        }
        std::vector<std::pair<std::string, std::vector<double>>> seeds{{"positive-neutral-base", base}};
        if (phase.support.size() > 1) for (auto j : phase.support) if (z[j] == 0.0) {
            auto seed = remainder;
            double sum = 0.0; for (auto i : phase.support) if (i != j) sum += seed[i];
            for (auto i : phase.support) seed[i] = i == j ? 0.99 : 0.01 * seed[i] / sum;
            seeds.emplace_back("near-pure-" + std::to_string(j), std::move(seed));
        }
        if (phase.kind == PhaseKind::Vapor) {
            auto seed = base;
            double maximum = -std::numeric_limits<double>::infinity();
            for (auto i : phase.support) maximum = std::max(maximum, -shift[i]);
            for (auto i : phase.support) seed[i] = std::exp(-shift[i] - maximum);
            seeds.emplace_back("ideal-gas", std::move(seed));
        }
        for (const auto& [previous, composition] : removed)
            if (previous.kind == phase.kind && previous.support == phase.support)
                seeds.emplace_back("removed-phase", composition);
        std::vector<std::pair<std::vector<double>, double>> rooted;
        for (auto& [seed_name, seed] : seeds) {
            StabilityTrial trial; trial.phase = phase; trial.seed = seed_name; trial.evidence = "unresolved";
            try {
                seed = reconstruct_trial_composition(std::move(seed), phase.support, z);
                ++report.work["reconstruction_evaluations"];
                const auto root = stability_root(*c.model, at.temperature, at.pressure, seed, phase.kind, report.work);
                const double eta = c.model->packing_fraction_coefficient(at.temperature, seed) * root.molar_density;
                packing_admission(eta);
                if (std::any_of(rooted.begin(), rooted.end(), [&](const auto& prior) {
                        return prior.first == seed && std::abs(std::log(prior.second / eta))
                            <= policy::DENSITY_ANCHOR_LOG_WINDOW.value; })) continue;
                rooted.emplace_back(seed, eta);
                if (options.trial_budget && report.work["trial_invocations"] >= *options.trial_budget) {
                    report.message = "search-incomplete: trial budget exhausted";
                    count_stability_work(report.work, trace, begin);
                    return report;
                }
                auto compiled = compile_stability_trial(*c.model, at.temperature, at.pressure, phase, seed, shift,
                    trial.row_transformation);
                trial.basis = compiled.B;
                core::Scope scope("tangent_plane_trial");
                const auto initial = trial_coordinates(compiled, seed, eta);
                ++report.work["trial_invocations"];
                EquilibriumResult audit;
                try {
                    trial.optimizer = run_driver(*c.model, compiled, initial, structural_result(compiled, scope),
                        scope, &trial.raw_coordinates, &audit);
                } catch (const std::exception& error) { trial.optimizer = audit; trial.optimizer.message = error.what(); }
                count_driver(report.work, trial.optimizer, "trial_driver_invocations");
                trial.message = trial.optimizer.message;
                if (trial.raw_coordinates.empty()) throw std::invalid_argument("no admitted nonfatal optimizer point");
                const auto raw = compiled.prepare(trial.raw_coordinates, true);
                report.work["energy_evaluations"] += 1.0;
                trial.raw_residuals = compiled.evaluate_at(raw).values;
                trial.composition = reconstruct_trial_composition(raw.composition[0], phase.support, z);
                ++report.work["reconstruction_evaluations"];
                const double volume = c.model->packing_fraction_coefficient(at.temperature, trial.composition) / raw.eta[0];
                trial.quantities = phase_energy(*c.model, at.temperature, at.pressure, volume, trial.composition,
                    report.standard_over_rt, report.work);
                auto& q = trial.quantities;
                double charge = 0.0, sum = -1.0, plane = 0.0, displacement = 0.0;
                for (auto i : phase.support) {
                    sum += trial.composition[i]; charge += z[i] * trial.composition[i];
                    plane += trial.composition[i] * report.plane_over_rt[i];
                    displacement += std::abs(trial.composition[i] - raw.composition[0][i]);
                }
                q["normalization_residual"] = sum; q["charge_residual"] = charge;
                q["reconstruction_l1"] = displacement;
                q["tpd_over_rt"] = q.at("g_over_rt") - plane;
                double dual = 0.0;
                for (std::size_t k = 0; k < compiled.rank; ++k)
                    dual += trial.raw_coordinates[compiled.primal_size + k] * compiled.balance[k];
                q["trial_dual_over_rt"] = dual; q["trial_identity_defect_over_rt"] = q["tpd_over_rt"] - dual;
                q["raw_objective_over_rt"] = trial.optimizer.objective;
                trial.evidence = "no_resolved_descent";
                if (q["tpd_over_rt"] < 0.0) donor_descent(c, candidate, trial, report, options);
            } catch (const std::exception& error) { trial.message += "; " + std::string(error.what()); }
            report.trials.push_back(std::move(trial));
        }
    }
    report.completed = true;
    const auto unresolved = std::count_if(report.trials.begin(), report.trials.end(), [](const auto& t) {
        return t.evidence == "unresolved" || t.evidence == "indeterminate" || t.evidence == "negative_plane_without_positive_donor"; });
    const auto descending = std::count_if(report.trials.begin(), report.trials.end(), [](const auto& t) { return t.resolved_descent; });
    report.message = descending ? "numerically resolved descending split with estimated uncertainty"
        : "no witness in the recorded trial set";
    report.message += "; unresolved trials " + std::to_string(unresolved);
    if (!options.trials && std::all_of(supports.front().support.begin(), supports.front().support.end(),
        [&](auto j) { return z[j] != 0; })) report.message += "; no neutral species";
    count_stability_work(report.work, trace, begin);
    return report;
}
} // namespace

StabilityReport stability_test(const CompiledSystem& c, const EquilibriumResult& result, const StabilityOptions& options) {
    core::Scope scope("stability_test");
    core::Capture capture;
    auto report = test_snapshot(c, result, options, {}, capture.trace());
    count_stability_work(report.work, capture.trace());
    return report;
}

// The event belongs solely to check_topology. Merge finite inventories and extensive volumes,
// then reclose the recipient on its declared density branch before the reduced solve.
std::pair<EquilibriumProblem, std::vector<double>> remove_discovery_phase(
    const CompiledSystem& c, const std::vector<double>& coordinates, Quantities* retained_work) {
    Quantities local_work;
    auto& work = retained_work ? *retained_work : local_work;
    const auto event = check_topology(c, coordinates);
    if (event.event == TopologyEvent::None) throw std::invalid_argument("no topology removal event");
    const auto at = c.prepare(coordinates, true);
    work["energy_evaluations"] += at.amount.size();
    std::size_t dropped = event.phase.value_or(event.pair ? event.pair->second : 0);
    std::optional<std::size_t> recipient;
    if (event.pair) recipient = event.pair->first;
    const auto contains = [&](std::size_t p) {
        const auto& support = c.problem.phases[p].support;
        return std::all_of(c.problem.phases[dropped].support.begin(), c.problem.phases[dropped].support.end(),
            [&](auto i) { return std::find(support.begin(), support.end(), i) != support.end(); });
    };
    if (!recipient) for (std::size_t p = 0; p < at.amount.size(); ++p)
        if (p != dropped && contains(p)) { recipient = p; break; }
    if (!recipient || !contains(*recipient)) throw std::invalid_argument("removal has no inventory recipient");
    if (!(at.amount[dropped] >= 0.0 && at.amount[*recipient] > 0.0))
        throw std::invalid_argument("removal attempt has invalid amounts");
    std::vector<std::vector<double>> n(at.amount.size(), std::vector<double>(c.n_feed.size(), 0.0));
    std::vector<double> before(c.n_feed.size(), 0.0), after(before);
    for (std::size_t p = 0; p < n.size(); ++p) for (std::size_t i = 0; i < before.size(); ++i) {
        n[p][i] = at.amount[p] * at.composition[p][i]; before[i] += n[p][i];
    }
    for (std::size_t i = 0; i < before.size(); ++i) n[*recipient][i] += n[dropped][i];
    const double merged_volume = at.amount[*recipient] * at.volume[*recipient] + at.amount[dropped] * at.volume[dropped];
    EquilibriumProblem reduced = c.problem;
    reduced.phases.erase(reduced.phases.begin() + static_cast<std::ptrdiff_t>(dropped));
    std::vector<double> seed;
    const auto& z = c.model->parameters().charge_number;
    for (std::size_t p = 0; p < n.size(); ++p) if (p != dropped) {
        const double amount = std::accumulate(n[p].begin(), n[p].end(), 0.0);
        std::vector<double> w(n[p].size(), 0.0);
        double charge = 0.0, sum = 0.0;
        for (auto i : c.problem.phases[p].support) {
            w[i] = n[p][i] / amount;
            if (!(w[i] > 0.0 && std::isfinite(w[i]))) throw std::invalid_argument("removal destroys support");
            sum += w[i]; charge += z[i] * w[i]; after[i] += amount * w[i]; seed.push_back(std::log(w[i]));
        }
        if (std::abs(sum - 1.0) > policy::COMPOSITION_SUM_TOLERANCE.value || std::abs(charge) > policy::COMPOSITION_SUM_TOLERANCE.value)
            throw std::invalid_argument("removal arithmetic feasibility failed");
        seed.push_back(amount);
        double volume = p == *recipient ? merged_volume / amount : at.volume[p];
        packing_admission(c.model->packing_fraction_coefficient(at.temperature, w) / volume);
        if (p == *recipient) {
            const auto root = stability_root(*c.model, at.temperature, at.pressure, w, c.problem.phases[p].kind, work);
            volume = 1.0 / root.molar_density;
        }
        seed.push_back(std::log(c.model->packing_fraction_coefficient(at.temperature, w) / volume));
    }
    std::vector<double> error(before.size());
    for (std::size_t i = 0; i < before.size(); ++i) error[i] = after[i] - before[i];
    work["removal_component_error_mol"] = maximum_abs(error);
    work["removal_moiety_error_mol"] = maximum_abs(matvec(c.B, error));
    if (maximum_abs(error) > policy::ADDITION_MATERIAL_BALANCE.value * std::max(1.0, c.feed_total)
        || maximum_abs(matvec(c.B, error)) > policy::ADDITION_MATERIAL_BALANCE.value * std::max(1.0, c.feed_total))
        throw std::invalid_argument("removal loses inventory");
    return {std::move(reduced), std::move(seed)};
}

EquilibriumResult discover_phases(const HelmholtzModel& model, const EquilibriumProblem& request) {
    if (!request.discovery || request.topology_declared || !request.phases.empty()
        || !request.temperature.is_positive() || !request.pressure.is_positive()
        || request.temperature.free || request.pressure.free)
        throw std::invalid_argument("discovery requires explicit unknown topology at fixed positive T/P");
    core::Capture capture;
    const auto& options = *request.discovery;
    const auto feed = resolve_feed(model, request.feed);
    const auto& z = model.parameters().charge_number;
    std::vector<std::size_t> united;
    bool has_charge = false;
    for (std::size_t i = 0; i < feed.size(); ++i) {
        bool participant = feed[i] > 0.0;
        for (const auto& reaction : request.reactions) {
            if (reaction.stoichiometry.size() != feed.size()) throw std::invalid_argument("stoichiometry size");
            participant |= reaction.stoichiometry[i] != 0.0;
        }
        if (participant) { united.push_back(i); has_charge |= z[i] != 0.0; }
    }
    auto record = std::make_shared<DiscoveryRecord>();
    EquilibriumProblem declaration = request;
    declaration.discovery.reset(); declaration.topology_declared = true;
    PhaseSpec initial; initial.name = "L1"; initial.support = united; initial.support_declared = true;
    if (!request.reactions.empty()) {
        const auto base = reconstruct_trial_composition(std::vector<double>(feed.size(), 1.0), united, z);
        for (auto i : united) initial.composition_guess.push_back(base[i]);
    }
    declaration.phases = {initial};
    auto compiled = compile_problem(model, declaration);
    // Admit the requested trial physics before any declaration solve.
    stability_supports(compiled, options);
    const std::size_t rank = compiled.rank;
    const std::size_t budget = options.declaration_budget.value_or(2 * rank);
    std::size_t calls = 0;
    std::vector<std::pair<PhaseSpec, std::vector<double>>> removed;
    auto invoke = [&](const CompiledSystem& c, std::vector<double> seed, const std::string& role,
                      std::vector<double>& attempt) {
        EquilibriumResult result;
        if (calls >= budget) { result.message = "declaration budget exhausted"; return result; }
        ++calls;
        core::Scope scope("phase_discovery_declaration");
        EquilibriumResult audit = structural_result(c, scope);
        bool invoked = false;
        try {
            if (seed.empty()) seed = boundary_start(c);
            else seed = add_multipliers(c, seed);
            invoked = true;
            result = run_driver(model, c, std::move(seed), audit, scope, &attempt, &audit);
            if (result.success) {
                const auto branch = check_density_branches(model, c, result, {}, true, scope);
                if (!branch.accepted) { result.success = false; result.message = branch.message; clear_state(result); }
            }
        } catch (const std::exception& error) { result = audit; result.message = error.what(); }
        record->declarations.push_back(result); record->declaration_roles.push_back(role);
        record->declaration_attempt_coordinates.push_back(attempt);
        if (invoked) count_driver(record->work, result, "declaration_invocations");
        else ++record->work["declaration_initialization_failures"];
        return result;
    };
    EquilibriumResult best;
    auto best_compiled = compiled;
    for (const auto kind : {PhaseKind::Liquid, PhaseKind::Vapor}) {
        if (has_charge && kind == PhaseKind::Vapor) continue;
        declaration.phases[0].kind = kind; declaration.phases[0].name = kind == PhaseKind::Liquid ? "L1" : "V1";
        compiled = compile_problem(model, declaration);
        std::vector<double> attempt;
        auto result = invoke(compiled, {}, "initial", attempt);
        if (result.success && (!best.success || result.objective < best.objective)) {
            best = std::move(result); best_compiled = compiled;
        }
    }
    const auto finish = [&](bool success, const std::string& message) {
        record->declaration = best_compiled.problem; record->candidate = best;
        // The accepted state also supplies the fixed declaration's ordinary cold
        // initializer, including same-kind and three-phase declarations.
        if (best.has_state()) for (std::size_t p = 0; p < best.phase_count; ++p) {
            auto& phase = record->declaration.phases[p];
            phase.composition_guess.clear();
            for (auto j : phase.support)
                phase.composition_guess.push_back(best.mole_fractions[p * feed.size() + j]);
            phase.packing_guess = best.packing_fractions[p];
        }
        count_stability_work(record->work, capture.trace());
        EquilibriumResult result = best; result.discovery = record;
        if (!success) { result.success = false; result.message = message; clear_state(result); }
        return result;
    };
    if (!best.success) return finish(false, "discovery has no accepted initial candidate");
    for (std::size_t round = 0; round < 2 * rank + 1; ++round) {
        auto trial_options = options;
        if (options.trial_budget) {
            const auto used = static_cast<std::size_t>(record->work["trial_invocations"]);
            trial_options.trial_budget = used < *options.trial_budget ? *options.trial_budget - used : 0;
        }
        auto report = test_snapshot(best_compiled, best, trial_options, removed, capture.trace());
        for (const auto& [name, value] : report.work) record->work[name] += value;
        record->rounds.push_back(std::move(report));
        auto& current = record->rounds.back();
        if (!current.completed) return finish(false, current.message);
        std::vector<std::size_t> descending;
        for (std::size_t i = 0; i < current.trials.size(); ++i) if (current.trials[i].resolved_descent) {
            const auto& trial = current.trials[i];
            const bool duplicate = std::any_of(descending.begin(), descending.end(), [&](auto j) {
                const auto& prior = current.trials[j];
                return prior.phase.support == trial.phase.support && prior.phase.kind == trial.phase.kind
                    && largest_difference(prior.composition, trial.composition) <= policy::PHASE_COINCIDENCE_DISTANCE.value
                    && std::abs(std::log(prior.quantities.at("volume_m3_per_mol") / trial.quantities.at("volume_m3_per_mol")))
                        <= policy::PHASE_COINCIDENCE_DISTANCE.value;
            });
            if (!duplicate) descending.push_back(i);
        }
        if (descending.empty()) return finish(true, current.message);
        if (best.phase_count >= rank) return finish(false, "resolved descent at balance-rank phase limit; replacement deferred");
        std::stable_sort(descending.begin(), descending.end(), [&](auto i, auto j) {
            return current.trials[i].quantities.at("donor_slope_over_rt") < current.trials[j].quantities.at("donor_slope_over_rt"); });
        bool added = false;
        for (auto i : descending) {
            if (calls >= budget) break;
            auto next = best_compiled.problem;
            auto phase = current.trials[i].phase;
            auto& trial = current.trials[i];
            // The optimizer may cross its starting root's label. Reuse the density
            // owner's two admitted end branches; a unique root keeps its source label.
            try {
                DensityRoot sibling;
                auto root = solve_density(model, best.temperature[0], best.pressure[0], trial.composition,
                    phase.kind == PhaseKind::Liquid ? DensityBranch::Liquid : DensityBranch::Vapor,
                    std::nullopt, &sibling);
                const double volume = trial.quantities.at("volume_m3_per_mol");
                if (std::abs(std::log(root.molar_density * volume)) > policy::DENSITY_ANCHOR_LOG_WINDOW.value)
                    root = sibling;
                trial.quantities["addition_root_log_displacement"] = std::abs(std::log(root.molar_density * volume));
                trial.quantities["addition_root_residual"] = root.convergence.residual;
                if (!std::isfinite(root.convergence.residual) || !std::isfinite(trial.quantities.at("addition_root_log_displacement"))
                    || root.convergence.residual > policy::DENSITY_ROOT_RESIDUAL.value
                    || trial.quantities.at("addition_root_log_displacement") > policy::DENSITY_ANCHOR_LOG_WINDOW.value)
                    throw std::runtime_error("descending trial has no admitted density branch");
                if (root.branch != DensityBranch::Unspecified)
                    phase.kind = root.branch == DensityBranch::Liquid ? PhaseKind::Liquid : PhaseKind::Vapor;
                if (phase.kind == PhaseKind::Vapor && std::any_of(phase.support.begin(), phase.support.end(),
                    [&](auto j) { return z[j] != 0; })) throw std::runtime_error("ionic descending trial left liquid branch");
            } catch (const std::exception& error) { trial.message += "; addition: " + std::string(error.what()); continue; }
            phase.name = (phase.kind == PhaseKind::Liquid ? "L" : "V") + std::to_string(next.phases.size() + 1);
            next.phases.push_back(phase);
            auto system = compile_problem(model, next);
            std::vector<double> attempt;
            auto result = invoke(system, current.trials[i].split_coordinates, "addition", attempt);
            while (!result.success && result.topology_event != TopologyEvent::None && !attempt.empty() && calls < budget) {
                const auto before = system.prepare(attempt, true);
                const auto event = check_topology(system, attempt);
                const auto dropped = event.phase.value_or(event.pair ? event.pair->second : 0);
                removed.emplace_back(system.problem.phases[dropped], before.composition[dropped]);
                try {
                    auto [reduced, seed] = remove_discovery_phase(system, attempt, &record->work);
                    system = compile_problem(model, reduced); attempt.clear();
                    result = invoke(system, std::move(seed), "removal", attempt);
                } catch (const std::exception& error) { result.message += "; " + std::string(error.what()); break; }
            }
            if (!result.success || result.convergence.residual > request.row_tolerance) continue;
            double error_estimate = 0.0;
            for (std::size_t p = 0; p < current.present_phases.size(); ++p)
                error_estimate += best.amounts[p] * energy_estimate(current.present_phases[p]);
            trial.quantities["addition_gibbs_estimate_over_rt_mol"] = error_estimate;
            if (!(result.objective < best.objective - error_estimate)) continue;
            best = std::move(result); best_compiled = std::move(system); added = true; break;
        }
        if (!added) return finish(false, calls >= budget ? "resolved descent with exhausted declaration budget"
            : "resolved descent without an admitted lower-G declaration");
    }
    return finish(false, "search-incomplete: trial-round budget exhausted");
}

}  // namespace epcsaft::equilibrium
