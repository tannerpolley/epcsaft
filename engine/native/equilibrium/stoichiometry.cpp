#include "equilibrium/stoichiometry.hpp"

#include "core/physical_constants.hpp"
#include "properties/thermochemistry.hpp"

#include <Eigen/Dense>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>

namespace epcsaft::equilibrium {
namespace {

void require_finite(const std::vector<double>& values, const char* name) {
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
    if (!(std::abs(vector_norm(values) - 1.0) <= policy::COMPOSITION_SUM_TOLERANCE.value)) {
        throw std::invalid_argument(std::string(name) + " must sum to one");
    }
}

Eigen::MatrixXd reaction_matrix(
    std::size_t components, const std::vector<ReactionSpec>& reactions
) {
    Eigen::MatrixXd matrix(
        static_cast<Eigen::Index>(reactions.size()), static_cast<Eigen::Index>(components)
    );
    for (std::size_t row = 0; row < reactions.size(); ++row) {
        if (reactions[row].stoichiometry.size() != components) {
            throw std::invalid_argument("stoichiometry size does not match components");
        }
        for (std::size_t column = 0; column < components; ++column) {
            const double value = reactions[row].stoichiometry[column];
            if (!std::isfinite(value)) throw std::invalid_argument("invalid stoichiometry");
            matrix(static_cast<Eigen::Index>(row), static_cast<Eigen::Index>(column)) = value;
        }
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

void validate_basis_shape(
    const std::vector<std::vector<double>>& basis, std::size_t components
) {
    for (const auto& row : basis) {
        if (row.size() != components) throw std::logic_error("balance basis width");
        if (!std::all_of(row.begin(), row.end(),
                         [](double value) { return std::isfinite(value); })) {
            throw std::logic_error("balance basis is not finite");
        }
    }
}

}  // namespace

std::vector<double> resolve_feed(const HelmholtzModel& model, const FeedSpec& feed) {
    const std::size_t components = model.parameters().component_count();
    if (feed.values.size() != components) {
        throw std::invalid_argument("feed size does not match components");
    }
    require_finite(feed.values, "feed");
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
        validate_basis_shape(result, components);
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
    validate_basis_shape(result, components);
    return result;
}

std::vector<std::vector<double>> reaction_target_series(
    double temperature, const std::vector<ReactionSpec>& reactions, std::size_t order
) {
    require_positive(temperature, "temperature");
    if (order > 3) throw std::invalid_argument("reaction target order");
    const double R = GAS_CONSTANT_J_PER_MOL_K;
    std::vector<std::vector<double>> result(order + 1, std::vector<double>(reactions.size(), 0.0));
    for (std::size_t row = 0; row < reactions.size(); ++row) {
        const auto& reaction = reactions[row];
        // l[k] = d^k ln K / dT^k; t = -R T ln K, so t^(k) = -R (k l[k-1] + T l[k]).
        std::array<double, 4> l{};
        if (reaction.affinity_j_per_mol.has_value()) {
            result[0][row] = *reaction.affinity_j_per_mol;
            continue;
        }
        if (reaction.correlation.has_value()) {
            if (reaction.correlation->standard_state_id != EOS_STANDARD_STATE_ID) {
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

std::vector<std::vector<double>> standard_state_series(
    double temperature, std::size_t components,
    const std::vector<ReactionSpec>& reactions,
    const std::vector<double>& supplied, std::size_t order
) {
    require_positive(temperature, "temperature");
    validate_reaction_set(components, reactions);
    std::vector<std::vector<double>> result(order + 1, std::vector<double>(components, 0.0));
    if (reactions.empty() && supplied.empty()) return result;
    const auto targets = reaction_target_series(temperature, reactions, order);
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
    const std::vector<double>& supplied
) {
    return standard_state_series(temperature, components, reactions, supplied, 0).front();
}

bool conserves_total_moles(const std::vector<std::vector<double>>& basis) {
    if (basis.empty()) throw std::invalid_argument("empty balance basis");
    const std::size_t components = basis.front().size();
    validate_basis_shape(basis, components);
    Eigen::MatrixXd extended(static_cast<Eigen::Index>(basis.size() + 1),
                             static_cast<Eigen::Index>(components));
    for (std::size_t row = 0; row < basis.size(); ++row) {
        for (std::size_t column = 0; column < components; ++column) {
            extended(static_cast<Eigen::Index>(row), static_cast<Eigen::Index>(column)) =
                basis[row][column];
        }
    }
    extended.row(static_cast<Eigen::Index>(basis.size())).setConstant(
        1.0 / std::sqrt(static_cast<double>(components)));
    return numerical_rank(extended) == numerical_rank(extended.topRows(
        static_cast<Eigen::Index>(basis.size())));
}

std::vector<double> matvec(
    const std::vector<std::vector<double>>& matrix, const std::vector<double>& vector
) {
    std::vector<double> result(matrix.size(), 0.0);
    for (std::size_t row = 0; row < matrix.size(); ++row) {
        if (matrix[row].size() != vector.size()) throw std::invalid_argument("matrix width");
        for (std::size_t column = 0; column < vector.size(); ++column) {
            result[row] += matrix[row][column] * vector[column];
        }
    }
    return result;
}

double vector_norm(const std::vector<double>& values) {
    double total = 0.0;
    for (const double value : values) total += value;
    return total;
}

double feed_total(const std::vector<double>& feed) {
    if (feed.empty()) throw std::invalid_argument("empty feed");
    if (!std::all_of(feed.begin(), feed.end(), [](double value) {
            return std::isfinite(value) && value >= 0.0;
        })) {
        throw std::invalid_argument("invalid feed");
    }
    const double total = vector_norm(feed);
    if (!(std::isfinite(total) && total > 0.0)) {
        throw std::invalid_argument("feed has zero total");
    }
    return total;
}

std::vector<double> normalized_feed(const std::vector<double>& feed) {
    const double total = feed_total(feed);
    std::vector<double> result = feed;
    for (double& value : result) value /= total;
    return result;
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
        if (targets != 1) {
            throw std::invalid_argument("inconsistent standard state: reaction " + label
                                        + " needs exactly one target");
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

std::vector<double> reaction_affinities(
    const std::vector<ReactionSpec>& reactions,
    const std::vector<double>& standard_gibbs
) {
    if (standard_gibbs.empty()) return std::vector<double>(reactions.size(), 0.0);
    std::vector<double> result;
    result.reserve(reactions.size());
    for (const auto& reaction : reactions) {
        if (reaction.stoichiometry.size() != standard_gibbs.size()) {
            throw std::invalid_argument("stoichiometry size does not match standard state");
        }
        result.push_back(std::inner_product(
            reaction.stoichiometry.begin(), reaction.stoichiometry.end(),
            standard_gibbs.begin(), 0.0));
    }
    return result;
}

bool all_numeric(const std::vector<double>& values) {
    return std::all_of(values.begin(), values.end(),
                       [](double value) { return std::isfinite(value); });
}

}  // namespace epcsaft::equilibrium

