#include "equilibrium/stoichiometry.hpp"

#include "core/physical_constants.hpp"

#include <Eigen/Dense>

#include <algorithm>
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
        const double total = vector_norm(feed.values);
        if (!(total > 0.0)) throw std::invalid_argument("feed has zero total");
        for (std::size_t i = 0; i < components; ++i) result[i] = feed.values[i] / total;
    } else if (feed.kind == FeedKind::MassFractions) {
        const auto& masses = model.parameters().molar_mass_kg_per_mol;
        require_positive(vector_norm(feed.values), "mass fractions");
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
    const auto svd = matrix.jacobiSvd(Eigen::ComputeFullV);
    const Eigen::Index rows = static_cast<Eigen::Index>(components) - rank;
    std::vector<std::vector<double>> result(
        static_cast<std::size_t>(rows), std::vector<double>(components, 0.0)
    );
    for (Eigen::Index row = 0; row < rows; ++row) {
        for (std::size_t column = 0; column < components; ++column) {
            result[static_cast<std::size_t>(row)][column] =
                svd.matrixV()(static_cast<Eigen::Index>(column), rank + row);
        }
    }
    validate_basis_shape(result, components);
    return result;
}

std::vector<double> standard_state(
    double temperature, std::size_t components,
    const std::vector<ReactionSpec>& reactions,
    const std::vector<double>& supplied
) {
    require_positive(temperature, "temperature");
    validate_reaction_set(components, reactions);
    if (!supplied.empty()) {
        if (supplied.size() != components) throw std::invalid_argument("standard state size");
        if (!std::all_of(supplied.begin(), supplied.end(), [](double value) {
                return std::isfinite(value);
            })) throw std::invalid_argument("invalid standard state");
        return supplied;
    }
    if (reactions.empty()) return std::vector<double>(components, 0.0);
    const auto matrix = reaction_matrix(components, reactions);
    Eigen::VectorXd target(static_cast<Eigen::Index>(reactions.size()));
    for (std::size_t row = 0; row < reactions.size(); ++row) {
        if (reactions[row].affinity_j_per_mol.has_value()) {
            target(static_cast<Eigen::Index>(row)) = *reactions[row].affinity_j_per_mol;
        } else {
            if (!(reactions[row].equilibrium_constant > 0.0)) {
                throw std::invalid_argument("equilibrium constant is not positive");
            }
            target(static_cast<Eigen::Index>(row)) = -GAS_CONSTANT_J_PER_MOL_K * temperature
                * std::log(reactions[row].equilibrium_constant);
        }
    }
    const Eigen::VectorXd solution = matrix.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV)
        .solve(target);
    const double error = (matrix * solution - target).norm();
    if (error > std::sqrt(std::numeric_limits<double>::epsilon())
        * std::max(1.0, target.norm())) {
        throw std::invalid_argument("inconsistent standard state");
    }
    return std::vector<double>(solution.data(), solution.data() + solution.size());
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
        if (reaction.stoichiometry.size() != components) {
            throw std::invalid_argument("stoichiometry size does not match components");
        }
        if (!std::all_of(reaction.stoichiometry.begin(), reaction.stoichiometry.end(),
                         [](double value) { return std::isfinite(value); })) {
            throw std::invalid_argument("invalid stoichiometry in reaction "
                                        + reaction_name(reaction, row));
        }
        if (reaction.affinity_j_per_mol.has_value()
            && !std::isfinite(*reaction.affinity_j_per_mol)) {
            throw std::invalid_argument("invalid affinity in reaction "
                                        + reaction_name(reaction, row));
        }
        if (!reaction.affinity_j_per_mol.has_value()
            && (!std::isfinite(reaction.equilibrium_constant)
                || reaction.equilibrium_constant <= 0.0)) {
            throw std::invalid_argument("equilibrium constant is not positive in reaction "
                                        + reaction_name(reaction, row));
        }
    }
}

std::size_t reaction_rank(
    std::size_t components, const std::vector<ReactionSpec>& reactions
) {
    validate_reaction_set(components, reactions);
    return static_cast<std::size_t>(numerical_rank(reaction_matrix(components, reactions)));
}

std::vector<std::vector<double>> reaction_projector(
    std::size_t components, const std::vector<ReactionSpec>& reactions
) {
    const auto basis = balance_basis(components, reactions);
    std::vector<std::vector<double>> result(
        components, std::vector<double>(components, 0.0));
    for (const auto& row : basis) {
        for (std::size_t i = 0; i < components; ++i) {
            for (std::size_t j = 0; j < components; ++j) {
                result[i][j] += row[i] * row[j];
            }
        }
    }
    return result;
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

StoichiometryReport inspect_stoichiometry(
    std::size_t components, const std::vector<ReactionSpec>& reactions
) {
    StoichiometryReport report;
    report.components = components;
    report.reactions = reactions.size();
    if (components == 0) {
        report.finite = false;
        report.independent = false;
        return report;
    }
    for (const auto& reaction : reactions) {
        report.finite = report.finite && reaction.stoichiometry.size() == components;
        report.finite = report.finite && std::all_of(
            reaction.stoichiometry.begin(), reaction.stoichiometry.end(),
            [](double value) { return std::isfinite(value); });
    }
    if (!report.finite) {
        report.independent = false;
        return report;
    }
    report.rank = reaction_rank(components, reactions);
    report.independent = report.rank == reactions.size();
    report.finite = report.finite && report.rank <= components;
    // A dependent reaction set is reported, not silently repaired.
    return report;
}

std::vector<double> reaction_residuals(
    const std::vector<ReactionSpec>& reactions,
    const std::vector<double>& standard_gibbs,
    const std::vector<double>& requested_affinities
) {
    const auto affinities = reaction_affinities(reactions, standard_gibbs);
    if (!requested_affinities.empty() && requested_affinities.size() != affinities.size()) {
        throw std::invalid_argument("requested affinity size");
    }
    std::vector<double> result = affinities;
    if (!requested_affinities.empty()) {
        for (std::size_t i = 0; i < result.size(); ++i) {
            result[i] -= requested_affinities[i];
        }
    }
    return result;
}

std::vector<double> balance_residual(
    const std::vector<std::vector<double>>& basis,
    const std::vector<double>& amounts,
    const std::vector<double>& target
) {
    if (basis.size() != amounts.size()) {
        throw std::invalid_argument("balance amount size");
    }
    const auto actual = matvec(basis, amounts);
    if (actual.size() != target.size()) throw std::invalid_argument("balance target size");
    std::vector<double> result(actual.size(), 0.0);
    for (std::size_t i = 0; i < actual.size(); ++i) result[i] = actual[i] - target[i];
    return result;
}

std::vector<std::vector<double>> gram_matrix(
    const std::vector<std::vector<double>>& basis
) {
    if (basis.empty()) return {};
    const std::size_t width = basis.front().size();
    for (const auto& row : basis) {
        if (row.size() != width) throw std::invalid_argument("basis width");
    }
    std::vector<std::vector<double>> result(
        basis.size(), std::vector<double>(basis.size(), 0.0));
    for (std::size_t i = 0; i < basis.size(); ++i) {
        for (std::size_t j = 0; j < basis.size(); ++j) {
            for (std::size_t column = 0; column < width; ++column) {
                result[i][j] += basis[i][column] * basis[j][column];
            }
        }
    }
    return result;
}

std::vector<double> transpose_matvec(
    const std::vector<std::vector<double>>& matrix,
    const std::vector<double>& vector
) {
    if (matrix.size() != vector.size()) throw std::invalid_argument("matrix height");
    if (matrix.empty()) return {};
    const std::size_t width = matrix.front().size();
    std::vector<double> result(width, 0.0);
    for (std::size_t row = 0; row < matrix.size(); ++row) {
        if (matrix[row].size() != width) throw std::invalid_argument("matrix width");
        for (std::size_t column = 0; column < width; ++column) {
            result[column] += matrix[row][column] * vector[row];
        }
    }
    return result;
}

bool is_orthonormal(
    const std::vector<std::vector<double>>& basis, double tolerance
) {
    if (!(std::isfinite(tolerance) && tolerance >= 0.0)) return false;
    const auto gram = gram_matrix(basis);
    for (std::size_t i = 0; i < gram.size(); ++i) {
        for (std::size_t j = 0; j < gram.size(); ++j) {
            const double expected = i == j ? 1.0 : 0.0;
            if (std::abs(gram[i][j] - expected) > tolerance) return false;
        }
    }
    return true;
}

std::vector<double> reaction_affinity_targets(
    double temperature, const std::vector<ReactionSpec>& reactions
) {
    require_positive(temperature, "temperature");
    std::vector<double> result;
    result.reserve(reactions.size());
    for (const auto& reaction : reactions) {
        if (reaction.affinity_j_per_mol.has_value()) {
            result.push_back(*reaction.affinity_j_per_mol);
            continue;
        }
        require_positive(reaction.equilibrium_constant, "equilibrium constant");
        result.push_back(-GAS_CONSTANT_J_PER_MOL_K * temperature
                         * std::log(reaction.equilibrium_constant));
    }
    return result;
}

std::vector<double> project_to_balance(
    const std::vector<std::vector<double>>& basis,
    const std::vector<double>& composition
) {
    if (basis.empty()) return {};
    return matvec(basis, composition);
}

std::vector<double> reconstruct_from_balance(
    const std::vector<std::vector<double>>& basis,
    const std::vector<double>& coordinates
) {
    if (basis.size() != coordinates.size()) {
        throw std::invalid_argument("balance coordinate size");
    }
    return transpose_matvec(basis, coordinates);
}

double l1_norm(const std::vector<double>& values) {
    double result = 0.0;
    for (const double value : values) {
        if (!std::isfinite(value)) return std::numeric_limits<double>::infinity();
        result += std::abs(value);
    }
    return result;
}

bool all_finite(const std::vector<double>& values) {
    return std::all_of(values.begin(), values.end(),
                       [](double value) { return std::isfinite(value); });
}

double feed_mass_kg(
    const HelmholtzModel& model, const std::vector<double>& amounts
) {
    const auto& masses = model.parameters().molar_mass_kg_per_mol;
    if (amounts.size() != masses.size()) throw std::invalid_argument("amount size");
    double result = 0.0;
    for (std::size_t i = 0; i < amounts.size(); ++i) {
        if (!std::isfinite(amounts[i]) || amounts[i] < 0.0) {
            throw std::invalid_argument("invalid amount");
        }
        require_positive(masses[i], "molar mass");
        result += amounts[i] * masses[i];
    }
    if (!(std::isfinite(result) && result > 0.0)) {
        throw std::invalid_argument("feed has zero mass");
    }
    return result;
}

std::vector<double> scale_amounts(
    const std::vector<double>& amounts, double total
) {
    require_positive(total, "amount total");
    const double source = feed_total(amounts);
    std::vector<double> result = amounts;
    const double factor = total / source;
    for (double& value : result) value *= factor;
    return result;
}

bool is_composition(const std::vector<double>& values, double tolerance) {
    if (!(std::isfinite(tolerance) && tolerance >= 0.0) || values.empty()) return false;
    double total = 0.0;
    for (const double value : values) {
        if (!std::isfinite(value) || value < 0.0) return false;
        total += value;
    }
    return std::abs(total - 1.0) <= tolerance;
}

double projection_error(
    const std::vector<std::vector<double>>& basis,
    const std::vector<double>& composition
) {
    const auto coordinates = project_to_balance(basis, composition);
    const auto reconstruction = reconstruct_from_balance(basis, coordinates);
    if (reconstruction.size() != composition.size()) {
        throw std::invalid_argument("basis width");
    }
    double result = 0.0;
    for (std::size_t i = 0; i < composition.size(); ++i) {
        result = std::max(result, std::abs(reconstruction[i] - composition[i]));
    }
    return result;
}

double maximum_abs(const std::vector<double>& values) {
    double result = 0.0;
    for (const double value : values) result = std::max(result, std::abs(value));
    return result;
}

double basis_trace(const std::vector<std::vector<double>>& basis) {
    const auto gram = gram_matrix(basis);
    double result = 0.0;
    for (std::size_t i = 0; i < gram.size(); ++i) result += gram[i][i];
    return result;
}

bool same_subspace(
    const std::vector<std::vector<double>>& left,
    const std::vector<std::vector<double>>& right,
    double tolerance
) {
    if (!(std::isfinite(tolerance) && tolerance >= 0.0)) return false;
    if (left.empty() || right.empty()) return left.empty() && right.empty();
    if (left.front().size() != right.front().size()) return false;
    if (left.size() != right.size()) return false;
    for (const auto& row : left) {
        if (projection_error(right, row) > tolerance) return false;
    }
    for (const auto& row : right) {
        if (projection_error(left, row) > tolerance) return false;
    }
    return true;
}

double balance_total(const std::vector<double>& amounts) {
    if (!all_finite(amounts)) throw std::invalid_argument("invalid balance amounts");
    return std::accumulate(amounts.begin(), amounts.end(), 0.0);
}

bool basis_matches_components(
    const std::vector<std::vector<double>>& basis, std::size_t components
) {
    if (components == 0) return false;
    for (const auto& row : basis) {
        if (row.size() != components || !all_finite(row)) return false;
    }
    return true;
}

bool contains_reaction_name(
    const std::vector<ReactionSpec>& reactions, const std::string& name
) {
    for (std::size_t index = 0; index < reactions.size(); ++index) {
        if (reaction_name(reactions[index], index) == name) return true;
    }
    return false;
}

}  // namespace epcsaft::equilibrium

