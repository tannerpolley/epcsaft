#include "chemical_equilibrium.hpp"

#include <Eigen/Dense>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

namespace epcsaft_engine::equilibrium {
namespace {

void require_finite(double value, const char* field) {
    if (!std::isfinite(value)) {
        throw std::invalid_argument(std::string(field) + " must be finite");
    }
}

void require_positive(double value, const char* field) {
    require_finite(value, field);
    if (value <= 0.0) {
        throw std::invalid_argument(std::string(field) + " must be positive");
    }
}

bool is_trimmed_nonempty(const std::string& value) {
    if (value.empty()) {
        return false;
    }
    const auto non_space = [](unsigned char character) {
        return std::isspace(character) == 0;
    };
    return non_space(static_cast<unsigned char>(value.front()))
        && non_space(static_cast<unsigned char>(value.back()));
}

void require_unique_ids(
    const std::vector<std::string>& values,
    const char* field
) {
    if (values.empty()) {
        throw std::invalid_argument(std::string(field) + " must not be empty");
    }
    std::unordered_set<std::string> identities;
    identities.reserve(values.size());
    for (std::size_t index = 0; index < values.size(); ++index) {
        if (!is_trimmed_nonempty(values[index])) {
            throw std::invalid_argument(
                std::string(field) + "[" + std::to_string(index)
                + "] must be a nonempty trimmed identity"
            );
        }
        if (!identities.insert(values[index]).second) {
            throw std::invalid_argument(
                std::string(field) + " must contain unique identities"
            );
        }
    }
}

void require_matrix_shape(
    const DenseMatrix& matrix,
    std::size_t expected_rows,
    std::size_t expected_columns,
    const char* field
) {
    if (matrix.rows != expected_rows || matrix.columns != expected_columns) {
        throw std::invalid_argument(
            std::string(field) + " has shape ("
            + std::to_string(matrix.rows) + ", "
            + std::to_string(matrix.columns) + "), expected ("
            + std::to_string(expected_rows) + ", "
            + std::to_string(expected_columns) + ")"
        );
    }
    if (matrix.values.size() != matrix.rows * matrix.columns) {
        throw std::invalid_argument(
            std::string(field) + " storage is inconsistent with its shape"
        );
    }
    for (std::size_t index = 0; index < matrix.values.size(); ++index) {
        require_finite(matrix.values[index], field);
    }
}

Eigen::MatrixXd to_eigen(const DenseMatrix& source) {
    Eigen::MatrixXd result(
        static_cast<Eigen::Index>(source.rows),
        static_cast<Eigen::Index>(source.columns)
    );
    for (std::size_t row = 0; row < source.rows; ++row) {
        for (std::size_t column = 0; column < source.columns; ++column) {
            result(
                static_cast<Eigen::Index>(row),
                static_cast<Eigen::Index>(column)
            ) = source.values[row * source.columns + column];
        }
    }
    return result;
}

DenseMatrix from_eigen(const Eigen::MatrixXd& source) {
    DenseMatrix result;
    result.rows = static_cast<std::size_t>(source.rows());
    result.columns = static_cast<std::size_t>(source.cols());
    result.values.reserve(result.rows * result.columns);
    for (Eigen::Index row = 0; row < source.rows(); ++row) {
        for (Eigen::Index column = 0; column < source.cols(); ++column) {
            result.values.push_back(source(row, column));
        }
    }
    return result;
}

double inf_norm(const Eigen::MatrixXd& matrix) {
    if (matrix.size() == 0) {
        return 0.0;
    }
    return matrix.cwiseAbs().maxCoeff();
}

double inf_norm(const Eigen::VectorXd& vector) {
    if (vector.size() == 0) {
        return 0.0;
    }
    return vector.cwiseAbs().maxCoeff();
}

// Eigen's SVD fixes the singular-vector order; this final sign convention
// removes the remaining harmless +/- ambiguity.  The row space, rather than
// the particular basis orientation, is the physical object and is therefore
// also checked by the compiler diagnostics.
void canonicalize_row_signs(Eigen::MatrixXd& rows) {
    for (Eigen::Index row = 0; row < rows.rows(); ++row) {
        Eigen::Index pivot = 0;
        double pivot_magnitude = 0.0;
        for (Eigen::Index column = 0; column < rows.cols(); ++column) {
            const double magnitude = std::abs(rows(row, column));
            if (magnitude > pivot_magnitude) {
                pivot = column;
                pivot_magnitude = magnitude;
            }
        }
        if (pivot_magnitude > 0.0 && rows(row, pivot) < 0.0) {
            rows.row(row) *= -1.0;
        }
    }
}

}  // namespace

PhaseDirectionCompilation compile_phase_direction_system(
    const std::vector<std::string>& phase_ids,
    const std::vector<std::string>& species_ids,
    const std::vector<double>& charges,
    const DenseMatrix& balance_matrix,
    const DenseMatrix& phase_reference_amounts,
    const DenseMatrix& direction_matrix,
    const std::vector<double>& driving_constants,
    double rank_tolerance,
    double conservation_tolerance
) {
    require_unique_ids(phase_ids, "phase_ids");
    require_unique_ids(species_ids, "species_ids");
    require_positive(rank_tolerance, "rank_tolerance");
    require_positive(conservation_tolerance, "conservation_tolerance");

    const std::size_t phase_count = phase_ids.size();
    const std::size_t species_count = species_ids.size();
    const std::size_t dimension = phase_count * species_count;
    if (charges.size() != species_count) {
        throw std::invalid_argument("charges must match species_ids");
    }
    for (double charge : charges) {
        require_finite(charge, "charges");
    }
    require_matrix_shape(
        balance_matrix,
        balance_matrix.rows,
        species_count,
        "balance_matrix"
    );
    require_matrix_shape(
        phase_reference_amounts,
        phase_count,
        species_count,
        "phase_reference_amounts"
    );
    require_matrix_shape(
        direction_matrix,
        direction_matrix.rows,
        dimension,
        "direction_matrix"
    );
    const std::size_t direction_count = direction_matrix.rows;
    if (driving_constants.size() != direction_count) {
        throw std::invalid_argument(
            "driving_constants must match direction_matrix rows"
        );
    }
    for (double constant : driving_constants) {
        require_finite(constant, "driving_constants");
    }
    for (double amount : phase_reference_amounts.values) {
        if (amount < 0.0) {
            throw std::invalid_argument(
                "phase_reference_amounts must be nonnegative"
            );
        }
    }

    const Eigen::MatrixXd directions = to_eigen(direction_matrix);
    const Eigen::MatrixXd references = to_eigen(phase_reference_amounts);
    const Eigen::MatrixXd balances = to_eigen(balance_matrix);
    const double direction_scale = directions.size() == 0
        ? 0.0
        : directions.cwiseAbs().maxCoeff();
    const double rank_threshold =
        rank_tolerance * std::max(1.0, direction_scale);

    std::size_t direction_rank = 0;
    Eigen::MatrixXd invariants = Eigen::MatrixXd::Identity(
        static_cast<Eigen::Index>(dimension),
        static_cast<Eigen::Index>(dimension)
    );
    Eigen::VectorXd reference_potentials = Eigen::VectorXd::Zero(
        static_cast<Eigen::Index>(dimension)
    );
    if (direction_count != 0) {
        const Eigen::JacobiSVD<Eigen::MatrixXd> svd(
            directions,
            Eigen::ComputeFullU | Eigen::ComputeFullV
        );
        const Eigen::VectorXd singular_values = svd.singularValues();
        for (Eigen::Index value = 0; value < singular_values.size(); ++value) {
            if (singular_values(value) > rank_threshold) {
                ++direction_rank;
            }
        }
        if (direction_rank != direction_count) {
            throw std::invalid_argument(
                "direction_matrix rows must have full row rank"
            );
        }
        const Eigen::Index invariant_count = static_cast<Eigen::Index>(
            dimension - direction_rank
        );
        if (invariant_count == 0) {
            invariants.resize(0, static_cast<Eigen::Index>(dimension));
        } else {
            invariants = svd.matrixV()
                .rightCols(invariant_count)
                .transpose();
            canonicalize_row_signs(invariants);
        }
        const Eigen::VectorXd rhs = -Eigen::Map<const Eigen::VectorXd>(
            driving_constants.data(),
            static_cast<Eigen::Index>(direction_count)
        );
        const Eigen::MatrixXd gram = directions * directions.transpose();
        const Eigen::LDLT<Eigen::MatrixXd> factor(gram);
        if (factor.info() != Eigen::Success) {
            throw std::invalid_argument(
                "direction_matrix normal equations are singular"
            );
        }
        reference_potentials = directions.transpose() * factor.solve(rhs);
        if (!reference_potentials.allFinite()) {
            throw std::invalid_argument(
                "reference-potential solve produced a nonfinite result"
            );
        }
    } else {
        // An empty direction declaration is an ordinary nonreactive problem:
        // every phase/species amount is an invariant and g is zero.
        direction_rank = 0;
    }

    const Eigen::VectorXd reference_vector = Eigen::Map<const Eigen::VectorXd>(
        phase_reference_amounts.values.data(),
        static_cast<Eigen::Index>(dimension)
    );
    const Eigen::VectorXd invariant_totals = invariants * reference_vector;
    const Eigen::MatrixXd direction_invariant_product =
        invariants * directions.transpose();
    const Eigen::MatrixXd invariant_gram =
        invariants * invariants.transpose();
    const Eigen::MatrixXd identity = Eigen::MatrixXd::Identity(
        invariants.rows(),
        invariants.rows()
    );

    std::vector<double> global_balance_residuals(balance_matrix.rows, 0.0);
    double global_balance_residual_inf = 0.0;
    for (std::size_t balance = 0; balance < balance_matrix.rows; ++balance) {
        double residual = 0.0;
        double scale = 0.0;
        for (std::size_t direction = 0; direction < direction_count; ++direction) {
            double change = 0.0;
            double row_scale = 0.0;
            for (std::size_t phase = 0; phase < phase_count; ++phase) {
                for (std::size_t species = 0; species < species_count; ++species) {
                    const double coefficient = direction_matrix(
                        direction,
                        phase * species_count + species
                    );
                    const double balance_value = balance_matrix(balance, species);
                    change += coefficient * balance_value;
                    row_scale += std::abs(coefficient * balance_value);
                }
            }
            residual = std::max(residual, std::abs(change));
            scale = std::max(scale, row_scale);
        }
        global_balance_residuals[balance] = residual;
        global_balance_residual_inf = std::max(global_balance_residual_inf, residual);
        if (residual > conservation_tolerance * std::max(1.0, scale)) {
            throw std::invalid_argument(
                "direction_matrix violates global balance conservation"
            );
        }
    }

    std::vector<double> phase_charge_residuals(phase_count, 0.0);
    double phase_charge_residual_inf = 0.0;
    for (std::size_t phase = 0; phase < phase_count; ++phase) {
        double residual = 0.0;
        double scale = 0.0;
        for (std::size_t direction = 0; direction < direction_count; ++direction) {
            double change = 0.0;
            double row_scale = 0.0;
            for (std::size_t species = 0; species < species_count; ++species) {
                const double contribution = direction_matrix(
                    direction,
                    phase * species_count + species
                ) * charges[species];
                change += contribution;
                row_scale += std::abs(contribution);
            }
            residual = std::max(residual, std::abs(change));
            scale = std::max(scale, row_scale);
        }
        phase_charge_residuals[phase] = residual;
        phase_charge_residual_inf = std::max(phase_charge_residual_inf, residual);
        if (residual > conservation_tolerance * std::max(1.0, scale)) {
            throw std::invalid_argument(
                "direction_matrix violates per-phase charge conservation"
            );
        }
    }

    std::vector<double> reference_phase_charge_residuals(phase_count, 0.0);
    double reference_neutrality_residual_inf = 0.0;
    for (std::size_t phase = 0; phase < phase_count; ++phase) {
        double residual = 0.0;
        double scale = 0.0;
        for (std::size_t species = 0; species < species_count; ++species) {
            const double contribution = phase_reference_amounts(
                phase,
                species
            ) * charges[species];
            residual += contribution;
            scale += std::abs(contribution);
        }
        reference_phase_charge_residuals[phase] = std::abs(residual);
        reference_neutrality_residual_inf = std::max(
            reference_neutrality_residual_inf,
            std::abs(residual)
        );
        if (std::abs(residual)
            > conservation_tolerance * std::max(1.0, scale)) {
            throw std::invalid_argument(
                "phase_reference_amounts must define neutral reference phases"
            );
        }
    }

    Eigen::VectorXd potential_residual;
    if (direction_count != 0) {
        potential_residual = directions * reference_potentials
            + Eigen::Map<const Eigen::VectorXd>(
                driving_constants.data(),
                static_cast<Eigen::Index>(direction_count)
            );
    }

    PhaseDirectionCompilation result;
    result.phase_ids = phase_ids;
    result.species_ids = species_ids;
    result.charges = charges;
    result.balance_matrix = balance_matrix;
    result.direction_matrix = direction_matrix;
    result.invariant_matrix = from_eigen(invariants);
    result.phase_reference_amounts = phase_reference_amounts.values;
    result.driving_constants = driving_constants;
    result.invariant_totals.assign(invariant_totals.data(), invariant_totals.data()
        + invariant_totals.size());
    result.reference_potentials.assign(
        reference_potentials.data(),
        reference_potentials.data() + reference_potentials.size()
    );
    result.global_balance_residuals = std::move(global_balance_residuals);
    result.phase_charge_residuals = std::move(phase_charge_residuals);
    result.reference_phase_charge_residuals =
        std::move(reference_phase_charge_residuals);
    result.direction_rank = direction_rank;
    result.invariant_rank = static_cast<std::size_t>(invariants.rows());
    result.direction_scale = direction_scale;
    result.rank_tolerance = rank_tolerance;
    result.conservation_tolerance = conservation_tolerance;
    const Eigen::MatrixXd invariant_orthonormal_residual =
        invariant_gram - identity;
    result.invariant_orthonormal_residual_inf =
        inf_norm(invariant_orthonormal_residual);
    result.invariant_direction_residual_inf = inf_norm(direction_invariant_product);
    result.potential_equation_residual_inf = inf_norm(potential_residual);
    result.minimum_norm_potential_norm = reference_potentials.norm();
    result.global_balance_residual_inf = global_balance_residual_inf;
    result.phase_charge_residual_inf = phase_charge_residual_inf;
    result.reference_neutrality_residual_inf = reference_neutrality_residual_inf;
    return result;
}

}  // namespace epcsaft_engine::equilibrium
