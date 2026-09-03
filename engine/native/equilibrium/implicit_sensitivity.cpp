#include "equilibrium/gibbs_minimizer.hpp"

#include "core/trace.hpp"
#include "core/numerical_policy.hpp"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace epcsaft::equilibrium {
namespace {

Eigen::MatrixXd dense_matrix(
    const std::vector<double>& values, std::size_t rows, std::size_t columns
) {
    if (values.size() != rows * columns) throw std::invalid_argument("matrix size");
    Eigen::MatrixXd result(static_cast<Eigen::Index>(rows),
                           static_cast<Eigen::Index>(columns));
    for (std::size_t row = 0; row < rows; ++row) {
        for (std::size_t column = 0; column < columns; ++column) {
            result(static_cast<Eigen::Index>(row), static_cast<Eigen::Index>(column)) =
                values[row * columns + column];
        }
    }
    return result;
}

double condition_estimate(const Eigen::MatrixXd& matrix) {
    if (matrix.rows() == 0 || matrix.cols() == 0) return 0.0;
    const auto singular = matrix.jacobiSvd().singularValues();
    if (singular.size() == 0 || singular(singular.size() - 1) <= 0.0) {
        return std::numeric_limits<double>::infinity();
    }
    return singular(0) / singular(singular.size() - 1);
}

std::size_t parameter_count(
    const std::vector<double>& parameter_rows, std::size_t rows
) {
    if (rows == 0 || parameter_rows.empty() || parameter_rows.size() % rows != 0) return 0;
    return parameter_rows.size() / rows;
}

void validate_parameter_rows(
    const std::vector<double>& parameter_rows, std::size_t rows
) {
    if (!std::all_of(parameter_rows.begin(), parameter_rows.end(),
                     [](double value) { return std::isfinite(value); })) {
        throw std::invalid_argument("invalid parameter derivative");
    }
    if (!parameter_rows.empty() && (rows == 0 || parameter_rows.size() % rows != 0)) {
        throw std::invalid_argument("parameter derivative shape");
    }
}

std::vector<double> flatten(const Eigen::MatrixXd& matrix) {
    std::vector<double> result(
        static_cast<std::size_t>(matrix.rows() * matrix.cols()), 0.0);
    for (Eigen::Index row = 0; row < matrix.rows(); ++row) {
        for (Eigen::Index column = 0; column < matrix.cols(); ++column) {
            result[static_cast<std::size_t>(row * matrix.cols() + column)] =
                matrix(row, column);
        }
    }
    return result;
}

bool valid_point(const std::vector<double>& point) {
    return !point.empty() && std::all_of(
        point.begin(), point.end(), [](double value) { return std::isfinite(value); });
}

double maximum_abs(const Eigen::MatrixXd& matrix) {
    double result = 0.0;
    for (Eigen::Index row = 0; row < matrix.rows(); ++row) {
        for (Eigen::Index column = 0; column < matrix.cols(); ++column) {
            result = std::max(result, std::abs(matrix(row, column)));
        }
    }
    return result;
}

bool residual_is_small(const Eigen::MatrixXd& matrix,
                       const Eigen::MatrixXd& solution,
                       const Eigen::MatrixXd& right) {
    const double scale = std::max(1.0, maximum_abs(right));
    return maximum_abs(matrix * solution - right)
        <= policy::EQUILIBRIUM_ROW_RESIDUAL.value * scale;
}

bool square_and_numeric(const Eigen::MatrixXd& matrix) {
    if (matrix.rows() != matrix.cols()) return false;
    for (Eigen::Index row = 0; row < matrix.rows(); ++row) {
        for (Eigen::Index column = 0; column < matrix.cols(); ++column) {
            if (!std::isfinite(matrix(row, column))) return false;
        }
    }
    return true;
}

double relative_residual(const Eigen::MatrixXd& matrix,
                         const Eigen::MatrixXd& solution,
                         const Eigen::MatrixXd& right) {
    const double scale = std::max(1.0, maximum_abs(right));
    return maximum_abs(matrix * solution - right) / scale;
}

}  // namespace

SensitivityInfo inspect_sensitivity(
    const CompiledSystem& compiled,
    const std::vector<double>& point,
    const std::vector<double>& parameter_rows
) {
    core::Scope scope("inspect_implicit_sensitivity");
    const std::size_t rows = compiled.rows.size();
    const std::size_t dimension = point.size();
    if (!valid_point(point)) throw std::invalid_argument("invalid sensitivity point");
    validate_parameter_rows(parameter_rows, rows);
    const auto jacobian = compiled.jacobian(point);
    if (jacobian.size() != rows * dimension) {
        throw std::invalid_argument("equilibrium Jacobian shape");
    }
    const auto matrix = dense_matrix(jacobian, rows, dimension);
    const auto lu = matrix.fullPivLu();
    SensitivityInfo result;
    result.equations = rows;
    result.parameters = parameter_count(parameter_rows, rows);
    result.rank = static_cast<std::size_t>(lu.rank());
    result.condition_estimate = condition_estimate(matrix);
    result.nonsingular = square_and_numeric(matrix) && result.rank == dimension;
    scope.probe("equations", static_cast<double>(rows));
    scope.probe("dimension", static_cast<double>(dimension));
    scope.probe("condition", result.condition_estimate);
    scope.probe("jacobian_entries", static_cast<double>(jacobian.size()));
    scope.probe("parameter_entries", static_cast<double>(parameter_rows.size()));
    scope.probe("parameters", static_cast<double>(result.parameters));
    scope.probe("nonsingular", result.nonsingular ? 1.0 : 0.0);
    scope.probe("rank", static_cast<double>(result.rank));
    scope.probe("condition", result.condition_estimate);
    return result;
}

std::vector<double> solve_implicit(
    const CompiledSystem& compiled,
    const std::vector<double>& point,
    const std::vector<double>& parameter_rows
) {
    core::Scope scope("implicit_sensitivity");
    const std::size_t rows = compiled.rows.size();
    const std::size_t dimension = point.size();
    if (!valid_point(point)) throw std::invalid_argument("invalid sensitivity point");
    validate_parameter_rows(parameter_rows, rows);
    const std::size_t parameters = parameter_count(parameter_rows, rows);
    if (parameters == 0) return {};
    const auto jacobian = compiled.jacobian(point);
    if (jacobian.size() != rows * dimension || rows != dimension) return {};
    const auto info = inspect_sensitivity(compiled, point, parameter_rows);
    if (!info.nonsingular) return {};
    const auto matrix = dense_matrix(jacobian, rows, dimension);
    Eigen::MatrixXd right(static_cast<Eigen::Index>(rows),
                          static_cast<Eigen::Index>(parameters));
    for (std::size_t row = 0; row < rows; ++row) {
        for (std::size_t parameter = 0; parameter < parameters; ++parameter) {
            right(static_cast<Eigen::Index>(row), static_cast<Eigen::Index>(parameter)) =
                -parameter_rows[row * parameters + parameter];
        }
    }
    const Eigen::MatrixXd result = matrix.fullPivLu().solve(right);
    if (!result.allFinite() || !residual_is_small(matrix, result, right)) return {};
    scope.probe("parameters", static_cast<double>(parameters));
    scope.probe("max_solution", maximum_abs(result));
    scope.probe("relative_residual", relative_residual(matrix, result, right));
    scope.probe("equations", static_cast<double>(rows));
    scope.probe("dimension", static_cast<double>(dimension));
    return flatten(result);
}

std::vector<double> sensitivity(
    const CompiledSystem& compiled,
    const std::vector<double>& point,
    const std::vector<double>& lever
) {
    core::Scope scope("sensitivity");
    const std::size_t levers = 2 + compiled.model->parameters().component_count();
    if (lever.size() != levers)
        throw std::invalid_argument("lever must contain T, P, and feed weights");
    const auto F_theta = compiled.lever_jacobian(point);
    if (F_theta.size() != compiled.rows.size() * levers) {
        throw std::invalid_argument("lever Jacobian shape");
    }
    std::vector<double> selected(compiled.rows.size(), 0.0);
    for (std::size_t row = 0; row < compiled.rows.size(); ++row) {
        double F_theta_value = 0.0;
        for (std::size_t column = 0; column < levers; ++column)
            F_theta_value += F_theta[levers * row + column] * lever[column];
        selected[row] = -(-F_theta_value);
    }
    return solve_implicit(compiled, point, selected);
}

}  // namespace epcsaft::equilibrium
