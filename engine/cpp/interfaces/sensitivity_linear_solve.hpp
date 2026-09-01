#pragma once
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>
#include <Eigen/Dense>
namespace epcsaft_engine::interfaces {
inline constexpr double kSolvedStateSensitivityConditionNumberMax = 1.0e8;
inline double solved_state_sensitivity_backward_error_max(std::size_t n) {
  return 8192.0 * std::numeric_limits<double>::epsilon() * static_cast<double>(n);
}
struct RankRevealingLinearSolve {
  std::size_t rank = 0;
  double condition_number = std::numeric_limits<double>::infinity();
  double backward_error = std::numeric_limits<double>::infinity();
  std::vector<double> solution;
};
inline RankRevealingLinearSolve solve_rank_revealing_linear_system(
    const std::vector<double>& f_y, const std::vector<double>& f_p,
    std::size_t dimension, std::size_t parameter_count) {
  RankRevealingLinearSolve result;
  const auto finite = [](const std::vector<double>& values) {
    return std::all_of(values.begin(), values.end(), [](double value) {
      return std::isfinite(value);
    });
  };
  if (dimension == 0 || parameter_count == 0 ||
      f_y.size() != dimension * dimension ||
      f_p.size() != dimension * parameter_count || !finite(f_y) || !finite(f_p))
    return result;
  using Matrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
  Eigen::Map<const Matrix> matrix(f_y.data(), dimension, dimension);
  Eigen::Map<const Matrix> partials(f_p.data(), dimension, parameter_count);
  Eigen::VectorXd row_scale = Eigen::VectorXd::Ones(dimension);
  Eigen::VectorXd column_scale = Eigen::VectorXd::Ones(dimension);
  Matrix scaled = matrix;
  for (int pass = 0; pass < 4; ++pass) {
    for (Eigen::Index row = 0; row < scaled.rows(); ++row) {
      const double norm = scaled.row(row).cwiseAbs().maxCoeff();
      if (norm > 0.0 && std::isfinite(norm)) { row_scale[row] /= norm; scaled.row(row) /= norm; }
    }
    for (Eigen::Index column = 0; column < scaled.cols(); ++column) {
      const double norm = scaled.col(column).cwiseAbs().maxCoeff();
      if (norm > 0.0 && std::isfinite(norm)) { column_scale[column] /= norm; scaled.col(column) /= norm; }
    }
  }
  Eigen::JacobiSVD<Matrix> decomposition(scaled, Eigen::ComputeFullU | Eigen::ComputeFullV);
  decomposition.setThreshold(4096.0 * std::numeric_limits<double>::epsilon() * static_cast<double>(dimension));
  result.rank = static_cast<std::size_t>(decomposition.rank());
  const auto singular_values = decomposition.singularValues();
  if (result.rank == dimension && singular_values.size() != 0 && singular_values[dimension - 1] > 0.0)
    result.condition_number = singular_values[0] / singular_values[dimension - 1];
  if (result.rank != dimension || !std::isfinite(result.condition_number) ||
      result.condition_number > kSolvedStateSensitivityConditionNumberMax)
    return result;
  const Matrix scaled_rhs = (-row_scale).asDiagonal() * partials;
  const Matrix scaled_solution = decomposition.solve(scaled_rhs);
  const Matrix solution = column_scale.asDiagonal() * scaled_solution;
  const Matrix residual = matrix * solution + partials;
  const double numerator = residual.cwiseAbs().maxCoeff();
  const double denominator = std::max(1.0, matrix.cwiseAbs().maxCoeff() * solution.cwiseAbs().maxCoeff() + partials.cwiseAbs().maxCoeff());
  result.backward_error = numerator / denominator;
  if (!std::isfinite(result.backward_error) ||
      result.backward_error > solved_state_sensitivity_backward_error_max(dimension) ||
      !solution.allFinite())
    return result;
  result.solution.assign(solution.data(), solution.data() + solution.size());
  return result;
}
}  // namespace epcsaft_engine::interfaces
