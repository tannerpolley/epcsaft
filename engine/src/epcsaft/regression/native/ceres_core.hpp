#ifndef EPCSAFT_ENGINE_REGRESSION_CERES_CORE_HPP
#define EPCSAFT_ENGINE_REGRESSION_CERES_CORE_HPP

#include <ceres/ceres.h>

#include <cstddef>
#include <functional>
#include <limits>
#include <string>
#include <vector>

namespace epcsaft_engine::regression::internal {

struct ProblemShape final {
    std::size_t fitted_count;
    std::size_t lifted_count;
    std::size_t residual_count;

    std::size_t variable_count() const noexcept {
        return fitted_count + lifted_count;
    }
};

struct CoordinateBound final {
    double lower;
    double upper;
};

struct SolverControls final {
    int maximum_iterations;
    double maximum_solver_time_seconds;
    double function_tolerance;
    double gradient_tolerance;
    double parameter_tolerance;
};

struct MatrixDiagnostics final {
    std::vector<double> singular_values;
    int rank{0};
    double condition_number{std::numeric_limits<double>::infinity()};
};

struct JacobianDiagnostics final {
    MatrixDiagnostics full;
    MatrixDiagnostics projected_parameters;
};

using ExactEvaluator = std::function<bool(
    const double* variables,
    std::size_t variable_count,
    bool jacobian_requested,
    double* residuals,
    double* row_major_jacobian,
    std::string& failure_reason
)>;

struct SolveResult final {
    ceres::Solver::Summary summary;
    std::vector<double> variables;
    std::vector<double> residuals;
    std::vector<double> jacobian;
    MatrixDiagnostics full_jacobian;
    MatrixDiagnostics projected_parameter_jacobian;
    std::string failure_reason;
};

SolveResult solve(
    const ProblemShape& shape,
    const std::vector<double>& start,
    const std::vector<CoordinateBound>& bounds,
    const SolverControls& controls,
    const ExactEvaluator& evaluator
);

JacobianDiagnostics diagnose_jacobian(
    const ProblemShape& shape,
    const std::vector<double>& row_major_jacobian
);

}  // namespace epcsaft_engine::regression::internal

#endif
