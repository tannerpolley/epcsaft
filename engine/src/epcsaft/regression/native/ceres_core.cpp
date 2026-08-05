#include "ceres_core.hpp"

#include <Eigen/Dense>
#include <Eigen/SVD>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace epcsaft_engine::regression::internal {
namespace {

class ExactCost final : public ceres::CostFunction {
public:
    ExactCost(const ProblemShape& shape, const ExactEvaluator& evaluator)
        : shape_(shape), evaluator_(evaluator) {
        set_num_residuals(static_cast<int>(shape_.residual_count));
        mutable_parameter_block_sizes()->push_back(
            static_cast<int>(shape_.variable_count())
        );
    }

    bool Evaluate(
        double const* const* parameters,
        double* residuals,
        double** jacobians
    ) const override {
        const bool jacobian_requested =
            jacobians != nullptr && jacobians[0] != nullptr;
        std::fill(
            residuals,
            residuals + shape_.residual_count,
            std::numeric_limits<double>::quiet_NaN()
        );
        if (jacobian_requested) {
            std::fill(
                jacobians[0],
                jacobians[0]
                    + shape_.residual_count * shape_.variable_count(),
                std::numeric_limits<double>::quiet_NaN()
            );
        }
        std::string failure;
        const bool usable = evaluator_(
            parameters[0],
            shape_.variable_count(),
            jacobian_requested,
            residuals,
            jacobian_requested ? jacobians[0] : nullptr,
            failure
        );
        if (!had_failure_) {
            failure_reason_ = std::move(failure);
        }
        if (!usable) {
            if (!had_failure_ && failure_reason_.empty()) {
                failure_reason_ = "exact evaluator callback rejected evaluation";
            }
            had_failure_ = true;
            return false;
        }
        const bool complete_residuals = std::all_of(
            residuals,
            residuals + shape_.residual_count,
            [](double value) { return std::isfinite(value); }
        );
        const bool complete_jacobian =
            !jacobian_requested
            || std::all_of(
                jacobians[0],
                jacobians[0]
                    + shape_.residual_count * shape_.variable_count(),
                [](double value) { return std::isfinite(value); }
            );
        if (!complete_residuals || !complete_jacobian) {
            if (!had_failure_) {
                failure_reason_ =
                    "exact evaluator left a nonfinite or incomplete residual "
                    "or Jacobian buffer";
            }
            had_failure_ = true;
            return false;
        }
        return true;
    }

    const std::string& failure_reason() const noexcept {
        return failure_reason_;
    }

    bool had_failure() const noexcept {
        return had_failure_;
    }

private:
    ProblemShape shape_;
    ExactEvaluator evaluator_;
    mutable std::string failure_reason_;
    mutable bool had_failure_{false};
};

MatrixDiagnostics diagnose(const Eigen::MatrixXd& matrix) {
    MatrixDiagnostics diagnostics{};
    if (matrix.rows() == 0 || matrix.cols() == 0) {
        diagnostics.condition_number = 1.0;
        return diagnostics;
    }
    const Eigen::JacobiSVD<Eigen::MatrixXd> svd(
        matrix, Eigen::ComputeThinU | Eigen::ComputeThinV
    );
    diagnostics.singular_values.assign(
        svd.singularValues().data(),
        svd.singularValues().data() + svd.singularValues().size()
    );
    if (diagnostics.singular_values.empty()) {
        return diagnostics;
    }
    const double threshold =
        100.0 * std::numeric_limits<double>::epsilon()
        * static_cast<double>(std::max(matrix.rows(), matrix.cols()))
        * diagnostics.singular_values.front();
    diagnostics.rank = static_cast<int>(std::count_if(
        diagnostics.singular_values.cbegin(),
        diagnostics.singular_values.cend(),
        [threshold](double value) { return value > threshold; }
    ));
    if (diagnostics.rank == matrix.cols()) {
        diagnostics.condition_number =
            diagnostics.singular_values.front()
            / diagnostics.singular_values[
                static_cast<std::size_t>(diagnostics.rank - 1)
            ];
    }
    return diagnostics;
}

void diagnose_jacobian_impl(
    const ProblemShape& shape,
    const std::vector<double>& row_major_jacobian,
    SolveResult& result
) {
    using RowMajorMatrix = Eigen::Matrix<
        double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor
    >;
    const Eigen::Map<const RowMajorMatrix> full(
        row_major_jacobian.data(),
        static_cast<Eigen::Index>(shape.residual_count),
        static_cast<Eigen::Index>(shape.variable_count())
    );
    result.full_jacobian = diagnose(full);
    Eigen::MatrixXd projected = full.leftCols(
        static_cast<Eigen::Index>(shape.fitted_count)
    );
    if (shape.lifted_count > 0) {
        const Eigen::MatrixXd lifted = full.rightCols(
            static_cast<Eigen::Index>(shape.lifted_count)
        );
        const Eigen::JacobiSVD<Eigen::MatrixXd> svd(
            lifted, Eigen::ComputeThinU
        );
        int nuisance_rank = 0;
        if (svd.singularValues().size() > 0) {
            const double threshold =
                100.0 * std::numeric_limits<double>::epsilon()
                * static_cast<double>(
                    std::max(
                        shape.residual_count, shape.variable_count()
                    )
                )
                * svd.singularValues()[0];
            nuisance_rank = static_cast<int>(std::count_if(
                svd.singularValues().data(),
                svd.singularValues().data() + svd.singularValues().size(),
                [threshold](double value) { return value > threshold; }
            ));
        }
        if (nuisance_rank > 0) {
            const Eigen::MatrixXd basis =
                svd.matrixU().leftCols(nuisance_rank);
            projected -= basis * (basis.transpose() * projected);
        }
    }
    result.projected_parameter_jacobian = diagnose(projected);
}

ceres::Solver::Options solver_options(const SolverControls& controls) {
    ceres::Solver::Options options;
    options.linear_solver_type = ceres::DENSE_QR;
    options.max_num_iterations = controls.maximum_iterations;
    options.max_solver_time_in_seconds =
        controls.maximum_solver_time_seconds;
    options.function_tolerance = controls.function_tolerance;
    options.gradient_tolerance = controls.gradient_tolerance;
    options.parameter_tolerance = controls.parameter_tolerance;
    options.logging_type = ceres::SILENT;
    options.num_threads = 1;
    return options;
}

}  // namespace

JacobianDiagnostics diagnose_jacobian(
    const ProblemShape& shape,
    const std::vector<double>& row_major_jacobian
) {
    if (shape.fitted_count == 0 || shape.residual_count == 0
        || row_major_jacobian.size()
            != shape.residual_count * shape.variable_count()) {
        throw std::invalid_argument(
            "Jacobian shape and row-major payload must be complete"
        );
    }
    SolveResult result{};
    diagnose_jacobian_impl(shape, row_major_jacobian, result);
    return {result.full_jacobian, result.projected_parameter_jacobian};
}

SolveResult solve(
    const ProblemShape& shape,
    const std::vector<double>& start,
    const std::vector<CoordinateBound>& bounds,
    const SolverControls& controls,
    const ExactEvaluator& evaluator
) {
    if (shape.fitted_count == 0 || shape.residual_count == 0) {
        throw std::invalid_argument(
            "fitted and residual counts must be positive"
        );
    }
    if (shape.residual_count < shape.variable_count()) {
        throw std::invalid_argument(
            "residual count must be at least the fitted-plus-lifted "
            "variable count"
        );
    }
    if (start.size() != shape.variable_count()
        || bounds.size() != shape.variable_count()) {
        throw std::invalid_argument(
            "start and bounds must match the fitted-plus-lifted "
            "variable count"
        );
    }
    for (std::size_t index = 0; index < bounds.size(); ++index) {
        if (!std::isfinite(start[index])
            || !std::isfinite(bounds[index].lower)
            || !std::isfinite(bounds[index].upper)
            || bounds[index].lower >= bounds[index].upper
            || start[index] < bounds[index].lower
            || start[index] > bounds[index].upper) {
            throw std::invalid_argument(
                "every start and ordered coordinate bound must be finite "
                "and consistent"
            );
        }
    }

    SolveResult result{};
    result.variables = start;
    ExactCost cost(shape, evaluator);
    ceres::Problem::Options problem_options;
    problem_options.cost_function_ownership = ceres::DO_NOT_TAKE_OWNERSHIP;
    ceres::Problem problem(problem_options);
    problem.AddResidualBlock(&cost, nullptr, result.variables.data());
    for (std::size_t index = 0; index < bounds.size(); ++index) {
        problem.SetParameterLowerBound(
            result.variables.data(),
            static_cast<int>(index),
            bounds[index].lower
        );
        problem.SetParameterUpperBound(
            result.variables.data(),
            static_cast<int>(index),
            bounds[index].upper
        );
    }
    ceres::Solve(
        solver_options(controls), &problem, &result.summary
    );
    result.failure_reason = cost.failure_reason();

    result.residuals.resize(shape.residual_count);
    result.jacobian.resize(
        shape.residual_count * shape.variable_count()
    );
    std::fill(
        result.residuals.begin(),
        result.residuals.end(),
        std::numeric_limits<double>::quiet_NaN()
    );
    std::fill(
        result.jacobian.begin(),
        result.jacobian.end(),
        std::numeric_limits<double>::quiet_NaN()
    );
    std::string final_failure;
    if (!evaluator(
            result.variables.data(),
            result.variables.size(),
            true,
            result.residuals.data(),
            result.jacobian.data(),
            final_failure
        )) {
        if (result.failure_reason.empty()) {
            result.failure_reason = std::move(final_failure);
        }
        return result;
    }
    if (!std::all_of(
            result.residuals.cbegin(),
            result.residuals.cend(),
            [](double value) { return std::isfinite(value); }
        )
        || !std::all_of(
            result.jacobian.cbegin(),
            result.jacobian.cend(),
            [](double value) { return std::isfinite(value); }
        )) {
        if (result.failure_reason.empty()) {
            result.failure_reason =
                "exact evaluator left a nonfinite or incomplete final "
                "residual or Jacobian buffer";
        }
        return result;
    }
    if (result.summary.IsSolutionUsable() && !cost.had_failure()) {
        result.failure_reason.clear();
    }
    diagnose_jacobian_impl(shape, result.jacobian, result);
    return result;
}

}  // namespace epcsaft_engine::regression::internal
