#include "ceres_core.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void require(bool condition) {
    if (!condition) {
        throw std::runtime_error("ceres core test assertion failed");
    }
}

}  // namespace

void test_evaluator_fit();

int main() {
    using epcsaft_engine::regression::internal::CoordinateBound;
    using epcsaft_engine::regression::internal::ProblemShape;
    using epcsaft_engine::regression::internal::SolverControls;

    constexpr std::array<std::size_t, 3> slot_to_parameter{{0, 1, 0}};
    constexpr std::array<std::array<double, 3>, 4> slot_jacobian{{
        {{0.5, 0.0, 0.5}},
        {{0.0, 1.0, 0.0}},
        {{0.5, 1.0, 0.5}},
        {{1.0, -1.0, 1.0}},
    }};
    constexpr std::array<double, 4> lifted_jacobian{{1.0, 2.0, 3.0, 1.0}};
    std::array<std::array<double, 3>, 4> jacobian{};
    for (std::size_t row = 0; row < jacobian.size(); ++row) {
        for (std::size_t slot = 0; slot < slot_to_parameter.size(); ++slot) {
            jacobian[row][slot_to_parameter[slot]] +=
                slot_jacobian[row][slot];
        }
        jacobian[row][2] = lifted_jacobian[row];
    }
    constexpr std::array<double, 3> expected{{1.0, -2.0, 0.5}};
    std::array<double, 4> target{};
    for (std::size_t row = 0; row < target.size(); ++row) {
        for (std::size_t column = 0; column < expected.size(); ++column) {
            target[row] += jacobian[row][column] * expected[column];
        }
    }
    std::size_t value_only_calls = 0;
    std::size_t jacobian_calls = 0;
    const auto evaluator = [&](
        const double* variables,
        std::size_t variable_count,
        bool jacobian_requested,
        double* residuals,
        double* row_major_jacobian,
        std::string& failure_reason
    ) {
        if (variable_count != expected.size()) {
            failure_reason = "unexpected variable count";
            return false;
        }
        jacobian_requested ? ++jacobian_calls : ++value_only_calls;
        for (std::size_t row = 0; row < target.size(); ++row) {
            residuals[row] = -target[row];
            for (std::size_t column = 0; column < expected.size(); ++column) {
                residuals[row] +=
                    jacobian[row][column] * variables[column];
                if (jacobian_requested) {
                    row_major_jacobian[
                        row * expected.size() + column
                    ] = jacobian[row][column];
                }
            }
        }
        failure_reason.clear();
        return true;
    };
    const ProblemShape shape{2, 1, 4};
    const std::vector<CoordinateBound> bounds(
        shape.variable_count(), CoordinateBound{-10.0, 10.0}
    );
    const SolverControls controls{50, 5.0, 1.0e-14, 1.0e-14, 1.0e-14};
    const auto primary = epcsaft_engine::regression::internal::solve(
        shape, {0.0, 0.0, 0.0}, bounds, controls, evaluator
    );
    const auto confirmation = epcsaft_engine::regression::internal::solve(
        shape, {-3.0, 4.0, -2.0}, bounds, controls, evaluator
    );
    require(primary.summary.IsSolutionUsable());
    require(confirmation.summary.IsSolutionUsable());
    for (std::size_t index = 0; index < expected.size(); ++index) {
        require(std::abs(primary.variables[index] - expected[index]) < 1.0e-10);
        require(
            std::abs(confirmation.variables[index] - expected[index])
            < 1.0e-10
        );
    }
    require(primary.full_jacobian.rank == 3);
    require(primary.projected_parameter_jacobian.rank == 2);
    const auto preflight = epcsaft_engine::regression::internal::diagnose_jacobian(
        shape, primary.jacobian
    );
    require(preflight.full.rank == primary.full_jacobian.rank);
    require(
        preflight.projected_parameters.rank
        == primary.projected_parameter_jacobian.rank
    );
    require(
        preflight.full.condition_number
        == primary.full_jacobian.condition_number
    );
    require(primary.jacobian.size() == 12);
    require(value_only_calls > 0);
    require(jacobian_calls > 0);

    bool structural_gate_failed = false;
    try {
        static_cast<void>(epcsaft_engine::regression::internal::solve(
            ProblemShape{2, 1, 2},
            {0.0, 0.0, 0.0},
            bounds,
            controls,
            evaluator
        ));
    } catch (const std::invalid_argument&) {
        structural_gate_failed = true;
    }
    require(structural_gate_failed);

    const auto missing_second_column = [&](
        const double* variables,
        std::size_t variable_count,
        bool jacobian_requested,
        double* residuals,
        double* row_major_jacobian,
        std::string& failure_reason
    ) {
        if (variable_count != 2) {
            failure_reason = "unexpected variable count";
            return false;
        }
        for (std::size_t row = 0; row < 2; ++row) {
            residuals[row] = variables[0] - 1.0;
            if (jacobian_requested) {
                row_major_jacobian[2 * row] = 1.0;
                row_major_jacobian[2 * row + 1] = 0.0;
            }
        }
        failure_reason.clear();
        return true;
    };
    const auto deficient = epcsaft_engine::regression::internal::solve(
        ProblemShape{2, 0, 2},
        {0.0, 0.0},
        std::vector<CoordinateBound>(2, {-10.0, 10.0}),
        controls,
        missing_second_column
    );
    require(deficient.full_jacobian.rank == 1);
    require(deficient.projected_parameter_jacobian.rank == 1);

    const auto incomplete = [](
        const double* variables,
        std::size_t variable_count,
        bool jacobian_requested,
        double* residuals,
        double* row_major_jacobian,
        std::string& failure_reason
    ) {
        if (variable_count != 1) {
            failure_reason = "unexpected variable count";
            return false;
        }
        residuals[0] = variables[0] - 1.0;
        if (jacobian_requested) {
            static_cast<void>(row_major_jacobian);
        }
        failure_reason.clear();
        return true;
    };
    const auto incomplete_result = epcsaft_engine::regression::internal::solve(
        ProblemShape{1, 0, 1},
        {0.0},
        std::vector<CoordinateBound>(1, {-10.0, 10.0}),
        controls,
        incomplete
    );
    require(!incomplete_result.summary.IsSolutionUsable());
    require(
        incomplete_result.failure_reason.find("incomplete")
        != std::string::npos
    );
    int transient_calls = 0;
    const auto transient_failure = [&](
        const double* variables,
        std::size_t variable_count,
        bool jacobian_requested,
        double* residuals,
        double* row_major_jacobian,
        std::string& failure_reason
    ) {
        if (variable_count != 1) {
            failure_reason = "unexpected variable count";
            return false;
        }
        if (transient_calls++ == 0) {
            failure_reason = "transient evaluator failure";
            return false;
        }
        residuals[0] = variables[0] - 1.0;
        if (jacobian_requested) {
            row_major_jacobian[0] = 1.0;
        }
        failure_reason.clear();
        return true;
    };
    const auto transient_result = epcsaft_engine::regression::internal::solve(
        ProblemShape{1, 0, 1},
        {0.0},
        std::vector<CoordinateBound>(1, {-10.0, 10.0}),
        controls,
        transient_failure
    );
    require(transient_calls > 1);
    require(!transient_result.summary.IsSolutionUsable());
    require(
        transient_result.failure_reason.find("transient evaluator failure")
        != std::string::npos
    );
    test_evaluator_fit();
    return 0;
}
