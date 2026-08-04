#include "evaluator_fit.hpp"

#include <algorithm>
#include <climits>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string_view>
#include <utility>

namespace epcsaft_engine::regression::evaluator {
namespace {

struct Batch final {
    std::vector<double> values;
    std::vector<double> jacobian;
    std::vector<epcsaft_engine_regression_evaluator_row_result_v1> rows;
};

std::string bounded_text(
    const char* value, std::size_t capacity, const char* field
);

std::string text(const char* value, const char* field) {
    const std::string result = bounded_text(
        value,
        EPCSAFT_ENGINE_REGRESSION_EVALUATOR_V1_TEXT_CAPACITY,
        field
    );
    if (result.empty()) {
        throw std::invalid_argument(std::string(field) + " is missing");
    }
    return result;
}

std::string bounded_text(
    const char* value, std::size_t capacity, const char* field
) {
    if (value == nullptr) {
        throw std::runtime_error(std::string(field) + " is missing");
    }
    const std::size_t length = strnlen(value, capacity);
    if (length == capacity) {
        throw std::runtime_error(std::string(field) + " is not terminated");
    }
    return std::string(value, length);
}

std::string required_bounded_text(
    const char* value, std::size_t capacity, const char* field
) {
    const std::string result = bounded_text(value, capacity, field);
    if (result.empty()) {
        throw std::runtime_error(std::string(field) + " is empty");
    }
    return result;
}

void require_equal(
    const char* actual, const std::string& expected, const char* field
) {
    if (text(actual, field) != expected) {
        throw std::invalid_argument(std::string(field) + " does not match");
    }
}

void require_table(
    const char* const* values, std::size_t count, const char* field
) {
    if (values == nullptr) {
        throw std::invalid_argument(std::string(field) + " table is missing");
    }
    for (std::size_t index = 0; index < count; ++index) {
        if (bounded_text(
                values[index],
                EPCSAFT_ENGINE_REGRESSION_EVALUATOR_V1_TEXT_CAPACITY,
                field
            ).empty()) {
            throw std::invalid_argument(
                std::string(field) + " table is incomplete"
            );
        }
    }
}

std::size_t checked_product(
    std::size_t left, std::size_t right, const char* field
) {
    if (right != 0 && left > std::numeric_limits<std::size_t>::max() / right) {
        throw std::invalid_argument(std::string(field) + " size overflows");
    }
    return left * right;
}

Batch evaluate(
    const epcsaft_engine_regression_evaluator_sdk_v1& sdk,
    const Problem& problem,
    const std::vector<double>& physical_parameters,
    bool with_jacobian
) {
    Batch batch;
    batch.values.assign(
        problem.rows.size(), std::numeric_limits<double>::quiet_NaN()
    );
    if (with_jacobian) {
        batch.jacobian.assign(
            checked_product(
                problem.rows.size(),
                problem.parameters.size(),
                "evaluator Jacobian"
            ),
            std::numeric_limits<double>::quiet_NaN()
        );
    }
    batch.rows.resize(problem.rows.size());
    for (auto& row : batch.rows) {
        row = {};
        row.struct_size = sizeof(row);
    }
    epcsaft_engine_regression_evaluator_result_v1 output{};
    output.struct_size = sizeof(output);
    output.status = EPCSAFT_ENGINE_REGRESSION_EVALUATOR_STATUS_UNAVAILABLE_V1;
    output.row_count = problem.rows.size();
    output.parameter_count = problem.parameters.size();
    output.value_capacity = batch.values.size();
    output.jacobian_capacity = batch.jacobian.size();
    output.row_result_capacity = batch.rows.size();
    output.request_mode = -1;
    output.values = batch.values.data();
    output.jacobian = with_jacobian ? batch.jacobian.data() : nullptr;
    output.row_results = batch.rows.data();
    const int request_mode =
        with_jacobian
        ? EPCSAFT_ENGINE_REGRESSION_EVALUATOR_REQUEST_VALUES_AND_JACOBIAN_V1
        : EPCSAFT_ENGINE_REGRESSION_EVALUATOR_REQUEST_VALUES_ONLY_V1;
    const int returned = sdk.evaluate(
        sdk.model_context,
        physical_parameters.data(),
        physical_parameters.size(),
        request_mode,
        &output
    );
    if (output.struct_size != sizeof(output)
        || output.row_count != problem.rows.size()
        || output.parameter_count != problem.parameters.size()
        || output.value_capacity != batch.values.size()
        || output.jacobian_capacity != batch.jacobian.size()
        || output.row_result_capacity != batch.rows.size()
        || output.values != batch.values.data()
        || output.jacobian
            != (with_jacobian ? batch.jacobian.data() : nullptr)
        || output.row_results != batch.rows.data()) {
        throw std::runtime_error(
            "evaluator callback mutated the result shape or buffers"
        );
    }
    if (returned != output.status) {
        throw std::runtime_error(
            "evaluator callback returned inconsistent aggregate status"
        );
    }
    if (output.request_mode != request_mode) {
        throw std::runtime_error(
            "evaluator callback returned an inconsistent request mode"
        );
    }
    const std::string aggregate_error = bounded_text(
        output.error, sizeof(output.error), "evaluator aggregate error"
    );
    if (returned != EPCSAFT_ENGINE_REGRESSION_EVALUATOR_STATUS_OK_V1) {
        std::string reason = aggregate_error.empty()
            ? "evaluator batch was unavailable"
            : aggregate_error;
        for (std::size_t index = 0; index < batch.rows.size(); ++index) {
            const std::string row_reason = bounded_text(
                batch.rows[index].reason,
                sizeof(batch.rows[index].reason),
                "evaluator row reason"
            );
            if (!row_reason.empty()) {
                reason += "; " + problem.rows[index].id + ": " + row_reason;
            }
        }
        throw std::runtime_error(reason);
    }
    if (!aggregate_error.empty()) {
        throw std::runtime_error(
            "successful evaluator batch returned an aggregate error"
        );
    }
    if (bounded_text(
            output.provider_parameter_fingerprint,
            sizeof(output.provider_parameter_fingerprint),
            "result Provider parameter fingerprint"
        ) != problem.metadata.provider_parameter_fingerprint
        || bounded_text(
               output.artifact_identity,
               sizeof(output.artifact_identity),
               "result artifact identity"
           ) != problem.metadata.artifact_identity) {
        throw std::runtime_error(
            "evaluator result identity changed during the solve"
        );
    }
    for (std::size_t index = 0; index < batch.rows.size(); ++index) {
        const auto& certificate = batch.rows[index];
        const auto& spec = problem.rows[index];
        if (certificate.struct_size != sizeof(certificate)
            || certificate.status
                != EPCSAFT_ENGINE_REGRESSION_EVALUATOR_STATUS_OK_V1) {
            throw std::runtime_error(
                spec.id + ": evaluator row was not accepted"
            );
        }
        const std::string reason = bounded_text(
            certificate.reason, sizeof(certificate.reason), "row reason"
        );
        if (!reason.empty()) {
            throw std::runtime_error(
                spec.id + ": successful evaluator row returned a reason"
            );
        }
        const std::string solver_status = required_bounded_text(
            certificate.solver_status,
            sizeof(certificate.solver_status),
            "row solver status"
        );
        const std::string numerical_status = required_bounded_text(
            certificate.numerical_status,
            sizeof(certificate.numerical_status),
            "row numerical status"
        );
        const std::string physical_status = required_bounded_text(
            certificate.physical_status,
            sizeof(certificate.physical_status),
            "row physical status"
        );
        const std::string derivative_status = required_bounded_text(
            certificate.derivative_status,
            sizeof(certificate.derivative_status),
            "row derivative status"
        );
        if (solver_status
                != EPCSAFT_ENGINE_REGRESSION_EVALUATOR_V1_SOLVER_STATUS_SOLVE_SUCCEEDED
            || numerical_status
                != EPCSAFT_ENGINE_REGRESSION_EVALUATOR_V1_NUMERICAL_STATUS_PASSED
            || physical_status
                != EPCSAFT_ENGINE_REGRESSION_EVALUATOR_V1_PHYSICAL_STATUS_PASSED
            || (
                with_jacobian
                && derivative_status
                    != EPCSAFT_ENGINE_REGRESSION_EVALUATOR_V1_DERIVATIVE_STATUS_AVAILABLE
            )) {
            throw std::runtime_error(
                spec.id + ": evaluator row acceptance status is not canonical"
            );
        }
        static_cast<void>(required_bounded_text(
            certificate.chart_topology,
            sizeof(certificate.chart_topology),
            "row chart topology"
        ));
        if (bounded_text(
                certificate.provider_topology_fingerprint,
                sizeof(certificate.provider_topology_fingerprint),
                "row Provider topology fingerprint"
            ) != problem.metadata.provider_topology_fingerprint) {
            throw std::runtime_error(
                spec.id + ": Provider topology changed during evaluation"
            );
        }
        if (certificate.kkt_dimension == 0
            || certificate.kkt_rank != certificate.kkt_dimension
            || !std::isfinite(certificate.kkt_condition_number_inf)
            || certificate.kkt_condition_number_inf <= 0.0
            || certificate.kkt_condition_number_inf
                > problem.maximum_condition_number) {
            throw std::runtime_error(
                spec.id + ": upstream sensitivity certificate is singular "
                "or ill-conditioned"
            );
        }
        if (!std::isfinite(batch.values[index])
            || batch.values[index] <= 0.0) {
            throw std::runtime_error(
                spec.id + ": positive evaluator returned an invalid value"
            );
        }
    }
    if (with_jacobian
        && !std::all_of(
            batch.jacobian.cbegin(),
            batch.jacobian.cend(),
            [](double value) { return std::isfinite(value); }
        )) {
        throw std::runtime_error(
            "evaluator returned an incomplete or nonfinite Jacobian"
        );
    }
    return batch;
}

std::vector<double> to_physical(
    const Problem& problem,
    const double* solver_parameters,
    std::size_t count
) {
    if (count != problem.parameters.size()) {
        throw std::invalid_argument(
            "solver parameter count does not match the evaluator contract"
        );
    }
    std::vector<double> physical(count);
    for (std::size_t index = 0; index < count; ++index) {
        physical[index] =
            problem.parameters[index].origin
            + problem.parameters[index].scale * solver_parameters[index];
        if (!std::isfinite(physical[index])) {
            throw std::runtime_error(
                "affine parameter transform returned a nonfinite value"
            );
        }
    }
    return physical;
}

std::vector<std::size_t> training_indices(const Problem& problem) {
    std::vector<std::size_t> indices;
    for (std::size_t index = 0; index < problem.rows.size(); ++index) {
        if (problem.rows[index].partition == "training") {
            indices.push_back(index);
        }
    }
    return indices;
}

bool evaluate_training(
    const epcsaft_engine_regression_evaluator_sdk_v1& sdk,
    const Problem& problem,
    const std::vector<std::size_t>& training,
    const double* variables,
    std::size_t count,
    bool jacobian_requested,
    double* residuals,
    double* jacobian,
    std::string& failure_reason
) {
    try {
        const auto physical = to_physical(problem, variables, count);
        const Batch batch = evaluate(
            sdk, problem, physical, jacobian_requested
        );
        for (std::size_t output = 0; output < training.size(); ++output) {
            const std::size_t row_index = training[output];
            const auto& row = problem.rows[row_index];
            const double value = batch.values[row_index];
            residuals[output] =
                row.transform == "natural_log"
                ? (std::log(value) - std::log(row.observed)) / row.scale
                : (value - row.observed) / row.scale;
            if (jacobian_requested) {
                for (std::size_t parameter = 0;
                     parameter < problem.parameters.size();
                     ++parameter) {
                    const double physical_derivative =
                        batch.jacobian[
                            row_index * problem.parameters.size()
                            + parameter
                        ];
                    jacobian[
                        output * problem.parameters.size() + parameter
                    ] =
                        physical_derivative
                        * problem.parameters[parameter].scale
                        / row.scale
                        / (
                            row.transform == "natural_log"
                            ? value
                            : 1.0
                        );
                }
            }
        }
        failure_reason.clear();
        return true;
    } catch (const std::exception& error) {
        failure_reason = error.what();
        return false;
    } catch (...) {
        failure_reason = "unknown evaluator adapter failure";
        return false;
    }
}

}  // namespace

void validate_contract(
    const epcsaft_engine_regression_evaluator_sdk_v1& sdk,
    const Problem& problem
) {
    if (sdk.abi_version != EPCSAFT_ENGINE_REGRESSION_EVALUATOR_V1_ABI_VERSION
        || sdk.table_size < sizeof(sdk)
        || sdk.result_size != sizeof(epcsaft_engine_regression_evaluator_result_v1)
        || sdk.row_result_size
            != sizeof(epcsaft_engine_regression_evaluator_row_result_v1)
        || sdk.model_context == nullptr || sdk.evaluate == nullptr) {
        throw std::invalid_argument(
            "evaluator handle ABI, sizes, context, or callback is invalid"
        );
    }
    if (sdk.row_count > static_cast<std::size_t>(INT_MAX)
        || sdk.parameter_count > static_cast<std::size_t>(INT_MAX)) {
        throw std::invalid_argument(
            "evaluator row or parameter count exceeds Ceres limits"
        );
    }
    if (problem.parameters.empty() || problem.rows.empty()
        || sdk.parameter_count != problem.parameters.size()
        || sdk.row_count != problem.rows.size()) {
        throw std::invalid_argument(
            "evaluator row or parameter shape does not match the problem"
        );
    }
    static_cast<void>(checked_product(
        sdk.row_count, sdk.parameter_count, "evaluator Jacobian"
    ));
    require_table(sdk.row_ids, sdk.row_count, "row identity");
    require_table(sdk.state_ids, sdk.row_count, "state identity");
    require_table(
        sdk.state_schema_ids, sdk.row_count, "state schema identity"
    );
    require_table(
        sdk.observation_source_ids, sdk.row_count, "observation source"
    );
    require_table(sdk.source_locators, sdk.row_count, "source locator");
    require_table(sdk.primitive_ids, sdk.row_count, "primitive identity");
    require_table(sdk.primitive_units, sdk.row_count, "primitive unit");
    require_table(sdk.transform_ids, sdk.row_count, "transform identity");
    require_table(sdk.reference_ids, sdk.row_count, "reference identity");
    require_table(
        sdk.reference_fingerprints, sdk.row_count, "reference fingerprint"
    );
    require_table(
        sdk.parameter_ids, sdk.parameter_count, "parameter identity"
    );
    require_table(
        sdk.parameter_units, sdk.parameter_count, "parameter unit"
    );
    require_equal(
        sdk.evaluator_identity,
        problem.metadata.evaluator_identity,
        "evaluator identity"
    );
    require_equal(
        sdk.capability_id,
        problem.metadata.capability_id,
        "capability identity"
    );
    require_equal(
        sdk.capability_fingerprint,
        problem.metadata.capability_fingerprint,
        "capability fingerprint"
    );
    require_equal(
        sdk.provider_artifact_identity,
        problem.metadata.provider_artifact_identity,
        "Provider artifact identity"
    );
    require_equal(
        sdk.owner_artifact_identity,
        problem.metadata.owner_artifact_identity,
        "owner artifact identity"
    );
    require_equal(
        sdk.contract_fingerprint,
        problem.metadata.contract_fingerprint,
        "contract fingerprint"
    );
    require_equal(
        sdk.model_fingerprint,
        problem.metadata.model_fingerprint,
        "model fingerprint"
    );
    require_equal(
        sdk.provider_parameter_fingerprint,
        problem.metadata.provider_parameter_fingerprint,
        "Provider parameter fingerprint"
    );
    require_equal(
        sdk.expected_provider_topology_fingerprint,
        problem.metadata.provider_topology_fingerprint,
        "Provider topology fingerprint"
    );
    require_equal(
        sdk.artifact_identity,
        problem.metadata.artifact_identity,
        "evaluator artifact identity"
    );
    if (text(sdk.provider_sdk_capsule_name, "Provider SDK capsule name")
            != "epcsaft.native_sdk.v1"
        || sdk.provider_sdk_abi_version != 1u
        || sdk.provider_sdk_table_size == 0
        || sdk.provider_sdk_result_size == 0
        || sdk.provider_sdk_mixture_result_size == 0
        || sdk.provider_sdk_neutral_reference_result_size == 0
        || sdk.provider_sdk_neutral_reference_derivative_result_size == 0
        || sdk.provider_sdk_reacting_phase_parameter_result_size == 0) {
        throw std::invalid_argument(
            "evaluator Provider SDK artifact contract is incomplete"
        );
    }
    if ((sdk.single_thread_non_reentrant != 0
         && sdk.single_thread_non_reentrant != 1)
        || (sdk.value_only_avoids_derivative_work != 0
            && sdk.value_only_avoids_derivative_work != 1)) {
        throw std::invalid_argument(
            "evaluator threading or request-mode declaration is invalid"
        );
    }
    for (std::size_t index = 0; index < problem.rows.size(); ++index) {
        const auto& row = problem.rows[index];
        if (sdk.row_ids[index] != row.id
            || sdk.state_ids[index] != row.state_id
            || sdk.state_schema_ids[index] != row.state_schema_id
            || sdk.observation_source_ids[index] != row.source_id
            || sdk.source_locators[index] != row.source_locator
            || sdk.primitive_ids[index] != row.primitive_id
            || sdk.primitive_units[index] != row.primitive_unit
            || sdk.transform_ids[index] != row.transform
            || sdk.reference_ids[index] != row.reference_id
            || sdk.reference_fingerprints[index]
                != row.reference_fingerprint) {
            throw std::invalid_argument(
                row.id + ": evaluator row metadata does not match"
            );
        }
        if (row.transform != "identity"
            && row.transform != "natural_log") {
            throw std::invalid_argument(
                row.id + ": evaluator transform is unsupported"
            );
        }
        if (row.partition != "training"
            || !std::isfinite(row.observed) || row.observed <= 0.0
            || !std::isfinite(row.scale) || row.scale <= 0.0) {
            throw std::invalid_argument(
                row.id + ": partition, observed value, or scale is invalid"
            );
        }
    }
    for (std::size_t index = 0; index < problem.parameters.size(); ++index) {
        if (sdk.parameter_ids[index] != problem.parameters[index].id
            || sdk.parameter_units[index] != problem.parameters[index].unit) {
            throw std::invalid_argument(
                "evaluator parameter identity, order, or unit does not match"
            );
        }
    }
}

FitResult solve(
    const epcsaft_engine_regression_evaluator_sdk_v1& sdk,
    const Problem& problem
) {
    validate_contract(sdk, problem);
    const auto training = training_indices(problem);
    if (training.size() < problem.parameters.size()) {
        throw std::invalid_argument(
            "training residual count must be at least the parameter count"
        );
    }
    const internal::ProblemShape shape{
        problem.parameters.size(), 0, training.size()
    };
    std::vector<internal::CoordinateBound> bounds;
    bounds.reserve(problem.parameters.size());
    for (const auto& parameter : problem.parameters) {
        if (!std::isfinite(parameter.origin)
            || !std::isfinite(parameter.scale) || parameter.scale == 0.0
            || !std::isfinite(parameter.lower)
            || !std::isfinite(parameter.upper)
            || parameter.lower >= parameter.upper) {
            throw std::invalid_argument(
                "parameter transform or bounds are invalid"
            );
        }
        double lower = (parameter.lower - parameter.origin) / parameter.scale;
        double upper = (parameter.upper - parameter.origin) / parameter.scale;
        if (lower > upper) {
            std::swap(lower, upper);
        }
        bounds.push_back({lower, upper});
    }
    if (!std::isfinite(problem.maximum_condition_number)
        || problem.maximum_condition_number <= 0.0
        || problem.controls.maximum_iterations <= 0
        || !std::isfinite(problem.controls.maximum_solver_time_seconds)
        || problem.controls.maximum_solver_time_seconds <= 0.0
        || !std::isfinite(problem.controls.function_tolerance)
        || problem.controls.function_tolerance <= 0.0
        || !std::isfinite(problem.controls.gradient_tolerance)
        || problem.controls.gradient_tolerance <= 0.0
        || !std::isfinite(problem.controls.parameter_tolerance)
        || problem.controls.parameter_tolerance <= 0.0
        || !std::isfinite(problem.confirmation_parameter_delta)
        || problem.confirmation_parameter_delta <= 0.0
        || !std::isfinite(problem.confirmation_cost_delta)
        || problem.confirmation_cost_delta <= 0.0) {
        throw std::invalid_argument(
            "evaluator solver, conditioning, or confirmation controls are invalid"
        );
    }
    const auto exact = [&](
        const double* variables,
        std::size_t count,
        bool jacobian_requested,
        double* residuals,
        double* jacobian,
        std::string& failure_reason
    ) {
        return evaluate_training(
            sdk,
            problem,
            training,
            variables,
            count,
            jacobian_requested,
            residuals,
            jacobian,
            failure_reason
        );
    };

    FitResult fit;
    fit.solves.reserve(problem.starts.size());
    for (const auto& physical_start : problem.starts) {
        if (physical_start.size() != problem.parameters.size()) {
            throw std::invalid_argument(
                "every start must match the evaluator parameter count"
            );
        }
        std::vector<double> solver_start(problem.parameters.size());
        for (std::size_t index = 0; index < solver_start.size(); ++index) {
            solver_start[index] =
                (physical_start[index] - problem.parameters[index].origin)
                / problem.parameters[index].scale;
        }
        fit.solves.push_back(internal::solve(
            shape, solver_start, bounds, problem.controls, exact
        ));
    }
    if (fit.solves.size() < 2) {
        throw std::invalid_argument(
            "a primary and at least one confirmation start are required"
        );
    }
    const auto& primary = fit.solves.front();
    fit.confirmations_usable = true;
    for (std::size_t index = 1; index < fit.solves.size(); ++index) {
        const auto& confirmation = fit.solves[index];
        fit.confirmations_usable =
            fit.confirmations_usable
            && confirmation.summary.termination_type == ceres::CONVERGENCE
            && confirmation.summary.IsSolutionUsable()
            && confirmation.failure_reason.empty();
        for (std::size_t parameter = 0;
             parameter < problem.parameters.size();
             ++parameter) {
            fit.confirmation_parameter_delta = std::max(
                fit.confirmation_parameter_delta,
                std::abs(
                    primary.variables[parameter]
                    - confirmation.variables[parameter]
                )
            );
        }
        fit.confirmation_cost_delta = std::max(
            fit.confirmation_cost_delta,
            std::abs(
                primary.summary.final_cost
                - confirmation.summary.final_cost
            )
                / std::max({
                    std::abs(primary.summary.final_cost),
                    std::abs(confirmation.summary.final_cost),
                    problem.controls.function_tolerance,
                })
        );
    }
    const auto physical = to_physical(
        problem, primary.variables.data(), primary.variables.size()
    );
    const Batch final = evaluate(sdk, problem, physical, true);
    fit.rows.reserve(problem.rows.size());
    for (std::size_t row_index = 0; row_index < problem.rows.size();
         ++row_index) {
        const auto& spec = problem.rows[row_index];
        const double value = final.values[row_index];
        RowEvaluation row{};
        row.value = value;
        row.residual =
            spec.transform == "natural_log"
            ? (std::log(value) - std::log(spec.observed)) / spec.scale
            : (value - spec.observed) / spec.scale;
        row.scaled_solver_jacobian.resize(problem.parameters.size());
        for (std::size_t parameter = 0;
             parameter < problem.parameters.size();
             ++parameter) {
            const double provider_derivative = final.jacobian[
                checked_product(
                    row_index,
                    problem.parameters.size(),
                    "evaluator row offset"
                ) + parameter
            ];
            row.scaled_solver_jacobian[parameter] =
                provider_derivative
                * problem.parameters[parameter].scale
                / spec.scale
                / (spec.transform == "natural_log" ? value : 1.0);
        }
        row.certificate = final.rows[row_index];
        fit.rows.push_back(std::move(row));
    }
    fit.provider_parameter_fingerprint =
        problem.metadata.provider_parameter_fingerprint;
    return fit;
}

#ifdef EPCSAFT_ENGINE_REGRESSION_EVALUATOR_CORE_ONLY
bool evaluate_at_for_test(
    const epcsaft_engine_regression_evaluator_sdk_v1& sdk,
    const Problem& problem,
    const std::vector<double>& solver_parameters,
    std::vector<double>& residuals,
    std::vector<double>& jacobian,
    std::string& failure_reason
) {
    validate_contract(sdk, problem);
    const auto training = training_indices(problem);
    residuals.assign(training.size(), std::numeric_limits<double>::quiet_NaN());
    jacobian.assign(
        checked_product(
            training.size(),
            problem.parameters.size(),
            "evaluator test Jacobian"
        ),
        std::numeric_limits<double>::quiet_NaN()
    );
    return evaluate_training(
        sdk,
        problem,
        training,
        solver_parameters.data(),
        solver_parameters.size(),
        true,
        residuals.data(),
        jacobian.data(),
        failure_reason
    );
}
#endif

}  // namespace epcsaft_engine::regression::evaluator

#ifndef EPCSAFT_ENGINE_REGRESSION_EVALUATOR_CORE_ONLY
namespace epcsaft_engine::regression {
namespace {

struct PyOwned final {
    PyObject* value;

    explicit PyOwned(PyObject* object) : value(object) {}
    ~PyOwned() { Py_XDECREF(value); }
    PyOwned(const PyOwned&) = delete;
    PyOwned& operator=(const PyOwned&) = delete;
};

PyObject* fast_sequence(PyObject* value, const char* field) {
    PyObject* result = PySequence_Fast(value, field);
    if (result == nullptr) {
        throw std::invalid_argument(field);
    }
    return result;
}

std::string py_text(PyObject* value, const char* field) {
    if (!PyUnicode_Check(value)) {
        throw std::invalid_argument(std::string(field) + " must be text");
    }
    const char* text_value = PyUnicode_AsUTF8(value);
    if (text_value == nullptr || text_value[0] == '\0') {
        throw std::invalid_argument(std::string(field) + " must be nonempty");
    }
    if (strnlen(
            text_value,
            EPCSAFT_ENGINE_REGRESSION_EVALUATOR_V1_TEXT_CAPACITY
        ) >= EPCSAFT_ENGINE_REGRESSION_EVALUATOR_V1_TEXT_CAPACITY) {
        throw std::invalid_argument(std::string(field) + " is too long");
    }
    return text_value;
}

double py_number(PyObject* value, const char* field) {
    const double result = PyFloat_AsDouble(value);
    if (PyErr_Occurred() != nullptr || !std::isfinite(result)) {
        throw std::invalid_argument(std::string(field) + " must be finite");
    }
    return result;
}

long py_integer(PyObject* value, const char* field) {
    const long result = PyLong_AsLong(value);
    if (PyErr_Occurred() != nullptr) {
        throw std::invalid_argument(std::string(field) + " must be an integer");
    }
    return result;
}

std::vector<double> py_doubles(PyObject* value, const char* field) {
    PyOwned sequence(fast_sequence(value, field));
    std::vector<double> result;
    const Py_ssize_t count = PySequence_Fast_GET_SIZE(sequence.value);
    result.reserve(static_cast<std::size_t>(count));
    for (Py_ssize_t index = 0; index < count; ++index) {
        result.push_back(py_number(
            PySequence_Fast_GET_ITEM(sequence.value, index), field
        ));
    }
    return result;
}

evaluator::Problem parse_problem(PyObject* payload) {
    PyOwned root(fast_sequence(payload, "evaluator payload"));
    evaluator::Problem problem;
    if (PySequence_Fast_GET_SIZE(root.value) != 5) {
        throw std::invalid_argument(
            "evaluator payload must contain five contract blocks"
        );
    }
    PyOwned metadata(fast_sequence(
        PySequence_Fast_GET_ITEM(root.value, 0), "evaluator metadata"
    ));
    {
        if (PySequence_Fast_GET_SIZE(metadata.value) != 10) {
            throw std::invalid_argument(
                "evaluator metadata must contain ten identities"
            );
        }
        auto meta = [&](Py_ssize_t index, const char* field) {
            return py_text(PySequence_Fast_GET_ITEM(metadata.value, index), field);
        };
        problem.metadata = {
            meta(0, "evaluator identity"),
            meta(1, "capability identity"),
            meta(2, "capability fingerprint"),
            meta(3, "Provider artifact identity"),
            meta(4, "owner artifact identity"),
            meta(5, "contract fingerprint"),
            meta(6, "model fingerprint"),
            meta(7, "Provider parameter fingerprint"),
            meta(8, "Provider topology fingerprint"),
            meta(9, "artifact identity"),
        };
    }

    PyOwned parameters(fast_sequence(
        PySequence_Fast_GET_ITEM(root.value, 1), "evaluator parameters"
    ));
    {
        const Py_ssize_t parameter_count =
            PySequence_Fast_GET_SIZE(parameters.value);
        problem.parameters.reserve(
            static_cast<std::size_t>(parameter_count)
        );
        for (Py_ssize_t index = 0; index < parameter_count; ++index) {
            PyOwned item(fast_sequence(
                PySequence_Fast_GET_ITEM(parameters.value, index),
                "evaluator parameter"
            ));
            if (PySequence_Fast_GET_SIZE(item.value) != 6) {
                throw std::invalid_argument(
                    "evaluator parameter must contain six fields"
                );
            }
            problem.parameters.push_back({
                py_text(PySequence_Fast_GET_ITEM(item.value, 0), "parameter id"),
                py_text(PySequence_Fast_GET_ITEM(item.value, 1), "parameter unit"),
                py_number(PySequence_Fast_GET_ITEM(item.value, 2), "parameter origin"),
                py_number(PySequence_Fast_GET_ITEM(item.value, 3), "parameter scale"),
                py_number(PySequence_Fast_GET_ITEM(item.value, 4), "parameter lower"),
                py_number(PySequence_Fast_GET_ITEM(item.value, 5), "parameter upper"),
            });
        }
    }

    PyOwned starts(fast_sequence(
        PySequence_Fast_GET_ITEM(root.value, 2), "evaluator starts"
    ));
    {
        const Py_ssize_t start_count = PySequence_Fast_GET_SIZE(starts.value);
        problem.starts.reserve(static_cast<std::size_t>(start_count));
        for (Py_ssize_t index = 0; index < start_count; ++index) {
            problem.starts.push_back(py_doubles(
                PySequence_Fast_GET_ITEM(starts.value, index), "parameter start"
            ));
        }
    }

    PyOwned rows(fast_sequence(
        PySequence_Fast_GET_ITEM(root.value, 3), "evaluator rows"
    ));
    {
        const Py_ssize_t row_count = PySequence_Fast_GET_SIZE(rows.value);
        problem.rows.reserve(static_cast<std::size_t>(row_count));
        for (Py_ssize_t index = 0; index < row_count; ++index) {
            PyOwned item(fast_sequence(
                PySequence_Fast_GET_ITEM(rows.value, index), "evaluator row"
            ));
            if (PySequence_Fast_GET_SIZE(item.value) != 13) {
                throw std::invalid_argument(
                    "evaluator row must contain thirteen fields"
                );
            }
            problem.rows.push_back({
                py_text(PySequence_Fast_GET_ITEM(item.value, 0), "row id"),
                py_text(PySequence_Fast_GET_ITEM(item.value, 1), "partition"),
                py_text(PySequence_Fast_GET_ITEM(item.value, 2), "state id"),
                py_text(PySequence_Fast_GET_ITEM(item.value, 3), "state schema id"),
                py_text(PySequence_Fast_GET_ITEM(item.value, 4), "source id"),
                py_text(PySequence_Fast_GET_ITEM(item.value, 5), "source locator"),
                py_text(PySequence_Fast_GET_ITEM(item.value, 6), "primitive id"),
                py_text(PySequence_Fast_GET_ITEM(item.value, 7), "primitive unit"),
                py_text(PySequence_Fast_GET_ITEM(item.value, 8), "transform"),
                py_text(PySequence_Fast_GET_ITEM(item.value, 9), "reference id"),
                py_text(
                    PySequence_Fast_GET_ITEM(item.value, 10),
                    "reference fingerprint"
                ),
                py_number(PySequence_Fast_GET_ITEM(item.value, 11), "observed value"),
                py_number(PySequence_Fast_GET_ITEM(item.value, 12), "residual scale"),
            });
        }
    }

    PyOwned controls(fast_sequence(
        PySequence_Fast_GET_ITEM(root.value, 4), "evaluator controls"
    ));
    {
        if (PySequence_Fast_GET_SIZE(controls.value) != 8) {
            throw std::invalid_argument(
                "evaluator controls must contain eight values"
            );
        }
        problem.maximum_condition_number = py_number(
            PySequence_Fast_GET_ITEM(controls.value, 0),
            "maximum condition number"
        );
        const long iterations = py_integer(
            PySequence_Fast_GET_ITEM(controls.value, 1), "maximum iterations"
        );
        if (iterations <= 0
            || iterations > std::numeric_limits<int>::max()) {
            throw std::invalid_argument(
                "maximum iterations is outside the supported range"
            );
        }
        problem.controls = {
            static_cast<int>(iterations),
            py_number(
                PySequence_Fast_GET_ITEM(controls.value, 2),
                "maximum solver time"
            ),
            py_number(
                PySequence_Fast_GET_ITEM(controls.value, 3),
                "function tolerance"
            ),
            py_number(
                PySequence_Fast_GET_ITEM(controls.value, 4),
                "gradient tolerance"
            ),
            py_number(
                PySequence_Fast_GET_ITEM(controls.value, 5),
                "parameter tolerance"
            ),
        };
        problem.confirmation_parameter_delta = py_number(
            PySequence_Fast_GET_ITEM(controls.value, 6),
            "confirmation parameter delta"
        );
        problem.confirmation_cost_delta = py_number(
            PySequence_Fast_GET_ITEM(controls.value, 7),
            "confirmation cost delta"
        );
    }
    return problem;
}

const epcsaft_engine_regression_evaluator_sdk_v1* checked_sdk(
    PyObject* capsule
) {
    if (!PyCapsule_CheckExact(capsule)) {
        throw std::invalid_argument(
            "evaluator handle must be an exact CPython capsule"
        );
    }
    void* pointer = PyCapsule_GetPointer(
        capsule, EPCSAFT_ENGINE_REGRESSION_EVALUATOR_V1_CAPSULE_NAME
    );
    if (pointer == nullptr) {
        throw std::invalid_argument(
            "evaluator capsule name or pointer is invalid"
        );
    }
    struct Prefix final {
        std::uint32_t abi_version;
        std::size_t table_size;
    };
    Prefix prefix{};
    std::memcpy(&prefix, pointer, sizeof(prefix));
    if (prefix.abi_version
            != EPCSAFT_ENGINE_REGRESSION_EVALUATOR_V1_ABI_VERSION
        || prefix.table_size
            < sizeof(epcsaft_engine_regression_evaluator_sdk_v1)) {
        throw std::invalid_argument(
            "evaluator capsule ABI or table size is invalid"
        );
    }
    return static_cast<const epcsaft_engine_regression_evaluator_sdk_v1*>(
        pointer
    );
}

std::string termination_name(ceres::TerminationType type) {
    switch (type) {
        case ceres::CONVERGENCE:
            return "CONVERGENCE";
        case ceres::NO_CONVERGENCE:
            return "NO_CONVERGENCE";
        case ceres::FAILURE:
            return "FAILURE";
        case ceres::USER_SUCCESS:
            return "USER_SUCCESS";
        case ceres::USER_FAILURE:
            return "USER_FAILURE";
    }
    return "UNKNOWN";
}

PyObject* doubles_tuple(const std::vector<double>& values) {
    PyObject* result = PyTuple_New(
        static_cast<Py_ssize_t>(values.size())
    );
    if (result == nullptr) {
        return nullptr;
    }
    for (std::size_t index = 0; index < values.size(); ++index) {
        PyObject* item = PyFloat_FromDouble(values[index]);
        if (item == nullptr) {
            Py_DECREF(result);
            return nullptr;
        }
        PyTuple_SET_ITEM(result, static_cast<Py_ssize_t>(index), item);
    }
    return result;
}

PyObject* strings_tuple(const std::vector<std::string>& values) {
    PyObject* result = PyTuple_New(
        static_cast<Py_ssize_t>(values.size())
    );
    if (result == nullptr) {
        return nullptr;
    }
    for (std::size_t index = 0; index < values.size(); ++index) {
        PyObject* item = PyUnicode_FromString(values[index].c_str());
        if (item == nullptr) {
            Py_DECREF(result);
            return nullptr;
        }
        PyTuple_SET_ITEM(result, static_cast<Py_ssize_t>(index), item);
    }
    return result;
}

PyObject* rows_tuple(const evaluator::FitResult& fit) {
    PyObject* rows = PyTuple_New(
        static_cast<Py_ssize_t>(fit.rows.size())
    );
    if (rows == nullptr) {
        return nullptr;
    }
    for (std::size_t index = 0; index < fit.rows.size(); ++index) {
        const auto& row = fit.rows[index];
        PyObject* item = PyTuple_New(13);
        PyObject* derivatives = doubles_tuple(row.scaled_solver_jacobian);
        if (item == nullptr || derivatives == nullptr) {
            Py_XDECREF(item);
            Py_XDECREF(derivatives);
            Py_DECREF(rows);
            return nullptr;
        }
        PyTuple_SET_ITEM(item, 0, PyFloat_FromDouble(row.value));
        PyTuple_SET_ITEM(item, 1, PyFloat_FromDouble(row.residual));
        PyTuple_SET_ITEM(item, 2, derivatives);
        PyTuple_SET_ITEM(
            item, 3, PyUnicode_FromString(row.certificate.solver_status)
        );
        PyTuple_SET_ITEM(
            item, 4, PyUnicode_FromString(row.certificate.numerical_status)
        );
        PyTuple_SET_ITEM(
            item, 5, PyUnicode_FromString(row.certificate.physical_status)
        );
        PyTuple_SET_ITEM(
            item, 6, PyUnicode_FromString(row.certificate.derivative_status)
        );
        PyTuple_SET_ITEM(
            item, 7, PyUnicode_FromString(row.certificate.chart_topology)
        );
        PyTuple_SET_ITEM(
            item, 8,
            PyUnicode_FromString(
                row.certificate.provider_topology_fingerprint
            )
        );
        PyTuple_SET_ITEM(
            item, 9, PyLong_FromSize_t(row.certificate.kkt_dimension)
        );
        PyTuple_SET_ITEM(
            item, 10, PyLong_FromSize_t(row.certificate.kkt_rank)
        );
        PyTuple_SET_ITEM(
            item, 11,
            PyFloat_FromDouble(row.certificate.kkt_condition_number_inf)
        );
        PyTuple_SET_ITEM(
            item, 12, PyUnicode_FromString(row.certificate.reason)
        );
        for (Py_ssize_t field = 0; field < 13; ++field) {
            if (PyTuple_GET_ITEM(item, field) == nullptr) {
                Py_DECREF(item);
                Py_DECREF(rows);
                return nullptr;
            }
        }
        PyTuple_SET_ITEM(rows, static_cast<Py_ssize_t>(index), item);
    }
    return rows;
}

PyObject* capability_tuple(
    const epcsaft_engine_regression_evaluator_sdk_v1& sdk
) {
    return Py_BuildValue(
        "(sKkkkkkkOO)",
        sdk.provider_sdk_capsule_name,
        static_cast<unsigned long long>(sdk.provider_sdk_abi_version),
        static_cast<unsigned long>(sdk.provider_sdk_table_size),
        static_cast<unsigned long>(sdk.provider_sdk_result_size),
        static_cast<unsigned long>(sdk.provider_sdk_mixture_result_size),
        static_cast<unsigned long>(
            sdk.provider_sdk_neutral_reference_result_size
        ),
        static_cast<unsigned long>(
            sdk.provider_sdk_neutral_reference_derivative_result_size
        ),
        static_cast<unsigned long>(
            sdk.provider_sdk_reacting_phase_parameter_result_size
        ),
        sdk.single_thread_non_reentrant ? Py_True : Py_False,
        sdk.value_only_avoids_derivative_work ? Py_True : Py_False
    );
}

struct ThreadRelease final {
    ThreadRelease() : state(PyEval_SaveThread()) {}
    ~ThreadRelease() {
        PyEval_RestoreThread(state);
    }
    PyThreadState* state;
};

}  // namespace

PyObject* solve_evaluator_python(
    PyObject* capsule, PyObject* payload
) {
    try {
        const auto* sdk = checked_sdk(capsule);
        const evaluator::Problem problem = parse_problem(payload);
        evaluator::FitResult fit;
        if (sdk->single_thread_non_reentrant != 0) {
            fit = evaluator::solve(*sdk, problem);
        } else {
            ThreadRelease release;
            fit = evaluator::solve(*sdk, problem);
        }
        const auto& primary = fit.solves.front();
        std::vector<double> physical;
        std::vector<double> distances;
        std::vector<std::string> active_bounds;
        physical.reserve(problem.parameters.size());
        distances.reserve(problem.parameters.size());
        active_bounds.reserve(problem.parameters.size());
        for (std::size_t index = 0; index < problem.parameters.size(); ++index) {
            const auto& parameter = problem.parameters[index];
            const double value =
                parameter.origin + parameter.scale * primary.variables[index];
            physical.push_back(value);
            distances.push_back(std::min(
                value - parameter.lower, parameter.upper - value
            ));
            const double tolerance =
                std::sqrt(std::numeric_limits<double>::epsilon())
                * std::max(1.0, parameter.upper - parameter.lower);
            active_bounds.push_back(
                std::abs(value - parameter.lower) <= tolerance
                ? "lower"
                : std::abs(value - parameter.upper) <= tolerance
                ? "upper"
                : ""
            );
        }
        PyObject* result = PyTuple_New(28);
        PyObject* physical_tuple = doubles_tuple(physical);
        PyObject* distance_tuple = doubles_tuple(distances);
        PyObject* active_tuple = strings_tuple(active_bounds);
        PyObject* residuals = doubles_tuple(primary.residuals);
        PyObject* jacobian = doubles_tuple(primary.jacobian);
        PyObject* full_singular = doubles_tuple(
            primary.full_jacobian.singular_values
        );
        PyObject* projected_singular = doubles_tuple(
            primary.projected_parameter_jacobian.singular_values
        );
        PyObject* row_records = rows_tuple(fit);
        PyObject* capability = capability_tuple(*sdk);
        if (result == nullptr || physical_tuple == nullptr
            || distance_tuple == nullptr || active_tuple == nullptr
            || residuals == nullptr || jacobian == nullptr
            || full_singular == nullptr || projected_singular == nullptr
            || row_records == nullptr || capability == nullptr) {
            Py_XDECREF(result);
            Py_XDECREF(physical_tuple);
            Py_XDECREF(distance_tuple);
            Py_XDECREF(active_tuple);
            Py_XDECREF(residuals);
            Py_XDECREF(jacobian);
            Py_XDECREF(full_singular);
            Py_XDECREF(projected_singular);
            Py_XDECREF(row_records);
            Py_XDECREF(capability);
            return nullptr;
        }
        PyTuple_SET_ITEM(
            result, 0,
            PyUnicode_FromString(
                termination_name(primary.summary.termination_type).c_str()
            )
        );
        PyTuple_SET_ITEM(
            result, 1,
            Py_NewRef(primary.summary.IsSolutionUsable() ? Py_True : Py_False)
        );
        PyTuple_SET_ITEM(
            result, 2, PyFloat_FromDouble(primary.summary.initial_cost)
        );
        PyTuple_SET_ITEM(
            result, 3, PyFloat_FromDouble(primary.summary.final_cost)
        );
        PyTuple_SET_ITEM(
            result, 4, PyLong_FromSize_t(primary.summary.iterations.size())
        );
        PyTuple_SET_ITEM(result, 5, physical_tuple);
        PyTuple_SET_ITEM(result, 6, distance_tuple);
        PyTuple_SET_ITEM(result, 7, active_tuple);
        PyTuple_SET_ITEM(result, 8, residuals);
        PyTuple_SET_ITEM(result, 9, jacobian);
        PyTuple_SET_ITEM(result, 10, full_singular);
        PyTuple_SET_ITEM(
            result, 11, PyLong_FromLong(primary.full_jacobian.rank)
        );
        PyTuple_SET_ITEM(
            result, 12,
            PyFloat_FromDouble(primary.full_jacobian.condition_number)
        );
        PyTuple_SET_ITEM(result, 13, projected_singular);
        PyTuple_SET_ITEM(
            result, 14,
            PyLong_FromLong(
                primary.projected_parameter_jacobian.rank
            )
        );
        PyTuple_SET_ITEM(
            result, 15,
            PyFloat_FromDouble(
                primary.projected_parameter_jacobian.condition_number
            )
        );
        PyTuple_SET_ITEM(
            result, 16, PyLong_FromSize_t(fit.solves.size() - 1)
        );
        PyTuple_SET_ITEM(
            result, 17,
            PyFloat_FromDouble(fit.confirmation_parameter_delta)
        );
        PyTuple_SET_ITEM(
            result, 18, PyFloat_FromDouble(fit.confirmation_cost_delta)
        );
        PyTuple_SET_ITEM(
            result, 19,
            Py_NewRef(fit.confirmations_usable ? Py_True : Py_False)
        );
        PyTuple_SET_ITEM(result, 20, row_records);
        PyTuple_SET_ITEM(
            result, 21, PyUnicode_FromString(primary.failure_reason.c_str())
        );
        PyTuple_SET_ITEM(
            result, 22, PyLong_FromSize_t(problem.parameters.size())
        );
        PyTuple_SET_ITEM(
            result, 23,
            PyLong_FromSize_t(
                static_cast<std::size_t>(std::count_if(
                    problem.rows.cbegin(),
                    problem.rows.cend(),
                    [](const evaluator::Row& row) {
                        return row.partition == "training";
                    }
                ))
            )
        );
        PyTuple_SET_ITEM(
            result, 24,
            PyLong_FromLong(primary.summary.num_residual_evaluations)
        );
        PyTuple_SET_ITEM(
            result, 25,
            PyLong_FromLong(primary.summary.num_jacobian_evaluations)
        );
        PyTuple_SET_ITEM(result, 26, capability);
        PyTuple_SET_ITEM(
            result, 27,
            PyUnicode_FromString(
                fit.provider_parameter_fingerprint.c_str()
            )
        );
        for (Py_ssize_t index = 0; index < 28; ++index) {
            if (PyTuple_GET_ITEM(result, index) == nullptr) {
                Py_DECREF(result);
                return nullptr;
            }
        }
        return result;
    } catch (const std::exception& error) {
        if (PyErr_Occurred() != nullptr) {
            PyErr_Clear();
        }
        PyErr_SetString(PyExc_RuntimeError, error.what());
        return nullptr;
    } catch (...) {
        if (PyErr_Occurred() != nullptr) {
            PyErr_Clear();
        }
        PyErr_SetString(PyExc_RuntimeError, "unknown evaluator fit failure");
        return nullptr;
    }
}

}  // namespace epcsaft_engine::regression
#endif
