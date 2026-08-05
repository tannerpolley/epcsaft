#include "figiel_water_factor_fit.hpp"

#include <epcsaft/native_sdk_v1.h>

#include <ceres/ceres.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace epcsaft_engine::regression {
namespace {

constexpr std::size_t row_count = 21;

struct Observation final {
    std::string row_id;
    double molality;
    double observed;
};

struct Start final {
    std::string name;
    double value;
};

struct Payload final {
    std::vector<Observation> observations;
    std::string expected_fingerprint;
    std::vector<Start> starts;
    double temperature;
    double pressure;
    double lower_bound;
    double upper_bound;
    int max_iterations;
    double max_solver_time_seconds;
    double function_tolerance;
    double gradient_tolerance;
    double parameter_tolerance;
    double rank_multiplier;
};

struct Row final {
    std::string row_id;
    double molality;
    double observed;
    double log_modeled;
    double modeled;
    double residual;
    double provider_log_derivative;
    double residual_derivative;
    double reference_molality;
    double reference_convergence;
    double derivative_convergence;
    std::string fingerprint;
};

struct Evaluation final {
    std::vector<double> residuals;
    std::vector<double> jacobian;
    std::vector<Row> rows;
};

struct SolveOutcome final {
    std::string name;
    ceres::Solver::Summary summary;
    double parameter;
    Evaluation evaluation;
    double singular_value;
    double rank_threshold;
    int rank;
    bool complete_column;
    std::string failure_reason;
};

std::string text(PyObject* object, const char* label) {
    if (!PyUnicode_Check(object)) {
        throw std::invalid_argument(std::string(label) + " must be text");
    }
    const char* value = PyUnicode_AsUTF8(object);
    if (value == nullptr) {
        throw std::invalid_argument(std::string(label) + " must be UTF-8");
    }
    return value;
}

double number(PyObject* object, const char* label) {
    const double value = PyFloat_AsDouble(object);
    if (PyErr_Occurred() != nullptr) {
        PyErr_Clear();
        throw std::invalid_argument(std::string(label) + " must be numeric");
    }
    if (!std::isfinite(value)) {
        throw std::invalid_argument(std::string(label) + " must be finite");
    }
    return value;
}

std::size_t nonnegative_integer(PyObject* object, const char* label) {
    const unsigned long long value = PyLong_AsUnsignedLongLong(object);
    if (PyErr_Occurred() != nullptr) {
        PyErr_Clear();
        throw std::invalid_argument(
            std::string(label) + " must be a nonnegative integer"
        );
    }
    return static_cast<std::size_t>(value);
}

Payload parse_payload(PyObject* object) {
    PyObject* sequence = PySequence_Fast(
        object, "water-factor payload must be a sequence"
    );
    if (sequence == nullptr || PySequence_Fast_GET_SIZE(sequence) != 12) {
        Py_XDECREF(sequence);
        throw std::invalid_argument("water-factor payload has the wrong length");
    }

    Payload payload{};
    PyObject* observations = PySequence_Fast(
        PySequence_Fast_GET_ITEM(sequence, 0),
        "water-factor observations must be a sequence"
    );
    if (
        observations == nullptr
        || PySequence_Fast_GET_SIZE(observations)
            != static_cast<Py_ssize_t>(row_count)
    ) {
        Py_XDECREF(observations);
        Py_DECREF(sequence);
        throw std::invalid_argument(
            "water-factor contract requires 21 observations"
        );
    }
    payload.observations.reserve(row_count);
    for (std::size_t index = 0; index < row_count; ++index) {
        PyObject* row = PySequence_Fast(
            PySequence_Fast_GET_ITEM(
                observations, static_cast<Py_ssize_t>(index)
            ),
            "water-factor observation must be a sequence"
        );
        if (row == nullptr || PySequence_Fast_GET_SIZE(row) != 3) {
            Py_XDECREF(row);
            Py_DECREF(observations);
            Py_DECREF(sequence);
            throw std::invalid_argument(
                "water-factor observation has the wrong length"
            );
        }
        Observation observation{
            text(PySequence_Fast_GET_ITEM(row, 0), "row id"),
            number(PySequence_Fast_GET_ITEM(row, 1), "molality"),
            number(PySequence_Fast_GET_ITEM(row, 2), "observed MIAC"),
        };
        Py_DECREF(row);
        if (
            observation.row_id.empty()
            || observation.molality <= 0.0
            || observation.observed <= 0.0
        ) {
            Py_DECREF(observations);
            Py_DECREF(sequence);
            throw std::invalid_argument(
                "water-factor observation is outside the frozen domain"
            );
        }
        payload.observations.push_back(std::move(observation));
    }
    Py_DECREF(observations);

    payload.expected_fingerprint = text(
        PySequence_Fast_GET_ITEM(sequence, 1), "expected fingerprint"
    );
    PyObject* starts = PySequence_Fast(
        PySequence_Fast_GET_ITEM(sequence, 2), "starts must be a sequence"
    );
    if (starts == nullptr || PySequence_Fast_GET_SIZE(starts) != 2) {
        Py_XDECREF(starts);
        Py_DECREF(sequence);
        throw std::invalid_argument(
            "water-factor contract requires two starts"
        );
    }
    for (Py_ssize_t index = 0; index < 2; ++index) {
        PyObject* start = PySequence_Fast(
            PySequence_Fast_GET_ITEM(starts, index), "start must be a pair"
        );
        if (start == nullptr || PySequence_Fast_GET_SIZE(start) != 2) {
            Py_XDECREF(start);
            Py_DECREF(starts);
            Py_DECREF(sequence);
            throw std::invalid_argument("water-factor start must be a pair");
        }
        payload.starts.push_back(Start{
            text(PySequence_Fast_GET_ITEM(start, 0), "start name"),
            number(PySequence_Fast_GET_ITEM(start, 1), "start value"),
        });
        Py_DECREF(start);
    }
    Py_DECREF(starts);

    payload.temperature = number(
        PySequence_Fast_GET_ITEM(sequence, 3), "temperature"
    );
    payload.pressure = number(
        PySequence_Fast_GET_ITEM(sequence, 4), "pressure"
    );
    PyObject* bounds = PySequence_Fast(
        PySequence_Fast_GET_ITEM(sequence, 5), "bounds must be a pair"
    );
    if (bounds == nullptr || PySequence_Fast_GET_SIZE(bounds) != 2) {
        Py_XDECREF(bounds);
        Py_DECREF(sequence);
        throw std::invalid_argument("water-factor bounds must be a pair");
    }
    payload.lower_bound = number(
        PySequence_Fast_GET_ITEM(bounds, 0), "lower bound"
    );
    payload.upper_bound = number(
        PySequence_Fast_GET_ITEM(bounds, 1), "upper bound"
    );
    Py_DECREF(bounds);
    payload.max_iterations = static_cast<int>(
        nonnegative_integer(
            PySequence_Fast_GET_ITEM(sequence, 6), "maximum iterations"
        )
    );
    payload.max_solver_time_seconds = number(
        PySequence_Fast_GET_ITEM(sequence, 7), "maximum solver time"
    );
    payload.function_tolerance = number(
        PySequence_Fast_GET_ITEM(sequence, 8), "function tolerance"
    );
    payload.gradient_tolerance = number(
        PySequence_Fast_GET_ITEM(sequence, 9), "gradient tolerance"
    );
    payload.parameter_tolerance = number(
        PySequence_Fast_GET_ITEM(sequence, 10), "parameter tolerance"
    );
    payload.rank_multiplier = number(
        PySequence_Fast_GET_ITEM(sequence, 11), "rank multiplier"
    );
    Py_DECREF(sequence);

    if (
        payload.expected_fingerprint.rfind("sha256:", 0) != 0
        || payload.temperature != 298.15
        || payload.pressure != 100000.0
        || payload.lower_bound != 1.0
        || payload.upper_bound != 2.0
        || payload.starts[0].name != "primary"
        || payload.starts[0].value != 1.2
        || payload.starts[1].name != "upper"
        || payload.starts[1].value != 1.8
        || payload.max_iterations != 500
        || payload.max_solver_time_seconds != 180.0
        || payload.function_tolerance != 1.0e-10
        || payload.gradient_tolerance != 1.0e-10
        || payload.parameter_tolerance != 1.0e-10
        || payload.rank_multiplier != 100.0
    ) {
        throw std::invalid_argument(
            "water-factor numerical contract does not match the frozen design"
        );
    }
    return payload;
}

const epcsaft_native_sdk_v1* checked_table(
    PyObject* capsule, const Payload& payload
) {
    if (!PyCapsule_CheckExact(capsule)) {
        throw std::invalid_argument(
            "provider transport must be an exact CPython capsule"
        );
    }
    void* pointer = PyCapsule_GetPointer(
        capsule, EPCSAFT_NATIVE_SDK_V1_CAPSULE_NAME
    );
    if (pointer == nullptr) {
        throw std::runtime_error("provider capability unavailable");
    }
    struct Prefix final {
        std::uint32_t abi_version;
        std::size_t table_size;
    } prefix{};
    std::memcpy(&prefix, pointer, sizeof(prefix));
    constexpr std::size_t minimum_size =
        offsetof(
            epcsaft_native_sdk_v1,
            evaluate_aqueous_miac_solvation_factor_batch
        )
        + sizeof(epcsaft_evaluate_aqueous_miac_solvation_factor_batch_v1);
    if (
        prefix.abi_version != EPCSAFT_NATIVE_SDK_V1_ABI_VERSION
        || prefix.table_size < minimum_size
    ) {
        throw std::runtime_error("provider capability unavailable");
    }
    const auto* table = static_cast<const epcsaft_native_sdk_v1*>(pointer);
    if (
        table->component_count != 3
        || table->component_ids == nullptr
        || table->component_charges == nullptr
        || table->component_ids[0] == nullptr
        || table->component_ids[1] == nullptr
        || table->component_ids[2] == nullptr
        || std::string(table->component_ids[0]) != "water"
        || std::string(table->component_ids[1]) != "sodium-cation"
        || std::string(table->component_ids[2]) != "bromide-anion"
        || table->component_charges[0] != 0
        || table->component_charges[1] != 1
        || table->component_charges[2] != -1
        || table->aqueous_miac_solvation_factor_result_size
            != sizeof(epcsaft_aqueous_miac_solvation_factor_result_v1)
        || table->evaluate_aqueous_miac_solvation_factor_batch == nullptr
        || table->model_context == nullptr
        || payload.expected_fingerprint.empty()
    ) {
        throw std::runtime_error("provider NaBr water-factor capability unavailable");
    }
    return table;
}

Evaluation evaluate(
    const epcsaft_native_sdk_v1* table,
    const Payload& payload,
    double parameter
) {
    std::vector<double> molalities;
    molalities.reserve(row_count);
    for (const auto& observation : payload.observations) {
        molalities.push_back(observation.molality);
    }
    std::vector<epcsaft_aqueous_miac_solvation_factor_result_v1> results(
        row_count
    );
    for (auto& result : results) {
        result.struct_size = sizeof(result);
    }
    const int status = table->evaluate_aqueous_miac_solvation_factor_batch(
        table->model_context,
        payload.expected_fingerprint.c_str(),
        payload.temperature,
        payload.pressure,
        molalities.data(),
        molalities.size(),
        parameter,
        results.data(),
        results.size()
    );
    if (status != EPCSAFT_NATIVE_STATUS_OK_V1) {
        throw std::runtime_error(
            std::string("provider water-factor batch failed: ")
            + results.front().error
        );
    }

    Evaluation evaluation{
        std::vector<double>(row_count),
        std::vector<double>(row_count),
        std::vector<Row>(),
    };
    evaluation.rows.reserve(row_count);
    for (std::size_t index = 0; index < row_count; ++index) {
        const auto& result = results[index];
        const auto& observation = payload.observations[index];
        const std::string fingerprint(
            result.parameter_fingerprint,
            strnlen(
                result.parameter_fingerprint,
                EPCSAFT_NATIVE_SDK_V1_FINGERPRINT_SIZE
            )
        );
        if (
            result.status != status
            || fingerprint != payload.expected_fingerprint
            || !std::isfinite(
                result.log_mean_ionic_activity_coefficient_molality
            )
            || !std::isfinite(result.derivative)
            || !std::isfinite(result.reference_molality_mol_per_kg)
            || !std::isfinite(result.reference_convergence_error)
            || !std::isfinite(result.reference_derivative_convergence_error)
        ) {
            throw std::runtime_error(
                std::string("provider water-factor row failed: ")
                + result.error
            );
        }
        const double modeled = std::exp(
            result.log_mean_ionic_activity_coefficient_molality
        );
        const double ratio = modeled / observation.observed;
        const double residual = 1.0 - ratio;
        const double derivative = -ratio * result.derivative;
        evaluation.residuals[index] = residual;
        evaluation.jacobian[index] = derivative;
        evaluation.rows.push_back(Row{
            observation.row_id,
            observation.molality,
            observation.observed,
            result.log_mean_ionic_activity_coefficient_molality,
            modeled,
            residual,
            result.derivative,
            derivative,
            result.reference_molality_mol_per_kg,
            result.reference_convergence_error,
            result.reference_derivative_convergence_error,
            fingerprint,
        });
    }
    return evaluation;
}

class WaterFactorCost final : public ceres::CostFunction {
public:
    WaterFactorCost(
        const epcsaft_native_sdk_v1* table, const Payload& payload
    ) : table_(table), payload_(payload) {
        set_num_residuals(static_cast<int>(row_count));
        mutable_parameter_block_sizes()->push_back(1);
    }

    bool Evaluate(
        double const* const* values, double* residuals, double** jacobians
    ) const override {
        try {
            if (!cached_parameter_.has_value() || values[0][0] != *cached_parameter_) {
                cached_evaluation_ = evaluate(table_, payload_, values[0][0]);
                cached_parameter_ = values[0][0];
            }
            const Evaluation& result = *cached_evaluation_;
            std::copy(
                result.residuals.begin(), result.residuals.end(), residuals
            );
            if (jacobians != nullptr && jacobians[0] != nullptr) {
                std::copy(
                    result.jacobian.begin(),
                    result.jacobian.end(),
                    jacobians[0]
                );
            }
            failure_reason_.clear();
            return true;
        } catch (const std::exception& error) {
            failure_reason_ = error.what();
            return false;
        }
    }

    const std::string& failure_reason() const noexcept {
        return failure_reason_;
    }

private:
    const epcsaft_native_sdk_v1* table_;
    const Payload& payload_;
    mutable std::string failure_reason_;
    mutable std::optional<double> cached_parameter_;
    mutable std::optional<Evaluation> cached_evaluation_;
};

std::string termination_name(ceres::TerminationType termination) {
    switch (termination) {
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

SolveOutcome solve_one(
    const epcsaft_native_sdk_v1* table,
    const Payload& payload,
    const Start& start
) {
    SolveOutcome outcome{};
    outcome.name = start.name;
    outcome.parameter = start.value;
    WaterFactorCost cost(table, payload);
    ceres::Problem::Options problem_options;
    problem_options.cost_function_ownership = ceres::DO_NOT_TAKE_OWNERSHIP;
    ceres::Problem problem(problem_options);
    problem.AddResidualBlock(&cost, nullptr, &outcome.parameter);
    problem.SetParameterLowerBound(&outcome.parameter, 0, payload.lower_bound);
    problem.SetParameterUpperBound(&outcome.parameter, 0, payload.upper_bound);

    ceres::Solver::Options options;
    options.linear_solver_type = ceres::DENSE_QR;
    options.max_num_iterations = payload.max_iterations;
    options.max_solver_time_in_seconds = payload.max_solver_time_seconds;
    options.function_tolerance = payload.function_tolerance;
    options.gradient_tolerance = payload.gradient_tolerance;
    options.parameter_tolerance = payload.parameter_tolerance;
    options.logging_type = ceres::SILENT;
    options.num_threads = 1;
    ceres::Solve(options, &problem, &outcome.summary);

    outcome.failure_reason = cost.failure_reason();
    outcome.evaluation = evaluate(table, payload, outcome.parameter);
    double squared_norm = 0.0;
    outcome.complete_column = false;
    for (const double derivative : outcome.evaluation.jacobian) {
        if (!std::isfinite(derivative)) {
            outcome.complete_column = false;
            squared_norm = std::numeric_limits<double>::quiet_NaN();
            break;
        }
        outcome.complete_column = outcome.complete_column || derivative != 0.0;
        squared_norm += derivative * derivative;
    }
    outcome.singular_value = std::sqrt(squared_norm);
    outcome.rank_threshold =
        outcome.singular_value
        * static_cast<double>(row_count)
        * std::numeric_limits<double>::epsilon()
        * payload.rank_multiplier;
    outcome.rank =
        std::isfinite(outcome.singular_value)
        && outcome.singular_value > outcome.rank_threshold
        ? 1
        : 0;
    return outcome;
}

PyObject* rows_to_tuple(const std::vector<Row>& rows) {
    PyObject* tuple = PyTuple_New(static_cast<Py_ssize_t>(rows.size()));
    if (tuple == nullptr) {
        return nullptr;
    }
    for (std::size_t index = 0; index < rows.size(); ++index) {
        const Row& row = rows[index];
        PyObject* item = Py_BuildValue(
            "(sdddddddddds)",
            row.row_id.c_str(),
            row.molality,
            row.observed,
            row.log_modeled,
            row.modeled,
            row.residual,
            row.provider_log_derivative,
            row.residual_derivative,
            row.reference_molality,
            row.reference_convergence,
            row.derivative_convergence,
            row.fingerprint.c_str()
        );
        if (item == nullptr) {
            Py_DECREF(tuple);
            return nullptr;
        }
        PyTuple_SET_ITEM(tuple, static_cast<Py_ssize_t>(index), item);
    }
    return tuple;
}

PyObject* outcome_to_tuple(const SolveOutcome& outcome) {
    PyObject* rows = rows_to_tuple(outcome.evaluation.rows);
    if (rows == nullptr) {
        return nullptr;
    }
    PyObject* result = PyTuple_New(13);
    if (result == nullptr) {
        Py_DECREF(rows);
        return nullptr;
    }
    PyTuple_SET_ITEM(
        result, 0, PyUnicode_FromString(outcome.name.c_str())
    );
    PyTuple_SET_ITEM(
        result,
        1,
        PyUnicode_FromString(
            termination_name(outcome.summary.termination_type).c_str()
        )
    );
    PyTuple_SET_ITEM(
        result,
        2,
        Py_NewRef(outcome.summary.IsSolutionUsable() ? Py_True : Py_False)
    );
    PyTuple_SET_ITEM(
        result, 3, PyFloat_FromDouble(outcome.summary.initial_cost)
    );
    PyTuple_SET_ITEM(
        result, 4, PyFloat_FromDouble(outcome.summary.final_cost)
    );
    PyTuple_SET_ITEM(
        result,
        5,
        PyLong_FromSize_t(outcome.summary.iterations.size())
    );
    PyTuple_SET_ITEM(result, 6, PyFloat_FromDouble(outcome.parameter));
    PyTuple_SET_ITEM(result, 7, rows);
    PyTuple_SET_ITEM(result, 8, PyFloat_FromDouble(outcome.singular_value));
    PyTuple_SET_ITEM(result, 9, PyFloat_FromDouble(outcome.rank_threshold));
    PyTuple_SET_ITEM(result, 10, PyLong_FromLong(outcome.rank));
    PyTuple_SET_ITEM(
        result, 11, Py_NewRef(outcome.complete_column ? Py_True : Py_False)
    );
    PyTuple_SET_ITEM(
        result,
        12,
        PyUnicode_FromString(outcome.failure_reason.c_str())
    );
    return result;
}

}  // namespace

PyObject* solve_figiel_water_factor_python(
    PyObject* capsule, PyObject* payload_object
) {
    try {
        const Payload payload = parse_payload(payload_object);
        const epcsaft_native_sdk_v1* table = checked_table(capsule, payload);
        PyObject* starts = PyTuple_New(
            static_cast<Py_ssize_t>(payload.starts.size())
        );
        if (starts == nullptr) {
            return nullptr;
        }
        for (std::size_t index = 0; index < payload.starts.size(); ++index) {
            PyObject* item = outcome_to_tuple(
                solve_one(table, payload, payload.starts[index])
            );
            if (item == nullptr) {
                Py_DECREF(starts);
                return nullptr;
            }
            PyTuple_SET_ITEM(
                starts, static_cast<Py_ssize_t>(index), item
            );
        }
        return starts;
    } catch (const std::exception& error) {
        if (PyErr_Occurred() != nullptr) {
            PyErr_Clear();
        }
        PyErr_SetString(PyExc_RuntimeError, error.what());
        return nullptr;
    }
}

PyObject* evaluate_figiel_water_factor_python(
    PyObject* capsule, PyObject* payload_object, double parameter
) {
    try {
        const Payload payload = parse_payload(payload_object);
        const epcsaft_native_sdk_v1* table = checked_table(capsule, payload);
        return rows_to_tuple(evaluate(table, payload, parameter).rows);
    } catch (const std::exception& error) {
        if (PyErr_Occurred() != nullptr) {
            PyErr_Clear();
        }
        PyErr_SetString(PyExc_RuntimeError, error.what());
        return nullptr;
    }
}

}  // namespace epcsaft_engine::regression
