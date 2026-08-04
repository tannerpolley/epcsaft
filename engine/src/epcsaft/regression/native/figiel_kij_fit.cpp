#include "figiel_kij_fit.hpp"

#include <epcsaft/native_sdk_v1.h>

#include <ceres/ceres.h>
#include <Eigen/Dense>
#include <Eigen/SVD>

#include <algorithm>
#include <array>
#include <chrono>
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

constexpr std::size_t salt_count = 6;
constexpr std::size_t parameter_count = 11;
constexpr std::size_t row_count = 164;

std::uint64_t deadline_after(double seconds) {
    const auto now = std::chrono::steady_clock::now().time_since_epoch();
    const auto delta = std::chrono::duration_cast<std::chrono::nanoseconds>(
        std::chrono::duration<double>(seconds)
    );
    return static_cast<std::uint64_t>(
        std::chrono::duration_cast<std::chrono::nanoseconds>(now + delta)
            .count()
    );
}

struct Observation final {
    std::string row_id;
    std::string salt;
    std::size_t salt_index;
    double molality;
    double observed;
    std::array<std::size_t, 3> columns;
};

struct Schedule final {
    std::string name;
    std::array<double, parameter_count> values;
    bool reverse;
};

struct Payload final {
    std::vector<Observation> observations;
    std::array<std::string, salt_count> expected_fingerprints;
    std::vector<Schedule> schedules;
    std::array<double, 2> bounds;
    std::array<double, parameter_count> published_parameters;
    double temperature;
    double pressure;
    int max_iterations;
    double max_solver_time_seconds;
    double function_tolerance;
    double gradient_tolerance;
    double parameter_tolerance;
    double rank_multiplier;
};

struct Row final {
    std::string row_id;
    std::string salt;
    double molality;
    double observed;
    double log_modeled;
    double modeled;
    double residual;
    std::array<double, 3> local_derivative;
    double reference_molality;
    double reference_convergence;
    double derivative_convergence;
    std::string fingerprint;
};

struct Evaluation final {
    std::array<double, row_count> residuals{};
    std::array<double, row_count * parameter_count> jacobian{};
    std::array<Row, row_count> rows{};
};

struct CoordinateSolve final {
    std::size_t coordinate;
    double parameter;
    double initial_cost;
    double final_cost;
    std::size_t iterations;
    std::string termination;
    bool solution_usable;
    std::string failure_reason;
};

struct SolveOutcome final {
    std::string name;
    std::string coordinate_order;
    std::string termination{"NO_CONVERGENCE"};
    bool solution_usable{false};
    double initial_cost{std::numeric_limits<double>::infinity()};
    double final_cost{std::numeric_limits<double>::infinity()};
    std::size_t iterations{0};
    std::array<double, parameter_count> parameters{};
    std::array<CoordinateSolve, parameter_count> coordinate_solves{};
    Evaluation evaluation;
    std::array<double, parameter_count> singular_values{};
    double rank_threshold{0.0};
    int rank{0};
    double condition_number{std::numeric_limits<double>::infinity()};
    std::array<double, parameter_count> least_sensitive_direction{};
    std::array<bool, parameter_count> complete_columns{};
    std::string failure_reason;
};

std::string text(PyObject* object, const char* label) {
    if (!PyUnicode_Check(object)) {
        throw std::invalid_argument(std::string(label) + " must be text");
    }
    Py_ssize_t size = 0;
    const char* value = PyUnicode_AsUTF8AndSize(object, &size);
    if (value == nullptr) {
        throw std::invalid_argument(std::string(label) + " must be UTF-8");
    }
    return std::string(value, static_cast<std::size_t>(size));
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

std::size_t index_value(PyObject* object, const char* label) {
    const unsigned long long value = PyLong_AsUnsignedLongLong(object);
    if (PyErr_Occurred() != nullptr) {
        PyErr_Clear();
        throw std::invalid_argument(
            std::string(label) + " must be a nonnegative integer"
        );
    }
    return static_cast<std::size_t>(value);
}

template <std::size_t size>
std::array<double, size> doubles(PyObject* object, const char* label) {
    PyObject* sequence = PySequence_Fast(object, label);
    if (
        sequence == nullptr
        || PySequence_Fast_GET_SIZE(sequence) != static_cast<Py_ssize_t>(size)
    ) {
        Py_XDECREF(sequence);
        throw std::invalid_argument(std::string(label) + " has the wrong length");
    }
    std::array<double, size> result{};
    for (std::size_t index = 0; index < size; ++index) {
        result[index] = number(
            PySequence_Fast_GET_ITEM(
                sequence, static_cast<Py_ssize_t>(index)
            ),
            label
        );
    }
    Py_DECREF(sequence);
    return result;
}

Payload parse_payload(PyObject* object) {
    PyObject* sequence = PySequence_Fast(
        object, "Figiel interaction payload must be a sequence"
    );
    if (sequence == nullptr || PySequence_Fast_GET_SIZE(sequence) != 14) {
        Py_XDECREF(sequence);
        throw std::invalid_argument(
            "Figiel interaction payload has the wrong length"
        );
    }
    Payload payload{};
    PyObject* observations = PySequence_Fast(
        PySequence_Fast_GET_ITEM(sequence, 0),
        "observations must be a sequence"
    );
    if (
        observations == nullptr
        || PySequence_Fast_GET_SIZE(observations)
            != static_cast<Py_ssize_t>(row_count)
    ) {
        Py_XDECREF(observations);
        Py_DECREF(sequence);
        throw std::invalid_argument(
            "Figiel interaction fit requires exactly 164 observations"
        );
    }
    payload.observations.reserve(row_count);
    std::array<std::size_t, salt_count> salt_rows{};
    for (std::size_t row_index = 0; row_index < row_count; ++row_index) {
        PyObject* row = PySequence_Fast(
            PySequence_Fast_GET_ITEM(
                observations, static_cast<Py_ssize_t>(row_index)
            ),
            "observation must be a sequence"
        );
        if (row == nullptr || PySequence_Fast_GET_SIZE(row) != 6) {
            Py_XDECREF(row);
            Py_DECREF(observations);
            Py_DECREF(sequence);
            throw std::invalid_argument("observation has the wrong length");
        }
        PyObject* columns_object = PySequence_Fast(
            PySequence_Fast_GET_ITEM(row, 5),
            "column map must be a sequence"
        );
        if (
            columns_object == nullptr
            || PySequence_Fast_GET_SIZE(columns_object) != 3
        ) {
            Py_XDECREF(columns_object);
            Py_DECREF(row);
            Py_DECREF(observations);
            Py_DECREF(sequence);
            throw std::invalid_argument("column map must contain three indices");
        }
        std::array<std::size_t, 3> columns{};
        for (std::size_t local = 0; local < columns.size(); ++local) {
            columns[local] = index_value(
                PySequence_Fast_GET_ITEM(
                    columns_object, static_cast<Py_ssize_t>(local)
                ),
                "column index"
            );
            if (columns[local] >= parameter_count) {
                Py_DECREF(columns_object);
                Py_DECREF(row);
                Py_DECREF(observations);
                Py_DECREF(sequence);
                throw std::invalid_argument(
                    "column index is outside the 11-parameter block"
                );
            }
        }
        Py_DECREF(columns_object);
        Observation observed{
            text(PySequence_Fast_GET_ITEM(row, 0), "row id"),
            text(PySequence_Fast_GET_ITEM(row, 1), "salt"),
            index_value(PySequence_Fast_GET_ITEM(row, 2), "salt index"),
            number(PySequence_Fast_GET_ITEM(row, 3), "molality"),
            number(PySequence_Fast_GET_ITEM(row, 4), "observed MIAC"),
            columns,
        };
        Py_DECREF(row);
        if (
            observed.salt_index >= salt_count || observed.molality <= 0.0
            || observed.observed <= 0.0
        ) {
            Py_DECREF(observations);
            Py_DECREF(sequence);
            throw std::invalid_argument(
                "observation is outside the frozen aqueous domain"
            );
        }
        ++salt_rows[observed.salt_index];
        payload.observations.push_back(std::move(observed));
    }
    Py_DECREF(observations);
    if (salt_rows != std::array<std::size_t, salt_count>{29, 29, 28, 29, 21, 28}) {
        Py_DECREF(sequence);
        throw std::invalid_argument("salt row counts do not match the frozen dataset");
    }
    PyObject* fingerprints = PySequence_Fast(
        PySequence_Fast_GET_ITEM(sequence, 1),
        "fingerprints must be a sequence"
    );
    if (
        fingerprints == nullptr
        || PySequence_Fast_GET_SIZE(fingerprints)
            != static_cast<Py_ssize_t>(salt_count)
    ) {
        Py_XDECREF(fingerprints);
        Py_DECREF(sequence);
        throw std::invalid_argument("fingerprints must contain six entries");
    }
    for (std::size_t index = 0; index < salt_count; ++index) {
        payload.expected_fingerprints[index] = text(
            PySequence_Fast_GET_ITEM(
                fingerprints, static_cast<Py_ssize_t>(index)
            ),
            "fingerprint"
        );
    }
    Py_DECREF(fingerprints);
    PyObject* schedules = PySequence_Fast(
        PySequence_Fast_GET_ITEM(sequence, 2),
        "coordinate schedules must be a sequence"
    );
    if (
        schedules == nullptr || PySequence_Fast_GET_SIZE(schedules) != 2
    ) {
        Py_XDECREF(schedules);
        Py_DECREF(sequence);
        throw std::invalid_argument(
            "conditional recovery requires two declared starts"
        );
    }
    for (
        Py_ssize_t index = 0;
        index < PySequence_Fast_GET_SIZE(schedules);
        ++index
    ) {
        PyObject* item = PySequence_Fast(
            PySequence_Fast_GET_ITEM(schedules, index),
            "coordinate schedule must be a triple"
        );
        if (item == nullptr || PySequence_Fast_GET_SIZE(item) != 3) {
            Py_XDECREF(item);
            Py_DECREF(schedules);
            Py_DECREF(sequence);
            throw std::invalid_argument(
                "coordinate schedule must contain name, values, and order"
            );
        }
        const std::string order = text(
            PySequence_Fast_GET_ITEM(item, 2), "coordinate order"
        );
        if (order != "forward" && order != "reverse") {
            Py_DECREF(item);
            Py_DECREF(schedules);
            Py_DECREF(sequence);
            throw std::invalid_argument(
                "coordinate order must be forward or reverse"
            );
        }
        payload.schedules.push_back(Schedule{
            text(PySequence_Fast_GET_ITEM(item, 0), "schedule name"),
            doubles<parameter_count>(
                PySequence_Fast_GET_ITEM(item, 1), "schedule values"
            ),
            order == "reverse",
        });
        Py_DECREF(item);
    }
    Py_DECREF(schedules);
    payload.temperature = number(
        PySequence_Fast_GET_ITEM(sequence, 3), "temperature"
    );
    payload.pressure = number(
        PySequence_Fast_GET_ITEM(sequence, 4), "pressure"
    );
    payload.bounds = doubles<2>(
        PySequence_Fast_GET_ITEM(sequence, 5), "bounds"
    );
    payload.published_parameters = doubles<parameter_count>(
        PySequence_Fast_GET_ITEM(sequence, 6), "published parameters"
    );
    payload.max_iterations = static_cast<int>(
        index_value(PySequence_Fast_GET_ITEM(sequence, 7), "maximum iterations")
    );
    payload.max_solver_time_seconds = number(
        PySequence_Fast_GET_ITEM(sequence, 8), "maximum solver time"
    );
    payload.function_tolerance = number(
        PySequence_Fast_GET_ITEM(sequence, 9), "function tolerance"
    );
    payload.gradient_tolerance = number(
        PySequence_Fast_GET_ITEM(sequence, 10), "gradient tolerance"
    );
    payload.parameter_tolerance = number(
        PySequence_Fast_GET_ITEM(sequence, 11), "parameter tolerance"
    );
    payload.rank_multiplier = number(
        PySequence_Fast_GET_ITEM(sequence, 12), "rank multiplier"
    );
    const std::string specification_id = text(
        PySequence_Fast_GET_ITEM(sequence, 13),
        "specification identity"
    );
    Py_DECREF(sequence);
    constexpr std::array<double, parameter_count> published_parameters{
        -0.4, -0.3, -0.1, -0.3, -0.3, 0.8, 0.8, 0.0, 0.5, 0.65, -0.35
    };
    if (
        payload.temperature != 298.15 || payload.pressure != 100000.0
        || payload.bounds != std::array<double, 2>{-1.0, 1.0}
        || payload.published_parameters != published_parameters
        || payload.schedules[0].name != "primary"
        || payload.schedules[0].values
            != std::array<double, parameter_count>{}
        || payload.schedules[0].reverse
        || payload.schedules[1].name != "confirmation"
        || payload.schedules[1].values
            != std::array<double, parameter_count>{
                0.25, 0.25, 0.25, 0.25, 0.25, 0.25,
                0.25, 0.25, 0.25, 0.25, 0.25
            }
        || !payload.schedules[1].reverse
        || payload.max_iterations != 50
        || payload.max_solver_time_seconds != 180.0
        || payload.function_tolerance != 1.0e-10
        || payload.gradient_tolerance != 1.0e-10
        || payload.parameter_tolerance != 1.0e-10
        || payload.rank_multiplier != 100.0
        || specification_id != "figiel-2025-aqueous-kij-v1"
    ) {
        throw std::invalid_argument(
            "Figiel interaction numerical contract does not match the design"
        );
    }
    return payload;
}

const epcsaft_native_sdk_v1* checked_table(
    PyObject* capsule,
    std::size_t salt_index
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
        throw std::runtime_error("provider aqueous-kij capability unavailable");
    }
    struct Prefix final {
        std::uint32_t abi_version;
        std::size_t table_size;
    } prefix{};
    std::memcpy(&prefix, pointer, sizeof(prefix));
    constexpr std::size_t minimum_size =
        offsetof(
            epcsaft_native_sdk_v1,
            evaluate_aqueous_miac_kij_batch_bounded
        )
        + sizeof(epcsaft_evaluate_aqueous_miac_kij_batch_bounded_v1);
    if (
        prefix.abi_version != EPCSAFT_NATIVE_SDK_V1_ABI_VERSION
        || prefix.table_size < minimum_size
    ) {
        throw std::runtime_error("provider aqueous-kij capability unavailable");
    }
    const auto* table = static_cast<const epcsaft_native_sdk_v1*>(pointer);
    constexpr std::array<const char*, salt_count> cations{
        "lithium-cation",
        "sodium-cation",
        "potassium-cation",
        "lithium-cation",
        "sodium-cation",
        "potassium-cation",
    };
    constexpr std::array<const char*, salt_count> anions{
        "chloride-anion",
        "chloride-anion",
        "chloride-anion",
        "bromide-anion",
        "bromide-anion",
        "bromide-anion",
    };
    if (
        table->component_count != 3 || table->component_ids == nullptr
        || table->component_charges == nullptr
        || table->component_ids[0] == nullptr
        || table->component_ids[1] == nullptr
        || table->component_ids[2] == nullptr
        || std::string(table->component_ids[0]) != "water"
        || std::string(table->component_ids[1]) != cations[salt_index]
        || std::string(table->component_ids[2]) != anions[salt_index]
        || table->component_charges[0] != 0
        || table->component_charges[1] != 1
        || table->component_charges[2] != -1
    ) {
        throw std::invalid_argument(
            "provider component order does not match the frozen aqueous salt"
        );
    }
    if (
        table->aqueous_miac_kij_result_size
            != sizeof(epcsaft_aqueous_miac_kij_result_v1)
        || table->evaluation_budget_size
            != sizeof(epcsaft_native_evaluation_budget_v1)
        || table->model_context == nullptr
        || table->evaluate_aqueous_miac_kij_batch_bounded == nullptr
    ) {
        throw std::runtime_error("provider aqueous-kij capability unavailable");
    }
    return table;
}

std::array<const epcsaft_native_sdk_v1*, salt_count> parse_tables(
    PyObject* capsules
) {
    PyObject* sequence = PySequence_Fast(
        capsules, "provider capsules must be a sequence"
    );
    if (
        sequence == nullptr
        || PySequence_Fast_GET_SIZE(sequence)
            != static_cast<Py_ssize_t>(salt_count)
    ) {
        Py_XDECREF(sequence);
        throw std::invalid_argument(
            "provider capsules must contain six salt models"
        );
    }
    std::array<const epcsaft_native_sdk_v1*, salt_count> tables{};
    for (std::size_t index = 0; index < salt_count; ++index) {
        tables[index] = checked_table(
            PySequence_Fast_GET_ITEM(
                sequence, static_cast<Py_ssize_t>(index)
            ),
            index
        );
    }
    Py_DECREF(sequence);
    return tables;
}

Evaluation evaluate(
    const std::array<const epcsaft_native_sdk_v1*, salt_count>& tables,
    const Payload& payload,
    const std::array<double, parameter_count>& parameters,
    std::uint64_t deadline_ns,
    std::optional<std::size_t> active_coordinate = std::nullopt
) {
    Evaluation evaluation{};
    const epcsaft_native_evaluation_budget_v1 budget{
        sizeof(epcsaft_native_evaluation_budget_v1),
        EPCSAFT_NATIVE_EVALUATION_BUDGET_UNLIMITED_V1,
        deadline_ns,
    };
    const auto evaluate_salt = [&](std::size_t salt_index) {
        std::vector<std::size_t> row_indices;
        std::vector<double> molalities;
        for (std::size_t row_index = 0; row_index < row_count; ++row_index) {
            if (payload.observations[row_index].salt_index == salt_index) {
                row_indices.push_back(row_index);
                molalities.push_back(payload.observations[row_index].molality);
            }
        }
        const auto& columns = payload.observations[row_indices.front()].columns;
        if (
            active_coordinate.has_value()
            && std::find(
                columns.begin(), columns.end(), *active_coordinate
            ) == columns.end()
        ) {
            return;
        }
        for (const std::size_t row_index : row_indices) {
            if (payload.observations[row_index].columns != columns) {
                throw std::invalid_argument(
                    "rows for one salt do not share a parameter block"
                );
            }
        }
        std::array<double, 3> local_parameters{};
        for (std::size_t local = 0; local < local_parameters.size(); ++local) {
            local_parameters[local] = parameters[columns[local]];
        }
        std::vector<epcsaft_aqueous_miac_kij_result_v1> results(
            row_indices.size()
        );
        for (auto& result : results) {
            result.struct_size = sizeof(result);
        }
        const auto* table = tables[salt_index];
        const int status = table->evaluate_aqueous_miac_kij_batch_bounded(
            table->model_context,
            payload.expected_fingerprints[salt_index].c_str(),
            payload.temperature,
            payload.pressure,
            molalities.data(),
            molalities.size(),
            local_parameters.data(),
            local_parameters.size(),
            results.data(),
            results.size(),
            &budget
        );
        if (status != EPCSAFT_NATIVE_STATUS_OK_V1) {
            throw std::runtime_error(
                std::string("provider bounded aqueous-kij evaluation failed: ")
                + results.front().error
            );
        }
        for (
            std::size_t local_row = 0;
            local_row < row_indices.size();
            ++local_row
        ) {
            const std::size_t row_index = row_indices[local_row];
            const Observation& observed = payload.observations[row_index];
            const auto& result = results[local_row];
            const std::string fingerprint(
                result.parameter_fingerprint,
                strnlen(
                    result.parameter_fingerprint,
                    EPCSAFT_NATIVE_SDK_V1_FINGERPRINT_SIZE
                )
            );
            if (
                result.status != status
                || fingerprint
                    != payload.expected_fingerprints[observed.salt_index]
                || !std::isfinite(
                    result.log_mean_ionic_activity_coefficient_molality
                )
                || !std::isfinite(result.reference_molality_mol_per_kg)
                || !std::isfinite(result.reference_convergence_error)
                || !std::isfinite(
                    result.reference_derivative_convergence_error
                )
                || !std::all_of(
                    result.derivative,
                    result.derivative + 3,
                    [](double value) { return std::isfinite(value); }
                )
            ) {
                throw std::runtime_error(
                    "provider aqueous-kij row diagnostics failed"
                );
            }
            const double log_modeled =
                result.log_mean_ionic_activity_coefficient_molality;
            const double modeled = std::exp(log_modeled);
            const double ratio = modeled / observed.observed;
            const double residual = 1.0 - ratio;
            evaluation.residuals[row_index] = residual;
            std::array<double, 3> local_derivative{};
            for (std::size_t local = 0; local < 3; ++local) {
                local_derivative[local] = result.derivative[local];
                evaluation.jacobian[
                    row_index * parameter_count + columns[local]
                ] = -ratio * result.derivative[local];
            }
            evaluation.rows[row_index] = Row{
                observed.row_id,
                observed.salt,
                observed.molality,
                observed.observed,
                log_modeled,
                modeled,
                residual,
                local_derivative,
                result.reference_molality_mol_per_kg,
                result.reference_convergence_error,
                result.reference_derivative_convergence_error,
                fingerprint,
            };
        }
    };
    for (std::size_t salt_index = 0; salt_index < salt_count; ++salt_index) {
        evaluate_salt(salt_index);
    }
    return evaluation;
}

class CoordinateCost final : public ceres::CostFunction {
public:
    CoordinateCost(
        const std::array<const epcsaft_native_sdk_v1*, salt_count>& tables,
        const Payload& payload,
        const std::array<double, parameter_count>& fixed_parameters,
        std::size_t coordinate,
        std::uint64_t deadline_ns
    ) : tables_(tables),
        payload_(payload),
        fixed_parameters_(fixed_parameters),
        coordinate_(coordinate),
        deadline_ns_(deadline_ns) {
        for (std::size_t row = 0; row < row_count; ++row) {
            const auto& columns = payload.observations[row].columns;
            if (
                std::find(columns.begin(), columns.end(), coordinate_)
                != columns.end()
            ) {
                row_indices_.push_back(row);
            }
        }
        if (row_indices_.empty()) {
            throw std::invalid_argument(
                "coordinate has no source-bound observations"
            );
        }
        set_num_residuals(static_cast<int>(row_indices_.size()));
        mutable_parameter_block_sizes()->push_back(1);
    }

    bool Evaluate(
        double const* const* values,
        double* residuals,
        double** jacobians
    ) const override {
        try {
            if (!cached_ || values[0][0] != cached_value_) {
                auto parameters = fixed_parameters_;
                parameters[coordinate_] = values[0][0];
                cached_evaluation_ = evaluate(
                    tables_,
                    payload_,
                    parameters,
                    deadline_ns_,
                    coordinate_
                );
                cached_value_ = values[0][0];
                cached_ = true;
            }
        } catch (const std::exception& error) {
            failure_reason_ = error.what();
            return false;
        }
        for (
            std::size_t local_row = 0;
            local_row < row_indices_.size();
            ++local_row
        ) {
            const std::size_t row = row_indices_[local_row];
            residuals[local_row] = cached_evaluation_.residuals[row];
            if (jacobians != nullptr && jacobians[0] != nullptr) {
                jacobians[0][local_row] = cached_evaluation_.jacobian[
                    row * parameter_count + coordinate_
                ];
            }
        }
        return true;
    }

    const std::string& failure_reason() const {
        return failure_reason_;
    }

private:
    const std::array<const epcsaft_native_sdk_v1*, salt_count>& tables_;
    const Payload& payload_;
    std::array<double, parameter_count> fixed_parameters_;
    std::size_t coordinate_;
    std::uint64_t deadline_ns_;
    std::vector<std::size_t> row_indices_;
    mutable bool cached_{false};
    mutable double cached_value_{0.0};
    mutable Evaluation cached_evaluation_{};
    mutable std::string failure_reason_;
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

CoordinateSolve solve_coordinate(
    const std::array<const epcsaft_native_sdk_v1*, salt_count>& tables,
    const Payload& payload,
    const Schedule& schedule,
    std::size_t coordinate,
    std::uint64_t deadline_ns
) {
    auto fixed_parameters = payload.published_parameters;
    double value = schedule.values[coordinate];
    CoordinateCost cost(
        tables,
        payload,
        fixed_parameters,
        coordinate,
        deadline_ns
    );
    ceres::Problem::Options problem_options;
    problem_options.cost_function_ownership = ceres::DO_NOT_TAKE_OWNERSHIP;
    ceres::Problem problem(problem_options);
    problem.AddResidualBlock(&cost, nullptr, &value);
    problem.SetParameterLowerBound(&value, 0, payload.bounds[0]);
    problem.SetParameterUpperBound(&value, 0, payload.bounds[1]);
    ceres::Solver::Options options;
    options.linear_solver_type = ceres::DENSE_QR;
    options.max_num_iterations = payload.max_iterations;
    options.max_solver_time_in_seconds = payload.max_solver_time_seconds;
    options.function_tolerance = payload.function_tolerance;
    options.gradient_tolerance = payload.gradient_tolerance;
    options.parameter_tolerance = payload.parameter_tolerance;
    options.logging_type = ceres::SILENT;
    options.num_threads = 1;
    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);
    std::string failure_reason = cost.failure_reason();
    if (failure_reason.empty() && (
            summary.termination_type != ceres::CONVERGENCE
            || !summary.IsSolutionUsable() || !std::isfinite(value)
        )) {
        failure_reason = "terminated " + termination_name(
            summary.termination_type
        );
    }
    return CoordinateSolve{
        coordinate,
        value,
        summary.initial_cost,
        summary.final_cost,
        summary.iterations.size(),
        termination_name(summary.termination_type),
        summary.IsSolutionUsable() && std::isfinite(value),
        std::move(failure_reason),
    };
}

SolveOutcome solve_one(
    const std::array<const epcsaft_native_sdk_v1*, salt_count>& tables,
    const Payload& payload,
    const Schedule& schedule
) {
    SolveOutcome outcome{};
    outcome.name = schedule.name;
    outcome.coordinate_order = schedule.reverse ? "reverse" : "forward";
    outcome.parameters = payload.published_parameters;
    bool schedule_usable = true;
    outcome.initial_cost = 0.0;
    outcome.final_cost = 0.0;
    for (
        std::size_t order_index = 0;
        order_index < parameter_count;
        ++order_index
    ) {
        const std::size_t coordinate = schedule.reverse
            ? parameter_count - 1 - order_index
            : order_index;
        const CoordinateSolve solved = solve_coordinate(
            tables,
            payload,
            schedule,
            coordinate,
            deadline_after(payload.max_solver_time_seconds)
        );
        outcome.coordinate_solves[solved.coordinate] = solved;
        outcome.parameters[solved.coordinate] = solved.parameter;
        outcome.initial_cost += solved.initial_cost;
        outcome.final_cost += solved.final_cost;
        outcome.iterations += solved.iterations;
        if (
            solved.termination != "CONVERGENCE"
            || !solved.solution_usable
            || !solved.failure_reason.empty()
        ) {
            schedule_usable = false;
            if (!outcome.failure_reason.empty()) {
                outcome.failure_reason += "; ";
            }
            outcome.failure_reason +=
                "coordinate " + std::to_string(solved.coordinate)
                + ": "
                + (
                    solved.failure_reason.empty()
                    ? solved.termination
                    : solved.failure_reason
                );
        }
    }
    try {
        outcome.evaluation = evaluate(
            tables,
            payload,
            outcome.parameters,
            deadline_after(payload.max_solver_time_seconds)
        );
    } catch (const std::exception& error) {
        schedule_usable = false;
        if (!outcome.failure_reason.empty()) {
            outcome.failure_reason += "; ";
        }
        outcome.failure_reason += "final reporting evaluation: ";
        outcome.failure_reason += error.what();
        outcome.termination = "NO_CONVERGENCE";
        outcome.solution_usable = false;
        return outcome;
    }
    outcome.termination = schedule_usable ? "CONVERGENCE" : "NO_CONVERGENCE";
    outcome.solution_usable = schedule_usable;
    for (std::size_t column = 0; column < parameter_count; ++column) {
        bool complete = true;
        bool nonzero = false;
        for (std::size_t row = 0; row < row_count; ++row) {
            const double value =
                outcome.evaluation.jacobian[row * parameter_count + column];
            complete = complete && std::isfinite(value);
            nonzero = nonzero || value != 0.0;
        }
        outcome.complete_columns[column] = complete && nonzero;
    }
    Eigen::Map<
        const Eigen::Matrix<
            double,
            Eigen::Dynamic,
            Eigen::Dynamic,
            Eigen::RowMajor
        >
    > row_major(
        outcome.evaluation.jacobian.data(),
        static_cast<Eigen::Index>(row_count),
        static_cast<Eigen::Index>(parameter_count)
    );
    const Eigen::JacobiSVD<Eigen::MatrixXd> svd(
        Eigen::MatrixXd(row_major), Eigen::ComputeFullV
    );
    for (std::size_t index = 0; index < parameter_count; ++index) {
        outcome.singular_values[index] =
            svd.singularValues()[static_cast<Eigen::Index>(index)];
    }
    outcome.rank_threshold =
        outcome.singular_values.front()
        * static_cast<double>(std::max(row_count, parameter_count))
        * std::numeric_limits<double>::epsilon() * payload.rank_multiplier;
    outcome.rank = static_cast<int>(std::count_if(
        outcome.singular_values.begin(),
        outcome.singular_values.end(),
        [&](double value) { return value > outcome.rank_threshold; }
    ));
    if (outcome.singular_values.back() > 0.0) {
        outcome.condition_number =
            outcome.singular_values.front() / outcome.singular_values.back();
    }
    const auto least =
        svd.matrixV().col(static_cast<Eigen::Index>(parameter_count - 1));
    for (std::size_t index = 0; index < parameter_count; ++index) {
        outcome.least_sensitive_direction[index] =
            least[static_cast<Eigen::Index>(index)];
    }
    return outcome;
}

template <typename Container>
PyObject* doubles_to_tuple(const Container& values) {
    PyObject* tuple = PyTuple_New(static_cast<Py_ssize_t>(values.size()));
    if (tuple == nullptr) {
        return nullptr;
    }
    for (std::size_t index = 0; index < values.size(); ++index) {
        PyObject* value = PyFloat_FromDouble(values[index]);
        if (value == nullptr) {
            Py_DECREF(tuple);
            return nullptr;
        }
        PyTuple_SET_ITEM(tuple, static_cast<Py_ssize_t>(index), value);
    }
    return tuple;
}

PyObject* bools_to_tuple(
    const std::array<bool, parameter_count>& values
) {
    PyObject* tuple = PyTuple_New(static_cast<Py_ssize_t>(values.size()));
    if (tuple == nullptr) {
        return nullptr;
    }
    for (std::size_t index = 0; index < values.size(); ++index) {
        PyTuple_SET_ITEM(
            tuple,
            static_cast<Py_ssize_t>(index),
            Py_NewRef(values[index] ? Py_True : Py_False)
        );
    }
    return tuple;
}

PyObject* rows_to_tuple(const std::array<Row, row_count>& rows) {
    PyObject* tuple = PyTuple_New(static_cast<Py_ssize_t>(rows.size()));
    if (tuple == nullptr) {
        return nullptr;
    }
    for (std::size_t index = 0; index < rows.size(); ++index) {
        const Row& row = rows[index];
        PyObject* derivative = doubles_to_tuple(row.local_derivative);
        PyObject* item = PyTuple_New(12);
        if (derivative == nullptr || item == nullptr) {
            Py_XDECREF(derivative);
            Py_XDECREF(item);
            Py_DECREF(tuple);
            return nullptr;
        }
        PyTuple_SET_ITEM(item, 0, PyUnicode_FromString(row.row_id.c_str()));
        PyTuple_SET_ITEM(item, 1, PyUnicode_FromString(row.salt.c_str()));
        PyTuple_SET_ITEM(item, 2, PyFloat_FromDouble(row.molality));
        PyTuple_SET_ITEM(item, 3, PyFloat_FromDouble(row.observed));
        PyTuple_SET_ITEM(item, 4, PyFloat_FromDouble(row.log_modeled));
        PyTuple_SET_ITEM(item, 5, PyFloat_FromDouble(row.modeled));
        PyTuple_SET_ITEM(item, 6, PyFloat_FromDouble(row.residual));
        PyTuple_SET_ITEM(item, 7, derivative);
        PyTuple_SET_ITEM(item, 8, PyFloat_FromDouble(row.reference_molality));
        PyTuple_SET_ITEM(
            item, 9, PyFloat_FromDouble(row.reference_convergence)
        );
        PyTuple_SET_ITEM(
            item, 10, PyFloat_FromDouble(row.derivative_convergence)
        );
        PyTuple_SET_ITEM(
            item, 11, PyUnicode_FromString(row.fingerprint.c_str())
        );
        PyTuple_SET_ITEM(tuple, static_cast<Py_ssize_t>(index), item);
    }
    return tuple;
}

PyObject* evaluation_to_tuple(const Evaluation& evaluation) {
    PyObject* residuals = doubles_to_tuple(evaluation.residuals);
    PyObject* jacobian = doubles_to_tuple(evaluation.jacobian);
    PyObject* rows = rows_to_tuple(evaluation.rows);
    if (residuals == nullptr || jacobian == nullptr || rows == nullptr) {
        Py_XDECREF(residuals);
        Py_XDECREF(jacobian);
        Py_XDECREF(rows);
        return nullptr;
    }
    PyObject* result = PyTuple_New(3);
    if (result == nullptr) {
        Py_DECREF(residuals);
        Py_DECREF(jacobian);
        Py_DECREF(rows);
        return nullptr;
    }
    PyTuple_SET_ITEM(result, 0, residuals);
    PyTuple_SET_ITEM(result, 1, jacobian);
    PyTuple_SET_ITEM(result, 2, rows);
    return result;
}

PyObject* coordinate_solutions_to_tuple(
    const std::array<CoordinateSolve, parameter_count>& outcomes
);

PyObject* solution_to_tuple(const SolveOutcome& outcome) {
    PyObject* parameters = doubles_to_tuple(outcome.parameters);
    PyObject* rows = rows_to_tuple(outcome.evaluation.rows);
    PyObject* singular = doubles_to_tuple(outcome.singular_values);
    PyObject* least = doubles_to_tuple(outcome.least_sensitive_direction);
    PyObject* complete = bools_to_tuple(outcome.complete_columns);
    PyObject* coordinate_solves = coordinate_solutions_to_tuple(
        outcome.coordinate_solves
    );
    if (
        parameters == nullptr || rows == nullptr || singular == nullptr
        || least == nullptr || complete == nullptr
        || coordinate_solves == nullptr
    ) {
        Py_XDECREF(parameters);
        Py_XDECREF(rows);
        Py_XDECREF(singular);
        Py_XDECREF(least);
        Py_XDECREF(complete);
        Py_XDECREF(coordinate_solves);
        return nullptr;
    }
    PyObject* result = PyTuple_New(17);
    if (result == nullptr) {
        Py_DECREF(parameters);
        Py_DECREF(rows);
        Py_DECREF(singular);
        Py_DECREF(least);
        Py_DECREF(complete);
        Py_DECREF(coordinate_solves);
        return nullptr;
    }
    PyTuple_SET_ITEM(result, 0, PyUnicode_FromString(outcome.name.c_str()));
    PyTuple_SET_ITEM(
        result,
        1,
        PyUnicode_FromString(outcome.coordinate_order.c_str())
    );
    PyTuple_SET_ITEM(
        result,
        2,
        PyUnicode_FromString(outcome.termination.c_str())
    );
    PyTuple_SET_ITEM(
        result,
        3,
        Py_NewRef(outcome.solution_usable ? Py_True : Py_False)
    );
    PyTuple_SET_ITEM(
        result, 4, PyFloat_FromDouble(outcome.initial_cost)
    );
    PyTuple_SET_ITEM(
        result, 5, PyFloat_FromDouble(outcome.final_cost)
    );
    PyTuple_SET_ITEM(result, 6, PyLong_FromSize_t(outcome.iterations));
    PyTuple_SET_ITEM(result, 7, parameters);
    PyTuple_SET_ITEM(result, 8, rows);
    PyTuple_SET_ITEM(result, 9, singular);
    PyTuple_SET_ITEM(result, 10, PyFloat_FromDouble(outcome.rank_threshold));
    PyTuple_SET_ITEM(result, 11, PyLong_FromLong(outcome.rank));
    PyTuple_SET_ITEM(
        result, 12, PyFloat_FromDouble(outcome.condition_number)
    );
    PyTuple_SET_ITEM(result, 13, least);
    PyTuple_SET_ITEM(result, 14, complete);
    PyTuple_SET_ITEM(result, 15, coordinate_solves);
    PyTuple_SET_ITEM(
        result, 16, PyUnicode_FromString(outcome.failure_reason.c_str())
    );
    return result;
}

PyObject* coordinate_solution_to_tuple(const CoordinateSolve& outcome) {
    PyObject* result = PyTuple_New(8);
    if (result == nullptr) {
        return nullptr;
    }
    PyTuple_SET_ITEM(result, 0, PyLong_FromSize_t(outcome.coordinate));
    PyTuple_SET_ITEM(result, 1, PyFloat_FromDouble(outcome.parameter));
    PyTuple_SET_ITEM(result, 2, PyFloat_FromDouble(outcome.initial_cost));
    PyTuple_SET_ITEM(result, 3, PyFloat_FromDouble(outcome.final_cost));
    PyTuple_SET_ITEM(result, 4, PyLong_FromSize_t(outcome.iterations));
    PyTuple_SET_ITEM(
        result, 5, PyUnicode_FromString(outcome.termination.c_str())
    );
    PyTuple_SET_ITEM(
        result,
        6,
        Py_NewRef(outcome.solution_usable ? Py_True : Py_False)
    );
    PyTuple_SET_ITEM(
        result, 7, PyUnicode_FromString(outcome.failure_reason.c_str())
    );
    return result;
}

PyObject* coordinate_solutions_to_tuple(
    const std::array<CoordinateSolve, parameter_count>& outcomes
) {
    PyObject* result = PyTuple_New(static_cast<Py_ssize_t>(outcomes.size()));
    if (result == nullptr) {
        return nullptr;
    }
    for (std::size_t index = 0; index < outcomes.size(); ++index) {
        PyObject* item = coordinate_solution_to_tuple(outcomes[index]);
        if (item == nullptr) {
            Py_DECREF(result);
            return nullptr;
        }
        PyTuple_SET_ITEM(result, static_cast<Py_ssize_t>(index), item);
    }
    return result;
}

}  // namespace

PyObject* evaluate_figiel_kij_python(
    PyObject* capsules,
    PyObject* payload_object,
    PyObject* parameters_object
) {
    try {
        const Payload payload = parse_payload(payload_object);
        const auto tables = parse_tables(capsules);
        const auto parameters = doubles<parameter_count>(
            parameters_object, "Figiel interaction parameters"
        );
        return evaluation_to_tuple(evaluate(
            tables,
            payload,
            parameters,
            deadline_after(payload.max_solver_time_seconds)
        ));
    } catch (const std::exception& error) {
        if (PyErr_Occurred() != nullptr) {
            PyErr_Clear();
        }
        PyErr_SetString(PyExc_RuntimeError, error.what());
        return nullptr;
    }
}

PyObject* solve_figiel_kij_python(
    PyObject* capsules,
    PyObject* payload_object,
    std::size_t schedule_index
) {
    try {
        const Payload payload = parse_payload(payload_object);
        const auto tables = parse_tables(capsules);
        if (schedule_index >= payload.schedules.size()) {
            throw std::invalid_argument(
                "coordinate schedule index is outside the frozen contract"
            );
        }
        return solution_to_tuple(
            solve_one(tables, payload, payload.schedules[schedule_index])
        );
    } catch (const std::exception& error) {
        if (PyErr_Occurred() != nullptr) {
            PyErr_Clear();
        }
        PyErr_SetString(PyExc_RuntimeError, error.what());
        return nullptr;
    }
}

}  // namespace epcsaft_engine::regression
