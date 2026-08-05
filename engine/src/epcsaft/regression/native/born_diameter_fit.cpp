#include "born_diameter_fit.hpp"

#include <epcsaft/native_sdk_v1.h>

#include <ceres/ceres.h>
#include <Eigen/Dense>
#include <Eigen/SVD>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace epcsaft_engine::regression {
namespace {

constexpr std::size_t ion_count = 5;

struct Target final {
    std::string target_id;
    std::string ion_label;
    std::string active_component_id;
    std::string counterion_component_id;
    double target_j_per_mol;
    double published_diameter_angstrom;
    std::string expected_fingerprint;
};

struct Payload final {
    std::vector<std::string> identity;
    std::array<Target, ion_count> targets;
    double temperature;
    double pressure;
    double reference_molality;
    double reference_convergence_max;
    double diameter_origin;
    double diameter_scale;
    std::array<double, 2> diameter_bounds;
    std::array<double, 2> scaled_bounds;
    std::array<std::array<double, ion_count>, 3> starts;
    int max_iterations;
    double function_tolerance;
    double gradient_tolerance;
    double parameter_tolerance;
    double rank_multiplier;
};

struct Row final {
    double value;
    double derivative;
    double raw_error;
    double scaled_residual;
    double scaled_jacobian;
    double reference_molality;
    double reference_convergence_error;
    std::string fingerprint;
};

struct Evaluation final {
    std::array<double, ion_count> residuals{};
    std::array<double, ion_count * ion_count> jacobian{};
    std::array<Row, ion_count> rows;
};

struct SolveOutcome final {
    std::string name;
    ceres::Solver::Summary summary;
    std::array<double, ion_count> transformed{};
    Evaluation evaluation;
    std::array<double, ion_count> singular_values{};
    int rank{0};
    double condition_number{std::numeric_limits<double>::infinity()};
    bool complete_columns{false};
    std::string failure_reason;
};

std::string text(PyObject* object, const char* label) {
    if (!PyUnicode_Check(object)) {
        throw std::invalid_argument(std::string(label) + " must be text");
    }
    Py_ssize_t size = 0;
    const char* value = PyUnicode_AsUTF8AndSize(object, &size);
    if (value == nullptr) {
        throw std::invalid_argument(std::string(label) + " must be UTF-8 text");
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

template <std::size_t Size>
std::array<double, Size> fixed_doubles(PyObject* object, const char* label) {
    PyObject* sequence = PySequence_Fast(object, label);
    if (sequence == nullptr) throw std::invalid_argument(label);
    if (PySequence_Fast_GET_SIZE(sequence) != static_cast<Py_ssize_t>(Size)) {
        Py_DECREF(sequence);
        throw std::invalid_argument(std::string(label) + " has the wrong length");
    }
    std::array<double, Size> result{};
    for (std::size_t index = 0; index < Size; ++index) {
        result[index] = number(
            PySequence_Fast_GET_ITEM(sequence, static_cast<Py_ssize_t>(index)), label
        );
    }
    Py_DECREF(sequence);
    return result;
}

Payload parse_payload(PyObject* object) {
    PyObject* sequence = PySequence_Fast(object, "Born payload must be a sequence");
    if (sequence == nullptr) throw std::invalid_argument("Born payload must be a sequence");
    if (PySequence_Fast_GET_SIZE(sequence) != 16) {
        Py_DECREF(sequence);
        throw std::invalid_argument("Born payload has the wrong length");
    }
    Payload payload{};
    PyObject* identity = PySequence_Fast(
        PySequence_Fast_GET_ITEM(sequence, 0), "Born identity must be a sequence"
    );
    if (identity == nullptr || PySequence_Fast_GET_SIZE(identity) == 0) {
        Py_XDECREF(identity);
        Py_DECREF(sequence);
        throw std::invalid_argument("Born identity must be a nonempty sequence");
    }
    const Py_ssize_t identity_size = PySequence_Fast_GET_SIZE(identity);
    payload.identity.reserve(static_cast<std::size_t>(identity_size));
    for (Py_ssize_t index = 0; index < identity_size; ++index) {
        payload.identity.push_back(text(PySequence_Fast_GET_ITEM(identity, index), "Born identity"));
    }
    Py_DECREF(identity);

    PyObject* targets = PySequence_Fast(
        PySequence_Fast_GET_ITEM(sequence, 1), "Born targets must be a sequence"
    );
    if (targets == nullptr || PySequence_Fast_GET_SIZE(targets) != 5) {
        Py_XDECREF(targets);
        Py_DECREF(sequence);
        throw std::invalid_argument("Born targets must contain exactly five rows");
    }
    for (std::size_t index = 0; index < ion_count; ++index) {
        PyObject* row = PySequence_Fast(
            PySequence_Fast_GET_ITEM(targets, static_cast<Py_ssize_t>(index)),
            "Born target row must be a sequence"
        );
        if (row == nullptr || PySequence_Fast_GET_SIZE(row) != 7) {
            Py_XDECREF(row);
            Py_DECREF(targets);
            Py_DECREF(sequence);
            throw std::invalid_argument("Born target row has the wrong length");
        }
        payload.targets[index] = Target{
            text(PySequence_Fast_GET_ITEM(row, 0), "target id"),
            text(PySequence_Fast_GET_ITEM(row, 1), "ion label"),
            text(PySequence_Fast_GET_ITEM(row, 2), "active component id"),
            text(PySequence_Fast_GET_ITEM(row, 3), "counterion component id"),
            number(PySequence_Fast_GET_ITEM(row, 4), "target"),
            number(PySequence_Fast_GET_ITEM(row, 5), "published diameter"),
            text(PySequence_Fast_GET_ITEM(row, 6), "expected fingerprint"),
        };
        Py_DECREF(row);
    }
    Py_DECREF(targets);
    payload.temperature = number(PySequence_Fast_GET_ITEM(sequence, 2), "temperature");
    payload.pressure = number(PySequence_Fast_GET_ITEM(sequence, 3), "pressure");
    payload.reference_molality = number(PySequence_Fast_GET_ITEM(sequence, 4), "reference molality");
    payload.reference_convergence_max = number(
        PySequence_Fast_GET_ITEM(sequence, 5), "reference convergence maximum"
    );
    payload.diameter_origin = number(PySequence_Fast_GET_ITEM(sequence, 6), "diameter origin");
    payload.diameter_scale = number(PySequence_Fast_GET_ITEM(sequence, 7), "diameter scale");
    payload.diameter_bounds = fixed_doubles<2>(
        PySequence_Fast_GET_ITEM(sequence, 8), "diameter bounds"
    );
    payload.scaled_bounds = fixed_doubles<2>(
        PySequence_Fast_GET_ITEM(sequence, 9), "scaled bounds"
    );
    PyObject* starts = PySequence_Fast(
        PySequence_Fast_GET_ITEM(sequence, 10), "Born starts must be a sequence"
    );
    if (starts == nullptr || PySequence_Fast_GET_SIZE(starts) != 3) {
        Py_XDECREF(starts);
        Py_DECREF(sequence);
        throw std::invalid_argument("Born starts must contain primary, lower, and upper");
    }
    for (std::size_t index = 0; index < 3; ++index) {
        payload.starts[index] = fixed_doubles<ion_count>(
            PySequence_Fast_GET_ITEM(starts, static_cast<Py_ssize_t>(index)), "Born start"
        );
    }
    Py_DECREF(starts);
    const long max_iterations = PyLong_AsLong(PySequence_Fast_GET_ITEM(sequence, 11));
    if (PyErr_Occurred() != nullptr) {
        PyErr_Clear();
        Py_DECREF(sequence);
        throw std::invalid_argument("maximum iterations must be an integer");
    }
    payload.max_iterations = static_cast<int>(max_iterations);
    payload.function_tolerance = number(PySequence_Fast_GET_ITEM(sequence, 12), "function tolerance");
    payload.gradient_tolerance = number(PySequence_Fast_GET_ITEM(sequence, 13), "gradient tolerance");
    payload.parameter_tolerance = number(PySequence_Fast_GET_ITEM(sequence, 14), "parameter tolerance");
    payload.rank_multiplier = number(PySequence_Fast_GET_ITEM(sequence, 15), "rank multiplier");
    Py_DECREF(sequence);

    constexpr std::array<const char*, ion_count> expected_target_ids{
        "figiel2025-s5-Lip-reported-average", "figiel2025-s5-Nap-reported-average",
        "figiel2025-s5-Kp-reported-average", "figiel2025-s5-Clm-reported-average",
        "figiel2025-s5-Brm-reported-average",
    };
    constexpr std::array<const char*, ion_count> expected_ions{"Li+", "Na+", "K+", "Cl-", "Br-"};
    constexpr std::array<const char*, ion_count> expected_active{
        "lithium-cation", "sodium-cation", "potassium-cation", "chloride-anion",
        "bromide-anion",
    };
    constexpr std::array<const char*, ion_count> expected_counter{
        "chloride-anion", "chloride-anion", "chloride-anion", "sodium-cation",
        "sodium-cation",
    };
    constexpr std::array<double, ion_count> expected_values{
        -486200.0, -381100.0, -309100.0, -314900.0, -290900.0,
    };
    constexpr std::array<double, ion_count> expected_published{2.784, 3.445, 4.150, 4.100, 4.480};
    constexpr std::array<const char*, ion_count> expected_fingerprints{
        "sha256:9cef6373572a729c7c3519ee04c261b4f9d821d37ae2b74dc916a4a02e863149",
        "sha256:e97c462061bf92d2d75622a8d025a4f7e94106edccf2c924c88142b1b6a03439",
        "sha256:1e885c57ebf7226d28a153ad4c3afb5e9ed8786b04aa0ff1b65e776b157d67d3",
        "sha256:ac5a3ee0bf044d9fc63086a41f04f7b9f34a506daf67ac9627fd082a519c62a3",
        "sha256:7db8316e0c2352d992c7482d4c1218463020d749d5365ec8a8f5e491643eb9b2",
    };
    for (std::size_t index = 0; index < ion_count; ++index) {
        const Target& target = payload.targets[index];
        if (target.target_id != expected_target_ids[index]
            || target.ion_label != expected_ions[index]
            || target.active_component_id != expected_active[index]
            || target.counterion_component_id != expected_counter[index]
            || target.target_j_per_mol != expected_values[index]
            || target.published_diameter_angstrom != expected_published[index]
            || target.expected_fingerprint != expected_fingerprints[index]) {
            throw std::invalid_argument("Born targets do not match the compiled five-row contract");
        }
    }
    if (payload.temperature != 298.15 || payload.pressure != 100000.0
        || payload.reference_molality != 1.0e-12
        || payload.reference_convergence_max != 5.0e-5
        || payload.diameter_origin != 3.0 || payload.diameter_scale != 1.0
        || payload.diameter_bounds != std::array<double, 2>{1.0, 6.0}
        || payload.scaled_bounds != std::array<double, 2>{-2.0, 3.0}
        || payload.starts[0] != std::array<double, ion_count>{3.0, 3.0, 3.0, 3.0, 3.0}
        || payload.starts[1] != std::array<double, ion_count>{2.0, 2.0, 2.0, 2.0, 2.0}
        || payload.starts[2] != std::array<double, ion_count>{5.0, 5.0, 5.0, 5.0, 5.0}
        || payload.max_iterations != 500 || payload.function_tolerance != 1.0e-10
        || payload.gradient_tolerance != 1.0e-10 || payload.parameter_tolerance != 1.0e-10
        || payload.rank_multiplier != 100.0) {
        throw std::invalid_argument("Born numerical contract does not match the frozen design");
    }
    return payload;
}

const epcsaft_native_sdk_v1* checked_born_table(PyObject* capsule, const Target& target) {
    if (!PyCapsule_CheckExact(capsule)) {
        throw std::invalid_argument("provider transport must be an exact CPython capsule");
    }
    void* pointer = PyCapsule_GetPointer(capsule, EPCSAFT_NATIVE_SDK_V1_CAPSULE_NAME);
    if (pointer == nullptr) throw std::runtime_error("provider capability unavailable");
    struct Prefix final {
        std::uint32_t abi_version;
        std::size_t table_size;
    } prefix{};
    std::memcpy(&prefix, pointer, sizeof(prefix));
    constexpr std::size_t minimum_size =
        offsetof(epcsaft_native_sdk_v1, evaluate_ion_solvation_born)
        + sizeof(epcsaft_evaluate_ion_solvation_born_v1);
    if (prefix.abi_version != EPCSAFT_NATIVE_SDK_V1_ABI_VERSION
        || prefix.table_size < minimum_size) {
        throw std::runtime_error("provider capability unavailable");
    }
    const auto* table = static_cast<const epcsaft_native_sdk_v1*>(pointer);
    if (table->ion_solvation_born_result_size != sizeof(epcsaft_ion_solvation_born_result_v1)
        || table->model_context == nullptr || table->evaluate_ion_solvation_born == nullptr) {
        throw std::runtime_error("provider capability unavailable");
    }
    if (table->component_count != 3 || table->component_ids == nullptr
        || table->component_charges == nullptr || table->component_ids[0] == nullptr
        || table->component_ids[1] == nullptr || table->component_ids[2] == nullptr
        || std::string(table->component_ids[0]) != "water"
        || std::string(table->component_ids[1]) != target.active_component_id
        || std::string(table->component_ids[2]) != target.counterion_component_id) {
        throw std::invalid_argument("provider component order does not match the Born target");
    }
    const int expected_active_charge = target.ion_label.back() == '+' ? 1 : -1;
    if (table->component_charges[0] != 0
        || table->component_charges[1] != expected_active_charge
        || table->component_charges[2] != -expected_active_charge) {
        throw std::invalid_argument("provider charge order does not match the Born target");
    }
    return table;
}

std::array<const epcsaft_native_sdk_v1*, ion_count> parse_tables(
    PyObject* capsules,
    const Payload& payload
) {
    PyObject* sequence = PySequence_Fast(capsules, "provider capsules must be a sequence");
    if (sequence == nullptr || PySequence_Fast_GET_SIZE(sequence) != 5) {
        Py_XDECREF(sequence);
        throw std::invalid_argument("provider capsules must contain exactly five models");
    }
    std::array<const epcsaft_native_sdk_v1*, ion_count> tables{};
    for (std::size_t index = 0; index < ion_count; ++index) {
        tables[index] = checked_born_table(
            PySequence_Fast_GET_ITEM(sequence, static_cast<Py_ssize_t>(index)),
            payload.targets[index]
        );
    }
    Py_DECREF(sequence);
    return tables;
}

Evaluation evaluate(
    const std::array<const epcsaft_native_sdk_v1*, ion_count>& tables,
    const Payload& payload,
    const std::array<double, ion_count>& diameters
) {
    Evaluation evaluation{};
    for (std::size_t index = 0; index < ion_count; ++index) {
        if (!std::isfinite(diameters[index]) || diameters[index] <= 0.0) {
            throw std::invalid_argument("Born diameter must be finite and positive");
        }
        epcsaft_ion_solvation_born_result_v1 result{};
        result.struct_size = sizeof(result);
        const int status = tables[index]->evaluate_ion_solvation_born(
            tables[index]->model_context,
            payload.temperature,
            payload.pressure,
            diameters[index],
            &result
        );
        if (status != EPCSAFT_NATIVE_STATUS_OK_V1 || result.status != status) {
            throw std::runtime_error(std::string("provider Born evaluation failed: ") + result.error);
        }
        const std::string fingerprint(
            result.parameter_fingerprint,
            strnlen(result.parameter_fingerprint, EPCSAFT_NATIVE_SDK_V1_FINGERPRINT_SIZE)
        );
        if (fingerprint != payload.targets[index].expected_fingerprint) {
            throw std::runtime_error("provider Born fingerprint mismatch");
        }
        if (!std::isfinite(result.solvation_gibbs_j_per_mol)
            || !std::isfinite(result.derivative_j_per_mol_per_angstrom)
            || result.reference_molality_mol_per_kg != payload.reference_molality
            || !std::isfinite(result.reference_convergence_error)
            || result.reference_convergence_error > payload.reference_convergence_max) {
            throw std::runtime_error("provider Born reference diagnostics failed");
        }
        const double scale = std::abs(payload.targets[index].target_j_per_mol);
        const double raw = result.solvation_gibbs_j_per_mol
            - payload.targets[index].target_j_per_mol;
        const double scaled = raw / scale;
        const double jacobian = result.derivative_j_per_mol_per_angstrom
            * payload.diameter_scale / scale;
        evaluation.residuals[index] = scaled;
        evaluation.jacobian[index * ion_count + index] = jacobian;
        evaluation.rows[index] = Row{
            result.solvation_gibbs_j_per_mol,
            result.derivative_j_per_mol_per_angstrom,
            raw,
            scaled,
            jacobian,
            result.reference_molality_mol_per_kg,
            result.reference_convergence_error,
            fingerprint,
        };
    }
    return evaluation;
}

class BornCost final : public ceres::CostFunction {
public:
    BornCost(
        const std::array<const epcsaft_native_sdk_v1*, ion_count>& tables,
        const Payload& payload
    ) : tables_(tables), payload_(payload) {
        set_num_residuals(static_cast<int>(ion_count));
        mutable_parameter_block_sizes()->push_back(static_cast<int>(ion_count));
    }

    bool Evaluate(
        double const* const* parameters,
        double* residuals,
        double** jacobians
    ) const override {
        try {
            std::array<double, ion_count> diameters{};
            for (std::size_t index = 0; index < ion_count; ++index) {
                diameters[index] = payload_.diameter_origin
                    + payload_.diameter_scale * parameters[0][index];
            }
            const Evaluation result = evaluate(tables_, payload_, diameters);
            std::copy(result.residuals.begin(), result.residuals.end(), residuals);
            if (jacobians != nullptr && jacobians[0] != nullptr) {
                std::copy(result.jacobian.begin(), result.jacobian.end(), jacobians[0]);
            }
            failure_reason_.clear();
            return true;
        } catch (const std::exception& error) {
            failure_reason_ = error.what();
            return false;
        }
    }

    const std::string& failure_reason() const noexcept { return failure_reason_; }

private:
    std::array<const epcsaft_native_sdk_v1*, ion_count> tables_;
    const Payload& payload_;
    mutable std::string failure_reason_;
};

std::string termination_name(ceres::TerminationType termination) {
    switch (termination) {
        case ceres::CONVERGENCE: return "CONVERGENCE";
        case ceres::NO_CONVERGENCE: return "NO_CONVERGENCE";
        case ceres::FAILURE: return "FAILURE";
        case ceres::USER_SUCCESS: return "USER_SUCCESS";
        case ceres::USER_FAILURE: return "USER_FAILURE";
    }
    return "UNKNOWN";
}

SolveOutcome solve_one(
    const std::string& name,
    const std::array<const epcsaft_native_sdk_v1*, ion_count>& tables,
    const Payload& payload,
    const std::array<double, ion_count>& start
) {
    SolveOutcome outcome{};
    outcome.name = name;
    for (std::size_t index = 0; index < ion_count; ++index) {
        outcome.transformed[index] = (start[index] - payload.diameter_origin)
            / payload.diameter_scale;
    }
    BornCost cost(tables, payload);
    ceres::Problem::Options problem_options;
    problem_options.cost_function_ownership = ceres::DO_NOT_TAKE_OWNERSHIP;
    ceres::Problem problem(problem_options);
    problem.AddResidualBlock(&cost, nullptr, outcome.transformed.data());
    for (std::size_t index = 0; index < ion_count; ++index) {
        problem.SetParameterLowerBound(outcome.transformed.data(), index, payload.scaled_bounds[0]);
        problem.SetParameterUpperBound(outcome.transformed.data(), index, payload.scaled_bounds[1]);
    }
    ceres::Solver::Options options;
    options.linear_solver_type = ceres::DENSE_QR;
    options.max_num_iterations = payload.max_iterations;
    options.function_tolerance = payload.function_tolerance;
    options.gradient_tolerance = payload.gradient_tolerance;
    options.parameter_tolerance = payload.parameter_tolerance;
    options.logging_type = ceres::SILENT;
    options.num_threads = 1;
    ceres::Solve(options, &problem, &outcome.summary);
    if (!cost.failure_reason().empty()) outcome.failure_reason = cost.failure_reason();
    std::array<double, ion_count> diameters{};
    for (std::size_t index = 0; index < ion_count; ++index) {
        diameters[index] = payload.diameter_origin
            + payload.diameter_scale * outcome.transformed[index];
    }
    outcome.evaluation = evaluate(tables, payload, diameters);
    outcome.complete_columns = std::all_of(
        outcome.evaluation.jacobian.begin(), outcome.evaluation.jacobian.end(),
        [](double value) { return std::isfinite(value); }
    );
    Eigen::Matrix<double, 5, 5, Eigen::RowMajor> matrix;
    std::copy(
        outcome.evaluation.jacobian.begin(), outcome.evaluation.jacobian.end(), matrix.data()
    );
    const Eigen::JacobiSVD<Eigen::Matrix<double, 5, 5, Eigen::RowMajor>> svd(
        matrix, Eigen::ComputeFullU | Eigen::ComputeFullV
    );
    for (std::size_t index = 0; index < ion_count; ++index) {
        outcome.singular_values[index] = svd.singularValues()[static_cast<Eigen::Index>(index)];
    }
    const double threshold = outcome.singular_values[0] * static_cast<double>(ion_count)
        * std::numeric_limits<double>::epsilon() * payload.rank_multiplier;
    outcome.rank = static_cast<int>(std::count_if(
        outcome.singular_values.begin(), outcome.singular_values.end(),
        [threshold](double value) { return value > threshold; }
    ));
    if (outcome.singular_values.back() > 0.0) {
        outcome.condition_number = outcome.singular_values.front()
            / outcome.singular_values.back();
    }
    return outcome;
}

PyObject* doubles_to_tuple(const double* values, std::size_t size) {
    PyObject* tuple = PyTuple_New(static_cast<Py_ssize_t>(size));
    if (tuple == nullptr) return nullptr;
    for (std::size_t index = 0; index < size; ++index) {
        PyObject* value = PyFloat_FromDouble(values[index]);
        if (value == nullptr) {
            Py_DECREF(tuple);
            return nullptr;
        }
        PyTuple_SET_ITEM(tuple, static_cast<Py_ssize_t>(index), value);
    }
    return tuple;
}

PyObject* strings_to_tuple(const std::vector<std::string>& values) {
    PyObject* tuple = PyTuple_New(static_cast<Py_ssize_t>(values.size()));
    if (tuple == nullptr) return nullptr;
    for (std::size_t index = 0; index < values.size(); ++index) {
        PyObject* value = PyUnicode_FromStringAndSize(
            values[index].data(), static_cast<Py_ssize_t>(values[index].size())
        );
        if (value == nullptr) {
            Py_DECREF(tuple);
            return nullptr;
        }
        PyTuple_SET_ITEM(tuple, static_cast<Py_ssize_t>(index), value);
    }
    return tuple;
}

PyObject* rows_to_tuple(const std::array<Row, ion_count>& rows) {
    PyObject* tuple = PyTuple_New(5);
    if (tuple == nullptr) return nullptr;
    for (std::size_t index = 0; index < ion_count; ++index) {
        const Row& row = rows[index];
        PyObject* item = Py_BuildValue(
            "(ddddddds)",
            row.value,
            row.derivative,
            row.raw_error,
            row.scaled_residual,
            row.scaled_jacobian,
            row.reference_molality,
            row.reference_convergence_error,
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

PyObject* fingerprints_to_tuple(const std::array<Row, ion_count>& rows) {
    PyObject* tuple = PyTuple_New(5);
    if (tuple == nullptr) return nullptr;
    for (std::size_t index = 0; index < ion_count; ++index) {
        PyObject* item = PyUnicode_FromString(rows[index].fingerprint.c_str());
        if (item == nullptr) {
            Py_DECREF(tuple);
            return nullptr;
        }
        PyTuple_SET_ITEM(tuple, static_cast<Py_ssize_t>(index), item);
    }
    return tuple;
}

PyObject* solution_to_tuple(const SolveOutcome& outcome) {
    PyObject* transformed = doubles_to_tuple(outcome.transformed.data(), ion_count);
    PyObject* residuals = doubles_to_tuple(outcome.evaluation.residuals.data(), ion_count);
    PyObject* jacobian = doubles_to_tuple(
        outcome.evaluation.jacobian.data(), ion_count * ion_count
    );
    PyObject* rows = rows_to_tuple(outcome.evaluation.rows);
    PyObject* singular = doubles_to_tuple(outcome.singular_values.data(), ion_count);
    if (transformed == nullptr || residuals == nullptr || jacobian == nullptr
        || rows == nullptr || singular == nullptr) {
        Py_XDECREF(transformed);
        Py_XDECREF(residuals);
        Py_XDECREF(jacobian);
        Py_XDECREF(rows);
        Py_XDECREF(singular);
        return nullptr;
    }
    PyObject* tuple = PyTuple_New(15);
    if (tuple == nullptr) {
        Py_DECREF(transformed);
        Py_DECREF(residuals);
        Py_DECREF(jacobian);
        Py_DECREF(rows);
        Py_DECREF(singular);
        return nullptr;
    }
    PyTuple_SET_ITEM(tuple, 0, PyUnicode_FromString(outcome.name.c_str()));
    PyTuple_SET_ITEM(tuple, 1, PyUnicode_FromString(termination_name(outcome.summary.termination_type).c_str()));
    PyTuple_SET_ITEM(tuple, 2, Py_NewRef(outcome.summary.IsSolutionUsable() ? Py_True : Py_False));
    PyTuple_SET_ITEM(tuple, 3, PyFloat_FromDouble(outcome.summary.initial_cost));
    PyTuple_SET_ITEM(tuple, 4, PyFloat_FromDouble(outcome.summary.final_cost));
    PyTuple_SET_ITEM(tuple, 5, PyLong_FromSize_t(outcome.summary.iterations.size()));
    PyTuple_SET_ITEM(tuple, 6, transformed);
    PyTuple_SET_ITEM(tuple, 7, residuals);
    PyTuple_SET_ITEM(tuple, 8, jacobian);
    PyTuple_SET_ITEM(tuple, 9, rows);
    PyTuple_SET_ITEM(tuple, 10, singular);
    PyTuple_SET_ITEM(tuple, 11, PyLong_FromLong(outcome.rank));
    PyTuple_SET_ITEM(tuple, 12, PyFloat_FromDouble(outcome.condition_number));
    PyTuple_SET_ITEM(tuple, 13, Py_NewRef(outcome.complete_columns ? Py_True : Py_False));
    PyTuple_SET_ITEM(tuple, 14, PyUnicode_FromString(outcome.failure_reason.c_str()));
    return tuple;
}

PyObject* evaluation_to_python(const Evaluation& evaluation, const Payload& payload) {
    PyObject* residuals = doubles_to_tuple(evaluation.residuals.data(), ion_count);
    PyObject* jacobian = doubles_to_tuple(evaluation.jacobian.data(), ion_count * ion_count);
    PyObject* rows = rows_to_tuple(evaluation.rows);
    PyObject* fingerprints = fingerprints_to_tuple(evaluation.rows);
    PyObject* identity = strings_to_tuple(payload.identity);
    if (residuals == nullptr || jacobian == nullptr || rows == nullptr
        || fingerprints == nullptr || identity == nullptr) {
        Py_XDECREF(residuals);
        Py_XDECREF(jacobian);
        Py_XDECREF(rows);
        Py_XDECREF(fingerprints);
        Py_XDECREF(identity);
        return nullptr;
    }
    PyObject* result = PyTuple_New(5);
    if (result == nullptr) {
        Py_DECREF(residuals);
        Py_DECREF(jacobian);
        Py_DECREF(rows);
        Py_DECREF(fingerprints);
        Py_DECREF(identity);
        return nullptr;
    }
    PyTuple_SET_ITEM(result, 0, residuals);
    PyTuple_SET_ITEM(result, 1, jacobian);
    PyTuple_SET_ITEM(result, 2, rows);
    PyTuple_SET_ITEM(result, 3, fingerprints);
    PyTuple_SET_ITEM(result, 4, identity);
    return result;
}

}  // namespace

PyObject* evaluate_born_python(
    PyObject* capsules,
    PyObject* payload_object,
    PyObject* diameters_object
) {
    try {
        const Payload payload = parse_payload(payload_object);
        const auto tables = parse_tables(capsules, payload);
        const auto diameters = fixed_doubles<ion_count>(diameters_object, "Born diameters");
        return evaluation_to_python(evaluate(tables, payload, diameters), payload);
    } catch (const std::exception& error) {
        if (PyErr_Occurred() != nullptr) PyErr_Clear();
        PyErr_SetString(PyExc_RuntimeError, error.what());
        return nullptr;
    }
}

PyObject* solve_born_python(PyObject* capsules, PyObject* payload_object) {
    try {
        const Payload payload = parse_payload(payload_object);
        const auto tables = parse_tables(capsules, payload);
        constexpr std::array<const char*, 3> names{"primary", "lower", "upper"};
        PyObject* solutions = PyTuple_New(3);
        if (solutions == nullptr) return nullptr;
        for (std::size_t index = 0; index < 3; ++index) {
            const SolveOutcome outcome = solve_one(names[index], tables, payload, payload.starts[index]);
            PyObject* item = solution_to_tuple(outcome);
            if (item == nullptr) {
                Py_DECREF(solutions);
                return nullptr;
            }
            PyTuple_SET_ITEM(solutions, static_cast<Py_ssize_t>(index), item);
        }
        PyObject* identity = strings_to_tuple(payload.identity);
        if (identity == nullptr) {
            Py_DECREF(solutions);
            return nullptr;
        }
        PyObject* result = PyTuple_New(2);
        if (result == nullptr) {
            Py_DECREF(solutions);
            Py_DECREF(identity);
            return nullptr;
        }
        PyTuple_SET_ITEM(result, 0, solutions);
        PyTuple_SET_ITEM(result, 1, identity);
        return result;
    } catch (const std::exception& error) {
        if (PyErr_Occurred() != nullptr) PyErr_Clear();
        PyErr_SetString(PyExc_RuntimeError, error.what());
        return nullptr;
    }
}

}  // namespace epcsaft_engine::regression
