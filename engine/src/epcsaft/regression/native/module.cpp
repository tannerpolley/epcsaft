#include "born_diameter_fit.hpp"
#include "ceres_core.hpp"
#include "evaluator_fit.hpp"
#include "figiel_kij_fit.hpp"
#include "figiel_water_factor_fit.hpp"
#include "general_fit.hpp"
#include "pure_saturation_fit.hpp"

#include <epcsaft/native_sdk_v1.h>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <exception>
#include <vector>

namespace epcsaft_engine::regression {
namespace {

struct NativeSdkPrefix final {
    std::uint32_t abi_version;
    std::size_t table_size;
    std::size_t result_size;
    const void* model_context;
    epcsaft_evaluate_pure_phase_v1 evaluate_pure_phase;
};

}  // namespace

const epcsaft_native_sdk_v1* checked_provider_table(PyObject* capsule) {
    if (!PyCapsule_CheckExact(capsule)) {
        PyErr_SetString(PyExc_TypeError, "provider transport must be an exact CPython capsule");
        return nullptr;
    }
    void* pointer = PyCapsule_GetPointer(capsule, EPCSAFT_NATIVE_SDK_V1_CAPSULE_NAME);
    if (pointer == nullptr) {
        return nullptr;
    }
    NativeSdkPrefix prefix{};
    std::memcpy(&prefix, pointer, sizeof(prefix));
    if (prefix.abi_version != EPCSAFT_NATIVE_SDK_V1_ABI_VERSION) {
        PyErr_SetString(PyExc_RuntimeError, "provider native SDK ABI version mismatch");
        return nullptr;
    }
    constexpr std::size_t minimum_size =
        offsetof(epcsaft_native_sdk_v1, evaluate_pure_phase_parameters)
        + sizeof(epcsaft_evaluate_pure_phase_parameters_v1);
    if (prefix.table_size < minimum_size) {
        PyErr_SetString(PyExc_RuntimeError, "provider native SDK table lacks the parameterized tail");
        return nullptr;
    }
    const auto* table = static_cast<const epcsaft_native_sdk_v1*>(pointer);
    if (table->parameterized_result_size != sizeof(epcsaft_parameterized_phase_block_result_v1)) {
        PyErr_SetString(PyExc_RuntimeError, "provider parameterized result size mismatch");
        return nullptr;
    }
    if (table->model_context == nullptr || table->evaluate_pure_phase_parameters == nullptr) {
        PyErr_SetString(PyExc_RuntimeError, "provider parameterized evaluator is unavailable");
        return nullptr;
    }
    return table;
}

Py_ssize_t minimum_parameterized_table_size() noexcept {
    return static_cast<Py_ssize_t>(
        offsetof(epcsaft_native_sdk_v1, evaluate_pure_phase_parameters)
        + sizeof(epcsaft_evaluate_pure_phase_parameters_v1)
    );
}

Py_ssize_t parameterized_result_size() noexcept {
    return static_cast<Py_ssize_t>(sizeof(epcsaft_parameterized_phase_block_result_v1));
}

PyObject* transport_info(PyObject* capsule) {
    const epcsaft_native_sdk_v1* table = checked_provider_table(capsule);
    if (table == nullptr) {
        return nullptr;
    }
    epcsaft_parameterized_phase_block_result_v1 result{};
    result.struct_size = sizeof(result);
    const int status = table->evaluate_pure_phase_parameters(
        table->model_context,
        130.0,
        1.0,
        1.0e-3,
        1.0,
        3.7039,
        150.03,
        &result
    );
    if (status != EPCSAFT_NATIVE_STATUS_OK_V1 || result.status != status) {
        PyErr_Format(
            PyExc_RuntimeError,
            "provider parameterized evaluator rejected the transport probe: %s",
            result.error
        );
        return nullptr;
    }
    const std::size_t fingerprint_length = strnlen(
        result.parameter_fingerprint, EPCSAFT_NATIVE_SDK_V1_FINGERPRINT_SIZE
    );
    return Py_BuildValue(
        "(knns#)",
        static_cast<unsigned long>(table->abi_version),
        static_cast<Py_ssize_t>(table->table_size),
        static_cast<Py_ssize_t>(table->parameterized_result_size),
        result.parameter_fingerprint,
        static_cast<Py_ssize_t>(fingerprint_length)
    );
}

}  // namespace epcsaft_engine::regression

namespace {

PyObject* py_transport_info(PyObject*, PyObject* capsule) {
    return epcsaft_engine::regression::transport_info(capsule);
}

PyObject* py_minimum_parameterized_table_size(PyObject*, PyObject*) {
    return PyLong_FromSsize_t(epcsaft_engine::regression::minimum_parameterized_table_size());
}

PyObject* py_parameterized_result_size(PyObject*, PyObject*) {
    return PyLong_FromSsize_t(epcsaft_engine::regression::parameterized_result_size());
}

PyObject* py_evaluate(PyObject*, PyObject* args) {
    PyObject* capsule = nullptr;
    PyObject* payload = nullptr;
    PyObject* variables = nullptr;
    if (!PyArg_ParseTuple(args, "OOO:evaluate", &capsule, &payload, &variables)) {
        return nullptr;
    }
    return epcsaft_engine::regression::evaluate_python(capsule, payload, variables);
}

PyObject* py_report_pure_saturation(PyObject*, PyObject* args) {
    PyObject* capsule = nullptr;
    PyObject* payload = nullptr;
    PyObject* reporting_rows = nullptr;
    PyObject* parameters = nullptr;
    if (!PyArg_ParseTuple(
            args,
            "OOOO:report_pure_saturation",
            &capsule,
            &payload,
            &reporting_rows,
            &parameters
        )) {
        return nullptr;
    }
    return epcsaft_engine::regression::report_python(
        capsule, payload, reporting_rows, parameters
    );
}

PyObject* py_evaluate_born(PyObject*, PyObject* args) {
    PyObject* capsules = nullptr;
    PyObject* payload = nullptr;
    PyObject* diameters = nullptr;
    if (!PyArg_ParseTuple(args, "OOO:evaluate_born", &capsules, &payload, &diameters)) {
        return nullptr;
    }
    return epcsaft_engine::regression::evaluate_born_python(capsules, payload, diameters);
}

PyObject* py_solve_born(PyObject*, PyObject* args) {
    PyObject* capsules = nullptr;
    PyObject* payload = nullptr;
    if (!PyArg_ParseTuple(args, "OO:solve_born", &capsules, &payload)) {
        return nullptr;
    }
    return epcsaft_engine::regression::solve_born_python(capsules, payload);
}

PyObject* py_solve_figiel_water_factor(PyObject*, PyObject* args) {
    PyObject* capsule = nullptr;
    PyObject* payload = nullptr;
    if (!PyArg_ParseTuple(
            args,
            "OO:solve_figiel_water_factor",
            &capsule,
            &payload
        )) {
        return nullptr;
    }
    return epcsaft_engine::regression::solve_figiel_water_factor_python(
        capsule, payload
    );
}

PyObject* py_evaluate_figiel_water_factor(PyObject*, PyObject* args) {
    PyObject* capsule = nullptr;
    PyObject* payload = nullptr;
    double parameter = 0.0;
    if (!PyArg_ParseTuple(
            args,
            "OOd:evaluate_figiel_water_factor",
            &capsule,
            &payload,
            &parameter
        )) {
        return nullptr;
    }
    return epcsaft_engine::regression::evaluate_figiel_water_factor_python(
        capsule, payload, parameter
    );
}

PyObject* py_solve_figiel_kij(PyObject*, PyObject* args) {
    PyObject* capsules = nullptr;
    PyObject* payload = nullptr;
    Py_ssize_t schedule_index = 0;
    if (!PyArg_ParseTuple(
            args,
            "OOn:solve_figiel_kij",
            &capsules,
            &payload,
            &schedule_index
        )) {
        return nullptr;
    }
    if (schedule_index < 0) {
        PyErr_SetString(PyExc_ValueError, "schedule index must be nonnegative");
        return nullptr;
    }
    return epcsaft_engine::regression::solve_figiel_kij_python(
        capsules,
        payload,
        static_cast<std::size_t>(schedule_index)
    );
}

PyObject* py_evaluate_figiel_kij(PyObject*, PyObject* args) {
    PyObject* capsules = nullptr;
    PyObject* payload = nullptr;
    PyObject* parameters = nullptr;
    if (!PyArg_ParseTuple(
            args,
            "OOO:evaluate_figiel_kij",
            &capsules,
            &payload,
            &parameters
        )) {
        return nullptr;
    }
    return epcsaft_engine::regression::evaluate_figiel_kij_python(
        capsules, payload, parameters
    );
}

PyObject* py_parameter_capabilities(PyObject*, PyObject* capsule) {
    return epcsaft_engine::regression::parameter_capabilities_python(capsule);
}

PyObject* py_evaluate_general(PyObject*, PyObject* args) {
    PyObject* capsule = nullptr;
    PyObject* payload = nullptr;
    PyObject* variables = nullptr;
    if (!PyArg_ParseTuple(
            args,
            "OOO:evaluate_general",
            &capsule,
            &payload,
            &variables
        )) {
        return nullptr;
    }
    return epcsaft_engine::regression::evaluate_general_python(
        capsule, payload, variables
    );
}

PyObject* py_solve_general(PyObject*, PyObject* args) {
    PyObject* capsule = nullptr;
    PyObject* payload = nullptr;
    if (!PyArg_ParseTuple(
            args, "OO:solve_general", &capsule, &payload
        )) {
        return nullptr;
    }
    return epcsaft_engine::regression::solve_general_python(capsule, payload);
}

PyObject* py_evaluate_general_start(PyObject*, PyObject* args) {
    PyObject* capsule = nullptr;
    PyObject* payload = nullptr;
    PyObject* physical_start = nullptr;
    if (!PyArg_ParseTuple(
            args,
            "OOO:evaluate_general_start",
            &capsule,
            &payload,
            &physical_start
        )) {
        return nullptr;
    }
    return epcsaft_engine::regression::evaluate_general_start_python(
        capsule, payload, physical_start
    );
}

PyObject* py_solve_evaluator(PyObject*, PyObject* args) {
    PyObject* capsule = nullptr;
    PyObject* payload = nullptr;
    if (!PyArg_ParseTuple(
            args, "OO:solve_evaluator", &capsule, &payload
        )) {
        return nullptr;
    }
    return epcsaft_engine::regression::solve_evaluator_python(capsule, payload);
}

PyObject* matrix_diagnostics_python(
    const epcsaft_engine::regression::internal::MatrixDiagnostics& diagnostics
) {
    PyObject* singular = PyTuple_New(
        static_cast<Py_ssize_t>(diagnostics.singular_values.size())
    );
    if (singular == nullptr) {
        return nullptr;
    }
    for (std::size_t index = 0; index < diagnostics.singular_values.size(); ++index) {
        PyObject* value = PyFloat_FromDouble(diagnostics.singular_values[index]);
        if (value == nullptr) {
            Py_DECREF(singular);
            return nullptr;
        }
        PyTuple_SET_ITEM(singular, static_cast<Py_ssize_t>(index), value);
    }
    PyObject* result = Py_BuildValue(
        "(Nid)", singular, diagnostics.rank, diagnostics.condition_number
    );
    return result;
}

PyObject* py_diagnose_jacobian(PyObject*, PyObject* args) {
    Py_ssize_t fitted = 0;
    Py_ssize_t lifted = 0;
    Py_ssize_t residuals = 0;
    PyObject* values = nullptr;
    if (!PyArg_ParseTuple(
            args, "nnnO:diagnose_jacobian",
            &fitted, &lifted, &residuals, &values
        )) {
        return nullptr;
    }
    if (fitted <= 0 || lifted < 0 || residuals <= 0) {
        PyErr_SetString(
            PyExc_ValueError,
            "fitted and residual counts must be positive and lifted count "
            "must be nonnegative"
        );
        return nullptr;
    }
    if (lifted > PY_SSIZE_T_MAX - fitted) {
        PyErr_SetString(PyExc_OverflowError, "Jacobian variable count overflows");
        return nullptr;
    }
    PyObject* sequence = PySequence_Fast(
        values, "Jacobian values must be a finite sequence"
    );
    if (sequence == nullptr) {
        return nullptr;
    }
    std::vector<double> jacobian(
        static_cast<std::size_t>(PySequence_Fast_GET_SIZE(sequence))
    );
    for (Py_ssize_t index = 0; index < PySequence_Fast_GET_SIZE(sequence); ++index) {
        jacobian[static_cast<std::size_t>(index)] = PyFloat_AsDouble(
            PySequence_Fast_GET_ITEM(sequence, index)
        );
        if (PyErr_Occurred()) {
            Py_DECREF(sequence);
            return nullptr;
        }
        if (!std::isfinite(jacobian[static_cast<std::size_t>(index)])) {
            Py_DECREF(sequence);
            PyErr_SetString(
                PyExc_ValueError,
                "Jacobian values must be finite"
            );
            return nullptr;
        }
    }
    Py_DECREF(sequence);
    try {
        const auto diagnostics = epcsaft_engine::regression::internal::diagnose_jacobian(
            {
                static_cast<std::size_t>(fitted),
                static_cast<std::size_t>(lifted),
                static_cast<std::size_t>(residuals),
            },
            jacobian
        );
        PyObject* full = matrix_diagnostics_python(diagnostics.full);
        PyObject* projected = matrix_diagnostics_python(
            diagnostics.projected_parameters
        );
        if (full == nullptr || projected == nullptr) {
            Py_XDECREF(full);
            Py_XDECREF(projected);
            return nullptr;
        }
        return Py_BuildValue("(NN)", full, projected);
    } catch (const std::exception& error) {
        PyErr_SetString(PyExc_ValueError, error.what());
        return nullptr;
    }
}

PyMethodDef methods[] = {
    {"transport_info", py_transport_info, METH_O, "Validate the Engine model capsule."},
    {
        "minimum_parameterized_table_size",
        py_minimum_parameterized_table_size,
        METH_NOARGS,
        "Return the minimum v1 parameterized table size."
    },
    {
        "parameterized_result_size",
        py_parameterized_result_size,
        METH_NOARGS,
        "Return the required parameterized result size."
    },
    {"evaluate", py_evaluate, METH_VARARGS, "Evaluate exact pure-saturation residuals and Jacobian."},
    {
        "report_pure_saturation",
        py_report_pure_saturation,
        METH_VARARGS,
        "Evaluate pure-saturation reporting closure at fixed fitted parameters."
    },
    {
        "evaluate_born",
        py_evaluate_born,
        METH_VARARGS,
        "Evaluate exact Figiel Born residuals and Jacobian."
    },
    {"solve_born", py_solve_born, METH_VARARGS, "Fit the five Figiel Born diameters."},
    {
        "solve_figiel_water_factor",
        py_solve_figiel_water_factor,
        METH_VARARGS,
        "Fit the Figiel NaBr water solvation factor."
    },
    {
        "evaluate_figiel_water_factor",
        py_evaluate_figiel_water_factor,
        METH_VARARGS,
        "Evaluate Figiel NaBr water-factor residuals and Jacobian."
    },
    {
        "solve_figiel_kij",
        py_solve_figiel_kij,
        METH_VARARGS,
        "Fit the eleven Figiel aqueous interaction parameters."
    },
    {
        "evaluate_figiel_kij",
        py_evaluate_figiel_kij,
        METH_VARARGS,
        "Evaluate exact Figiel aqueous interaction residuals and Jacobian."
    },
    {
        "parameter_capabilities",
        py_parameter_capabilities,
        METH_O,
        "Read exact model-bound Provider parameter capabilities."
    },
    {
        "evaluate_general",
        py_evaluate_general,
        METH_VARARGS,
        "Evaluate exact general neutral-binary residuals and Jacobian."
    },
    {
        "solve_general",
        py_solve_general,
        METH_VARARGS,
        "Fit one shared neutral-binary interaction parameter."
    },
    {
        "evaluate_general_start",
        py_evaluate_general_start,
        METH_VARARGS,
        "Evaluate a joint-pure problem at its exact resolved solver start."
    },
    {
        "diagnose_jacobian",
        py_diagnose_jacobian,
        METH_VARARGS,
        "Diagnose an exact general Jacobian without running Ceres."
    },
    {
        "solve_evaluator",
        py_solve_evaluator,
        METH_VARARGS,
        "Fit source-bound positive observations from an exact native evaluator."
    },
    {nullptr, nullptr, 0, nullptr},
};

PyModuleDef module = {
    PyModuleDef_HEAD_INIT,
    "_native",
    "Native Ceres regression transport.",
    -1,
    methods,
    nullptr,
    nullptr,
    nullptr,
    nullptr,
};

}  // namespace

PyMODINIT_FUNC PyInit__native() {
    return PyModule_Create(&module);
}
