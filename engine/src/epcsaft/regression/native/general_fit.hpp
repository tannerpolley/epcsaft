#ifndef EPCSAFT_ENGINE_REGRESSION_GENERAL_FIT_HPP
#define EPCSAFT_ENGINE_REGRESSION_GENERAL_FIT_HPP

#include <Python.h>

namespace epcsaft_engine::regression {

PyObject* parameter_capabilities_python(PyObject* capsule);
PyObject* evaluate_general_python(
    PyObject* capsule, PyObject* payload, PyObject* variables
);
PyObject* evaluate_general_start_python(
    PyObject* capsule, PyObject* payload, PyObject* physical_start
);
PyObject* solve_general_python(PyObject* capsule, PyObject* payload);

}  // namespace epcsaft_engine::regression

#endif
