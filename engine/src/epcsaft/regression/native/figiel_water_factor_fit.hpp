#pragma once

#include <Python.h>

namespace epcsaft_engine::regression {

PyObject* solve_figiel_water_factor_python(PyObject* capsule, PyObject* payload);
PyObject* evaluate_figiel_water_factor_python(
    PyObject* capsule, PyObject* payload, double parameter
);

}  // namespace epcsaft_engine::regression
