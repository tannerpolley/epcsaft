#pragma once

#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include <cstddef>

namespace epcsaft_engine::regression {

PyObject* solve_figiel_kij_python(
    PyObject* capsules,
    PyObject* payload,
    std::size_t schedule_index
);
PyObject* evaluate_figiel_kij_python(
    PyObject* capsules,
    PyObject* payload,
    PyObject* parameters
);

}  // namespace epcsaft_engine::regression
