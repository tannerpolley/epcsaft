#ifndef EPCSAFT_ENGINE_REGRESSION_BORN_DIAMETER_FIT_HPP
#define EPCSAFT_ENGINE_REGRESSION_BORN_DIAMETER_FIT_HPP

#include <Python.h>

namespace epcsaft_engine::regression {

PyObject* evaluate_born_python(
    PyObject* capsules,
    PyObject* payload,
    PyObject* diameters
);
PyObject* solve_born_python(PyObject* capsules, PyObject* payload);

}  // namespace epcsaft_engine::regression

#endif
