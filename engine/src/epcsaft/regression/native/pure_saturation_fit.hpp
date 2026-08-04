#ifndef EPCSAFT_ENGINE_REGRESSION_PURE_SATURATION_FIT_HPP
#define EPCSAFT_ENGINE_REGRESSION_PURE_SATURATION_FIT_HPP

#include <Python.h>
#include <epcsaft/native_sdk_v1.h>

namespace epcsaft_engine::regression {

PyObject* transport_info(PyObject* capsule);
const epcsaft_native_sdk_v1* checked_provider_table(PyObject* capsule);
PyObject* evaluate_python(PyObject* capsule, PyObject* payload, PyObject* variables);
PyObject* report_python(
    PyObject* capsule,
    PyObject* payload,
    PyObject* reporting_rows,
    PyObject* parameters
);
Py_ssize_t minimum_parameterized_table_size() noexcept;
Py_ssize_t parameterized_result_size() noexcept;

}  // namespace epcsaft_engine::regression

#endif
