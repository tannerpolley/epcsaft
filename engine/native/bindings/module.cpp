#include <pybind11/pybind11.h>

#include <stdexcept>

namespace py = pybind11;

void register_eos(py::module_& module);

PYBIND11_MODULE(_core, module) {
    py::register_exception<std::invalid_argument>(module, "ParameterError", PyExc_ValueError);
    py::register_exception<std::runtime_error>(module, "EosError", PyExc_RuntimeError);
    register_eos(module);
}
