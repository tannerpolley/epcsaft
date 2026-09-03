#include "core/convergence.hpp"
#include "core/declared.hpp"

#include <pybind11/pybind11.h>

#include <stdexcept>

namespace py = pybind11;

void register_eos(py::module_& module);
void register_diagnostics(py::module_& module);
namespace epcsaft::equilibrium {
void register_equilibrium(py::module_& module);
}
PYBIND11_MODULE(_core, module) {
    py::register_exception<std::invalid_argument>(module, "ParameterError", PyExc_ValueError);
    py::register_exception<std::runtime_error>(module, "EosError", PyExc_RuntimeError);
    py::register_exception<epcsaft::core::ConvergenceFailure>(module, "ConvergenceError", module.attr("EosError").ptr());
    py::register_exception<epcsaft::core::MissingDeclaration>(module, "DeclarationError", module.attr("ParameterError").ptr());
    register_eos(module);
    register_diagnostics(module);
    epcsaft::equilibrium::register_equilibrium(module);
}
