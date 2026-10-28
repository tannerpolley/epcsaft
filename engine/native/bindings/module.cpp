#include "core/convergence.hpp"
#include "core/declared.hpp"
#include "properties/helmholtz_model.hpp"

#include <nanobind/nanobind.h>
#include <nanobind/stl/optional.h>

#include <stdexcept>

namespace nb = nanobind;

void register_eos(nb::module_& module);
void register_diagnostics(nb::module_& module);
namespace epcsaft::equilibrium {
void register_equilibrium(nb::module_& module);
}
namespace epcsaft::regression {
void register_regression(nb::module_& module);
}
NB_MODULE(_core, module) {
    nb::exception<std::invalid_argument>(module, "ParameterError", PyExc_ValueError);
    nb::exception<std::runtime_error>(module, "EosError", PyExc_RuntimeError);
    nb::exception<epcsaft::DensityStructureError>(module, "DensityStructureError", module.attr("EosError").ptr());
    nb::exception<epcsaft::core::ConvergenceFailure>(module, "ConvergenceError", module.attr("EosError").ptr());
    nb::exception<epcsaft::core::DeclarationError>(module, "DeclarationError", module.attr("ParameterError").ptr());
    nb::exception<epcsaft::IonicReferenceBackground>(module, "IonicReferenceBackground", module.attr("ParameterError").ptr());
    register_eos(module);
    register_diagnostics(module);
    epcsaft::equilibrium::register_equilibrium(module);
    epcsaft::regression::register_regression(module);
    nb::register_exception_translator([](const std::exception_ptr& p, void* payload) {
        const nb::handle owner(static_cast<PyObject*>(payload));
        const auto publish = [&](const char* name, const std::exception& error, const auto& record,
                                 const std::optional<epcsaft::core::Convergence>& first) {
            auto type = owner.attr(name); auto instance = type(error.what());
            instance.attr("record") = nb::cast(record);
            instance.attr("first_numerical_failure") = nb::cast(first);
            PyErr_SetObject(type.ptr(), instance.ptr());
        };
        try { std::rethrow_exception(p); }
        catch (const epcsaft::core::ConvergenceFailure& e) { publish("ConvergenceError", e, e.record, e.first_numerical_failure); }
        catch (const epcsaft::core::DeclarationError& e) { publish("DeclarationError", e, e.record, std::nullopt); }
        catch (const epcsaft::DensityStructureError& e) { publish("DensityStructureError", e, e.first_numerical_failure, e.first_numerical_failure); }
    }, module.ptr());
}
