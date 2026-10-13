#include "core/convergence.hpp"
#include "core/numerical_policy.hpp"
#include "core/trace.hpp"

#include "fields.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <optional>

namespace py = pybind11;
namespace {

struct CaptureWrapper final {
    std::optional<epcsaft::core::Capture> capture;
    std::optional<epcsaft::core::Trace> captured;
    bool diagnostics = false;
    explicit CaptureWrapper(bool detail) : diagnostics(detail) {}
};

py::object tolerance_reference(const epcsaft::policy::Tolerance& tolerance) {
    if (tolerance.upstream == nullptr) return py::none();
    return py::cast(tolerance.upstream, py::return_value_policy::reference);
}

}  // namespace

void register_diagnostics(py::module_& module) {
    using namespace epcsaft;
    auto stops = py::enum_<core::Stop>(module, "Stop");
    for (int value = static_cast<int>(core::Stop::converged);
         value <= static_cast<int>(core::Stop::nonfinite); ++value) {
        const auto stop = static_cast<core::Stop>(value);
        stops.value(core::stop_name(stop), stop);
    }
    py::enum_<policy::Category>(module, "Category")
        EPCSAFT_VALUES(policy::Category, ThermodynamicConstraint, AlgorithmTolerance);
    py::class_<policy::Tolerance>(module, "Tolerance")
        EPCSAFT_READONLY(policy::Tolerance, name, value, category, source, factor)
        .def_property_readonly("upstream", tolerance_reference);
    module.def("tolerances", [] {
        py::list entries;
        for (const auto* entry : policy::ENTRIES) entries.append(py::cast(entry, py::return_value_policy::reference));
        return entries;
    }, "Every named numerical policy entry (read-only); values change only in numerical_policy.hpp.");
    py::class_<core::Convergence>(module, "Convergence")
        .def_property_readonly("criterion", [](const core::Convergence& record) -> py::object {
            if (record.criterion == nullptr) return py::none();
            // A value copy outlives the record that may own a caller-valued criterion.
            return py::cast(*record.criterion);
        })
        EPCSAFT_READONLY(core::Convergence, residual, floor, amplification, contraction, iterations, evaluations, stop)
        .def_property_readonly("residual_over_criterion", [](const core::Convergence& record) -> py::object {
            if (record.criterion == nullptr) return py::none();
            return py::float_(record.residual / record.criterion->value);
        })
        .def_property_readonly("residual_over_floor", [](const core::Convergence& record) -> py::object {
            if (!record.floor.has_value()) return py::none();
            return py::float_(record.residual / *record.floor);
        })
        .def_property_readonly("amplification_over_declared", [](const core::Convergence& record) -> py::object {
            if (!record.amplification.has_value() || record.criterion == nullptr) return py::none();
            return py::float_(*record.amplification / record.criterion->factor);
        })
        .def("__repr__", [](const core::Convergence& record) {
            return core::describe(record);
        });
    py::class_<core::Event>(module, "Event")
        EPCSAFT_READONLY(core::Event, name, span, iteration, values, reason);
    py::class_<core::Span>(module, "Span")
        EPCSAFT_READONLY(core::Span, name, parent, convergence, wall_us);
    py::class_<core::Trace>(module, "Trace")
        EPCSAFT_READONLY(core::Trace, spans, events, current, incomplete);
    py::class_<CaptureWrapper>(module, "Capture")
        .def(py::init<bool>(),
             py::arg("diagnostics") = false)
        .def("__enter__", [](CaptureWrapper& wrapper) -> CaptureWrapper& {
            wrapper.capture.emplace(wrapper.diagnostics);
            wrapper.captured.reset();
            return wrapper;
        }, py::return_value_policy::reference_internal)
        .def("__exit__", [](CaptureWrapper& wrapper, py::object, py::object, py::object) {
            wrapper.captured.emplace(wrapper.capture->take());
            wrapper.capture.reset();
            return false;
        })
        .def_property_readonly("trace", [](const CaptureWrapper& wrapper) -> py::object {
            if (!wrapper.captured.has_value()) return py::none();
            return py::cast(&*wrapper.captured, py::return_value_policy::reference_internal, py::cast(wrapper));
        })
        .def_property_readonly("spans", [](const CaptureWrapper& wrapper) -> const std::vector<core::Span>& {
            if (!wrapper.captured.has_value())
                throw std::runtime_error("trace is not available until the capture closes");
            return wrapper.captured->spans;
        }, py::return_value_policy::reference_internal)
        .def_property_readonly("events", [](const CaptureWrapper& wrapper) -> const std::vector<core::Event>& {
            if (!wrapper.captured.has_value())
                throw std::runtime_error("trace is not available until the capture closes");
            return wrapper.captured->events;
        }, py::return_value_policy::reference_internal);
}
