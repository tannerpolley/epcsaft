#include "core/convergence.hpp"
#include "core/numerical_policy.hpp"
#include "core/trace.hpp"

#include "fields.hpp"
#include "solver_options.hpp"

#include <nanobind/nanobind.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <optional>

namespace nb = nanobind;
namespace {

struct CaptureWrapper final {
    std::optional<epcsaft::core::Capture> capture;
    std::optional<epcsaft::core::Trace> captured;
    bool diagnostics = false;
    explicit CaptureWrapper(bool detail) : diagnostics(detail) {}
};

nb::object tolerance_reference(const epcsaft::policy::Tolerance& tolerance) {
    if (tolerance.upstream == nullptr) return nb::none();
    return nb::cast(tolerance.upstream, nb::rv_policy::reference);
}

}  // namespace

void register_diagnostics(nb::module_& module) {
    using namespace epcsaft;
    auto stops = nb::enum_<core::Stop>(module, "Stop");
    for (int value = static_cast<int>(core::Stop::converged);
         value <= static_cast<int>(core::Stop::nonfinite); ++value) {
        const auto stop = static_cast<core::Stop>(value);
        stops.value(core::stop_name(stop), stop);
    }
    nb::enum_<policy::Category>(module, "Category")
        EPCSAFT_VALUES(policy::Category, ThermodynamicConstraint, AlgorithmTolerance);
    nb::class_<policy::Tolerance>(module, "Tolerance")
        EPCSAFT_READONLY(policy::Tolerance, name, value, category, source, factor)
        .def_prop_ro("upstream", tolerance_reference);
    module.def("tolerances", [] {
        nb::list entries;
        for (const auto* entry : policy::ENTRIES) entries.append(nb::cast(entry, nb::rv_policy::reference));
        return entries;
    }, "Immutable source defaults for every numerical policy; result.solver_options reports effective stopping settings.");
    nb::class_<core::Convergence>(module, "Convergence")
        .def_prop_ro("criterion", [](const core::Convergence& record) -> nb::object {
            if (record.criterion == nullptr) return nb::none();
            // A value copy outlives the record that may own a caller-valued criterion.
            return nb::cast(*record.criterion);
        })
        EPCSAFT_READONLY(core::Convergence, residual, floor, amplification, contraction, iterations, evaluations, stop)
        .def_prop_ro("residual_over_criterion", [](const core::Convergence& record) -> nb::object {
            if (record.criterion == nullptr) return nb::none();
            return nb::float_(record.residual / record.criterion->value);
        })
        .def_prop_ro("residual_over_floor", [](const core::Convergence& record) -> nb::object {
            if (!record.floor.has_value()) return nb::none();
            return nb::float_(record.residual / *record.floor);
        })
        .def("__repr__", [](const core::Convergence& record) {
            return core::describe(record);
        });
    nb::class_<core::Event>(module, "Event")
        EPCSAFT_READONLY(core::Event, name, span, iteration, values, reason);
    nb::class_<core::Span>(module, "Span")
        EPCSAFT_READONLY(core::Span, name, parent, convergence, wall_us);
    nb::class_<core::Trace>(module, "Trace")
        EPCSAFT_READONLY(core::Trace, spans, events, current, incomplete);
    nb::class_<CaptureWrapper>(module, "Capture")
        .def(nb::init<bool>(),
             nb::arg("diagnostics") = false)
        .def("__enter__", [](CaptureWrapper& wrapper) -> CaptureWrapper& {
            wrapper.capture.emplace(wrapper.diagnostics);
            wrapper.captured.reset();
            return wrapper;
        }, nb::rv_policy::reference_internal)
        .def("__exit__", [](CaptureWrapper& wrapper, nb::object, nb::object, nb::object) {
            wrapper.captured.emplace(wrapper.capture->take());
            wrapper.capture.reset();
            return false;
        }, nb::arg("exc_type").none(), nb::arg("exc_value").none(), nb::arg("traceback").none())
        .def_prop_ro("trace", [](const CaptureWrapper& wrapper) -> nb::object {
            if (!wrapper.captured.has_value()) return nb::none();
            return nb::cast(&*wrapper.captured, nb::rv_policy::reference_internal, nb::find(&wrapper));
        })
        .def_prop_ro("spans", [](const CaptureWrapper& wrapper) -> const std::vector<core::Span>& {
            if (!wrapper.captured.has_value())
                throw std::runtime_error("trace is not available until the capture closes");
            return wrapper.captured->spans;
        }, nb::rv_policy::reference_internal)
        .def_prop_ro("events", [](const CaptureWrapper& wrapper) -> const std::vector<core::Event>& {
            if (!wrapper.captured.has_value())
                throw std::runtime_error("trace is not available until the capture closes");
            return wrapper.captured->events;
        }, nb::rv_policy::reference_internal);
}
