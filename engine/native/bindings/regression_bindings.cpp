#include "regression/fit.hpp"

#include "core/trace.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
namespace epcsaft::regression {

void register_regression(py::module_& module) {
    py::enum_<ObservationKind>(module, "ObservationKind")
        .value("Density", ObservationKind::Density)
        .value("PhaseDensity", ObservationKind::PhaseDensity)
        .value("EquilibriumPressure", ObservationKind::EquilibriumPressure)
        .value("EquilibriumTemperature", ObservationKind::EquilibriumTemperature)
        .value("PhaseMoleFraction", ObservationKind::PhaseMoleFraction)
        .value("PhaseFraction", ObservationKind::PhaseFraction);
    py::enum_<ResidualForm>(module, "ResidualForm")
        .value("Difference", ResidualForm::Difference)
        .value("LogRatio", ResidualForm::LogRatio);
    py::enum_<ObservationStatus>(module, "ObservationStatus")
        .value("Available", ObservationStatus::Available)
        .value("DensityUnavailable", ObservationStatus::DensityUnavailable)
        .value("CentralSolveUnavailable", ObservationStatus::CentralSolveUnavailable)
        .value("TopologyEvent", ObservationStatus::TopologyEvent)
        .value("MissingMolarMass", ObservationStatus::MissingMolarMass)
        .value("NonpositiveLogArgument", ObservationStatus::NonpositiveLogArgument)
        .value("ActionUnavailable", ObservationStatus::ActionUnavailable)
        .value("Nonfinite", ObservationStatus::Nonfinite)
        .value("InadmissibleTrial", ObservationStatus::InadmissibleTrial);
    py::class_<Observation>(module, "Observation")
        .def(py::init<>())
        .def_readwrite("kind", &Observation::kind)
        .def_readwrite("temperature", &Observation::temperature)
        .def_readwrite("pressure", &Observation::pressure)
        .def_readwrite("composition", &Observation::composition)
        .def_readwrite("branch", &Observation::branch)
        .def_readwrite("problem", &Observation::problem)
        .def_readwrite("phase", &Observation::phase)
        .def_readwrite("species", &Observation::species)
        .def_readwrite("mass_basis", &Observation::mass_basis)
        .def_readwrite("observed", &Observation::observed)
        .def_readwrite("form", &Observation::form)
        .def_readwrite("scale", &Observation::scale);
    py::class_<Prediction>(module, "Prediction")
        .def_readonly("value", &Prediction::value)
        .def_readonly("unit", &Prediction::unit)
        .def_readonly("status", &Prediction::status)
        .def_readonly("message", &Prediction::message)
        .def_readonly("residual", &Prediction::residual)
        .def_readonly("convergence", &Prediction::convergence)
        .def_readonly("branch", &Prediction::branch)
        .def_readonly("stable_root_count", &Prediction::stable_root_count)
        .def_readonly("requested_tolerance_met", &Prediction::requested_tolerance_met);
    py::class_<Coordinate>(module, "Coordinate")
        .def(py::init<>())
        .def_readwrite("family", &Coordinate::family)
        .def_readwrite("slots", &Coordinate::slots)
        .def_readwrite("origin", &Coordinate::origin)
        .def_readwrite("scale", &Coordinate::scale)
        .def_readwrite("lower", &Coordinate::lower)
        .def_readwrite("upper", &Coordinate::upper);
    py::class_<CoordinateSet>(module, "CoordinateSet")
        .def(py::init<>())
        .def_readwrite("coordinates", &CoordinateSet::coordinates)
        .def("active", [](const CoordinateSet& set) { return set.active_parameters(); })
        .def("tie", [](const CoordinateSet& set) { return set.tie(); });
    py::class_<FitControls>(module, "FitControls")
        .def(py::init<>())
        .def_readwrite("function_tolerance", &FitControls::function_tolerance)
        .def_readwrite("gradient_tolerance", &FitControls::gradient_tolerance)
        .def_readwrite("parameter_tolerance", &FitControls::parameter_tolerance)
        .def_readwrite("maximum_iterations", &FitControls::maximum_iterations)
        .def_readwrite("maximum_elapsed_time_seconds", &FitControls::maximum_elapsed_time_seconds);
    py::class_<TrialFailure>(module, "TrialFailure")
        .def_readonly("trial", &TrialFailure::trial)
        .def_readonly("observation", &TrialFailure::observation)
        .def_readonly("optimizer", &TrialFailure::optimizer)
        .def_readonly("status", &TrialFailure::status)
        .def_readonly("message", &TrialFailure::message);
    py::class_<FitCovariance>(module, "FitCovariance")
        .def_readonly("status", &FitCovariance::status)
        .def_readonly("assumption", &FitCovariance::assumption)
        .def_readonly("rank", &FitCovariance::rank)
        .def_readonly("singular_values", &FitCovariance::singular_values)
        .def_readonly("active_bounds", &FitCovariance::active_bounds)
        .def_readonly("optimizer", &FitCovariance::optimizer)
        .def_readonly("physical", &FitCovariance::physical)
        .def_readonly("variance_factor", &FitCovariance::variance_factor);
    py::class_<FitResult>(module, "FitResult")
        .def_readonly("status", &FitResult::status)
        .def_readonly("message", &FitResult::message)
        .def_readonly("usable", &FitResult::usable)
        .def_readonly("physical", &FitResult::physical)
        .def_readonly("parameters", &FitResult::parameters)
        .def_readonly("predictions", &FitResult::predictions)
        .def_readonly("residuals", &FitResult::residuals)
        .def_readonly("weighted_residuals", &FitResult::weighted_residuals)
        .def_readonly("physical_jacobian", &FitResult::physical_jacobian)
        .def_readonly("optimizer_jacobian", &FitResult::optimizer_jacobian)
        .def_readonly("statuses", &FitResult::statuses)
        .def_readonly("failures", &FitResult::failures)
        .def_readonly("iterations", &FitResult::iterations)
        .def_readonly("residual_evaluations", &FitResult::residual_evaluations)
        .def_readonly("jacobian_evaluations", &FitResult::jacobian_evaluations)
        .def_readonly("trial_failures", &FitResult::trial_failures)
        .def_readonly("initial_cost", &FitResult::initial_cost)
        .def_readonly("final_cost", &FitResult::final_cost)
        .def_readonly("covariance", &FitResult::covariance);
    module.def("predict", [](const HelmholtzModel& model, const Observation& observation) {
        core::Scope scope("predict");
        validate(observation, model);
        return predict(model, observation);
    }, py::arg("model"), py::arg("observation"));
    module.def("fit", [](const Parameters& parameters, std::optional<ThermochemistryRecord> thermochemistry,
                         std::vector<Coordinate> coordinates, std::vector<Observation> observations,
                         std::vector<double> weights, const FitControls& controls) {
        core::Scope scope("fit");
        return fit({parameters, std::move(thermochemistry), {std::move(coordinates)}, std::move(observations),
                    std::move(weights), controls});
    }, py::arg("parameters"), py::arg("thermochemistry"), py::arg("coordinates"), py::arg("observations"),
       py::arg("weights"), py::arg("controls"));
}

}  // namespace epcsaft::regression
