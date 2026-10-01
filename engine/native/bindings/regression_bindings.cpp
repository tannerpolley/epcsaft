#include "regression/fit.hpp"

#include "core/trace.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
namespace epcsaft::regression {

void register_regression(py::module_& module) {
    py::enum_<ObservationKind>(module, "ObservationKind")
        .value("MeanIonicActivityCoefficient", ObservationKind::MeanIonicActivityCoefficient)
        .value("OsmoticCoefficient", ObservationKind::OsmoticCoefficient)
        .value("SolvationGibbsEnergy", ObservationKind::SolvationGibbsEnergy)
        .value("TransferGibbsEnergy", ObservationKind::TransferGibbsEnergy)
        .value("BulkRelativePermittivity", ObservationKind::BulkRelativePermittivity)
        .value("IsobaricHeatCapacity", ObservationKind::IsobaricHeatCapacity)
        .value("LogEquilibriumConstant", ObservationKind::LogEquilibriumConstant)
        .value("Density", ObservationKind::Density)
        .value("PhaseDensity", ObservationKind::PhaseDensity)
        .value("EquilibriumPressure", ObservationKind::EquilibriumPressure)
        .value("EquilibriumTemperature", ObservationKind::EquilibriumTemperature)
        .value("PhaseMoleFraction", ObservationKind::PhaseMoleFraction)
        .value("PhaseFraction", ObservationKind::PhaseFraction)
        .value("PartialPressure", ObservationKind::PartialPressure)
        .value("PhaseAmountRatio", ObservationKind::PhaseAmountRatio)
        .value("AdditionEnthalpy", ObservationKind::AdditionEnthalpy);
    py::enum_<ResidualForm>(module, "ResidualForm")
        .value("Difference", ResidualForm::Difference)
        .value("LogRatio", ResidualForm::LogRatio);
    py::enum_<ObservationStatus>(module, "ObservationStatus")
        .value("NotEvaluated", ObservationStatus::NotEvaluated)
        .value("ReferenceUnavailable", ObservationStatus::ReferenceUnavailable)
        .value("Available", ObservationStatus::Available)
        .value("ThermochemistryUnavailable", ObservationStatus::ThermochemistryUnavailable)
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
        .def_readwrite("ion_indices", &Observation::ion_indices)
        .def_readwrite("stoichiometry", &Observation::stoichiometry)
        .def_readwrite("molality", &Observation::molality)
        .def_readwrite("reference_composition", &Observation::reference_composition)
        .def_readwrite("branch", &Observation::branch)
        .def_readwrite("problem", &Observation::problem)
        .def_readwrite("final_problem", &Observation::final_problem)
        .def_readwrite("feed_amount", &Observation::feed_amount)
        .def_readwrite("reaction", &Observation::reaction)
        .def_readwrite("phase", &Observation::phase)
        .def_readwrite("species", &Observation::species)
        .def_readwrite("denominator_phase", &Observation::denominator_phase)
        .def_readwrite("denominator_species", &Observation::denominator_species)
        .def_readwrite("mass_basis", &Observation::mass_basis)
        .def_readwrite("observed", &Observation::observed)
        .def_readwrite("form", &Observation::form)
        .def_readwrite("scale", &Observation::scale);
    py::enum_<CheckMetric>(module, "CheckMetric")
        .value("AbsoluteDifference", CheckMetric::AbsoluteDifference).value("AbsoluteLog10Ratio", CheckMetric::AbsoluteLog10Ratio);
    py::class_<HydrationHold>(module, "HydrationHold")
        .def(py::init<>())
        .def_readwrite("reference", &HydrationHold::reference)
        .def_readwrite("lower", &HydrationHold::lower)
        .def_readwrite("upper", &HydrationHold::upper)
        .def_readwrite("energy_tolerance", &HydrationHold::energy_tolerance)
        .def_readwrite("diameter_tolerance", &HydrationHold::diameter_tolerance)
        .def_readwrite("maximum_iterations", &HydrationHold::maximum_iterations);
    py::class_<HoldResult>(module, "HoldResult")
        .def_readonly("energy_error", &HoldResult::energy_error)
        .def_readonly("correction", &HoldResult::correction)
        .def_readonly("uncertainty", &HoldResult::uncertainty)
        .def_readonly("derivative", &HoldResult::derivative)
        .def_readonly("parameter_derivative", &HoldResult::parameter_derivative)
        .def_readonly("iterations", &HoldResult::iterations)
        .def_readonly("cause", &HoldResult::cause);
    py::class_<EvaluatedCoordinates>(module, "EvaluatedCoordinates")
        .def_readonly("parameters", &EvaluatedCoordinates::parameters)
        .def_readonly("active", &EvaluatedCoordinates::active)
        .def_readonly("physical", &EvaluatedCoordinates::physical)
        .def_readonly("map", &EvaluatedCoordinates::map)
        .def_readonly("holds", &EvaluatedCoordinates::holds);
    py::class_<CheckGroup>(module, "CheckGroup")
        .def(py::init<>())
        .def_readwrite("name", &CheckGroup::name)
        .def_readwrite("rows", &CheckGroup::rows)
        .def_readwrite("metric", &CheckGroup::metric)
        .def_readwrite("maximum", &CheckGroup::maximum)
        .def_readwrite("mean", &CheckGroup::mean)
        .def_readwrite("baseline", &CheckGroup::baseline)
        .def_readwrite("report_only", &CheckGroup::report_only);
    py::class_<CheckResult>(module, "CheckResult")
        .def_readonly("name", &CheckResult::name)
        .def_readonly("unit", &CheckResult::unit)
        .def_readonly("status", &CheckResult::status)
        .def_readonly("rows", &CheckResult::rows)
        .def_readonly("failed_rows", &CheckResult::failed_rows)
        .def_readonly("errors", &CheckResult::errors)
        .def_readonly("slopes", &CheckResult::slopes)
        .def_readonly("maximum", &CheckResult::maximum)
        .def_readonly("mean", &CheckResult::mean)
        .def_readonly("argmax", &CheckResult::argmax)
        .def_readonly("maximum_threshold", &CheckResult::maximum_threshold)
        .def_readonly("mean_threshold", &CheckResult::mean_threshold)
        .def_readonly("baseline_maximum", &CheckResult::baseline_maximum)
        .def_readonly("baseline_mean", &CheckResult::baseline_mean)
        .def_readonly("report_only", &CheckResult::report_only);
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
    py::class_<PredictionGradient>(module, "PredictionGradient")
        .def_readonly("value", &PredictionGradient::value)
        .def_readonly("status", &PredictionGradient::status)
        .def_readonly("message", &PredictionGradient::message);
    py::enum_<ReactionCoefficient>(module, "ReactionCoefficient")
        .value("A", ReactionCoefficient::A).value("B", ReactionCoefficient::B)
        .value("C", ReactionCoefficient::C).value("D", ReactionCoefficient::D);
    py::class_<ReactionCoordinate>(module, "ReactionCoordinate")
        .def(py::init<>())
        .def_readwrite("name", &ReactionCoordinate::name)
        .def_readwrite("coefficient", &ReactionCoordinate::coefficient);
    py::class_<Coordinate>(module, "Coordinate")
        .def(py::init<>())
        .def_readwrite("family", &Coordinate::family)
        .def_readwrite("slots", &Coordinate::slots)
        .def_readwrite("reaction", &Coordinate::reaction)
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
        .def_readonly("message", &TrialFailure::message)
        .def_readonly("jacobian", &TrialFailure::jacobian);
    py::class_<FitCovariance>(module, "FitCovariance")
        .def_readonly("status", &FitCovariance::status)
        .def_readonly("assumption", &FitCovariance::assumption)
        .def_readonly("rank", &FitCovariance::rank)
        .def_readonly("singular_values", &FitCovariance::singular_values)
        .def_readonly("left_vectors", &FitCovariance::left_vectors)
        .def_readonly("right_vectors", &FitCovariance::right_vectors)
        .def_readonly("projected_gradient_norm", &FitCovariance::projected_gradient_norm)
        .def_readonly("objective_rows", &FitCovariance::objective_rows)
        .def_readonly("active_bounds", &FitCovariance::active_bounds)
        .def_readonly("optimizer", &FitCovariance::optimizer)
        .def_readonly("physical", &FitCovariance::physical)
        .def_readonly("variance_factor", &FitCovariance::variance_factor);
    py::class_<FitResult>(module, "FitResult")
        .def_readonly("status", &FitResult::status)
        .def_readonly("chart", &FitResult::chart)
        .def_readonly("checks", &FitResult::checks)
        .def_readonly("acceptance", &FitResult::acceptance)
        .def_readonly("optimizer", &FitResult::optimizer)
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
        .def_readonly("iteration_seconds", &FitResult::iteration_seconds)
        .def_readonly("iteration_costs", &FitResult::iteration_costs)
        .def_readonly("initial_cost", &FitResult::initial_cost)
        .def_readonly("final_cost", &FitResult::final_cost)
        .def_readonly("covariance", &FitResult::covariance);
    module.def("predict", [](const HelmholtzModel& model, const Observation& observation) {
        core::Scope scope("predict");
        validate(observation, model);
        return predict(model, observation);
    }, py::arg("model"), py::arg("observation"));
    module.def("differentiate", [](const HelmholtzModel& model, const Observation& observation) {
        core::Scope scope("differentiate");
        validate(observation, model);
        return differentiate(model, observation, predict(model, observation));
    }, py::arg("model"), py::arg("observation"));
    module.def("fit", [](const Parameters& parameters, std::optional<ThermochemistryRecord> thermochemistry,
                         std::vector<Coordinate> coordinates, std::vector<Observation> observations,
                         std::vector<double> weights, const FitControls& controls, std::vector<HydrationHold> holds,
                         std::vector<CheckGroup> checks, std::optional<Parameters> baseline, bool singular_vectors,
                         const std::optional<std::vector<std::vector<double>>>& starts, std::optional<int> iterations, std::optional<double> seconds) -> py::object {
        core::Scope scope("fit");
        FitRequest request{parameters, std::move(thermochemistry), {std::move(coordinates)}, std::move(observations),
                    std::move(weights), controls, std::move(holds), std::move(checks), std::move(baseline), singular_vectors};
        return starts ? py::cast(fit_sequence(request, *starts, iterations.value_or(controls.maximum_iterations), seconds.value_or(controls.maximum_elapsed_time_seconds))) : py::cast(fit(request));
    }, py::arg("parameters"), py::arg("thermochemistry"), py::arg("coordinates"), py::arg("observations"),
       py::arg("weights"), py::arg("controls"), py::arg("holds") = std::vector<HydrationHold>{},
       py::arg("checks") = std::vector<CheckGroup>{}, py::arg("baseline") = std::nullopt, py::arg("singular_vectors") = false,
       py::arg("starts") = std::nullopt, py::arg("iterations") = std::nullopt, py::arg("seconds") = std::nullopt);
    module.def("observations_at", [](const Parameters& parameters, std::vector<Coordinate> coordinates,
        const std::vector<double>& optimizer, std::vector<HydrationHold> holds, std::vector<Observation> observations, const std::vector<double>& errors) {
        FitRequest request; request.parameters = parameters; request.coordinates.coordinates = std::move(coordinates);
        request.holds = std::move(holds); request.observations = std::move(observations);
        return observations_at(request, optimizer, errors);
    });
    module.def("evaluate_coordinates", &evaluate_coordinates, py::arg("parameters"), py::arg("set"), py::arg("optimizer"),
        py::arg("holds"), py::arg("thermochemistry") = std::nullopt);
    module.def("evaluate_checks", [](const Parameters& stored, std::vector<Observation> observations,
        std::vector<CheckGroup> checks, std::optional<Parameters> baseline) {
        FitRequest request; request.observations = std::move(observations); request.checks = std::move(checks); request.baseline = std::move(baseline);
        return evaluate_checks(request, stored);
    }, py::arg("stored"), py::arg("observations"), py::arg("checks"), py::arg("baseline") = std::nullopt);
}

}  // namespace epcsaft::regression
