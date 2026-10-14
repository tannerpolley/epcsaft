#include "regression/fit.hpp"

#include "core/trace.hpp"

#include "fields.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <numeric>

namespace py = pybind11;
namespace epcsaft::regression {
namespace {

// The point readers' request: fit's model inputs at a declared point, without optimizer controls.
FitRequest point_request(const Parameters& parameters, std::vector<Coordinate> coordinates, std::vector<HydrationHold> holds,
                         std::vector<Observation> observations, std::optional<ThermochemistryRecord> thermochemistry) {
    FitRequest request; request.parameters = parameters; request.thermochemistry = std::move(thermochemistry);
    request.coordinates.coordinates = std::move(coordinates); request.holds = std::move(holds); request.observations = std::move(observations);
    return request;
}

}  // namespace

void register_regression(py::module_& module) {
    // Fit workers rerun this interpreter on the extension image they inherit as fd 3 (#380).
    if (const auto python = py::module_::import("sys").attr("executable"); !python.is_none() && py::len(python) > 0)
        fit_host.worker = {python.cast<std::string>(), "-I", "-c", "import importlib.machinery as i, importlib.util as u; "
            "s = u.spec_from_loader('epcsaft._core', i.ExtensionFileLoader('epcsaft._core', '/proc/self/fd/3')); "
            "m = u.module_from_spec(s); s.loader.exec_module(m); m._serve_fit()"};
    fit_host.interrupted = [] { return PyErr_CheckSignals() != 0; };
    py::register_exception_translator([](std::exception_ptr error) {
        try { if (error) std::rethrow_exception(error); } catch (const FitInterrupted&) {}  // the poll raised KeyboardInterrupt
    });
    module.def("_serve_fit", [] { serve_fit(0); });
    module.def("formula_composition", [](const std::vector<double>& solvent, const Parameters& parameters,
        const std::array<std::size_t, 2>& ions, const std::array<int, 2>& stoichiometry, double molality) {
        core::Scope scope("formula_composition");
        return formula_composition(solvent, parameters, ions, stoichiometry, molality);
    });
    py::enum_<ObservationKind>(module, "ObservationKind")
        EPCSAFT_VALUES(ObservationKind, MeanIonicActivityCoefficient, OsmoticCoefficient, SolvationGibbsEnergy,
            TransferGibbsEnergy, BulkRelativePermittivity, IsobaricHeatCapacity, EquilibriumIsobaricHeatCapacity,
            LogEquilibriumConstant, Density, PhaseDensity, EquilibriumPressure, EquilibriumTemperature,
            PhaseMoleFraction, PhaseFraction, PartialPressure, PhaseAmountRatio, AdditionEnthalpy);
    py::enum_<ResidualForm>(module, "ResidualForm")
        EPCSAFT_VALUES(ResidualForm, Difference, LogRatio);
    py::enum_<ObservationStatus>(module, "ObservationStatus")
        EPCSAFT_VALUES(ObservationStatus, NotEvaluated, ReferenceUnavailable, Available, ThermochemistryUnavailable,
            DensityUnavailable, CentralSolveUnavailable, TopologyEvent, MissingMolarMass, NonpositiveLogArgument,
            ActionUnavailable, Nonfinite, InadmissibleTrial);
    py::class_<Observation>(module, "Observation")
        .def(py::init<>())
        EPCSAFT_READWRITE(Observation, kind, temperature, pressure, composition, ion_indices, stoichiometry, molality,
            reference_composition, branch, problem, final_problem, feed_amount, reaction, phase, species,
            denominator_phase, denominator_species, mass_basis, observed, form, scale);
    py::enum_<CheckMetric>(module, "CheckMetric")
        EPCSAFT_VALUES(CheckMetric, AbsoluteDifference, AbsoluteLog10Ratio);
    py::class_<HydrationHold>(module, "HydrationHold")
        .def(py::init<>())
        EPCSAFT_READWRITE(HydrationHold, reference, lower, upper, energy_tolerance, diameter_tolerance,
            maximum_iterations);
    py::class_<HoldResult>(module, "HoldResult")
        EPCSAFT_READONLY(HoldResult, energy_error, correction, uncertainty, derivative, parameter_derivative,
            iterations);
    py::class_<EvaluatedCoordinates>(module, "EvaluatedCoordinates")
        EPCSAFT_READONLY(EvaluatedCoordinates, parameters, active, physical, map, holds, reaction_laws);
    module.def("reaction_log10_reciprocal", &reaction_log10_reciprocal, py::arg("law"));
    py::class_<CheckGroup>(module, "CheckGroup")
        .def(py::init<>())
        EPCSAFT_READWRITE(CheckGroup, name, rows, metric, maximum, mean, baseline, report_only);
    py::class_<CheckResult>(module, "CheckResult")
        EPCSAFT_READONLY(CheckResult, name, unit, status, rows, failed_rows, errors, slopes, maximum, mean, argmax,
            maximum_threshold, mean_threshold, baseline_maximum, baseline_mean, report_only);
    py::class_<Prediction>(module, "Prediction")
        EPCSAFT_READONLY(Prediction, value, unit, status, message, residual, convergence, branch, stable_root_count,
            requested_tolerance_met);
    py::class_<PredictionGradient>(module, "PredictionGradient")
        EPCSAFT_READONLY(PredictionGradient, value, status, message);
    py::enum_<ReactionCoefficient>(module, "ReactionCoefficient")
        EPCSAFT_VALUES(ReactionCoefficient, A, B, C, D, EnthalpyShift);
    py::class_<ReactionCoordinate>(module, "ReactionCoordinate")
        .def(py::init<>())
        EPCSAFT_READWRITE(ReactionCoordinate, name, coefficient, pivot_temperature);
    py::class_<CoordinateTarget>(module, "CoordinateTarget")
        .def(py::init<>())
        EPCSAFT_READWRITE(CoordinateTarget, family, indices, coefficient);
    py::class_<Coordinate>(module, "Coordinate")
        .def(py::init<>())
        EPCSAFT_READWRITE(Coordinate, targets, reaction, origin, scale, lower, upper);
    py::class_<CoordinateSet>(module, "CoordinateSet")
        .def(py::init<>())
        .def_readwrite("coordinates", &CoordinateSet::coordinates)
        .def("active", [](const CoordinateSet& set) { return set.active_parameters(); })
        .def("tie", [](const CoordinateSet& set) { return set.tie(); })
        .def("start", &CoordinateSet::start, py::arg("parameters"), py::arg("observations") = std::vector<Observation>{},
             py::arg("thermochemistry") = std::optional<ThermochemistryRecord>{});
    // The two budgets are scientific inputs (ALG-FIT); zero tolerances select the policy defaults.
    py::class_<FitControls>(module, "FitControls")
        .def(py::init([](int maximum_iterations, double maximum_elapsed_time_seconds, double function_tolerance,
                         double gradient_tolerance, double parameter_tolerance, int workers) {
            return FitControls{function_tolerance, gradient_tolerance, parameter_tolerance, maximum_iterations, maximum_elapsed_time_seconds, workers};
        }), py::kw_only(), py::arg("maximum_iterations"), py::arg("maximum_elapsed_time_seconds"),
            py::arg("function_tolerance") = 0.0, py::arg("gradient_tolerance") = 0.0, py::arg("parameter_tolerance") = 0.0,
            py::arg("workers") = 0)
        EPCSAFT_READWRITE(FitControls, function_tolerance, gradient_tolerance, parameter_tolerance, maximum_iterations,
            maximum_elapsed_time_seconds, workers);
    py::class_<TrialFailure>(module, "TrialFailure")
        EPCSAFT_READONLY(TrialFailure, trial, observation, optimizer, status, message, jacobian);
    py::class_<FitCovariance>(module, "FitCovariance")
        EPCSAFT_READONLY(FitCovariance, status, assumption, rank, singular_values, left_vectors, right_vectors,
            projected_gradient_norm, objective_rows, active_bounds, optimizer, physical, variance_factor);
    py::class_<MinimaxRequest>(module, "MinimaxRequest").def(py::init<>())
        EPCSAFT_READWRITE(MinimaxRequest, targets, group_scales, q, tolerance, activity_distance, margin);
    py::class_<MinimaxRow>(module, "MinimaxRow")
        EPCSAFT_READONLY(MinimaxRow, kind, row, group, metric, sign, threshold, scale, auxiliary, members);
    py::class_<MinimaxResult>(module, "MinimaxResult")
        EPCSAFT_READONLY(MinimaxResult, availability, application_status, solver_status, feasible, stationary, t,
            target_maximum, gap, primal, stationarity, dual_sign, complementarity, point, lower, upper, values,
            physical_values, jacobian, multipliers, physical_multipliers, z_lower, z_upper, physical_z_lower,
            physical_z_upper, lower_distances, upper_distances, hold_multipliers, full_stationarity, root_errors,
            slack_groups, slack_rows, active_lower, active_upper, hold_causes, rows, frozen_checks);
    py::class_<FitResult>(module, "FitResult")
        EPCSAFT_READONLY(FitResult, minimax, status, chart, checks, acceptance, optimizer, message, usable, physical,
            parameters, predictions, residuals, weighted_residuals, physical_jacobian, optimizer_jacobian, statuses,
            failures, iterations, residual_evaluations, jacobian_evaluations, trial_failures, iteration_seconds,
            iteration_costs, initial_cost, final_cost, covariance);
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
                         const std::optional<std::vector<std::vector<double>>>& starts, std::optional<int> iterations, std::optional<double> seconds, std::optional<MinimaxRequest> minimax) -> py::object {
        core::Scope scope("fit");
        FitRequest request{parameters, std::move(thermochemistry), {std::move(coordinates)}, std::move(observations),
                    std::move(weights), controls, std::move(holds), std::move(checks), std::move(baseline), singular_vectors, std::move(minimax)};
        return starts ? py::cast(fit_sequence(request, *starts, iterations.value_or(controls.maximum_iterations), seconds.value_or(controls.maximum_elapsed_time_seconds))) : py::cast(fit(request));
    }, py::arg("parameters"), py::arg("thermochemistry"), py::arg("coordinates"), py::arg("observations"),
       py::arg("weights"), py::arg("controls"), py::arg("holds") = std::vector<HydrationHold>{},
       py::arg("checks") = std::vector<CheckGroup>{}, py::arg("baseline") = std::nullopt, py::arg("singular_vectors") = false,
       py::arg("starts") = std::nullopt, py::arg("iterations") = std::nullopt, py::arg("seconds") = std::nullopt, py::arg("minimax") = std::nullopt);
    module.def("observations_at", [](const Parameters& parameters, std::vector<Coordinate> coordinates,
        const std::vector<double>& optimizer, std::vector<HydrationHold> holds, std::vector<Observation> observations,
        const std::vector<double>& errors, std::optional<ThermochemistryRecord> thermochemistry) {
        return observations_at(point_request(parameters, std::move(coordinates), std::move(holds), std::move(observations),
            std::move(thermochemistry)), optimizer, errors);
    }, py::arg("parameters"), py::arg("coordinates"), py::arg("optimizer"), py::arg("holds"), py::arg("observations"),
       py::arg("errors"), py::arg("thermochemistry") = std::nullopt);
    py::class_<ObservationEvaluation>(module, "ObservationEvaluation")
        EPCSAFT_READONLY(ObservationEvaluation, predictions, gradient, optimizer_gradient, chart,
            residual_optimizer_gradient, statuses);
    module.def("evaluate_observations", [](const Parameters& parameters, std::vector<Coordinate> coordinates,
        const std::vector<double>& optimizer, std::vector<HydrationHold> holds, std::vector<Observation> observations,
        std::optional<ThermochemistryRecord> thermochemistry) {
        core::Scope scope("evaluate_observations");
        const auto request = point_request(parameters, std::move(coordinates), std::move(holds), std::move(observations),
            std::move(thermochemistry));
        std::vector<std::size_t> rows(request.observations.size()); std::iota(rows.begin(), rows.end(), 0);
        FitResult failures; ObservationEvaluator evaluator(request, failures, rows); ObservationEvaluation out;
        if (!evaluator.evaluate(optimizer, true, out)) throw std::invalid_argument(failures.failures.back().message);
        return out;
    }, py::arg("parameters"), py::arg("coordinates"), py::arg("optimizer"), py::arg("holds"), py::arg("observations"),
       py::arg("thermochemistry") = std::nullopt);
    module.def("evaluate_minimax", [](const Parameters& parameters, std::vector<Coordinate> coordinates,
        const std::vector<double>& optimizer, std::vector<HydrationHold> holds, std::vector<Observation> observations,
        std::vector<CheckGroup> checks, std::optional<Parameters> baseline, const MinimaxRequest& minimax,
        std::optional<ThermochemistryRecord> thermochemistry) {
        core::Scope scope("evaluate_minimax");
        auto request = point_request(parameters, std::move(coordinates), std::move(holds), std::move(observations),
            std::move(thermochemistry));
        request.checks = std::move(checks); request.baseline = std::move(baseline); request.minimax = minimax;
        return fit_minimax(request, optimizer);
    }, py::arg("parameters"), py::arg("coordinates"), py::arg("optimizer"), py::arg("holds"), py::arg("observations"),
       py::arg("checks"), py::arg("baseline"), py::arg("minimax"), py::arg("thermochemistry") = std::nullopt);
    module.def("evaluate_coordinates", &evaluate_coordinates, py::arg("parameters"), py::arg("set"), py::arg("optimizer"),
        py::arg("holds"), py::arg("thermochemistry") = std::nullopt);
    module.def("evaluate_checks", [](const Parameters& stored, std::vector<Observation> observations,
        std::vector<CheckGroup> checks, std::optional<Parameters> baseline) {
        FitRequest request; request.observations = std::move(observations); request.checks = std::move(checks); request.baseline = std::move(baseline);
        return evaluate_checks(request, stored);
    }, py::arg("stored"), py::arg("observations"), py::arg("checks"), py::arg("baseline") = std::nullopt);
}

}  // namespace epcsaft::regression
