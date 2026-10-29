#include "regression/fit.hpp"

#include "core/trace.hpp"

#include "fields.hpp"
#include "solver_options.hpp"

#include <nanobind/nanobind.h>
#include <nanobind/stl/array.h>
#include <nanobind/stl/map.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <numeric>

namespace nb = nanobind;
namespace epcsaft::regression {
namespace {

thread_local std::optional<nb::python_error> interrupted_error;

template <typename Function>
auto fit_entry(Function&& function) {
    interrupted_error.reset();
    try {
        auto result = function();
        interrupted_error.reset();
        return result;
    } catch (const FitInterrupted&) {
        throw;  // The translator consumes the saved exception after native cleanup.
    } catch (...) {
        interrupted_error.reset();
        throw;
    }
}

// The point readers' request: fit's model inputs at a declared point, without optimizer controls.
FitRequest point_request(const Parameters& parameters, std::vector<Coordinate> coordinates, std::vector<HydrationHold> holds,
                         std::vector<Observation> observations, std::optional<ThermochemistryRecord> thermochemistry, const nb::object& solver_options) {
    FitRequest request; request.controls.solver_options = bindings::solver_options(solver_options, 7).resolved({}); request.parameters = parameters; request.thermochemistry = std::move(thermochemistry);
    request.coordinates.coordinates = std::move(coordinates); request.holds = std::move(holds); request.observations = std::move(observations);
    return request;
}

}  // namespace

void register_regression(nb::module_& module) {
    FitHost host;
    if (const nb::object python = nb::module_::import_("sys").attr("executable"); !python.is_none() && nb::len(python) > 0)
        // The native library loads by its installed path, so its $ORIGIN search path finds bundled libraries; a path that no
        // longer names the held file (an editable relink replaced it) refuses the worker.
        host.worker = {nb::cast<std::string>(python), "-I", "-c", "import ctypes, os, importlib.machinery as i, importlib.util as u; "
            "p = os.readlink('/proc/self/fd/4'); s, f = os.stat(p), os.fstat(4); "
            "(s.st_dev, s.st_ino) == (f.st_dev, f.st_ino) or os._exit(1); "
            "native = ctypes.CDLL(p, mode=ctypes.RTLD_GLOBAL); "
            "s = u.spec_from_loader('epcsaft._core', i.ExtensionFileLoader('epcsaft._core', '/proc/self/fd/3')); "
            "m = u.module_from_spec(s); s.loader.exec_module(m); m._serve_fit()"};
    host.interrupted = [] {
        if (PyErr_CheckSignals() == 0) return false;
        interrupted_error.emplace();  // Capture and clear before native worker cleanup.
        return true;
    };
    set_fit_host(std::move(host), reinterpret_cast<const void*>(&register_regression));
    nb::register_exception_translator([](const std::exception_ptr& error, void*) {
        try { if (error) std::rethrow_exception(error); } catch (const FitInterrupted&) {
            if (interrupted_error) {
                interrupted_error->restore();
                interrupted_error.reset();
            } else PyErr_SetString(PyExc_RuntimeError, "fit interruption has no saved Python exception");
        }
    });
    module.def("_serve_fit", [] { serve_fit(0); });
    module.def("formula_composition", [](const std::vector<double>& solvent, const Parameters& parameters,
        const std::array<std::size_t, 2>& ions, const std::array<int, 2>& stoichiometry, double molality) {
        core::Scope scope("formula_composition");
        return formula_composition(solvent, parameters, ions, stoichiometry, molality);
    });
    nb::enum_<ObservationKind>(module, "ObservationKind")
        EPCSAFT_VALUES(ObservationKind, MeanIonicActivityCoefficient, OsmoticCoefficient, MixtureOsmoticCoefficient, SolvationGibbsEnergy,
            TransferGibbsEnergy, BulkRelativePermittivity, IsobaricHeatCapacity, EquilibriumIsobaricHeatCapacity,
            LogEquilibriumConstant, Density, PhaseDensity, EquilibriumPressure, EquilibriumTemperature,
            PhaseMoleFraction, PhaseFraction, PartialPressure, PhaseAmountRatio, AdditionEnthalpy, SaturationCappedAdditionEnthalpy);
    nb::enum_<ResidualForm>(module, "ResidualForm")
        EPCSAFT_VALUES(ResidualForm, Difference, LogRatio);
    nb::enum_<ObservationStatus>(module, "ObservationStatus")
        EPCSAFT_VALUES(ObservationStatus, NotEvaluated, ReferenceUnavailable, Available, ThermochemistryUnavailable,
            DensityUnavailable, CentralSolveUnavailable, TopologyEvent, MissingMolarMass, NonpositiveLogArgument,
            ActionUnavailable, Nonfinite, InadmissibleTrial);
    nb::class_<Observation>(module, "Observation")
        .def(nb::init<>())
        EPCSAFT_READWRITE(Observation, kind, temperature, pressure, composition, ion_indices, stoichiometry, molality, finite_solute_molalities,
            reference_composition, branch, problem, final_problem, feed_amount, basis_amount, loading_cap, loading_bounds, loading_grid, reaction, phase, species,
            denominator_phase, denominator_species, mass_basis, observed, form, scale);
    nb::enum_<CheckMetric>(module, "CheckMetric")
        EPCSAFT_VALUES(CheckMetric, AbsoluteDifference, AbsoluteLog10Ratio);
    nb::class_<HydrationHold>(module, "HydrationHold")
        .def(nb::init<>())
        EPCSAFT_READWRITE(HydrationHold, reference, lower, upper, energy_tolerance, diameter_tolerance,
            maximum_iterations);
    nb::class_<HoldResult>(module, "HoldResult")
        EPCSAFT_READONLY(HoldResult, energy_error, correction, uncertainty, derivative, parameter_derivative,
            iterations, rank);
    nb::class_<EvaluatedCoordinates>(module, "EvaluatedCoordinates")
        EPCSAFT_READONLY(EvaluatedCoordinates, parameters, active, physical, map, holds, reaction_laws);
    module.def("reaction_log10_reciprocal", &reaction_log10_reciprocal, nb::arg("law"));
    nb::class_<CheckGroup>(module, "CheckGroup")
        .def(nb::init<>())
        EPCSAFT_READWRITE(CheckGroup, name, rows, metric, maximum, mean, baseline, report_only);
    nb::class_<CheckResult>(module, "CheckResult")
        EPCSAFT_READONLY(CheckResult, name, unit, status, rows, failed_rows, errors, slopes, maximum, mean, argmax,
            maximum_threshold, mean_threshold, baseline_maximum, baseline_mean, report_only);
    nb::class_<equilibrium::AffineLoading>(module, "AffineLoading")
        EPCSAFT_READONLY(equilibrium::AffineLoading, initial, direction, lower, upper, seed, pressure);
    nb::class_<equilibrium::SaturationOperand>(module, "SaturationOperand")
        EPCSAFT_READONLY(equilibrium::SaturationOperand, quantity, value, error);
    nb::class_<equilibrium::SaturationCandidate>(module, "SaturationCandidate")
        EPCSAFT_READONLY(equilibrium::SaturationCandidate, state, operands);
    nb::class_<equilibrium::SaturationComparison>(module, "SaturationComparison")
        EPCSAFT_READONLY(equilibrium::SaturationComparison, left, right, identity, operands);
    nb::class_<equilibrium::SaturationNode>(module, "SaturationNode")
        EPCSAFT_READONLY(equilibrium::SaturationNode, loading, mismatch, slope, state, estimates, action_cause, action_seconds, estimate_seconds, local_estimates, action_diagnostics);
    nb::class_<equilibrium::SaturationSearch>(module, "SaturationSearch")
        EPCSAFT_READONLY(equilibrium::SaturationSearch, grid, nodes, events, candidates, comparisons, causes, groups,
            accepted, first_numerical_failure, brackets, iterations, evaluations, hessian_calls, attempts, complete, declaration, work_at_entry, uniqueness_scope, detection_limits);
    nb::class_<Prediction>(module, "Prediction")
        EPCSAFT_SOLVER_OPTIONS(Prediction, 7)
        EPCSAFT_READONLY(Prediction, value, unit, status, message, residual, convergence, branch, stable_root_count,
            requested_tolerance_met, first_numerical_failure, saturation, saturation_loading, endpoint_loading, added_amount, cap_branch, saturation_loading_gradient, endpoint_loading_gradient);
    nb::class_<PredictionGradient>(module, "PredictionGradient")
        EPCSAFT_READONLY(PredictionGradient, value, status, message, parameter_actions, first_numerical_failure, saturation_loading_gradient, endpoint_loading_gradient);
    nb::enum_<ReactionCoefficient>(module, "ReactionCoefficient")
        EPCSAFT_VALUES(ReactionCoefficient, A, B, C, D, EnthalpyShift);
    nb::class_<ReactionCoordinate>(module, "ReactionCoordinate")
        .def(nb::init<>())
        EPCSAFT_READWRITE(ReactionCoordinate, name, coefficient, pivot_temperature);
    nb::class_<CoordinateTarget>(module, "CoordinateTarget")
        .def(nb::init<>())
        EPCSAFT_READWRITE(CoordinateTarget, family, indices, coefficient);
    nb::class_<Coordinate>(module, "Coordinate")
        .def(nb::init<>())
        EPCSAFT_READWRITE(Coordinate, targets, reaction, origin, scale, lower, upper);
    nb::class_<CoordinateSet>(module, "CoordinateSet")
        .def(nb::init<>())
        .def_rw("coordinates", &CoordinateSet::coordinates)
        .def("active", [](const CoordinateSet& set) { return set.active_parameters(); })
        .def("tie", [](const CoordinateSet& set) { return set.tie(); })
        .def("start", &CoordinateSet::start, nb::arg("parameters"), nb::arg("observations") = std::vector<Observation>{},
             nb::arg("thermochemistry") = std::optional<ThermochemistryRecord>{});
    // The two budgets are scientific inputs (ALG-FIT); zero tolerances select the policy defaults.
    auto controls = nb::class_<FitControls>(module, "FitControls")
        .def("__init__", [](FitControls* self, int maximum_iterations, double maximum_elapsed_time_seconds, nb::object function_tolerance,
                         nb::object gradient_tolerance, nb::object parameter_tolerance, int workers, nb::object solver_options) {
            FitControls c{bindings::solver_options(solver_options), maximum_iterations, maximum_elapsed_time_seconds, workers};
            for (const auto& [i, value] : {std::pair{policy::function, function_tolerance}, {policy::gradient, gradient_tolerance}, {policy::parameter, parameter_tolerance}})
                bindings::solver_option(c.solver_options, policy::SOLVER_POLICIES[i]->name, value, 63, true);
            new (self) FitControls(std::move(c));
        }, nb::kw_only(), nb::arg("maximum_iterations"), nb::arg("maximum_elapsed_time_seconds"),
            nb::arg("function_tolerance").none() = 0.0, nb::arg("gradient_tolerance").none() = 0.0, nb::arg("parameter_tolerance").none() = 0.0,
            nb::arg("workers") = 0, nb::arg("solver_options") = nb::none())
        EPCSAFT_SOLVER_OPTIONS(FitControls, 63)
        EPCSAFT_READWRITE(FitControls, maximum_iterations, maximum_elapsed_time_seconds, workers);
    for (const auto& [name, i] : {std::pair{"function_tolerance", policy::function}, {"gradient_tolerance", policy::gradient}, {"parameter_tolerance", policy::parameter}})
        controls.def_prop_rw(name, [i](const FitControls& c) { return c.solver_options.value(i); }, [i](FitControls& c, nb::object v) {
            bindings::solver_option(c.solver_options, policy::SOLVER_POLICIES[i]->name, v, 63, true);
        }, nb::for_setter(nb::arg("v").none()));
    nb::class_<TrialFailure>(module, "TrialFailure")
        EPCSAFT_READONLY(TrialFailure, trial, observation, optimizer, status, message, jacobian, parameter_actions, first_numerical_failure);
    nb::class_<FitCovariance>(module, "FitCovariance")
        EPCSAFT_READONLY(FitCovariance, status, assumption, rank, singular_values, left_vectors, right_vectors,
            projected_gradient_norm, objective_rows, active_bounds, optimizer, physical, variance_factor);
    nb::enum_<MeasuredInputKind>(module, "MeasuredInputKind")
        EPCSAFT_VALUES(MeasuredInputKind, Temperature, Pressure, MoleFraction, Molality);
    nb::class_<MeasuredInput>(module, "MeasuredInput").def(nb::init<>())
        EPCSAFT_READWRITE(MeasuredInput, kind, species, sigma, lower, upper, unit);
    nb::class_<MeasurementGroup>(module, "MeasurementGroup").def(nb::init<>())
        EPCSAFT_READWRITE(MeasurementGroup, rows, inputs, dependent_species, feed_total_mol);
    nb::class_<MeasurementUncertainty>(module, "MeasurementUncertainty").def(nb::init<>())
        EPCSAFT_READWRITE(MeasurementUncertainty, groups, output_sigma, covariance, correlation);
    nb::class_<InputReconciliation>(module, "InputReconciliation")
        EPCSAFT_READONLY(InputReconciliation, values, corrections, penalties, jacobian, row_groups, column_groups,
            column_inputs, output_cost, input_cost);
    nb::class_<MinimaxRequest>(module, "MinimaxRequest").def(nb::init<>())
        EPCSAFT_READWRITE(MinimaxRequest, targets, group_scales, q, tolerance, activity_distance, margin);
    nb::class_<MinimaxRow>(module, "MinimaxRow")
        EPCSAFT_READONLY(MinimaxRow, kind, row, group, metric, sign, threshold, scale, auxiliary, members);
    nb::class_<MinimaxResult>(module, "MinimaxResult")
        EPCSAFT_READONLY(MinimaxResult, availability, application_status, solver_status, linear_solver, feasible, stationary, t,
            target_maximum, gap, primal, stationarity, dual_sign, complementarity, point, lower, upper, values,
            physical_values, jacobian, multipliers, physical_multipliers, z_lower, z_upper, physical_z_lower,
            physical_z_upper, lower_distances, upper_distances, hold_multipliers, full_stationarity, root_errors,
            slack_groups, slack_rows, active_lower, active_upper, hold_causes, rows, frozen_checks, ipopt_diagnostics);
    nb::class_<FitResult>(module, "FitResult")
        EPCSAFT_SOLVER_OPTIONS(FitResult, 63)
        EPCSAFT_READONLY(FitResult, saturation_predictions, minimax, status, chart, checks, acceptance, optimizer, message, usable, physical,
            parameters, predictions, residuals, weighted_residuals, physical_jacobian, optimizer_jacobian, statuses,
            failures, iterations, residual_evaluations, jacobian_evaluations, trial_failures, iteration_seconds,
            iteration_costs, initial_cost, final_cost, covariance, reconciliation, rank_records, first_numerical_failure);
    module.def("predict", [](const HelmholtzModel& model, const Observation& observation) {
        core::Scope scope("predict");
        validate(observation, model);
        return predict(model, observation);
    }, nb::arg("model"), nb::arg("observation"));
    module.def("differentiate", [](const HelmholtzModel& model, const Observation& observation) {
        core::Scope scope("differentiate");
        validate(observation, model);
        return differentiate(model, observation, predict(model, observation));
    }, nb::arg("model"), nb::arg("observation"));
    module.def("fit", [](const Parameters& parameters, std::optional<ThermochemistryRecord> thermochemistry,
                         std::vector<Coordinate> coordinates, std::vector<Observation> observations,
                         std::vector<double> weights, const FitControls& controls, std::vector<HydrationHold> holds,
                         std::vector<CheckGroup> checks, std::optional<Parameters> baseline, bool singular_vectors,
                         const std::optional<std::vector<std::vector<double>>>& starts, std::optional<int> iterations, std::optional<double> seconds,
                         std::optional<MinimaxRequest> minimax, std::optional<MeasurementUncertainty> uncertainty) -> nb::object {
        return fit_entry([&]() -> nb::object {
            core::Scope scope("fit");
            FitRequest request{parameters, std::move(thermochemistry), {std::move(coordinates)}, std::move(observations),
                        std::move(weights), controls, std::move(holds), std::move(checks), std::move(baseline), singular_vectors,
                        std::move(minimax), std::move(uncertainty)};
            return starts ? nb::cast(fit_sequence(request, *starts, iterations.value_or(controls.maximum_iterations), seconds.value_or(controls.maximum_elapsed_time_seconds))) : nb::cast(fit(request));
        });
    }, nb::arg("parameters"), nb::arg("thermochemistry"), nb::arg("coordinates"), nb::arg("observations"),
       nb::arg("weights"), nb::arg("controls"), nb::arg("holds") = std::vector<HydrationHold>{},
       nb::arg("checks") = std::vector<CheckGroup>{}, nb::arg("baseline") = std::nullopt, nb::arg("singular_vectors") = false,
       nb::arg("starts") = std::nullopt, nb::arg("iterations") = std::nullopt, nb::arg("seconds") = std::nullopt,
       nb::arg("minimax") = std::nullopt, nb::arg("uncertainty") = std::nullopt);
    module.def("observations_at", [](const Parameters& parameters, std::vector<Coordinate> coordinates,
        const std::vector<double>& optimizer, std::vector<HydrationHold> holds, std::vector<Observation> observations,
        const std::vector<double>& errors, std::optional<ThermochemistryRecord> thermochemistry, nb::object solver_options) {
        return observations_at(point_request(parameters, std::move(coordinates), std::move(holds), std::move(observations),
            std::move(thermochemistry), solver_options), optimizer, errors);
    }, nb::arg("parameters"), nb::arg("coordinates"), nb::arg("optimizer"), nb::arg("holds"), nb::arg("observations"),
       nb::arg("errors"), nb::arg("thermochemistry") = std::nullopt, nb::arg("solver_options") = nb::none());
    nb::class_<ObservationEvaluation>(module, "ObservationEvaluation")
        EPCSAFT_SOLVER_OPTIONS(ObservationEvaluation, 63)
        EPCSAFT_READONLY(ObservationEvaluation, saturation_predictions, predictions, gradient, optimizer_gradient, chart,
            residual_optimizer_gradient, statuses, reconciliation, state_densities, reference_values, first_numerical_failure, first_numerical_failure_observation);
    module.def("evaluate_observations", [](const Parameters& parameters, std::vector<Coordinate> coordinates,
        const std::vector<double>& optimizer, std::vector<HydrationHold> holds, std::vector<Observation> observations,
        std::optional<ThermochemistryRecord> thermochemistry, std::optional<MeasurementUncertainty> uncertainty, nb::object solver_options) {
        core::Scope scope("evaluate_observations");
        auto request = point_request(parameters, std::move(coordinates), std::move(holds), std::move(observations),
            std::move(thermochemistry), solver_options);
        request.uncertainty = std::move(uncertainty);
        std::vector<std::size_t> rows(request.observations.size()); std::iota(rows.begin(), rows.end(), 0);
        FitResult failures; ObservationEvaluator evaluator(request, failures, rows); ObservationEvaluation out;
        if (!evaluator.evaluate(optimizer, true, out)) throw std::invalid_argument(failures.failures.back().message);
        return out;
    }, nb::arg("parameters"), nb::arg("coordinates"), nb::arg("optimizer"), nb::arg("holds"), nb::arg("observations"),
       nb::arg("thermochemistry") = std::nullopt, nb::arg("uncertainty") = std::nullopt, nb::arg("solver_options") = nb::none());
    module.def("evaluate_minimax", [](const Parameters& parameters, std::vector<Coordinate> coordinates,
        const std::vector<double>& optimizer, std::vector<HydrationHold> holds, std::vector<Observation> observations,
        std::vector<CheckGroup> checks, std::optional<Parameters> baseline, const MinimaxRequest& minimax,
        std::optional<ThermochemistryRecord> thermochemistry, nb::object solver_options) {
        return fit_entry([&] {
            core::Scope scope("evaluate_minimax");
            auto request = point_request(parameters, std::move(coordinates), std::move(holds), std::move(observations),
                std::move(thermochemistry), solver_options);
            request.checks = std::move(checks); request.baseline = std::move(baseline); request.minimax = minimax;
            return fit_minimax(request, optimizer);
        });
    }, nb::arg("parameters"), nb::arg("coordinates"), nb::arg("optimizer"), nb::arg("holds"), nb::arg("observations"),
       nb::arg("checks"), nb::arg("baseline"), nb::arg("minimax"), nb::arg("thermochemistry") = std::nullopt, nb::arg("solver_options") = nb::none());
    module.def("evaluate_coordinates", [](const Parameters& p, const CoordinateSet& c, const std::vector<double>& u,
        const std::vector<HydrationHold>& h, const std::optional<ThermochemistryRecord>& t, nb::object options) {
        return evaluate_coordinates(p, c, u, h, t, bindings::solver_options(options, 7).resolved({}));
    }, nb::arg("parameters"), nb::arg("set"), nb::arg("optimizer"),
        nb::arg("holds"), nb::arg("thermochemistry") = std::nullopt, nb::arg("solver_options") = nb::none());
    module.def("evaluate_checks", [](const Parameters& stored, std::vector<Observation> observations,
        std::vector<CheckGroup> checks, std::optional<Parameters> baseline, nb::object solver_options) {
        FitRequest request; request.controls.solver_options = bindings::solver_options(solver_options, 7).resolved({}); request.observations = std::move(observations); request.checks = std::move(checks); request.baseline = std::move(baseline);
        return evaluate_checks(request, stored);
    }, nb::arg("stored"), nb::arg("observations"), nb::arg("checks"), nb::arg("baseline") = std::nullopt, nb::arg("solver_options") = nb::none());
}

}  // namespace epcsaft::regression
