#include "equilibrium/problem.hpp"
#include "equilibrium/gibbs_minimizer.hpp"
#include "equilibrium_record.hpp"

#include "fields.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <optional>

namespace py = pybind11;
namespace epcsaft::equilibrium {

void register_equilibrium(py::module_& module) {
    module.def("property_asl_record", [](py::object model) {
        const py::object native_owner = model.attr("_model"), document = model.attr("parameters").attr("_document");
        const auto& native = native_owner.cast<const HelmholtzModel&>();
        if (document.is_none() || !native.active_parameters().empty() || !native.thermochemistry())
            throw std::invalid_argument("ASL property record requires document parameters, inactive and with thermochemistry");
        return asl_property_record(document.cast<ParameterDocument>(), *native.thermochemistry());
    }, py::arg("model"));
    module.attr("EOS_STANDARD_STATE_ID") = EOS_STANDARD_STATE_ID;
    py::enum_<PhaseKind>(module, "PhaseKind")
        EPCSAFT_VALUES(PhaseKind, Liquid, Vapor);
    py::enum_<AmountKind>(module, "AmountKind")
        EPCSAFT_VALUES(AmountKind, Free, Pinned);
    py::enum_<FeedKind>(module, "FeedKind")
        EPCSAFT_VALUES(FeedKind, MoleFractions, Amounts, MassFractions, Molalities);
    py::class_<AmountSpec>(module, "AmountSpec")
        .def(py::init<>())
        EPCSAFT_READWRITE(AmountSpec, kind, fraction)
        .def_static("free", &AmountSpec::free)
        .def_static("pinned", &AmountSpec::pinned);
    py::class_<FreeValue>(module, "FreeValue")
        .def(py::init<>())
        EPCSAFT_READWRITE(FreeValue, free, value)
        .def_static("fixed", &FreeValue::fixed)
        .def_static("variable", [](std::optional<double> guess) {
            return guess.has_value() ? FreeValue::variable(*guess) : FreeValue::variable();
        }, py::arg("guess") = py::none());
    py::class_<PhaseSpec>(module, "PhaseSpec")
        .def(py::init<>())
        EPCSAFT_READWRITE(PhaseSpec, name, kind, support, support_declared, amount, composition_guess, packing_guess);
    py::class_<FeedSpec>(module, "FeedSpec")
        .def(py::init<>())
        EPCSAFT_READWRITE(FeedSpec, kind, values, solvent, solvent_mass_kg);
    py::enum_<SpeciesReference>(module, "SpeciesReference")
        EPCSAFT_VALUES(SpeciesReference, InfiniteDilution, PureLiquid, PureGas);
    py::class_<NeutralReferenceDeclaration>(module, "NeutralReferenceDeclaration")
        .def(py::init<>())
        EPCSAFT_READWRITE(NeutralReferenceDeclaration, solvent_index, solute_molality_ratios,
            standard_molality_mol_per_kg, species_reference);
    py::enum_<ReferenceSourceBasis>(module, "ReferenceSourceBasis")
        EPCSAFT_VALUES(ReferenceSourceBasis, CommonMolalityInfiniteDilution, RawMoleFractionInfiniteDilution);
    py::enum_<NeutralReferenceStatus>(module, "NeutralReferenceStatus")
        EPCSAFT_VALUES(NeutralReferenceStatus, Available, InvalidDeclaration, NonNeutralDeclaration,
            UnsupportedSourceBasis, FreeTemperature, FreePressure, DensityUnavailable, BranchChange, Nonfinite,
            TerminalChangeExceeded);
    py::class_<NeutralReferenceReactionDiagnostic>(module, "NeutralReferenceReactionDiagnostic")
        EPCSAFT_READONLY(NeutralReferenceReactionDiagnostic, reaction_index, terminal_contraction,
            observed_terminal_change, terminal_contraction_pressure_derivative,
            terminal_contraction_pressure_second_derivative, reference_pressure_pa, status, message);
    py::class_<ReactionReference>(module, "ReactionReference")
        .def(py::init([](ReferenceSourceBasis basis, std::optional<double> pressure) {
            return ReactionReference{pressure, basis};
        }), py::arg("source_basis"), py::arg("reference_pressure_pa") = py::none())
        EPCSAFT_READWRITE(ReactionReference, source_basis, reference_pressure_pa);
    py::class_<ReactionSpec>(module, "ReactionSpec")
        .def(py::init<>())
        EPCSAFT_READWRITE(ReactionSpec, name, stoichiometry, equilibrium_constant, affinity_j_per_mol, correlation,
            reference);
    py::class_<StabilityOptions>(module, "StabilityOptions")
        .def(py::init<>())
        EPCSAFT_READWRITE(StabilityOptions, evaluate_supplied_trials, trials, declaration_budget, trial_budget,
            standard_uncertainty_over_rt);
    py::class_<CriticalDeclaration>(module, "CriticalDeclaration")
        .def(py::init<>())
        EPCSAFT_READWRITE(CriticalDeclaration, density_guess_mol_per_m3, mode_guess);
    py::enum_<EnvelopeCoordinate>(module, "EnvelopeCoordinate")
        EPCSAFT_VALUES(EnvelopeCoordinate, LogTemperature, LogPressure, LogK);
    py::class_<EnvelopeCondition>(module, "EnvelopeCondition")
        .def(py::init<>())
        EPCSAFT_READWRITE(EnvelopeCondition, kind, component, target);
    py::class_<EquilibriumProblem>(module, "EquilibriumProblem")
        .def(py::init<>())
        EPCSAFT_READWRITE(EquilibriumProblem, phases, temperature, pressure, feed, reactions)
        .def_readwrite("neutral_reference", &EquilibriumProblem::source_reference)
        EPCSAFT_READWRITE(EquilibriumProblem, standard_gibbs_j_per_mol, topology_declared, row_tolerance, enthalpy_j,
            entropy_j_per_k, critical, envelope, discovery);
    py::enum_<RowKind>(module, "RowKind")
        EPCSAFT_VALUES(RowKind, Chemical, Pressure, Moiety, Electro, Normalization, Enthalpy, Entropy, CriticalMode,
            CriticalNormalization, CriticalCubic, EnvelopeCondition);
    py::class_<Row>(module, "Row")
        EPCSAFT_READONLY(Row, kind, phase, index)
        .def_property_readonly("eqid", [](const Row& row) { return row_kind_name(row.kind); });
    py::class_<PointEvaluation>(module, "PointEvaluation")
        .def_readonly("values", &PointEvaluation::values)
        .def_readonly("jacobian", &PointEvaluation::jacobian)
        EPCSAFT_READONLY(PointEvaluation, objective, floor);
    py::class_<CompiledSystem>(module, "CompiledSystem")
        .def_property_readonly("rows", [](const CompiledSystem& system) { return system.rows; })
        EPCSAFT_READONLY(CompiledSystem, pinned_count, free_intensive_count, rank, charge_row_count, primal_size,
            generic_nullity, feed_total, has_ions, has_reaction, conserves_total_moles, jacobian_rank, removed_species)
        .def("generic_point", &CompiledSystem::generic_point)
        .def("evaluate", [](const CompiledSystem& system, const std::vector<double>& point) {
            return system.evaluate(point); })
        .def("jacobian", [](const CompiledSystem& system, const std::vector<double>& point) {
            return system.jacobian(point); })
        .def("lever_jacobian", &CompiledSystem::lever_jacobian)
        .def("objective_gradient", [](const CompiledSystem& system, const std::vector<double>& point) {
            return system.objective_gradient(point); })
        .def("hessian", [](const CompiledSystem& system, const std::vector<double>& point,
                           const std::vector<double>& multipliers) {
            return system.hessian(point, multipliers); },
            py::arg("point"), py::arg("constraint_multipliers") = std::vector<double>{});
    py::class_<EquilibriumResult>(module, "EquilibriumResult", "Phases in declared order; K, Pa, mol, mol/m^3; per-phase-per-component lists are phase-major in component_ids order.")
        .def_readonly("neutral_reference_diagnostics", &EquilibriumResult::reference_diagnostics)
        EPCSAFT_READONLY(EquilibriumResult, discovery, success, message, solver_status, temperature, pressure, amounts,
            mole_fractions, molar_densities, packing_fractions,
            // Dimensionless reduced-energy balance/charge stationarity coordinates of the solution.
            multipliers,
            // Extensive-objective multipliers of the non-chemical, non-pressure `rows` in order, formed at
            // the returned point from `multipliers`; raw optimizer duals only when a phase amount is pinned.
            constraint_multipliers, residuals, classified_residuals)
        .def_readonly("sensitivity", &EquilibriumResult::sensitivity)
        EPCSAFT_READONLY(EquilibriumResult, coordinates, phase_count, iterations, evaluations, hessian_calls, objective,
            reduced_hessian_minimum, convergence, rows, unknown_count, jacobian_rank, generic_nullity,
            requested_tolerance_met, topology_event, topology_phase, topology_pair, critical);
    py::class_<CriticalRecord>(module, "CriticalRecord")
        EPCSAFT_READONLY(CriticalRecord, eigenvalues, mode, cubic, zero_band, condition_estimate, backward_error, rank);
    py::class_<StabilityTrial>(module, "StabilityTrial")
        EPCSAFT_READONLY(StabilityTrial, phase, seed, evidence, message, optimizer, raw_coordinates, raw_residuals,
            composition, split_coordinates, basis, row_transformation, quantities, split_checks, donor,
            resolved_descent);
    py::class_<StabilityReport>(module, "StabilityReport")
        EPCSAFT_READONLY(StabilityReport, completed, message, trials, candidate_coordinates, standard_over_rt,
            plane_over_rt, present_phases, work);
    py::class_<DiscoveryRecord, std::shared_ptr<DiscoveryRecord>>(module, "DiscoveryRecord")
        EPCSAFT_READONLY(DiscoveryRecord, declaration, candidate, rounds, declarations, declaration_roles,
            declaration_attempt_coordinates, work);
    module.def("cooling_onset", &cooling_onset, py::arg("compiled"));
    module.def("stability_test", &stability_test, py::arg("compiled"), py::arg("central"),
        py::arg("options") = StabilityOptions{});
    py::enum_<TopologyEvent>(module, "TopologyEvent")
        .value("None_", TopologyEvent::None)
        EPCSAFT_VALUES(TopologyEvent, PhaseDisappeared, PhasesCoalesced);
    auto continuation_stops = py::enum_<ContinuationStop>(module, "ContinuationStop");
    for (int value = static_cast<int>(ContinuationStop::reached_target);
         value <= static_cast<int>(ContinuationStop::critical_endpoint); ++value) {
        const auto stop = static_cast<ContinuationStop>(value);
        continuation_stops.value(continuation_stop_name(stop), stop);
    }
    py::class_<ContinuationOptions>(module, "ContinuationOptions", "Ordinary path: path-fraction steps only halve, so initial_step is the largest. Envelope: steps are capped at envelope_max_log_step; a rejected crossing or closing step ends the trace.")
        .def(py::init<>())
        .def(py::init([](double initial_step, double minimum_step, std::size_t max_steps, std::size_t max_retries_per_step,
                         double envelope_max_log_step, double envelope_straddle) {
            return ContinuationOptions{initial_step, minimum_step, max_steps, max_retries_per_step, envelope_max_log_step, envelope_straddle};
        }), py::kw_only(), py::arg("initial_step"), py::arg("minimum_step"), py::arg("max_steps"), py::arg("max_retries_per_step"),
            py::arg("envelope_max_log_step") = 0.0, py::arg("envelope_straddle") = 0.0)
        EPCSAFT_READWRITE(ContinuationOptions, initial_step, minimum_step, max_steps, max_retries_per_step,
            envelope_max_log_step, envelope_straddle);
    py::class_<ContinuationRequest>(module, "ContinuationRequest")
        .def(py::init<>())
        EPCSAFT_READWRITE(ContinuationRequest, start_problem, target_problem, start_result, options, parameter_targets);
    py::class_<ContinuationResult>(module, "ContinuationResult")
        EPCSAFT_READONLY(ContinuationResult, success, message, stop, path_values, parameter_values, states,
            reached_path, attempted_path, event_phase, event_pair, attempts, accepted_steps, retries, last_result,
            cause_stop, cause_convergence, cause_message, conditions, critical_points, critical_seeds);
    py::class_<EquilibriumStart>(module, "EquilibriumStart")
        .def(py::init([](std::vector<double> coordinates, std::vector<double> amounts, std::vector<std::vector<double>> compositions,
                         std::vector<std::optional<double>> packing_fractions, std::optional<double> temperature,
                         std::optional<double> pressure, std::optional<EquilibriumResult> previous) {
            return EquilibriumStart{std::move(coordinates), std::move(amounts), std::move(compositions),
                                    std::move(packing_fractions), temperature, pressure, std::move(previous)};
        }), py::kw_only(), py::arg("coordinates") = std::vector<double>{}, py::arg("amounts") = std::vector<double>{},
            py::arg("compositions") = std::vector<std::vector<double>>{}, py::arg("packing_fractions") = std::vector<std::optional<double>>{},
            py::arg("temperature") = py::none(), py::arg("pressure") = py::none(), py::arg("previous") = py::none())
        EPCSAFT_READWRITE(EquilibriumStart, coordinates, amounts, compositions, packing_fractions, temperature,
            pressure, previous);
    py::class_<EquilibriumWorkLimits>(module, "EquilibriumWorkLimits")
        .def(py::init([](std::optional<std::size_t> attempts, std::optional<std::size_t> iterations,
                         std::optional<std::size_t> attempt_iterations, std::optional<double> seconds, std::optional<double> attempt_seconds) {
            return EquilibriumWorkLimits{attempts, iterations, attempt_iterations, seconds, attempt_seconds};
        }), py::kw_only(), py::arg("attempts") = py::none(), py::arg("iterations") = py::none(),
            py::arg("attempt_iterations") = py::none(), py::arg("seconds") = py::none(), py::arg("attempt_seconds") = py::none())
        EPCSAFT_READWRITE(EquilibriumWorkLimits, attempts, iterations, attempt_iterations, seconds, attempt_seconds);
    py::class_<EquilibriumStartsResult>(module, "EquilibriumStartsResult")
        EPCSAFT_READONLY(EquilibriumStartsResult, selected, results, causes, attempt_seconds, iterations, seconds,
            overrun_seconds, budget_exhausted);
    module.def("solve_equilibrium_starts", &solve_equilibrium_starts);
    module.def("continue_equilibrium", &continue_equilibrium, py::arg("model"), py::arg("request"));
    // Solved-state actions: the existing native records and owner, transported unchanged.
    py::enum_<SolvedStateObservableKind>(module, "SolvedStateObservableKind")
        EPCSAFT_VALUES(SolvedStateObservableKind, Temperature, Pressure, PhaseAmount, PhaseComponentAmount,
            PhaseMoleFraction, PhaseVolume, PhaseMolarDensity, PhasePackingFraction, PhasePressure,
            PhaseResidualChemicalPotentialOverRt, PhaseProperty, PhaseLogActivity, PhaseGibbsMolarVolume, SystemEnthalpy,
            EquilibriumIsobaricHeatCapacity);
    py::class_<SolvedStateObservable>(module, "SolvedStateObservable")
        .def(py::init<SolvedStateObservableKind, std::size_t, std::size_t, PropertyObservable>(), py::arg("kind"), py::arg("phase") = 0, py::arg("component") = 0,
            py::arg("property") = PropertyObservable::IdealEnthalpy)
        EPCSAFT_READONLY(SolvedStateObservable, kind, phase, component, property);
    py::class_<SolvedStateActionDirection>(module, "SolvedStateActionDirection")
        .def(py::init<double, double, std::vector<double>, std::vector<double>, std::vector<std::array<double, 4>>, double, double>(),
             py::arg("temperature_k"), py::arg("pressure_pa"), py::arg("feed_amounts_mol"),
             py::arg("active_parameters"),
             py::arg("reaction_coefficients") = std::vector<std::array<double, 4>>{},
             py::arg("enthalpy_j") = 0.0, py::arg("entropy_j_per_k") = 0.0)
        EPCSAFT_READONLY(SolvedStateActionDirection, temperature_k, pressure_pa, feed_amounts_mol, active_parameters,
            reaction_coefficients, enthalpy_j, entropy_j_per_k);
    py::class_<SolvedStateActionRequest>(module, "SolvedStateActionRequest")
        .def(py::init<std::vector<SolvedStateObservable>, std::vector<SolvedStateActionDirection>>(), py::arg("observables"), py::arg("directions"))
        EPCSAFT_READONLY(SolvedStateActionRequest, observables, directions);
    py::enum_<SolvedStateActionStatus>(module, "SolvedStateActionStatus")
        EPCSAFT_VALUES(SolvedStateActionStatus, Available, CentralSolveUnavailable, DefiningEquationsNotSatisfied,
            PhaseBoundaryNotDifferentiable, NonuniqueCoalescedPhases, RankDeficient, NonfiniteJacobian,
            LinearSolveResidual, SupportChange, IonicInsertion, EnergyOrderUnavailable, DensityUnavailable,
            PressureSingular, AssociationSingular, AssociationIllConditioned, AssociationResidual,
            MissingThermochemistry, ReferenceUnavailable, OutsideInterval, ThermalBoundary, CentralPropertyUnavailable,
            Nonfinite, IllConditioned);
    py::class_<SolvedStateActionResult>(module, "SolvedStateActionResult")
        EPCSAFT_READONLY(SolvedStateActionResult, observable, unit, order, value, action, status, property_diagnostics,
            diagnostic_message);
    py::class_<SolvedStateActionDiagnostics>(module, "SolvedStateActionDiagnostics")
        EPCSAFT_READONLY(SolvedStateActionDiagnostics, equation_count, state_dimension, maximum_energy_order,
            state_factorizations, linear_solves, rank, condition_estimate, maximum_linear_residual, central_residual,
            central_floor);
    py::class_<SolvedStateActionBatch>(module, "SolvedStateActionBatch")
        EPCSAFT_READONLY(SolvedStateActionBatch, results, diagnostics);
    module.def("solved_state_actions", &solved_state_actions,
               py::arg("compiled"), py::arg("central"), py::arg("request"),
               "Exact order 1-3 actions at `central`, the accepted solve of this same compiled problem.\n\n"
               "Directions carry K, Pa, absolute feed mol (component width) and active-parameter increments\n"
               "(active width, native units); a free T or P entry must be zero. No solve is repeated.");
    module.def("solved_state_first_actions", &solved_state_first_actions);
    py::class_<InventoryBound>(module, "InventoryBound")
        EPCSAFT_READONLY(InventoryBound, species, status, termination, maximum, feed_normalized_amount, log_fraction,
            balance_residual, normalization_residual, lower_bound_residual, witness, amounts_mol, finite_witness,
            limiting, supported, floor, lp_error, finite_residual, feasibility_termination);
    module.def("inventory_bounds", &inventory_bounds, py::arg("compiled"), py::arg("species"),
        py::arg("fraction"), py::arg("floor"));
    module.def("equilibrium_asl_record", [](py::object model, py::object problem, const ContinuationOptions& continuation,
            std::vector<double> offset, std::vector<std::vector<double>> matrix, std::vector<double> anchor,
            std::size_t identity_component, double minimum_gap) {
        const py::object document = model.attr("parameters").attr("_document");
        if (document.is_none())  // fitted parameters carry no source document to export (#29)
            throw std::invalid_argument("ASL record requires parameters loaded from a parameter document; fitted parameters have none");
        const py::object native_owner = model.attr("_model");
        const auto& native = native_owner.cast<const HelmholtzModel&>();
        if (!native.active_parameters().empty() || native.thermochemistry())
            throw std::invalid_argument("ASL record requires inactive parameters and explicit problem reaction references");
        const auto p = problem.attr("_native")(model).cast<EquilibriumProblem>();
        AslFeedMap map{std::move(offset), std::move(anchor), std::move(matrix), identity_component, minimum_gap};
        map.validate(p, native.parameters());
        compile_problem(native, p);
        return asl_record(document.cast<ParameterDocument>(), p, continuation, map);
    }, py::arg("model"), py::arg("problem"), py::arg("continuation"), py::arg("feed_offset"),
       py::arg("feed_matrix"), py::arg("anchor_controls"), py::arg("identity_component"), py::arg("minimum_gap"));
    module.def("compile_problem", &compile_problem, py::keep_alive<0, 1>());
    module.def("sensitivity", &sensitivity);
    module.def("solve_equilibrium", &solve_equilibrium, py::arg("model"), py::arg("problem"),
               py::arg("start") = std::nullopt);
}

}  // namespace epcsaft::equilibrium
