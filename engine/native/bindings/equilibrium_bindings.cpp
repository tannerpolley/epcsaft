#include "equilibrium/problem.hpp"
#include "equilibrium/gibbs_minimizer.hpp"
#include "equilibrium_record.hpp"

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

#include <optional>

namespace nb = nanobind;
namespace epcsaft::equilibrium {

// A fit the start document cannot represent carries the refusal message in place of the document.
static ParameterDocument record_document(const nb::object& model) {
    const nb::object document = model.attr("parameters").attr("_document");
    if (document.is_none()) throw std::invalid_argument("ASL record: parameters built without a parameter document");
    if (nb::isinstance<nb::str>(document)) throw std::invalid_argument("ASL record: " + nb::cast<std::string>(document));
    return nb::cast<ParameterDocument>(document);
}

void register_equilibrium(nb::module_& module) {
    module.def("property_asl_record", [](nb::object model) {
        const auto document = record_document(model);
        const nb::object native_owner = model.attr("_model");
        const auto& native = nb::cast<const HelmholtzModel&>(native_owner);
        if (!native.active_parameters().empty() || !native.thermochemistry())
            throw std::invalid_argument("ASL property record requires inactive parameters with thermochemistry");
        return asl_property_record(document, *native.thermochemistry(), native.solver_options());
    }, nb::arg("model"));
    module.attr("IDEAL_MOLE_FRACTION_STANDARD_STATE_ID") = IDEAL_MOLE_FRACTION_STANDARD_STATE_ID;
    module.attr("IDEAL_MOLARITY_STANDARD_STATE_ID") = IDEAL_MOLARITY_STANDARD_STATE_ID;
    nb::class_<IdealMoleFraction>(module, "IdealMoleFraction")
        .def("__init__", [](IdealMoleFraction* self, std::vector<std::string> ids, std::vector<int> charges, nb::object options,
            std::optional<std::string> basis, std::optional<double> volume, std::optional<std::string> solvent) {
            new (self) IdealMoleFraction(std::move(ids), std::move(charges), bindings::solver_options(options, 4),
                std::move(basis), volume, std::move(solvent));
        }, nb::arg("component_ids"), nb::arg("charge_number"), nb::arg("solver_options").none() = nb::none(),
           nb::arg("activity_basis").none() = nb::none(), nb::arg("solution_volume_m3").none() = nb::none(),
           nb::arg("unit_activity_solvent").none() = nb::none())
        EPCSAFT_READONLY(IdealMoleFraction, component_ids, charge_number, activity_basis, solution_volume_m3, unit_activity_solvent)
        EPCSAFT_SOLVER_OPTIONS(IdealMoleFraction, 4);
    module.def("effective_water_amounts", &effective_water_amounts, nb::arg("model"), nb::arg("result"),
        nb::arg("water"), nb::arg("ions"), nb::arg("molecular_species"));
    nb::class_<SpeciesActivities>(module, "SpeciesActivities")
        EPCSAFT_READONLY(SpeciesActivities, log_activities, activities, causes, solvent_index, standard_molality_mol_per_kg);
    nb::class_<ReactionActivities>(module, "ReactionActivities")
        EPCSAFT_READONLY(ReactionActivities, phase_index, reaction_index, basis, cause, reference_pressure_pa,
            log_quotient, species);
    module.def("reaction_activities", nb::overload_cast<const HelmholtzModel&, const EquilibriumResult&>(&reaction_activities), nb::arg("model"), nb::arg("result"));
    module.def("reaction_activities", nb::overload_cast<const IdealMoleFraction&, const EquilibriumResult&>(&reaction_activities), nb::arg("model"), nb::arg("result"));
    module.attr("EOS_STANDARD_STATE_ID") = EOS_STANDARD_STATE_ID;
    nb::enum_<PhaseKind>(module, "PhaseKind")
        EPCSAFT_VALUES(PhaseKind, Liquid, Vapor);
    nb::enum_<AmountKind>(module, "AmountKind")
        EPCSAFT_VALUES(AmountKind, Free, Pinned);
    nb::enum_<FeedKind>(module, "FeedKind")
        EPCSAFT_VALUES(FeedKind, MoleFractions, Amounts, MassFractions, Molalities);
    nb::class_<AmountSpec>(module, "AmountSpec")
        .def(nb::init<>())
        EPCSAFT_READWRITE(AmountSpec, kind, fraction)
        .def_static("free", &AmountSpec::free)
        .def_static("pinned", &AmountSpec::pinned);
    nb::class_<FreeValue>(module, "FreeValue")
        .def(nb::init<>())
        EPCSAFT_READWRITE(FreeValue, free, value)
        .def_static("fixed", &FreeValue::fixed)
        .def_static("variable", [](std::optional<double> guess) {
            return guess.has_value() ? FreeValue::variable(*guess) : FreeValue::variable();
        }, nb::arg("guess") = nb::none());
    nb::class_<PhaseSpec>(module, "PhaseSpec")
        .def(nb::init<>())
        EPCSAFT_READWRITE(PhaseSpec, name, kind, support, support_declared, amount, composition_guess, packing_guess);
    nb::class_<FeedSpec>(module, "FeedSpec")
        .def(nb::init<>())
        EPCSAFT_READWRITE(FeedSpec, kind, values, solvent, solvent_mass_kg);
    nb::enum_<SpeciesReference>(module, "SpeciesReference")
        EPCSAFT_VALUES(SpeciesReference, InfiniteDilution, PureLiquid, PureGas);
    nb::class_<NeutralReferenceDeclaration>(module, "NeutralReferenceDeclaration")
        .def(nb::init<>())
        EPCSAFT_READWRITE(NeutralReferenceDeclaration, solvent_index, solute_molality_ratios,
            standard_molality_mol_per_kg, species_reference);
    nb::enum_<ReferenceSourceBasis>(module, "ReferenceSourceBasis")
        EPCSAFT_VALUES(ReferenceSourceBasis, CommonMolalityInfiniteDilution, RawMoleFractionInfiniteDilution);
    nb::enum_<NeutralReferenceStatus>(module, "NeutralReferenceStatus")
        EPCSAFT_VALUES(NeutralReferenceStatus, Available, InvalidDeclaration, NonNeutralDeclaration,
            UnsupportedSourceBasis, FreeTemperature, FreePressure, DensityUnavailable, BranchChange, Nonfinite);
    nb::class_<NeutralReferenceReactionDiagnostic>(module, "NeutralReferenceReactionDiagnostic")
        EPCSAFT_SOLVER_OPTIONS(NeutralReferenceReactionDiagnostic, 7)
        EPCSAFT_READONLY(NeutralReferenceReactionDiagnostic, reaction_index, terminal_contraction,
            observed_terminal_change, diagnostics, density_convergence, terminal_contraction_pressure_derivative,
            terminal_contraction_pressure_second_derivative, reference_pressure_pa, status, message);
    nb::class_<ReactionReference>(module, "ReactionReference")
        .def("__init__", [](ReactionReference* self, ReferenceSourceBasis basis, std::optional<double> pressure) {
            new (self) ReactionReference{pressure, basis};
        }, nb::arg("source_basis"), nb::arg("reference_pressure_pa") = nb::none())
        EPCSAFT_READWRITE(ReactionReference, source_basis, reference_pressure_pa);
    nb::class_<ReactionSpec>(module, "ReactionSpec")
        .def(nb::init<>())
        EPCSAFT_READWRITE(ReactionSpec, name, stoichiometry, equilibrium_constant, affinity_j_per_mol, correlation,
            reference);
    nb::class_<SolidSpec>(module, "SolidSpec")
        .def(nb::init<>())
        EPCSAFT_READWRITE(SolidSpec, name, stoichiometry, equilibrium_constant, correlation, reference,
            molar_volume_m3_per_mol, amount);
    nb::class_<StabilityOptions>(module, "StabilityOptions")
        .def(nb::init<>())
        EPCSAFT_READWRITE(StabilityOptions, evaluate_supplied_trials, trials, declaration_budget, trial_budget,
            standard_uncertainty_over_rt);
    nb::class_<CriticalDeclaration>(module, "CriticalDeclaration")
        .def(nb::init<>())
        EPCSAFT_READWRITE(CriticalDeclaration, density_guess_mol_per_m3, mode_guess);
    nb::enum_<EnvelopeCoordinate>(module, "EnvelopeCoordinate")
        EPCSAFT_VALUES(EnvelopeCoordinate, LogTemperature, LogPressure, LogK);
    nb::class_<EnvelopeCondition>(module, "EnvelopeCondition")
        .def(nb::init<>())
        EPCSAFT_READWRITE(EnvelopeCondition, kind, component, target);
    nb::class_<EquilibriumProblem>(module, "EquilibriumProblem")
        .def_prop_rw("solver_options", [](const EquilibriumProblem& v) { return bindings::solver_options_view(v.solver_options, 4); },
            [](EquilibriumProblem& v, nb::object mapping) { v.solver_options = bindings::solver_options(mapping, 4); }, nb::for_setter(nb::arg("mapping").none()))
        .def_prop_rw("row_tolerance", [](const EquilibriumProblem& v) { return v.solver_options.value(policy::row); },
            [](EquilibriumProblem& v, nb::object value) { bindings::solver_option(v.solver_options, policy::EQUILIBRIUM_ROW_RESIDUAL.name, value, 4); }, nb::for_setter(nb::arg("value").none()))
        .def(nb::init<>())
        EPCSAFT_READWRITE(EquilibriumProblem, phases, temperature, pressure, feed, reactions)
        .def_rw("neutral_reference", &EquilibriumProblem::source_reference)
        EPCSAFT_READWRITE(EquilibriumProblem, standard_gibbs_j_per_mol, topology_declared, enthalpy_j,
            entropy_j_per_k, critical, envelope, azeotrope, discovery, solids);
    nb::enum_<RowKind>(module, "RowKind")
        EPCSAFT_VALUES(RowKind, Chemical, Pressure, Moiety, Electro, Normalization, Enthalpy, Entropy, CriticalMode,
            CriticalNormalization, CriticalCubic, EnvelopeCondition, AzeotropeCondition, SolidAffinity, ReactionMassAction, ReactionMolarity);
    nb::class_<Row>(module, "Row")
        EPCSAFT_READONLY(Row, kind, phase, index)
        .def_prop_ro("eqid", [](const Row& row) { return row_kind_name(row.kind); });
    nb::class_<PointEvaluation>(module, "PointEvaluation")
        EPCSAFT_READONLY(PointEvaluation, values, objective, floor);
    nb::class_<CompiledSystem>(module, "CompiledSystem")
        .def_prop_ro("rows", [](const CompiledSystem& system) { return system.rows; })
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
            nb::arg("point"), nb::arg("constraint_multipliers") = std::vector<double>{});
    nb::class_<EquilibriumResult>(module, "EquilibriumResult", "Phases in declared order; K, Pa, mol, mol/m^3; per-phase-per-component lists are phase-major in component_ids order.")
        .def_prop_ro("solver_options", [](const EquilibriumResult& v) { return v.declaration ? bindings::solver_options_view(v.declaration->solver_options, 7) : nb::list(); })
        EPCSAFT_READONLY(EquilibriumResult, solver_settings, diagnostic_options_file, linear_solver)
        .def_ro("neutral_reference_diagnostics", &EquilibriumResult::reference_diagnostics)
        EPCSAFT_READONLY(EquilibriumResult, discovery, success, message, solver_status, temperature, pressure, amounts,
            mole_fractions, molar_densities, packing_fractions,
            // Dimensionless reduced-energy balance/charge stationarity coordinates of the solution.
            multipliers,
            // Extensive-objective multipliers of the non-chemical, non-pressure `rows` in order, formed at
            // the returned point from `multipliers`; raw optimizer duals only when a phase amount is pinned.
            constraint_multipliers, residuals, classified_residuals, coordinates, phase_count, iterations, evaluations,
            hessian_calls, objective, reduced_hessian_minimum, convergence, rows, unknown_count, jacobian_rank,
            generic_nullity, requested_tolerance_met, topology_event, topology_phase, topology_pair, critical,
            // Per declared solid: amount in mol and saturation index SI_s in decades.
            solid_amounts, solid_affinity_decades);
    nb::class_<CriticalRecord>(module, "CriticalRecord")
        EPCSAFT_READONLY(CriticalRecord, eigenvalues, mode, cubic, zero_band, condition_estimate, backward_error, rank);
    nb::class_<StabilityTrial>(module, "StabilityTrial")
        EPCSAFT_READONLY(StabilityTrial, phase, seed, evidence, message, optimizer, raw_coordinates, raw_residuals,
            composition, split_coordinates, basis, row_transformation, quantities, split_checks, donor,
            resolved_descent);
    nb::class_<StabilityReport>(module, "StabilityReport")
        EPCSAFT_READONLY(StabilityReport, completed, message, trials, candidate_coordinates, standard_over_rt,
            plane_over_rt, present_phases, work);
    nb::class_<DiscoveryRecord>(module, "DiscoveryRecord")
        EPCSAFT_READONLY(DiscoveryRecord, declaration, candidate, rounds, declarations, declaration_roles,
            declaration_attempt_coordinates, work);
    module.def("cooling_onset", &cooling_onset, nb::arg("compiled"));
    module.def("stability_test", &stability_test, nb::arg("compiled"), nb::arg("central"),
        nb::arg("options") = StabilityOptions{});
    nb::enum_<TopologyEvent>(module, "TopologyEvent")
        .value("None_", TopologyEvent::None)
        EPCSAFT_VALUES(TopologyEvent, PhaseDisappeared, PhasesCoalesced);
    auto continuation_stops = nb::enum_<ContinuationStop>(module, "ContinuationStop");
    for (int value = static_cast<int>(ContinuationStop::reached_target);
         value <= static_cast<int>(ContinuationStop::critical_endpoint); ++value) {
        const auto stop = static_cast<ContinuationStop>(value);
        continuation_stops.value(continuation_stop_name(stop), stop);
    }
    nb::class_<ContinuationOptions>(module, "ContinuationOptions", "Ordinary path: path-fraction steps only halve, so initial_step is the largest. Envelope: steps are capped at envelope_max_log_step; a rejected crossing or closing step ends the trace.")
        .def(nb::init<>())
        .def("__init__", [](ContinuationOptions* self, double initial_step, double minimum_step, std::size_t max_steps, std::size_t max_retries_per_step,
                         double envelope_max_log_step, double envelope_straddle) {
            new (self) ContinuationOptions{initial_step, minimum_step, max_steps, max_retries_per_step, envelope_max_log_step, envelope_straddle};
        }, nb::kw_only(), nb::arg("initial_step"), nb::arg("minimum_step"), nb::arg("max_steps"), nb::arg("max_retries_per_step"),
            nb::arg("envelope_max_log_step") = 0.0, nb::arg("envelope_straddle") = 0.0)
        EPCSAFT_READWRITE(ContinuationOptions, initial_step, minimum_step, max_steps, max_retries_per_step,
            envelope_max_log_step, envelope_straddle);
    nb::class_<ContinuationRequest>(module, "ContinuationRequest")
        .def(nb::init<>())
        EPCSAFT_READWRITE(ContinuationRequest, start_problem, target_problem, start_result, options, parameter_targets);
    nb::class_<ContinuationResult>(module, "ContinuationResult")
        EPCSAFT_READONLY(ContinuationResult, success, message, stop, path_values, parameter_values, states,
            reached_path, attempted_path, event_phase, event_pair, attempts, accepted_steps, retries, last_result,
            cause_stop, cause_convergence, cause_message, conditions, critical_points, critical_seeds);
    nb::class_<EquilibriumStart>(module, "EquilibriumStart")
        .def("__init__", [](EquilibriumStart* self, std::vector<double> coordinates, std::vector<double> amounts, std::vector<std::vector<double>> compositions,
                         std::vector<std::optional<double>> packing_fractions, std::optional<double> temperature,
                         std::optional<double> pressure, std::optional<EquilibriumResult> previous) {
            new (self) EquilibriumStart{std::move(coordinates), std::move(amounts), std::move(compositions),
                                    std::move(packing_fractions), temperature, pressure, std::move(previous)};
        }, nb::kw_only(), nb::arg("coordinates") = std::vector<double>{}, nb::arg("amounts") = std::vector<double>{},
            nb::arg("compositions") = std::vector<std::vector<double>>{}, nb::arg("packing_fractions") = std::vector<std::optional<double>>{},
            nb::arg("temperature") = nb::none(), nb::arg("pressure") = nb::none(), nb::arg("previous") = nb::none())
        EPCSAFT_READWRITE(EquilibriumStart, coordinates, amounts, compositions, packing_fractions, temperature,
            pressure, previous);
    nb::class_<EquilibriumWorkLimits>(module, "EquilibriumWorkLimits")
        .def("__init__", [](EquilibriumWorkLimits* self, std::optional<std::size_t> attempts, std::optional<std::size_t> iterations,
                         std::optional<std::size_t> attempt_iterations, std::optional<double> seconds, std::optional<double> attempt_seconds) {
            new (self) EquilibriumWorkLimits{attempts, iterations, attempt_iterations, seconds, attempt_seconds};
        }, nb::kw_only(), nb::arg("attempts") = nb::none(), nb::arg("iterations") = nb::none(),
            nb::arg("attempt_iterations") = nb::none(), nb::arg("seconds") = nb::none(), nb::arg("attempt_seconds") = nb::none())
        EPCSAFT_READWRITE(EquilibriumWorkLimits, attempts, iterations, attempt_iterations, seconds, attempt_seconds);
    nb::class_<EquilibriumStartsResult>(module, "EquilibriumStartsResult")
        EPCSAFT_READONLY(EquilibriumStartsResult, selected, results, causes, attempt_seconds, iterations, seconds,
            overrun_seconds, budget_exhausted);
    module.def("solve_equilibrium_starts", &solve_equilibrium_starts);
    module.def("continue_equilibrium", &continue_equilibrium, nb::arg("model"), nb::arg("request"));
    // Solved-state actions: the existing native records and owner, transported unchanged.
    nb::enum_<SolvedStateObservableKind>(module, "SolvedStateObservableKind")
        EPCSAFT_VALUES(SolvedStateObservableKind, Temperature, Pressure, PhaseAmount, PhaseComponentAmount,
            PhaseMoleFraction, PhaseVolume, PhaseMolarDensity, PhasePackingFraction, PhasePressure,
            PhaseResidualChemicalPotentialOverRt, PhaseProperty, PhaseLogActivity, PhaseGibbsMolarVolume, SystemEnthalpy,
            EquilibriumIsobaricHeatCapacity, SolidAmount, SolidSaturationIndex);
    nb::class_<SolvedStateObservable>(module, "SolvedStateObservable")
        .def(nb::init<SolvedStateObservableKind, std::size_t, std::size_t, PropertyObservable>(), nb::arg("kind"), nb::arg("phase") = 0, nb::arg("component") = 0,
            nb::arg("property") = PropertyObservable::IdealEnthalpy)
        EPCSAFT_READONLY(SolvedStateObservable, kind, phase, component, property);
    nb::class_<SolvedStateActionDirection>(module, "SolvedStateActionDirection")
        .def(nb::init<double, double, std::vector<double>, std::vector<double>, std::vector<std::array<double, 4>>, double, double>(),
             nb::arg("temperature_k"), nb::arg("pressure_pa"), nb::arg("feed_amounts_mol"),
             nb::arg("active_parameters"),
             nb::arg("reaction_coefficients") = std::vector<std::array<double, 4>>{},
             nb::arg("enthalpy_j") = 0.0, nb::arg("entropy_j_per_k") = 0.0)
        EPCSAFT_READONLY(SolvedStateActionDirection, temperature_k, pressure_pa, feed_amounts_mol, active_parameters,
            reaction_coefficients, enthalpy_j, entropy_j_per_k);
    nb::class_<SolvedStateActionRequest>(module, "SolvedStateActionRequest")
        .def(nb::init<std::vector<SolvedStateObservable>, std::vector<SolvedStateActionDirection>>(), nb::arg("observables"), nb::arg("directions"))
        EPCSAFT_READONLY(SolvedStateActionRequest, observables, directions);
    nb::enum_<SolvedStateActionStatus>(module, "SolvedStateActionStatus")
        EPCSAFT_VALUES(SolvedStateActionStatus, Available, CentralSolveUnavailable, DefiningEquationsNotSatisfied,
            PhaseBoundaryNotDifferentiable, NonuniqueCoalescedPhases, RankDeficient, NonfiniteJacobian,
            LinearSolveResidual, SupportChange, IonicInsertion, EnergyOrderUnavailable, DensityUnavailable,
            PressureSingular, AssociationSingular, AssociationIllConditioned, AssociationResidual,
            MissingThermochemistry, ReferenceUnavailable, OutsideInterval, ThermalBoundary, CentralPropertyUnavailable,
            Nonfinite, IllConditioned);
    nb::class_<SolvedStateActionResult>(module, "SolvedStateActionResult")
        EPCSAFT_READONLY(SolvedStateActionResult, observable, unit, order, value, action, status, property_diagnostics,
            diagnostic_message);
    nb::class_<SolvedStateActionDiagnostics>(module, "SolvedStateActionDiagnostics")
        EPCSAFT_READONLY(SolvedStateActionDiagnostics, equation_count, state_dimension, maximum_energy_order,
            state_factorizations, linear_solves, rank, condition_estimate, maximum_linear_residual, central_residual,
            central_floor);
    nb::class_<SolvedStateActionBatch>(module, "SolvedStateActionBatch")
        EPCSAFT_READONLY(SolvedStateActionBatch, results, diagnostics);
    module.def("solved_state_actions", &solved_state_actions,
               nb::arg("compiled"), nb::arg("central"), nb::arg("request"),
               "Exact order 1-3 actions at `central`, the accepted solve of this same compiled problem.\n\n"
               "Directions carry K, Pa, absolute feed mol (component width) and active-parameter increments\n"
               "(active width, native units); a free T or P entry must be zero. No solve is repeated.");
    module.def("solved_state_first_actions", &solved_state_first_actions);
    nb::class_<InventoryBound>(module, "InventoryBound")
        EPCSAFT_READONLY(InventoryBound, species, status, termination, maximum, feed_normalized_amount, log_fraction,
            balance_residual, normalization_residual, lower_bound_residual, witness, amounts_mol, finite_witness,
            limiting, supported, floor, lp_error, finite_residual, feasibility_termination);
    module.def("inventory_bounds", &inventory_bounds, nb::arg("compiled"), nb::arg("species"),
        nb::arg("fraction"), nb::arg("floor"));
    module.def("equilibrium_asl_record", [](nb::object model, nb::object problem, const ContinuationOptions& continuation,
            std::vector<double> offset, std::vector<std::vector<double>> matrix, std::vector<double> anchor,
            std::size_t identity_component, double minimum_gap) {
        const auto document = record_document(model);
        const nb::object native_owner = model.attr("_model");
        const auto& native = nb::cast<const HelmholtzModel&>(native_owner);
        if (!native.active_parameters().empty() || native.thermochemistry())
            throw std::invalid_argument("ASL record requires inactive parameters and explicit problem reaction references");
        const auto p = nb::cast<EquilibriumProblem>(problem.attr("_native")(model));
        if (!p.solids.empty()) throw std::invalid_argument("ASL record does not carry solids");
        AslFeedMap map{std::move(offset), std::move(anchor), std::move(matrix), identity_component, minimum_gap};
        map.validate(p, native.parameters());
        compile_problem(native, p);
        return asl_record(document, p, continuation, map, native.solver_options());
    }, nb::arg("model"), nb::arg("problem"), nb::arg("continuation"), nb::arg("feed_offset"),
       nb::arg("feed_matrix"), nb::arg("anchor_controls"), nb::arg("identity_component"), nb::arg("minimum_gap"));
    module.def("linear_solver", [] { return linear_solver(); });  // the process-level record
    module.def("compile_problem", nb::overload_cast<const HelmholtzModel&, const EquilibriumProblem&>(&compile_problem), nb::keep_alive<0, 1>());
    module.def("compile_problem", nb::overload_cast<const IdealMoleFraction&, const EquilibriumProblem&>(&compile_problem), nb::keep_alive<0, 1>());
    module.def("sensitivity", &sensitivity);
    module.def("solve_equilibrium", nb::overload_cast<const HelmholtzModel&, const EquilibriumProblem&, const std::optional<EquilibriumResult>&>(&solve_equilibrium), nb::arg("model"), nb::arg("problem"),
               nb::arg("start") = std::nullopt);
    module.def("solve_equilibrium", nb::overload_cast<const IdealMoleFraction&, const EquilibriumProblem&, const std::optional<EquilibriumResult>&>(&solve_equilibrium),
        nb::arg("model"), nb::arg("problem"), nb::arg("start") = std::nullopt, nb::keep_alive<0, 1>());
}

}  // namespace epcsaft::equilibrium
