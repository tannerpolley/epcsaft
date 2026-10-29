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

static nb::dict convergence_values(const core::Convergence& c) {
    nb::dict values;
    const bool measured = c.criterion || !std::isnan(c.residual);
    values["residual"] = std::isnan(c.residual) ? nb::none() : nb::cast(c.residual);
    values["floor"] = nb::cast(c.floor); values["amplification"] = nb::cast(c.amplification);
    values["contraction"] = nb::cast(c.contraction);
    values["iterations"] = measured ? nb::cast(c.iterations) : nb::none();
    values["evaluations"] = measured ? nb::cast(c.evaluations) : nb::none();
    values["stop"] = measured ? nb::cast(core::stop_name(c.stop)) : nb::none();
    values["criterion"] = c.criterion ? nb::cast(c.criterion->value) : nb::none();
    return values;
}

static nb::object failure_values(const std::optional<core::Convergence>& record) {
    if (!record) return nb::none();
    auto values = convergence_values(*record);
    values["cause"] = record->cause; values["rank"] = nb::cast(record->rank);
    values["backward_error"] = nb::cast(record->backward_error);
    return values;
}

static nb::dict saturation_values(const EquilibriumResult& r) {
    nb::dict out;
    const auto finite = [](double value) -> nb::object { return std::isfinite(value) ? nb::cast(value) : nb::none(); };
    out["convergence"] = std::isnan(r.convergence.residual) && !r.convergence.criterion
        ? nb::object(nb::none()) : nb::object(convergence_values(r.convergence));
    out["P"] = nb::cast(r.pressure); out["T"] = nb::cast(r.temperature);
    out["success"] = r.success; out["requested_tolerance_met"] = nb::cast(r.requested_tolerance_met);
    out["message"] = r.message; out["solver_status"] = r.solver_status;
    out["max_residual"] = r.residuals.empty() ? nb::none() : finite(r.residual_norm());
    out["iterations"] = r.iterations; out["evaluations"] = r.evaluations; out["hessian_calls"] = r.hessian_calls;
    out["first_numerical_failure"] = failure_values(r.first_numerical_failure);
    out["terminal_numerical_failure"] = failure_values(r.terminal_numerical_failure);
    out["last_full_evaluation"] = failure_values(r.last_full_evaluation);
    nb::object rank = nb::none();
    if (r.last_rank) {
        nb::dict values;
        values["basis"] = r.last_rank->basis; values["dimension"] = r.last_rank->dimension;
        values["rank"] = r.last_rank->rank; values["threshold"] = r.last_rank->threshold;
        rank = std::move(values);
    }
    out["rank"] = std::move(rank);
    nb::object linear = nb::none();
    if (r.linear_solve_diagnostics) {
        const auto& d = *r.linear_solve_diagnostics;
        nb::dict values;
        values["method"] = d.system; values["reason"] = d.reason;
        values["residual"] = finite(d.residual_norm); values["evaluation_floor"] = finite(d.evaluation_floor);
        values["inverse_norm"] = finite(d.jacobian_inverse_norm); values["condition_estimate"] = finite(d.condition_estimate);
        values["forward_error_estimate"] = finite(d.forward_error_estimate);
        values["working_precision_bits"] = d.precision_bits;
        values["refinement_iterations"] = d.refinement_iterations;
        linear = std::move(values);
    }
    out["linear"] = std::move(linear);
    nb::list phases;
    for (std::size_t p = 0; r.has_state() && p < r.phase_count; ++p) {
        nb::dict phase;
        phase["phase"] = r.declaration->phases[p].kind == PhaseKind::Liquid ? "liquid" : "vapor";
        phase["rho"] = r.molar_densities[p];
        const auto n = r.mole_fractions.size() / r.phase_count;
        phase["x"] = nb::cast(std::vector<double>(r.mole_fractions.begin() + p * n, r.mole_fractions.begin() + (p + 1) * n));
        phase["P"] = p < r.phase_pressures.size() ? nb::cast(r.phase_pressures[p]) : nb::none();
        phases.append(std::move(phase));
    }
    out["phases"] = std::move(phases);
    out["rejected_reduction"] = r.rejected_reduction ? nb::object(saturation_values(*r.rejected_reduction)) : nb::object(nb::none());
    return out;
}

void register_equilibrium(nb::module_& module) {
    nb::class_<core::NumericalRank>(module, "NumericalRank")
        EPCSAFT_READONLY(core::NumericalRank, basis, dimension, rank, threshold);
    using core::Quantity; using core::ComparisonMode; using core::NumericalEstimate;
    using core::ComparisonEligibility; using core::EquilibriumComparison;
    using core::compare_equilibrium_quantity;
    using core::project_numerical_estimate;
    nb::enum_<Quantity>(module, "Quantity") EPCSAFT_VALUES(Quantity, Linear, Positive, LogPositive);
    nb::enum_<ComparisonMode>(module, "ComparisonMode") EPCSAFT_VALUES(ComparisonMode, Strict, EstimatedError);
    nb::class_<NumericalEstimate>(module, "NumericalEstimate").def(nb::init<>())
        EPCSAFT_READWRITE(NumericalEstimate, absolute, chart, available, resolved, cause);
    nb::class_<ComparisonEligibility>(module, "ComparisonEligibility").def(nb::init<>())
        EPCSAFT_READWRITE(ComparisonEligibility, declaration, branch, phase_matching, inventory,
            residual_a, residual_b, admission_a, admission_b);
    nb::class_<EquilibriumComparison>(module, "EquilibriumComparison")
        EPCSAFT_READONLY(EquilibriumComparison, judgment, cause, difference, chart_difference, radius_sum);
    module.def("compare_equilibrium_quantity", &compare_equilibrium_quantity,
        nb::arg("quantity"), nb::arg("a"), nb::arg("b"), nb::arg("ordinary_allowance"),
        nb::arg("estimate_a"), nb::arg("estimate_b"), nb::arg("eligibility"), nb::arg("mode"));
    module.def("project_numerical_estimate", &project_numerical_estimate<double>,
        nb::arg("value"), nb::arg("mapped"), nb::arg("radius"), nb::arg("positive"), nb::arg("constant"));
    nb::class_<EquilibriumNumericalEstimates>(module, "EquilibriumNumericalEstimates")
        EPCSAFT_READONLY(EquilibriumNumericalEstimates, basis, cause, source_reference_method, rank, dimension, precision_bits, threshold,
            condition, inverse_norm, backward_error, residuals, row_terms, row_inner, row_floor, physical_scale,
            temperature, pressure, loading, densities, packing, composition, amounts);
    module.def("equilibrium_numerical_estimates", nb::overload_cast<const HelmholtzModel&, const EquilibriumResult&>(
        &equilibrium_numerical_estimates), nb::arg("model"), nb::arg("result"));
    module.def("equilibrium_numerical_estimates", nb::overload_cast<const IdealMoleFraction&, const EquilibriumResult&>(
        &equilibrium_numerical_estimates), nb::arg("model"), nb::arg("result"));
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
            CriticalNormalization, CriticalCubic, EnvelopeCondition, AzeotropeCondition, SolidAffinity, ReactionMassAction, ReactionMolarity, SaturationPressure);
    nb::class_<Row>(module, "Row")
        EPCSAFT_READONLY(Row, kind, phase, index)
        .def_prop_ro("eqid", [](const Row& row) { return row_kind_name(row.kind); });
    nb::class_<PointEvaluation>(module, "PointEvaluation")
        EPCSAFT_READONLY(PointEvaluation, values, objective, floor);
    nb::class_<CompiledSystem>(module, "CompiledSystem")
        .def_prop_ro("rows", [](const CompiledSystem& system) { return system.rows; })
        EPCSAFT_READONLY(CompiledSystem, first_numerical_failure, rank_records, last_rank, pinned_count, free_intensive_count, rank, charge_row_count, primal_size,
            generic_nullity, feed_total, has_ions, has_reaction, conserves_total_moles, jacobian_rank, removed_species)
        .def("solve", [](const CompiledSystem& system) { return solve_compiled(system); }, nb::keep_alive<0, 1>())
        .def("solve_saturation", &solve_saturation, nb::keep_alive<0, 1>())
        .def("generic_point", &CompiledSystem::generic_point)
        .def("numerical_estimates", [](const CompiledSystem& system, std::vector<double> coordinates,
                double temperature, double pressure, std::vector<double> amounts, std::vector<double> composition,
                std::vector<double> densities, std::vector<double> packing) {
            const auto phases = system.problem.phases.size(), components = system.n_feed.size();
            if (coordinates.size() != system.unknown_count()
                || amounts.size() != phases || composition.size() != phases * components
                || (!system.ideal && (densities.size() != phases || packing.size() != phases)))
                throw std::invalid_argument("numerical estimate snapshot shape");
            EquilibriumResult result;
            result.success = true; result.temperature = temperature; result.pressure = pressure;
            result.declaration = std::make_shared<const EquilibriumProblem>(system.problem);
            result.coordinates = std::move(coordinates); result.amounts = std::move(amounts);
            result.mole_fractions = std::move(composition); result.molar_densities = std::move(densities);
            result.packing_fractions = std::move(packing);
            return equilibrium_numerical_estimates(system, result);
        }, nb::arg("coordinates"), nb::arg("temperature"), nb::arg("pressure"), nb::arg("amounts"),
           nb::arg("composition"), nb::arg("densities"), nb::arg("packing"))
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
    nb::class_<LinearSolveDiagnostics>(module, "LinearSolveDiagnostics")
        EPCSAFT_READONLY(LinearSolveDiagnostics, system, reason, residual_norm, evaluation_floor, jacobian_inverse_norm,
            condition_estimate, forward_error_estimate, precision_bits, refinement_iterations);
    nb::class_<IpoptIteration>(module, "IpoptIteration")
        EPCSAFT_READONLY(IpoptIteration, mode, iteration, line_search_trials, objective, inf_pr, inf_du, mu, step_norm, regularization, alpha_du, alpha_pr);
    nb::class_<IpoptPoint>(module, "IpoptPoint")
        EPCSAFT_READONLY(IpoptPoint, tnlp_coordinates, coordinate_scale, tnlp_lower_multipliers, tnlp_upper_multipliers,
            tnlp_constraint_multipliers, declaration, rows, constraint_rows, objective_factor);
    nb::class_<IpoptDiagnostics>(module, "IpoptDiagnostics")
        EPCSAFT_READONLY(IpoptDiagnostics, last, first_restoration_iteration, initialization_return, application_return,
            restoration_entries, callbacks, warm_initialized, predicted, physical_amounts, failed_point);
    nb::class_<EquilibriumResult>(module, "EquilibriumResult", "Phases in declared order; K, Pa, mol, mol/m^3; per-phase-per-component lists are phase-major in component_ids order.")
        .def("saturation_values", [](nb::object self) {
            auto values = saturation_values(nb::cast<const EquilibriumResult&>(self));
            values["diagnostics"] = self;
            return values;
        })
        .def_prop_ro("rejected_reduction", [](const EquilibriumResult& r) { return r.rejected_reduction ? nb::cast(*r.rejected_reduction) : nb::none(); })
        .def_prop_ro("solver_options", [](const EquilibriumResult& v) { return v.declaration ? bindings::solver_options_view(v.declaration->solver_options, 7) : nb::list(); })
        EPCSAFT_READONLY(EquilibriumResult, solver_settings, diagnostic_options_file, linear_solver, first_numerical_failure, rank_records, ipopt_diagnostics, recovered_ipopt_diagnostics)
        .def_ro("neutral_reference_diagnostics", &EquilibriumResult::reference_diagnostics)
        EPCSAFT_READONLY(EquilibriumResult, terminal_numerical_failure, last_full_evaluation, phase_pressures, discovery, success, message, solver_status, temperature, pressure, amounts,
            mole_fractions, molar_densities, packing_fractions,
            // Dimensionless reduced-energy balance/charge stationarity coordinates of the solution.
            multipliers,
            // Extensive-objective multipliers of the non-chemical, non-pressure `rows` in order, formed at
            // the returned point from `multipliers`; raw optimizer duals only when a phase amount is pinned.
            constraint_multipliers, residuals, classified_residuals, coordinates, phase_count, iterations, evaluations,
            hessian_calls, objective, reduced_hessian_minimum, convergence, last_rank, linear_solve_diagnostics, numerical_estimates, rows, unknown_count, jacobian_rank,
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
        EPCSAFT_READONLY(StabilityReport, first_numerical_failure, last_rank, completed, message, trials, candidate_coordinates, standard_over_rt,
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
        .def_prop_ro("first_numerical_failure", [](const ContinuationResult& r) { return r.observations.first_numerical_failure; })
        .def_prop_ro("rank_records", [](const ContinuationResult& r) { return r.observations.rank_records; })
        .def_prop_ro("ipopt_diagnostics", [](const ContinuationResult& r) { return r.observations.ipopt_diagnostics; })
        .def_prop_ro("recovered_ipopt_diagnostics", [](const ContinuationResult& r) { return r.observations.recovered_ipopt_diagnostics; })
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
    module.def("sensitivity", [](const CompiledSystem& compiled, const std::vector<double>& point, const std::vector<double>& lever) {
        std::optional<core::Convergence> failure;
        auto result = sensitivity(compiled, point, lever, &failure);
        if (failure) throw core::ConvergenceFailure(*failure);
        return result;
    }, nb::arg("compiled"), nb::arg("point"), nb::arg("lever"));
    module.def("solve_equilibrium", nb::overload_cast<const HelmholtzModel&, const EquilibriumProblem&, const std::optional<EquilibriumResult>&>(&solve_equilibrium), nb::arg("model"), nb::arg("problem"),
               nb::arg("start") = std::nullopt);
    module.def("solve_equilibrium", nb::overload_cast<const IdealMoleFraction&, const EquilibriumProblem&, const std::optional<EquilibriumResult>&>(&solve_equilibrium),
        nb::arg("model"), nb::arg("problem"), nb::arg("start") = std::nullopt, nb::keep_alive<0, 1>());
}

}  // namespace epcsaft::equilibrium
