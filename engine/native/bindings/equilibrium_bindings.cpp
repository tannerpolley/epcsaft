#include "equilibrium/problem.hpp"
#include "equilibrium/gibbs_minimizer.hpp"
#include "equilibrium/stoichiometry.hpp"
#include "equilibrium_record.hpp"

#include "domain_reports.hpp"
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
#include <utility>

namespace nb = nanobind;
namespace epcsaft::equilibrium {

// Python opens this only while the epcsaft logger is at DEBUG; it collects the calling thread's Ipopt iterations.
struct IpoptIterations final {
    enum class State { not_entered, active, closed };
    State state = State::not_entered;
    std::vector<IpoptIterationRecord> records;
    std::vector<IpoptIterationRecord>* previous = nullptr;
    IpoptIterations& enter();
    bool exit(const nb::args&);
};
IpoptIterations& IpoptIterations::enter() {
    core::Scope scope("ipopt_iterations");
    if (state == State::active) throw std::runtime_error("capture is already active");
    records.clear();
    state = State::active;
    previous = std::exchange(ipopt_iteration_log(), &records);
    return *this;
}
bool IpoptIterations::exit(const nb::args&) {
    core::Scope scope("ipopt_iterations");
    if (state != State::active) throw std::runtime_error("capture is not active");
    ipopt_iteration_log() = previous;
    state = State::closed;
    return false;
}

static void phase_domains(std::vector<OutsideDomain>& rows, const HelmholtzModel& model,
                          const EquilibriumResult& result, std::size_t row = 0) {
    if (result.has_state())
        for (std::size_t i = 0; i < result.phase_count; ++i)
            outside_domains(rows, model.parameters(), model.domains.get(), result.mole_fractions, result.temperature, result.pressure, row, i);
}

static auto equilibrium_domains(EquilibriumResult&& result, const HelmholtzModel* model) {
    return bindings::domain_result(std::move(result), model,
        [&](auto& rows, const auto& r) { phase_domains(rows, *model, r); });
}

void register_equilibrium(nb::module_& module) {
    nb::class_<core::ConditioningObservation>(module, "ConditioningObservation")
        EPCSAFT_READONLY(core::ConditioningObservation, operation, basis, norm, method, cause, point_kind, criterion,
            dimension, rank, precision_bits, attempt, stage, iteration, threshold, condition, matrix_norm, error,
            criterion_value, coordinates, row_scale, column_scale, controls, phase, component, export_displacement, available)
        EPCSAFT_REPR(ConditioningObservation, "operation", "basis", "rank", "dimension", "condition", "cause", "available");
    nb::class_<core::NumericalRank>(module, "NumericalRank")
        EPCSAFT_READONLY(core::NumericalRank, basis, dimension, rank, threshold, observation, first_conditioning_warning)
        EPCSAFT_REPR(NumericalRank, "basis", "dimension", "rank", "threshold", "first_conditioning_warning");
    using equilibrium::Quantity; using equilibrium::ComparisonMode; using equilibrium::NumericalEstimate;
    using equilibrium::ComparisonEligibility; using equilibrium::EquilibriumComparison;
    using equilibrium::compare_equilibrium_quantity;
    using equilibrium::project_numerical_estimate;
    nb::enum_<Quantity>(module, "Quantity") EPCSAFT_VALUES(Quantity, Linear, Positive, LogPositive);
    nb::enum_<ComparisonMode>(module, "ComparisonMode") EPCSAFT_VALUES(ComparisonMode, Strict, EstimatedError);
    nb::class_<NumericalEstimate>(module, "NumericalEstimate").def(nb::init<>())
        EPCSAFT_READWRITE(NumericalEstimate, absolute, chart, available, resolved, cause)
        EPCSAFT_REPR(NumericalEstimate, "absolute", "chart", "available", "resolved", "cause");
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
    module.def("property_asl_record", [](const ParameterDocument& document, const HelmholtzModel& model, bool require_phase_identity,
            const std::vector<SolidSpec>& solids, const std::optional<NeutralReferenceDeclaration>& neutral_reference) {
        return asl_property_record(document, model.active_parameters(), model.thermochemistry(), model.solver_options(),
            require_phase_identity, solids, neutral_reference);
    }, nb::arg("document"), nb::arg("model"), nb::arg("require_phase_identity"), nb::arg("solids"),
       nb::arg("neutral_reference").none());
    module.attr("IDEAL_MOLE_FRACTION_STANDARD_STATE_ID") = IDEAL_MOLE_FRACTION_STANDARD_STATE_ID;
    module.attr("IDEAL_MOLARITY_STANDARD_STATE_ID") = IDEAL_MOLARITY_STANDARD_STATE_ID;
    nb::class_<IdealMoleFraction>(module, "IdealMoleFraction")
        .def("__init__", [](IdealMoleFraction* self, std::vector<std::string> ids, std::vector<int> charges, nb::object options,
            std::optional<std::string> basis, std::optional<double> volume, std::optional<std::string> solvent) {
            new (self) IdealMoleFraction(std::move(ids), std::move(charges), bindings::solver_options(options, ROW_SCOPE),
                std::move(basis), volume, std::move(solvent));
        }, nb::arg("component_ids"), nb::arg("charge_number"), nb::arg("solver_options").none() = nb::none(),
           nb::arg("activity_basis").none() = nb::none(), nb::arg("solution_volume_m3").none() = nb::none(),
           nb::arg("unit_activity_solvent").none() = nb::none())
        EPCSAFT_READONLY(IdealMoleFraction, component_ids, charge_number, activity_basis, solution_volume_m3, unit_activity_solvent)
        EPCSAFT_SOLVER_OPTIONS(IdealMoleFraction, ROW_SCOPE);
    nb::class_<SpeciesActivities>(module, "SpeciesActivities")
        EPCSAFT_READONLY(SpeciesActivities, log_activities, activities, causes, solvent_index, standard_molality_mol_per_kg);
    nb::class_<ReactionActivities>(module, "ReactionActivities")
        EPCSAFT_READONLY(ReactionActivities, phase_index, reaction_index, basis, cause, reference_pressure_pa,
            log_quotient, standard_pressure_pa, log_pressure_quotient, log_pressure_quotient_actions, pressure_cause, species);
    module.def("solid_activities", &solid_activities, nb::arg("result"));
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
        .def_static("pinned", &AmountSpec::pinned, nb::arg("value"));
    nb::class_<FreeValue>(module, "FreeValue")
        .def(nb::init<>())
        EPCSAFT_READWRITE(FreeValue, free, value)
        .def_static("fixed", &FreeValue::fixed, nb::arg("value"))
        .def_static("variable", [](std::optional<double> guess) {
            return guess.has_value() ? FreeValue::variable(*guess) : FreeValue::variable();
        }, nb::arg("guess") = nb::none());
    nb::class_<PhaseSpec>(module, "PhaseSpec")
        .def(nb::init<>())
        EPCSAFT_READWRITE(PhaseSpec, name, kind, support, support_declared, amount, composition_guess, packing_guess);
    nb::class_<FeedSpec>(module, "FeedSpec")
        .def(nb::init<>())
        EPCSAFT_READWRITE(FeedSpec, kind, values, solvent, solvent_mass_kg);
    module.def("salt_solute_molalities", &salt_solute_molalities,
        nb::arg("model"), nb::arg("solvent"), nb::arg("salts"));
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
        EPCSAFT_SOLVER_OPTIONS(NeutralReferenceReactionDiagnostic, PROPERTY_SCOPE)
        EPCSAFT_READONLY(NeutralReferenceReactionDiagnostic, reaction_index, terminal_contraction,
            diagnostics, density_convergence, terminal_contraction_pressure_derivative,
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
    nb::enum_<DiscoveryMode>(module, "DiscoveryMode")
        EPCSAFT_VALUES(DiscoveryMode, Full, Fast, FastTested);
    nb::class_<PreparedSuppliedTrials>(module, "PreparedSuppliedTrials")
        .def_ro("temperature", &PreparedSuppliedTrials::temperature)
        .def_prop_ro("work", [](const PreparedSuppliedTrials& s) { return s.work; })
        .def("__len__", [](const PreparedSuppliedTrials& s) { return s.trials.size(); });
    module.def("prepare_supplied_trials", &prepare_supplied_trials, nb::arg("model"), nb::arg("temperature"),
        nb::arg("trials"), nb::keep_alive<0, 1>());
    nb::class_<StabilityOptions>(module, "StabilityOptions")
        .def(nb::init<>())
        EPCSAFT_READWRITE(StabilityOptions, mode, evaluate_supplied_trials, prepared_trials, declaration_budget, trial_budget,
            standard_uncertainty_over_rt, retry_condensed_refusals)
        EPCSAFT_COPIED(StabilityOptions, trials);
    nb::class_<CriticalDeclaration>(module, "CriticalDeclaration")
        .def(nb::init<>())
        EPCSAFT_READWRITE(CriticalDeclaration, density_guess_mol_per_m3, mode_guess);
    nb::class_<SpinodalDeclaration, CriticalDeclaration>(module, "SpinodalDeclaration")
        .def(nb::init<>()) EPCSAFT_READWRITE(SpinodalDeclaration, density_nodes);
    nb::enum_<EnvelopeCoordinate>(module, "EnvelopeCoordinate")
        EPCSAFT_VALUES(EnvelopeCoordinate, LogTemperature, LogPressure, LogK);
    nb::class_<EnvelopeCondition>(module, "EnvelopeCondition")
        .def(nb::init<>())
        EPCSAFT_READWRITE(EnvelopeCondition, kind, component, target);
    nb::class_<EquilibriumProblem>(module, "EquilibriumProblem")
        .def_prop_rw("solver_options", [](const EquilibriumProblem& v) { return bindings::solver_options_view(v.solver_options, ROW_SCOPE); },
            [](EquilibriumProblem& v, nb::object mapping) { v.solver_options = bindings::solver_options(mapping, ROW_SCOPE); }, nb::for_setter(nb::arg("mapping").none()))
        .def_prop_rw("row_tolerance", [](const EquilibriumProblem& v) { return v.solver_options.value(policy::row); },
            [](EquilibriumProblem& v, nb::object value) { bindings::solver_option(v.solver_options, policy::EQUILIBRIUM_ROW_RESIDUAL.name, value, ROW_SCOPE); }, nb::for_setter(nb::arg("value").none()))
        .def(nb::init<>())
        EPCSAFT_READWRITE(EquilibriumProblem, temperature, pressure, feed) EPCSAFT_COPIED(EquilibriumProblem, phases, reactions, solids)
        .def_rw("neutral_reference", &EquilibriumProblem::source_reference)
        EPCSAFT_READWRITE(EquilibriumProblem, standard_gibbs_j_per_mol, topology_declared, enthalpy_j,
            entropy_j_per_k, volume_m3, internal_energy_j, critical, spinodal, envelope, azeotrope, discovery, prescribed_liquid_composition);
    nb::enum_<RowKind>(module, "RowKind")
        EPCSAFT_VALUES(RowKind, Chemical, Pressure, Moiety, Electro, Normalization, PrescribedComposition, Enthalpy, Entropy, InternalEnergy, Volume, CriticalMode,
            CriticalNormalization, CriticalCubic, SpinodalMode, SpinodalNormalization, EnvelopeCondition, AzeotropeCondition, SolidAffinity, ReactionMassAction, ReactionMolarity, SaturationPressure);
    nb::class_<Row>(module, "Row")
        EPCSAFT_READONLY(Row, kind, phase, index)
        EPCSAFT_REPR(Row, "kind", "phase", "index")
        .def_prop_ro("eqid", [](const Row& row) { return row_kind_name(row.kind); });
    nb::class_<PointEvaluation>(module, "PointEvaluation")
        EPCSAFT_READONLY(PointEvaluation, values, objective, floor);
    nb::class_<CompiledSystem>(module, "CompiledSystem")
        .def_prop_ro("rows", [](const CompiledSystem& system) { return system.rows; })
        EPCSAFT_READONLY(CompiledSystem, rank_records, last_rank, pinned_count, free_intensive_count, rank, charge_row_count, primal_size,
            generic_nullity, feed_total, has_ions, has_reaction, conserves_total_moles, jacobian_rank, removed_species)
        .def("solve", [](const CompiledSystem& system) {
            return equilibrium_domains(solve_compiled(system), system.model);
        }, nb::keep_alive<0, 1>())
        .def("solve_saturation", [](const CompiledSystem& system) {
            return equilibrium_domains(solve_saturation(system), system.model);
        }, nb::keep_alive<0, 1>())
        .def("generic_point", &CompiledSystem::generic_point)
        .def("coordinate_labels", &CompiledSystem::coordinate_labels)
        .def("jacobian_difference", &CompiledSystem::jacobian_difference, nb::arg("point"), nb::arg("relative_step"),
            nb::arg("tau"), nb::arg("kappa"))
        .def("numerical_estimates", [](const CompiledSystem& system, std::vector<double> coordinates,
                double temperature, double pressure, std::vector<double> amounts, std::vector<double> composition,
                std::vector<double> densities, std::vector<double> packing, const std::vector<double>& prescribed_errors) {
            const auto phases = system.problem.phases.size(), components = system.component_count();
            if (coordinates.size() != system.unknown_count()
                || amounts.size() != (system.problem.prescribed_liquid_composition ? 0 : phases) || composition.size() != phases * components
                || (!system.ideal && (densities.size() != phases || packing.size() != phases)))
                throw std::invalid_argument("numerical estimate snapshot shape: coordinates " + std::to_string(coordinates.size())
                    + " of " + std::to_string(system.unknown_count()) + ", amounts " + std::to_string(amounts.size())
                    + ", composition " + std::to_string(composition.size()) + " of phases x components "
                    + std::to_string(phases) + " x " + std::to_string(components) + ", densities "
                    + std::to_string(densities.size()) + ", packing " + std::to_string(packing.size()));
            EquilibriumResult result;
            result.success = true; result.temperature = temperature; result.pressure = pressure;
            result.declaration = std::make_shared<const EquilibriumProblem>(system.problem);
            result.coordinates = std::move(coordinates); result.amounts = std::move(amounts);
            result.mole_fractions = std::move(composition); result.molar_densities = std::move(densities);
            result.packing_fractions = std::move(packing);
            return equilibrium_numerical_estimates(system, result, nullptr, false, false, prescribed_errors);
        }, nb::arg("coordinates"), nb::arg("temperature"), nb::arg("pressure"), nb::arg("amounts"),
           nb::arg("composition"), nb::arg("densities"), nb::arg("packing"),
           nb::arg("prescribed_errors") = std::vector<double>{})
        .def("evaluate", [](const CompiledSystem& system, const std::vector<double>& point) {
            return system.evaluate(point); }, nb::arg("point"))
        .def("jacobian", [](const CompiledSystem& system, const std::vector<double>& point) {
            return system.jacobian(point); }, nb::arg("point"))
        .def("lever_jacobian", &CompiledSystem::lever_jacobian, nb::arg("point"))
        .def("objective_gradient", [](const CompiledSystem& system, const std::vector<double>& point) {
            return system.objective_gradient(point); }, nb::arg("point"))
        .def("hessian", [](const CompiledSystem& system, const std::vector<double>& point,
                           const std::vector<double>& multipliers) {
            return system.hessian(point, multipliers); },
            nb::arg("point"), nb::arg("constraint_multipliers") = std::vector<double>{});
    nb::class_<LinearSolveDiagnostics>(module, "LinearSolveDiagnostics")
        EPCSAFT_READONLY(LinearSolveDiagnostics, point_kind, system, reason, residual_norm, evaluation_floor, jacobian_inverse_norm,
            condition_estimate, forward_error_estimate, precision_bits, refinement_iterations);
    nb::class_<IpoptIteration>(module, "IpoptIteration")
        EPCSAFT_READONLY(IpoptIteration, mode, iteration, line_search_trials, objective, inf_pr, inf_du, mu, step_norm, regularization, alpha_du, alpha_pr);
    nb::class_<IpoptIterationRecord>(module, "IpoptIterationRecord") EPCSAFT_READONLY(IpoptIterationRecord, run, iteration);
    nb::class_<IpoptIterations>(module, "IpoptIterations")
        .def(nb::init<>())
        .def("__enter__", &IpoptIterations::enter, nb::rv_policy::reference_internal)
        .def("__exit__", &IpoptIterations::exit)
        .def_ro("records", &IpoptIterations::records);
    nb::class_<IpoptPoint>(module, "IpoptPoint")
        EPCSAFT_READONLY(IpoptPoint, tnlp_coordinates, coordinate_scale, tnlp_lower_multipliers, tnlp_upper_multipliers,
            tnlp_constraint_multipliers, declaration, rows, constraint_rows, objective_factor);
    nb::class_<IpoptDiagnostics>(module, "IpoptDiagnostics")
        EPCSAFT_READONLY(IpoptDiagnostics, last, first_restoration_iteration, initialization_return, application_return,
            restoration_entries, callbacks, warm_initialized, predicted, physical_amounts, failed_point)
        EPCSAFT_REPR(IpoptDiagnostics, "initialization_return", "application_return", "restoration_entries", "callbacks", "warm_initialized");
    nb::enum_<BranchEndCondition>(module, "BranchEndCondition")
        EPCSAFT_VALUES(BranchEndCondition, LoadingFold, VaporMechanicalLimit, VaporKindRefused);
    nb::class_<BubbleBranchEnd>(module, "BubbleBranchEnd")
        EPCSAFT_READONLY(BubbleBranchEnd, condition, feed_bracket_mol, limiting_feed_mol, limiting_pressure_pa);
    nb::class_<EquilibriumResult>(module, "EquilibriumResult", "Phases in declared order; K, Pa, mol, mol/m^3; per-phase-per-component lists are phase-major in component_ids order.", nb::dynamic_attr())
        EPCSAFT_DOMAIN_VIOLATIONS
        .def_prop_ro("critical_search", [](const EquilibriumResult& r) { return r.critical_search ? nb::cast(*r.critical_search) : nb::none(); })
        .def_prop_ro("rejected_reduction", [](const EquilibriumResult& r) { return r.rejected_reduction ? nb::cast(*r.rejected_reduction) : nb::none(); })
        .def_prop_ro("recovery_stability", [](const EquilibriumResult& r) { return r.recovery_stability ? nb::cast(*r.recovery_stability) : nb::none(); })
        .def_prop_ro("vanished_phase_stability", [](const EquilibriumResult& r) { return r.vanished_phase_stability ? nb::cast(*r.vanished_phase_stability) : nb::none(); })
        .def_prop_ro("solver_options", [](const EquilibriumResult& v) { return v.declaration ? bindings::solver_options_view(v.declaration->solver_options, PROPERTY_SCOPE) : nb::list(); })
        EPCSAFT_READONLY(EquilibriumResult, solver_settings, diagnostic_options_file, linear_solver, first_numerical_failure, rank_records, ipopt_diagnostics, recovered_ipopt_diagnostics)
        .def_ro("neutral_reference_diagnostics", &EquilibriumResult::reference_diagnostics)
        .def_ro("branch_end", &EquilibriumResult::branch_end)
        .def_ro("electric_potential_difference_v", &EquilibriumResult::electric_potential_difference_v)
        EPCSAFT_READONLY(EquilibriumResult, terminal_numerical_failure, last_full_evaluation, phase_pressures, discovery, success, message, solver_status, temperature, pressure, amounts,
            mole_fractions, molar_densities, packing_fractions,
            // Dimensionless reduced-energy balance/charge stationarity coordinates of the solution.
            multipliers,
            // Extensive-objective multipliers of the non-chemical, non-pressure `rows` in order, formed at
            // the returned point from `multipliers`; raw optimizer duals only when a phase amount is pinned.
            constraint_multipliers, residuals, classified_residuals, coordinates, phase_count, iterations, evaluations,
            hessian_calls, objective, reduced_hessian_minimum, convergence, last_rank, linear_solve_diagnostics, numerical_estimates, rows, unknown_count, jacobian_rank,
            generic_nullity, requested_tolerance_met, topology_event, topology_phase, topology_pair, critical, spinodal,
            // Per declared solid: amount in mol and saturation index SI_s in decades.
            solid_amounts, solid_affinity_decades)
        EPCSAFT_REPR(EquilibriumResult, "success", "temperature", "pressure", "amounts", "mole_fractions")
        .def("phase_mole_fractions", [](const EquilibriumResult& r, std::ptrdiff_t phase) {  // Python sequence indexing
            const auto count = static_cast<std::ptrdiff_t>(r.phase_count);
            if (phase < -count || phase >= count) throw nb::index_error("phase outside the declared phases");
            const auto n = static_cast<std::ptrdiff_t>(r.mole_fractions.size()) / count;
            const auto first = (phase < 0 ? phase + count : phase) * n;
            return std::vector<double>(r.mole_fractions.begin() + first, r.mole_fractions.begin() + first + n);
        }, nb::arg("phase"));
    nb::class_<HelmholtzModeRecord>(module, "HelmholtzModeRecord")
        EPCSAFT_READONLY(HelmholtzModeRecord, eigenvalues, mode, zero_band, condition_estimate, backward_error, rank);
    nb::class_<CriticalRecord, HelmholtzModeRecord>(module, "CriticalRecord")
        EPCSAFT_READONLY(CriticalRecord, cubic);
    nb::class_<SpinodalBracket>(module, "SpinodalBracket")
        EPCSAFT_READONLY(SpinodalBracket, lower, upper, lower_eigenvalue, upper_eigenvalue);
    nb::class_<SpinodalSearch>(module, "SpinodalSearch")
        EPCSAFT_READONLY(SpinodalSearch, brackets, finite_limits, unresolved_densities, lower, upper, seed, geometric_cap, evaluation_order, evaluations);
    nb::class_<SpinodalRecord, HelmholtzModeRecord>(module, "SpinodalRecord")
        EPCSAFT_READONLY(SpinodalRecord, search, normalized_null_residual);
    nb::class_<CriticalSearch>(module, "CriticalSearch")
        EPCSAFT_READONLY(CriticalSearch, seeds, results, added_work, scheduled, budget_exhausted, selected);
    nb::class_<StabilityPhaseComparison>(module, "StabilityPhaseComparison")
        EPCSAFT_READONLY(StabilityPhaseComparison, phase, added_species, missing_species, composition_max_distance, volume_log_distance);
    nb::class_<StabilityTrial>(module, "StabilityTrial")
        EPCSAFT_READONLY(StabilityTrial, phase, seed, evidence, message, optimizer, raw_coordinates, raw_residuals,
            composition, split_coordinates, basis, row_transformation, quantities, split_checks, present_phase_comparisons, diagnostic_cause, donor,
            resolved_descent);
    nb::class_<StabilityReport>(module, "StabilityReport")
        EPCSAFT_READONLY(StabilityReport, first_numerical_failure, last_rank, completed, message, trials, candidate_coordinates, standard_over_rt,
            plane_over_rt, present_phases, work, no_descent_in_trial_set);
    nb::class_<StabilityMeasure>(module, "StabilityMeasure")
        EPCSAFT_READONLY(StabilityMeasure, report, raw_minimum, value, uncertainty, minimizing_trial, minima, cause, event, available);
    nb::class_<DiscoveryRecord>(module, "DiscoveryRecord")
        EPCSAFT_READONLY(DiscoveryRecord, mode, termination, unsearched, declaration, candidate, rounds, declarations, declaration_roles,
            declaration_attempt_coordinates, work);
    module.def("stability_measure", &stability_measure, nb::arg("compiled"), nb::arg("result"), nb::arg("options") = StabilityOptions{});
    module.def("cooling_onset", &cooling_onset, nb::arg("compiled"));
    module.def("stability_test", &stability_test, nb::arg("compiled"), nb::arg("central"),
        nb::arg("options") = StabilityOptions{});
    nb::enum_<TopologyEvent>(module, "TopologyEvent")
        .value("None_", TopologyEvent::None)
        EPCSAFT_VALUES(TopologyEvent, PhaseDisappeared, PhasesCoalesced, BranchEnd);
    auto continuation_stops = nb::enum_<ContinuationStop>(module, "ContinuationStop");
    for (int value = static_cast<int>(ContinuationStop::reached_target);
         value <= static_cast<int>(ContinuationStop::turning_point); ++value) {
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
    nb::class_<ContinuationResult>(module, "ContinuationResult", nb::dynamic_attr())
        EPCSAFT_DOMAIN_VIOLATIONS
        .def_prop_ro("first_numerical_failure", [](const ContinuationResult& r) { return r.observations.first_numerical_failure; })
        .def_prop_ro("rank_records", [](const ContinuationResult& r) { return r.observations.rank_records; })
        .def_prop_ro("ipopt_diagnostics", [](const ContinuationResult& r) { return r.observations.ipopt_diagnostics; })
        .def_prop_ro("recovered_ipopt_diagnostics", [](const ContinuationResult& r) { return r.observations.recovered_ipopt_diagnostics; })
        EPCSAFT_READONLY(ContinuationResult, success, message, stop, path_values, parameter_values, states,
            reached_path, attempted_path, event_phase, event_pair, attempts, accepted_steps, retries, last_result,
            cause_stop, cause_convergence, cause_message, conditions, critical_points, critical_seeds)
        EPCSAFT_REPR(ContinuationResult, "success", "stop", "reached_path", "accepted_steps", "message");
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
    nb::class_<EquilibriumStartsResult>(module, "EquilibriumStartsResult", nb::dynamic_attr())
        EPCSAFT_DOMAIN_VIOLATIONS
        EPCSAFT_READONLY(EquilibriumStartsResult, selected, results, causes, attempt_seconds, iterations, seconds,
            overrun_seconds, budget_exhausted);
    module.def("solve_equilibrium_starts", [](const HelmholtzModel& model, const EquilibriumProblem& problem,
        const std::vector<EquilibriumStart>& starts, const EquilibriumWorkLimits& limits) {
        return bindings::domain_result(solve_equilibrium_starts(model, problem, starts, limits), &model,
            [&](auto& rows, const EquilibriumStartsResult& r) {
                for (std::size_t i = 0; i < r.results.size(); ++i) if (r.results[i]) phase_domains(rows, model, *r.results[i], i);
            });
    }, nb::arg("model"), nb::arg("problem"), nb::arg("starts"), nb::arg("limits"));
    module.def("continue_equilibrium", [](const HelmholtzModel& model, const ContinuationRequest& request) {
        return bindings::domain_result(continue_equilibrium(model, request), &model, [&](auto& rows, const ContinuationResult& r) {
            for (std::size_t i = 0; i < r.states.size(); ++i) phase_domains(rows, model, r.states[i], i);
            for (std::size_t i = 0; i < r.critical_points.size(); ++i) phase_domains(rows, model, r.critical_points[i], r.states.size() + i);
            if (r.last_result) phase_domains(rows, model, *r.last_result, r.states.size() + r.critical_points.size());
        });
    }, nb::arg("model"), nb::arg("request"));
    // Solved-state actions: the existing native records and owner, transported unchanged.
    nb::enum_<SolvedStateObservableKind>(module, "SolvedStateObservableKind")
        EPCSAFT_VALUES(SolvedStateObservableKind, Temperature, Pressure, PhaseAmount, PhaseComponentAmount,
            PhaseMoleFraction, PhaseVolume, PhaseMolarDensity, PhasePackingFraction, PhasePressure,
            PhaseResidualChemicalPotentialOverRt, PhaseProperty, PhaseLogActivity, PhaseGibbsMolarVolume, SystemEnthalpy,
            EquilibriumIsobaricHeatCapacity, SolidAmount, SolidSaturationIndex);
    nb::class_<SolvedStateObservable>(module, "SolvedStateObservable")
        .def(nb::init<SolvedStateObservableKind, std::size_t, std::size_t, PropertyObservable>(), nb::arg("kind"), nb::arg("phase") = 0, nb::arg("component") = 0,
            nb::arg("property") = PropertyObservable::IdealEnthalpy)
        EPCSAFT_READONLY(SolvedStateObservable, kind, phase, component, property)
        EPCSAFT_REPR(SolvedStateObservable, "kind", "phase", "component", "property");
    nb::class_<SolvedStateActionDirection>(module, "SolvedStateActionDirection")
        .def(nb::init<double, double, std::vector<double>, std::vector<double>, std::vector<std::array<double, 4>>, double, double, double, double>(),
             nb::arg("temperature_k"), nb::arg("pressure_pa"), nb::arg("feed_amounts_mol"),
             nb::arg("active_parameters"),
             nb::arg("reaction_coefficients") = std::vector<std::array<double, 4>>{},
             nb::arg("enthalpy_j") = 0.0, nb::arg("entropy_j_per_k") = 0.0, nb::arg("volume_m3") = 0.0, nb::arg("internal_energy_j") = 0.0)
        EPCSAFT_READONLY(SolvedStateActionDirection, temperature_k, pressure_pa, feed_amounts_mol, active_parameters,
            reaction_coefficients, enthalpy_j, entropy_j_per_k, volume_m3, internal_energy_j);
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
    nb::class_<SolvedStateActionBatch>(module, "SolvedStateActionBatch", nb::dynamic_attr())
        EPCSAFT_DOMAIN_VIOLATIONS
        EPCSAFT_READONLY(SolvedStateActionBatch, results, diagnostics);
    module.def("solved_state_difference", &solved_state_difference, nb::arg("compiled"), nb::arg("central"),
        nb::arg("observables"), nb::arg("direction"), nb::arg("step"), nb::arg("tau"), nb::arg("kappa"));
    module.def("solved_state_actions", [](const CompiledSystem& system, const EquilibriumResult& central,
        const SolvedStateActionRequest& request) {
        return bindings::domain_result(solved_state_actions(system, central, request), system.model,
            [&](auto& rows, const auto&) { phase_domains(rows, *system.model, central); });
    },
               nb::arg("compiled"), nb::arg("central"), nb::arg("request"),
               "Exact order 1-3 actions at `central`, the accepted solve of this same compiled problem.\n\n"
               "Directions carry K, Pa, absolute feed mol (component width) and active-parameter increments\n"
               "(active width, native units); a free T or P entry must be zero. No solve is repeated.");
    module.def("solved_state_first_actions", [](const CompiledSystem& system, const EquilibriumResult& central,
        const std::vector<SolvedStateObservable>& observables, const std::vector<SolvedStateActionDirection>& directions) {
        return bindings::domain_result(solved_state_first_actions(system, central, observables, directions), system.model,
            [&](auto& rows, const auto& r) { if (!r.empty()) phase_domains(rows, *system.model, central); });
    }, nb::arg("compiled"), nb::arg("central"), nb::arg("observables"), nb::arg("directions"));
    nb::class_<InventoryBound>(module, "InventoryBound")
        EPCSAFT_READONLY(InventoryBound, species, status, termination, maximum, feed_normalized_amount, log_fraction,
            balance_residual, normalization_residual, lower_bound_residual, witness, amounts_mol, finite_witness,
            limiting, supported, floor, lp_error, finite_residual, feasibility_termination);
    module.def("inventory_bounds", &inventory_bounds, nb::arg("compiled"), nb::arg("species"),
        nb::arg("fraction"), nb::arg("floor"));
    module.def("equilibrium_asl_record", [](const ParameterDocument& document, const HelmholtzModel& model, const EquilibriumProblem& problem,
            const ContinuationOptions& continuation, std::vector<double> offset, std::vector<std::vector<double>> matrix,
            std::vector<double> anchor, std::size_t identity_component, double minimum_gap) {
        AslFeedMap map{std::move(offset), std::move(anchor), std::move(matrix), identity_component, minimum_gap};
        map.validate(problem, model.parameters());
        compile_problem(model, problem);
        return asl_record(document, model.active_parameters(), model.thermochemistry(), problem, continuation, map, model.solver_options());
    }, nb::arg("document"), nb::arg("model"), nb::arg("problem"), nb::arg("continuation"), nb::arg("feed_offset"),
       nb::arg("feed_matrix"), nb::arg("anchor_controls"), nb::arg("identity_component"), nb::arg("minimum_gap"));
    module.def("linear_solver", [] { return linear_solver(); });  // the process-level record
    module.def("compile_problem", nb::overload_cast<const HelmholtzModel&, const EquilibriumProblem&, const std::vector<std::string>&>(&compile_problem),
        nb::arg("model"), nb::arg("problem"), nb::arg("component_ids") = std::vector<std::string>{}, nb::keep_alive<0, 1>());
    module.def("compile_problem", nb::overload_cast<const IdealMoleFraction&, const EquilibriumProblem&>(&compile_problem), nb::arg("model"),
        nb::arg("problem"), nb::keep_alive<0, 1>());
    module.def("sensitivity", [](const CompiledSystem& compiled, const std::vector<double>& point, const std::vector<double>& lever) {
        std::optional<core::Convergence> failure;
        auto result = sensitivity(compiled, point, lever, &failure);
        if (failure) throw core::ConvergenceFailure(*failure);
        return result;
    }, nb::arg("compiled"), nb::arg("point"), nb::arg("lever"));
    module.def("solve_equilibrium", [](const HelmholtzModel& model, const EquilibriumProblem& problem,
        const std::optional<EquilibriumResult>& start, const std::vector<std::string>& ids,
        std::optional<std::size_t> work_limit, std::optional<std::size_t> branch_end_work_limit,
        const CompiledSystem* prepared) {
        auto result = prepared && !work_limit && !branch_end_work_limit
            ? solve_equilibrium(model, problem, start, ids, *prepared)
            : solve_equilibrium(model, problem, start, ids, work_limit, branch_end_work_limit);
        return equilibrium_domains(std::move(result), &model);
    },
        nb::arg("model"), nb::arg("problem"), nb::arg("start") = std::nullopt, nb::arg("component_ids") = std::vector<std::string>{},
        nb::arg("work_limit") = std::nullopt, nb::kw_only(), nb::arg("branch_end_work_limit") = std::nullopt,
        nb::arg("prepared").none() = nb::none(), nb::keep_alive<0, 1>());
    module.def("solve_equilibrium", nb::overload_cast<const IdealMoleFraction&, const EquilibriumProblem&, const std::optional<EquilibriumResult>&>(&solve_equilibrium),
        nb::arg("model"), nb::arg("problem"), nb::arg("start") = std::nullopt, nb::keep_alive<0, 1>());
}

}  // namespace epcsaft::equilibrium
