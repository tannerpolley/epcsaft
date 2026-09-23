#include "equilibrium/problem.hpp"
#include "equilibrium/gibbs_minimizer.hpp"
#include "equilibrium_record.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <optional>

namespace py = pybind11;
namespace epcsaft::equilibrium {

void register_equilibrium(py::module_& module) {
    module.def("equilibrium_asl_record", [](py::object model, py::object problem, py::object continuation,
            std::vector<double> offset, std::vector<std::vector<double>> matrix, std::vector<double> anchor,
            std::size_t identity_component, double minimum_gap) {
        const py::object native_owner = model.attr("_model");
        const auto& native = native_owner.cast<const HelmholtzModel&>();
        if (!native.active_parameters().empty() || native.thermochemistry())
            throw std::invalid_argument("ASL record requires inactive parameters and explicit problem reaction references");
        const auto document = model.attr("parameters").attr("_document").cast<ParameterDocument>();
        const auto p = problem.attr("_native")(model).cast<EquilibriumProblem>();
        AslFeedMap map{std::move(offset), std::move(anchor), std::move(matrix), identity_component, minimum_gap};
        map.validate(p, native.parameters());
        compile_problem(native, p);
        ContinuationOptions options;
        options.initial_step = continuation.attr("initial_step").cast<double>();
        options.minimum_step = continuation.attr("minimum_step").cast<double>();
        options.max_steps = continuation.attr("max_steps").cast<std::size_t>();
        options.max_retries_per_step = continuation.attr("max_retries_per_step").cast<std::size_t>();
        return asl_record(document, p, options, map);
    }, py::arg("model"), py::arg("problem"), py::arg("continuation"), py::arg("feed_offset"),
       py::arg("feed_matrix"), py::arg("anchor_controls"), py::arg("identity_component"), py::arg("minimum_gap"));
    module.attr("EOS_STANDARD_STATE_ID") = EOS_STANDARD_STATE_ID;
    py::enum_<PhaseKind>(module, "PhaseKind")
        .value("Liquid", PhaseKind::Liquid)
        .value("Vapor", PhaseKind::Vapor);
    py::enum_<AmountKind>(module, "AmountKind")
        .value("Free", AmountKind::Free)
        .value("Pinned", AmountKind::Pinned);
    py::enum_<FeedKind>(module, "FeedKind")
        .value("MoleFractions", FeedKind::MoleFractions)
        .value("Amounts", FeedKind::Amounts)
        .value("MassFractions", FeedKind::MassFractions)
        .value("Molalities", FeedKind::Molalities);
    py::class_<AmountSpec>(module, "AmountSpec")
        .def(py::init<>())
        .def_readwrite("kind", &AmountSpec::kind)
        .def_readwrite("fraction", &AmountSpec::fraction)
        .def_static("free", &AmountSpec::free)
        .def_static("pinned", &AmountSpec::pinned);
    py::class_<FreeValue>(module, "FreeValue")
        .def(py::init<>())
        .def_readwrite("free", &FreeValue::free)
        .def_readwrite("value", &FreeValue::value)
        .def_static("fixed", &FreeValue::fixed)
        .def_static("variable", [](std::optional<double> guess) {
            return guess.has_value() ? FreeValue::variable(*guess) : FreeValue::variable();
        }, py::arg("guess") = py::none());
    py::class_<PhaseSpec>(module, "PhaseSpec")
        .def(py::init<>())
        .def_readwrite("name", &PhaseSpec::name)
        .def_readwrite("kind", &PhaseSpec::kind)
        .def_readwrite("support", &PhaseSpec::support)
        .def_readwrite("support_declared", &PhaseSpec::support_declared)
        .def_readwrite("amount", &PhaseSpec::amount)
        .def_readwrite("composition_guess", &PhaseSpec::composition_guess)
        .def_readwrite("packing_guess", &PhaseSpec::packing_guess);
    py::class_<FeedSpec>(module, "FeedSpec")
        .def(py::init<>())
        .def_readwrite("kind", &FeedSpec::kind)
        .def_readwrite("values", &FeedSpec::values)
        .def_readwrite("solvent", &FeedSpec::solvent)
        .def_readwrite("solvent_mass_kg", &FeedSpec::solvent_mass_kg);
    py::class_<NeutralReferenceDeclaration>(module, "NeutralReferenceDeclaration")
        .def(py::init<>())
        .def_readwrite("solvent_index", &NeutralReferenceDeclaration::solvent_index)
        .def_readwrite("solute_molality_ratios", &NeutralReferenceDeclaration::solute_molality_ratios)
        .def_readwrite("standard_molality_mol_per_kg", &NeutralReferenceDeclaration::standard_molality_mol_per_kg);
    py::enum_<ReferenceSourceBasis>(module, "ReferenceSourceBasis")
        .value("CommonMolalityInfiniteDilution", ReferenceSourceBasis::CommonMolalityInfiniteDilution)
        .value("RawMoleFractionInfiniteDilution", ReferenceSourceBasis::RawMoleFractionInfiniteDilution);
    py::enum_<NeutralReferenceStatus>(module, "NeutralReferenceStatus")
        .value("Available", NeutralReferenceStatus::Available)
        .value("InvalidDeclaration", NeutralReferenceStatus::InvalidDeclaration)
        .value("NonNeutralDeclaration", NeutralReferenceStatus::NonNeutralDeclaration)
        .value("UnsupportedSourceBasis", NeutralReferenceStatus::UnsupportedSourceBasis)
        .value("FreeTemperature", NeutralReferenceStatus::FreeTemperature)
        .value("FreePressure", NeutralReferenceStatus::FreePressure)
        .value("DensityUnavailable", NeutralReferenceStatus::DensityUnavailable)
        .value("BranchChange", NeutralReferenceStatus::BranchChange)
        .value("Nonfinite", NeutralReferenceStatus::Nonfinite)
        .value("TerminalChangeExceeded", NeutralReferenceStatus::TerminalChangeExceeded);
    py::class_<NeutralReferenceReactionDiagnostic>(module, "NeutralReferenceReactionDiagnostic")
        .def_readonly("reaction_index", &NeutralReferenceReactionDiagnostic::reaction_index)
        .def_readonly("terminal_contraction", &NeutralReferenceReactionDiagnostic::terminal_contraction)
        .def_readonly("observed_terminal_change", &NeutralReferenceReactionDiagnostic::observed_terminal_change)
        .def_readonly("terminal_contraction_pressure_derivative",
                      &NeutralReferenceReactionDiagnostic::terminal_contraction_pressure_derivative)
        .def_readonly("terminal_contraction_pressure_second_derivative",
                      &NeutralReferenceReactionDiagnostic::terminal_contraction_pressure_second_derivative)
        .def_readonly("reference_pressure_pa", &NeutralReferenceReactionDiagnostic::reference_pressure_pa)
        .def_readonly("status", &NeutralReferenceReactionDiagnostic::status)
        .def_readonly("message", &NeutralReferenceReactionDiagnostic::message);
    py::class_<ReactionReference>(module, "ReactionReference")
        .def(py::init([](ReferenceSourceBasis basis, std::optional<double> pressure) {
            return ReactionReference{pressure, basis};
        }), py::arg("source_basis"), py::arg("reference_pressure_pa") = py::none())
        .def_readwrite("source_basis", &ReactionReference::source_basis)
        .def_readwrite("reference_pressure_pa", &ReactionReference::reference_pressure_pa);
    py::class_<ReactionSpec>(module, "ReactionSpec")
        .def(py::init<>())
        .def_readwrite("name", &ReactionSpec::name)
        .def_readwrite("stoichiometry", &ReactionSpec::stoichiometry)
        .def_readwrite("equilibrium_constant", &ReactionSpec::equilibrium_constant)
        .def_readwrite("affinity_j_per_mol", &ReactionSpec::affinity_j_per_mol)
        .def_readwrite("correlation", &ReactionSpec::correlation)
        .def_readwrite("reference", &ReactionSpec::reference);
    py::class_<EquilibriumProblem>(module, "EquilibriumProblem")
        .def(py::init<>())
        .def_readwrite("phases", &EquilibriumProblem::phases)
        .def_readwrite("temperature", &EquilibriumProblem::temperature)
        .def_readwrite("pressure", &EquilibriumProblem::pressure)
        .def_readwrite("feed", &EquilibriumProblem::feed)
        .def_readwrite("reactions", &EquilibriumProblem::reactions)
        .def_readwrite("neutral_reference", &EquilibriumProblem::source_reference)
        .def_readwrite("standard_gibbs_j_per_mol", &EquilibriumProblem::standard_gibbs_j_per_mol)
        .def_readwrite("topology_declared", &EquilibriumProblem::topology_declared)
        .def_readwrite("row_tolerance", &EquilibriumProblem::row_tolerance);
    py::enum_<RowKind>(module, "RowKind")
        .value("Chemical", RowKind::Chemical)
        .value("Pressure", RowKind::Pressure)
        .value("Moiety", RowKind::Moiety)
        .value("Electro", RowKind::Electro)
        .value("Normalization", RowKind::Normalization);
    py::class_<Row>(module, "Row")
        .def_readonly("kind", &Row::kind)
        .def_readonly("phase", &Row::phase)
        .def_readonly("index", &Row::index)
        .def_readonly("eqid", &Row::eqid);
    py::class_<PointEvaluation>(module, "PointEvaluation")
        .def_readonly("values", &PointEvaluation::values)
        .def_readonly("jacobian", &PointEvaluation::jacobian)
        .def_readonly("objective", &PointEvaluation::objective)
        .def_readonly("floor", &PointEvaluation::floor);
    py::class_<CompiledSystem>(module, "CompiledSystem")
        .def_property_readonly("rows", [](const CompiledSystem& system) { return system.rows; })
        .def_readonly("pinned_count", &CompiledSystem::pinned_count)
        .def_readonly("free_intensive_count", &CompiledSystem::free_intensive_count)
        .def_readonly("rank", &CompiledSystem::rank)
        .def_readonly("charge_row_count", &CompiledSystem::charge_row_count)
        .def_readonly("primal_size", &CompiledSystem::primal_size)
        .def_readonly("rows_minus_unknowns", &CompiledSystem::rows_minus_unknowns)
        .def_readonly("generic_nullity", &CompiledSystem::generic_nullity)
        .def_readonly("feed_total", &CompiledSystem::feed_total)
        .def_readonly("has_ions", &CompiledSystem::has_ions)
        .def_readonly("has_reaction", &CompiledSystem::has_reaction)
        .def_readonly("conserves_total_moles", &CompiledSystem::conserves_total_moles)
        .def_readonly("jacobian_rank", &CompiledSystem::jacobian_rank)
        .def_readonly("removed_species", &CompiledSystem::removed_species)
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
    py::class_<EquilibriumResult>(module, "EquilibriumResult")
        .def_readonly("neutral_reference_diagnostics", &EquilibriumResult::reference_diagnostics)
        .def_readonly("success", &EquilibriumResult::success)
        .def_readonly("message", &EquilibriumResult::message)
        .def_readonly("solver_status", &EquilibriumResult::solver_status)
        .def_readonly("temperature", &EquilibriumResult::temperature)
        .def_readonly("pressure", &EquilibriumResult::pressure)
        .def_readonly("amounts", &EquilibriumResult::amounts)
        .def_readonly("mole_fractions", &EquilibriumResult::mole_fractions)
        .def_readonly("molar_densities", &EquilibriumResult::molar_densities)
        .def_readonly("packing_fractions", &EquilibriumResult::packing_fractions)
        // Dimensionless reduced-energy balance/charge stationarity coordinates of the solution.
        .def_readonly("multipliers", &EquilibriumResult::multipliers)
        // Extensive-objective multipliers of the non-chemical, non-pressure `rows` in order, formed at
        // the returned point from `multipliers`; raw optimizer duals only when a phase amount is pinned.
        .def_readonly("constraint_multipliers", &EquilibriumResult::constraint_multipliers)
        .def_readonly("residuals", &EquilibriumResult::residuals)
        .def_readonly("classified_residuals", &EquilibriumResult::classified_residuals)
        .def_readonly("sensitivity", &EquilibriumResult::sensitivity)
        .def_readonly("coordinates", &EquilibriumResult::coordinates)
        .def_readonly("phase_count", &EquilibriumResult::phase_count)
        .def_readonly("iterations", &EquilibriumResult::iterations)
        .def_readonly("evaluations", &EquilibriumResult::evaluations)
        .def_readonly("hessian_calls", &EquilibriumResult::hessian_calls)
        .def_readonly("objective", &EquilibriumResult::objective)
        .def_readonly("reduced_hessian_minimum", &EquilibriumResult::reduced_hessian_minimum)
        .def_readonly("convergence", &EquilibriumResult::convergence)
        .def_readonly("rows", &EquilibriumResult::rows)
        .def_readonly("unknown_count", &EquilibriumResult::unknown_count)
        .def_readonly("jacobian_rank", &EquilibriumResult::jacobian_rank)
        .def_readonly("generic_nullity", &EquilibriumResult::generic_nullity)
        .def_readonly("requested_tolerance_met", &EquilibriumResult::requested_tolerance_met)
        .def_readonly("topology_event", &EquilibriumResult::topology_event)
        .def_readonly("topology_phase", &EquilibriumResult::topology_phase)
        .def_readonly("topology_pair", &EquilibriumResult::topology_pair);
    py::enum_<TopologyEvent>(module, "TopologyEvent")
        .value("None_", TopologyEvent::None)
        .value("PhaseDisappeared", TopologyEvent::PhaseDisappeared)
        .value("PhasesCoalesced", TopologyEvent::PhasesCoalesced);
    py::enum_<ContinuationStop>(module, "ContinuationStop")
        .value("reached_target", ContinuationStop::reached_target)
        .value("initial_failure", ContinuationStop::initial_failure)
        .value("invalid_path", ContinuationStop::invalid_path)
        .value("singular_jacobian", ContinuationStop::singular_jacobian)
        .value("branch_lost", ContinuationStop::branch_lost)
        .value("phase_disappeared", ContinuationStop::phase_disappeared)
        .value("phases_coalesced", ContinuationStop::phases_coalesced)
        .value("left_domain", ContinuationStop::left_domain)
        .value("infeasible", ContinuationStop::infeasible)
        .value("nonfinite", ContinuationStop::nonfinite)
        .value("stalled", ContinuationStop::stalled)
        .value("step_too_small", ContinuationStop::step_too_small)
        .value("budget_exhausted", ContinuationStop::budget_exhausted);
    py::class_<ContinuationOptions>(module, "ContinuationOptions")
        .def(py::init<>())
        .def_readwrite("initial_step", &ContinuationOptions::initial_step)
        .def_readwrite("minimum_step", &ContinuationOptions::minimum_step)
        .def_readwrite("max_steps", &ContinuationOptions::max_steps)
        .def_readwrite("max_retries_per_step", &ContinuationOptions::max_retries_per_step);
    py::class_<ContinuationRequest>(module, "ContinuationRequest")
        .def(py::init<>())
        .def_readwrite("start_problem", &ContinuationRequest::start_problem)
        .def_readwrite("target_problem", &ContinuationRequest::target_problem)
        .def_readwrite("start_result", &ContinuationRequest::start_result)
        .def_readwrite("options", &ContinuationRequest::options);
    py::class_<ContinuationResult>(module, "ContinuationResult")
        .def_readonly("success", &ContinuationResult::success)
        .def_readonly("message", &ContinuationResult::message)
        .def_readonly("stop", &ContinuationResult::stop)
        .def_readonly("path_values", &ContinuationResult::path_values)
        .def_readonly("states", &ContinuationResult::states)
        .def_readonly("reached_path", &ContinuationResult::reached_path)
        .def_readonly("attempted_path", &ContinuationResult::attempted_path)
        .def_readonly("event_phase", &ContinuationResult::event_phase)
        .def_readonly("event_pair", &ContinuationResult::event_pair)
        .def_readonly("attempts", &ContinuationResult::attempts)
        .def_readonly("accepted_steps", &ContinuationResult::accepted_steps)
        .def_readonly("retries", &ContinuationResult::retries)
        .def_readonly("last_result", &ContinuationResult::last_result)
        .def_readonly("cause_stop", &ContinuationResult::cause_stop)
        .def_readonly("cause_convergence", &ContinuationResult::cause_convergence)
        .def_readonly("cause_message", &ContinuationResult::cause_message);
    module.def("continue_equilibrium", &continue_equilibrium, py::arg("model"), py::arg("request"));
    py::class_<SensitivityInfo>(module, "SensitivityInfo")
        .def_readonly("equations", &SensitivityInfo::equations)
        .def_readonly("parameters", &SensitivityInfo::parameters)
        .def_readonly("rank", &SensitivityInfo::rank)
        .def_readonly("condition_estimate", &SensitivityInfo::condition_estimate)
        .def_readonly("nonsingular", &SensitivityInfo::nonsingular);
    // E5 solved-state actions: the existing native records and owner, transported unchanged.
    py::enum_<SolvedStateObservableKind>(module, "SolvedStateObservableKind")
        .value("Temperature", SolvedStateObservableKind::Temperature)
        .value("Pressure", SolvedStateObservableKind::Pressure)
        .value("PhaseAmount", SolvedStateObservableKind::PhaseAmount)
        .value("PhaseComponentAmount", SolvedStateObservableKind::PhaseComponentAmount)
        .value("PhaseMoleFraction", SolvedStateObservableKind::PhaseMoleFraction)
        .value("PhaseVolume", SolvedStateObservableKind::PhaseVolume)
        .value("PhaseMolarDensity", SolvedStateObservableKind::PhaseMolarDensity)
        .value("PhasePackingFraction", SolvedStateObservableKind::PhasePackingFraction)
        .value("PhasePressure", SolvedStateObservableKind::PhasePressure)
        .value("PhaseResidualChemicalPotentialOverRt",
               SolvedStateObservableKind::PhaseResidualChemicalPotentialOverRt)
        .value("PhaseProperty", SolvedStateObservableKind::PhaseProperty)
        .value("PhaseLogActivity", SolvedStateObservableKind::PhaseLogActivity);
    py::class_<SolvedStateObservable>(module, "SolvedStateObservable")
        .def(py::init([](SolvedStateObservableKind kind, std::size_t phase, std::size_t component,
                         PropertyObservable property) {
            return SolvedStateObservable{kind, phase, component, property};
        }), py::arg("kind"), py::arg("phase") = 0, py::arg("component") = 0,
            py::arg("property") = PropertyObservable::IdealEnthalpy)
        .def_readonly("kind", &SolvedStateObservable::kind)
        .def_readonly("phase", &SolvedStateObservable::phase)
        .def_readonly("component", &SolvedStateObservable::component)
        .def_readonly("property", &SolvedStateObservable::property);
    py::class_<SolvedStateActionDirection>(module, "SolvedStateActionDirection")
        .def(py::init([](double temperature_k, double pressure_pa, std::vector<double> feed_amounts_mol,
                         std::vector<double> active_parameters) {
            return SolvedStateActionDirection{temperature_k, pressure_pa, std::move(feed_amounts_mol),
                                              std::move(active_parameters)};
        }), py::arg("temperature_k"), py::arg("pressure_pa"), py::arg("feed_amounts_mol"),
            py::arg("active_parameters"))
        .def_readonly("temperature_k", &SolvedStateActionDirection::temperature_k)
        .def_readonly("pressure_pa", &SolvedStateActionDirection::pressure_pa)
        .def_readonly("feed_amounts_mol", &SolvedStateActionDirection::feed_amounts_mol)
        .def_readonly("active_parameters", &SolvedStateActionDirection::active_parameters);
    py::class_<SolvedStateActionRequest>(module, "SolvedStateActionRequest")
        .def(py::init([](std::vector<SolvedStateObservable> observables,
                         std::vector<SolvedStateActionDirection> directions) {
            return SolvedStateActionRequest{std::move(observables), std::move(directions)};
        }), py::arg("observables"), py::arg("directions"))
        .def_readonly("observables", &SolvedStateActionRequest::observables)
        .def_readonly("directions", &SolvedStateActionRequest::directions);
    py::enum_<SolvedStateActionStatus>(module, "SolvedStateActionStatus")
        .value("Available", SolvedStateActionStatus::Available)
        .value("CentralSolveUnavailable", SolvedStateActionStatus::CentralSolveUnavailable)
        .value("DefiningEquationsNotSatisfied", SolvedStateActionStatus::DefiningEquationsNotSatisfied)
        .value("PhaseBoundaryNotDifferentiable", SolvedStateActionStatus::PhaseBoundaryNotDifferentiable)
        .value("NonuniqueCoalescedPhases", SolvedStateActionStatus::NonuniqueCoalescedPhases)
        .value("RankDeficient", SolvedStateActionStatus::RankDeficient)
        .value("NonfiniteJacobian", SolvedStateActionStatus::NonfiniteJacobian)
        .value("LinearSolveResidual", SolvedStateActionStatus::LinearSolveResidual)
        .value("SupportChange", SolvedStateActionStatus::SupportChange)
        .value("IonicInsertion", SolvedStateActionStatus::IonicInsertion)
        .value("EnergyOrderUnavailable", SolvedStateActionStatus::EnergyOrderUnavailable)
        .value("DensityUnavailable", SolvedStateActionStatus::DensityUnavailable)
        .value("PressureSingular", SolvedStateActionStatus::PressureSingular)
        .value("AssociationSingular", SolvedStateActionStatus::AssociationSingular)
        .value("AssociationIllConditioned", SolvedStateActionStatus::AssociationIllConditioned)
        .value("AssociationResidual", SolvedStateActionStatus::AssociationResidual)
        .value("MissingThermochemistry", SolvedStateActionStatus::MissingThermochemistry)
        .value("ReferenceUnavailable", SolvedStateActionStatus::ReferenceUnavailable)
        .value("OutsideInterval", SolvedStateActionStatus::OutsideInterval)
        .value("ThermalBoundary", SolvedStateActionStatus::ThermalBoundary)
        .value("CentralPropertyUnavailable", SolvedStateActionStatus::CentralPropertyUnavailable)
        .value("Nonfinite", SolvedStateActionStatus::Nonfinite)
        .value("IllConditioned", SolvedStateActionStatus::IllConditioned);
    py::class_<SolvedStateActionResult>(module, "SolvedStateActionResult")
        .def_readonly("observable", &SolvedStateActionResult::observable)
        .def_readonly("unit", &SolvedStateActionResult::unit)
        .def_readonly("order", &SolvedStateActionResult::order)
        .def_readonly("value", &SolvedStateActionResult::value)
        .def_readonly("action", &SolvedStateActionResult::action)
        .def_readonly("status", &SolvedStateActionResult::status)
        .def_readonly("property_diagnostics", &SolvedStateActionResult::property_diagnostics)
        .def_readonly("diagnostic_message", &SolvedStateActionResult::diagnostic_message);
    py::class_<SolvedStateActionDiagnostics>(module, "SolvedStateActionDiagnostics")
        .def_readonly("equation_count", &SolvedStateActionDiagnostics::equation_count)
        .def_readonly("state_dimension", &SolvedStateActionDiagnostics::state_dimension)
        .def_readonly("maximum_energy_order", &SolvedStateActionDiagnostics::maximum_energy_order)
        .def_readonly("state_factorizations", &SolvedStateActionDiagnostics::state_factorizations)
        .def_readonly("linear_solves", &SolvedStateActionDiagnostics::linear_solves)
        .def_readonly("rank", &SolvedStateActionDiagnostics::rank)
        .def_readonly("condition_estimate", &SolvedStateActionDiagnostics::condition_estimate)
        .def_readonly("maximum_linear_residual", &SolvedStateActionDiagnostics::maximum_linear_residual)
        .def_readonly("central_residual", &SolvedStateActionDiagnostics::central_residual)
        .def_readonly("central_floor", &SolvedStateActionDiagnostics::central_floor);
    py::class_<SolvedStateActionBatch>(module, "SolvedStateActionBatch")
        .def_readonly("results", &SolvedStateActionBatch::results)
        .def_readonly("diagnostics", &SolvedStateActionBatch::diagnostics);
    module.def("solved_state_actions", &solved_state_actions,
               py::arg("compiled"), py::arg("central"), py::arg("request"));
    module.def("solved_state_first_actions", &solved_state_first_actions);
    module.def("compile_problem", &compile_problem, py::keep_alive<0, 1>());
    module.def("sensitivity", &sensitivity);
    module.def("solve_equilibrium", &solve_equilibrium, py::arg("model"), py::arg("problem"),
               py::arg("start") = std::nullopt);
}

}  // namespace epcsaft::equilibrium
