#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "native_transport.hpp"

#include "flash.hpp"
#include "held2_certificates.hpp"
#include "phase_equilibrium.hpp"
#include "solved_state_sensitivity.hpp"

namespace py = pybind11;

namespace epcsaft_engine::equilibrium {
void bind_chemical_equilibrium(py::module_& module);
}

namespace {

const epcsaft_native_transport& checked_transport(const py::capsule& capsule) {
    const char* name = capsule.name();
    if (name == nullptr
        || std::string_view(name) != EPCSAFT_NATIVE_TRANSPORT_CAPSULE_NAME) {
        throw py::value_error("expected private Engine native transport");
    }
    const auto* transport = capsule.get_pointer<epcsaft_native_transport>();
    if (transport == nullptr || transport->model_context == nullptr) {
        throw py::value_error("native transport is incomplete");
    }
    return *transport;
}

const epcsaft_native_transport& checked_mixture_transport(
    const py::capsule& capsule
) {
    const epcsaft_native_transport& transport = checked_transport(capsule);
    if (transport.component_count < 2) {
        throw py::value_error("provider capsule must contain at least two components");
    }
    if (transport.evaluate_mixture_phase == nullptr
        && transport.evaluate_electrolyte_phase == nullptr) {
        throw py::value_error("provider capsule is missing its phase evaluator");
    }
    return transport;
}

py::dict evaluate_parameterized_phase(
    const py::capsule& capsule,
    double temperature_k,
    const std::vector<double>& amounts_mol,
    double volume_m3,
    double packing_fraction_min,
    double packing_fraction_max
) {
    const epcsaft_native_transport& transport = checked_transport(capsule);
    const epcsaft_engine::equilibrium::ProviderContext provider(transport);
    const auto active_parameters = provider.active_parameters();
    const auto evaluation = provider.evaluate_parameterized_phase(
        temperature_k,
        amounts_mol,
        volume_m3,
        packing_fraction_min,
        packing_fraction_max,
        active_parameters
    );
    py::dict result;
    result["active_parameter_identities"] = active_parameters.identities;
    result["value"] = evaluation.phase.value;
    result["state_gradient"] = evaluation.phase.gradient;
    result["state_hessian"] = evaluation.phase.hessian;
    result["parameter_derivatives"] = evaluation.parameter_derivatives;
    result["state_parameter_derivatives"] =
        evaluation.state_parameter_derivatives;
    result["pressure_pa"] = evaluation.phase.pressure_pa;
    result["pressure_parameter_derivatives_pa"] =
        evaluation.pressure_parameter_derivatives_pa;
    result["chemical_potential_parameter_derivatives_over_rt"] =
        evaluation.chemical_potential_parameter_derivatives_over_rt;
    result["packing_fraction"] = evaluation.packing.value;
    result["packing_state_gradient"] = evaluation.packing.gradient;
    result["packing_state_hessian"] = evaluation.packing.hessian;
    result["packing_parameter_derivatives"] =
        evaluation.packing_parameter_derivatives;
    result["molar_volume_bounds_m3_per_mol"] =
        evaluation.molar_volume_bounds_m3_per_mol;
    result["molar_volume_bound_parameter_derivatives_m3_per_mol"] =
        evaluation.molar_volume_bound_parameter_derivatives_m3_per_mol;
    result["topology_fingerprint"] = evaluation.topology_fingerprint;
    return result;
}

py::dict evaluate_solved_state_sensitivity(
    const py::capsule& capsule,
    double temperature_k,
    double pressure_pa,
    const std::vector<double>& feed_amounts_mol,
    const std::vector<std::vector<double>>& phase_amounts_mol,
    const std::vector<double>& phase_volumes_m3
) {
    const epcsaft_native_transport& transport = checked_transport(capsule);
    if (phase_amounts_mol.size() != phase_volumes_m3.size()) {
        throw py::value_error("phase amount and volume counts must match");
    }
    std::vector<std::string> component_ids;
    component_ids.reserve(transport.component_count);
    for (
        std::size_t component = 0;
        component < transport.component_count;
        ++component
    ) {
        const char* identity = transport.component_ids[component];
        if (identity == nullptr || identity[0] == '\0') {
            throw py::value_error("Provider component identity is incomplete");
        }
        component_ids.emplace_back(identity);
    }
    std::vector<epcsaft_engine::equilibrium::SolvedPhaseState> phases;
    phases.reserve(phase_amounts_mol.size());
    for (std::size_t phase = 0; phase < phase_amounts_mol.size(); ++phase) {
        phases.push_back({phase_amounts_mol[phase], phase_volumes_m3[phase]});
    }
    const epcsaft_engine::equilibrium::ProviderContext provider(transport);
    const auto sensitivity =
        epcsaft_engine::equilibrium::evaluate_solved_state_sensitivity(
            provider,
            temperature_k,
            pressure_pa,
            component_ids,
            feed_amounts_mol,
            phases
        );
    py::dict result;
    result["status"] = sensitivity.status;
    result["failure_reason"] = sensitivity.failure_reason;
    result["derivative_order"] = 1;
    result["higher_order_status"] = "unavailable";
    result["higher_order_reason"] =
        "higher_order_total_derivatives_not_implemented";
    result["tensor_layout"] = "row-major";
    result["residual_identities"] = sensitivity.residual_identities;
    result["state_identities"] = sensitivity.state_identities;
    result["output_identities"] = sensitivity.output_identities;
    result["parameter_identities"] = sensitivity.parameter_identities;
    result["residual_state_jacobian"] =
        sensitivity.residual_state_jacobian;
    result["residual_parameter_jacobian"] =
        sensitivity.residual_parameter_jacobian;
    result["state_parameter_jacobian"] =
        sensitivity.state_parameter_jacobian;
    result["output_parameter_jacobian"] =
        sensitivity.output_parameter_jacobian;
    result["rank"] = sensitivity.rank;
    result["dimension"] = sensitivity.dimension;
    result["parameter_count"] = sensitivity.parameter_identities.size();
    result["output_count"] = sensitivity.output_identities.size();
    result["scaled_condition_number_2"] =
        sensitivity.scaled_condition_number_2;
    result["linear_solve_backward_error"] =
        sensitivity.linear_solve_backward_error;
    result["reference_phase_index"] = sensitivity.reference_phase_index;
    return result;
}


epcsaft_engine::equilibrium::PhaseRole parse_phase_role(
    const std::string& value
) {
    if (value == "liquid") {
        return epcsaft_engine::equilibrium::PhaseRole::Liquid;
    }
    if (value == "vapor") {
        return epcsaft_engine::equilibrium::PhaseRole::Vapor;
    }
    throw py::value_error("phase role must be 'liquid' or 'vapor'");
}

epcsaft_engine::equilibrium::PhaseSupport parse_phase_support(
    const std::string& value
) {
    if (value == "all_components") {
        return epcsaft_engine::equilibrium::PhaseSupport::AllComponents;
    }
    if (value == "neutral_components") {
        return epcsaft_engine::equilibrium::PhaseSupport::NeutralComponents;
    }
    throw py::value_error(
        "phase support must be 'all_components' or 'neutral_components'"
    );
}

std::vector<epcsaft_engine::equilibrium::DeclaredPhaseSpecification>
parse_declared_phases(
    const std::vector<std::string>& roles,
    const std::vector<std::string>& supports
) {
    if (roles.size() != supports.size()) {
        throw py::value_error("phase role and support counts must match");
    }
    std::vector<
        epcsaft_engine::equilibrium::DeclaredPhaseSpecification
    > result;
    result.reserve(roles.size());
    for (std::size_t phase = 0; phase < roles.size(); ++phase) {
        result.push_back({
            parse_phase_role(roles[phase]),
            parse_phase_support(supports[phase]),
        });
    }
    return result;
}

std::vector<epcsaft_engine::equilibrium::PhaseHint> parse_phase_hints(
    const std::vector<std::size_t>& indices,
    const std::vector<std::vector<double>>& compositions,
    const std::vector<double>& densities
) {
    if (indices.size() != compositions.size()
        || indices.size() != densities.size()) {
        throw py::value_error("phase hint arrays must have equal length");
    }
    std::vector<epcsaft_engine::equilibrium::PhaseHint> result;
    result.reserve(indices.size());
    for (std::size_t hint = 0; hint < indices.size(); ++hint) {
        result.push_back({
            indices[hint],
            compositions[hint],
            densities[hint],
        });
    }
    return result;
}

py::dict phase_diagnostics_to_dict(
    const epcsaft_engine::equilibrium::PhaseEquilibriumDiagnostics& value
) {
    py::dict result;
    result["route"] = value.route;
    result["solver_status"] = value.solver_status;
    result["numerical_status"] = value.numerical_status;
    result["physical_status"] = value.physical_status;
    result["failure_stage"] = value.failure_stage.empty()
        ? py::none()
        : py::cast(value.failure_stage);
    result["failure_reason"] = value.failure_reason.empty()
        ? py::none()
        : py::cast(value.failure_reason);
    result["certificate_scope"] = value.certificate_scope;
    result["phase_selection_claim"] = value.phase_selection_claim;
    result["globality_claim"] = value.globality_claim;
    result["search_scope"] = value.search_scope;
    result["search_budget"] = value.search_budget;
    result["search_exclusions"] = value.search_exclusions;
    result["completeness_criterion"] = value.completeness_criterion;
    result["role_evidence"] = value.role_evidence;
    result["optimizer_solve_count"] = value.optimizer_solve_count;
    result["optimizer_iteration_count"] = value.optimizer_iteration_count;
    result["provider_evaluation_count"] = value.provider_evaluation_count;
    result["hint_count"] = value.hint_count;
    result["hint_reordered_start_count"] =
        value.hint_reordered_start_count;
    result["material_balance_inf"] = value.material_balance_inf;
    result["electroneutrality_inf"] = value.electroneutrality_inf;
    result["pressure_residual_inf"] = value.pressure_residual_inf;
    result["chemical_potential_residual_inf"] =
        value.chemical_potential_residual_inf;
    result["kkt_stationarity_inf"] = value.kkt_stationarity_inf;
    result["dual_sign_violation_inf"] = value.dual_sign_violation_inf;
    result["complementarity_inf"] = value.complementarity_inf;
    result["phase_distinctness"] = value.phase_distinctness;
    return result;
}

py::dict phase_state_to_dict(
    const epcsaft_engine::equilibrium::PhaseEquilibriumState& value
) {
    py::dict result;
    result["role"] =
        epcsaft_engine::equilibrium::phase_role_name(value.role);
    result["support"] =
        epcsaft_engine::equilibrium::phase_support_name(value.support);
    result["supported_component_ids"] = value.supported_component_ids;
    result["mole_fractions"] = value.mole_fractions;
    result["molar_volume_m3_per_mol"] = value.molar_volume_m3_per_mol;
    result["molar_density_mol_m3"] = value.molar_density_mol_m3;
    result["packing_fraction"] = value.packing_fraction;
    result["pressure_pa"] = value.pressure_pa;
    result["chemical_potential_component_ids"] =
        value.chemical_potential_component_ids;
    result["chemical_potential_over_rt"] =
        value.chemical_potential_over_rt;
    result["mechanical_class"] = value.mechanical_class;
    return result;
}

py::dict solved_sensitivity_to_dict(
    const epcsaft_engine::equilibrium::SolvedStateSensitivity& value
) {
    py::dict result;
    result["status"] = value.status;
    result["failure_reason"] = value.failure_reason.empty()
        ? py::none()
        : py::cast(value.failure_reason);
    result["derivative_order"] = 1;
    result["higher_order_status"] = "unavailable";
    result["higher_order_reason"] =
        "higher_order_total_derivatives_not_implemented";
    result["tensor_layout"] = "row-major";
    result["residual_identities"] = value.residual_identities;
    result["state_identities"] = value.state_identities;
    result["output_identities"] = value.output_identities;
    result["parameter_identities"] = value.parameter_identities;
    result["residual_state_jacobian"] =
        value.residual_state_jacobian;
    result["residual_parameter_jacobian"] =
        value.residual_parameter_jacobian;
    result["state_parameter_jacobian"] =
        value.state_parameter_jacobian;
    result["output_parameter_jacobian"] =
        value.output_parameter_jacobian;
    result["rank"] = value.rank;
    result["dimension"] = value.dimension;
    result["parameter_count"] = value.parameter_identities.size();
    result["output_count"] = value.output_identities.size();
    if (std::isfinite(value.scaled_condition_number_2)) {
        result["scaled_condition_number_2"] = value.scaled_condition_number_2;
    } else {
        result["scaled_condition_number_2"] = py::none();
    }
    if (std::isfinite(value.linear_solve_backward_error)) {
        result["linear_solve_backward_error"] = value.linear_solve_backward_error;
    } else {
        result["linear_solve_backward_error"] = py::none();
    }
    result["reference_phase_index"] = value.reference_phase_index;
    return result;
}

void set_optional_double(
    py::dict& result,
    const char* name,
    const std::optional<double>& value
) {
    result[name] = value ? py::cast(*value) : py::none();
}

py::dict neutral_phase_to_dict(
    const epcsaft_engine::equilibrium::HeldStageIIIPhaseEvaluation& phase
) {
    const double amount_mol = phase.amounts_mol[0] + phase.amounts_mol[1];
    py::dict result;
    result["amount_mol"] = amount_mol;
    result["mole_fractions"] = std::array<double, 2>{
        phase.amounts_mol[0] / amount_mol,
        phase.amounts_mol[1] / amount_mol,
    };
    result["volume_m3"] = phase.volume_m3;
    result["molar_density_mol_m3"] = amount_mol / phase.volume_m3;
    result["pressure_pa"] = phase.provider.pressure_pa;
    result["chemical_potential_over_rt"] = std::array<double, 2>{
        phase.provider.gradient[0],
        phase.provider.gradient[1],
    };
    return result;
}

py::dict neutral_discovery_to_dict(
    const epcsaft_engine::equilibrium::FlashResult& flash,
    int provider_evaluation_count
) {
    using epcsaft_engine::equilibrium::HeldResult;
    const HeldResult& solve = std::get<HeldResult>(flash.solve);
    int attempts = static_cast<int>(
        solve.stage_i.reference_attempts.size()
        + solve.stage_i.attempt_log.size()
        + solve.stage_ii.endpoint_attempts.size()
    );
    for (const auto& trace : solve.stage_ii.trace) {
        attempts += trace.lower_starts_completed;
    }
    for (const auto& refinement : solve.stage_iii_attempts) {
        attempts += refinement.result.local.attempts;
    }

    py::list profiles;
    profiles.append(solve.stage_i.search_profile);
    if (solve.stage_i.outcome == "negative_tpd") {
        profiles.append(solve.stage_ii.search_profile);
    }
    if (!solve.stage_iii_attempts.empty()) {
        profiles.append(solve.stage_iii_attempts.front().result.search_profile);
    }

    py::dict result;
    result["controller"] = "held_binary_v1";
    result["outcome"] = solve.outcome;
    result["search_status"] = solve.search_status;
    result["failure_reason"] = solve.failure_reason;
    result["attempts"] = attempts;
    result["major_iterations"] = solve.stage_ii.major_iterations;
    result["material_balance_max_abs"] = py::none();
    result["pressure_stationarity_max_relative"] = py::none();
    result["kkt_stationarity_max_abs"] = py::none();
    result["chemical_potential_max_relative"] = py::none();
    result["search_profiles"] = std::move(profiles);
    result["globality_certificate"] = "not_guaranteed";
    result["provider_evaluation_count"] = provider_evaluation_count;
    result["solver_status"] = "not_adjudicated";
    result["numerical_status"] = "not_adjudicated";
    result["physical_status"] = "not_adjudicated";

    const bool accepted_reference = solve.stage_i.has_reference
        && std::any_of(
            solve.stage_i.reference_attempts.begin(),
            solve.stage_i.reference_attempts.end(),
            [](const auto& attempt) { return attempt.accepted; }
        );
    const bool reference_solver_failed = std::any_of(
        solve.stage_i.reference_attempts.begin(),
        solve.stage_i.reference_attempts.end(),
        [](const auto& attempt) {
            return !attempt.solver_converged || !attempt.callback_error.empty();
        }
    );
    const bool reference_failed = std::any_of(
        solve.stage_i.reference_attempts.begin(),
        solve.stage_i.reference_attempts.end(),
        [](const auto& attempt) { return !attempt.accepted; }
    );
    const bool stage_i_solver_failed = std::any_of(
        solve.stage_i.attempt_log.begin(),
        solve.stage_i.attempt_log.end(),
        [](const auto& attempt) {
            return !attempt.solver_converged || !attempt.callback_error.empty();
        }
    );
    const bool stage_i_search_failed = std::any_of(
        solve.stage_i.attempt_log.begin(),
        solve.stage_i.attempt_log.end(),
        [](const auto& attempt) {
            return !attempt.solver_converged || !attempt.callback_error.empty()
                || !attempt.accepted;
        }
    );
    const bool completed_one_phase = accepted_reference
        && solve.stage_i.planned_starts.size() == 20
        && solve.stage_i.starts_completed == 20
        && solve.stage_i.negative_confirmations == 0
        && std::none_of(
            solve.stage_i.attempt_log.begin(),
            solve.stage_i.attempt_log.end(),
            [](const auto& attempt) { return attempt.kind == "confirmation"; }
        )
        && std::isfinite(solve.stage_i.best_tpd)
        && solve.stage_i.best_tpd >= -1.0e-8;
    if (completed_one_phase) {
        const bool solver_passed = !reference_solver_failed
            && !stage_i_solver_failed;
        const bool numerical_passed = solver_passed
            && !reference_failed && !stage_i_search_failed;
        result["solver_status"] = solver_passed ? "passed" : "failed";
        result["numerical_status"] = numerical_passed ? "passed" : "failed";
        result["physical_status"] =
            numerical_passed ? "passed" : "not_adjudicated";
    } else if (solve.outcome == "accepted" && !solve.stage_iii_attempts.empty()) {
        const auto& local = solve.stage_iii_attempts.back().result.local;
        result["solver_status"] = local.solver_converged ? "passed" : "failed";
        result["numerical_status"] = local.numerical_converged ? "passed" : "failed";
        result["physical_status"] = local.physical_accepted ? "passed" : "failed";
    }

    py::list phases;
    py::list phase_fractions;
    if (solve.outcome == "one_phase" && solve.stage_i.has_reference) {
        const auto& state = solve.stage_i.reference;
        py::dict phase;
        phase["amount_mol"] = 1.0;
        phase["mole_fractions"] = state.amounts_mol;
        phase["volume_m3"] = state.volume_m3;
        phase["molar_density_mol_m3"] = 1.0 / state.volume_m3;
        phase["pressure_pa"] = state.provider.pressure_pa;
        phase["chemical_potential_over_rt"] = std::array<double, 2>{
            state.provider.gradient[0], state.provider.gradient[1]
        };
        phases.append(std::move(phase));
        phase_fractions.append(1.0);
        result["total_free_energy_over_rt"] = state.g_bar;
        result["pressure_stationarity_max_relative"] =
            std::abs(state.pressure_stationarity_relative);
    } else if (!solve.stage_iii_attempts.empty()) {
        const auto& refinement = solve.stage_iii_attempts.back().result;
        const auto& local = refinement.local;
        if (local.solver_converged) {
            result["material_balance_max_abs"] = local.material_balance_max_abs;
            result["pressure_stationarity_max_relative"] =
                local.pressure_stationarity_max_relative;
            result["kkt_stationarity_max_abs"] = local.kkt_stationarity_max_abs;
            result["chemical_potential_max_relative"] =
                refinement.chemical_potential_max_relative;
        }
        if (solve.outcome == "accepted") {
            phases.append(neutral_phase_to_dict(local.evaluation.liquid));
            phases.append(neutral_phase_to_dict(local.evaluation.vapor));
            phase_fractions.append(
                local.evaluation.liquid.amounts_mol[0]
                + local.evaluation.liquid.amounts_mol[1]
            );
            phase_fractions.append(
                local.evaluation.vapor.amounts_mol[0]
                + local.evaluation.vapor.amounts_mol[1]
            );
            result["total_free_energy_over_rt"] = local.evaluation.objective;
        }
    }
    result["phases"] = std::move(phases);
    result["phase_fractions"] = std::move(phase_fractions);
    result["parameter_fingerprint"] = flash.parameter_fingerprint;
    if (flash.sensitivity) {
        result["sensitivity"] = solved_sensitivity_to_dict(*flash.sensitivity);
    } else {
        result["sensitivity"] = py::none();
    }
    return result;
}

py::dict held2_phase_to_dict(
    const epcsaft_engine::equilibrium::Held2Phase& phase
) {
    py::dict result;
    result["phase_fraction"] = phase.phase_fraction;
    result["mole_fractions"] = phase.physical_fractions_provider_order;
    result["molar_volume_m3_mol"] = phase.volume;
    result["pressure_pa"] = phase.pressure_pa;
    result["chemical_potential_over_rt"] = phase.chemical_potentials_over_rt;
    return result;
}

py::dict held2_discovery_to_dict(
    const epcsaft_engine::equilibrium::FlashResult& flash,
    int provider_evaluation_count
) {
    using epcsaft_engine::equilibrium::Held2AlgorithmResult;
    const Held2AlgorithmResult& solve =
        std::get<Held2AlgorithmResult>(flash.solve);
    py::dict result;
    result["controller"] = "held2_steps_1_to_10_v1";
    result["parameter_fingerprint"] = flash.parameter_fingerprint;
    result["outcome"] = solve.outcome;
    result["failure_stage"] = solve.failure_stage.empty()
        ? py::none() : py::cast(solve.failure_stage);
    result["failure_reason"] = solve.failure_reason.empty()
        ? py::none() : py::cast(solve.failure_reason);
    result["globality_certificate"] = solve.globality_certificate;
    result["provider_evaluation_count"] = provider_evaluation_count;
    std::uint64_t optimizer_solves = 0;
    std::uint64_t optimizer_iterations = 0;
    for (const auto& timing : solve.step_timings) {
        optimizer_solves += timing.optimizer_solves;
        optimizer_iterations += timing.optimizer_iterations;
    }
    result["optimizer_solves"] = optimizer_solves;
    result["optimizer_iterations"] = optimizer_iterations;
    int major_iterations = 0;
    for (const auto& step : solve.step4_history) {
        major_iterations = std::max(major_iterations, step.major_iteration + 1);
    }
    result["major_iterations"] = major_iterations;
    result["material_balance_max_abs"] = py::none();
    result["electroneutrality_inf"] = py::none();
    result["pressure_stationarity_max_relative"] = py::none();
    result["kkt_stationarity_max_abs"] = py::none();
    if (solve.step10 && solve.step10->final_certificate) {
        const auto& certificate = *solve.step10->final_certificate;
        result["material_balance_max_abs"] = std::max(
            certificate.modified_balance_inf,
            certificate.ordinary_balance_inf
        );
        result["electroneutrality_inf"] = certificate.electroneutrality_inf;
        result["pressure_stationarity_max_relative"] =
            certificate.pressure_residual_inf;
        result["kkt_stationarity_max_abs"] = certificate.kkt_residual_inf;
    }
    py::list phases;
    for (const auto& phase : solve.phases) {
        phases.append(held2_phase_to_dict(phase));
    }
    result["phases"] = std::move(phases);
    set_optional_double(
        result,
        "total_free_energy_over_rt",
        solve.total_free_energy_over_rt
    );
    if (flash.sensitivity) {
        result["sensitivity"] = solved_sensitivity_to_dict(*flash.sensitivity);
    } else {
        result["sensitivity"] = py::none();
    }
    return result;
}

py::dict solve_declared_fixed_tp(
    const py::capsule& capsule,
    double temperature_k,
    double pressure_pa,
    const std::vector<double>& feed_amounts_mol,
    const std::vector<std::string>& roles,
    const std::vector<std::string>& supports,
    const std::vector<std::size_t>& hint_indices,
    const std::vector<std::vector<double>>& hint_compositions,
    const std::vector<double>& hint_densities
) {
    const epcsaft_engine::equilibrium::ProviderContext provider(
        checked_mixture_transport(capsule)
    );
    const auto declared_phases = parse_declared_phases(roles, supports);
    const auto phase_hints = parse_phase_hints(
        hint_indices, hint_compositions, hint_densities
    );
    epcsaft_engine::equilibrium::FixedTpPhaseResult solved;
    {
        py::gil_scoped_release release;
        solved = epcsaft_engine::equilibrium::solve_declared_fixed_tp(
            provider,
            {
                temperature_k,
                pressure_pa,
                feed_amounts_mol,
                declared_phases,
                phase_hints,
            }
        );
    }
    py::dict result;
    py::list phases;
    for (const auto& phase : solved.phases) {
        phases.append(phase_state_to_dict(phase));
    }
    result["phases"] = std::move(phases);
    result["phase_amounts_mol"] = solved.phase_amounts_mol;
    result["phase_fractions"] = solved.phase_fractions;
    result["total_free_energy_over_rt"] =
        solved.total_free_energy_over_rt;
    result["parameter_fingerprint"] = solved.parameter_fingerprint;
    result["diagnostics"] = phase_diagnostics_to_dict(
        solved.diagnostics
    );
    if (solved.sensitivity.has_value()) {
        result["sensitivity"] = solved_sensitivity_to_dict(
            *solved.sensitivity
        );
    } else {
        result["sensitivity"] = py::none();
    }
    return result;
}

py::dict audit_phase_state(
    const py::capsule& capsule,
    const std::string& role,
    double temperature_k,
    const std::vector<double>& mole_fractions,
    double molar_volume_m3_per_mol
) {
    const epcsaft_engine::equilibrium::ProviderContext provider(
        checked_mixture_transport(capsule)
    );
    const auto parsed_role = parse_phase_role(role);
    epcsaft_engine::equilibrium::PhaseEquilibriumState audited;
    {
        py::gil_scoped_release release;
        audited = epcsaft_engine::equilibrium::audit_phase_state(
            provider,
            parsed_role,
            temperature_k,
            mole_fractions,
            molar_volume_m3_per_mol
        );
    }
    return phase_state_to_dict(audited);
}

py::dict solve_fixed_t_coexistence(
    const py::capsule& reference_capsule,
    const py::capsule& partner_capsule,
    double temperature_k,
    const std::vector<std::string>& roles,
    const std::vector<std::string>& supports,
    std::size_t reference_phase_index,
    const std::vector<double>& reference_mole_fractions,
    const std::vector<std::size_t>& hint_indices,
    const std::vector<std::vector<double>>& hint_compositions,
    const std::vector<double>& hint_densities
) {
    const epcsaft_engine::equilibrium::ProviderContext reference_provider(
        checked_mixture_transport(reference_capsule)
    );
    const epcsaft_engine::equilibrium::ProviderContext partner_provider(
        checked_transport(partner_capsule)
    );
    const auto declared_phases = parse_declared_phases(roles, supports);
    const auto phase_hints = parse_phase_hints(
        hint_indices, hint_compositions, hint_densities
    );
    epcsaft_engine::equilibrium::FixedTCoexistenceResult solved;
    {
        py::gil_scoped_release release;
        solved = epcsaft_engine::equilibrium::solve_fixed_t_coexistence(
            reference_provider,
            partner_provider,
            {
                temperature_k,
                declared_phases,
                reference_phase_index,
                reference_mole_fractions,
                phase_hints,
            }
        );
    }
    py::dict result;
    py::list phases;
    for (const auto& phase : solved.phases) {
        phases.append(phase_state_to_dict(phase));
    }
    result["phases"] = std::move(phases);
    result["reference_phase_index"] = solved.reference_phase_index;
    result["pressure_pa"] = solved.pressure_pa;
    result["parameter_fingerprint"] = solved.parameter_fingerprint;
    result["diagnostics"] = phase_diagnostics_to_dict(
        solved.diagnostics
    );
    if (solved.sensitivity.has_value()) {
        result["sensitivity"] = solved_sensitivity_to_dict(
            *solved.sensitivity
        );
    } else {
        result["sensitivity"] = py::none();
    }
    return result;
}
py::dict solve_held_discovery(
    const py::capsule& capsule,
    double temperature_k,
    double pressure_pa,
    const std::vector<double>& overall_mole_fractions
) {
    const epcsaft_engine::equilibrium::ProviderContext provider(
        checked_transport(capsule)
    );
    epcsaft_engine::equilibrium::FlashResult flash;
    {
        py::gil_scoped_release release;
        flash = epcsaft_engine::equilibrium::solve_tp_flash(
            provider,
            {
                temperature_k,
                pressure_pa,
                overall_mole_fractions,
            }
        );
    }
    const int provider_evaluation_count = provider.evaluation_count();
    if (std::holds_alternative<
            epcsaft_engine::equilibrium::Held2AlgorithmResult
        >(flash.solve)) {
        return held2_discovery_to_dict(flash, provider_evaluation_count);
    }
    return neutral_discovery_to_dict(flash, provider_evaluation_count);
}

py::dict held2_audit_farkas_certificate(
    const std::vector<std::vector<double>>& matrix,
    const std::vector<double>& row_lower,
    const std::vector<double>& row_upper,
    const std::vector<double>& column_lower,
    const std::vector<double>& column_upper,
    const std::vector<double>& row_ray
) {
    const auto audit =
        epcsaft_engine::equilibrium::audit_held2_farkas_certificate(
            matrix,
            row_lower,
            row_upper,
            column_lower,
            column_upper,
            row_ray
        );
    py::dict result;
    result["reason"] = audit.reason;
    result["row_multipliers"] = audit.row_multipliers;
    result["normalization"] = audit.normalization;
    result["row_sign_violation_inf"] = audit.row_sign_violation_inf;
    result["dual_feasibility_violation_inf"] =
        audit.dual_feasibility_violation_inf;
    result["contradiction_margin"] = audit.contradiction_margin;
    result["contradiction_scale"] = audit.contradiction_scale;
    result["contradiction_threshold"] = audit.contradiction_threshold;
    result["solver_ray_recovered_without_presolve"] =
        audit.solver_ray_recovered_without_presolve;
    result["accepted"] = audit.accepted;
    return result;
}

py::dict held2_audit_step5_kkt(
    const std::vector<double>& variables,
    const std::vector<double>& lower,
    const std::vector<double>& upper,
    const std::vector<double>& objective_gradient,
    const std::vector<std::vector<double>>& constraints,
    const std::vector<double>& constraint_upper,
    const std::vector<double>& lower_bound_multipliers,
    const std::vector<double>& upper_bound_multipliers,
    const std::vector<double>& constraint_multipliers,
    double pullback_residual,
    double pullback_scale,
    double pressure_residual,
    bool same_major_iteration,
    bool step4_binding_valid,
    bool pressure_branch_valid
) {
    const auto audit = epcsaft_engine::equilibrium::audit_held2_step5_kkt(
        variables,
        lower,
        upper,
        objective_gradient,
        constraints,
        constraint_upper,
        lower_bound_multipliers,
        upper_bound_multipliers,
        constraint_multipliers,
        pullback_residual,
        pullback_scale,
        pressure_residual,
        same_major_iteration,
        step4_binding_valid,
        pressure_branch_valid
    );
    py::dict result;
    result["reason"] = audit.reason;
    result["primal_residual_inf"] = audit.primal_residual_inf;
    result["dual_sign_violation_inf"] = audit.dual_sign_violation_inf;
    result["pullback_residual_inf"] = audit.pullback_residual_inf;
    result["pullback_scale"] = audit.pullback_scale;
    result["pressure_residual"] = audit.pressure_residual;
    result["stationarity_residual_inf"] =
        audit.stationarity_residual_inf;
    result["complementarity_inf"] = audit.complementarity_inf;
    result["active_constraint_count"] = audit.active_constraint_count;
    result["active_jacobian_rank"] = audit.active_jacobian_rank;
    result["accepted"] = audit.accepted;
    return result;
}

}  // namespace

PYBIND11_MODULE(_equilibrium, module) {
    module.doc() = "Native local equilibrium formulation";
    epcsaft_engine::equilibrium::bind_chemical_equilibrium(module);
    module.def(
        "evaluate_parameterized_phase",
        &evaluate_parameterized_phase,
        py::arg("capsule"),
        py::arg("temperature_k"),
        py::arg("amounts_mol"),
        py::arg("volume_m3"),
        py::arg("packing_fraction_min"),
        py::arg("packing_fraction_max")
    );
    module.def(
        "evaluate_solved_state_sensitivity",
        &evaluate_solved_state_sensitivity,
        py::arg("capsule"),
        py::arg("temperature_k"),
        py::arg("pressure_pa"),
        py::arg("feed_amounts_mol"),
        py::arg("phase_amounts_mol"),
        py::arg("phase_volumes_m3")
    );
    module.def(
        "_audit_phase_state",
        &audit_phase_state,
        py::arg("capsule"),
        py::arg("role"),
        py::arg("temperature_k"),
        py::arg("mole_fractions"),
        py::arg("molar_volume_m3_per_mol")
    );
    module.def(
        "_solve_declared_fixed_tp",
        &solve_declared_fixed_tp,
        py::arg("capsule"),
        py::arg("temperature_k"),
        py::arg("pressure_pa"),
        py::arg("feed_amounts_mol"),
        py::arg("roles"),
        py::arg("supports"),
        py::arg("hint_indices"),
        py::arg("hint_compositions"),
        py::arg("hint_densities_mol_m3")
    );
    module.def(
        "_solve_fixed_t_coexistence",
        &solve_fixed_t_coexistence,
        py::arg("reference_capsule"),
        py::arg("partner_capsule"),
        py::arg("temperature_k"),
        py::arg("roles"),
        py::arg("supports"),
        py::arg("reference_phase_index"),
        py::arg("reference_mole_fractions"),
        py::arg("hint_indices"),
        py::arg("hint_compositions"),
        py::arg("hint_densities_mol_m3")
    );
    module.def(
        "_solve_held_discovery",
        &solve_held_discovery,
        py::arg("capsule"),
        py::arg("temperature_k"),
        py::arg("pressure_pa"),
        py::arg("overall_mole_fractions")
    );
    module.def(
        "_held2_audit_farkas_certificate",
        &held2_audit_farkas_certificate,
        py::arg("matrix"),
        py::arg("row_lower"),
        py::arg("row_upper"),
        py::arg("column_lower"),
        py::arg("column_upper"),
        py::arg("row_ray")
    );
    module.def(
        "_held2_audit_step5_kkt",
        &held2_audit_step5_kkt,
        py::arg("variables"),
        py::arg("lower"),
        py::arg("upper"),
        py::arg("objective_gradient"),
        py::arg("constraints"),
        py::arg("constraint_upper"),
        py::arg("lower_bound_multipliers"),
        py::arg("upper_bound_multipliers"),
        py::arg("constraint_multipliers"),
        py::arg("pullback_residual"),
        py::arg("pullback_scale"),
        py::arg("pressure_residual"),
        py::arg("same_major_iteration"),
        py::arg("step4_binding_valid"),
        py::arg("pressure_branch_valid")
    );
}
