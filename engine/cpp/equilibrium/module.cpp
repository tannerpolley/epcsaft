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

#include "flash.hpp"
#include "phase_stability_certificates.hpp"
#include "phase_equilibrium.hpp"
#include "solved_state_sensitivity.hpp"

namespace py = pybind11;

namespace epcsaft_engine::equilibrium {
void bind_chemical_equilibrium(py::module_& module);
}

namespace {

py::dict source_reference_transfer(
    const epcsaft::NativeModel& model,
    double temperature_k,
    double pressure_pa,
    const std::vector<std::string>& component_ids,
    const std::vector<double>& source_solvent_mole_fractions,
    const std::vector<std::array<std::size_t, 2>>& ion_pair_indices,
    const std::vector<std::array<int, 2>>& ion_pair_stoichiometry,
    const std::vector<std::size_t>& pure_component_indices,
    const std::vector<std::string>& pure_component_phases,
    const std::string& phase,
    const std::string& reference_convention,
    const std::string& reference_state_id,
    const std::string& activity_convention_id,
    const std::string& reference_state_fingerprint,
    double standard_molality_mol_per_kg,
    double source_reference_pressure_pa,
    std::uint32_t required_derivatives
) {
    const epcsaft_engine::equilibrium::ProviderContext provider(model);
    const auto result = provider.evaluate_source_reference_transfer(
        temperature_k,
        pressure_pa,
        {
            component_ids,
            source_solvent_mole_fractions,
            ion_pair_indices,
            ion_pair_stoichiometry,
            pure_component_indices,
            pure_component_phases,
            phase,
            reference_convention,
            reference_state_id,
            activity_convention_id,
            reference_state_fingerprint,
            standard_molality_mol_per_kg,
            source_reference_pressure_pa,
            required_derivatives,
        }
    );
    py::dict value;
    value["status"] = result.status;
    value["error"] = result.error;
    if (result.status != epcsaft_engine::equilibrium::SOURCE_REFERENCE_OK) {
        return value;
    }
    value["component_count"] = result.component_count;
    value["neutral_basis_row_count"] = result.neutral_basis_row_count;
    value["neutral_basis"] = result.neutral_basis;
    value["log_fugacity_contractions"] = result.log_fugacity_contractions;
    value["activity_scale_log_contractions"] =
        result.activity_scale_log_contractions;
    value["transfer_log_contractions"] = result.transfer_log_contractions;
    value["source_composition"] = result.source_composition;
    value["derivative_availability"] = result.derivative_availability;
    value["pressure_derivatives_per_pa"] = result.pressure_derivatives_per_pa;
    value["parameter_derivatives"] = result.parameter_derivatives;
    value["temperature_k"] = result.temperature_k;
    value["pressure_pa"] = result.pressure_pa;
    value["source_reference_pressure_pa"] = result.source_reference_pressure_pa;
    value["standard_molality_mol_per_kg"] = result.standard_molality_mol_per_kg;
    value["solvent_molar_mass_kg_per_mol"] =
        result.solvent_molar_mass_kg_per_mol;
    value["reference_limit_molality_mol_per_kg"] =
        result.reference_limit_molality_mol_per_kg;
    value["reference_convergence_error"] = result.reference_convergence_error;
    value["source_temperature_min_k"] = result.source_temperature_min_k;
    value["source_temperature_max_k"] = result.source_temperature_max_k;
    value["source_pressure_min_pa"] = result.source_pressure_min_pa;
    value["source_pressure_max_pa"] = result.source_pressure_max_pa;
    value["candidate_temperature_min_k"] = result.candidate_temperature_min_k;
    value["candidate_temperature_max_k"] = result.candidate_temperature_max_k;
    value["candidate_pressure_min_pa"] = result.candidate_pressure_min_pa;
    value["candidate_pressure_max_pa"] = result.candidate_pressure_max_pa;
    value["candidate_qualification"] = result.candidate_qualification;
    value["candidate_source_fingerprint"] =
        result.candidate_source_fingerprint;
    value["parameter_fingerprint"] = result.parameter_fingerprint;
    value["topology_fingerprint"] = result.topology_fingerprint;
    value["component_order_fingerprint"] = result.component_order_fingerprint;
    value["reference_state_fingerprint"] = result.reference_state_fingerprint;
    value["domain_fingerprint"] = result.domain_fingerprint;
    value["artifact_fingerprint"] = result.artifact_fingerprint;
    value["helmholtz_basis_id"] = result.helmholtz_basis_id;
    return value;
}

py::dict evaluate_parameterized_phase(
    const epcsaft::NativeModel& model,
    double temperature_k,
    const std::vector<double>& amounts_mol,
    double volume_m3,
    double packing_fraction_min,
    double packing_fraction_max
) {
    const epcsaft_engine::equilibrium::ProviderContext provider(model);
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
    const epcsaft::NativeModel& model,
    double temperature_k,
    double pressure_pa,
    const std::vector<double>& feed_amounts_mol,
    const std::vector<std::vector<double>>& phase_amounts_mol,
    const std::vector<double>& phase_volumes_m3
) {
    const epcsaft_engine::equilibrium::ProviderContext provider(model);
    if (phase_amounts_mol.size() != phase_volumes_m3.size()) {
        throw py::value_error("phase amount and volume counts must match");
    }
    const std::vector<std::string> component_ids = provider.component_ids();
    std::vector<epcsaft_engine::equilibrium::SolvedPhaseState> phases;
    phases.reserve(phase_amounts_mol.size());
    for (std::size_t phase = 0; phase < phase_amounts_mol.size(); ++phase) {
        phases.push_back({phase_amounts_mol[phase], phase_volumes_m3[phase]});
    }
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
    if (value == "declared_neutral_components") {
        return epcsaft_engine::equilibrium::PhaseSupport::DeclaredNeutralComponents;
    }
    throw py::value_error(
        "phase support must be 'all_components', 'neutral_components', or 'declared_neutral_components'"
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

std::vector<epcsaft_engine::equilibrium::DeclaredPhaseInitialGuess>
parse_declared_phase_initial_guesses(
    const std::vector<std::size_t>& indices,
    const std::vector<std::vector<double>>& compositions,
    const std::vector<double>& densities,
    const std::vector<double>& phase_fractions
) {
    if (indices.size() != compositions.size()
        || indices.size() != densities.size()
        || indices.size() != phase_fractions.size()) {
        throw py::value_error(
            "declared phase initial-guess arrays must have equal length"
        );
    }
    std::vector<
        epcsaft_engine::equilibrium::DeclaredPhaseInitialGuess
    > result;
    result.reserve(indices.size());
    for (std::size_t guess = 0; guess < indices.size(); ++guess) {
        result.push_back({
            indices[guess],
            compositions[guess],
            densities[guess],
            phase_fractions[guess],
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
    result["initial_guess_count"] = value.initial_guess_count;
    result["accepted_start_count"] = value.accepted_start_count;
    result["distinct_root_count"] = value.distinct_root_count;
    result["final_correction_attempted"] =
        value.final_correction_attempted;
    result["final_correction_accepted"] =
        value.final_correction_accepted;
    result["final_correction_iteration_count"] =
        value.final_correction_iteration_count;
    result["pre_correction_residual_inf"] =
        value.pre_correction_residual_inf;
    result["post_correction_residual_inf"] =
        value.post_correction_residual_inf;
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

py::dict phase_stability_phase_to_dict(
    const epcsaft_engine::equilibrium::PhaseStabilityPhase& phase
) {
    py::dict result;
    result["phase_fraction"] = phase.phase_fraction;
    result["mole_fractions"] = phase.physical_fractions_provider_order;
    result["molar_volume_m3_mol"] = phase.volume;
    result["pressure_pa"] = phase.pressure_pa;
    result["chemical_potential_over_rt"] = phase.chemical_potentials_over_rt;
    return result;
}

py::dict phase_stability_discovery_to_dict(
    const epcsaft_engine::equilibrium::FlashResult& flash,
    int provider_evaluation_count
) {
    using epcsaft_engine::equilibrium::PhaseStabilityDiscoveryResult;
    const PhaseStabilityDiscoveryResult& solve = flash.solve;
    py::dict result;
    result["controller"] = "phase_stability_stages_1_to_7_grepe_v1";
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
    py::list stage_evidence;
    for (const auto& timing : solve.step_timings) {
        stage_evidence.append(
            "stage_" + std::to_string(timing.step) + ":"
            + timing.terminal_status + ":" + timing.terminal_reason
        );
    }
    result["stage_evidence"] = std::move(stage_evidence);
    result["optimizer_solves"] = optimizer_solves;
    result["optimizer_iterations"] = optimizer_iterations;
    int major_iterations = 0;
    for (const auto& step : solve.step4_history) {
        major_iterations = std::max(major_iterations, step.major_iteration + 1);
    }
    result["major_iterations"] = major_iterations;
    set_optional_double(
        result, "material_balance_max_abs", solve.material_balance_inf
    );
    set_optional_double(
        result, "electroneutrality_inf", solve.electroneutrality_inf
    );
    set_optional_double(
        result,
        "pressure_stationarity_max_relative",
        solve.pressure_residual_inf
    );
    set_optional_double(
        result,
        "chemical_potential_residual_inf",
        solve.chemical_potential_residual_inf
    );
    set_optional_double(
        result, "kkt_stationarity_max_abs", solve.kkt_stationarity_inf
    );
    set_optional_double(
        result, "phase_distinctness", solve.phase_distinctness
    );
    py::list phases;
    for (const auto& phase : solve.phases) {
        phases.append(phase_stability_phase_to_dict(phase));
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
    const epcsaft::NativeModel& model,
    double temperature_k,
    double pressure_pa,
    const std::vector<double>& feed_amounts_mol,
    const std::vector<std::string>& roles,
    const std::vector<std::string>& supports,
    const std::vector<std::size_t>& hint_indices,
    const std::vector<std::vector<double>>& hint_compositions,
    const std::vector<double>& hint_densities,
    const std::vector<std::size_t>& initial_guess_indices,
    const std::vector<std::vector<double>>& initial_guess_compositions,
    const std::vector<double>& initial_guess_densities,
    const std::vector<double>& initial_guess_phase_fractions
) {
    const epcsaft_engine::equilibrium::ProviderContext provider(model);
    const auto declared_phases = parse_declared_phases(roles, supports);
    const auto phase_hints = parse_phase_hints(
        hint_indices, hint_compositions, hint_densities
    );
    const auto initial_guesses = parse_declared_phase_initial_guesses(
        initial_guess_indices,
        initial_guess_compositions,
        initial_guess_densities,
        initial_guess_phase_fractions
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
                initial_guesses,
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
    const epcsaft::NativeModel& model,
    const std::string& role,
    double temperature_k,
    const std::vector<double>& mole_fractions,
    double molar_volume_m3_per_mol
) {
    const epcsaft_engine::equilibrium::ProviderContext provider(model);
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
    const epcsaft::NativeModel& reference_model,
    const epcsaft::NativeModel& partner_model,
    double temperature_k,
    const std::vector<std::string>& roles,
    const std::vector<std::string>& supports,
    std::size_t reference_phase_index,
    const std::vector<double>& reference_mole_fractions,
    const std::vector<std::string>& reference_direction_identities,
    const std::vector<std::vector<double>>& reference_composition_directions,
    const std::vector<std::size_t>& hint_indices,
    const std::vector<std::vector<double>>& hint_compositions,
    const std::vector<double>& hint_densities
) {
    const epcsaft_engine::equilibrium::ProviderContext reference_provider(
        reference_model
    );
    const epcsaft_engine::equilibrium::ProviderContext partner_provider(
        partner_model
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
                reference_direction_identities,
                reference_composition_directions,
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
    result["partner_parameter_fingerprint"] =
        solved.partner_parameter_fingerprint;
    result["partner_component_ids"] = solved.partner_component_ids;
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
    if (solved.reference_direction_sensitivity.has_value()) {
        result["reference_direction_sensitivity"] =
            solved_sensitivity_to_dict(
                *solved.reference_direction_sensitivity
            );
    } else {
        result["reference_direction_sensitivity"] = py::none();
    }
    if (solved.temperature_sensitivity.has_value()) {
        result["temperature_sensitivity"] = solved_sensitivity_to_dict(
            *solved.temperature_sensitivity
        );
    } else {
        result["temperature_sensitivity"] = py::none();
    }
    return result;
}
py::dict solve_phase_stability_discovery(
    const epcsaft::NativeModel& model,
    double temperature_k,
    double pressure_pa,
    const std::vector<double>& overall_mole_fractions
) {
    const epcsaft_engine::equilibrium::ProviderContext provider(model);
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
    return phase_stability_discovery_to_dict(flash, provider_evaluation_count);
}

py::dict phase_stability_audit_farkas_certificate(
    const std::vector<std::vector<double>>& matrix,
    const std::vector<double>& row_lower,
    const std::vector<double>& row_upper,
    const std::vector<double>& column_lower,
    const std::vector<double>& column_upper,
    const std::vector<double>& row_ray
) {
    const auto audit =
        epcsaft_engine::equilibrium::audit_phase_stability_farkas_certificate(
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

py::dict phase_stability_audit_step5_kkt(
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
    const auto audit = epcsaft_engine::equilibrium::audit_phase_stability_stage5_kkt(
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

void bind_equilibrium(py::module_& module) {
    epcsaft_engine::equilibrium::bind_chemical_equilibrium(module);
    module.def(
        "source_reference_transfer",
        &source_reference_transfer,
        py::arg("model"),
        py::arg("temperature_k"),
        py::arg("pressure_pa"),
        py::arg("component_ids"),
        py::arg("source_solvent_mole_fractions"),
        py::arg("ion_pair_indices"),
        py::arg("ion_pair_stoichiometry"),
        py::arg("pure_component_indices"),
        py::arg("pure_component_phases"),
        py::arg("phase"),
        py::arg("reference_convention"),
        py::arg("reference_state_id"),
        py::arg("activity_convention_id"),
        py::arg("reference_state_fingerprint"),
        py::arg("standard_molality_mol_per_kg"),
        py::arg("source_reference_pressure_pa"),
        py::arg("required_derivatives")
    );
    module.def(
        "evaluate_parameterized_phase",
        &evaluate_parameterized_phase,
        py::arg("model"),
        py::arg("temperature_k"),
        py::arg("amounts_mol"),
        py::arg("volume_m3"),
        py::arg("packing_fraction_min"),
        py::arg("packing_fraction_max")
    );
    module.def(
        "evaluate_solved_state_sensitivity",
        &evaluate_solved_state_sensitivity,
        py::arg("model"),
        py::arg("temperature_k"),
        py::arg("pressure_pa"),
        py::arg("feed_amounts_mol"),
        py::arg("phase_amounts_mol"),
        py::arg("phase_volumes_m3")
    );
    module.def(
        "_audit_phase_state",
        &audit_phase_state,
        py::arg("model"),
        py::arg("role"),
        py::arg("temperature_k"),
        py::arg("mole_fractions"),
        py::arg("molar_volume_m3_per_mol")
    );
    module.def(
        "_solve_declared_fixed_tp",
        &solve_declared_fixed_tp,
        py::arg("model"),
        py::arg("temperature_k"),
        py::arg("pressure_pa"),
        py::arg("feed_amounts_mol"),
        py::arg("roles"),
        py::arg("supports"),
        py::arg("hint_indices"),
        py::arg("hint_compositions"),
        py::arg("hint_densities_mol_m3"),
        py::arg("initial_guess_indices"),
        py::arg("initial_guess_compositions"),
        py::arg("initial_guess_densities_mol_m3"),
        py::arg("initial_guess_phase_fractions")
    );
    module.def(
        "_solve_fixed_t_coexistence",
        &solve_fixed_t_coexistence,
        py::arg("reference_model"),
        py::arg("partner_model"),
        py::arg("temperature_k"),
        py::arg("roles"),
        py::arg("supports"),
        py::arg("reference_phase_index"),
        py::arg("reference_mole_fractions"),
        py::arg("reference_direction_identities"),
        py::arg("reference_composition_directions"),
        py::arg("hint_indices"),
        py::arg("hint_compositions"),
        py::arg("hint_densities_mol_m3")
    );
    module.def(
        "_solve_phase_stability_discovery",
        &solve_phase_stability_discovery,
        py::arg("model"),
        py::arg("temperature_k"),
        py::arg("pressure_pa"),
        py::arg("overall_mole_fractions")
    );
    module.def(
        "_phase_stability_audit_farkas_certificate",
        &phase_stability_audit_farkas_certificate,
        py::arg("matrix"),
        py::arg("row_lower"),
        py::arg("row_upper"),
        py::arg("column_lower"),
        py::arg("column_upper"),
        py::arg("row_ray")
    );
    module.def(
        "_phase_stability_audit_step5_kkt",
        &phase_stability_audit_step5_kkt,
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
