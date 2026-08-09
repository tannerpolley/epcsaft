#include "result_json.hpp"

#include <cmath>
#include <optional>
#include <string>

#include <nlohmann/json.hpp>

namespace epcsaft_engine::equilibrium {
namespace {

using Json = nlohmann::json;

Json finite_or_null(double value) {
    return std::isfinite(value) ? Json(value) : Json(nullptr);
}

Json phase_to_json(const PhaseEquilibriumState& phase) {
    return {
        {"role", phase_role_name(phase.role)},
        {"amount_mol", phase.amount_mol},
        {"amounts_mol", phase.amounts_mol},
        {"mole_fractions", phase.mole_fractions},
        {"volume_m3", phase.volume_m3},
        {"molar_density_mol_m3", phase.molar_density_mol_m3},
        {"packing_fraction", phase.packing_fraction},
        {"pressure_pa", phase.pressure_pa},
        {"chemical_potential_over_rt", phase.chemical_potential_over_rt},
    };
}

Json phases_to_json(const std::vector<PhaseEquilibriumState>& phases) {
    Json result = Json::array();
    for (const auto& phase : phases) {
        result.push_back(phase_to_json(phase));
    }
    return result;
}

Json diagnostics_to_json(const PhaseEquilibriumDiagnostics& value) {
    return {
        {"solver_status", value.solver_status},
        {"numerical_status", value.numerical_status},
        {"physical_status", value.physical_status},
        {"failure_stage", value.failure_stage.empty() ? Json(nullptr) : Json(value.failure_stage)},
        {"failure_reason", value.failure_reason.empty() ? Json(nullptr) : Json(value.failure_reason)},
        {"failure_detail", value.failure_detail.empty() ? Json(nullptr) : Json(value.failure_detail)},
        {"certificate_scope", value.certificate_scope},
        {"phase_selection_claim", value.phase_selection_claim},
        {"globality_claim", value.globality_claim},
        {"completeness_criterion", value.completeness_criterion},
        {"searched_topologies", value.searched_topologies},
        {"excluded_topologies", value.excluded_topologies},
        {"equation_state_identities", value.equation_state_identities},
        {"equation_residual_identities", value.equation_residual_identities},
        {"tensor_layout", value.tensor_layout},
        {"solver_derivative_contract", value.solver_derivative_contract},
        {"sensitivity_derivative_contract", value.sensitivity_derivative_contract},
        {"parameter_parameter_requirement", value.parameter_parameter_requirement},
        {"topology_budget", value.topology_budget},
        {"topology_count", value.topology_count},
        {"start_budget_per_topology", value.start_budget_per_topology},
        {"attempted_start_count", value.attempted_start_count},
        {"certified_candidate_count", value.certified_candidate_count},
        {"optimizer_iteration_count", value.optimizer_iteration_count},
        {"hint_count", value.hint_count},
        {"hint_reordered_start_count", value.hint_reordered_start_count},
        {"material_balance_inf", value.material_balance_inf},
        {"electroneutrality_inf", value.electroneutrality_inf},
        {"pressure_residual_inf", value.pressure_residual_inf},
        {"chemical_potential_residual_inf", value.chemical_potential_residual_inf},
        {"kkt_stationarity_inf", value.kkt_stationarity_inf},
        {"minimum_reduced_hessian_eigenvalue", finite_or_null(value.minimum_reduced_hessian_eigenvalue)},
        {"kkt_dimension", value.kkt_dimension},
        {"kkt_rank", value.kkt_rank},
        {"scaled_kkt_condition_number_2", finite_or_null(value.scaled_kkt_condition_number_2)},
        {"material_balance_tolerance_mol", value.material_balance_tolerance_mol},
        {"electroneutrality_tolerance_mol", value.electroneutrality_tolerance_mol},
        {"pressure_relative_tolerance", value.pressure_relative_tolerance},
        {"chemical_potential_tolerance_over_rt", value.chemical_potential_tolerance_over_rt},
        {"kkt_stationarity_tolerance", value.kkt_stationarity_tolerance},
        {"phase_distinct_tolerance", value.phase_distinct_tolerance},
        {"liquid_compressibility_factor_maximum", value.liquid_compressibility_factor_maximum},
        {"vapor_compressibility_factor_minimum", value.vapor_compressibility_factor_minimum},
        {"reduced_hessian_relative_tolerance", value.reduced_hessian_relative_tolerance},
        {"kkt_condition_number_maximum", value.kkt_condition_number_maximum},
    };
}

Json sensitivity_to_json(const SolvedStateSensitivity& value) {
    return {
        {"status", value.status},
        {"failure_reason", value.failure_reason.empty() ? Json(nullptr) : Json(value.failure_reason)},
        {"derivative_order", 1},
        {"higher_order_status", "unavailable"},
        {"higher_order_reason", "higher_order_total_derivatives_not_implemented"},
        {"tensor_layout", "row-major"},
        {"residual_identities", value.residual_identities},
        {"state_identities", value.state_identities},
        {"output_identities", value.output_identities},
        {"parameter_identities", value.parameter_identities},
        {"residual_state_jacobian", value.residual_state_jacobian},
        {"residual_parameter_jacobian", value.residual_parameter_jacobian},
        {"state_parameter_jacobian", value.state_parameter_jacobian},
        {"output_parameter_jacobian", value.output_parameter_jacobian},
        {"dimension", value.dimension},
        {"parameter_count", value.parameter_identities.size()},
        {"output_count", value.output_identities.size()},
        {"rank", value.rank},
        {"scaled_condition_number_2", finite_or_null(value.scaled_condition_number_2)},
        {"linear_solve_backward_error", finite_or_null(value.linear_solve_backward_error)},
        {"reference_phase_index", value.reference_phase_index},
    };
}

}  // namespace

std::string phase_equilibrium_result_to_json(
    const PhaseEquilibriumResult& result
) {
    Json roles = Json::array();
    for (PhaseRole role : result.problem.constrained_roles) {
        roles.push_back(phase_role_name(role));
    }
    Json payload = {
        {"mode", phase_mode_name(result.problem.mode)},
        {"temperature_k", result.problem.temperature_k},
        {"pressure_pa", result.problem.pressure_pa},
        {"feed_amounts_mol", result.problem.feed_amounts_mol},
        {"constrained_roles", std::move(roles)},
        {"parameter_fingerprint", result.parameter_fingerprint},
        {"phases", phases_to_json(result.phases)},
        {"best_candidate", phases_to_json(result.best_candidate)},
        {"total_free_energy_over_rt", finite_or_null(result.total_free_energy_over_rt)},
        {"diagnostics", diagnostics_to_json(result.diagnostics)},
        {"sensitivity", result.sensitivity ? sensitivity_to_json(*result.sensitivity) : Json(nullptr)},
    };
    return payload.dump() + "\n";
}

}  // namespace epcsaft_engine::equilibrium
