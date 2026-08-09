#pragma once

#include <cstddef>
#include <limits>
#include <optional>
#include <string>
#include <vector>

#include "provider.hpp"
#include "solved_state_sensitivity.hpp"

namespace epcsaft_engine::equilibrium {

enum class PhaseEquilibriumMode {
    Discovery,
    Constrained,
};

enum class PhaseRole {
    Liquid,
    Vapor,
};

struct PhaseHint {
    PhaseRole role = PhaseRole::Liquid;
    std::vector<double> mole_fractions;
    double molar_density_mol_m3 = 0.0;
};

struct PhaseEquilibriumProblem {
    double temperature_k = 0.0;
    double pressure_pa = 0.0;
    std::vector<double> feed_amounts_mol;
    PhaseEquilibriumMode mode = PhaseEquilibriumMode::Discovery;
    std::vector<PhaseRole> constrained_roles;
    std::vector<PhaseHint> hints;
};

struct PhaseEquilibriumState {
    PhaseRole role = PhaseRole::Liquid;
    double amount_mol = 0.0;
    std::vector<double> amounts_mol;
    std::vector<double> mole_fractions;
    double volume_m3 = 0.0;
    double molar_density_mol_m3 = 0.0;
    double packing_fraction = 0.0;
    double pressure_pa = 0.0;
    std::vector<double> chemical_potential_over_rt;
};

struct PhaseEquilibriumDiagnostics {
    std::string solver_status = "not_evaluated";
    std::string numerical_status = "not_evaluated";
    std::string physical_status = "not_evaluated";
    std::string failure_stage;
    std::string failure_reason;
    std::string failure_detail;
    std::string certificate_scope = "none";
    std::string phase_selection_claim = "none";
    std::string globality_claim = "not_established";
    std::string completeness_criterion;
    std::vector<std::string> searched_topologies;
    std::vector<std::string> excluded_topologies;
    std::vector<std::string> equation_state_identities;
    std::vector<std::string> equation_residual_identities;
    std::string tensor_layout = "row-major";
    std::string solver_derivative_contract =
        "CppAD state gradient and state-state Hessian";
    std::string sensitivity_derivative_contract =
        "CppAD state-parameter partials plus exact implicit linear algebra";
    std::string parameter_parameter_requirement = "not_consumed";
    int topology_budget = 0;
    int topology_count = 0;
    int start_budget_per_topology = 0;
    int attempted_start_count = 0;
    int certified_candidate_count = 0;
    int optimizer_iteration_count = 0;
    int hint_count = 0;
    int hint_reordered_start_count = 0;
    double material_balance_inf = 0.0;
    double electroneutrality_inf = 0.0;
    double pressure_residual_inf = 0.0;
    double chemical_potential_residual_inf = 0.0;
    double kkt_stationarity_inf = 0.0;
    double minimum_reduced_hessian_eigenvalue = 0.0;
    std::size_t kkt_dimension = 0;
    std::size_t kkt_rank = 0;
    double scaled_kkt_condition_number_2 = 0.0;
    double material_balance_tolerance_mol = 0.0;
    double electroneutrality_tolerance_mol = 0.0;
    double pressure_relative_tolerance = 0.0;
    double chemical_potential_tolerance_over_rt = 0.0;
    double kkt_stationarity_tolerance = 0.0;
    double phase_distinct_tolerance = 0.0;
    double liquid_compressibility_factor_maximum = 0.0;
    double vapor_compressibility_factor_minimum = 0.0;
    double reduced_hessian_relative_tolerance = 0.0;
    double kkt_condition_number_maximum = 0.0;
};

struct PhaseEquilibriumResult {
    PhaseEquilibriumProblem problem;
    std::string parameter_fingerprint;
    std::vector<PhaseEquilibriumState> phases;
    std::vector<PhaseEquilibriumState> best_candidate;
    double total_free_energy_over_rt = std::numeric_limits<double>::quiet_NaN();
    PhaseEquilibriumDiagnostics diagnostics;
    std::optional<SolvedStateSensitivity> sensitivity;
};

[[nodiscard]] PhaseEquilibriumResult solve_phase_equilibrium(
    const ProviderContext& provider,
    const PhaseEquilibriumProblem& problem
);

[[nodiscard]] const char* phase_role_name(PhaseRole role) noexcept;
[[nodiscard]] const char* phase_mode_name(PhaseEquilibriumMode mode) noexcept;

}  // namespace epcsaft_engine::equilibrium
