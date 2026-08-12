#pragma once

#include <cstddef>
#include <optional>
#include <string>
#include <vector>

#include "provider.hpp"
#include "solved_state_sensitivity.hpp"

namespace epcsaft_engine::equilibrium {

enum class PhaseRole {
    Liquid,
    Vapor,
};

enum class PhaseSupport {
    AllComponents,
    NeutralComponents,
};

struct DeclaredPhaseSpecification {
    PhaseRole role = PhaseRole::Liquid;
    PhaseSupport support = PhaseSupport::AllComponents;
};

struct PhaseHint {
    std::size_t phase_index = 0;
    std::vector<double> mole_fractions;
    double molar_density_mol_m3 = 0.0;
};

struct FixedTpPhaseProblem {
    double temperature_k = 0.0;
    double pressure_pa = 0.0;
    std::vector<double> feed_amounts_mol;
    std::vector<DeclaredPhaseSpecification> phases;
    std::vector<PhaseHint> hints;
};

struct FixedTCoexistenceProblem {
    double temperature_k = 0.0;
    std::vector<DeclaredPhaseSpecification> phases;
    std::size_t reference_phase_index = 0;
    std::vector<double> reference_mole_fractions;
    std::vector<PhaseHint> hints;
};

struct PhaseEquilibriumState {
    PhaseRole role = PhaseRole::Liquid;
    PhaseSupport support = PhaseSupport::AllComponents;
    std::vector<std::string> supported_component_ids;
    std::vector<double> mole_fractions;
    double molar_volume_m3_per_mol = 0.0;
    double molar_density_mol_m3 = 0.0;
    double packing_fraction = 0.0;
    double pressure_pa = 0.0;
    std::vector<std::string> chemical_potential_component_ids;
    std::vector<double> chemical_potential_over_rt;
    std::string mechanical_class = "unclassified";
};

struct PhaseEquilibriumDiagnostics {
    std::string route;
    std::string solver_status = "not_run";
    std::string numerical_status = "not_adjudicated";
    std::string physical_status = "not_adjudicated";
    std::string failure_stage;
    std::string failure_reason;
    std::string certificate_scope;
    std::string phase_selection_claim;
    std::string globality_claim = "not_established";
    std::string search_scope;
    std::string search_budget;
    std::vector<std::string> search_exclusions;
    std::string completeness_criterion;
    std::vector<std::string> role_evidence;
    int optimizer_solve_count = 0;
    int optimizer_iteration_count = 0;
    int provider_evaluation_count = 0;
    int hint_count = 0;
    int hint_reordered_start_count = 0;
    double material_balance_inf = 0.0;
    double electroneutrality_inf = 0.0;
    double pressure_residual_inf = 0.0;
    double chemical_potential_residual_inf = 0.0;
    double kkt_stationarity_inf = 0.0;
    double dual_sign_violation_inf = 0.0;
    double complementarity_inf = 0.0;
    double phase_distinctness = 0.0;
};

struct FixedTpPhaseResult {
    std::vector<PhaseEquilibriumState> phases;
    std::vector<double> phase_amounts_mol;
    std::vector<double> phase_fractions;
    double total_free_energy_over_rt = 0.0;
    std::string parameter_fingerprint;
    PhaseEquilibriumDiagnostics diagnostics;
    std::optional<SolvedStateSensitivity> sensitivity;
};

struct FixedTCoexistenceResult {
    std::vector<PhaseEquilibriumState> phases;
    std::size_t reference_phase_index = 0;
    double pressure_pa = 0.0;
    std::string parameter_fingerprint;
    PhaseEquilibriumDiagnostics diagnostics;
    std::optional<SolvedStateSensitivity> sensitivity;
};

[[nodiscard]] std::vector<std::string> phase_provider_component_ids(
    const ProviderContext& provider
);

[[nodiscard]] std::vector<double> phase_provider_component_charges(
    const ProviderContext& provider
);

[[nodiscard]] bool phase_topology_is_charged(
    const std::vector<double>& charges
) noexcept;

[[nodiscard]] MixturePhaseEvaluation evaluate_nonreactive_phase(
    const ProviderContext& provider,
    const std::vector<double>& charges,
    double temperature_k,
    const std::vector<double>& mole_fractions,
    double molar_volume_m3_per_mol
);

[[nodiscard]] std::string phase_mechanical_class(
    const MixturePhaseEvaluation& phase,
    double temperature_k,
    double molar_volume_m3_per_mol
);

[[nodiscard]] PhaseEquilibriumState make_phase_equilibrium_state(
    const ProviderContext& provider,
    PhaseRole role,
    PhaseSupport support,
    double temperature_k,
    const std::vector<double>& mole_fractions,
    double molar_volume_m3_per_mol,
    const MixturePhaseEvaluation& phase
);

[[nodiscard]] FixedTpPhaseResult solve_declared_fixed_tp(
    const ProviderContext& provider,
    const FixedTpPhaseProblem& problem
);

[[nodiscard]] FixedTCoexistenceResult solve_fixed_t_coexistence(
    const ProviderContext& reference_provider,
    const ProviderContext& partner_provider,
    const FixedTCoexistenceProblem& problem
);

[[nodiscard]] PhaseEquilibriumState audit_phase_state(
    const ProviderContext& provider,
    PhaseRole role,
    double temperature_k,
    const std::vector<double>& mole_fractions,
    double molar_volume_m3_per_mol
);

[[nodiscard]] const char* phase_role_name(PhaseRole role) noexcept;
[[nodiscard]] const char* phase_support_name(PhaseSupport support) noexcept;

}  // namespace epcsaft_engine::equilibrium
