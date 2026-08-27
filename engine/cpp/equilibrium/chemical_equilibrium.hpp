#pragma once

#include <array>
#include <cstddef>
#include <limits>
#include <string>
#include <vector>

namespace epcsaft_engine::equilibrium {

inline constexpr double kChemicalBalanceTolerance = 1.0e-9;
inline constexpr double kChemicalPressureTolerance = 1.0e-8;
// A 1e-6 dimensionless affinity is about 2.48e-3 J/mol at 298.15 K.
// That is below the meaningful resolution of configured reactive phase
// calculations while remaining substantially tighter than process data.
inline constexpr double kChemicalAffinityTolerance = 1.0e-6;
inline constexpr double kChemicalKktTolerance = 1.0e-7;
inline constexpr double kChemicalChartPullbackTolerance = 1.0e-10;
inline constexpr double kChemicalObjectiveDerivativeTolerance = 1.0e-5;
inline constexpr double kChemicalConstraintDerivativeTolerance = 1.0e-5;
inline constexpr double kChemicalHessianDerivativeTolerance = 2.0e-4;
inline constexpr double kChemicalConditionNumberMax = 1.0e6;
inline constexpr double kPhaseDirectionRankTolerance = 1.0e-12;
inline constexpr double kPhaseDirectionConservationTolerance = 1.0e-12;

struct DenseMatrix {
    std::size_t rows = 0;
    std::size_t columns = 0;
    std::vector<double> values;

    [[nodiscard]] double operator()(std::size_t row, std::size_t column) const;
    double& operator()(std::size_t row, std::size_t column);
};

// A source-neutral, phase/species direction compilation.  The columns of
// ``direction_matrix`` are ordered phase-major/species-minor.  ``invariant_matrix``
// contains an orthonormal row basis for the complement of the declared
// directions; its rows therefore represent the conserved coordinates left
// after the caller has declared reactions and transfers.  This is deliberately
// a compiler seam only: it does not own an equilibrium solve.
struct PhaseDirectionCompilation {
    std::vector<std::string> phase_ids;
    std::vector<std::string> species_ids;
    std::vector<double> charges;
    DenseMatrix balance_matrix;
    DenseMatrix direction_matrix;
    DenseMatrix invariant_matrix;
    std::vector<double> phase_reference_amounts;
    std::vector<double> driving_constants;
    std::vector<double> invariant_totals;
    std::vector<double> reference_potentials;
    std::vector<double> global_balance_residuals;
    std::vector<double> phase_charge_residuals;
    std::vector<double> reference_phase_charge_residuals;
    std::size_t direction_rank = 0;
    std::size_t invariant_rank = 0;
    double direction_scale = 0.0;
    double rank_tolerance = 0.0;
    double conservation_tolerance = 0.0;
    double invariant_orthonormal_residual_inf = 0.0;
    double invariant_direction_residual_inf = 0.0;
    double potential_equation_residual_inf = 0.0;
    double minimum_norm_potential_norm = 0.0;
    double global_balance_residual_inf = 0.0;
    double phase_charge_residual_inf = 0.0;
    double reference_neutrality_residual_inf = 0.0;
};

struct EquilibriumConstantRecord {
    std::string source_id;
    std::string reference_id;
    std::string reaction_orientation;
    std::string conversion_id;
    bool dimensionless = false;
    double temperature_k = 0.0;
    double pressure_pa = 0.0;
};

struct ReactionSystemInput {
    std::vector<std::string> species_ids;
    std::vector<int> charges;
    std::string provider_fingerprint;
    std::vector<double> molar_masses_kg_per_mol;
    DenseMatrix balance_matrix;
    std::vector<double> conserved_totals;
    DenseMatrix reaction_matrix;
    std::vector<double> feed_amounts;
    std::vector<double> ln_k;
    std::vector<EquilibriumConstantRecord> equilibrium_constant_records;
    double temperature_k = 0.0;
    double pressure_pa = 0.0;
};

struct CompiledReactionSystem {
    bool ideal_gas = false;
    std::string provider_fingerprint;
    std::size_t original_species_count = 0;
    std::vector<int> original_charges;
    std::vector<double> original_molar_masses_kg_per_mol;
    std::vector<double> original_feed_amounts;
    std::size_t species_count = 0;
    std::vector<int> charges;
    std::vector<double> molar_masses_kg_per_mol;
    DenseMatrix supplied_balance_matrix;
    std::vector<std::size_t> retained_species_indices;
    std::vector<std::size_t> removed_species_indices;
    DenseMatrix balance_matrix;
    DenseMatrix reaction_matrix;
    DenseMatrix supplied_reaction_transform;
    std::vector<double> balance_totals;
    std::vector<double> feed_amounts;
    std::vector<double> g_ref;
    std::vector<double> g_ref_pressure_derivatives_per_pa;
    std::vector<double> g_ref_parameter_derivatives;
};

struct AmountChart {
    std::vector<int> charges;
    std::vector<std::size_t> cation_indices;
    std::vector<std::size_t> anion_indices;
    std::vector<std::size_t> neutral_indices;

    [[nodiscard]] std::size_t coordinate_count() const;
    [[nodiscard]] bool ionic() const;
};

struct AmountChartEvaluation {
    std::vector<double> amounts;
    std::vector<double> jacobian;
    std::vector<double> amount_hessians;
    double minimum_amount = 0.0;
    double charge_residual = 0.0;
};

struct ReactiveFinitePhaseStart {
    std::vector<double> amounts_mol;
    double molar_volume_m3_per_mol = std::numeric_limits<double>::quiet_NaN();
};

struct ReactiveFinitePhaseOutput {
    std::string selector;
    std::size_t phase_index = 0;
    std::vector<double> coefficients;
    std::vector<double> solvent_mass_coefficients_kg_per_mol;
};

struct ReactiveFinitePhaseResult {
    bool accepted = false;
    std::string solver_status;
    std::string callback_error;
    std::string callback_failure_kind = "none";
    std::string failure_reason;
    std::string numerical_status = "not_adjudicated";
    std::string physical_status = "not_adjudicated";
    std::string provider_domain_status = "not_adjudicated";
    std::vector<std::pair<std::string, double>> resolved_policy;
    std::vector<std::vector<double>> phase_amounts_mol;
    std::vector<double> phase_volumes_m3;
    std::vector<double> phase_pressures_pa;
    std::vector<double> phase_packing_fractions;
    std::vector<std::vector<double>> phase_chemical_potentials_over_rt;
    std::vector<std::vector<double>> phase_log_fugacities;
    std::vector<double> local_variables;
    std::vector<double> local_lower_bounds;
    std::vector<double> local_upper_bounds;
    std::vector<std::string> input_order;
    std::vector<double> kkt_input_jacobian;
    std::vector<double> input_amount_derivatives;
    std::vector<double> input_volume_derivatives;
    std::vector<double> input_pressure_derivatives;
    std::vector<double> input_chemical_potential_derivatives;
    std::vector<double> input_log_fugacity_derivatives;
    std::vector<double> kkt_root_values;
    std::vector<double> kkt_root_jacobian;
    std::size_t kkt_dimension = 0;
    std::size_t kkt_rank = 0;
    double condition_number_inf = std::numeric_limits<double>::infinity();
    std::vector<double> output_values;
    std::vector<double> output_input_derivatives;
    std::vector<double> output_local_derivatives;
    std::vector<double> output_parameter_partials;
    double objective_over_rt = std::numeric_limits<double>::quiet_NaN();
    double balance_inf_norm = std::numeric_limits<double>::quiet_NaN();
    double pressure_relative_inf_norm = std::numeric_limits<double>::quiet_NaN();
    double reaction_affinity_inf_norm = std::numeric_limits<double>::quiet_NaN();
    double transfer_chemical_potential_inf_norm =
        std::numeric_limits<double>::quiet_NaN();
    double kkt_stationarity_inf_norm = std::numeric_limits<double>::quiet_NaN();
    double phase_distinctness = std::numeric_limits<double>::quiet_NaN();
    double minimum_amount_mol = std::numeric_limits<double>::quiet_NaN();
    double minimum_bound_margin = std::numeric_limits<double>::quiet_NaN();
    double reduced_hessian_min_eigenvalue =
        std::numeric_limits<double>::quiet_NaN();
    double reduced_hessian_scale = std::numeric_limits<double>::quiet_NaN();
    double reduced_hessian_threshold = std::numeric_limits<double>::quiet_NaN();
    std::size_t declared_direction_rank = 0;
    std::size_t admissible_nullity = 0;
    std::size_t optimizer_iterations = 0;
    std::size_t provider_evaluation_count = 0;
    bool derivative_audits_performed = false;
    std::string derivative_audit_status = "not_requested";
    double objective_gradient_audit_max_relative_error =
        std::numeric_limits<double>::quiet_NaN();
    double constraint_jacobian_audit_max_relative_error =
        std::numeric_limits<double>::quiet_NaN();
    double lagrangian_hessian_audit_max_relative_error =
        std::numeric_limits<double>::quiet_NaN();
    std::vector<std::string> parameter_order;
    std::string sensitivity_status = "not_requested";
    std::string sensitivity_failure_reason;
    std::vector<double> amount_derivatives;
    std::vector<double> volume_derivatives;
    std::vector<double> pressure_derivatives;
    std::vector<double> chemical_potential_derivatives;
    std::vector<double> log_fugacity_derivatives;
};


class ProviderContext;
struct ProviderActiveParameterSet;

struct SourceReferenceTransferEvidence {
    std::string status;
    std::string domain_status;
    std::string convergence_status;
    std::string reference_state_id;
    std::string activity_convention_id;
    std::vector<std::string> component_ids;
    std::size_t neutral_basis_row_count = 0;
    std::vector<double> neutral_basis;
    std::vector<double> log_fugacity_contractions;
    std::vector<double> activity_scale_log_contractions;
    std::vector<double> transfer_log_contractions;
    std::vector<double> source_composition;
    std::vector<std::string> derivative_availability;
    std::vector<double> pressure_derivatives_per_pa;
    std::vector<std::string> active_parameter_identities;
    std::vector<double> parameter_derivatives;
    double temperature_k = 0.0;
    double pressure_pa = 0.0;
    double source_reference_pressure_pa = 0.0;
    double standard_molality_mol_per_kg = 0.0;
    double reference_convergence_error = 0.0;
    double source_temperature_min_k = 0.0;
    double source_temperature_max_k = 0.0;
    double source_pressure_min_pa = 0.0;
    double source_pressure_max_pa = 0.0;
    std::string parameter_fingerprint;
    std::string topology_fingerprint;
    std::string component_order_fingerprint;
    std::string reference_state_fingerprint;
    std::string domain_fingerprint;
    std::string artifact_fingerprint;
    std::string helmholtz_basis_id;
};

struct SourceStandardStateResult {
    std::vector<double> standard_offsets;
    std::vector<double> ln_k_provider_basis;
    std::vector<double> ln_k_pressure_derivatives_per_pa;
    std::vector<double> ln_k_parameter_derivatives;
    double representation_residual_inf_norm = 0.0;
};

[[nodiscard]] CompiledReactionSystem compile_reaction_system(
    const ReactionSystemInput& input,
    bool require_complete_species_basis = true
);

[[nodiscard]] std::vector<double> reconstruct_compiled_reference(
    const CompiledReactionSystem& system,
    const std::vector<double>& supplied_ln_k
);

[[nodiscard]] PhaseDirectionCompilation compile_phase_direction_system(
    const std::vector<std::string>& phase_ids,
    const std::vector<std::string>& species_ids,
    const std::vector<double>& charges,
    const DenseMatrix& balance_matrix,
    const DenseMatrix& phase_reference_amounts,
    const DenseMatrix& direction_matrix,
    const std::vector<double>& driving_constants,
    double rank_tolerance = kPhaseDirectionRankTolerance,
    double conservation_tolerance = kPhaseDirectionConservationTolerance
);

[[nodiscard]] AmountChart make_amount_chart(const std::vector<int>& charges);

[[nodiscard]] AmountChartEvaluation evaluate_amount_chart(
    const AmountChart& chart,
    const std::vector<double>& coordinates
);

[[nodiscard]] std::vector<double> invert_amount_chart(
    const AmountChart& chart,
    const std::vector<double>& amounts
);

[[nodiscard]] std::vector<std::size_t> homogeneous_structural_zeros(
    const DenseMatrix& balance_matrix,
    const std::vector<double>& balance_totals,
    const std::vector<int>& charges
);


[[nodiscard]] SourceStandardStateResult transform_source_standard_state(
    const DenseMatrix& reaction_matrix,
    const std::vector<double>& source_ln_k,
    const std::vector<double>& log_activity_scale_factors,
    const std::vector<std::string>& species_ids,
    const std::vector<int>& charges,
    const std::string& provider_fingerprint,
    double temperature_k,
    double pressure_pa,
    const SourceReferenceTransferEvidence& reference
);


[[nodiscard]] ReactiveFinitePhaseResult solve_provider_reactive_finite_phases(
    const CompiledReactionSystem& system,
    const ProviderContext& provider,
    double temperature_k,
    double pressure_pa,
    const std::vector<std::array<double, 2>>& packing_fraction_intervals,
    const std::vector<ReactiveFinitePhaseStart>& starts,
    double trace_floor,
    const ProviderActiveParameterSet* active_parameters = nullptr,
    int max_iterations = 500,
    bool perform_derivative_audits = false,
    const PhaseDirectionCompilation* phase_directions = nullptr,
    const std::vector<std::vector<std::size_t>>& phase_support_indices = {},
    const std::vector<ReactiveFinitePhaseOutput>& outputs = {}
);

[[nodiscard]] ReactiveFinitePhaseResult evaluate_provider_reactive_finite_phases_local(
    const CompiledReactionSystem& system,
    const ProviderContext& provider,
    double temperature_k,
    double pressure_pa,
    const std::vector<std::array<double, 2>>& packing_fraction_intervals,
    const std::vector<double>& local_variables,
    const std::vector<double>& local_lower_bounds,
    const std::vector<double>& local_upper_bounds,
    const ProviderActiveParameterSet* active_parameters = nullptr,
    const PhaseDirectionCompilation* phase_directions = nullptr,
    const std::vector<std::vector<std::size_t>>& phase_support_indices = {},
    const std::vector<ReactiveFinitePhaseOutput>& outputs = {}
);

[[nodiscard]] ReactiveFinitePhaseResult finalize_provider_reactive_finite_phases_local(
    const CompiledReactionSystem& system,
    const ProviderContext& provider,
    double temperature_k,
    double pressure_pa,
    const std::vector<std::array<double, 2>>& packing_fraction_intervals,
    const std::vector<double>& local_variables,
    const std::vector<double>& local_lower_bounds,
    const std::vector<double>& local_upper_bounds,
    const ProviderActiveParameterSet* active_parameters = nullptr,
    const PhaseDirectionCompilation* phase_directions = nullptr,
    const std::vector<std::vector<std::size_t>>& phase_support_indices = {},
    const std::vector<ReactiveFinitePhaseOutput>& outputs = {},
    double trace_floor = 0.0,
    bool perform_derivative_audits = false
);

}  // namespace epcsaft_engine::equilibrium
