#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <optional>
#include <string>
#include <vector>

#include "phase_stability_certificates.hpp"
#include "phase_stability_tolerances.hpp"

namespace epcsaft_engine::equilibrium {

inline constexpr double kPhaseStabilityPackingFractionMinimum = 1.0e-6;
inline constexpr double kPhaseStabilityPackingFractionMaximum = 0.74;
inline constexpr double kPhaseStabilityModifiedLowerScale = 1.0e-10;

enum class PhaseStabilityCompositionDomain {
    FiniteSearch,
    TraceRefinement,
};

struct PhaseStabilityPolytopeConstraint {
    std::string name;
    std::vector<double> coefficients;
    double upper_bound = 0.0;
};

struct PhaseStabilityCoordinates {
    std::vector<double> charges;
    std::size_t eliminated_index = 0;
    std::size_t dependent_index = 0;
    std::vector<std::size_t> paper_to_provider_indices;
    std::vector<std::size_t> provider_to_paper_indices;
    std::vector<std::size_t> compact_to_paper_indices;
    std::vector<std::size_t> retained_indices;
    std::vector<std::size_t> independent_indices;
    std::vector<double> modified_factors;
    std::vector<double> independent_lower_bounds;
    std::vector<double> independent_upper_bounds;
    std::vector<PhaseStabilityPolytopeConstraint> polytope_constraints;
};

struct PhaseStabilityStageTiming {
    int step = 0;
    int invocation_count = 0;
    double wall_seconds = 0.0;
    double cpu_seconds = 0.0;
    // Provider-backed thermodynamic callback invocations.
    std::uint64_t provider_evaluations = 0;
    std::uint64_t optimizer_solves = 0;
    std::uint64_t optimizer_iterations = 0;
    std::string terminal_status = "not_run";
    std::string terminal_reason = "not_run";
    int next_step = 0;
};

struct PhaseStabilityStateEvaluation;

[[nodiscard]] PhaseStabilityCoordinates make_phase_stability_coordinates(
    const std::vector<double>& charges,
    const std::vector<std::string>& component_ids = {}
);

[[nodiscard]] std::vector<double> phase_stability_transform_physical_fractions(
    const PhaseStabilityCoordinates& coordinates,
    const std::vector<double>& physical_fractions
);

[[nodiscard]] std::vector<double> phase_stability_lift_independent_fractions(
    const PhaseStabilityCoordinates& coordinates,
    const std::vector<double>& independent_modified_fractions,
    PhaseStabilityCompositionDomain domain = PhaseStabilityCompositionDomain::FiniteSearch
);

[[nodiscard]] std::vector<double> phase_stability_physical_fraction_jacobian(
    const PhaseStabilityCoordinates& coordinates
);

[[nodiscard]] std::vector<double> phase_stability_map_unit_cube_to_independent_fractions(
    const PhaseStabilityCoordinates& coordinates,
    const std::vector<double>& unit_cube_coordinates,
    double total_ion_mole_fraction_max
);

[[nodiscard]] std::vector<double> phase_stability_transform_modified_potentials(
    const PhaseStabilityCoordinates& coordinates,
    const std::vector<double>& chemical_potentials
);

struct PhaseStabilityPhysicalPhaseBlock {
    double helmholtz_over_rt = 0.0;
    std::vector<double> gradient;
    std::vector<double> hessian;
    double pressure_pa = 0.0;
};

using PhaseStabilityPhysicalVolumeBoundsEvaluator =
    std::function<std::array<double, 2>(const std::vector<double>&)>;

struct PhaseStabilityStateEvaluation {
    std::vector<double> modified_fractions;
    std::vector<double> physical_amounts;
    double volume = 0.0;
    double objective = 0.0;
    std::vector<double> gradient;
    std::vector<double> hessian;
    std::vector<double> modified_potentials;
    double pressure_stationarity_relative = 0.0;
    double pressure_stationarity_derivative_log_volume = 0.0;
    double helmholtz_over_rt_reference_amount = 0.0;
    double pressure_pa = 0.0;
    std::vector<double> chemical_potentials_over_rt;
};

using PhaseStabilityStateEvaluator = std::function<PhaseStabilityStateEvaluation(
    const std::vector<double>&,
    double
)>;

using PhaseStabilityVolumeBoundsEvaluator = std::function<std::array<double, 2>(
    const std::vector<double>&
)>;

struct PhaseStabilityPressureScanPoint {
    double log_volume = 0.0;
    double volume = 0.0;
    double pressure_residual = 0.0;
    double pressure_derivative_log_volume = 0.0;
    double objective = 0.0;
    bool valid = false;
    std::string failure;
};

struct PhaseStabilityPressureScanInterval {
    double lower_log_volume = 0.0;
    double upper_log_volume = 0.0;
    int depth = 0;
    std::string status;
};

struct PhaseStabilityPressureRoot {
    double log_volume = 0.0;
    double volume = 0.0;
    double objective = 0.0;
    double pressure_residual = 0.0;
    double pressure_derivative_log_volume = 0.0;
    double objective_curvature_log_volume = 0.0;
    std::string mechanical_class = "unclassified";
    std::string origin;
    bool boundary = false;
    PhaseStabilityStateEvaluation state;
};

struct PhaseStabilityPressureEnvelopeResult {
    std::string outcome = "indeterminate";
    std::string failure_reason;
    std::string root_completeness = "not_proven";
    std::string selection_scope =
        "lowest_among_detected_strict_stable_roots";
    int selected_root_index = -1;
    int evaluation_failure_count = 0;
    int refinement_failure_count = 0;
    int stationary_point_count = 0;
    int tangential_root_count = 0;
    int marginal_root_count = 0;
    int boundary_root_count = 0;
    int objective_tie_count = 0;
    int deduplicated_root_count = 0;
    double lower_log_volume = 0.0;
    double upper_log_volume = 0.0;
    std::vector<PhaseStabilityPressureScanPoint> scan_points;
    std::vector<PhaseStabilityPressureScanInterval> intervals;
    std::vector<PhaseStabilityPressureRoot> roots;
};

[[nodiscard]] PhaseStabilityPressureEnvelopeResult evaluate_phase_stability_pressure_envelope(
    const std::vector<double>& independent_modified_fractions,
    const std::array<double, 2>& molar_volume_bounds,
    const PhaseStabilityStateEvaluator& evaluator,
    int initial_interval_count,
    int maximum_subdivision_depth = 8
);

struct PhaseStabilityStageICandidate {
    std::vector<double> modified_fractions;
    double volume = 0.0;
    double tpd = 0.0;
    double reduced_gibbs = 0.0;
    std::vector<double> reduced_gibbs_gradient;
};

struct PhaseStabilityStageIICandidate {
    std::vector<double> independent_modified_fractions;
    double volume = 0.0;
    double phase_coordinate = 0.0;
};

[[nodiscard]] PhaseStabilityStateEvaluation evaluate_phase_stability_phase_block(
    const PhaseStabilityCoordinates& coordinates,
    const std::vector<double>& independent_modified_fractions,
    double log_volume,
    double pressure_over_rt,
    double target_pressure_pa,
    const PhaseStabilityPhysicalPhaseBlock& block,
    PhaseStabilityCompositionDomain domain = PhaseStabilityCompositionDomain::FiniteSearch
);

}  // namespace epcsaft_engine::equilibrium
