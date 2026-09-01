#pragma once

#include <functional>
#include <cstdint>
#include <optional>
#include <string>
#include <vector>

#include "phase_stability_stage4.hpp"
#include "phase_stability_certificates.hpp"

namespace epcsaft_engine::equilibrium {

using PhaseStabilityPackingFractionEvaluator = std::function<double(
    const std::vector<double>&,
    double
)>;

struct PhaseStabilityResourceProfile {
    int step2_search_work_budget = 6500;
    int step5_start_cap = 128;
    int step7_major_iteration_cap = 128;
    int step5_recovery_start_cap = 2048;
};

struct PhaseStabilityDiluteFaceRestartEvidence {
    double face_radius = 0.0;
    double bound_activity_atol = 0.0;
    double dual_sign_atol = 0.0;
    double complementarity_atol = 0.0;
    std::vector<std::size_t> coordinate_indices;
    std::vector<std::size_t> provider_component_indices;
    std::vector<double> lower_bound_distances;
    std::vector<double> lower_bound_multipliers;
    std::vector<double> complementarity_products;
    std::string solver_status = "not_run";
    bool attempted = false;
    bool accepted = false;
};

struct PhaseStabilityLocalCertificate {
    std::uint64_t start_ordinal = 0;
    std::string start_family;
    std::string solver_status;
    std::optional<double> local_value;
    bool finite_and_in_domain = false;
    std::optional<PhaseStabilityStage5KktCertificate> kkt;
    std::optional<PhaseStabilityDiluteFaceRestartEvidence> dilute_face_restart;
    bool accepted = false;
};

struct PhaseStabilityStage5Result {
    std::string status = "indeterminate";
    std::string reason = "not_run";
    std::optional<double> lower_value;
    std::optional<PhaseStabilityMPoint> terminal;
    std::uint64_t starts_consumed = 0;
    std::vector<PhaseStabilityLocalCertificate> attempts;
    PhaseStabilityStageTiming timing;
};

[[nodiscard]] bool retain_phase_stability_m_point(
    PhaseStabilityPersistentState& state,
    PhaseStabilityMPoint& point
);

[[nodiscard]] PhaseStabilityStage5Result run_phase_stability_stage5(
    const PhaseStabilityStage1Result& step1,
    const PhaseStabilityStage4Result& step4,
    PhaseStabilityPersistentState& state,
    const PhaseStabilityStateEvaluator& evaluator,
    const PhaseStabilityResourceProfile& resources
);

}  // namespace epcsaft_engine::equilibrium
