#pragma once

#include "phase_stability_stage5.hpp"

namespace epcsaft_engine::equilibrium {

struct PhaseStabilityCandidateDecision {
    std::uint64_t insertion_id = 0;
    bool gap_passed = false;
    bool derivative_passed = false;
    bool pairwise_distinct = false;
    bool retained = false;
    std::string reason;
};

struct PhaseStabilityStage6Result {
    std::string status = "indeterminate";
    std::string reason = "not_run";
    std::vector<PhaseStabilityMPoint> candidates;
    std::vector<PhaseStabilityCandidateDecision> decisions;
    PhaseStabilityStageTiming timing;
};

[[nodiscard]] PhaseStabilityStage6Result run_phase_stability_stage6(
    const PhaseStabilityStage1Result& step1,
    const PhaseStabilityStage4Result& step4,
    const PhaseStabilityPersistentState& state,
    const PhaseStabilityPackingFractionEvaluator& packing_fraction
);

}  // namespace epcsaft_engine::equilibrium
