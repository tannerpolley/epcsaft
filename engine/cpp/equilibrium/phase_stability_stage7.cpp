#include "phase_stability_stage7.hpp"

namespace epcsaft_engine::equilibrium {

PhaseStabilityStage7Result run_phase_stability_stage7(
    PhaseStabilityPersistentState& state,
    const PhaseStabilityStage5Result& step5,
    const PhaseStabilityStage6Result& step6,
    const PhaseStabilityResourceProfile& resources
) {
    PhaseStabilityStage7Result result;
    result.timing.invocation_count = 1;
    if (step5.status != "complete" || step6.status != "complete") {
        result.reason = "invalid_step7_input";
        return result;
    }
    ++state.major_iteration;
    if (state.major_iteration >= resources.step7_major_iteration_cap) {
        result.reason = "major_iteration_limit";
    } else {
        result.status = "complete";
        result.reason = "continue_candidate_search";
    }
    result.timing.terminal_status = result.status;
    result.timing.terminal_reason = result.reason;
    result.timing.next_step = result.status == "complete" ? 4 : 0;
    return result;
}

}  // namespace epcsaft_engine::equilibrium
