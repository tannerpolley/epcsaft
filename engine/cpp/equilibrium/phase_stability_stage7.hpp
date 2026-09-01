#pragma once

#include "phase_stability_stage6.hpp"

namespace epcsaft_engine::equilibrium {

struct PhaseStabilityStage7Result {
    std::string status = "indeterminate";
    std::string reason = "not_run";
    PhaseStabilityStageTiming timing;
};

[[nodiscard]] PhaseStabilityStage7Result run_phase_stability_stage7(
    PhaseStabilityPersistentState& state,
    const PhaseStabilityStage5Result& step5,
    const PhaseStabilityStage6Result& step6,
    const PhaseStabilityResourceProfile& resources
);

}  // namespace epcsaft_engine::equilibrium
