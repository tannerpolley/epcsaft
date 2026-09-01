#pragma once

#include "phase_stability.hpp"
#include "solved_state_sensitivity.hpp"

#include <optional>
#include <string>

namespace epcsaft_engine::equilibrium {

struct FlashResult {
    FlashInput input;
    std::string parameter_fingerprint;
    PhaseStabilityDiscoveryResult solve;
    std::optional<SolvedStateSensitivity> sensitivity;
};

[[nodiscard]] FlashResult solve_tp_flash(
    const ProviderContext& provider,
    const FlashInput& input
);

}  // namespace epcsaft_engine::equilibrium
