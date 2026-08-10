#pragma once

#include <optional>
#include <string>
#include <variant>
#include <vector>

#include "held.hpp"
#include "held2_algorithm.hpp"
#include "held2_progress.hpp"
#include "provider.hpp"
#include "solved_state_sensitivity.hpp"

namespace epcsaft_engine::equilibrium {

struct FlashResult {
    FlashInput input;
    std::string parameter_fingerprint;
    std::variant<HeldResult, Held2AlgorithmResult> solve;
    std::optional<SolvedStateSensitivity> sensitivity;
};

[[nodiscard]] FlashResult solve_tp_flash(
    const ProviderContext& provider,
    const FlashInput& input,
    Held2ProgressObserver* observer = nullptr
);

}  // namespace epcsaft_engine::equilibrium
