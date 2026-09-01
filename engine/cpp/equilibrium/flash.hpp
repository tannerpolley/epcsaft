#pragma once

#include <optional>
#include <string>
#include <variant>
#include <vector>

#include "held.hpp"
#include "held2_algorithm.hpp"
#include "eos.hpp"
#include "solved_state_sensitivity.hpp"

namespace epcsaft_engine::equilibrium {

struct FlashResult {
    FlashInput input;
    std::string parameter_fingerprint;
    std::variant<HeldResult, Held2AlgorithmResult> solve;
    std::optional<SolvedStateSensitivity> sensitivity;
    std::size_t eos_evaluation_count = 0;
};

[[nodiscard]] FlashResult solve_tp_flash(
    const epcsaft::NativeModel& eos,
    const FlashInput& input
);

}  // namespace epcsaft_engine::equilibrium
