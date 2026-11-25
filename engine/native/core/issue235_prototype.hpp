#pragma once

#include "core/export.hpp"

#include <array>
#include <chrono>
#include <cstdint>

namespace epcsaft::issue235 {

enum Lever : unsigned { curvature = 1, eos = 2, refusal = 4, reduced = 8, shadow = 16, meter = 32 };
enum Counter {
    hessian_trial_calls, hessian_trial_ns, hessian_pinned_calls, hessian_pinned_ns,
    hessian_entries, hessian_differences, eos_probes, eos_hits, eos_ineligible,
    eos_checks, eos_differences, eos_bytes, refusal_count, refusal_before_ns,
    refusal_publication_ns, trial_eligible, trial_unrepresentable, trial_fallback,
    trial_ns, trial_full_ns, trial_reduced_iterations, trial_full_iterations,
    hessian_check_refusals, eos_evictions
};

using Clock = std::chrono::steady_clock;
inline std::uint64_t elapsed(Clock::time_point begin) {
    return std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - begin).count();
}
struct Row {
    unsigned mode = 0;
    std::uint64_t generation = 0;
    std::array<std::uint64_t, 24> counts{};
    Clock::time_point decision{};
};
extern thread_local Row row;
inline bool enabled(Lever lever) { return row.mode & lever; }
}
