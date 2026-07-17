#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace epcsaft::detail {

struct DensityRoot final {
    double molar_density;
    double bracket_lower;
    double bracket_upper;
    double pressure_residual;
    double pressure_density_derivative;
    std::size_t iterations;
    bool stable;
    std::string branch;
};

template <typename Evaluate>
DensityRoot refine_density_root(
    double pressure,
    double bracket_lower,
    double bracket_upper,
    double residual_lower,
    Evaluate&& evaluate
) {
    constexpr std::size_t max_iterations = 100;
    const double pressure_tolerance = std::max(1.0e-6, 1.0e-10 * pressure);
    std::size_t iteration = 0;
    for (; iteration < max_iterations; ++iteration) {
        const double midpoint = 0.5 * (bracket_lower + bracket_upper);
        const auto state = evaluate(midpoint);
        const double residual = state.pressure - pressure;
        if (!std::isfinite(residual)) {
            throw std::runtime_error("density residual is nonfinite");
        }
        if (residual_lower * residual <= 0.0) {
            bracket_upper = midpoint;
        } else {
            bracket_lower = midpoint;
            residual_lower = residual;
        }
        if (std::abs(residual) <= pressure_tolerance
            && (bracket_upper - bracket_lower)
                <= 1.0e-10 * std::max(1.0, midpoint)) {
            break;
        }
        const double slope = state.pressure_density_derivative;
        if (std::isfinite(slope) && slope != 0.0) {
            const double candidate = midpoint - residual / slope;
            if (candidate > bracket_lower && candidate < bracket_upper) {
                const auto candidate_state = evaluate(candidate);
                const double candidate_residual = candidate_state.pressure - pressure;
                if (std::isfinite(candidate_residual)) {
                    if (residual_lower * candidate_residual <= 0.0) {
                        bracket_upper = candidate;
                    } else {
                        bracket_lower = candidate;
                        residual_lower = candidate_residual;
                    }
                }
            }
        }
    }
    if (iteration == max_iterations) {
        throw std::runtime_error("density root refinement did not converge");
    }
    const double density = 0.5 * (bracket_lower + bracket_upper);
    const auto state = evaluate(density);
    const double residual = state.pressure - pressure;
    if (std::abs(residual) > 10.0 * pressure_tolerance) {
        throw std::runtime_error("density root pressure residual is too large");
    }
    const double slope = state.pressure_density_derivative;
    return {
        density,
        bracket_lower,
        bracket_upper,
        residual,
        slope,
        iteration + 1,
        std::isfinite(slope) && slope > 0.0,
        "",
    };
}

template <typename Evaluate>
DensityRoot solve_density(
    double pressure,
    double ideal_density,
    double maximum_density,
    const std::string& phase,
    Evaluate&& evaluate
) {
    if (!(pressure > 0.0) || !(maximum_density > 0.0)
        || !std::isfinite(pressure) || !std::isfinite(maximum_density)) {
        throw std::invalid_argument("invalid density solve bounds");
    }
    if (phase != "" && phase != "vapor" && phase != "liquid") {
        throw std::invalid_argument("invalid density phase");
    }
    const double lower = std::max(1.0e-12, ideal_density * 1.0e-8);
    const double upper = maximum_density * 0.999;
    if (!(lower < upper)) {
        throw std::runtime_error("density interval is empty");
    }

    constexpr std::size_t scan_intervals = 360;
    const double log_lower = std::log(lower);
    const double log_upper = std::log(upper);
    std::vector<DensityRoot> roots;
    double previous_density = lower;
    auto previous_state = evaluate(previous_density);
    double previous_residual = previous_state.pressure - pressure;
    for (std::size_t index = 1; index <= scan_intervals; ++index) {
        const double fraction = static_cast<double>(index) / scan_intervals;
        const double density = std::exp(
            log_lower + fraction * (log_upper - log_lower)
        );
        const auto state = evaluate(density);
        const double residual = state.pressure - pressure;
        if (std::isfinite(previous_residual) && std::isfinite(residual)
            && previous_residual * residual <= 0.0) {
            const auto root = refine_density_root(
                pressure,
                previous_density,
                density,
                previous_residual,
                evaluate
            );
            if (roots.empty()
                || std::abs(root.molar_density - roots.back().molar_density)
                    > 1.0e-8 * std::max(1.0, root.molar_density)) {
                roots.push_back(root);
            }
        }
        previous_density = density;
        previous_residual = residual;
    }
    std::vector<DensityRoot> stable;
    std::copy_if(
        roots.begin(), roots.end(), std::back_inserter(stable),
        [](const DensityRoot& root) { return root.stable; }
    );
    if (stable.empty()) {
        throw std::runtime_error("no certified stable density root");
    }
    if (stable.size() == 1) {
        stable.front().branch = "single";
        return stable.front();
    }
    if (phase.empty()) {
        throw std::runtime_error("density root is ambiguous");
    }
    if (phase == "vapor") {
        stable.front().branch = "vapor";
        return stable.front();
    }
    stable.back().branch = "liquid";
    return stable.back();
}

}  // namespace epcsaft::detail
