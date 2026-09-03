#pragma once

#include "core/trace.hpp"
#include "equilibrium/problem.hpp"
#include "model_support.hpp"

#include <cstddef>
#include <algorithm>
#include <array>
#include <cstdio>
#include <exception>
#include <string>
#include <vector>

namespace equilibrium_support {

inline int failures = 0;

inline void fail(const char* file, int line, const char* expression) {
    std::fprintf(stderr, "%s:%d: failed: %s\n", file, line, expression);
    ++failures;
}

inline void check(bool condition, const char* file, int line, const char* expression) {
    if (!condition) fail(file, line, expression);
}

inline void print_trace(const epcsaft::core::Trace& trace) {
    for (std::size_t i = 0; i < trace.spans.size(); ++i) {
        const auto& span = trace.spans[i];
        std::fprintf(stderr, "  span[%zu] %s %.0fus\n", i, span.name, span.wall_us);
    }
    for (const auto& event : trace.events) {
        std::fprintf(stderr, "  event span=%u iteration=%u %s", event.span,
                     event.iteration, event.name);
        for (const double value : event.values) std::fprintf(stderr, " %.17g", value);
        std::fputc('\n', stderr);
    }
}

template <class Function>
void run_check(const char* name, Function function) {
    epcsaft::core::Capture capture;
    const int before = failures;
    try {
        function();
    } catch (const std::exception& error) {
        std::fprintf(stderr, "%s: exception: %s\n", name, error.what());
        ++failures;
    }
    if (failures == before) {
        std::printf("%s passed\n", name);
    } else {
        std::fprintf(stderr, "%s failed\n", name);
        print_trace(capture.trace());
    }
}

inline epcsaft::equilibrium::PhaseSpec phase(
    const char* name, epcsaft::equilibrium::PhaseKind kind,
    epcsaft::equilibrium::AmountSpec amount = {}
) {
    epcsaft::equilibrium::PhaseSpec result;
    result.name = name;
    result.kind = kind;
    result.amount = amount;
    return result;
}

inline epcsaft::equilibrium::EquilibriumProblem flash(std::size_t components) {
    epcsaft::equilibrium::EquilibriumProblem result;
    result.phases = {phase("L", epcsaft::equilibrium::PhaseKind::Liquid)};
    result.temperature = epcsaft::equilibrium::FreeValue::fixed(300.0);
    result.pressure = epcsaft::equilibrium::FreeValue::fixed(1.0e5);
    result.feed.kind = epcsaft::equilibrium::FeedKind::Amounts;
    result.feed.values.assign(components, 1.0 / static_cast<double>(components));
    return result;
}

inline epcsaft::equilibrium::EquilibriumProblem boundary() {
    auto result = flash(2);
    result.phases = {
        phase("L", epcsaft::equilibrium::PhaseKind::Liquid),
        phase("V", epcsaft::equilibrium::PhaseKind::Vapor,
              epcsaft::equilibrium::AmountSpec::pinned(0.0)),
    };
    result.pressure = epcsaft::equilibrium::FreeValue::variable(1.0e5);
    return result;
}

inline epcsaft::Parameters ascani_reactive_ternary() {
    auto result = model_support::parameters(3);
    result.segment_count = {2.40, 1.08, 2.80};
    result.sigma_constant = {3.20, 3.00, 3.80};
    result.dispersion_energy_over_k = {200.0, 400.0, 280.0};
    result.molar_mass_kg_per_mol = {0.030, 0.031, 0.059};
    result.k_ij[1] = result.k_ij[3] = -0.045;
    result.k_ij[2] = result.k_ij[6] = 0.045;
    result.k_ij[5] = result.k_ij[7] = -0.025;
    result.site_component = {0, 0, 1, 1};
    result.site_multiplicity = {1, 1, 1, 1};
    result.association_active.assign(16, 0);
    result.association_energy_over_k.assign(16, 0.0);
    result.association_volume.assign(16, 0.0);
    result.association_active[0 * 4 + 1] = result.association_active[1 * 4 + 0] = 1;
    result.association_active[2 * 4 + 3] = result.association_active[3 * 4 + 2] = 1;
    result.association_energy_over_k[0 * 4 + 1] = result.association_energy_over_k[1 * 4 + 0] = 2500.0;
    result.association_energy_over_k[2 * 4 + 3] = result.association_energy_over_k[3 * 4 + 2] = 2500.0;
    result.association_volume[0 * 4 + 1] = result.association_volume[1 * 4 + 0] = 0.05;
    result.association_volume[2 * 4 + 3] = result.association_volume[3 * 4 + 2] = 0.05;
    return result;
}

inline epcsaft::Parameters relabel_neutral_binary() {
    auto result = model_support::neutral_binary();
    std::swap(result.segment_count[0], result.segment_count[1]);
    std::swap(result.sigma_constant[0], result.sigma_constant[1]);
    std::swap(result.dispersion_energy_over_k[0], result.dispersion_energy_over_k[1]);
    std::swap(result.molar_mass_kg_per_mol[0], result.molar_mass_kg_per_mol[1]);
    std::swap(result.k_ij[1], result.k_ij[2]);
    return result;
}

inline epcsaft::equilibrium::EquilibriumProblem ascani_reactive_problem() {
    auto result = flash(3);
    result.phases = {
        phase("L1", epcsaft::equilibrium::PhaseKind::Liquid),
        phase("L2", epcsaft::equilibrium::PhaseKind::Liquid),
    };
    result.phases[0].composition_guess = {0.05, 0.90, 0.05};
    result.phases[1].composition_guess = {0.05, 0.05, 0.90};
    result.reactions = {{"A+B=C", {-1.0, -1.0, 1.0}, 2.25, std::nullopt}};
    result.standard_gibbs_j_per_mol = {0.0, 0.0, 0.0};
    result.feed.values = {0.05, 0.475, 0.475};
    result.temperature = epcsaft::equilibrium::FreeValue::fixed(298.15);
    result.pressure = epcsaft::equilibrium::FreeValue::fixed(1.0e5);
    return result;
}

enum class PhaseSet { OneLiquid, OneVapor, Two, TwoVaporLiquid, Three, ThreeVaporLiquid };
enum class AmountPattern { AllFree, FirstPinned, LastPinned, SplitPinned };
enum class IntensivePattern { BothFixed, TemperatureFree, PressureFree, BothFree };
enum class Chemistry { Neutral, Isomerisation, Electrolyte, TwoReactions };

inline epcsaft::equilibrium::EquilibriumProblem lattice_problem(
    PhaseSet phase_set,
    AmountPattern amounts,
    IntensivePattern intensives,
    Chemistry chemistry
) {
    const std::size_t count = phase_set == PhaseSet::OneLiquid
        || phase_set == PhaseSet::OneVapor ? 1
        : phase_set == PhaseSet::Two || phase_set == PhaseSet::TwoVaporLiquid ? 2 : 3;
    std::size_t components = chemistry == Chemistry::Neutral ? 2 : 3;
    auto result = flash(components);
    result.phases.clear();
    for (std::size_t phase_index = 0; phase_index < count; ++phase_index) {
        const bool vapor = phase_set == PhaseSet::OneVapor
            || (phase_set == PhaseSet::TwoVaporLiquid && phase_index == 0)
            || (phase_set == PhaseSet::ThreeVaporLiquid && phase_index == 0);
        const auto kind = vapor ? epcsaft::equilibrium::PhaseKind::Vapor
                                : epcsaft::equilibrium::PhaseKind::Liquid;
        auto item = phase(
            (std::string("p") + std::to_string(phase_index)).c_str(), kind
        );
        item.name = std::string("p") + std::to_string(phase_index);
        if (amounts == AmountPattern::FirstPinned && phase_index == 0)
            item.amount = epcsaft::equilibrium::AmountSpec::pinned(0.0);
        if (amounts == AmountPattern::LastPinned && phase_index + 1 == count)
            item.amount = epcsaft::equilibrium::AmountSpec::pinned(0.0);
        if (amounts == AmountPattern::SplitPinned && phase_index > 0)
            item.amount = epcsaft::equilibrium::AmountSpec::pinned(0.25);
        result.phases.push_back(item);
    }
    result.temperature = intensives == IntensivePattern::TemperatureFree
        || intensives == IntensivePattern::BothFree
        ? epcsaft::equilibrium::FreeValue::variable(300.0)
        : epcsaft::equilibrium::FreeValue::fixed(300.0);
    result.pressure = intensives == IntensivePattern::PressureFree
        || intensives == IntensivePattern::BothFree
        ? epcsaft::equilibrium::FreeValue::variable(1.0e5)
        : epcsaft::equilibrium::FreeValue::fixed(1.0e5);
    if (chemistry == Chemistry::Isomerisation) {
        result.reactions = {{"isomerisation", {1.0, -1.0, 0.0}, 2.0, std::nullopt}};
    } else if (chemistry == Chemistry::TwoReactions) {
        result.reactions = {
            {"reaction_a", {1.0, -1.0, 0.0}, 1.0, std::nullopt},
            {"reaction_b", {0.0, 1.0, -1.0}, 1.0, std::nullopt},
        };
    } else if (chemistry == Chemistry::Electrolyte) {
        result.feed.values = {0.98, 0.01, 0.01};
        result.feed.kind = epcsaft::equilibrium::FeedKind::Amounts;
    }
    if (chemistry != Chemistry::Electrolyte) {
        result.feed.values.assign(components, 1.0 / static_cast<double>(components));
    }
    return result;
}

inline std::array<PhaseSet, 6> phase_sets() {
    return {PhaseSet::OneLiquid, PhaseSet::OneVapor, PhaseSet::Two,
            PhaseSet::TwoVaporLiquid, PhaseSet::Three, PhaseSet::ThreeVaporLiquid};
}

inline std::array<AmountPattern, 4> amount_patterns() {
    return {AmountPattern::AllFree, AmountPattern::FirstPinned,
            AmountPattern::LastPinned, AmountPattern::SplitPinned};
}

inline std::array<IntensivePattern, 4> intensive_patterns() {
    return {IntensivePattern::BothFixed, IntensivePattern::TemperatureFree,
            IntensivePattern::PressureFree, IntensivePattern::BothFree};
}

inline std::array<Chemistry, 4> chemistry_values() {
    return {Chemistry::Neutral, Chemistry::Isomerisation,
            Chemistry::Electrolyte, Chemistry::TwoReactions};
}

inline int free_intensives(IntensivePattern pattern) {
    return pattern == IntensivePattern::BothFree ? 2
        : pattern == IntensivePattern::BothFixed ? 0 : 1;
}

inline int pinned_phases(
    PhaseSet phase_set, AmountPattern pattern
) {
    const int count = phase_set == PhaseSet::OneLiquid
        || phase_set == PhaseSet::OneVapor ? 1
        : phase_set == PhaseSet::Two || phase_set == PhaseSet::TwoVaporLiquid ? 2 : 3;
    if (pattern == AmountPattern::AllFree) return 0;
    if (pattern == AmountPattern::SplitPinned) return count > 1 ? count - 1 : 0;
    return 1;
}

}  // namespace equilibrium_support

#define CHECK(condition) \
    equilibrium_support::check(static_cast<bool>(condition), __FILE__, __LINE__, #condition)
