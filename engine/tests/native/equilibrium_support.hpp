#pragma once

#include "core/trace.hpp"
#include "equilibrium/problem.hpp"
#include "properties_support.hpp"
#include "recording_support.hpp"

#include <cstddef>
#include <algorithm>
#include <array>
#include <cstdio>
#include <exception>
#include <string>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace equilibrium_support {

inline int failures = 0;

// Record labels of declaration enums.
inline const char* phase_kind_name(epcsaft::equilibrium::PhaseKind kind) noexcept {
    return kind == epcsaft::equilibrium::PhaseKind::Liquid ? "liquid" : "vapor";
}

inline const char* amount_kind_name(epcsaft::equilibrium::AmountKind kind) noexcept {
    return kind == epcsaft::equilibrium::AmountKind::Free ? "free" : "pinned";
}

inline const char* feed_kind_name(epcsaft::equilibrium::FeedKind kind) noexcept {
    using epcsaft::equilibrium::FeedKind;
    switch (kind) {
    case FeedKind::MoleFractions: return "mole_fractions";
    case FeedKind::Amounts: return "amounts";
    case FeedKind::MassFractions: return "mass_fractions";
    case FeedKind::Molalities: return "molalities";
    }
    return "unknown";
}

// nu_r . g for each reaction; zero affinities without a standard state.
inline std::vector<double> reaction_affinities(
    const std::vector<epcsaft::equilibrium::ReactionSpec>& reactions, const std::vector<double>& standard_gibbs
) {
    std::vector<double> result(reactions.size(), 0.0);
    if (standard_gibbs.empty()) return result;
    for (std::size_t r = 0; r < reactions.size(); ++r) {
        if (reactions[r].stoichiometry.size() != standard_gibbs.size())
            throw std::invalid_argument("stoichiometry size does not match standard state");
        result[r] = std::inner_product(reactions[r].stoichiometry.begin(), reactions[r].stoichiometry.end(),
                                       standard_gibbs.begin(), 0.0);
    }
    return result;
}

inline void fail(const char* file, int line, const char* expression) {
    std::fprintf(stderr, "%s:%d: failed: %s\n", file, line, expression);
    ++failures;
}

inline void check(bool condition, const char* file, int line, const char* expression) {
    if (!condition) {
        fail(file, line, expression);
        throw recording_support::Failure(
            std::string(file) + ":" + std::to_string(line)
            + ": failed: " + expression
        );
    }
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

// Synthetic three-species associating parameters; manufactured numerical evidence, not
// literature data.
inline epcsaft::Parameters synthetic_reactive_ternary() {
    auto result = properties_support::parameters(3);
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
    auto result = properties_support::neutral_binary();
    std::swap(result.segment_count[0], result.segment_count[1]);
    std::swap(result.sigma_constant[0], result.sigma_constant[1]);
    std::swap(result.dispersion_energy_over_k[0], result.dispersion_energy_over_k[1]);
    std::swap(result.molar_mass_kg_per_mol[0], result.molar_mass_kg_per_mol[1]);
    std::swap(result.k_ij[1], result.k_ij[2]);
    return result;
}

inline epcsaft::equilibrium::EquilibriumProblem synthetic_reactive_problem() {
    auto result = flash(3);
    result.phases = {
        phase("L1", epcsaft::equilibrium::PhaseKind::Liquid),
        phase("L2", epcsaft::equilibrium::PhaseKind::Liquid),
    };
    result.phases[0].composition_guess = {0.05, 0.90, 0.05};
    result.phases[1].composition_guess = {0.05, 0.05, 0.90};
    // Manufactured numerical case: supplied g=0 governs, so the consistent target is K=1.
    result.reactions = {{"A+B=C", {-1.0, -1.0, 1.0}, 1.0, std::nullopt}};
    result.standard_gibbs_j_per_mol = {0.0, 0.0, 0.0};
    result.feed.values = {0.05, 0.475, 0.475};
    result.temperature = epcsaft::equilibrium::FreeValue::fixed(298.15);
    result.pressure = epcsaft::equilibrium::FreeValue::fixed(1.0e5);
    return result;
}

enum class PhaseSet { OneLiquid, OneVapor, Two, TwoVaporLiquid, Three, ThreeVaporLiquid };
enum class AmountPattern { AllFree, FirstPinned, LastPinned, SplitPinned, AllPinned };
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
        // A single declared phase pinned alone is an all-pinned row; its admitted fraction is 1.
        const double lone_pin = count == 1 ? 1.0 : 0.0;
        if (amounts == AmountPattern::FirstPinned && phase_index == 0)
            item.amount = epcsaft::equilibrium::AmountSpec::pinned(lone_pin);
        if (amounts == AmountPattern::LastPinned && phase_index + 1 == count)
            item.amount = epcsaft::equilibrium::AmountSpec::pinned(lone_pin);
        if (amounts == AmountPattern::SplitPinned && phase_index > 0)
            item.amount = epcsaft::equilibrium::AmountSpec::pinned(0.25);
        if (amounts == AmountPattern::AllPinned)
            item.amount = epcsaft::equilibrium::AmountSpec::pinned(1.0 / static_cast<double>(count));
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

inline std::array<AmountPattern, 5> amount_patterns() {
    return {AmountPattern::AllFree, AmountPattern::FirstPinned,
            AmountPattern::LastPinned, AmountPattern::SplitPinned, AmountPattern::AllPinned};
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
    if (pattern == AmountPattern::AllPinned) return count;
    if (pattern == AmountPattern::SplitPinned) return count > 1 ? count - 1 : 0;
    return 1;
}

// Reorder every per-component record so species identity, not vector position, is tested.
inline epcsaft::Parameters permute_parameters(
    const epcsaft::Parameters& source, const std::vector<std::size_t>& order
) {
    const std::size_t n = order.size();
    auto result = source;
    const auto permute = [&](auto& target, const auto& origin) {
        for (std::size_t i = 0; i < n; ++i) target[i] = origin[order[i]];
    };
    const auto permute_matrix = [&](auto& target, const auto& origin) {
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < n; ++j) target[i * n + j] = origin[order[i] * n + order[j]];
    };
    permute(result.segment_count, source.segment_count);
    permute(result.sigma_constant, source.sigma_constant);
    permute(result.dispersion_energy_over_k, source.dispersion_energy_over_k);
    permute(result.molar_mass_kg_per_mol, source.molar_mass_kg_per_mol);
    permute(result.charge_number, source.charge_number);
    permute(result.packing_diameter_angstrom, source.packing_diameter_angstrom);
    permute(result.debye_huckel_diameter_angstrom, source.debye_huckel_diameter_angstrom);
    permute(result.born_diameter_angstrom, source.born_diameter_angstrom);
    permute(result.solvation_factor, source.solvation_factor);
    permute(result.relative_permittivity, source.relative_permittivity);
    permute(result.permittivity_reference_temperature_k, source.permittivity_reference_temperature_k);
    permute(result.ion_suppression_coefficient, source.ion_suppression_coefficient);
    permute(result.dielectric_volume_m3_per_mol, source.dielectric_volume_m3_per_mol);
    permute(result.dielectric_temperature_k, source.dielectric_temperature_k);
    permute(result.dipole_scaling_debye2, source.dipole_scaling_debye2);
    permute(result.polarizability_scaling_angstrom3, source.polarizability_scaling_angstrom3);
    permute(result.correlation_integral_parameter, source.correlation_integral_parameter);
    permute_matrix(result.k_ij, source.k_ij);
    permute_matrix(result.k_ij_reciprocal_temperature_slope, source.k_ij_reciprocal_temperature_slope);
    permute_matrix(result.k_ij_reference_temperature_k, source.k_ij_reference_temperature_k);
    permute_matrix(result.l_ij, source.l_ij);
    std::vector<int> inverse(n);
    for (std::size_t i = 0; i < n; ++i) inverse[order[i]] = static_cast<int>(i);
    for (auto& component : result.sigma_term_component) component = inverse[component];
    for (auto& component : result.site_component) component = inverse[component];
    for (auto& component : result.permittivity_term_component) component = inverse[component];
    return result;
}

template <class Vector>
Vector permute_vector(const Vector& source, const std::vector<std::size_t>& order) {
    Vector result = source;
    for (std::size_t i = 0; i < order.size(); ++i) result[i] = source[order[i]];
    return result;
}

}  // namespace equilibrium_support

#define CHECK(condition) \
    equilibrium_support::check(static_cast<bool>(condition), __FILE__, __LINE__, #condition)
