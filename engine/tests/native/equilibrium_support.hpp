#pragma once

#include "core/trace.hpp"
#include "equilibrium/gibbs_minimizer.hpp"
#include "equilibrium/problem.hpp"
#include "properties_support.hpp"
#include "recording_support.hpp"

#include <cstddef>
#include <algorithm>
#include <cstdio>
#include <exception>
#include <string>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace equilibrium_support {

// Record labels of declaration enums.
inline const char* phase_kind_name(epcsaft::equilibrium::PhaseKind kind) noexcept {
    return kind == epcsaft::equilibrium::PhaseKind::Liquid ? "liquid" : "vapor";
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

// The neutral-binary bubble problem at T: L free, V pinned 0, equimolar feed, P free from 1e5 Pa.
inline epcsaft::equilibrium::EquilibriumProblem bubble_at(double temperature) {
    auto problem = boundary();
    problem.feed.values = {0.5, 0.5};
    problem.temperature = epcsaft::equilibrium::FreeValue::fixed(temperature);
    return problem;
}

// The E4-B interior tie line: neutral binary at 180 K, feed at the mean of the bubble
// compositions, both amounts free at half the bubble pressure; the recorded coordinate seed
// is used only if the cold solve of the same problem is not accepted.
struct Interior final {
    epcsaft::equilibrium::EquilibriumProblem problem;
    epcsaft::equilibrium::EquilibriumResult start;
    double bubble_pressure = 0.0;
    std::vector<double> feed;
};

inline Interior interior(const epcsaft::HelmholtzModel& model) {
    Interior result;
    const auto bubble = epcsaft::equilibrium::solve_equilibrium(model, bubble_at(180.0));
    recording_support::check(bubble.success, __FILE__, __LINE__, "bubble.success");
    result.bubble_pressure = *bubble.pressure;
    result.feed = {0.5 * (bubble.mole_fractions[0] + bubble.mole_fractions[2]),
                   0.5 * (bubble.mole_fractions[1] + bubble.mole_fractions[3])};
    result.problem = bubble_at(180.0);
    result.problem.feed.values = result.feed;
    result.problem.phases[1].amount = epcsaft::equilibrium::AmountSpec::free();
    result.problem.pressure = epcsaft::equilibrium::FreeValue::fixed(0.5 * result.bubble_pressure);
    result.start = epcsaft::equilibrium::solve_equilibrium(model, result.problem);
    if (!(result.start.success && result.start.requested_tolerance_met.value_or(false))) {
        const std::vector<double> seed = {-1.4805904749125276, -0.25812755981945695, 0.25979825064587636,
            -0.9468581761373848, -0.10102638484899645, -2.342461524948757, 0.7402017493541235,
            -4.6620373689869785, 6.106829374848107, 3.7251308960023377};
        result.start = epcsaft::equilibrium::solve_equilibrium_from_coordinates(model, result.problem, seed);
    }
    recording_support::check(result.start.success, __FILE__, __LINE__, "result.start.success");
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

// The declared state an equilibrium row is evaluated at.
inline std::string state_inputs(const epcsaft::equilibrium::CompiledSystem& compiled) {
    const auto& problem = compiled.problem;
    return recording_support::fields({recording_support::field("temperature", problem.temperature.value),
        recording_support::field("pressure", problem.pressure.value), recording_support::field("feed", problem.feed.values)});
}

}  // namespace equilibrium_support

#define CHECK(condition) \
    recording_support::check(static_cast<bool>(condition), __FILE__, __LINE__, #condition)
