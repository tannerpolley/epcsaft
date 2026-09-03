#include "equilibrium/gibbs_minimizer.hpp"

#include "core/trace.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

namespace epcsaft::equilibrium {
namespace {

std::vector<std::size_t> offsets(const CompiledSystem& compiled) {
    std::vector<std::size_t> result;
    std::size_t offset = 0;
    for (const auto& phase : compiled.problem.phases) {
        result.push_back(offset);
        offset += phase.support.size();
        if (phase.amount.kind == AmountKind::Free) ++offset;
        ++offset;
    }
    return result;
}

void seed_composition(const CompiledSystem& compiled, std::size_t phase,
                      const std::vector<double>& composition,
                      std::vector<double>& point, std::size_t offset) {
    const auto& spec = compiled.problem.phases[phase];
    double total = 0.0;
    for (std::size_t i = 0; i < spec.support.size(); ++i) {
        const std::size_t species = spec.support[i];
        double value = composition[species];
        if (!(value > 0.0)) value = 1.0 / static_cast<double>(spec.support.size());
        point[offset + i] = value;
        total += value;
    }
    for (std::size_t i = 0; i < spec.support.size(); ++i)
        point[offset + i] = std::log(point[offset + i] / total);
}

double seed_temperature(const CompiledSystem& compiled) {
    return compiled.problem.temperature.is_positive()
        ? compiled.problem.temperature.value : 300.0;
}

double seed_pressure(const CompiledSystem& compiled) {
    return compiled.problem.pressure.is_positive()
        ? compiled.problem.pressure.value : 1.0e5;
}

void validate_density_seed(const DensityRoot& root) {
    if (!(std::isfinite(root.molar_density) && root.molar_density > 0.0)) {
        throw std::runtime_error("boundary seed has invalid density");
    }
}

std::vector<std::vector<double>> start_compositions(
    const CompiledSystem& compiled, core::Scope& scope
) {
    const std::size_t phases = compiled.problem.phases.size();
    const std::size_t components = compiled.model->parameters().component_count();
    std::vector<std::vector<double>> result(phases, compiled.feed_composition);
    bool guessed = true;
    for (std::size_t phase = 0; phase < phases; ++phase) {
        const auto& spec = compiled.problem.phases[phase];
        guessed = guessed && spec.composition_guess.size() == spec.support.size();
        if (spec.composition_guess.size() == spec.support.size()) {
            result[phase].assign(components, 0.0);
            for (std::size_t i = 0; i < spec.support.size(); ++i)
                result[phase][spec.support[i]] = spec.composition_guess[i];
        }
    }
    if (phases == 1) return result;
    if (guessed) {
        scope.probe("start_rule", 1.0);
        return result;
    }
    if (phases != 2 || compiled.problem.phases[0].kind == compiled.problem.phases[1].kind) {
        scope.probe("start_rule", 3.0);
        throw std::runtime_error("start_symmetric: phases "
            + compiled.problem.phases[0].name + " and "
            + compiled.problem.phases[1].name + " require composition guesses");
    }
    const double temperature = seed_temperature(compiled);
    const double pressure = seed_pressure(compiled);
    const auto vapor = compiled.model->state_at_pressure(
        temperature, pressure, compiled.feed_composition, DensityBranch::Vapor);
    const auto liquid = compiled.model->state_at_pressure(
        temperature, pressure, compiled.feed_composition, DensityBranch::Liquid);
    if (!vapor.log_fugacity_coefficient || !liquid.log_fugacity_coefficient)
        throw std::runtime_error("boundary start has no fugacity coefficients");
    const auto& ln_phi_vapor = *vapor.log_fugacity_coefficient;
    const auto& ln_phi_liquid = *liquid.log_fugacity_coefficient;
    std::vector<double> x(components), y(components);
    double x_total = 0.0, y_total = 0.0;
    for (std::size_t i = 0; i < components; ++i) {
        const double k = std::exp(ln_phi_liquid[i] - ln_phi_vapor[i]);
        scope.probe("k_factor", k, static_cast<std::uint32_t>(i));
        x[i] = compiled.feed_composition[i] / k;
        y[i] = compiled.feed_composition[i] * k;
        x_total += x[i]; y_total += y[i];
    }
    for (std::size_t i = 0; i < components; ++i) { x[i] /= x_total; y[i] /= y_total; }
    for (std::size_t phase = 0; phase < phases; ++phase)
        result[phase] = compiled.problem.phases[phase].kind == PhaseKind::Vapor ? y : x;
    scope.probe("start_rule", 2.0);
    return result;
}

}  // namespace

std::vector<double> boundary_start(const CompiledSystem& compiled) {
    core::Scope scope("boundary_start");
    std::vector<double> point = compiled.generic_point();
    const auto phase_offsets = offsets(compiled);
    const auto compositions = start_compositions(compiled, scope);
    std::size_t free_phases = 0;
    for (const auto& phase : compiled.problem.phases) {
        if (phase.amount.kind == AmountKind::Free) ++free_phases;
    }
    double feed_total = 0.0;
    for (const double value : compiled.n_feed) feed_total += value;
    const double free_amount = compiled.n_feed.empty() ? 1.0
        : feed_total / static_cast<double>(std::max<std::size_t>(1, free_phases));
    scope.probe("free_phases", static_cast<double>(free_phases));
    scope.probe("feed_total", feed_total);
    for (std::size_t phase = 0; phase < compiled.problem.phases.size(); ++phase) {
        seed_composition(compiled, phase, compositions[phase], point, phase_offsets[phase]);
        std::size_t cursor = phase_offsets[phase] + compiled.problem.phases[phase].support.size();
        if (compiled.problem.phases[phase].amount.kind == AmountKind::Free) point[cursor++] = free_amount;
        const double t = seed_temperature(compiled);
        const double p = seed_pressure(compiled);
        std::vector<double> composition(compiled.model->parameters().component_count(), 0.0);
        for (std::size_t i = 0; i < compiled.problem.phases[phase].support.size(); ++i) {
            composition[compiled.problem.phases[phase].support[i]] =
                std::exp(point[phase_offsets[phase] + i]);
        }
        const DensityBranch branch = compiled.problem.phases[phase].kind == PhaseKind::Vapor
            ? DensityBranch::Vapor : DensityBranch::Liquid;
        DensityRoot root;
        try {
            root = solve_density(*compiled.model, t, p, composition, branch, std::nullopt);
        } catch (const std::runtime_error& error) {
            throw std::runtime_error(
                std::string("boundary seed has no density root: ") + error.what()
            );
        }
        if (compiled.model->parameters().component_count() == 1
            && (!root.stable_root_count.has_value() || *root.stable_root_count < 2)) {
            throw std::runtime_error("boundary seed has no coexisting density roots");
        }
        validate_density_seed(root);
        scope.probe("phase", static_cast<double>(phase));
        scope.probe("branch_density", root.molar_density,
                    static_cast<std::uint32_t>(phase));
        scope.probe("pressure", p);
        scope.probe("temperature", t);
        point[cursor] = std::log(compiled.model->packing_fraction_coefficient(t, composition)
                                 * root.molar_density);
    }
    return point;
}

// Boundary seeds remain ordinary equilibrium coordinates; no second driver is introduced.

}  // namespace epcsaft::equilibrium
