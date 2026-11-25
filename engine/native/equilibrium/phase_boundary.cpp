#include "equilibrium/gibbs_minimizer.hpp"

#include "core/bitwise_equal.hpp"
#include "core/numerical_policy.hpp"
#include "core/trace.hpp"
#include "equilibrium/stoichiometry.hpp"

#include <Eigen/Dense>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <functional>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>

namespace epcsaft {
struct PressureDensityDerivative final {
    static double at(const HelmholtzModel& model, const equilibrium::PhasePoint& point, std::size_t phase,
        double density, const std::vector<double>& composition) {
        const double volume = 1.0 / density;
        const auto* completed = point.state_model == model.identity()
            && core::bitwise_equal(volume, point.volume[phase])
            && core::bitwise_equal(composition, point.composition[phase])
            ? &point.derivative[phase] : nullptr;
        return model.state(point.temperature, density, composition, false, false, false, completed).pressure_density_derivative;
    }
};
}
namespace epcsaft::equilibrium {

bool supports_split_seed(const CompiledSystem& compiled) {
    const auto& p = compiled.problem;
    const auto components = compiled.model->parameters().component_count();
    return components > 1 && p.phases.size() == 2
        && compiled.model->parameters().site_component.empty()
        && p.phases[0].kind != p.phases[1].kind
        && !p.temperature.free && !p.pressure.free
        && !compiled.has_ions && !compiled.has_reaction && !compiled.feasibility_cell()
        && std::all_of(compiled.feed_composition.begin(), compiled.feed_composition.end(),
            [](double z) { return std::isfinite(z) && z > 0.0; })
        && p.phases[0].support.size() == components && p.phases[1].support.size() == components;
}

namespace {

// The one rethrow of a continuation corrector's numerical failure, after keeping its observations.
EquilibriumResult unless_failed(EquilibriumResult solved, StageWork* observations) {
    if (!solved.failure) return solved;
    if (observations) observations->observe(solved);
    std::rethrow_exception(solved.failure);
}

std::vector<double> seed_composition(const CompiledSystem& compiled, std::size_t phase,
                      std::vector<double> composition,
                      std::vector<double>& point, std::size_t offset, bool declared = false) {
    const auto& spec = compiled.problem.phases[phase];
    const bool guessed = declared && spec.declared_composition_guess && !compiled.problem.prescribed_liquid_composition;
    double total = 0.0;
    for (std::size_t i = 0; i < spec.support.size(); ++i) {
        const std::size_t species = spec.support[i];
        double value = composition[species];
        if (!(value > 0.0)) value = 1.0 / static_cast<double>(spec.support.size());
        // ALG-EQUILIBRIUM trace floor: a guessed entry on the composition box feels only the barrier.
        else if (guessed) value = std::max(value, policy::EQUILIBRIUM_TRACE_MOIETY_FRACTION.value);
        point[offset + i] = value;
        total += value;
    }
    if (guessed) {
        const std::size_t components = compiled.feed_composition.size();
        bool identity = compiled.B.size() == components;
        for (std::size_t i = 0; i < compiled.B.size(); ++i)
            for (std::size_t j = 0; j < components; ++j)
                identity = identity && compiled.B[i][j] == static_cast<double>(i == j);
        std::vector<double> upper(spec.support.size());
        double removed = 0.0, recipients = 0.0;
        for (std::size_t i = 0; i < spec.support.size(); ++i) {
            const auto species = spec.support[i];
            double u = identity ? compiled.feed_composition[species] : compiled.inventory_upper_bounds[species];
            if (!(u > 0.0 && u < policy::EQUILIBRIUM_TRACE_MOIETY_FRACTION.value))
                u = std::numeric_limits<double>::infinity();
            upper[i] = u;
            const double g = point[offset + i] / total;
            removed += g - std::min(g, u);
            if (!std::isfinite(u)) recipients += g;
        }
        if (removed > 0.0 && recipients > 0.0) {
            for (std::size_t i = 0; i < spec.support.size(); ++i) {
                const double g = point[offset + i] / total;
                const double value = std::isfinite(upper[i]) ? std::min(g, upper[i]) : g + removed * g / recipients;
                composition[spec.support[i]] = value;
                point[offset + i] = std::log(value);
            }
            return composition;
        }
    }
    for (std::size_t i = 0; i < spec.support.size(); ++i)
        point[offset + i] = std::log(point[offset + i] / total);
    return composition;
}

// A start's compositions in declared phase order and, after a split, its free amounts; vapor_density anchors the
// true-species split's vapor root; a positive pressure replaces a pure fluid's free-pressure seed (NaN: the loop rule found none).
struct PhaseStart final {
    std::vector<std::vector<double>> compositions;
    std::vector<double> amounts{};
    double vapor_density = 0.0, pressure = 0.0, loop_pressure = 0.0;
    std::array<std::optional<State>, 2> feed_states{};
    std::vector<double> densities{};
    bool declared_guesses = false;
};

// The Rachford-Rice split of `total` mol of composition z: beta is the interior root of
// sum_i z_i (K_i - 1)/(1 - beta + beta K_i) by bisection; x_i = z_i/(1 - beta + beta K_i), y_i = K_i x_i,
// N_L = total (1 - beta) and N_V = total beta. K_i = 0 marks a liquid-confined species; its positive inventory makes
// the right end -infinity. Empty without a bracket, below PHASE_DISAPPEARANCE_AMOUNT of the feed or at a
// nonpositive fraction.
std::optional<PhaseStart> rachford_rice(const CompiledSystem& compiled, const std::vector<double>& z,
                                        const std::vector<double>& k, double total, core::Scope& scope) {
    const auto balance = [&](double beta) {
        double value = 0.0;
        for (std::size_t i = 0; i < z.size(); ++i)
            if (z[i] > 0.0) value += z[i] * (k[i] - 1.0) / ((1.0 - beta) + beta * k[i]);
        return value;
    };
    const double left = balance(0.0), right = balance(1.0);
    if (!(std::isfinite(left) && left > 0.0 && right < 0.0)) return std::nullopt;
    double lower = 0.0, upper = 1.0;
    for (int step = 0; step < std::numeric_limits<double>::digits; ++step) {
        const double beta = 0.5 * (lower + upper);
        if (balance(beta) > 0.0) lower = beta;
        else upper = beta;
    }
    const double beta = 0.5 * (lower + upper), share = total / compiled.feed_total;
    if (!(share * beta > policy::PHASE_DISAPPEARANCE_AMOUNT.value
          && share * (1.0 - beta) > policy::PHASE_DISAPPEARANCE_AMOUNT.value)) return std::nullopt;
    std::vector<double> x(z.size(), 0.0), y(z.size(), 0.0);
    for (std::size_t i = 0; i < z.size(); ++i) {
        if (!(z[i] > 0.0)) continue;
        x[i] = z[i] / ((1.0 - beta) + beta * k[i]);
        y[i] = k[i] * x[i];
        if (!(std::isfinite(x[i]) && x[i] > 0.0 && std::isfinite(y[i]) && (k[i] == 0.0 || y[i] > 0.0))) return std::nullopt;
    }
    PhaseStart start;
    for (const auto& phase : compiled.problem.phases) {
        const bool vapor = phase.kind == PhaseKind::Vapor;
        start.compositions.push_back(vapor ? y : x);
        start.amounts.push_back(total * (vapor ? beta : 1.0 - beta));
    }
    scope.probe("vapor_fraction_seed", beta);
    return start;
}

// The feed root resolves its liquid side: two stable roots, or a sole root denser than its negative-slope density.
template <class Root> bool liquid_feed(const Root& r) {
    return r.stable_root_count.value_or(0) > 1 || (r.negative_slope_density && r.molar_density > *r.negative_slope_density);
}

bool pressure_loop_eligible(const CompiledSystem& c) {
    return c.problem.pressure.free && !c.problem.temperature.free && newton_eligible(c) && recovery_eligible(c)
        && (c.model->parameters().component_count() == 1 || !c.problem.pressure.is_positive());
}

PhaseStart start_compositions(const CompiledSystem& compiled, core::Scope& scope) {
    const std::size_t phases = compiled.problem.phases.size();
    const std::size_t components = compiled.model->parameters().component_count();
    std::vector<std::vector<double>> result(phases, compiled.feed_composition);
    if (pressure_loop_eligible(compiled)) if (auto rho = superancillary_densities(*compiled.model, compiled.seed_temperature())) {
        const double pressure = compiled.model->pure_isothermal(compiled.seed_temperature(), 1.0 / (*rho)[1]).pressure;
        scope.probe("superancillary_start", pressure);
        if (compiled.problem.phases[0].kind == PhaseKind::Vapor) std::swap((*rho)[0], (*rho)[1]);
        return {result, {}, 0.0, pressure, 0.0, {}, {(*rho)[0], (*rho)[1]}};
    }
    std::string unguessed;
    for (std::size_t phase = 0; phase < phases; ++phase) {
        const auto& spec = compiled.problem.phases[phase];
        if (spec.composition_guess.size() != spec.support.size()) {
            unguessed += (unguessed.empty() ? "" : ", ") + spec.name;
            continue;
        }
        result[phase].assign(components, 0.0);
        for (std::size_t i = 0; i < spec.support.size(); ++i)
            result[phase][spec.support[i]] = spec.composition_guess[i];
    }
    if (phases == 1) return {result, {}, 0.0, 0.0, 0.0, {}, {}, unguessed.empty()};
    if (unguessed.empty()) {
        scope.probe("start_rule", 1.0);
        return {result, {}, 0.0, 0.0, 0.0, {}, {}, true};
    }
    if (phases != 2 || compiled.problem.phases[0].kind == compiled.problem.phases[1].kind) {
        scope.probe("start_rule", 3.0);
        throw core::DeclarationError("phases without a composition_guess: " + unguessed
            + "; a start without guesses covers one liquid and one vapor only;"
              " give each phase a composition_guess or use phases=Discover()");
    }
    const double temperature = compiled.seed_temperature();
    const double pressure = compiled.seed_pressure();
    DensityRoot liquid_root;
    const auto vapor_root = solve_density(
        *compiled.model, temperature, pressure, compiled.feed_composition,
        DensityBranch::Vapor, std::nullopt, &liquid_root);
    auto vapor = compiled.model->state(
        temperature, vapor_root.molar_density, compiled.feed_composition);
    auto liquid = compiled.model->state(
        temperature, liquid_root.molar_density, compiled.feed_composition);
    if (!vapor.log_fugacity_coefficient || !liquid.log_fugacity_coefficient)
        throw std::runtime_error("boundary start has no fugacity coefficients");
    // A sole unlabelled feed root cannot supply both phases' fugacity estimates:
    // their self-ratio is one. Use the ideal-vapor limit for this start only.
    const auto ln_phi_vapor = vapor_root.stable_root_count == 1
        ? std::vector<std::optional<double>>(components, 0.0) : *vapor.log_fugacity_coefficient;
    const auto& ln_phi_liquid = *liquid.log_fugacity_coefficient;
    // Pure fixed-T, free-P saturation with a sole stable root at P0: P at the scan's negative-slope density lies
    // between the single mechanical loop's spinodal pressures, so it seeds the pressure and both densities.
    double witness = 0.0;
    if (components == 1 && compiled.problem.pressure.free && !compiled.problem.temperature.free
        && vapor_root.stable_root_count == 1) {
        witness = vapor_root.negative_slope_density ? compiled.model->state(
            temperature, *vapor_root.negative_slope_density, compiled.feed_composition).pressure : NAN;
        if (std::isfinite(witness) && witness > 0.0) { scope.probe("loop_pressure_seed", witness); return {result, {}, 0.0, witness}; }
        witness = NAN;
    }
    const double loop_pressure = pressure_loop_eligible(compiled) && !liquid_feed(liquid_root)
        && vapor_root.negative_slope_density ? compiled.model->state(
            temperature, *vapor_root.negative_slope_density, compiled.feed_composition).pressure : 0.0;
    liquid.stable_root_count = liquid_root.stable_root_count; liquid.negative_slope_density = liquid_root.negative_slope_density;
    vapor.stable_root_count = vapor_root.stable_root_count; vapor.negative_slope_density = vapor_root.negative_slope_density;
    std::vector<double> x(components), y(components), k_values(components);
    double x_total = 0.0, y_total = 0.0;
    for (std::size_t i = 0; i < components; ++i) {
        if (compiled.feed_composition[i] == 0.0) continue;
        if (std::find(compiled.removed_species.begin(), compiled.removed_species.end(), i)
            != compiled.removed_species.end()) continue;
        const double k = std::exp(require_component_property(ln_phi_liquid[i]) - require_component_property(ln_phi_vapor[i]));
        k_values[i] = k;
        scope.probe("k_factor", k, static_cast<std::uint32_t>(i));
        x[i] = compiled.feed_composition[i] / k;
        y[i] = compiled.feed_composition[i] * k;
        x_total += x[i]; y_total += y[i];
    }
    for (std::size_t i = 0; i < components; ++i) { x[i] /= x_total; y[i] /= y_total; }
    scope.probe("start_rule", 2.0);
    if (supports_split_seed(compiled)
        && std::all_of(k_values.begin(), k_values.end(), [](double k) { return std::isfinite(k) && k > 0.0; }))
        if (auto split = rachford_rice(compiled, compiled.feed_composition, k_values, compiled.feed_total, scope))
            return std::move(*split);
    // Use the conserved finite liquid only when the feed root resolves its liquid side.
    for (std::size_t phase = 0; phase < phases; ++phase) {
        const auto& own = compiled.problem.phases[phase].amount, & other = compiled.problem.phases[1 - phase].amount;
        result[phase] = !compiled.has_reaction && own.kind == AmountKind::Free && other.kind == AmountKind::Pinned
            && other.fraction == 0.0 && ((compiled.problem.phases[phase].kind == PhaseKind::Liquid && liquid_feed(liquid_root)) || [&] {
                const auto& support = compiled.problem.phases[1 - phase].support;
                for (std::size_t i = 0; i < components; ++i)
                    if (compiled.feed_composition[i] > 0.0 && std::find(support.begin(), support.end(), i) == support.end()) return true;
                return false;
            }()) ? compiled.feed_composition : compiled.problem.phases[phase].kind == PhaseKind::Vapor ? y : x;
    }
    return {result, {}, 0.0, witness, loop_pressure, {std::move(liquid), std::move(vapor)}};
}

// ALG-EQUILIBRIUM: a converged ln P start, otherwise the existing one-step bubble seed.
double fugacity_pressure_seed(const CompiledSystem& compiled, std::vector<double>& point,
                            const std::vector<std::size_t>& phase_offsets, const PhaseStart& start, core::Scope& scope, bool& corrected, bool correct_dense, const DriverLimits* limits, StageWork* stages, bool allow_loop = true) {
    const auto& problem = compiled.problem;
    const double seed = compiled.seed_pressure(), fallback = start.pressure > 0.0 ? start.pressure : seed;
    const bool loop = allow_loop && pressure_loop_eligible(compiled);
    if (!problem.pressure.free || problem.temperature.free || problem.phases.size() != 2
        || problem.phases[0].kind == problem.phases[1].kind
        || (!loop && compiled.model->parameters().component_count() < 2)) return fallback;
    const std::size_t vapor = problem.phases[0].kind == PhaseKind::Vapor ? 0 : 1;
    const bool zero_v = problem.phases[vapor].is_zero_pinned();
    if (!loop && !zero_v) return fallback;
    const std::size_t pinned = zero_v ? vapor : 1 - vapor;
    const auto& support = problem.phases[pinned].support;
    const std::size_t components = compiled.model->parameters().component_count();
    std::vector<std::vector<double>> x(2, std::vector<double>(components, 0.0));
    for (std::size_t phase = 0; phase < 2; ++phase) {
        const auto& own = problem.phases[phase];
        for (std::size_t i = 0; i < own.support.size(); ++i)
            x[phase][own.support[i]] = loop ? (phase != pinned || own.composition_guess.empty()
                ? compiled.feed_composition[own.support[i]] : own.composition_guess[i]) : std::exp(point[phase_offsets[phase] + i]);
    }
    const double temperature = compiled.seed_temperature();
    double pressure = loop && std::isfinite(start.loop_pressure) && start.loop_pressure > 0.0 ? start.loop_pressure : fallback;
    try {
        for (int pass = 0; pass < (loop ? policy::EQUILIBRIUM_START_PRESSURE_PASSES.value : 1.0); ++pass) {
            std::array<State, 2> states;
            for (std::size_t side = 0; side < 2; ++side) {
                const auto& composition = x[side == 1 ? vapor : 1 - vapor];
                states[side] = loop && pass == 0 && core::bitwise_equal(pressure, seed) && core::bitwise_equal(composition, compiled.feed_composition) && start.feed_states[side]
                    ? *start.feed_states[side] : compiled.model->state_at_pressure(temperature, pressure, composition,
                        side == 0 ? DensityBranch::Liquid : DensityBranch::Vapor);
            }
            if (correct_dense && states[1].stable_root_count == 1 && liquid_feed(states[1])) {
                corrected = true;
                scope.probe("dense_vapor_start", pressure);
                scope.probe("dense_vapor_roots", static_cast<double>(states[1].stable_root_count.value_or(0)));
                scope.probe("dense_vapor_witness", states[1].negative_slope_density.value_or(NAN));
                for (int reduction = 0; reduction < policy::EQUILIBRIUM_START_PRESSURE_PASSES.value; ++reduction) {
                    if (limits) limits->check();
                    pressure *= 0.5;
                    scope.probe("dense_vapor_pressure", pressure);
                    states[1] = compiled.model->state_at_pressure(temperature, pressure, x[vapor], DensityBranch::Vapor);
                    scope.probe("dense_vapor_density", states[1].molar_density);
                    if (states[1].returned_density_branch() == DensityBranch::Vapor) break;
                    if (!(states[1].stable_root_count == 1 && liquid_feed(states[1])))
                        throw std::runtime_error("dense vapor start: reduced-pressure root side is unresolved");
                }
                if (states[1].returned_density_branch() != DensityBranch::Vapor)
                    throw std::runtime_error("dense vapor start: eight pressure reductions did not establish vapor");
                states[0] = compiled.model->state_at_pressure(temperature, pressure, x[1-vapor], DensityBranch::Liquid);
                if (states[0].returned_density_branch() == DensityBranch::Vapor)
                    throw std::runtime_error("dense vapor start: reduced-pressure liquid lies on vapor side");
            }
            if (!states[0].log_fugacity_coefficient || !states[1].log_fugacity_coefficient) {
                if (corrected) throw std::runtime_error("dense vapor start: fugacity unavailable");
                break;
            }
            if (loop && states[0].stable_root_count == 1 && !liquid_feed(states[0])) break;
            const bool ideal = loop && states[1].stable_root_count == 1 && liquid_feed(states[1]);
            const double dz = (ideal ? 1.0 : states[1].compressibility_factor) - states[0].compressibility_factor;
            if (loop && !(dz > 0.0)) break;
            std::vector<double> next(support.size());
            double sum = 0.0, change = 0.0;
            for (std::size_t k = 0; k < support.size(); ++k) {
                const std::size_t i = support[k];
                const double log_k = require_component_property((*states[0].log_fugacity_coefficient)[i])
                    - (ideal ? 0.0 : require_component_property((*states[1].log_fugacity_coefficient)[i]));
                next[k] = x[1 - pinned][i] * std::exp(zero_v ? log_k : -log_k);
                sum += next[k];
            }
            const double step = (zero_v ? std::log(sum) : -std::log(sum)) / dz;
            const double stepped = loop ? pressure * std::exp(step) : (corrected ? pressure : seed) * sum;
            if (!(std::isfinite(stepped) && stepped > 0.0)) {
                if (corrected) throw std::runtime_error("dense vapor start: fugacity pressure is invalid");
                break;
            }
            for (std::size_t k = 0; k < support.size(); ++k) {
                next[k] /= sum; change = std::max(change, std::abs(next[k] - x[pinned][support[k]]));
                x[pinned][support[k]] = next[k];
            }
            if (!loop || (std::abs(step) < policy::EQUILIBRIUM_START_PRESSURE_STEP.value
                && change < policy::EQUILIBRIUM_START_PRESSURE_STEP.value)) {
                if (ideal) break;
                if (loop || std::all_of(problem.phases.begin(), problem.phases.end(),
                        [](const PhaseSpec& phase) { return phase.composition_guess.size() == phase.support.size(); }))
                    for (std::size_t k = 0; k < support.size(); ++k)
                        if (next[k] > 0.0) point[phase_offsets[pinned] + k] = std::log(next[k]);
                if (loop) seed_composition(compiled, 1 - pinned, compiled.feed_composition, point, phase_offsets[1 - pinned]);
                scope.probe("fugacity_pressure_seed", stepped); return stepped;
            }
            pressure = stepped;
        }
    } catch (const std::exception& error) {
        if (corrected) throw;
        if (stages) retain_density_failure(stages->first_numerical_failure, error, "pressure seed");
    }
    return loop ? fugacity_pressure_seed(compiled, point, phase_offsets, start, scope, corrected, false, limits, stages, false) : fallback;
}

// Ideal-solution chemical equilibrium with frozen c_j: minimize sum_j n_j (c_j + ln n_j/n) subject to B n = b
// by the Gordon-McBride element-potential Newton iteration and step control (NASA RP-1311, 1994, eqs
// 2.24-2.26, 3.1-3.3), the box separating major from trace species; the trace limit applies only to a
// species growing faster than the total.
bool ideal_speciation(const std::vector<std::vector<double>>& B, const std::vector<double>& b,
                      const std::vector<double>& c, std::vector<double>& n, StageWork* observations) {
    const std::size_t r = B.size(), m = c.size();
    const double trace = -std::log(policy::EQUILIBRIUM_COMPOSITION_BOX.value);
    std::vector<double> ln_n(m), g(m);
    for (std::size_t j = 0; j < m; ++j) ln_n[j] = std::log(n[j]);
    double ln_total = std::log(std::accumulate(n.begin(), n.end(), 0.0));
    for (int iteration = 0; iteration < 100; ++iteration) {
        Eigen::MatrixXd A = Eigen::MatrixXd::Zero(r + 1, r + 1);
        Eigen::VectorXd rhs = Eigen::VectorXd::Zero(r + 1);
        rhs(r) = std::exp(ln_total);
        A(r, r) = -rhs(r);
        for (std::size_t k = 0; k < r; ++k) rhs(k) = b[k];
        for (std::size_t j = 0; j < m; ++j) {
            n[j] = std::exp(ln_n[j]);
            g[j] = c[j] + ln_n[j] - ln_total;
            for (std::size_t k = 0; k <= r; ++k) {
                const double a = k < r ? B[k][j] : 1.0;
                for (std::size_t i = 0; i < r; ++i) A(k, i) += a * B[i][j] * n[j];
                A(k, r) += a * n[j];
                rhs(k) += a * n[j] * (g[j] - 1.0);
            }
        }
        const auto lu = A.fullPivLu();
        if (observations) core::retain_rank(observations->rank_records, core::numerical_rank(lu, "ideal_speciation"));
        const Eigen::VectorXd x = lu.solve(rhs);
        double change = std::abs(x(r)), major = 5.0 * change, control = 1.0;
        for (std::size_t j = 0; j < m; ++j) {
            double step = x(r) - g[j];
            for (std::size_t i = 0; i < r; ++i) step += B[i][j] * x(i);
            change = std::max(change, std::abs(step));
            if (ln_n[j] - ln_total > -trace) major = std::max(major, std::abs(step));
            else if (step > std::max(0.0, x(r)))
                control = std::min(control, (ln_total - ln_n[j] - 0.5 * trace) / (step - x(r)));
            g[j] = step;
        }
        if (!std::isfinite(change)) { if (observations) core::retain_failure(observations->first_numerical_failure, "nonfinite ideal speciation change", "seed"); return false; }
        if (change <= policy::EQUILIBRIUM_START_SPECIATION_STEP.value) return true;
        control = std::max(0.0, std::min(control, 2.0 / std::max(major, 2.0)));
        for (std::size_t j = 0; j < m; ++j) ln_n[j] += control * g[j];
        ln_total += control * x(r);
    }
    if (observations) core::retain_failure(observations->first_numerical_failure, "ideal speciation iteration cap", "seed");
    return false;
}

// A reacting declaration of one free phase starts from the ideal-solution speciation of the feed balances,
// with the residual, volume and standard terms of its chemical rows frozen at the seeded composition and
// refreshed over sweeps (#154, #157): an unspeciated start leaves the log chart's trace
// species without restoring force and the optimizer cycles. Zero guesses enter at the optimizer box; a
// failed sweep keeps the last speciation.
// Returns the speciated amount and the last sweep's density root, which seeds the start's packing.
// `index` selects the speciated phase; its chemical rows alone supply c_j, so other phases' coordinates stay as given.
std::optional<std::pair<double, double>> speciate_single_phase(const CompiledSystem& compiled, std::vector<double>& point,
                                                               std::vector<double> composition, core::Scope& scope,
                                                               std::size_t index = 0, StageWork* observations = nullptr, bool selected_chemical_phase = false) {
    const auto& phase = compiled.problem.phases[index];
    const std::size_t m = phase.support.size(), offset = compiled.phase_offsets()[index];
    std::vector<std::size_t> chemical;  // phase `index`'s chemical rows, in support order
    for (std::size_t row = 0; row < compiled.rows.size(); ++row)
        if (compiled.rows[row].kind == RowKind::Chemical && compiled.rows[row].phase == index) chemical.push_back(row);
    const double t = compiled.seed_temperature(), p = compiled.seed_pressure();
    std::vector<std::vector<double>> B(compiled.rank, std::vector<double>(m));
    for (std::size_t k = 0; k < compiled.rank; ++k)
        for (std::size_t j = 0; j < m; ++j) B[k][j] = compiled.B[k][phase.support[j]];
    for (const std::size_t i : phase.support) composition[i] = std::max(composition[i], policy::EQUILIBRIUM_COMPOSITION_BOX.value);
    std::vector<double> n(m), speciated, trial = point;
    for (std::size_t j = 0; j < m; ++j) n[j] = compiled.feed_total * composition[phase.support[j]];
    std::optional<double> density, speciated_density;
    try {
        const auto root = [&](std::optional<double> anchor) {
            const auto root = solve_density(*compiled.model, t, p, composition,
                branch_of(phase.kind), anchor);
            if (observations) core::retain_first(observations->first_numerical_failure, root.first_numerical_failure);
            return root.molar_density;
        };
        for (int sweep = 0; sweep < 2; ++sweep) {
            const double total = std::accumulate(composition.begin(), composition.end(), 0.0);
            for (double& value : composition) value /= total;
            try {  // continue the previous sweep's root; a speciation jump beyond its window solves it anew
                density = root(density);
            } catch (const std::exception& error) {
                if (observations) retain_density_failure(observations->first_numerical_failure, error, "speciation density");
                if (!density) throw;
                density = root(std::nullopt);
            }
            for (std::size_t j = 0; j < m; ++j) trial[offset + j] = std::log(composition[phase.support[j]]);
            trial[offset + m] = compiled.feed_total;
            trial[offset + m + 1] = std::log((compiled.model->parameters().ideal_gas ? 1.0 : compiled.model->packing_fraction_coefficient(t, composition)) * *density);
            std::fill(trial.begin() + static_cast<std::ptrdiff_t>(compiled.primal_size), trial.end(), 0.0);
            const auto selected = selected_chemical_phase ? std::optional<std::size_t>(index) : std::nullopt;
            const auto& parameters = compiled.model->parameters();
            const bool gradient_only = selected && parameters.site_pathway.empty() && parameters.site_correction.empty();
            const auto prepared = compiled.prepare_values(trial, true, observations ? &observations->first_numerical_failure : nullptr,
                selected, {}, false, gradient_only);
            const auto rows = compiled.evaluate_at(prepared, selected).values;
            std::vector<double> c(m);
            for (std::size_t j = 0; j < m; ++j) c[j] = rows[chemical[j]] - trial[offset + j];
            if (!ideal_speciation(B, compiled.balance, c, n, observations)) break;
            speciated = n;
            speciated_density = density;
            for (std::size_t j = 0; j < m; ++j) composition[phase.support[j]] = n[j];
        }
    } catch (const std::exception& error) {
        if (observations) retain_density_failure(observations->first_numerical_failure, error, "speciation seed");
    }
    if (speciated.empty()) return std::nullopt;
    n = speciated;
    const double total = std::accumulate(n.begin(), n.end(), 0.0);
    for (std::size_t j = 0; j < m; ++j) point[offset + j] = std::log(n[j] / total);
    scope.probe("speciated_amount", total);
    return std::make_pair(total, *speciated_density);
}

// ALG-EQUILIBRIUM true-species split: a reacting flash of one free liquid and one free uncharged vapor at fixed T and P
// starts from the ideal speciation (z*, N*) of the feed as that liquid, frozen and split by Rachford-Rice with
// K_i = phi_i^L(z*)/phi_i^V(y) on the vapor support and K_i = 0 on liquid-confined species. y only estimates phi^V:
// two substitution passes y ~ z* phi^L/phi^V; the split's vapor is K_i x_i. The split conserves N* z*, so its charge
// and, to the speciation's accuracy, B n = b. Empty outside that domain or without an interior split: the
// apparent-species start follows, before the one solve.
std::optional<PhaseStart> true_species_split(
    const CompiledSystem& compiled, core::Scope& scope, const DriverLimits* limits, StageWork* stages) {
    const auto& problem = compiled.problem;
    if (!compiled.has_reaction || compiled.feasibility_cell() || problem.temperature.free || problem.pressure.free
        || problem.phases.size() != 2 || problem.phases[0].kind == problem.phases[1].kind) return std::nullopt;
    const std::size_t liquid = problem.phases[0].kind == PhaseKind::Liquid ? 0 : 1, vapor = 1 - liquid;
    for (const auto& phase : problem.phases)
        if (!phase.is_free() || !phase.composition_guess.empty()) return std::nullopt;
    const auto& ls = problem.phases[liquid].support;
    const auto& vs = problem.phases[vapor].support;
    // Compilation puts every fed or reacting, non-absent species on an effective support and only
    // zero-charge species on a vapor support; vapor inside liquid then puts each on the liquid support.
    for (const std::size_t i : vs)
        if (std::find(ls.begin(), ls.end(), i) == ls.end()) return std::nullopt;
    const std::size_t components = compiled.model->parameters().component_count();
    std::vector<double> z(components, 0.0), k(components, 0.0), y(components, 0.0), ln_phi_vapor(components, 0.0);
    double total = 0.0, vapor_density = 0.0;
    try {
        auto point = compiled.generic_point();
        const auto speciated = speciate_single_phase(compiled, point, compiled.feed_composition, scope, liquid, stages, true);
        if (limits) limits->check();
        if (!speciated) return std::nullopt;
        total = speciated->first;
        for (std::size_t j = 0; j < ls.size(); ++j) z[ls[j]] = std::exp(point[compiled.phase_offsets()[liquid] + j]);
        const double t = compiled.seed_temperature(), p = compiled.seed_pressure();
        const auto liquid_state = compiled.model->state_at_pressure(t, p, z, DensityBranch::Liquid);
        if (stages) core::retain_first(stages->first_numerical_failure, liquid_state.first_numerical_failure);
        if (!liquid_state.log_fugacity_coefficient) return std::nullopt;
        const auto& ln_phi_liquid = *liquid_state.log_fugacity_coefficient;
        for (int pass = 0; pass < 2; ++pass) {
            double sum = 0.0;
            for (const std::size_t i : vs) sum += y[i] = z[i] * std::exp(require_component_property(ln_phi_liquid[i]) - ln_phi_vapor[i]);
            for (const std::size_t i : vs) y[i] /= sum;
            const auto gas = [&] {  // the first pass's root continues to the second y; outside its window, a new scan
                if (pass > 0) try { return compiled.model->state_at_pressure(t, p, y, DensityBranch::Vapor, vapor_density); }
                    catch (const std::exception& error) { if (stages) retain_density_failure(stages->first_numerical_failure, error, "split vapor anchor"); }
                return compiled.model->state_at_pressure(t, p, y, DensityBranch::Vapor);
            }();
            if (stages) core::retain_first(stages->first_numerical_failure, gas.first_numerical_failure);
            if (!gas.log_fugacity_coefficient) return std::nullopt;
            vapor_density = gas.molar_density;
            for (const std::size_t i : vs) ln_phi_vapor[i] = require_component_property((*gas.log_fugacity_coefficient)[i]);
        }
        for (const std::size_t i : vs) {
            k[i] = std::exp(require_component_property(ln_phi_liquid[i]) - ln_phi_vapor[i]);
            if (!(std::isfinite(k[i]) && k[i] > 0.0)) return std::nullopt;
        }
    } catch (const std::exception& error) {  // an exhausted budget stops the start; any other failure keeps the fallback
        if (stages) retain_density_failure(stages->first_numerical_failure, error, "split seed");
        if (limits) limits->check();
        return std::nullopt;
    }
    auto split = rachford_rice(compiled, z, k, total, scope);
    if (!split) return std::nullopt;
    split->vapor_density = vapor_density;
    scope.probe("start_rule", 4.0);
    return split;
}

// ---------------------------------------------------------------- continuation


struct BoundaryAttempt final {
    StageWork observations;
    bool accepted = false;
    bool terminal = false;
    ContinuationStop stop = ContinuationStop::stalled;
    std::string message;
    std::optional<EquilibriumResult> result;
    std::optional<core::Convergence> convergence;
};

bool same_intensive(const FreeValue& a, const FreeValue& b, bool endpoint) {
    if (a.free != b.free) return false;
    return !a.free || !endpoint || (a.absent == b.absent && a.value == b.value);
}

std::string layout_difference(const CompiledSystem& a, const CompiledSystem& b) {
    for (std::size_t p = 0; p < a.problem.phases.size(); ++p)
        if (a.problem.phases[p].support != b.problem.phases[p].support)
            return "admitted support of phase " + a.problem.phases[p].name + " differs";
    if (a.rows.size() != b.rows.size()) return "row count differs";
    for (std::size_t r = 0; r < a.rows.size(); ++r)
        if (a.rows[r].kind != b.rows[r].kind || a.rows[r].phase != b.rows[r].phase
            || a.rows[r].index != b.rows[r].index)
            return "row " + std::to_string(r) + " identity differs";
    if (a.primal_size != b.primal_size || a.rank != b.rank
        || a.charge_row_count != b.charge_row_count || a.pinned_count != b.pinned_count
        || a.free_intensive_count != b.free_intensive_count)
        return "coordinate layout differs";
    return "";
}

struct Path final {
    EquilibriumProblem base;
    std::vector<double> feed_start, feed_target;
    std::vector<double> parameter_start, parameter_target;
    double t1 = 0.0, p1 = 0.0, h1 = 0.0, s1 = 0.0, v1 = 0.0, u1 = 0.0;  // the start values are base's

    EquilibriumProblem at(double q) const {
        EquilibriumProblem problem = base;
        problem.feed.kind = FeedKind::Amounts;
        problem.feed.values.resize(feed_start.size());
        for (std::size_t i = 0; i < feed_start.size(); ++i)
            problem.feed.values[i] = (1.0 - q) * feed_start[i] + q * feed_target[i];
        if (!problem.temperature.free) problem.temperature = FreeValue::fixed((1.0 - q) * base.temperature.value + q * t1);
        if (!problem.pressure.free) problem.pressure = FreeValue::fixed((1.0 - q) * base.pressure.value + q * p1);
        if (base.internal_energy_j) problem.internal_energy_j = (1.0 - q) * *base.internal_energy_j + q * u1;
        if (base.volume_m3) problem.volume_m3 = (1.0 - q) * *base.volume_m3 + q * v1;
        if (base.enthalpy_j) problem.enthalpy_j = (1.0 - q) * *base.enthalpy_j + q * h1;
        if (base.entropy_j_per_k) problem.entropy_j_per_k = (1.0 - q) * *base.entropy_j_per_k + q * s1;
        return problem;
    }

    std::vector<double> lever(const CompiledSystem& compiled) const {
        std::vector<double> result(compiled.control_width(), 0.0);
        result[0] = base.temperature.free ? 0.0 : t1 - base.temperature.value;
        result[1] = base.pressure.free ? 0.0 : p1 - base.pressure.value;
        for (std::size_t i = 0; i < feed_start.size(); ++i) result[2 + i] = feed_target[i] - feed_start[i];
        for (std::size_t i = 0; i < parameter_target.size(); ++i)
            result[2 + feed_start.size() + i] = parameter_target[i] - parameter_start[i];
        if (!base.enthalpy_j && !base.entropy_j_per_k && !base.volume_m3 && !base.internal_energy_j) return result;
        result.resize(compiled.extended_control_width(), 0.0);
        const auto target = result.size() - 2 - base.volume_m3.has_value() - base.internal_energy_j.has_value();
        if (base.volume_m3) result[target+2] = v1 - *base.volume_m3;
        if (base.internal_energy_j) result.back() = u1 - *base.internal_energy_j;
        result[target] = base.enthalpy_j ? h1 - *base.enthalpy_j : 0.0;
        result[target + 1] = base.entropy_j_per_k ? s1 - *base.entropy_j_per_k : 0.0;
        return result;
    }
};

// Molar density of a prepared phase on the driver's log-composition chart.
double chart_density(const HelmholtzModel& model, const PhasePoint& prepared, std::size_t phase) {
    const auto& x = prepared.composition[phase];
    return (model.parameters().ideal_gas ? std::accumulate(x.begin(), x.end(), 0.0) : 1.0) / prepared.volume[phase];
}

// A predicted free amount below zero is projected onto its bound rather than rejected: the
// amount bound is the disappearance boundary itself, and only the bound-constrained corrector
// can report the typed event there (counted in `projected`). A log mole fraction above zero is
// projected onto its corrector bound zero. Every other coordinate outside its domain is rejected.
bool inside_domain(const CompiledSystem& compiled, std::vector<double>& y, std::size_t& projected, bool affine = false) {
    if (!std::all_of(y.begin(), y.end(), [](double v) { return std::isfinite(v); })) return false;
    std::size_t offset = 0;
    for (const auto& phase : compiled.problem.phases) {
        for (std::size_t i = 0; i < phase.support.size(); ++i, ++offset) {
            if (y[offset] > 0.0) { y[offset] = 0.0; ++projected; }
        }
        if (phase.amount.kind == AmountKind::Free) {
            if (y[offset] < 0.0) { y[offset] = 0.0; ++projected; }
            ++offset;
        }
        const double coordinate = std::exp(y[offset++]);
        if (compiled.model->parameters().ideal_gas) {
            if (!(coordinate > 0.0) || !std::isfinite(coordinate)
                || !std::isfinite(1.0 / coordinate)) return false;
        } else if (!admitted_packing_fraction(coordinate)) return false;
    }
    offset += compiled.mode_count;
    if (affine && (y[offset] < compiled.problem.loading->lower || y[offset++] > compiled.problem.loading->upper)) return false;
    if (compiled.problem.temperature.free && !(y[offset++] > 0.0)) return false;
    if (compiled.problem.pressure.free && !(y[offset++] > 0.0)) return false;
    return true;
}

// Passage through a critical point between consecutive accepted states swaps a pair's
// identities: every nonzero component of the inter-phase difference vector (delta ln eta and
// each delta x_i) changes sign while both states pass the coincidence predicate. That is the
// pair's coalescence (owner decision 2026-09-16); zero components carry no sign information.
std::optional<std::pair<std::size_t, std::size_t>> swapped_pair(
    const CompiledSystem& compiled, const std::vector<double>& before, const std::vector<double>& after
) {
    // Accepted states only: every chart composition has a positive total.
    const auto a = check_topology(compiled, before);
    const auto b = check_topology(compiled, after);
    for (std::size_t left = 0; left + 1 < a.compositions.size(); ++left) {
        for (std::size_t right = left + 1; right < a.compositions.size(); ++right) {
            std::vector<double> before_delta{a.log_packing[left] - a.log_packing[right]};
            std::vector<double> after_delta{b.log_packing[left] - b.log_packing[right]};
            for (std::size_t i = 0; i < a.compositions[left].size(); ++i) {
                before_delta.push_back(a.compositions[left][i] - a.compositions[right][i]);
                after_delta.push_back(b.compositions[left][i] - b.compositions[right][i]);
            }
            bool informative = false, swapped = true;
            for (std::size_t i = 0; i < before_delta.size(); ++i) {
                if (before_delta[i] == 0.0 || after_delta[i] == 0.0) continue;
                informative = true;
                swapped = swapped && (before_delta[i] > 0.0) != (after_delta[i] > 0.0);
            }
            if (informative && swapped) return std::make_pair(left, right);
        }
    }
    return std::nullopt;
}



ContinuationStop stop_of(const EquilibriumResult& result, bool& terminal) {
    if (result.topology_event == TopologyEvent::PhaseDisappeared) { terminal = true; return ContinuationStop::phase_disappeared; }
    if (result.topology_event == TopologyEvent::PhasesCoalesced) { terminal = true; return ContinuationStop::phases_coalesced; }
    if (result.spinodal && !result.success && result.convergence.stop == core::Stop::left_domain) return ContinuationStop::branch_lost;
    switch (result.convergence.stop) {
    case core::Stop::nonfinite: terminal = true; return ContinuationStop::nonfinite;
    case core::Stop::infeasible: terminal = true; return ContinuationStop::infeasible;
    case core::Stop::left_domain: return ContinuationStop::left_domain;
    case core::Stop::singular: return ContinuationStop::singular_jacobian;
    default: return ContinuationStop::stalled;
    }
}

// `supplied`: a caller's result; the driver's own results already carry this problem's classified rows.
BoundaryAttempt admit_result(
    const CompiledSystem& compiled, EquilibriumResult result,
    const std::vector<double>& anchors, bool initial, core::Scope& scope, bool supplied = false
) {
    BoundaryAttempt outcome;
    outcome.convergence = result.convergence;
    const auto refuse = [&](ContinuationStop stop, std::string message) {
        outcome.stop = stop;
        outcome.message = std::move(message);
        if (!core::CalculationWork::stopped())
            core::retain_failure(result.first_numerical_failure, outcome.message, "continuation admission");
        outcome.observations.observe(result);
        outcome.result = std::move(result);
        return std::move(outcome);
    };
    scope.probe("driver_result_present", 1.0);
    if (std::isfinite(result.convergence.residual)) scope.probe("row_residual", result.convergence.residual);
    if (!result.has_state()) return refuse(stop_of(result, outcome.terminal), result.message);
    const double recomputed = supplied ? compiled.classified_residual_norm(result.coordinates)
        : maximum_abs(result.classified_residuals);
    if (!(result.requested_tolerance_met.value_or(false)) || !(recomputed <= compiled.problem.solver_options.value(policy::row)))
        return refuse(ContinuationStop::stalled, "available state misses the caller row tolerance: " + result.message);
    const auto topology = check_topology(compiled, result.coordinates);
    if (topology.event != TopologyEvent::None) {
        result.topology_event = topology.event;
        result.topology_phase = topology.phase;
        result.topology_pair = topology.pair;
        return refuse(stop_of(result, outcome.terminal), "supplied state fails the candidate topology predicate");
    }
    BranchCheck branch;
    try { branch = compiled.mode_count ? admit_null_mode(compiled, result) : check_density_branches(compiled, result, anchors, initial, scope); }
    catch (const std::exception& error) {
        if (compiled.problem.spinodal && compiled.spinodal_diagnostics() && !std::isfinite(compiled.spinodal_diagnostics()->seed))
            result = spinodal_search_failure(compiled, std::move(result), error.what());
        branch.message = error.what();
    }
    outcome.terminal = branch.terminal;
    if (!branch.accepted) return refuse(ContinuationStop::branch_lost, branch.message + "; driver: " + result.message);
    outcome.accepted = true;
    outcome.message = branch.message;
    outcome.observations.observe(result);
    outcome.result = std::move(result);
    return outcome;
}

ContinuationResult finish(ContinuationResult result, ContinuationStop stop, std::string message) {
    result.stop = stop;
    result.success = stop == ContinuationStop::reached_target || stop == ContinuationStop::critical_endpoint;
    result.message = std::string(continuation_stop_name(stop)) + (message.empty() ? "" : ": " + message);
    return result;
}

void reject(ContinuationResult& result, const BoundaryAttempt& attempt, double q) {
    result.attempted_path = q;
    result.cause_stop = attempt.stop;
    result.cause_convergence = attempt.convergence;
    result.cause_message = attempt.message;
    result.last_result = attempt.result;
    if (attempt.result.has_value()) {
        result.last_result->optimizer_state.reset();
        result.event_phase = attempt.result->topology_phase;
        result.event_pair = attempt.result->topology_pair;
    }
}

// One predictor-corrector attempt shared by the straight path and the envelope: generated layout and
// rank, the Euler predictor over `delta` along the implicit tangent at the accepted state, the admitted
// domain and the same-driver corrector. Continuation uses the previous or predicted density as a branch witness.
// A declared target independently selects its density roots.
enum class DensityAdmission { Previous, Predicted, Declared };
BoundaryAttempt advance(const CompiledSystem& start, const CompiledSystem& current,
                const std::vector<double>& y, const CompiledSystem& compiled, const EquilibriumProblem& problem,
                const std::function<std::vector<double>(std::optional<core::Convergence>*)>& tangent_at, double delta,
                std::vector<double> anchors, DensityAdmission admission, core::Scope& scope, std::shared_ptr<const OptimizerState> warm,
                const DriverLimits* limits = nullptr, StageWork* observations = nullptr, const SpinodalSearch* carried = nullptr) {
    BoundaryAttempt attempt;
    attempt.terminal = true;
    std::vector<double> tangent;
    if (const std::string generated = layout_difference(start, compiled); !generated.empty()) {
        attempt.stop = ContinuationStop::invalid_path;
        attempt.message = generated;
    } else if (!compiled.feasible) {
        // Same layout, generic-point rank loss at the generated problem.
        attempt.stop = ContinuationStop::singular_jacobian;
        attempt.observations.observe(structural_result(compiled, scope));
        attempt.message = "generated problem loses generic rank (nullity " + std::to_string(compiled.generic_nullity) + ")";
    } else if ((tangent = tangent_at(&attempt.convergence)).empty()) {
        // A state-local tangent cannot change by halving the step: singular is terminal
        // (owner decision 2026-09-16), the branch end at that accepted state.
        const auto jacobian = current.jacobian(y);
        const bool numeric = std::all_of(jacobian.begin(), jacobian.end(), [](double v) { return std::isfinite(v); });
        const auto info = inspect_sensitivity(current, y);
        attempt.stop = numeric ? ContinuationStop::singular_jacobian : ContinuationStop::nonfinite;
        attempt.message = numeric
            ? "implicit tangent is unavailable at the accepted state (rank " + std::to_string(info.rank) + " of "
                + std::to_string(y.size()) + ", condition estimate " + core::optional_number(info.condition_estimate) + ")"
            : "nonfinite operands in the implicit tangent";
    } else {
        auto predicted = shifted_coordinates(y, tangent, delta);
        std::size_t projected = 0;
        const bool inside = inside_domain(compiled, predicted, projected);
        scope.probe("predictor_projected", static_cast<double>(projected));
        if (inside) {
            std::optional<CompiledSystem> selected;
            if (compiled.problem.spinodal) {
                SpinodalSearch scan;
                try {
                    const auto offset = compiled.phase_offsets()[0] + compiled.feed_composition.size() + 1;
                    selected = spinodal_corrector(compiled, limits, anchors[0], carried, scan);
                    predicted[offset] = std::log(scan.seed * compiled.model->packing_fraction_coefficient(compiled.seed_temperature(), compiled.feed_composition));
                } catch (const std::bad_alloc&) { throw; }
                catch (const std::exception& error) {
                    return admit_result(compiled, spinodal_search_failure(compiled, structural_result(compiled, scope), error.what(), carried ? &scan : nullptr), anchors, false, scope);
                }
            }
            const auto& corrector = selected ? *selected : compiled;
            auto solved = unless_failed(solve_equilibrium_from_coordinates(corrector, predicted, nullptr, limits, std::move(warm)), observations);
            if (admission != DensityAdmission::Predicted || compiled.mode_count)
                return admit_result(corrector, std::move(solved), anchors, admission == DensityAdmission::Declared, scope);
            auto previous = admit_result(compiled, solved, anchors, false, scope);
            if (previous.accepted || previous.stop != ContinuationStop::branch_lost) return previous;
            const auto at = compiled.prepare(predicted, true);
            for (std::size_t p = 0; p < anchors.size(); ++p) anchors[p] = chart_density(*compiled.model, at, p);
            auto witnessed = admit_result(compiled, std::move(solved), anchors, false, scope);
            return witnessed.accepted ? witnessed : previous;
        }
        attempt.terminal = false;
        attempt.stop = ContinuationStop::left_domain;
        attempt.message = "predicted coordinates leave the admitted domain";
    }
    core::retain_first(attempt.observations.first_numerical_failure, attempt.convergence);
    core::retain_failure(attempt.observations.first_numerical_failure, attempt.message, "continuation predictor");
    return attempt;
}

// An identity swap between accepted states without a requested ln K sign change is the pair's coalescence.
void coalesced(BoundaryAttempt& attempt, const CompiledSystem& compiled, std::pair<std::size_t, std::size_t> pair) {
    attempt.accepted = false;
    attempt.terminal = true;
    attempt.stop = ContinuationStop::phases_coalesced;
    attempt.result->topology_event = TopologyEvent::PhasesCoalesced;
    attempt.result->topology_pair = pair;
    attempt.message = "phases " + compiled.problem.phases[pair.first].name + " and " + compiled.problem.phases[pair.second].name
        + " swapped identity between accepted states: passage through a critical point";
}

// J dy/ds = -dF/ds of an envelope cell, dF/ds = -1 on its condition row (the last row).
std::vector<double> envelope_tangent(const CompiledSystem& compiled, const PhasePoint& prepared, ContinuationResult& result) {
    std::vector<double> rows(compiled.rows.size(), 0.0);
    rows.back() = -1.0;
    std::optional<core::Convergence> failure;
    auto tangent = solve_implicit(compiled, prepared, rows, &failure);
    core::retain_first(result.observations.first_numerical_failure, failure);
    if (tangent.empty()) result.cause_convergence = std::move(failure);
    return tangent;
}

// Value and d/ds of ln T, ln P, ln rho of phases 0 and 1, then ln K_i = ln x_i^inc - ln x_i^feed on the
// shared support, from the geometry owner's first Taylor coefficient along the tangent.
std::vector<std::array<double, 2>> envelope_coordinates(const CompiledSystem& compiled, const PhasePoint& prepared,
                                                        const std::vector<double>& tangent) {
    const auto& y = prepared.point;
    const auto series = compiled.row_series(prepared, {y, tangent},
        std::vector<double>(compiled.control_width(), 0.0), 1);
    std::vector<std::array<double, 2>> result{{std::log(series.temperature[0]), series.temperature[1] / series.temperature[0]},
                                              {std::log(series.pressure[0]), series.pressure[1] / series.pressure[0]}};
    for (const auto& phase : series.phase) {
        double n0 = 0.0, n1 = 0.0;  // rho = sum(x)/V on the unnormalized log-composition chart
        for (std::size_t i = 0; i < compiled.n_feed.size(); ++i) { n0 += phase.point[0][2 + i]; n1 += phase.point[1][2 + i]; }
        result.push_back({std::log(n0 / phase.point[0][1]), n1 / n0 - phase.point[1][1] / phase.point[0][1]});
    }
    const std::size_t pin = compiled.problem.phases[0].is_pinned() ? 0 : 1;  // the zero-pinned phase
    const auto &pinned = prepared.shape[pin].x, &feed = prepared.shape[1 - pin].x;
    for (std::size_t i = 0; i < feed.size(); ++i)
        result.push_back({y[pinned[i]] - y[feed[i]], tangent[pinned[i]] - tangent[feed[i]]});
    return result;
}

// Cubic Hermite value at s = 0 from (value, d/ds) at s0 and s1 in one coordinate.
double hermite(double s0, double s1, const std::array<double, 2>& a, const std::array<double, 2>& b) {
    const double h = s1 - s0, t = -s0 / h;
    return (2 * t * t * t - 3 * t * t + 1) * a[0] + (t * t * t - 2 * t * t + t) * h * a[1]
        + (-2 * t * t * t + 3 * t * t) * b[0] + (t * t * t - t * t) * h * b[1];
}

// The feed's critical cell seeded at (ln T, ln P, ln rho*) with u = D^-1 dn of the common-volume displacement
// dn = V*(rho_inc x_inc - rho_feed x_feed) at `state` (V* scales out; orientation kept).
EquilibriumResult critical_candidate(const HelmholtzModel& model, const CompiledSystem& compiled,
                                     const EquilibriumResult& state, const std::array<double, 3>& seed) {
    const std::size_t c = compiled.n_feed.size(), pin = compiled.problem.phases[0].is_pinned() ? 0 : 1;
    EquilibriumProblem cell = compiled.problem;
    cell.envelope.reset();
    cell.phases = {compiled.problem.phases[1 - pin]};
    cell.temperature = FreeValue::variable(std::exp(seed[0]));
    cell.pressure = FreeValue::variable(std::exp(seed[1]));
    cell.critical = CriticalDeclaration{std::exp(seed[2]), std::vector<double>(c)};
    for (std::size_t i = 0; i < c; ++i)
        cell.critical->mode_guess[i] = (state.molar_densities[pin] * state.mole_fractions[pin * c + i]
            - state.molar_densities[1 - pin] * state.mole_fractions[(1 - pin) * c + i])
            / std::sqrt(compiled.feed_composition[i]);
    return solve_compiled(compile_problem(model, cell));
}

// Envelope mode of continue_equilibrium (ALG-EQUILIBRIUM; Michelsen 1980 local parametrization in the raw T/P
// chart): the condition target s is the path parameter and the trace first increases the declared
// coordinate. A mixture closes at the initial pressure on the other branch after one accepted ln K straddle;
// a pure fluid stops at its separately solved critical candidate. Anything else is a diagnosed prefix.
ContinuationResult trace_envelope(const HelmholtzModel& model, const ContinuationRequest& request) {
    core::Scope scope("equilibrium_envelope");
    ContinuationResult result;
    auto& observations = result.observations;
    const auto& options = request.options;
    const double window = policy::DENSITY_ANCHOR_LOG_WINDOW.value, straddle_band = options.envelope_straddle;
    // A pure trace ends within 4 W of coalescence; below about 2 W the window refuses both roots.
    const double coalescence_gap = 4.0 * window;
    if (!(options.envelope_max_log_step > 0.0) || !(straddle_band > 0.0)
        || !(straddle_band + straddle_band <= options.envelope_max_log_step)  // the straddle jump fits one step
        || !(options.minimum_step > 0.0) || options.max_steps == 0)
        throw std::invalid_argument("invalid envelope options: need envelope_straddle > 0, 2 envelope_straddle <= "
            "envelope_max_log_step, minimum_step > 0 and max_steps > 0");
    EquilibriumProblem base = request.start_problem;
    const CompiledSystem start = compile_problem(model, base);
    if (!start.feasible) { observations.observe(structural_result(start, scope)); return finish(std::move(result), ContinuationStop::invalid_path, "envelope cell is not full rank"); }
    const bool pure = start.n_feed.size() == 1;
    auto admitted = admit_result(start, unless_failed(solve_equilibrium_from_coordinates(model, base,
        request.start_result ? request.start_result->coordinates : boundary_start(start, true, nullptr, &observations)),
        &observations), {}, true, scope);
    observations.observe(admitted.observations);
    if (!admitted.accepted) {
        reject(result, admitted, base.envelope->target);
        return finish(std::move(result), ContinuationStop::initial_failure, admitted.message);
    }
    CompiledSystem current = start;
    std::vector<double> y = admitted.result->coordinates, motion;
    const double closing_pressure = std::log(*admitted.result->pressure);
    result.states.push_back(std::move(*admitted.result));
    result.path_values.push_back(base.envelope->target);
    result.conditions.push_back(*base.envelope);
    double step = 0.0;  // zero until a step's predictor error is measured
    bool passed = false, closing = false;
    std::size_t local_retries = 0;
    bool transfer = true;
    while (result.accepted_steps < options.max_steps) {
        const auto at = current.prepare(y);
        auto tangent = envelope_tangent(current, at, result);
        if (tangent.empty()) {
            return finish(std::move(result), ContinuationStop::singular_jacobian, "envelope tangent is unavailable at the accepted state");
        }
        auto c = envelope_coordinates(current, at, tangent);
        // Orientation by continuity of the dimensionless (d ln T, d ln P, d ln K)/ds.
        std::vector<double> heading = {c[0][1], c[1][1]};
        for (std::size_t k = 4; k < c.size(); ++k) heading.push_back(c[k][1]);
        if (motion.empty()) motion = heading;
        double sign = std::inner_product(heading.begin(), heading.end(), motion.begin(), 0.0) < 0.0 ? -1.0 : 1.0;
        for (double& v : heading) v *= sign;
        // The largest dimensionless tangent component selects the coordinate (Michelsen 1980, Eq. 18); the
        // closing corrector holds ln P. A switch rescales the same tangent and restarts s at the current value.
        std::size_t selected = 1;
        for (std::size_t k = 0; !closing && k < c.size(); ++k)
            if (k != 2 && k != 3 && std::abs(c[k][1]) > std::abs(c[selected][1])) selected = k;
        // The accepted state's own condition: a rejected attempt leaves it, and the tangent, unchanged.
        EnvelopeCondition active = *current.problem.envelope;
        const EnvelopeCondition condition{selected == 0 ? EnvelopeCoordinate::LogTemperature
            : selected == 1 ? EnvelopeCoordinate::LogPressure : EnvelopeCoordinate::LogK,
            selected < 4 ? 0 : current.problem.phases[0].support[selected - 4], c[selected][0]};
        if (condition.kind != active.kind || condition.component != active.component) {
            const double slope = c[selected][1];
            for (double& v : tangent) v /= slope;
            for (auto& v : c) v[1] /= slope;
            sign *= slope < 0.0 ? -1.0 : 1.0;
            active = condition;
            scope.probe("condition_switch", static_cast<double>(selected));
        }
        const double s = active.target, gap = c[2][0] - c[3][0];
        if (pure && std::abs(gap) <= coalescence_gap) {
            // The first accepted state within 4 W of coalescence ends the saturation trace at the separately
            // solved critical candidate seeded there; the last distinct saturation state is retained.
            const std::array<double, 3> seed{c[0][0], c[1][0], 0.5 * (c[2][0] + c[3][0])};
            result.critical_seeds.push_back({seed.begin(), seed.end()});
            auto candidate = critical_candidate(model, current, result.states.back(), seed);
            if (!candidate.success) return finish(std::move(result), ContinuationStop::stalled, "pure critical candidate: " + candidate.message);
            result.critical_points.push_back(std::move(candidate));
            return finish(std::move(result), ContinuationStop::critical_endpoint, "last distinct saturation state retained");
        }
        // A small ln K decreasing in magnitude (s sign < 0) is approached at the same-sign threshold; only the
        // straddle to -s crosses zero, and s = 0 is never a target.
        const bool toward_zero = active.kind == EnvelopeCoordinate::LogK && !passed && s * sign < 0.0;
        const bool straddle = toward_zero && std::abs(s) <= straddle_band;
        double s_try = closing ? closing_pressure : -s;
        if (!closing && !straddle) {
            // First step: density change at most 0.9 W; later: the predictor error; a pure fluid: 3/8 of the gap.
            double delta = step > 0.0 ? step
                : std::min(options.envelope_max_log_step, 0.9 * window / std::max(std::abs(c[2][1]), std::abs(c[3][1])));
            if (pure) delta = std::min(delta, 0.375 * std::abs(gap / (c[2][1] - c[3][1])));
            s_try = s + sign * delta;
            if (toward_zero && (std::abs(s_try) < straddle_band || s_try * s <= 0.0)) s_try = std::copysign(straddle_band, s);
        }
        auto problem = base;
        problem.envelope = active;
        problem.envelope->target = s_try;
        // Each anchor is the Euler-predicted density, also through the straddle's crossing.
        std::vector<double> expected{std::exp(c[2][0] + (s_try - s) * c[2][1]), std::exp(c[3][0] + (s_try - s) * c[3][1])};
        if (straddle) {  // the straddle attempt itself expects the feed on the other density branch
            for (auto& phase : problem.phases) phase.kind = phase.kind == PhaseKind::Liquid ? PhaseKind::Vapor : PhaseKind::Liquid;
        }
        ++result.attempts;
        core::Scope attempt_scope("equilibrium_envelope_attempt");
        attempt_scope.probe("path_value", s_try);
        CompiledSystem compiled = compile_problem(model, problem);
        auto attempt = advance(start, current, y, compiled, problem, [&](auto*) { return tangent; }, s_try - s, expected, DensityAdmission::Previous, attempt_scope, transfer ? result.states.back().optimizer_state : nullptr, nullptr, &observations);
        observations.observe(attempt.observations);
        if (attempt.accepted) {
            const auto pair = swapped_pair(compiled, y, attempt.result->coordinates);
            if (pair && !straddle) coalesced(attempt, compiled, *pair);
            else if (straddle && !pair) {
                attempt.accepted = false;
                attempt.stop = ContinuationStop::stalled;
                attempt.message = "requested ln K sign change produced no identity swap";
            } else if (straddle) {
                // Both endpoint tangents in the same s; cubic Hermite values at s = 0 seed the critical cell.
                const auto passed_at = compiled.prepare(attempt.result->coordinates);
                const auto far = envelope_tangent(compiled, passed_at, result);
                if (far.empty()) return finish(std::move(result), ContinuationStop::singular_jacobian, "passage tangent unavailable");
                const auto there = envelope_coordinates(compiled, passed_at, far);
                const std::array<double, 3> seed{hermite(s, s_try, c[0], there[0]), hermite(s, s_try, c[1], there[1]),
                    0.5 * (hermite(s, s_try, c[2], there[2]) + hermite(s, s_try, c[3], there[3]))};
                result.critical_seeds.push_back({seed.begin(), seed.end()});
                auto candidate = critical_candidate(model, compiled, *attempt.result, seed);
                if (!candidate.success) return finish(std::move(result), ContinuationStop::stalled, "critical candidate: " + candidate.message);
                result.critical_points.push_back(std::move(candidate));
                passed = true;
            }
        }
        attempt_scope.probe("accepted", attempt.accepted ? 1.0 : 0.0);
        if (attempt.convergence.has_value()) attempt_scope.finish(*attempt.convergence);
        if (attempt.accepted) {
            const bool closed = closing;
            // The mixture closes when an accepted step after the passage brackets the initial pressure.
            closing = passed && (c[1][0] - closing_pressure) * (std::log(*attempt.result->pressure) - closing_pressure) <= 0.0;
            ++result.accepted_steps;
            local_retries = 0;
            transfer = true;
            // Euler error grows as h^2: the next step targets a predictor density error of W/2.
            const auto& rho = attempt.result->molar_densities;
            const double e = std::max(std::abs(std::log(rho[0] / expected[0])), std::abs(std::log(rho[1] / expected[1])));
            step = std::min(options.envelope_max_log_step, std::abs(s_try - s) * std::min(2.0, std::sqrt(window / (2.0 * e))));
            motion = heading;
            y = attempt.result->coordinates;
            base = problem;
            current = std::move(compiled);
            result.path_values.push_back(s_try);
            result.conditions.push_back(*base.envelope);
            result.reached_path = s_try;
            result.states.push_back(std::move(*attempt.result));
            if (closed) return finish(std::move(result), ContinuationStop::reached_target, "");
            continue;
        }
        transfer = false;
        core::retain_failure(observations.first_numerical_failure, attempt.message, "continuation admission");
        reject(result, attempt, s_try);
        // A straddle or closing target cannot change by halving: its refusal ends the trace as a prefix.
        // An ordinary step that leaves the branch window or whose corrector ends on coincident phases halves.
        attempt.terminal &= attempt.stop != ContinuationStop::branch_lost
            && (attempt.stop != ContinuationStop::phases_coalesced || attempt.result->has_state());
        if (attempt.terminal || straddle || closing) return finish(std::move(result), attempt.stop, attempt.message);
        step = std::abs(s_try - s) / 2.0;
        if (step < options.minimum_step) return finish(std::move(result), ContinuationStop::step_too_small, attempt.message);
        if (local_retries++ >= options.max_retries_per_step)
            return finish(std::move(result), ContinuationStop::budget_exhausted, attempt.message);
        ++result.retries;
    }
    return finish(std::move(result), ContinuationStop::budget_exhausted, "max_steps reached before the envelope closed");
}

}  // namespace

struct SpinodalKey final {
    struct Request { const CompiledSystem& compiled; const DriverLimits* limits; bool corrector; };
    using Entry = core::Retained<SpinodalKey, SpinodalSearch>;
    std::shared_ptr<const char> producer;
    double temperature, guess;
    std::vector<double> feed, nodes;
    std::vector<std::size_t> support;
    bool corrector;

    bool matches(const Request& r) const {
        const auto& c = r.compiled;
        return c.model && producer == c.model->identity() && c.problem.spinodal && corrector == r.corrector
            && core::bitwise_equal(temperature, c.seed_temperature())
            && core::bitwise_equal(guess, c.problem.spinodal->density_guess_mol_per_m3)
            && core::bitwise_equal(feed, c.feed_composition) && core::bitwise_equal(nodes, c.problem.spinodal->density_nodes)
            && support == c.problem.phases[0].support;
    }
    static bool admits(const Request& r, const SpinodalSearch& result) {
        return std::isfinite(result.seed) && !(r.limits && r.limits->expired());
    }
    static auto resolve(const CompiledSystem& c, const DriverLimits* limits) {
        return c.spinodal_selection_.resolve(Request{c, limits, c.spinodal_corrector_});
    }
    static void publish(const CompiledSystem& c, const SpinodalSearch& report, bool corrector) {
        c.spinodal_selection_ = Entry::publish({c.model->identity(), c.seed_temperature(), c.problem.spinodal->density_guess_mol_per_m3,
            c.feed_composition, c.problem.spinodal->density_nodes, c.problem.phases[0].support, corrector}, report);
        c.spinodal_corrector_ = corrector;
    }
    static auto& begin_diagnostics(const CompiledSystem& c) { return c.spinodal_diagnostics_.emplace(); }
};

SpinodalSearch spinodal_search(const CompiledSystem& compiled, const DriverLimits* limits, double anchor, const SpinodalSearch* carried, SpinodalSearch* attempted) {
    if (!carried) if (const auto admitted = SpinodalKey::resolve(compiled, limits)) return *admitted;
    const auto previous_selection = carried && carried == attempted ? std::optional(*carried) : std::nullopt;
    if (previous_selection) carried = &*previous_selection;
    SpinodalSearch local;
    auto& scan = carried ? (attempted ? *attempted : local) : SpinodalKey::begin_diagnostics(compiled);
    if (carried && attempted) scan = {};
    const auto& declaration = *compiled.problem.spinodal;
    if (!carried) anchor = declaration.density_guess_mol_per_m3;
    const auto& x = compiled.feed_composition;
    const double T = compiled.seed_temperature();
    scan.geometric_cap = compiled.model->molar_density_upper_bound(T, x);
    const auto sample = [&](double rho) {
        if (limits) limits->check();
        ++scan.evaluations;
        if (!admitted_packing_fraction(compiled.model->packing_fraction_coefficient(T, x) * rho))
            throw std::domain_error("density outside packing domain (Hessian request)");
        return helmholtz_spectrum(compiled.model->evaluate(T, 1 / rho, x, DerivativeOrder::Hessian), x,
            compiled.problem.phases[0].support);
    };
    std::optional<HelmholtzSpectrum> previous;
    double last = NAN;
    for (double rho : declaration.density_nodes) {
        std::optional<HelmholtzSpectrum> at;
        try { at = sample(rho); }
        catch (const std::bad_alloc&) { throw; }
        catch (const BudgetExhausted&) { throw; }
        catch (const std::exception& error) {
            DensityDiagnosticPoint limit;
            limit.lower = previous ? last : NAN; limit.upper = rho;
            limit.cause = error.what(); limit.stop = "hessian_unavailable";
            scan.finite_limits.push_back(std::move(limit));
        }
        if (at && std::abs(at->values[0]) <= at->band) { scan.unresolved_densities.push_back(rho); continue; }
        if (previous && at && ((previous->values[0] > previous->band && at->values[0] < -at->band)
            || (previous->values[0] < -previous->band && at->values[0] > at->band))) {
            scan.brackets.push_back({last, rho, previous->values[0], at->values[0]});
        }
        previous = std::move(at); last = rho;
    }
    double distance = INFINITY;
    for (const auto& bracket : scan.brackets) {
        if (carried) {
            const auto prior = std::find_if(carried->brackets.begin(), carried->brackets.end(),
                [&](const auto& b) { return b.lower == carried->lower; });
            if (prior == carried->brackets.end() || std::signbit(prior->lower_eigenvalue) != std::signbit(bracket.lower_eigenvalue)
                || bracket.lower > carried->upper || bracket.upper < carried->lower) continue;
            const double middle = std::sqrt(bracket.lower * bracket.upper);
            if (std::any_of(scan.finite_limits.begin(), scan.finite_limits.end(), [&](const auto& gap) {
                return gap.upper >= std::min(anchor, middle) && (!std::isfinite(gap.lower) || gap.lower <= std::max(anchor, middle));
            })) continue;
        }
        const double d = std::abs(std::log(std::sqrt(bracket.lower * bracket.upper) / anchor));
        if (d < distance) { distance = d; scan.lower = bracket.lower; scan.upper = bracket.upper; }
    }
    if (!std::isfinite(scan.lower)) throw std::domain_error(carried ? "spinodal search lost the carried crossing" : "no resolved connected spinodal bracket in declared scan");
    double a = scan.lower, b = scan.upper, fa = sample(a).values[0];
    for (std::size_t k = 0; k < static_cast<std::size_t>(policy::DENSITY_MAX_ITERATIONS.value); ++k) {
        const double m = std::sqrt(a * b);
        HelmholtzSpectrum at;
        try { at = sample(m); }
        catch (const std::bad_alloc&) { throw; }
        catch (const BudgetExhausted&) { scan.seed = NAN; throw; }
        catch (const std::exception& error) {
            DensityDiagnosticPoint limit;
            limit.lower = a; limit.upper = m; limit.cause = error.what(); limit.stop = "hessian_unavailable";
            scan.finite_limits.push_back(std::move(limit)); scan.seed = NAN; throw;
        }
        scan.seed = m;
        if (std::abs(at.values[0]) <= at.band || m == a || m == b) break;
        if (std::signbit(fa) != std::signbit(at.values[0])) b = m;
        else { a = m; fa = at.values[0]; }
    }
    if (limits && limits->expired()) { scan.seed = NAN; limits->check(); }
    if (!carried) SpinodalKey::publish(compiled, scan, false);
    return scan;
}

CompiledSystem spinodal_corrector(const CompiledSystem& compiled, const DriverLimits* limits, double anchor,
    const SpinodalSearch* carried, SpinodalSearch& attempted) {
    attempted = spinodal_search(compiled, limits, anchor, carried, &attempted);
    auto selected = compiled;
    SpinodalKey::begin_diagnostics(selected) = attempted;
    SpinodalKey::publish(selected, attempted, true);
    return selected;
}

std::vector<double> boundary_start(const CompiledSystem& compiled, bool speciate, const DriverLimits* limits, StageWork* stages,
                                   const std::vector<bool>& packed, bool correct_dense) {
    core::Scope scope("boundary_start");
    if (limits) limits->check();
    std::vector<double> point = compiled.generic_point();
    const auto phase_offsets = compiled.phase_offsets();
    if (compiled.mode_count) {
        // The declared T and density at the feed, P there, and u declared or the lowest mode of M (the
        // mode rows' u columns at this point).
        const auto& z = compiled.feed_composition;
        const double t = compiled.seed_temperature(), rho = compiled.problem.spinodal ? spinodal_search(compiled, limits).seed : compiled.mode_seed().density_guess_mol_per_m3;
        const std::size_t c = z.size(), mode = compiled.mode_offset();
        seed_composition(compiled, 0, z, point, phase_offsets[0]);
        point[phase_offsets[0] + c] = compiled.feed_total;
        point[phase_offsets[0] + c + 1] = std::log(compiled.model->packing_fraction_coefficient(t, z) * rho);
        if (compiled.problem.temperature.free) point[compiled.primal_size - 2] = t;
        point[compiled.primal_size - 1] = compiled.model->state(t, rho, z).pressure;
        auto u = compiled.mode_seed().mode_guess;
        if (u.empty()) {
            const auto eigen = helmholtz_spectrum(compiled.model->evaluate(t, 1 / rho, z, DerivativeOrder::Hessian),
                z, compiled.problem.phases[0].support);
            u.assign(eigen.mode.data(), eigen.mode.data() + c);
        }
        const double norm = std::sqrt(std::inner_product(u.begin(), u.end(), u.begin(), 0.0));
        for (std::size_t i = 0; i < c; ++i) point[mode + i] = u[i] / norm;
        return point;
    }
    if (compiled.problem.volume_m3 && compiled.problem.phases.size() == 1) {
        const auto at = compiled.prepare(point);
        const auto& a = at.shape[0];
        seed_composition(compiled, 0, compiled.feed_composition, point, phase_offsets[0]);
        point[*a.amount] = compiled.feed_total;
        const auto rho = compiled.feed_total / *compiled.problem.volume_m3;
        point[a.volume] = std::log(compiled.model->packing_fraction_coefficient(compiled.seed_temperature(), compiled.feed_composition) * rho);
        point[compiled.primal_size - 1] = compiled.model->state(compiled.seed_temperature(), rho, compiled.feed_composition).pressure;
        return point;
    }
    if (compiled.has_caloric_condition() || compiled.problem.volume_m3) {
        auto seed_problem = compiled.problem;
        seed_problem.volume_m3.reset();
        seed_problem.internal_energy_j.reset();
        if (compiled.problem.internal_energy_j) seed_problem.temperature = FreeValue::fixed(compiled.seed_temperature());
        seed_problem.enthalpy_j.reset();
        seed_problem.entropy_j_per_k.reset();
        const bool pure = compiled.model->parameters().component_count() == 1 && seed_problem.phases.size() == 2;
        if (pure) seed_problem.phases[1].amount = AmountSpec::pinned(0.0);
        else seed_problem.temperature = FreeValue::fixed(compiled.seed_temperature());
        if (compiled.problem.volume_m3) {
            seed_problem.pressure = pure ? FreeValue::variable(compiled.seed_pressure()) : FreeValue::fixed(compiled.seed_pressure());
        }
        const auto source = compile_problem(*compiled.model, seed_problem);
        // The seed alone stops at its own budget; a seed without a state is discarded either way (#579).
        auto seed_initial = boundary_start(source, true, limits, stages);
        DriverLimits budget = limits ? *limits : DriverLimits{std::numeric_limits<std::size_t>::max(), std::chrono::steady_clock::time_point::max()};
        budget.iterations = std::min(budget.iterations, budget.consumed + static_cast<std::size_t>(policy::EQUILIBRIUM_CALORIC_SEED_ITERATIONS.value));
        auto seed = structural_result(source, scope);
        std::exception_ptr failure;
        try {
            seed = solve_equilibrium_from_coordinates(*compiled.model, seed_problem, seed_initial, nullptr, &budget);
            failure = seed.failure;
        } catch (...) {
            failure = std::current_exception();
        }
        if (limits) {
            limits->consumed = budget.consumed;
            limits->evaluations = budget.evaluations;
            limits->hessian_calls = budget.hessian_calls;
            limits->exhausted = budget.exhausted;
        }
        if (stages) stages->add(seed);
        if (failure) std::rethrow_exception(failure);
        scope.probe("caloric_seed_cause", {}, 0, seed.message.c_str());
        if (compiled.problem.volume_m3) {
            if (!seed.has_state()) return point;
            point[compiled.primal_size - 1] = *seed.pressure;
            const auto a = source.prepare(seed.coordinates), b = compiled.prepare(point);
            for (std::size_t p = 0; p < a.shape.size(); ++p) {
                for (std::size_t i = 0; i < a.shape[p].x.size(); ++i) point[b.shape[p].x[i]] = seed.coordinates[a.shape[p].x[i]];
                point[b.shape[p].volume] = seed.coordinates[a.shape[p].volume];
                if (!pure || !compiled.problem.internal_energy_j) point[*b.shape[p].amount] = a.amount[p];
            }
            if (pure) {
                const double split = (*compiled.problem.volume_m3 / compiled.feed_total - a.volume[0]) / (a.volume[1] - a.volume[0]);
                if (split > 0.0 && split < 1.0) {
                    point[*b.shape[0].amount] = compiled.feed_total * (1.0 - split);
                    point[*b.shape[1].amount] = compiled.feed_total * split;
                } else if (!compiled.problem.internal_energy_j) throw std::domain_error("volume target outside declared coexistence interval");
            }
            return point;
        }
        const bool anchored = compiled.anchored_calorics();
        const auto inside = [&](double t) {  // every declared species' caloric record admits t
            if (anchored) try {  // records, reaction laws and completion
                return reference_calorics(compiled, t, compiled.seed_pressure(), {}).status == SolvedStateActionStatus::Available;
            } catch (const std::invalid_argument&) { return false; }
            for (const auto& phase : compiled.problem.phases)
                for (auto i : phase.support) {
                    const auto& intervals = compiled.model->thermochemistry()->component_intervals[i];
                    if (std::none_of(intervals.begin(), intervals.end(), [&](const auto& v) { return ideal_interval_contains(v, t); })) return false;
                }
            return true;
        };
        // #610: a mixture seed without a state starts at the first liquid-sided T_s r^j in the records, else at T_s.
        std::optional<CompiledSystem> cooled;
        try {
            const auto liquid_sided = [&](double t) {
                DensityRoot root;
                solve_density(*compiled.model, t, compiled.seed_pressure(), compiled.feed_composition, DensityBranch::Vapor, std::nullopt, &root);
                return liquid_feed(root);
            };
            double t = compiled.seed_temperature();
            while (!seed.has_state() && !pure && !liquid_sided(t) && inside(t)) t *= policy::COOLING_RATIO.value;
            seed_problem.temperature = FreeValue::fixed(t);
            if (t != compiled.seed_temperature() && inside(t)) cooled = compile_problem(*compiled.model, seed_problem);
        } catch (const std::exception& error) {
            if (stages) retain_density_failure(stages->first_numerical_failure, error, "cooled seed");
        }
        if (seed.has_state() || cooled) {
            const auto coordinates = cooled ? boundary_start(*cooled, true, limits, stages) : seed.coordinates;
            double T = cooled ? cooled->seed_temperature() : *seed.temperature;
            const auto prepared = (cooled ? *cooled : source).prepare(coordinates);
            point[compiled.primal_size - 1] = T;  // admits the records at the seed T, or names the refusal
            const auto destination = compiled.prepare(point);
            std::vector<double> values, capacities;
            for (std::size_t p = 0; p < prepared.shape.size(); ++p) {
                const auto& a = prepared.shape[p];
                const auto& b = destination.shape[p];
                for (std::size_t i = 0; i < a.x.size(); ++i) point[b.x[i]] = coordinates[a.x[i]];
                point[b.volume] = coordinates[a.volume];
                point[*b.amount] = prepared.amount[p];
                if (cooled) continue;
                const auto state = compiled.model->state(T, seed.molar_densities[p], prepared.composition[p]);
                if (destination.calorics) {  // anchored ideal plus the residual parts of the frozen phase
                    const auto ideal = anchored_ideal(*destination.calorics, compiled.problem.phases[p].support,
                        prepared.composition[p], 1.0 / seed.molar_densities[p], T, 0.0);
                    values.push_back(compiled.problem.enthalpy_j ? ideal.enthalpy + state.residual_enthalpy : ideal.entropy + state.residual_entropy);
                    capacities.push_back(ideal.isobaric_heat_capacity + state.residual_isobaric_heat_capacity.value_or(std::numeric_limits<double>::quiet_NaN()));
                    continue;
                }
                const auto& value = compiled.problem.enthalpy_j ? state.total_enthalpy : state.total_entropy;
                if (!value) throw std::invalid_argument("caloric seed state has no finite total caloric value");
                values.push_back(*value);
                capacities.push_back(state.total_isobaric_heat_capacity.value_or(std::numeric_limits<double>::quiet_NaN()));
            }
            const double target = compiled.problem.enthalpy_j ? *compiled.problem.enthalpy_j : *compiled.problem.entropy_j_per_k;
            if (pure) {
                const double split = (target / compiled.feed_total - values[0]) / (values[1] - values[0]);
                if (!(split > 0.0 && split < 1.0)) throw std::domain_error("caloric target outside declared coexistence interval");
                point[*destination.shape[0].amount] = compiled.feed_total * (1.0 - split);
                point[*destination.shape[1].amount] = compiled.feed_total * split;
            } else if (!cooled) {
                const double value = std::inner_product(values.begin(), values.end(), prepared.amount.begin(), 0.0);
                const double cp = std::inner_product(capacities.begin(), capacities.end(), prepared.amount.begin(), 0.0);
                const double next = compiled.problem.enthalpy_j ? T + (target - value) / cp : T * std::exp((target - value) / cp);
                if (cp > 0.0 && std::isfinite(next) && next > 0.0 && inside(next)) T = next;
            }
            point[compiled.primal_size - 1] = T;
            return point;
        }
    }
    auto start = true_species_split(compiled, scope, limits, stages);
    if (!start) start = start_compositions(compiled, scope);
    auto& compositions = start->compositions;
    const std::size_t free_phases = compiled.problem.phases.size() - compiled.pinned_count;
    const double feed_total = compiled.feed_total;
    std::optional<double> speciated_density;
    double free_amount = feed_total;  // ALG-EQUILIBRIUM default free amount (#407)
    scope.probe("free_phases", static_cast<double>(free_phases));
    scope.probe("feed_total", feed_total);
    for (std::size_t phase = 0; phase < compiled.problem.phases.size(); ++phase)
        compositions[phase] = seed_composition(compiled, phase, compositions[phase], point, phase_offsets[phase], start->declared_guesses);
    if (speciate && compiled.has_reaction && compiled.problem.phases.size() == 1 && free_phases == 1)
        if (const auto speciated = speciate_single_phase(compiled, point, compositions.front(), scope, 0, stages))
            std::tie(free_amount, speciated_density) = *speciated;
    bool corrected = false;
    const auto& problem = compiled.problem;
    correct_dense = correct_dense && compiled.has_reaction && !compiled.model->parameters().ideal_gas
        && problem.solids.empty() && !problem.discovery && !problem.enthalpy_j && !problem.entropy_j_per_k && !problem.volume_m3
        && !problem.azeotrope && !problem.critical && !problem.spinodal && !problem.envelope && !problem.loading
        && problem.phases.size() == 2 && free_phases == 1;
    const double p = !start->densities.empty() ? start->pressure : fugacity_pressure_seed(compiled, point, phase_offsets, *start, scope, corrected, correct_dense, limits, stages);
    if (compiled.problem.pressure.free) point[compiled.primal_size - 1] = p;
    compiled.pin_liquid_composition(point);
    for (std::size_t phase = 0; phase < compiled.problem.phases.size(); ++phase) {
        if (limits) limits->check();
        std::size_t cursor = phase_offsets[phase] + compiled.problem.phases[phase].support.size();
        if (compiled.problem.phases[phase].amount.kind == AmountKind::Free) {
            point[cursor++] = start->amounts.empty() ? free_amount : start->amounts[phase];
        }
        const double t = compiled.seed_temperature();
        std::vector<double> composition(compiled.model->parameters().component_count(), 0.0);
        for (std::size_t i = 0; i < compiled.problem.phases[phase].support.size(); ++i) {
            composition[compiled.problem.phases[phase].support[i]] =
                std::exp(point[phase_offsets[phase] + i]);
        }
        if ((supports_split_seed(compiled) || (phase < packed.size() && packed[phase]))
            && compiled.problem.phases[phase].packing_guess.has_value()) {
            const double packing = *compiled.problem.phases[phase].packing_guess;
            if (!admitted_packing_fraction(packing)) {
                throw std::invalid_argument("packing guess is outside the admitted range");
            }
            point[cursor] = std::log(packing);
            scope.probe("packing_guess", packing, static_cast<std::uint32_t>(phase));
            continue;
        }
        const auto density = start->densities.empty() ? speciated_density : std::optional(start->densities[phase]);
        if (density.has_value()) {
            point[cursor] = std::log((compiled.model->parameters().ideal_gas ? 1.0 : compiled.model->packing_fraction_coefficient(t, composition)) * *density);
            continue;
        }
        const DensityBranch branch = branch_of(compiled.problem.phases[phase].kind);
        DensityRoot root;
        try {
            const auto unanchored = [&] { return corrected && branch == DensityBranch::Vapor
                ? declared_root(*compiled.model, t, p, composition, compiled.problem.phases[phase].kind) : solve_density(*compiled.model, t, p, composition, branch, std::nullopt); };
            if (!corrected && start->vapor_density > 0.0 && branch == DensityBranch::Vapor) try {  // the phi^V estimate's root continues
                root = solve_density(*compiled.model, t, p, composition, branch, start->vapor_density);
            } catch (const DensityStructureError&) { throw;
            } catch (const std::exception& error) { if (stages) retain_density_failure(stages->first_numerical_failure, error, "boundary anchor"); root = unanchored(); }
            else root = unanchored();
        } catch (const std::runtime_error& error) {
            std::optional<core::Convergence> local_first;
            auto& first = stages ? stages->first_numerical_failure : local_first;
            retain_density_failure(first, error, "boundary density");
            throw DensityStructureError(std::string("boundary seed has no density root: ") + error.what(), first);
        }
        if (stages) { core::retain_first(stages->first_numerical_failure, root.first_numerical_failure);
            for (const auto& rank : root.rank_records) core::retain_rank(stages->rank_records, rank); }
        // Coexistence needs two stable pure roots only when more than one phase is declared;
        // a single declared pure phase uses the existing unique-root branch rules.
        if (compiled.model->parameters().component_count() == 1
            && compiled.problem.phases.size() > 1
            && (!root.stable_root_count.has_value() || *root.stable_root_count < 2)) {
            throw std::runtime_error(std::string("boundary seed has no coexisting density roots")
                + (std::isnan(start->pressure) ? " and no positive mechanical-loop witness pressure" : "")
                + " at T = " + core::optional_number(t) + " K; T may exceed the model critical temperature");
        }
        if (limits) limits->check();
        if (!(std::isfinite(root.molar_density) && root.molar_density > 0.0))
            throw std::runtime_error("boundary seed has invalid density");
        scope.probe("phase", static_cast<double>(phase));
        scope.probe("branch_density", root.molar_density,
                    static_cast<std::uint32_t>(phase));
        scope.probe("pressure", p);
        scope.probe("temperature", t);
        point[cursor] = std::log((compiled.model->parameters().ideal_gas ? 1.0 : compiled.model->packing_fraction_coefficient(t, composition))
                                 * root.molar_density);
    }
    return point;
}

// Every field that is not a path lever must be identical at both endpoints. Without `endpoint`, free T/P
// seeds and the labels that leave the state and its rows unchanged (phase names, row_tolerance, the
// support_declared and topology_declared flags) are not compared.
std::string declaration_difference(const EquilibriumProblem& a, const EquilibriumProblem& b, bool endpoint) {
    if (a.prescribed_liquid_composition != b.prescribed_liquid_composition) return "prescribed liquid composition differs";
    if (a.critical.has_value() != b.critical.has_value() || a.spinodal.has_value() != b.spinodal.has_value())
        return "null-mode declaration differs";
    if (a.spinodal && (a.spinodal->density_nodes != b.spinodal->density_nodes
        || a.spinodal->density_guess_mol_per_m3 != b.spinodal->density_guess_mol_per_m3)) return "spinodal search differs";
    if (!(a.source_reference == b.source_reference)) return "source reference declaration differs";
    if (a.internal_energy_j.has_value() != b.internal_energy_j.has_value() || a.volume_m3.has_value() != b.volume_m3.has_value()) return "volume condition kind differs";
    if (a.enthalpy_j.has_value() != b.enthalpy_j.has_value()
        || a.entropy_j_per_k.has_value() != b.entropy_j_per_k.has_value()) return "caloric condition kind differs";
    if (a.phases.empty() || b.phases.empty()) return "endpoint declares no phases";
    if (a.phases.size() != b.phases.size()) return "phase count differs";
    for (std::size_t p = 0; p < a.phases.size(); ++p) {
        const auto& x = a.phases[p];
        const auto& y = b.phases[p];
        if ((endpoint && x.name != y.name) || x.kind != y.kind) return "phase " + x.name + " identity differs";
        if (x.support != y.support || (endpoint && x.support_declared != y.support_declared))
            return "phase " + x.name + " support differs";
        if (x.amount.kind != y.amount.kind || (x.is_pinned() && x.amount.fraction != y.amount.fraction))
            return "phase " + x.name + " amount declaration differs";
    }
    if (a.reactions.size() != b.reactions.size()) return "reaction count differs";
    for (std::size_t r = 0; r < a.reactions.size(); ++r) {
        const auto& x = a.reactions[r];
        const auto& y = b.reactions[r];
        if (!same_reaction(x, y)) return "reaction " + reaction_name(x, r) + " differs";
    }
    if (a.standard_gibbs_j_per_mol != b.standard_gibbs_j_per_mol) return "standard Gibbs input differs";
    if (a.solids.size() != b.solids.size()) return "solid count differs";
    for (std::size_t s = 0; s < a.solids.size(); ++s)
        if (a.solids[s].stoichiometry != b.solids[s].stoichiometry || a.solids[s].amount.kind != b.solids[s].amount.kind)
            return "solid " + a.solids[s].name + " differs";
    if (endpoint && a.topology_declared != b.topology_declared) return "topology declaration differs";
    if (a.feed.kind != b.feed.kind || a.feed.solvent != b.feed.solvent
        || a.feed.solvent_mass_kg != b.feed.solvent_mass_kg || a.feed.values.size() != b.feed.values.size())
        return "feed basis, solvent set or component count differs";
    if (endpoint && a.solver_options.value(policy::row) != b.solver_options.value(policy::row)) return "row_tolerance differs";
    if (!same_intensive(a.temperature, b.temperature, endpoint)) return "temperature free/fixed declaration differs";
    if (!same_intensive(a.pressure, b.pressure, endpoint)) return "pressure free/fixed declaration differs";
    return "";
}

std::string declared_state_difference(const EquilibriumProblem& a, const EquilibriumProblem& b) {
    auto difference = declaration_difference(a, b, false);
    if (a.feed.values != b.feed.values) difference = "feed values differ";
    if (!b.temperature.free && a.temperature.value != b.temperature.value) difference = "fixed temperature differs";
    if (!b.pressure.free && a.pressure.value != b.pressure.value) difference = "fixed pressure differs";
    return difference;
}

namespace {

bool check_declared_side(const DensityRoot& root, DensityBranch branch, bool witness = false) {
    const auto side = properties::density_side(root.branch, root.molar_density, root.negative_slope_density);
    if (side == DensityBranch::Unspecified || side == branch) return true;
    if (witness && root.negative_slope_density) return false;
    const bool liquid = side == DensityBranch::Liquid;
    throw std::runtime_error(std::string(liquid ? "vapor" : "liquid") + " branch not admitted by the finite search: the sole stable root "
        + core::optional_number(root.molar_density) + " mol/m3 lies " + (liquid ? "above" : "below")
        + " the negative-slope sample " + core::optional_number(root.negative_slope_density) + " mol/m3");
}

}  // namespace

DensityRoot declared_root(const HelmholtzModel& model, double T, double P, const std::vector<double>& x,
    PhaseKind kind, bool* same_side) {
    const auto branch = branch_of(kind);
    auto root = solve_density(model, T, P, x, branch, std::nullopt);
    if (!check_declared_side(root, branch, same_side != nullptr)) *same_side = false;
    return root;
}

namespace {

// Corrected-phase branch identity: the corrected density and the selected stable
// root must share the anchor's window, with a positive corrected pressure slope. The physical
// state is recomputed from the coordinates, never taken from the supplied result's arrays.
BranchCheck density_branches_at(
    const CompiledSystem& compiled, const PhasePoint& prepared,
    const std::vector<double>& anchors, bool initial, core::Scope& scope, const DriverLimits* limits,
    const std::array<DensityRoot, 2>* roots, std::optional<BranchEndCondition>* witness = nullptr
) {
    const auto& model = *compiled.model;
    BranchCheck outcome;
    const auto packing_order = [&] {
        // Declared kinds: a liquid less densely packed than a vapor is the other boundary with its labels swapped (#599).
        for (std::size_t l = 0; initial && l < prepared.eta.size(); ++l)
            for (std::size_t v = 0; v < prepared.eta.size(); ++v)
                if (compiled.problem.phases[l].kind == PhaseKind::Liquid && compiled.problem.phases[v].kind == PhaseKind::Vapor
                    && !(prepared.eta[l] > prepared.eta[v])) {
                    outcome.message = "declared liquid " + compiled.problem.phases[l].name + " (packing fraction " + core::optional_number(prepared.eta[l])
                        + ") is not more densely packed than declared vapor " + compiled.problem.phases[v].name + " ("
                        + core::optional_number(prepared.eta[v]) + "): swapped phase labels";
                    return false;
                }
        return true;
    };
    if (witness && !packing_order()) return outcome;
    const bool pure = one_component_two_density_eligible(compiled);
    for (std::size_t index = 0; index < compiled.problem.phases.size(); ++index) {
        // For branch-end evidence, establish the incumbent before examining a refused vapor.
        const auto phase = witness && compiled.problem.phases[0].kind == PhaseKind::Vapor ? 1 - index : index;
        if (witness && core::CalculationWork::stop_at("branch-end density")) return outcome;
        if (limits) limits->check();
        const auto& spec = compiled.problem.phases[phase];
        std::vector<double> x = prepared.composition[phase];
        const double total = std::accumulate(x.begin(), x.end(), 0.0);
        for (double& value : x) value /= total;  // the driver's own normalisation of the log chart
        const double corrected = chart_density(model, prepared, phase);
        const auto branch = branch_of(spec.kind);
        const double slope = pure ? prepared.pure[phase].pressure_density_derivative
            : PressureDensityDerivative::at(model, prepared, phase, corrected, x);
        if (limits) limits->check();
        if (witness && spec.kind == PhaseKind::Vapor && !(std::isfinite(slope) && slope > 0.0)) {
            if (std::isfinite(slope)) *witness = BranchEndCondition::VaporMechanicalLimit;
            scope.probe("corrected_slope", slope, phase);
            return outcome;
        }
        double anchor = anchors.empty() ? 0.0 : anchors[phase];
        DensityRoot root;
        bool same_side = true;
        try {
            if (initial) {
                if (witness && !roots && core::CalculationWork::stop_at("branch-end unanchored density")) return outcome;
                root = roots ? (*roots)[spec.kind == PhaseKind::Vapor ? 1 : 0]
                    : declared_root(model, prepared.temperature, prepared.pressure, x, spec.kind,
                        witness && spec.kind == PhaseKind::Vapor ? &same_side : nullptr);
                if (roots) check_declared_side(root, branch);
                if (pure && root.stable_root_count.value_or(0) < 2) {
                    outcome.message = "one-component two-density requires two distinct stable density branches"; return outcome;
                }
                anchor = root.molar_density;
                if (limits) limits->check();
            }
            if (same_side && (!initial || !roots)) {
                if (witness && core::CalculationWork::stop_at("branch-end anchored density")) return outcome;
                root = solve_density(model, prepared.temperature, prepared.pressure, x, branch, anchor);
            }
            if (limits) limits->check();
        } catch (const std::exception& error) {
            outcome.message = "phase " + spec.name + " has no locally unique anchored root: " + error.what();
            return outcome;
        }
        const double displacement = std::abs(std::log(corrected / anchor));
        const double selected_displacement = std::abs(std::log(root.molar_density / anchor));
        scope.probe("selected_density", root.molar_density, static_cast<std::uint32_t>(phase));
        scope.probe("corrected_density", corrected, static_cast<std::uint32_t>(phase));
        scope.probe("anchor_density", anchor, static_cast<std::uint32_t>(phase));
        scope.probe("corrected_displacement", displacement, static_cast<std::uint32_t>(phase));
        scope.probe("corrected_slope", slope, static_cast<std::uint32_t>(phase));
        scope.probe("density_residual", root.convergence.residual, static_cast<std::uint32_t>(phase));
        const bool ok = std::isfinite(corrected) && corrected > 0.0
            && std::isfinite(slope) && slope > 0.0
            && root.convergence.residual <= root.convergence.criterion->value
            && displacement <= policy::DENSITY_ANCHOR_LOG_WINDOW.value
            && selected_displacement <= policy::DENSITY_ANCHOR_LOG_WINDOW.value;
        if (!ok) {
            outcome.terminal = true;
            outcome.message = "phase " + spec.name + " corrected density left the anchored branch window";
            return outcome;
        }
        if (!same_side) {
            *witness = BranchEndCondition::VaporKindRefused;
            scope.probe("branch_end_wrong_side", {root.molar_density, *root.negative_slope_density});
            return outcome;
        }
    }
    if (!witness && !packing_order()) return outcome;
    outcome.accepted = true;
    return outcome;
}

}  // namespace

BranchCheck check_density_branches(
    const CompiledSystem& compiled, const EquilibriumResult& result,
    const std::vector<double>& anchors, bool initial, core::Scope& scope, const DriverLimits* limits
) {
    if (core::CalculationWork::active_limit())
        scope.probe("density_admission_work", static_cast<double>(core::CalculationWork::current()->units()));
    if (core::CalculationWork::stop_at("density admission")) return {false, false, "work allowance exhausted"};
    const auto& model = *compiled.model;
    if (limits) limits->check();
    // Admission reads the driver's published physical point when it matches the current chart.
    const bool pure = one_component_two_density_eligible(compiled);
    const auto* supplied = result.has_state() && result.prepared_endpoint && same_result_origin(compiled, result)
        && core::bitwise_equal(result.prepared_endpoint->point, result.coordinates)
        && (!pure || !result.prepared_endpoint->pure.empty()) ? result.prepared_endpoint.get() : nullptr;
    const auto owned = supplied ? std::optional<PhasePoint>{} : std::make_optional(compiled.prepare(result.coordinates, true));
    const auto& prepared = supplied ? *supplied : *owned;
    if (initial && !limits && pure) {
        std::array<DensityRoot, 2> roots;
        bool complete = false;
        try {
            auto x = prepared.composition.front();
            const double total = std::accumulate(x.begin(), x.end(), 0.0);
            for (double& value : x) value /= total;
            roots[0] = solve_density(model, prepared.temperature, prepared.pressure, x,
                DensityBranch::Liquid, std::nullopt, &roots[1]);
            complete = roots[0].branch == DensityBranch::Liquid && roots[1].branch == DensityBranch::Vapor
                && roots[0].stable_root_count.value_or(0) >= 2 && roots[1].stable_root_count.value_or(0) >= 2
                && roots[0].molar_density > roots[1].molar_density && roots[1].molar_density > 0.0;
        } catch (const std::exception&) {
            // Ordinary admission owns the phase and first cause of any failed shared search.
        }
        if (complete) {
            const auto shared = density_branches_at(compiled, prepared, anchors, initial, scope, limits, &roots);
            if (shared.accepted) return shared;
        }
    }
    return density_branches_at(compiled, prepared, anchors, initial, scope, limits, nullptr);
}

const char* continuation_stop_name(ContinuationStop stop) noexcept {
    switch (stop) {
    case ContinuationStop::reached_target: return "reached_target";
    case ContinuationStop::initial_failure: return "initial_failure";
    case ContinuationStop::invalid_path: return "invalid_path";
    case ContinuationStop::singular_jacobian: return "singular_jacobian";
    case ContinuationStop::branch_lost: return "branch_lost";
    case ContinuationStop::phase_disappeared: return "phase_disappeared";
    case ContinuationStop::phases_coalesced: return "phases_coalesced";
    case ContinuationStop::left_domain: return "left_domain";
    case ContinuationStop::infeasible: return "infeasible";
    case ContinuationStop::nonfinite: return "nonfinite";
    case ContinuationStop::stalled: return "stalled";
    case ContinuationStop::step_too_small: return "step_too_small";
    case ContinuationStop::budget_exhausted: return "budget_exhausted";
    case ContinuationStop::critical_endpoint: return "critical_endpoint";
    case ContinuationStop::turning_point: return "turning_point";
    }
    return "unknown";
}

// `limits` (recovery, #565) also caps attempted steps at max_steps; a supplied start with `anchors` is admitted on them.
static ContinuationResult continue_path(const HelmholtzModel& model, const ContinuationRequest& request,
                                        const DriverLimits* limits, std::vector<double> anchors, StageWork* output = nullptr) {
    if (!request.start_problem.solids.empty() || !request.target_problem.solids.empty())
        throw std::invalid_argument("continuation is not admitted with solids");
    if (request.start_problem.envelope) return trace_envelope(model, request);
    ContinuationResult result;
    auto& observations = output ? *output : result.observations;
    const auto& options = request.options;
    if (!std::isfinite(options.initial_step) || options.initial_step <= 0.0 || options.initial_step > 1.0
        || !std::isfinite(options.minimum_step) || options.minimum_step <= 0.0
        || options.minimum_step > options.initial_step || options.max_steps == 0) {
        throw std::invalid_argument("invalid continuation options");
    }
    auto a = request.start_problem, b = request.target_problem;
    a.solver_options = a.solver_options.resolved(model.solver_options());
    b.solver_options = b.solver_options.resolved(model.solver_options());
    const std::string declared = declaration_difference(a, b);
    if (!declared.empty()) return finish(std::move(result), ContinuationStop::invalid_path, declared);
    CompiledSystem start = compile_problem(model, request.start_problem);
    const bool parameter_path = !request.parameter_targets.empty();
    if (start.problem.spinodal && (parameter_path || request.start_problem.feed.values != request.target_problem.feed.values))
        return finish(std::move(result), ContinuationStop::invalid_path, "spinodal continuation holds composition and parameters fixed");
    const auto parameter_start = model.active_parameter_values();
    if (parameter_path && (request.parameter_targets.size() != parameter_start.size()
        || !std::all_of(request.parameter_targets.begin(), request.parameter_targets.end(), [](double v) { return std::isfinite(v); })))
        return finish(std::move(result), ContinuationStop::invalid_path, "invalid physical parameter endpoints");
    const auto parameter_model = [&](double s) {
        if (!parameter_path) return std::unique_ptr<HelmholtzModel>{};
        auto parameters = model.parameters();
        for (std::size_t i = 0; i < parameter_start.size(); ++i) {
            const auto& coordinate = model.active_parameters()[i];
            for (const auto index : coordinate.indices)
                detail::family_values(parameters, coordinate.family).first[index] = (1.0 - s) * parameter_start[i] + s * request.parameter_targets[i];
        }
        return std::make_unique<HelmholtzModel>(std::move(parameters), model.active_parameters(), model.thermochemistry(), model.solver_options());
    };
    auto target_model = parameter_model(1.0);
    CompiledSystem target = compile_problem(target_model ? *target_model : model, request.target_problem);
    const std::string layout = layout_difference(start, target);
    if (!layout.empty()) return finish(std::move(result), ContinuationStop::invalid_path, layout);
    if (!start.feasible || !target.feasible) {
        core::Scope rank_scope("continuation_endpoint_rank");
        observations.observe(structural_result(start, rank_scope)); observations.observe(structural_result(target, rank_scope));
        return finish(std::move(result), ContinuationStop::invalid_path, "endpoint is not full rank at its generic point");
    }
    Path path;
    path.base = request.start_problem;
    path.feed_start = start.n_feed;
    path.feed_target = target.n_feed;
    path.t1 = request.target_problem.temperature.value;
    path.p1 = request.target_problem.pressure.value;
    path.u1 = request.target_problem.internal_energy_j.value_or(0.0);
    path.v1 = request.target_problem.volume_m3.value_or(0.0);
    path.h1 = request.target_problem.enthalpy_j.value_or(0.0);
    path.s1 = request.target_problem.entropy_j_per_k.value_or(0.0);
    const bool moves = request.start_problem.internal_energy_j != request.target_problem.internal_energy_j || request.start_problem.volume_m3 != request.target_problem.volume_m3
        || request.start_problem.enthalpy_j != request.target_problem.enthalpy_j
        || request.start_problem.entropy_j_per_k != request.target_problem.entropy_j_per_k || path.feed_start != path.feed_target
        || (!path.base.temperature.free && path.base.temperature.value != path.t1)
        || (!path.base.pressure.free && path.base.pressure.value != path.p1);
    if (parameter_path && moves)
        return finish(std::move(result), ContinuationStop::invalid_path, "parameter path must hold intensives and resolved feed fixed");
    if (!moves && (!parameter_path || parameter_start == request.parameter_targets))
        return finish(std::move(result), ContinuationStop::invalid_path,
            parameter_path ? "no active parameter differs" : "no fixed intensive or resolved feed differs");
    if (parameter_path) { path.parameter_start = parameter_start; path.parameter_target = request.parameter_targets; }

    // Initial admission at q = 0.
    std::vector<double> y_current;
    CompiledSystem current;
    std::unique_ptr<HelmholtzModel> current_model;
    {
        core::Scope initial("equilibrium_continuation_initial");
        initial.probe("path_value", 0.0);
        BoundaryAttempt admitted;
        const auto solve_initial = [&](const std::vector<double>& coordinates) {
            return unless_failed(solve_equilibrium_from_coordinates(start, coordinates), &observations);
        };
        if (request.start_result.has_value()) {
            const auto& supplied = *request.start_result;
            const bool layout_ok = supplied.coordinates.size() == start.unknown_count()
                && supplied.phase_count == start.problem.phases.size()
                && supplied.mole_fractions.size() == supplied.phase_count * model.parameters().component_count()
                && supplied.rows.size() == start.rows.size();
            if (!layout_ok) {
                initial.probe("accepted", 0.0);
                initial.probe(continuation_stop_name(ContinuationStop::invalid_path), 1.0);
                return finish(std::move(result), ContinuationStop::invalid_path, "supplied start result layout differs");
            }
            admitted = admit_result(start, parameter_path
                ? solve_initial(supplied.coordinates) : supplied, anchors, anchors.empty(), initial,
                !parameter_path);
        } else {
            EquilibriumResult seed;
            try { seed = solve_initial(boundary_start(start, true, nullptr, &observations)); }
            catch (const std::bad_alloc&) { throw; }
            catch (const std::exception& error) {
                if (!start.problem.spinodal || !start.spinodal_diagnostics() || std::isfinite(start.spinodal_diagnostics()->seed)) throw;
                seed = spinodal_search_failure(start, structural_result(start, initial), error.what());
            }
            admitted = admit_result(start, std::move(seed), anchors, true, initial);
        }
        if (!output || !request.start_result || parameter_path) observations.observe(admitted.observations);
        else core::retain_first(observations.first_numerical_failure, admitted.observations.first_numerical_failure);
        initial.probe("accepted", admitted.accepted ? 1.0 : 0.0);
        initial.probe(continuation_stop_name(admitted.accepted ? ContinuationStop::reached_target : admitted.stop), 1.0);
        if (admitted.convergence.has_value()) initial.finish(*admitted.convergence);
        if (!admitted.accepted) {
            reject(result, admitted, 0.0);
            return finish(std::move(result), ContinuationStop::initial_failure, admitted.message);
        }
        result.path_values.push_back(0.0);
        if (parameter_path) result.parameter_values.push_back(parameter_start);
        result.states.push_back(*admitted.result);
        y_current = admitted.result->coordinates;
        anchors = admitted.result->molar_densities;
        current = start;
    }

    double q = 0.0;
    double dq = options.initial_step;
    std::size_t local_retries = 0;
    bool transfer = true;
    double w = NAN, w_before = NAN, q_before = NAN, q_failed = 1.0;  // turning-point record
    while (true) {
        if (q >= 1.0) return finish(std::move(result), ContinuationStop::reached_target, "");
        if (result.accepted_steps >= options.max_steps)
            return finish(std::move(result), ContinuationStop::budget_exhausted, "max_steps reached before the target");
        if (limits && (limits->expired() || result.attempts >= options.max_steps))
            return finish(std::move(result), ContinuationStop::budget_exhausted, "recovery work or attempts exhausted");
        const double q_try = std::min(q + dq, 1.0);
        ++result.attempts;
        core::Scope attempt_scope("equilibrium_continuation_attempt");
        attempt_scope.probe("path_value", q_try);
        const EquilibriumProblem problem = path.at(q_try);
        auto trial_model = parameter_model(q_try);
        const auto& trial = trial_model ? *trial_model : model;
        CompiledSystem compiled = compile_problem(trial, problem);
        BoundaryAttempt attempt = advance(start, current, y_current, compiled, problem, [&](auto* failure) {
            auto tangent = sensitivity(current, y_current, path.lever(current), failure);
            w = 1.0 / (1.0 + std::inner_product(tangent.begin(), tangent.end(), tangent.begin(), 0.0));
            return tangent;
        }, q_try - q, anchors, parameter_path ? DensityAdmission::Predicted : DensityAdmission::Previous, attempt_scope,
            transfer ? result.states.back().optimizer_state : nullptr, limits, &observations,
            result.states.back().spinodal ? &result.states.back().spinodal->search : nullptr);
        observations.observe(attempt.observations);
        if (attempt.accepted)
            if (const auto pair = swapped_pair(compiled, y_current, attempt.result->coordinates)) coalesced(attempt, compiled, *pair);
        if (!attempt.result.has_value()) attempt_scope.probe("driver_result_present", 0.0);
        attempt_scope.probe("accepted", attempt.accepted ? 1.0 : 0.0);
        attempt_scope.probe(continuation_stop_name(attempt.accepted ? ContinuationStop::reached_target : attempt.stop), 1.0);
        if (attempt.convergence.has_value()) attempt_scope.finish(*attempt.convergence);
        if (attempt.accepted) {
            if (current.problem.spinodal) {
                const auto first = current.mode_offset();
                double overlap = 0;
                for (std::size_t i = 0; i < current.mode_count; ++i) overlap += y_current[first + i] * attempt.result->coordinates[first + i];
                if (overlap < 0) for (std::size_t i = 0; i < current.mode_count; ++i) {
                    attempt.result->coordinates[first + i] *= -1;
                    if (attempt.result->spinodal) attempt.result->spinodal->mode[i] *= -1;
                    attempt.result->residuals[current.mode_row()+i] *= -1;
                    attempt.result->classified_residuals[current.mode_row()+i] *= -1;
                    attempt.result->prepared_endpoint.reset();
                    attempt.result->optimizer_state.reset();
                }
            }
            w_before = w;
            q_before = q;
            q = q_try;
            ++result.accepted_steps;
            local_retries = 0;
            transfer = true;
            result.path_values.push_back(q);
            y_current = attempt.result->coordinates;
            anchors = attempt.result->molar_densities;
            if (parameter_path) result.parameter_values.push_back(trial.active_parameter_values());
            current_model = std::move(trial_model);
            current = std::move(compiled);
            result.reached_path = q;
            result.states.push_back(std::move(*attempt.result));
            continue;
        }
        transfer = false;
        core::retain_failure(observations.first_numerical_failure, attempt.message, "continuation admission");
        reject(result, attempt, q_try);
        if (attempt.terminal) return finish(std::move(result), attempt.stop, attempt.message);
        if (local_retries == 0) q_failed = q_try;
        const double halved = dq / 2.0;
        if (halved < options.minimum_step || local_retries >= options.max_retries_per_step) {
            // ALG-EQUILIBRIUM turning point: w = 1/(1 + |dy/ds|^2) of the last two accepted states vanishes by the first failed s.
            const double fold = q + (q - q_before) * w / (w_before - w);
            const auto stop = halved < options.minimum_step ? ContinuationStop::step_too_small : ContinuationStop::budget_exhausted;
            if (!(w_before > w && w > 0.0 && fold <= q_failed)) return finish(std::move(result), stop, attempt.message);
            return finish(std::move(result), ContinuationStop::turning_point, "path component of the tangent vanishes at s "
                + core::optional_number(fold) + "; " + continuation_stop_name(stop) + ": " + attempt.message);
        }
        dq = halved;
        ++local_retries;
        ++result.retries;
    }
}

ContinuationResult continue_equilibrium(const HelmholtzModel& model, const ContinuationRequest& request) {
    core::Scope scope("continue_equilibrium");
    if (request.start_problem.prescribed_liquid_composition || request.target_problem.prescribed_liquid_composition)
        throw std::invalid_argument("continuation unavailable for prescribed liquid composition");
    return continue_path(model, request, nullptr, {});
}

bool supports_work_limit(const EquilibriumProblem& p) {
    if (p.phase_calculation != PhaseCalculation::Declared || p.phases.size() != 2 || p.temperature.free || !p.pressure.free || p.feed.kind != FeedKind::Amounts
        || p.discovery || p.critical || p.spinodal || p.envelope || p.azeotrope || p.loading
        || !p.solids.empty() || p.volume_m3 || p.internal_energy_j || p.enthalpy_j || p.entropy_j_per_k)
        return false;
    const auto& liquid = p.phases[p.phases[0].kind == PhaseKind::Liquid ? 0 : 1];
    const auto& vapor = p.phases[p.phases[0].kind == PhaseKind::Liquid ? 1 : 0];
    return liquid.kind == PhaseKind::Liquid && liquid.is_free() && vapor.kind == PhaseKind::Vapor && vapor.is_zero_pinned();
}

namespace {

struct BubbleBranchSample final {
    std::vector<double> coordinates, tangent;
    double log_density;
    std::optional<BranchEndCondition> refusal;
};

std::optional<BubbleBranchSample> correct_branch_density(
    const CompiledSystem& c, std::vector<double> y, std::size_t v, double lambda,
    core::Scope& scope, bool admitted_seed = false) try {
    const auto tolerance = c.problem.solver_options.value(policy::row);
    const auto admissible = [&](auto& point) {
        std::size_t projected = 0;
        return inside_domain(c, point, projected, true) && !projected;
    };
    const auto evaluate = [&](const auto& point) {
        auto prepared = c.prepare(point, true);
        auto rows = c.evaluate_at(prepared).values;
        rows.back() = std::log(chart_density(*c.model, prepared, v)) - lambda;
        return std::pair{std::move(prepared), std::move(rows)};
    };
    if (core::CalculationWork::stop_at("branch-end preparation") || !admissible(y)) return {};
    auto [at, f] = evaluate(y);
    const auto predicted = at.volume;
    while (true) {
        const auto norm = maximum_abs(f);
        if (!std::isfinite(norm) || core::CalculationWork::stop_at("branch-end Jacobian")) return {};
        auto jac = c.jacobian_at(at);
        auto row = jac.end() - at.point.size();
        std::fill(row, jac.end(), 0.0);
        row[at.shape[v].volume] = 1.0;
        for (std::size_t j = 0; j < at.shape[v].x.size(); ++j)
            row[at.shape[v].x[j]] = -at.packing[v].component[c.problem.phases[v].support[j]] / at.packing[v].coefficient;
        const auto error = std::max(norm, maximum_abs(c.classified_residuals(at, f)));
        const bool corrected = error <= tolerance;
        if (corrected) { f.assign(at.point.size(), 0.0); f.back() = -1.0; }
        else if (admitted_seed) return {};
        auto direction = solve_implicit(c, at, f, nullptr, std::move(jac));
        if (direction.empty()) return {};
        if (corrected) {
            if (check_topology(c, at.point).event != TopologyEvent::None) return {};
            for (std::size_t phase = 0; phase < 2; ++phase)
                if (std::abs(std::log(predicted[phase] / at.volume[phase])) > policy::DENSITY_ANCHOR_LOG_WINDOW.value) return {};
            std::optional<BranchEndCondition> witness;
            if (!admitted_seed && !density_branches_at(c, at, {}, true, scope, nullptr, nullptr, &witness).accepted && !witness) return {};
            scope.probe("branch_end_sample", {at.point[c.loading_column()], at.pressure, chart_density(*c.model, at, v),
                direction[c.loading_column()], direction[c.primal_size - 1], error});
            return BubbleBranchSample{std::move(at.point), std::move(direction), lambda, witness};
        }
        const auto trial = [&](double damping, std::uint32_t, const auto& coordinates, const auto& contracts) {
            if (core::CalculationWork::stop_at("branch-end damping")) return NewtonTrial::Stop;
            auto candidate = coordinates(damping);
            if (!admissible(candidate) || candidate == at.point) return NewtonTrial::Refused;
            try {
                auto [prepared, rows] = evaluate(candidate);
                if (!contracts(maximum_abs(rows))) return NewtonTrial::Refused;
                at = std::move(prepared);
                f = std::move(rows);
                return NewtonTrial::Accepted;
            } catch (const std::bad_alloc&) { throw; }
            catch (const std::exception&) { return NewtonTrial::Refused; }
        };
        if (damped_newton_step(at.point, direction, norm, trial) != NewtonTrial::Accepted) return {};
    }
} catch (const std::bad_alloc&) { throw; }
catch (const std::exception&) { return {}; }

std::optional<BubbleBranchEnd> classify_branch_end(
    const CompiledSystem& source, const CompiledSystem& target, const EquilibriumResult& seed,
    std::size_t component, std::size_t v, core::Scope& scope) {
    const auto q0 = source.n_feed[component], q1 = target.n_feed[component];
    if (source.model->parameters().ideal_gas || !(q0 > 0)) return {};
    auto request = source.problem;
    request.loading = AffineLoading{source.n_feed, std::vector<double>(source.n_feed.size()), q0, q1, q0, *seed.pressure};
    request.loading->initial[component] = 0.0;
    request.loading->direction[component] = 1.0;
    const auto c = compile_problem(*source.model, request);
    const auto column = c.loading_column();
    auto y = seed.coordinates;
    y.insert(y.begin() + column, q0);
    auto left = correct_branch_density(c, std::move(y), v, std::log(seed.molar_densities[v]), scope, true);
    if (!left || left->tangent[column] == 0) return {};
    const auto sign = std::copysign(1.0, left->tangent[column]);
    const auto resolution = policy::DENSITY_ANCHOR_LOG_WINDOW.value / 4;
    const double step = sign * 0.05; // r1 MEA trace spacing.
    std::optional<double> upper;
    while (!core::CalculationWork::stop_at("branch-end bracket")) {
        const auto delta = upper ? (*upper - left->log_density) / 2 : step;
        const auto lambda = left->log_density + delta;
        if (lambda == left->log_density || (upper && lambda == *upper)) return {};
        auto predicted = shifted_coordinates(left->coordinates, left->tangent, delta);
        predicted[column] = std::clamp(predicted[column], q0, q1);
        auto right = correct_branch_density(c, std::move(predicted), v, lambda, scope);
        if (right && right->coordinates[column] < q1 && !swapped_pair(c, left->coordinates, right->coordinates)) {
            const auto directed = sign * right->tangent[column];
            const bool admitted = !right->refusal;
            if (admitted && directed > 0) {
                left = std::move(right);
                continue;
            }
            if (std::abs(delta) <= resolution && (admitted ? directed < 0 : directed > 0)) {
                const auto bracket = std::pair{left->coordinates[column], right->coordinates[column]};
                const auto& estimate = admitted && bracket.second > bracket.first ? *right : *left;
                const auto variation = std::abs(delta) * (std::abs(left->tangent[column]) + std::abs(right->tangent[column]));
                if (admitted && !(estimate.coordinates[column] + variation < q1)) return {};
                auto feed = source.n_feed;
                feed[component] = estimate.coordinates[column];
                return BubbleBranchEnd{right->refusal.value_or(BranchEndCondition::LoadingFold), bracket,
                    std::move(feed), estimate.coordinates[column + 1]};
            }
        } else if (std::abs(delta) <= resolution) return {};
        upper = lambda;
    }
    return {};
}

}  // namespace

EquilibriumResult solve_from_feed_increment(const CompiledSystem& target, const EquilibriumResult& seed, core::Scope& scope,
    std::optional<std::size_t> branch_end_work_limit) {
    if (seed.model_identity != target.model->identity() || !seed.declaration || !seed.has_state()
        || !seed.requested_tolerance_met.value_or(false) || !supports_work_limit(*seed.declaration))
        throw std::invalid_argument("bubble source must be admitted on the target model");
    if (const auto difference = declaration_difference(*seed.declaration, target.problem, true); !difference.empty())
        throw std::invalid_argument("bubble source: " + difference);
    if (seed.declaration->temperature.value != target.problem.temperature.value)
        throw std::invalid_argument("bubble source temperature differs");
    const auto source = compile_problem(*target.model, *seed.declaration);
    if (const auto difference = layout_difference(source, target); !difference.empty())
        throw std::invalid_argument("bubble source: " + difference);
    if (seed.coordinates.size() != source.unknown_count())
        throw std::invalid_argument("bubble source coordinates differ from compiled layout");
    const auto& charge = target.model->parameters().charge_number;
    std::size_t additions = 0, component = 0;
    for (std::size_t i = 0; i < target.n_feed.size(); ++i) {
        if (target.n_feed[i] == source.n_feed[i]) continue;
        if (target.n_feed[i] < source.n_feed[i] || (!charge.empty() && charge[i] != 0))
            throw std::invalid_argument("bubble source differs by more than a neutral feed addition");
        component = i;
        ++additions;
    }
    if (additions != 1) throw std::invalid_argument("bubble source needs one increased neutral feed amount");
    if (core::CalculationWork::stop_at("source admission")) return structural_result(target, scope);
    const auto admitted = admit_result(source, seed, {}, true, scope, true);
    if (!admitted.accepted) {
        auto result = structural_result(target, scope);
        result.message = "bubble source: " + admitted.message;
        result.first_numerical_failure = admitted.observations.first_numerical_failure;
        return result;
    }
    if (core::CalculationWork::stop_at("feed prediction")) return structural_result(target, scope);
    Path path{source.problem, source.n_feed, target.n_feed, {}, {}};
    path.t1 = target.problem.temperature.value;
    auto attempt = advance(source, source, seed.coordinates, target, target.problem, [&](auto* failure) {
        return sensitivity(source, seed.coordinates, path.lever(source), failure);
    }, 1.0, {}, DensityAdmission::Declared, scope, seed.optimizer_state);
    auto result = attempt.result ? std::move(*attempt.result) : structural_result(target, scope);
    if (attempt.accepted) {
        result = refine_selected_endpoint(target, std::move(result), scope);
        if (!branch_end_work_limit || (result.has_state() && !core::CalculationWork::stopped()
            && core::CalculationWork::current()->units() <= core::CalculationWork::active_limit()->units)) return result;
    } else {
        result.message = std::move(attempt.message);
        core::retain_first(result.first_numerical_failure, attempt.observations.first_numerical_failure);
    }
    if (branch_end_work_limit && attempt.result && result.topology_event == TopologyEvent::None) {
        scope.probe("branch_end_entry", static_cast<double>(core::CalculationWork::current()->units()));
        core::CalculationWork::Limit limit{*branch_end_work_limit};
        core::CalculationWork work(true, &limit);
        const auto vapor = source.problem.phases.front().kind == PhaseKind::Vapor ? 0 : 1;
        try { result.branch_end = classify_branch_end(source, target, seed, component, vapor, scope); }
        catch (const std::bad_alloc&) { throw; } catch (const std::exception&) {}
        if (work.units() > limit.units) core::CalculationWork::stop_at("branch-end return");
        if (limit.boundary) result.branch_end.reset();
        if (const auto& end = result.branch_end) {
            result.topology_event = TopologyEvent::BranchEnd;
            result.topology_phase = vapor;
            result.solver_status = "BRANCH_END";
            const char* condition = end->condition == BranchEndCondition::LoadingFold ? "loading fold"
                : end->condition == BranchEndCondition::VaporMechanicalLimit ? "vapor mechanical limit" : "loss of declared vapor admission";
            result.message.insert(0, std::string("branch end: ") + condition + " at "
                + core::optional_number(end->limiting_feed_mol[component]) + " mol of increased feed and "
                + core::optional_number(end->limiting_pressure_pa) + " Pa" + (result.message.empty() ? "" : "; "));
        }
        scope.probe("branch_end_work", {static_cast<double>(work.units()), static_cast<double>(limit.units)}, 0,
            limit.boundary ? limit.boundary : "");
    }
    refuse(result);
    return result;
}

// ALG-EQUILIBRIUM: pure saturation uses four neighboring sources and direct target transport; mixtures use
// the first admitted source at T_seed r^j, then fixed-T continuation, or for a
// free-T boundary a same-source chart correction at the source pressure and fixed-P continuation. Internal
// declarations carry no guesses; boundary endpoints share one free-intensive seed. No second source after a path stop.
EquilibriumResult recover_from_source(const CompiledSystem& target, const DriverLimits& limits) {
    core::Scope scope("equilibrium_recovery");
    const auto& model = *target.model;
    const double tolerance = target.problem.solver_options.value(policy::row);
    const bool free_temperature = target.problem.temperature.free;
    EquilibriumResult refused, found;
    StageWork observations;
    const auto publish = [&](EquilibriumResult result) { std::move(observations).publish_observations(result); return result; };
    auto declared = target.problem;
    declared.phase_calculation = PhaseCalculation::Discovery;
    for (auto& phase : declared.phases) { phase.composition_guess.clear(); phase.packing_guess.reset(); }
    if (declared.pressure.free) declared.pressure = FreeValue::variable();
    auto source = declared;
    if (free_temperature) source.pressure = FreeValue::variable();
    try {
        if (one_component_two_density_eligible(target)) {
            for (double ratio : {policy::COOLING_RATIO.value, 1.0 / policy::COOLING_RATIO.value, 1.25, 1.5}) {
                limits.check();
                source.temperature = FreeValue::fixed(target.seed_temperature() * ratio);
                const auto failed = [&](const EquilibriumResult& solved) {
                    observations.observe(solved);
                    with_failure(solved, [&](const std::exception& error) { retain_density_failure(observations.first_numerical_failure, error, "saturation source"); });
                    refused.message = solved.message;
                };
                std::optional<CompiledSystem> compiled;
                try { compiled.emplace(compile_problem(model, source)); }
                catch (const std::bad_alloc&) { throw; }
                catch (const std::exception& error) {  // a source that does not compile is skipped like one that fails
                    ++observations.stage;
                    retain_density_failure(observations.first_numerical_failure, error, "saturation source");
                    refused.message = error.what();
                    continue;
                }
                found = try_solve_compiled(*compiled, std::nullopt, nullptr, &limits);
                if (found.failure) { failed(found); continue; }
                observations.observe(found);
                scope.probe("recovery_source", {source.temperature.value}, 0, found.message.c_str());
                if (!found.success) continue;
                auto reached = try_solve_compiled(target, found, nullptr, &limits);
                if (reached.failure) { failed(reached); continue; }
                observations.observe(reached);
                if (reached.success) return publish(std::move(reached));
                refused.message = reached.message;
            }
            if (refused.message.empty()) refused.message = "no admitted saturation source: " + found.message;
            return publish(std::move(refused));
        }
        for (int j = 1; j <= static_cast<int>(policy::EQUILIBRIUM_RECOVERY_SOURCES.value) && !found.success; ++j) {
            limits.check();  // the slot is spent before compilation
            source.temperature = FreeValue::fixed(target.seed_temperature() * std::pow(policy::COOLING_RATIO.value, j));
            found = unless_failed(try_solve_compiled(compile_problem(model, source), std::nullopt, nullptr, &limits), &observations);
            observations.observe(found);
            scope.probe("recovery_source", {source.temperature.value}, static_cast<std::uint32_t>(j), found.message.c_str());
        }
        if (!found.success) { refused.message = "no admitted lower-temperature source: " + found.message; return publish(std::move(refused)); }
        const auto work = static_cast<std::size_t>(policy::EQUILIBRIUM_RECOVERY_WORK.value);
        ContinuationRequest request{source, declared, found, {0.05, 1.0e-6, work, 10}, {}};
        if (free_temperature) {
            request.start_problem = declared;
            request.start_problem.pressure = FreeValue::fixed(*found.pressure);
            auto chart = found.coordinates;
            chart[target.primal_size - 1] = *found.temperature;
            auto corrected = unless_failed(solve_equilibrium_from_coordinates(model, request.start_problem, chart, nullptr, &limits), &observations);
            observations.observe(corrected);
            double drift = corrected.success ? std::max(std::abs(*corrected.temperature / *found.temperature - 1.0),
                largest_difference(corrected.mole_fractions, found.mole_fractions)) : INFINITY;
            for (std::size_t p = 0; corrected.success && p < found.amounts.size(); ++p)
                drift = std::max({drift, std::abs(corrected.amounts[p] - found.amounts[p]) / target.feed_total,
                    std::abs(corrected.molar_densities[p] / found.molar_densities[p] - 1.0)});
            if (!(drift <= tolerance)) {
                refused.message = "free-temperature chart correction left the source state: " + corrected.message;
                core::retain_failure(observations.first_numerical_failure, refused.message, "source recovery admission");
                return publish(std::move(refused));
            }
            request.start_result = std::move(corrected);
        }
        auto path = [&] { core::Scope span("continue_equilibrium"); return continue_path(model, request, &limits, found.molar_densities, &observations); }();
        if (!path.success) {
            refused.message = "from source at " + core::optional_number(*found.temperature) + " K: " + path.message;
            return publish(std::move(refused));
        }
        auto reached = std::move(path.states.back());
        const auto admitted = check_density_branches(target, reached, {}, true, scope, &limits);
        if (!admitted.accepted || !(target.classified_residual_norm(reached.coordinates) <= tolerance)) {
            refused.message = "reached target refused: " + admitted.message;
                core::retain_failure(observations.first_numerical_failure, refused.message, "source recovery admission");
            return publish(std::move(refused));
        }
        reached.declaration = std::make_shared<const EquilibriumProblem>(target.problem);
        return publish(std::move(reached));
    } catch (const std::bad_alloc&) { throw; }
    catch (const std::exception& error) {
        retain_density_failure(observations.first_numerical_failure, error, "source recovery");
        refused.message = error.what();
        return publish(std::move(refused));
    }
}

EquilibriumResult recover_from_source(const CompiledSystem& target, EquilibriumResult failed) {
    if (failed.success || !one_component_two_density_eligible(target) || std::getenv("EPCSAFT_IPOPT_OPTIONS_FILE")) return failed;
    const auto cap = static_cast<std::size_t>(policy::EQUILIBRIUM_RECOVERY_WORK.value);
    auto work = DriverLimits::recovery_allowance(std::chrono::steady_clock::time_point::max(), cap);
    auto recovered = recover_from_source(target, work);
    StageWork observations; observations.observe(failed); observations.observe(recovered);
    std::move(observations).publish_observations(recovered.success ? recovered : failed);
    failed.iterations += work.consumed; failed.evaluations += work.evaluations; failed.hessian_calls += work.hessian_calls;
    if (!recovered.success) { failed.message += "; recovery: " + recovered.message; return failed; }
    recovered.message += "; recovered after: " + failed.message;
    recovered.rejected_reduction = failed.rejected_reduction;
    recovered.iterations = failed.iterations; recovered.evaluations = failed.evaluations; recovered.hessian_calls = failed.hessian_calls;
    return recovered;
}


void compare_loading_candidates(SaturationSearch& out, ComparisonEligibility eligibility) {
    const auto usable = [](const NumericalEstimate& e) {
        return e.available && e.resolved && std::isfinite(e.absolute) && e.absolute >= 0
            && std::isfinite(e.chart) && e.chart >= 0;
    };
    for (std::size_t i = 0; i < out.candidates.size(); ++i) {
        if (!out.candidates[i].admitted) continue;
        for (std::size_t j = 0; j < i; ++j) {
            if (!out.candidates[j].admitted) continue;
            SaturationComparison comparison{i, j, {}, "duplicate"};
            const auto& a = out.candidates[i].operands;
            const auto& b = out.candidates[j].operands;
            if (!eligibility.valid() || a.empty() || a.size() != b.size()) comparison.identity = "unresolved";
            for (std::size_t k = 0; k < a.size() && k < b.size(); ++k) {
                auto e = eligibility;
                e.admission_a = e.admission_a && usable(a[k].error) && a[k].quantity == b[k].quantity;
                e.admission_b = e.admission_b && usable(b[k].error);
                const auto c = compare_equilibrium_quantity(a[k].quantity, a[k].value, b[k].value, 0,
                    a[k].error, b[k].error, e, ComparisonMode::EstimatedError);
                comparison.operands.push_back(c);
                if (c.judgment == "failed") {
                    if (e.valid() && (c.cause == "outside_estimates" || c.cause == "ordinary_achievable")) {
                        if (comparison.identity != "unresolved") comparison.identity = "distinct";
                    } else comparison.identity = "unresolved";
                }
            }
            out.comparisons.push_back(std::move(comparison));
        }
        std::vector<std::size_t> duplicates;
        for (std::size_t g = 0; g < out.groups.size(); ++g) {
            std::size_t count = 0;
            for (auto member : out.groups[g]) for (const auto& c : out.comparisons)
                if (c.left == i && c.right == member) count += c.identity == "duplicate";
            if (count && count != out.groups[g].size()) out.causes.push_back("root_comparison_unresolved");
            if (count) duplicates.push_back(g);
        }
        if (duplicates.empty()) out.groups.push_back({i});
        else if (duplicates.size() == 1) out.groups[duplicates[0]].push_back(i);
        else out.causes.push_back("root_comparison_unresolved");
    }
    for (const auto& c : out.comparisons)
        if (c.identity == "unresolved") out.causes.push_back("root_comparison_unresolved");
    if (out.groups.size() > 1 || out.brackets > 1) out.causes.push_back("multiple_roots");
}

static unsigned loading_slope_reversals(const SaturationNode& node, const SaturationNode* previous) {
    if (!std::isfinite(node.mismatch)) return 0;
    const auto sign = [](double value) { return (value > 0) - (value < 0); };
    unsigned count = !std::isfinite(node.slope) || node.slope == 0;
    if (previous) {
        const double delta = node.mismatch - previous->mismatch;
        if (sign(node.slope) != sign(previous->slope) || (std::abs(delta) > node.estimates.pressure.absolute + previous->estimates.pressure.absolute
            && sign(delta) != sign(node.slope))) ++count;
    }
    return count;
}

void admit_loading_search(SaturationSearch& out) {
    compare_loading_candidates(out, {});
    for (std::size_t i = 0; i < out.nodes.size(); ++i)
        out.causes.insert(out.causes.end(), loading_slope_reversals(out.nodes[i], i ? &out.nodes[i - 1] : nullptr), "unresolved_slope_reversal");

    if (out.groups.empty() && out.causes.empty()) out.causes.push_back("no_resolved_root");
    if (out.complete && out.causes.empty() && out.groups.size() == 1) out.accepted = out.groups[0][0];
}

SaturationSearch search_loading(const HelmholtzModel& model, EquilibriumProblem declaration, const AffineLoading& affine,
    const std::vector<double>& grid, const DriverLimits* supplied, std::chrono::steady_clock::time_point deadline, const SaturationSearch* previous) {
    using Clock = std::chrono::steady_clock;
    if (grid.size() < 2 || grid.front() != affine.lower || grid.back() != affine.upper)
        throw std::invalid_argument("saturation grid must include both declared bounds");
    for (std::size_t i = 0; i < grid.size(); ++i)
        if (!std::isfinite(grid[i]) || grid[i] <= 0 || (i && grid[i] <= grid[i - 1]))
            throw std::invalid_argument("saturation grid must be positive and strictly increasing");
    const auto cap = static_cast<std::size_t>(policy::EQUILIBRIUM_RECOVERY_WORK.value);
    auto own = DriverLimits::recovery_allowance(deadline, cap);
    own.attempt_cap = cap;
    const auto& work = supplied ? *supplied : own;
    SaturationSearch out; out.grid = grid; out.declaration = affine;
    out.work_at_entry = {work.consumed, work.evaluations, work.hessian_calls, work.attempts};
    const auto fail = [&](const std::string& cause) { out.causes.push_back(cause); };
    const auto fail_at = [&](double alpha, const std::string& cause) { fail("loading " + core::optional_number(alpha) + ": " + cause); };
    const auto solve = [&](double loading, double lower, double upper, double target,
                           const std::optional<EquilibriumResult>& start, bool composed) {
        if (work.attempt_cap && work.attempts >= work.attempt_cap) work.exhausted = true;
        work.check();
        auto problem = declaration;
        problem.feed.kind = FeedKind::Amounts;
        problem.feed.values = affine.initial;
        for (std::size_t i = 0; i < problem.feed.values.size(); ++i)
            problem.feed.values[i] += loading * affine.direction[i];
        problem.loading.reset();
        if (composed) { problem.loading = affine; problem.loading->seed = loading;
            problem.loading->lower = lower; problem.loading->upper = upper; problem.loading->pressure = target; }
        auto compiled = compile_problem(model, problem);
        auto state = try_solve_compiled(compiled, start, nullptr, &work);
        if (composed && state.failure) std::rethrow_exception(state.failure);
        if (!out.first_numerical_failure) out.first_numerical_failure = state.first_numerical_failure;
        return std::make_pair(std::move(compiled), std::move(state));
    };
    const auto usable = [](const NumericalEstimate& e) {
        return e.available && e.resolved && std::isfinite(e.absolute) && e.absolute >= 0;
    };
    const auto sign = [](double value) { return (value > 0) - (value < 0); };
    const auto evaluate_node = [&](double alpha, std::optional<EquilibriumResult> start, std::string prior_start, SaturationNode* retry = nullptr) -> std::optional<std::string> {
        const auto before = work.attempts;
        if (start) start->multipliers.clear();
        auto [c, state] = solve(alpha, affine.lower, affine.upper, affine.pressure,
            start, false);
        auto& n = retry ? *retry : out.nodes.emplace_back();
        n.driver_attempts = work.attempts - before;
        n.state = std::move(state); n.loading = alpha; n.prior_start = std::move(prior_start);
        if (n.state.failure) std::rethrow_exception(n.state.failure);
        if (!n.state.success || !n.state.has_state()) { return n.state.message; }
        n.mismatch = *n.state.pressure - affine.pressure;
        const auto estimate = [](const CompiledSystem& system, const EquilibriumResult& result, const PhasePoint* prepared = nullptr) {
            core::Scope scope("loading_node_uncertainty");
            return equilibrium_numerical_estimates(system, result, prepared, false, false, {}, true);
        };
        auto begun = Clock::now();
        n.estimates = estimate(c, n.state, n.state.prepared_endpoint.get());
        n.estimate_seconds = std::chrono::duration<double>(Clock::now() - begun).count();
        SolvedStateActionDirection direction; direction.feed_amounts_mol = affine.direction;
        direction.active_parameters.resize(model.active_parameters().size());
        begun = Clock::now();
        const auto action = solved_state_first_actions(c, n.state,
            {{SolvedStateObservableKind::Pressure}}, {direction});
        n.action_seconds = std::chrono::duration<double>(Clock::now() - begun).count();
        if (!action.empty()) n.action_diagnostics = std::make_shared<SolvedStateActionDiagnostics>(action[0].diagnostics);
        if (action.empty() || action[0].results.empty() || !action[0].results[0].action || action[0].results[0].status != SolvedStateActionStatus::Available) {
            n.action_cause = action.empty() || action[0].results.empty() ? "action_unavailable" : solved_state_status_name(action[0].results[0].status);
            if (!action.empty() && !action[0].results.empty())
                n.action_cause += ":" + action[0].results[0].diagnostic_message.value_or("");
            return n.action_cause;
        }
        n.slope = *action[0].results[0].action;
        auto local = c.problem;
        local.loading = affine; local.loading->seed = alpha; local.loading->pressure = *n.state.pressure;
        const auto augmented = compile_problem(model, local);
        auto mapped = n.state;
        mapped.coordinates.insert(mapped.coordinates.begin() + augmented.loading_column(), alpha);
        mapped.declaration = std::make_shared<EquilibriumProblem>(augmented.problem);
        mapped.prepared_endpoint.reset();
        mapped.numerical_estimates.reset();
        const auto local_errors = estimate(augmented, mapped);
        n.estimates.loading = local_errors.loading;
        n.local_estimates = local_errors;
        if (!usable(local_errors.loading)) { return local_errors.cause; }
        if (!usable(n.estimates.pressure)) { return n.estimates.cause; }
        return {};
    };
    const bool reuse = previous && previous->accepted && previous->grid == grid && previous->nodes.size() == grid.size()
        && previous->declaration.initial == affine.initial && previous->declaration.direction == affine.direction
        && previous->declaration.pressure == affine.pressure;
    bool bracket_found = false;
    try {
        for (double alpha : grid) {
            const auto i = out.nodes.size();
            std::string prior_start = reuse ? "bracket_continuation" : "not_supplied";
            SaturationNode* retry = nullptr;
            if (reuse && bracket_found) {
                const auto before = work.attempts;
                try {
                    const auto cause = evaluate_node(alpha, previous->nodes[i].state, "accepted");
                    if (!cause && !loading_slope_reversals(out.nodes.back(), &out.nodes[i - 1])) continue;
                    prior_start = "fallback:" + cause.value_or("unresolved_slope_reversal");
                } catch (const BudgetExhausted&) { throw; }
                catch (const std::bad_alloc&) { throw; }
                catch (const std::exception& error) {
                    core::retain_failure(out.first_numerical_failure, error, "prior loading node");
                    prior_start = std::string("fallback:") + error.what();
                }
                // Retain the failed start even if the shared budget stops its ordinary retry.
                if (out.nodes.size() == i) out.nodes.emplace_back();
                auto& rejected = out.nodes.back();
                rejected.loading = alpha; rejected.prior_start = prior_start; rejected.driver_attempts = work.attempts - before;
                out.rejected_starts.push_back(std::move(rejected));
                retry = &out.nodes.back();
                *retry = {}; retry->loading = alpha; retry->prior_start = prior_start;
            }
            auto cause = evaluate_node(alpha, i ? std::make_optional(out.nodes[i - 1].state) : std::nullopt, std::move(prior_start), retry);
            if (cause) { fail_at(alpha, *cause); fail("unexamined_gap"); break; }
            // Preserve the correction seed reached by the ordinary loading continuation.
            const auto& n = out.nodes.back();
            if (i && std::abs(n.mismatch) > n.estimates.pressure.absolute &&
                (sign(n.mismatch) != sign(out.nodes[i - 1].mismatch) ||
                 std::abs(out.nodes[i - 1].mismatch) <= out.nodes[i - 1].estimates.pressure.absolute))
                bracket_found = true;
        }
        if (out.nodes.size() == grid.size() && out.causes.empty()) {
            for (std::size_t i = 0; i < out.nodes.size(); ++i) {
                const auto& a = out.nodes[i];
                const bool event = std::abs(a.mismatch) <= a.estimates.pressure.absolute;
                if (event) out.events.push_back({i ? i - 1 : i, std::min(i + 1, grid.size() - 1)});
                if (!i || event) continue;
                const auto& b = out.nodes[i - 1];
                if (std::abs(b.mismatch) > b.estimates.pressure.absolute && sign(a.mismatch) != sign(b.mismatch)) {
                    ++out.brackets; out.events.push_back({i - 1, i});
                }
            }
            for (const auto& event : out.events) {
                const auto lo = grid[event[0]], hi = grid[event[1]];
                double seed = (lo + hi) / 2;
                auto start_node = event[0];
                for (std::size_t i = event[0]; i <= event[1]; ++i)
                    if (std::abs(out.nodes[i].mismatch) <= out.nodes[i].estimates.pressure.absolute) { seed = grid[i]; start_node = i; }
                auto [c, state] = solve(seed, lo, hi, affine.pressure, out.nodes[start_node].state, true);
                SaturationCandidate candidate; candidate.state = std::move(state);
                if (!candidate.state.success || !candidate.state.has_state()) {
                    fail(candidate.state.message); fail("unexamined_gap"); out.candidates.push_back(std::move(candidate)); continue;
                }
                const auto estimates = equilibrium_numerical_estimates(c, candidate.state);
                candidate.state.numerical_estimates = estimates;
                const auto alpha = candidate.state.coordinates[c.loading_column()];
                if (!usable(estimates.loading) || estimates.densities.size() != candidate.state.phase_count
                    || estimates.composition.size() != candidate.state.mole_fractions.size()) {
                    fail(estimates.cause); fail("unexamined_gap"); out.candidates.push_back(std::move(candidate)); continue;
                }
                candidate.admitted = true;
                const auto append = [&](double value, const NumericalEstimate& error, bool positive = false) {
                    candidate.operands.push_back({positive ? Quantity::Positive : Quantity::Linear, value, error});
                    if (!usable(error)) { candidate.admitted = false; fail("root_comparison_unresolved"); }
                };
                append(alpha, estimates.loading); append(*candidate.state.temperature, estimates.temperature);
                append(*candidate.state.pressure, estimates.pressure);
                for (std::size_t p = 0; p < candidate.state.phase_count; ++p) {
                    append(candidate.state.molar_densities[p], estimates.densities[p], true);
                    append(candidate.state.packing_fractions[p], estimates.packing[p], true);
                    append(candidate.state.amounts[p], estimates.amounts[p]);
                    for (std::size_t j = 0; j < c.n_feed.size(); ++j)
                        append(candidate.state.mole_fractions[p * c.n_feed.size() + j], estimates.composition[p * c.n_feed.size() + j]);
                }
                if (alpha < lo || alpha > hi || c.classified_residual_norm(candidate.state.coordinates) > policy::EQUILIBRIUM_ROW_RESIDUAL.value
                    || std::abs(*candidate.state.pressure / affine.pressure - 1) > policy::EQUILIBRIUM_ROW_RESIDUAL.value)
                    { candidate.admitted = false; fail("central_residual_refused"); fail("unexamined_gap"); }
                out.candidates.push_back(std::move(candidate));
            }
            out.complete = out.candidates.size() == out.events.size()
                && std::find(out.causes.begin(), out.causes.end(), "unexamined_gap") == out.causes.end();
        }
    } catch (const std::exception& error) {
        core::retain_failure(out.first_numerical_failure, error, "saturation search");
        fail(dynamic_cast<const BudgetExhausted*>(&error) ? "budget_exhausted" : error.what()); fail("unexamined_gap");
    }
    admit_loading_search(out);
    out.iterations = work.consumed; out.evaluations = work.evaluations;
    out.hessian_calls = work.hessian_calls; out.attempts = work.attempts;
    return out;
}

}  // namespace epcsaft::equilibrium
