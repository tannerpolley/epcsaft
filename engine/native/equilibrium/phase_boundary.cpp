#include "equilibrium/gibbs_minimizer.hpp"

#include "core/numerical_policy.hpp"
#include "core/trace.hpp"
#include "equilibrium/stoichiometry.hpp"

#include <Eigen/Dense>

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>

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

std::vector<double> seed_composition(const CompiledSystem& compiled, std::size_t phase,
                      std::vector<double> composition,
                      std::vector<double>& point, std::size_t offset, bool declared = false) {
    const auto& spec = compiled.problem.phases[phase];
    double total = 0.0;
    for (std::size_t i = 0; i < spec.support.size(); ++i) {
        const std::size_t species = spec.support[i];
        double value = composition[species];
        if (!(value > 0.0)) value = 1.0 / static_cast<double>(spec.support.size());
        point[offset + i] = value;
        total += value;
    }
    if (declared && spec.declared_composition_guess) {
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

void validate_density_seed(const DensityRoot& root) {
    if (!(std::isfinite(root.molar_density) && root.molar_density > 0.0)) {
        throw std::runtime_error("boundary seed has invalid density");
    }
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
                            const std::vector<std::size_t>& phase_offsets, const PhaseStart& start, core::Scope& scope, bool allow_loop = true) {
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
                states[side] = loop && pass == 0 && pressure == seed && composition == compiled.feed_composition && start.feed_states[side]
                    ? *start.feed_states[side] : compiled.model->state_at_pressure(temperature, pressure, composition,
                        side == 0 ? DensityBranch::Liquid : DensityBranch::Vapor);
            }
            if (!states[0].log_fugacity_coefficient || !states[1].log_fugacity_coefficient) break;
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
            const double stepped = loop ? pressure * std::exp(step) : seed * sum;
            if (!(std::isfinite(stepped) && stepped > 0.0)) break;
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
    } catch (const std::exception&) {}  // A failed start leaves the original path unchanged.
    return loop ? fugacity_pressure_seed(compiled, point, phase_offsets, start, scope, false) : fallback;
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
        if (observations) core::retain_rank(observations->rank_records, lu, "ideal_speciation");
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
                                                               std::size_t index = 0, StageWork* observations = nullptr) {
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
                phase.kind == PhaseKind::Vapor ? DensityBranch::Vapor : DensityBranch::Liquid, anchor);
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
            const auto rows = compiled.evaluate_at(compiled.prepare(trial, true, observations ? &observations->first_numerical_failure : nullptr)).values;
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
        const auto speciated = speciate_single_phase(compiled, point, compiled.feed_composition, scope, liquid, stages);
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


struct Attempt final {
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
    double t1 = 0.0, p1 = 0.0, h1 = 0.0, s1 = 0.0;  // the start values are base's

    EquilibriumProblem at(double q) const {
        EquilibriumProblem problem = base;
        problem.feed.kind = FeedKind::Amounts;
        problem.feed.values.resize(feed_start.size());
        for (std::size_t i = 0; i < feed_start.size(); ++i)
            problem.feed.values[i] = (1.0 - q) * feed_start[i] + q * feed_target[i];
        if (!problem.temperature.free) problem.temperature = FreeValue::fixed((1.0 - q) * base.temperature.value + q * t1);
        if (!problem.pressure.free) problem.pressure = FreeValue::fixed((1.0 - q) * base.pressure.value + q * p1);
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
        if (!base.enthalpy_j && !base.entropy_j_per_k) return result;
        result.resize(compiled.extended_control_width(), 0.0);
        result[result.size() - 2] = base.enthalpy_j ? h1 - *base.enthalpy_j : 0.0;
        result.back() = base.entropy_j_per_k ? s1 - *base.entropy_j_per_k : 0.0;
        return result;
    }
};

// A predicted free amount below zero is projected onto its bound rather than rejected: the
// amount bound is the disappearance boundary itself, and only the bound-constrained corrector
// can report the typed event there (counted in `projected`). A log mole fraction above zero is
// projected onto its corrector bound zero. Every other coordinate outside its domain is rejected.
bool inside_domain(const CompiledSystem& compiled, std::vector<double>& y, std::size_t& projected) {
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
    switch (result.convergence.stop) {
    case core::Stop::nonfinite: terminal = true; return ContinuationStop::nonfinite;
    case core::Stop::infeasible: terminal = true; return ContinuationStop::infeasible;
    case core::Stop::left_domain: return ContinuationStop::left_domain;
    case core::Stop::singular: return ContinuationStop::singular_jacobian;
    default: return ContinuationStop::stalled;
    }
}

// `supplied`: a caller's result; the driver's own results already carry this problem's classified rows.
Attempt admit_result(
    const HelmholtzModel& model, const CompiledSystem& compiled, EquilibriumResult result,
    const std::vector<double>& anchors, bool initial, core::Scope& scope, bool supplied = false
) {
    Attempt outcome;
    outcome.convergence = result.convergence;
    const auto refuse = [&](ContinuationStop stop, std::string message) {
        outcome.stop = stop;
        outcome.message = std::move(message);
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
    BranchCheck branches;
    try { branches = check_density_branches(model, compiled, result, anchors, initial, scope); }
    catch (const std::exception& error) { branches.message = error.what(); }
    if (!branches.accepted) branches.message += "; driver: " + result.message;
    outcome.accepted = branches.accepted;
    outcome.terminal = branches.terminal;
    outcome.message = branches.message;
    if (!branches.accepted) return refuse(ContinuationStop::branch_lost, branches.message);
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

void reject(ContinuationResult& result, const Attempt& attempt, double q) {
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
// domain and the same-driver corrector with carried anchors.
Attempt advance(const HelmholtzModel& trial, const CompiledSystem& start, const CompiledSystem& current,
                const std::vector<double>& y, const CompiledSystem& compiled, const EquilibriumProblem& problem,
                const std::function<std::vector<double>(std::optional<core::Convergence>*)>& tangent_at, double delta,
                const std::vector<double>& anchors, core::Scope& scope, std::shared_ptr<const OptimizerState> warm,
                const DriverLimits* limits = nullptr, StageWork* observations = nullptr) {
    Attempt attempt;
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
                + std::to_string(y.size()) + ", condition estimate " + std::to_string(info.condition_estimate) + ")"
            : "nonfinite operands in the implicit tangent";
    } else {
        std::vector<double> predicted = y;
        for (std::size_t i = 0; i < predicted.size(); ++i) predicted[i] += delta * tangent[i];
        std::size_t projected = 0;
        const bool inside = inside_domain(compiled, predicted, projected);
        scope.probe("predictor_projected", static_cast<double>(projected));
        if (inside) {
            EquilibriumResult audit;
            try { return admit_result(trial, compiled, solve_equilibrium_from_coordinates(trial, problem, predicted, nullptr, limits, std::move(warm), &audit), anchors, false, scope); }
            catch (...) { if (observations) observations->observe(audit); throw; }
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
void coalesced(Attempt& attempt, const CompiledSystem& compiled, std::pair<std::size_t, std::size_t> pair) {
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
    return solve_equilibrium(model, cell);
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
    EquilibriumResult audit;
    auto admitted = [&] {
        try { return admit_result(model, start, solve_equilibrium_from_coordinates(model, base,
            request.start_result ? request.start_result->coordinates : boundary_start(start, true, nullptr, &observations),
            nullptr, nullptr, {}, &audit), {}, true, scope); }
        catch (...) { observations.observe(audit); throw; }
    }();
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
        auto attempt = advance(model, start, current, y, compiled, problem, [&](auto*) { return tangent; }, s_try - s, expected, attempt_scope, transfer ? result.states.back().optimizer_state : nullptr, nullptr, &observations);
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

std::vector<double> boundary_start(const CompiledSystem& compiled, bool speciate, const DriverLimits* limits, StageWork* stages,
                                   const std::vector<bool>& packed) {
    core::Scope scope("boundary_start");
    if (limits) limits->check();
    std::vector<double> point = compiled.generic_point();
    const auto phase_offsets = compiled.phase_offsets();
    if (compiled.problem.critical) {
        // The declared T and density at the feed, P there, and u declared or the lowest mode of M (the
        // mode rows' u columns at this point).
        const auto& z = compiled.feed_composition;
        const double t = compiled.seed_temperature(), rho = compiled.problem.critical->density_guess_mol_per_m3;
        const std::size_t c = z.size(), mode = compiled.mode_offset();
        seed_composition(compiled, 0, z, point, phase_offsets[0]);
        point[phase_offsets[0] + c] = compiled.feed_total;
        point[phase_offsets[0] + c + 1] = std::log(compiled.model->packing_fraction_coefficient(t, z) * rho);
        point[compiled.primal_size - 2] = t;
        point[compiled.primal_size - 1] = compiled.model->state(t, rho, z).pressure;
        auto u = compiled.problem.critical->mode_guess;
        if (u.empty()) {
            const auto jacobian = compiled.jacobian(point);
            Eigen::MatrixXd M(c, c);
            for (std::size_t i = 0; i < c; ++i)
                for (std::size_t j = 0; j < c; ++j) M(i, j) = jacobian[(compiled.rows.size() - c - 2 + i) * point.size() + mode + j];
            const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigen(M);
            u.assign(eigen.eigenvectors().col(0).data(), eigen.eigenvectors().col(0).data() + c);
        }
        const double norm = std::sqrt(std::inner_product(u.begin(), u.end(), u.begin(), 0.0));
        for (std::size_t i = 0; i < c; ++i) point[mode + i] = u[i] / norm;
        return point;
    }
    if (compiled.has_caloric_condition()) {
        auto seed_problem = compiled.problem;
        seed_problem.enthalpy_j.reset();
        seed_problem.entropy_j_per_k.reset();
        const bool pure = compiled.model->parameters().component_count() == 1 && seed_problem.phases.size() == 2;
        if (pure) seed_problem.phases[1].amount = AmountSpec::pinned(0.0);
        else seed_problem.temperature = FreeValue::fixed(compiled.seed_temperature());
        const auto source = compile_problem(*compiled.model, seed_problem);
        // The seed alone stops at its own budget; a seed without a state is discarded either way (#579).
        auto seed_initial = boundary_start(source, true, limits, stages);
        DriverLimits budget = limits ? *limits : DriverLimits{std::numeric_limits<std::size_t>::max(), std::chrono::steady_clock::time_point::max()};
        budget.iterations = std::min(budget.iterations, budget.consumed + static_cast<std::size_t>(policy::EQUILIBRIUM_CALORIC_SEED_ITERATIONS.value));
        auto seed = structural_result(source, scope);
        std::exception_ptr failure;
        try {
            seed = solve_equilibrium_from_coordinates(*compiled.model, seed_problem, seed_initial, nullptr, &budget, {}, &seed);
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
        const auto inside = [&](double t) {  // every declared species' caloric record admits t
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
        } catch (const std::exception&) {}
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
    const double p = !start->densities.empty() ? start->pressure : fugacity_pressure_seed(compiled, point, phase_offsets, *start, scope);
    if (compiled.problem.pressure.free) point[compiled.primal_size - 1] = p;
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
        const DensityBranch branch = compiled.problem.phases[phase].kind == PhaseKind::Vapor
            ? DensityBranch::Vapor : DensityBranch::Liquid;
        DensityRoot root;
        try {
            const auto unanchored = [&] { return solve_density(*compiled.model, t, p, composition, branch, std::nullopt); };
            if (start->vapor_density > 0.0 && branch == DensityBranch::Vapor) try {  // the phi^V estimate's root continues
                root = solve_density(*compiled.model, t, p, composition, branch, start->vapor_density);
            } catch (const DensityStructureError&) { throw;
            } catch (const std::exception& error) { if (stages) retain_density_failure(stages->first_numerical_failure, error, "boundary anchor"); root = unanchored(); }
            else root = unanchored();
        } catch (const std::runtime_error& error) {
            if (stages) retain_density_failure(stages->first_numerical_failure, error, "boundary density");
            throw std::runtime_error(
                std::string("boundary seed has no density root: ") + error.what()
            );
        }
        if (stages) core::retain_first(stages->first_numerical_failure, root.first_numerical_failure);
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
        validate_density_seed(root);
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
    if (!(a.source_reference == b.source_reference)) return "source reference declaration differs";
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

void check_declared_side(const DensityRoot& root, DensityBranch branch) {
    const auto side = properties::density_side(root.branch, root.molar_density, root.negative_slope_density);
    const bool liquid = side == DensityBranch::Liquid;
    if (side != DensityBranch::Unspecified && side != branch)
        throw std::runtime_error(std::string(liquid ? "vapor" : "liquid") + " branch not admitted by the finite search: the sole stable root "
            + core::optional_number(root.molar_density) + " mol/m3 lies " + (liquid ? "above" : "below")
            + " the negative-slope sample " + core::optional_number(root.negative_slope_density) + " mol/m3");
}

}  // namespace

DensityRoot declared_root(const HelmholtzModel& model, double T, double P, const std::vector<double>& x, PhaseKind kind) {
    const auto branch = kind == PhaseKind::Vapor ? DensityBranch::Vapor : DensityBranch::Liquid;
    auto root = solve_density(model, T, P, x, branch, std::nullopt);
    check_declared_side(root, branch);
    return root;
}

namespace {

// Corrected-phase branch identity: the corrected density and the anchored locally unique
// root must share the anchor's window, with a positive corrected pressure slope. The physical
// state is recomputed from the coordinates, never taken from the supplied result's arrays.
BranchCheck density_branches_at(
    const HelmholtzModel& model, const CompiledSystem& compiled, const PhasePoint& prepared,
    const std::vector<double>& anchors, bool initial, core::Scope& scope, const DriverLimits* limits,
    const std::array<DensityRoot, 2>* roots
) {
    BranchCheck outcome;
    const bool pure = one_component_two_density_eligible(compiled);
    for (std::size_t phase = 0; phase < compiled.problem.phases.size(); ++phase) {
        if (limits) limits->check();
        const auto& spec = compiled.problem.phases[phase];
        std::vector<double> x = prepared.composition[phase];
        const double total = std::accumulate(x.begin(), x.end(), 0.0);
        for (double& value : x) value /= total;  // the driver's own normalisation of the log chart
        const double corrected = (model.parameters().ideal_gas ? total : 1.0) / prepared.volume[phase];
        const auto branch = spec.kind == PhaseKind::Vapor ? DensityBranch::Vapor : DensityBranch::Liquid;
        const double slope = pure ? prepared.pure[phase].pressure_density_derivative
            : model.state(prepared.temperature, corrected, x).pressure_density_derivative;
        if (limits) limits->check();
        double anchor = anchors.empty() ? 0.0 : anchors[phase];
        DensityRoot root;
        try {
            if (initial) {
                if (roots) {
                    root = (*roots)[spec.kind == PhaseKind::Vapor ? 1 : 0];
                    check_declared_side(root, branch);
                } else {
                    root = declared_root(model, prepared.temperature, prepared.pressure, x, spec.kind);
                }
                if (pure && root.stable_root_count.value_or(0) < 2) {
                    outcome.message = "one-component two-density requires two distinct stable density branches"; return outcome;
                }
                anchor = root.molar_density;
                if (limits) limits->check();
            }
            root = solve_density(model, prepared.temperature, prepared.pressure, x, branch, anchor);
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
    }
    // Declared kinds: a liquid less densely packed than a vapor is the other boundary with its labels swapped (#599).
    for (std::size_t l = 0; initial && l < prepared.eta.size(); ++l)
        for (std::size_t v = 0; v < prepared.eta.size(); ++v)
            if (compiled.problem.phases[l].kind == PhaseKind::Liquid && compiled.problem.phases[v].kind == PhaseKind::Vapor
                && !(prepared.eta[l] > prepared.eta[v])) {
                outcome.message = "declared liquid " + compiled.problem.phases[l].name + " (packing fraction " + core::optional_number(prepared.eta[l])
                    + ") is not more densely packed than declared vapor " + compiled.problem.phases[v].name + " ("
                    + core::optional_number(prepared.eta[v]) + "): swapped phase labels";
                return outcome;
            }
    outcome.accepted = true;
    return outcome;
}

}  // namespace

BranchCheck check_density_branches(
    const HelmholtzModel& model, const CompiledSystem& compiled, const EquilibriumResult& result,
    const std::vector<double>& anchors, bool initial, core::Scope& scope, const DriverLimits* limits,
    const PhasePoint* supplied
) {
    if (limits) limits->check();
    const auto owned = supplied ? std::optional<PhasePoint>{} : std::make_optional(compiled.prepare(result.coordinates, true));
    const auto& prepared = supplied ? *supplied : *owned;
    if (initial && !limits && one_component_two_density_eligible(compiled)) {
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
            const auto shared = density_branches_at(model, compiled, prepared, anchors, initial, scope, limits, &roots);
            if (shared.accepted) return shared;
        }
    }
    return density_branches_at(model, compiled, prepared, anchors, initial, scope, limits, nullptr);
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
    path.h1 = request.target_problem.enthalpy_j.value_or(0.0);
    path.s1 = request.target_problem.entropy_j_per_k.value_or(0.0);
    const bool moves = request.start_problem.enthalpy_j != request.target_problem.enthalpy_j
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
    CompiledSystem current = start;
    std::unique_ptr<HelmholtzModel> current_model;
    {
        core::Scope initial("equilibrium_continuation_initial");
        initial.probe("path_value", 0.0);
        Attempt admitted;
        const auto solve_initial = [&](const std::vector<double>& coordinates) {
            EquilibriumResult audit;
            try { return solve_equilibrium_from_coordinates(model, request.start_problem, coordinates, nullptr, nullptr, {}, &audit); }
            catch (...) { observations.observe(audit); throw; }
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
            admitted = admit_result(model, start, parameter_path
                ? solve_initial(supplied.coordinates) : supplied, anchors, anchors.empty(), initial,
                !parameter_path);
        } else {
            admitted = admit_result(model, start,
                solve_initial(boundary_start(start, true, nullptr, &observations)),
                anchors, true, initial);
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
    }

    double q = 0.0;
    double dq = options.initial_step;
    std::size_t local_retries = 0;
    bool transfer = true;
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
        Attempt attempt = advance(trial, start, current, y_current, compiled, problem, [&](auto* failure) {
            return sensitivity(current, y_current, path.lever(current), failure);
        }, q_try - q, anchors, attempt_scope, transfer ? result.states.back().optimizer_state : nullptr, limits, &observations);
        observations.observe(attempt.observations);
        if (attempt.accepted)
            if (const auto pair = swapped_pair(compiled, y_current, attempt.result->coordinates)) coalesced(attempt, compiled, *pair);
        if (!attempt.result.has_value()) attempt_scope.probe("driver_result_present", 0.0);
        attempt_scope.probe("accepted", attempt.accepted ? 1.0 : 0.0);
        attempt_scope.probe(continuation_stop_name(attempt.accepted ? ContinuationStop::reached_target : attempt.stop), 1.0);
        if (attempt.convergence.has_value()) attempt_scope.finish(*attempt.convergence);
        if (attempt.accepted) {
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
        const double halved = dq / 2.0;
        if (halved < options.minimum_step)
            return finish(std::move(result), ContinuationStop::step_too_small, attempt.message);
        if (local_retries >= options.max_retries_per_step)
            return finish(std::move(result), ContinuationStop::budget_exhausted, attempt.message);
        dq = halved;
        ++local_retries;
        ++result.retries;
    }
}

ContinuationResult continue_equilibrium(const HelmholtzModel& model, const ContinuationRequest& request) {
    core::Scope scope("continue_equilibrium");
    return continue_path(model, request, nullptr, {});
}

// ALG-EQUILIBRIUM neutral recovery: the first admitted source at T_seed r^j, then fixed-T continuation, or for a
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
    for (auto& phase : declared.phases) { phase.composition_guess.clear(); phase.packing_guess.reset(); }
    if (declared.pressure.free) declared.pressure = FreeValue::variable();
    auto source = declared;
    if (free_temperature) source.pressure = FreeValue::variable();
    try {
        for (int j = 1; j <= static_cast<int>(policy::EQUILIBRIUM_RECOVERY_SOURCES.value) && !found.success; ++j) {
            limits.check();  // the slot is spent before compilation
            source.temperature = FreeValue::fixed(target.seed_temperature() * std::pow(policy::COOLING_RATIO.value, j));
            EquilibriumResult audit;
            try { found = solve_compiled(compile_problem(model, source), std::nullopt, nullptr, &limits, {}, false, &audit); }
            catch (...) { observations.observe(audit); throw; }
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
            EquilibriumResult corrected;
            try { corrected = solve_equilibrium_from_coordinates(model, request.start_problem, chart, nullptr, &limits, {}, &corrected); }
            catch (...) { observations.observe(corrected); throw; }
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
        const auto admitted = check_density_branches(model, target, reached, {}, true, scope, &limits);
        if (!admitted.accepted || !(target.classified_residual_norm(reached.coordinates) <= tolerance)) {
            refused.message = "reached target refused: " + admitted.message;
                core::retain_failure(observations.first_numerical_failure, refused.message, "source recovery admission");
            return publish(std::move(refused));
        }
        reached.declaration = std::make_shared<const EquilibriumProblem>(target.problem);
        return publish(std::move(reached));
    } catch (const std::exception& error) {
        retain_density_failure(observations.first_numerical_failure, error, "source recovery");
        refused.message = error.what();
        return publish(std::move(refused));
    }
}


void compare_loading_candidates(SaturationSearch& out, core::ComparisonEligibility eligibility) {
    const auto usable = [](const core::NumericalEstimate& e) {
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
                const auto c = core::compare_equilibrium_quantity(a[k].quantity, a[k].value, b[k].value, 0,
                    a[k].error, b[k].error, e, core::ComparisonMode::EstimatedError);
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

void admit_loading_search(SaturationSearch& out) {
    const auto sign = [](double value) { return (value > 0) - (value < 0); };
    compare_loading_candidates(out, {});
    for (std::size_t i = 0; i < out.nodes.size(); ++i) {
        const auto& a = out.nodes[i];
        if (!std::isfinite(a.slope) || a.slope == 0) out.causes.push_back("unresolved_slope_reversal");
        if (!i) continue;
        const auto& b = out.nodes[i - 1]; const auto delta = a.mismatch - b.mismatch;
        if (sign(a.slope) != sign(b.slope) || (std::abs(delta) > a.estimates.pressure.absolute + b.estimates.pressure.absolute
            && sign(delta) != sign(a.slope))) out.causes.push_back("unresolved_slope_reversal");
    }

    if (out.groups.empty() && out.causes.empty()) out.causes.push_back("no_resolved_root");
    if (out.complete && out.causes.empty() && out.groups.size() == 1) out.accepted = out.groups[0][0];
}

SaturationSearch search_loading(const HelmholtzModel& model, EquilibriumProblem declaration, const AffineLoading& affine,
    const std::vector<double>& grid, const DriverLimits* supplied, std::chrono::steady_clock::time_point deadline) {
    using Clock = std::chrono::steady_clock;
    if (grid.size() < 2 || grid.front() != affine.lower || grid.back() != affine.upper)
        throw std::invalid_argument("saturation grid must include both declared bounds");
    for (std::size_t i = 0; i < grid.size(); ++i)
        if (!std::isfinite(grid[i]) || grid[i] <= 0 || (i && grid[i] <= grid[i - 1]))
            throw std::invalid_argument("saturation grid must be positive and strictly increasing");
    const auto cap = static_cast<std::size_t>(policy::EQUILIBRIUM_RECOVERY_WORK.value);
    DriverLimits own{cap, deadline, 0, true, 4 * cap, 2 * cap};
    own.attempt_cap = cap;
    const auto& work = supplied ? *supplied : own;
    SaturationSearch out; out.grid = grid; out.declaration = affine;
    out.work_at_entry = {work.consumed, work.evaluations, work.hessian_calls, work.attempts};
    const auto fail = [&](const std::string& cause) { out.causes.push_back(cause); };
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
        auto state = solve_compiled(compiled, start, nullptr, &work);
        if (!out.first_numerical_failure) out.first_numerical_failure = state.first_numerical_failure;
        return std::make_pair(std::move(compiled), std::move(state));
    };
    const auto usable = [](const core::NumericalEstimate& e) {
        return e.available && e.resolved && std::isfinite(e.absolute) && e.absolute >= 0;
    };
    const auto sign = [](double value) { return (value > 0) - (value < 0); };
    try {
        for (double alpha : grid) {
            auto [c, state] = solve(alpha, affine.lower, affine.upper, affine.pressure,
                out.nodes.empty() ? std::nullopt : std::make_optional(out.nodes.back().state), false);
            SaturationNode node; node.loading = alpha; node.state = std::move(state);
            out.nodes.push_back(std::move(node)); auto& n = out.nodes.back();
            if (!n.state.success || !n.state.has_state()) { fail(n.state.message); fail("unexamined_gap"); break; }
            n.mismatch = *n.state.pressure - affine.pressure;
            auto begun = Clock::now();
            n.estimates = equilibrium_numerical_estimates(c, n.state);
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
                fail(n.action_cause); fail("unexamined_gap"); break;
            }
            n.slope = *action[0].results[0].action;
            auto local = c.problem;
            local.loading = affine; local.loading->seed = alpha; local.loading->pressure = *n.state.pressure;
            const auto augmented = compile_problem(model, local);
            auto mapped = n.state;
            mapped.coordinates.insert(mapped.coordinates.begin() + augmented.loading_column(), alpha);
            mapped.declaration = std::make_shared<EquilibriumProblem>(augmented.problem);
            const auto local_errors = equilibrium_numerical_estimates(augmented, mapped);
            n.estimates.loading = local_errors.loading;
            n.local_estimates = local_errors;
            if (!usable(local_errors.loading)) { fail(local_errors.cause); fail("unexamined_gap"); break; }
            if (!usable(n.estimates.pressure)) { fail(n.estimates.cause); fail("unexamined_gap"); break; }
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
                const auto append = [&](double value, const core::NumericalEstimate& error, bool positive = false) {
                    candidate.operands.push_back({positive ? core::Quantity::Positive : core::Quantity::Linear, value, error});
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
        fail(work.expired() ? "budget_exhausted" : error.what()); fail("unexamined_gap");
    }
    admit_loading_search(out);
    out.iterations = work.consumed; out.evaluations = work.evaluations;
    out.hessian_calls = work.hessian_calls; out.attempts = work.attempts;
    return out;
}

}  // namespace epcsaft::equilibrium
