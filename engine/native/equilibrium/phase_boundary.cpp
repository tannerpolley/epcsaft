#include "equilibrium/gibbs_minimizer.hpp"

#include "core/numerical_policy.hpp"
#include "core/trace.hpp"
#include "equilibrium/stoichiometry.hpp"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
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
        && !compiled.has_ions && !compiled.has_reaction && compiled.pinned_count == 0
        && std::all_of(compiled.feed_composition.begin(), compiled.feed_composition.end(),
            [](double z) { return std::isfinite(z) && z > 0.0; })
        && p.phases[0].support.size() == components && p.phases[1].support.size() == components;
}

namespace {

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
    const CompiledSystem& compiled, core::Scope& scope, std::optional<double>& vapor_fraction
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
    DensityRoot liquid_root;
    const auto vapor_root = solve_density(
        *compiled.model, temperature, pressure, compiled.feed_composition,
        DensityBranch::Vapor, std::nullopt, &liquid_root);
    const auto vapor = compiled.model->state(
        temperature, vapor_root.molar_density, compiled.feed_composition);
    const auto liquid = compiled.model->state(
        temperature, liquid_root.molar_density, compiled.feed_composition);
    if (!vapor.log_fugacity_coefficient || !liquid.log_fugacity_coefficient)
        throw std::runtime_error("boundary start has no fugacity coefficients");
    // A sole unlabelled feed root cannot supply both phases' fugacity estimates:
    // their self-ratio is one. Use the ideal-vapor limit for this start only.
    const auto ln_phi_vapor = vapor_root.stable_root_count == 1
        ? std::vector<std::optional<double>>(components, 0.0) : *vapor.log_fugacity_coefficient;
    const auto& ln_phi_liquid = *liquid.log_fugacity_coefficient;
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
    if (supports_split_seed(compiled)
        && std::all_of(k_values.begin(), k_values.end(), [](double k) { return std::isfinite(k) && k > 0.0; })) {
        const auto balance = [&](double beta) {
            double value = 0.0;
            for (std::size_t i = 0; i < components; ++i)
                value += compiled.feed_composition[i] * (k_values[i] - 1.0)
                    / ((1.0 - beta) + beta * k_values[i]);
            return value;
        };
        const double left = balance(0.0), right = balance(1.0);
        if (std::isfinite(left) && std::isfinite(right) && left > 0.0 && right < 0.0) {
            double lower = 0.0, upper = 1.0;
            for (int step = 0; step < std::numeric_limits<double>::digits; ++step) {
                const double beta = 0.5 * (lower + upper);
                if (balance(beta) > 0.0) lower = beta;
                else upper = beta;
            }
            const double beta = 0.5 * (lower + upper);
            std::vector<double> trial_x(components), trial_y(components);
            bool usable = beta > policy::STABILITY_PHASE_AMOUNT_FLOOR.value
                && 1.0 - beta > policy::STABILITY_PHASE_AMOUNT_FLOOR.value;
            for (std::size_t i = 0; i < components; ++i) {
                trial_x[i] = compiled.feed_composition[i] / ((1.0 - beta) + beta * k_values[i]);
                trial_y[i] = k_values[i] * trial_x[i];
                usable = usable && std::isfinite(trial_x[i]) && trial_x[i] > 0.0
                    && std::isfinite(trial_y[i]) && trial_y[i] > 0.0;
            }
            if (usable) {
                vapor_fraction = beta;
                x = std::move(trial_x); y = std::move(trial_y);
                scope.probe("vapor_fraction_seed", beta);
            }
        }
    }
    for (std::size_t phase = 0; phase < phases; ++phase)
        result[phase] = compiled.problem.phases[phase].kind == PhaseKind::Vapor ? y : x;
    scope.probe("start_rule", 2.0);
    return result;
}

// A mixture liquid with a vapor pinned at zero amount, at free pressure and fixed temperature:
// move the caller's pressure seed by one fugacity step, P1 = P0 sum_{i in V} x_i phi_i^L/phi_i^V,
// from the seeded compositions (phi^L ~ 1/P: the seeded liquid's saturation pressure) and the same step's
// y_i = x_i (phi_i^L/phi_i^V) P0/P1. A failed state keeps the caller's seeds; equations and tolerances do not move.
double fugacity_pressure_seed(const CompiledSystem& compiled, std::vector<double>& point,
                            const std::vector<std::size_t>& phase_offsets, core::Scope& scope) {
    const auto& problem = compiled.problem;
    const double seed = seed_pressure(compiled);
    // Pure saturation keeps its coexistence-root start.
    if (!problem.pressure.free || problem.temperature.free || problem.phases.size() != 2
        || problem.phases[0].kind == problem.phases[1].kind
        || compiled.model->parameters().component_count() < 2) return seed;
    const std::size_t vapor = problem.phases[0].kind == PhaseKind::Vapor ? 0 : 1;
    const auto& pinned = problem.phases[vapor].amount;
    if (pinned.kind != AmountKind::Pinned || pinned.fraction != 0.0) return seed;
    const std::size_t components = compiled.model->parameters().component_count();
    std::vector<std::vector<double>> x(2, std::vector<double>(components, 0.0));
    for (std::size_t phase = 0; phase < 2; ++phase) {
        const auto& support = problem.phases[phase].support;
        for (std::size_t i = 0; i < support.size(); ++i)
            x[phase][support[i]] = std::exp(point[phase_offsets[phase] + i]);
    }
    const double temperature = seed_temperature(compiled);
    try {
        const auto liquid = compiled.model->state_at_pressure(
            temperature, seed, x[1 - vapor], DensityBranch::Liquid);
        const auto gas = compiled.model->state_at_pressure(
            temperature, seed, x[vapor], DensityBranch::Vapor);
        if (!liquid.log_fugacity_coefficient || !gas.log_fugacity_coefficient) return seed;
        const auto& support = problem.phases[vapor].support;
        std::vector<double> y(support.size());
        double sum = 0.0;
        for (std::size_t k = 0; k < support.size(); ++k) {
            const std::size_t i = support[k];
            y[k] = x[1 - vapor][i] * std::exp(require_component_property((*liquid.log_fugacity_coefficient)[i])
                                              - require_component_property((*gas.log_fugacity_coefficient)[i]));
            sum += y[k];
        }
        const double stepped = seed * sum;
        if (!(std::isfinite(stepped) && stepped > 0.0)) return seed;
        // y is as good as the seeded liquid: only the all-guessed caller start (start_compositions) carries one;
        // the default z/K liquid would return y = z, the trivial solution.
        if (std::all_of(problem.phases.begin(), problem.phases.end(),
                [](const PhaseSpec& phase) { return phase.composition_guess.size() == phase.support.size(); }))
            for (std::size_t k = 0; k < support.size(); ++k)
                if (y[k] > 0.0) point[phase_offsets[vapor] + k] = std::log(y[k] / sum);
        scope.probe("fugacity_pressure_seed", stepped);
        return stepped;
    } catch (const std::exception&) {
        return seed;
    }
}

// Ideal-solution chemical equilibrium with frozen c_j: minimize sum_j n_j (c_j + ln n_j/n) subject to B n = b
// by the Gordon-McBride element-potential Newton iteration and step control (NASA RP-1311, 1994, eqs
// 2.24-2.26, 3.1-3.3), the box separating major from trace species; the trace limit applies only to a
// species growing faster than the total.
bool ideal_speciation(const std::vector<std::vector<double>>& B, const std::vector<double>& b,
                      const std::vector<double>& c, std::vector<double>& n) {
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
        const Eigen::VectorXd x = A.fullPivLu().solve(rhs);
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
        if (!std::isfinite(change)) return false;
        if (change <= policy::EQUILIBRIUM_START_SPECIATION_STEP.value) return true;
        control = std::max(0.0, std::min(control, 2.0 / std::max(major, 2.0)));
        for (std::size_t j = 0; j < m; ++j) ln_n[j] += control * g[j];
        ln_total += control * x(r);
    }
    return false;
}

// A reacting declaration of one free phase starts from the ideal-solution speciation of the feed balances,
// with the residual, volume and standard terms of its chemical rows frozen at the seeded composition and
// refreshed over sweeps (#154, #157): an unspeciated start leaves the log chart's trace
// species without restoring force and the optimizer cycles. Zero guesses enter at the optimizer box; a
// failed sweep keeps the last speciation.
// Returns the speciated amount and the last sweep's density root, which seeds the start's packing.
std::optional<std::pair<double, double>> speciate_single_phase(const CompiledSystem& compiled, std::vector<double>& point,
                                                               std::vector<double> composition, core::Scope& scope) {
    const auto& phase = compiled.problem.phases.front();
    const std::size_t m = phase.support.size();
    const double t = seed_temperature(compiled), p = seed_pressure(compiled);
    std::vector<std::vector<double>> B(compiled.rank, std::vector<double>(m));
    for (std::size_t k = 0; k < compiled.rank; ++k)
        for (std::size_t j = 0; j < m; ++j) B[k][j] = compiled.B[k][phase.support[j]];
    for (const std::size_t i : phase.support) composition[i] = std::max(composition[i], policy::EQUILIBRIUM_COMPOSITION_BOX.value);
    std::vector<double> n(m), speciated, trial = point;
    for (std::size_t j = 0; j < m; ++j) n[j] = compiled.feed_total * composition[phase.support[j]];
    std::optional<double> density, speciated_density;
    try {
        const auto root = [&](std::optional<double> anchor) {
            return solve_density(*compiled.model, t, p, composition,
                phase.kind == PhaseKind::Vapor ? DensityBranch::Vapor : DensityBranch::Liquid, anchor).molar_density;
        };
        for (int sweep = 0; sweep < 2; ++sweep) {
            const double total = std::accumulate(composition.begin(), composition.end(), 0.0);
            for (double& value : composition) value /= total;
            try {  // continue the previous sweep's root; a speciation jump beyond its window solves it anew
                density = root(density);
            } catch (const std::exception&) {
                if (!density) throw;
                density = root(std::nullopt);
            }
            for (std::size_t j = 0; j < m; ++j) trial[j] = std::log(composition[phase.support[j]]);
            trial[m] = compiled.feed_total;
            trial[m + 1] = std::log(compiled.model->packing_fraction_coefficient(t, composition) * *density);
            std::fill(trial.begin() + static_cast<std::ptrdiff_t>(compiled.primal_size), trial.end(), 0.0);
            const auto rows = compiled.evaluate_at(compiled.prepare(trial, true)).values;
            std::vector<double> c(m);
            for (std::size_t j = 0; j < m; ++j) c[j] = rows[j] - trial[j];  // chemical rows lead, in support order
            if (!ideal_speciation(B, compiled.balance, c, n)) break;
            speciated = n;
            speciated_density = density;
            for (std::size_t j = 0; j < m; ++j) composition[phase.support[j]] = n[j];
        }
    } catch (const std::exception&) {
    }
    if (speciated.empty()) return std::nullopt;
    n = speciated;
    const double total = std::accumulate(n.begin(), n.end(), 0.0);
    for (std::size_t j = 0; j < m; ++j) point[j] = std::log(n[j] / total);
    scope.probe("speciated_amount", total);
    return std::make_pair(total, *speciated_density);
}

}  // namespace

std::vector<double> boundary_start(const CompiledSystem& compiled, bool speciate) {
    core::Scope scope("boundary_start");
    std::vector<double> point = compiled.generic_point();
    const auto phase_offsets = compiled.phase_offsets();
    std::optional<double> vapor_fraction;
    const auto compositions = start_compositions(compiled, scope, vapor_fraction);
    std::size_t free_phases = 0;
    for (const auto& phase : compiled.problem.phases) {
        if (phase.amount.kind == AmountKind::Free) ++free_phases;
    }
    double feed_total = 0.0;
    for (const double value : compiled.n_feed) feed_total += value;
    std::optional<double> speciated_density;
    double free_amount = compiled.n_feed.empty() ? 1.0
        : feed_total / static_cast<double>(std::max<std::size_t>(1, free_phases));
    scope.probe("free_phases", static_cast<double>(free_phases));
    scope.probe("feed_total", feed_total);
    for (std::size_t phase = 0; phase < compiled.problem.phases.size(); ++phase)
        seed_composition(compiled, phase, compositions[phase], point, phase_offsets[phase]);
    if (speciate && compiled.has_reaction && compiled.problem.phases.size() == 1 && free_phases == 1)
        if (const auto speciated = speciate_single_phase(compiled, point, compositions.front(), scope))
            std::tie(free_amount, speciated_density) = *speciated;
    const double p = fugacity_pressure_seed(compiled, point, phase_offsets, scope);
    if (compiled.problem.pressure.free) point[compiled.primal_size - 1] = p;
    for (std::size_t phase = 0; phase < compiled.problem.phases.size(); ++phase) {
        std::size_t cursor = phase_offsets[phase] + compiled.problem.phases[phase].support.size();
        if (compiled.problem.phases[phase].amount.kind == AmountKind::Free) {
            point[cursor++] = vapor_fraction.has_value()
                ? feed_total * (compiled.problem.phases[phase].kind == PhaseKind::Vapor
                    ? *vapor_fraction : 1.0 - *vapor_fraction) : free_amount;
        }
        const double t = seed_temperature(compiled);
        std::vector<double> composition(compiled.model->parameters().component_count(), 0.0);
        for (std::size_t i = 0; i < compiled.problem.phases[phase].support.size(); ++i) {
            composition[compiled.problem.phases[phase].support[i]] =
                std::exp(point[phase_offsets[phase] + i]);
        }
        if (supports_split_seed(compiled) && compiled.problem.phases[phase].packing_guess.has_value()) {
            const double packing = *compiled.problem.phases[phase].packing_guess;
            if (!(std::isfinite(packing)
                  && packing >= policy::PACKING_FRACTION_FLOOR.value
                  && packing <= policy::PACKING_FRACTION_CEILING.value)) {
                throw std::invalid_argument("packing guess is outside the admitted range");
            }
            point[cursor] = std::log(packing);
            scope.probe("packing_guess", packing, static_cast<std::uint32_t>(phase));
            continue;
        }
        if (speciated_density.has_value()) {  // the last sweep's root, at nearly the speciated composition
            point[cursor] = std::log(compiled.model->packing_fraction_coefficient(t, composition) * *speciated_density);
            continue;
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
        // Coexistence needs two stable pure roots only when more than one phase is declared;
        // a single declared pure phase uses the existing unique-root branch rules.
        if (compiled.model->parameters().component_count() == 1
            && compiled.problem.phases.size() > 1
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

// ---------------------------------------------------------------- continuation

namespace epcsaft::equilibrium {
namespace {

struct Attempt final {
    bool accepted = false;
    bool terminal = false;
    ContinuationStop stop = ContinuationStop::stalled;
    std::string message;
    std::optional<EquilibriumResult> result;
    std::optional<core::Convergence> convergence;
};

bool same_intensive(const FreeValue& a, const FreeValue& b) {
    if (a.free != b.free) return false;
    return !a.free || (a.absent == b.absent && a.value == b.value);
}

// Every field that is not a path lever must be identical at both endpoints.
std::string declaration_difference(const EquilibriumProblem& a, const EquilibriumProblem& b) {
    const auto &u = a.source_reference, &v = b.source_reference;
    if (u.has_value() != v.has_value() || (u && (u->solvent_index != v->solvent_index
        || u->solute_molality_ratios != v->solute_molality_ratios || u->standard_molality_mol_per_kg != v->standard_molality_mol_per_kg)))
        return "source reference declaration differs";
    if (a.phases.empty() || b.phases.empty()) return "endpoint declares no phases";
    if (a.phases.size() != b.phases.size()) return "phase count differs";
    for (std::size_t p = 0; p < a.phases.size(); ++p) {
        const auto& x = a.phases[p];
        const auto& y = b.phases[p];
        if (x.name != y.name || x.kind != y.kind) return "phase " + x.name + " identity differs";
        if (x.support != y.support || x.support_declared != y.support_declared)
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
    if (a.topology_declared != b.topology_declared) return "topology declaration differs";
    if (a.feed.kind != b.feed.kind || a.feed.solvent != b.feed.solvent
        || a.feed.solvent_mass_kg != b.feed.solvent_mass_kg || a.feed.values.size() != b.feed.values.size())
        return "feed basis, solvent set or component count differs";
    if (a.row_tolerance != b.row_tolerance) return "row_tolerance differs";
    if (!same_intensive(a.temperature, b.temperature)) return "temperature free/fixed declaration differs";
    if (!same_intensive(a.pressure, b.pressure)) return "pressure free/fixed declaration differs";
    return "";
}

std::string layout_difference(const CompiledSystem& a, const CompiledSystem& b) {
    for (std::size_t p = 0; p < a.problem.phases.size(); ++p)
        if (a.problem.phases[p].support != b.problem.phases[p].support)
            return "admitted support of phase " + a.problem.phases[p].name + " differs";
    if (a.rows.size() != b.rows.size()) return "row count differs";
    for (std::size_t r = 0; r < a.rows.size(); ++r)
        if (a.rows[r].kind != b.rows[r].kind || a.rows[r].phase != b.rows[r].phase
            || a.rows[r].index != b.rows[r].index || std::string(a.rows[r].eqid) != b.rows[r].eqid)
            return "row " + std::to_string(r) + " identity differs";
    if (a.primal_size != b.primal_size || a.unknown_count() != b.unknown_count() || a.rank != b.rank
        || a.charge_row_count != b.charge_row_count || a.pinned_count != b.pinned_count
        || a.free_intensive_count != b.free_intensive_count)
        return "coordinate layout differs";
    return "";
}

struct Path final {
    EquilibriumProblem base;
    std::vector<double> feed_start, feed_target;
    std::vector<double> parameter_start, parameter_target;
    double t0 = 0.0, t1 = 0.0, p0 = 0.0, p1 = 0.0;

    EquilibriumProblem at(double q) const {
        EquilibriumProblem problem = base;
        problem.feed.kind = FeedKind::Amounts;
        problem.feed.values.resize(feed_start.size());
        for (std::size_t i = 0; i < feed_start.size(); ++i)
            problem.feed.values[i] = (1.0 - q) * feed_start[i] + q * feed_target[i];
        if (!problem.temperature.free) problem.temperature = FreeValue::fixed((1.0 - q) * t0 + q * t1);
        if (!problem.pressure.free) problem.pressure = FreeValue::fixed((1.0 - q) * p0 + q * p1);
        return problem;
    }

    std::vector<double> lever(std::size_t width) const {
        std::vector<double> result(width, 0.0);
        result[0] = base.temperature.free ? 0.0 : t1 - t0;
        result[1] = base.pressure.free ? 0.0 : p1 - p0;
        for (std::size_t i = 0; i < feed_start.size(); ++i) result[2 + i] = feed_target[i] - feed_start[i];
        for (std::size_t i = 0; i < parameter_target.size(); ++i)
            result[2 + feed_start.size() + i] = parameter_target[i] - parameter_start[i];
        return result;
    }
};

// A predicted free amount below zero is projected onto its bound rather than rejected: the
// amount bound is the disappearance boundary itself, and only the bound-constrained corrector
// can report the typed event there (counted in `projected`). Every other coordinate outside
// its domain is rejected.
bool inside_domain(const CompiledSystem& compiled, std::vector<double>& y, std::size_t& projected) {
    if (!std::all_of(y.begin(), y.end(), [](double v) { return std::isfinite(v); })) return false;
    std::size_t offset = 0;
    for (const auto& phase : compiled.problem.phases) {
        for (std::size_t i = 0; i < phase.support.size(); ++i, ++offset)
            if (y[offset] < std::log(policy::STABILITY_COMPOSITION_FLOOR.value) || y[offset] > 0.0) return false;
        if (phase.amount.kind == AmountKind::Free) {
            if (y[offset] < 0.0) { y[offset] = 0.0; ++projected; }
            ++offset;
        }
        if (y[offset] < std::log(policy::PACKING_FRACTION_FLOOR.value)
            || y[offset] > std::log(policy::PACKING_FRACTION_CEILING.value)) return false;
        ++offset;
    }
    if (compiled.problem.temperature.free && !(y[offset++] > 0.0)) return false;
    if (compiled.problem.pressure.free && !(y[offset++] > 0.0)) return false;
    return true;
}

// Per-phase normalised composition and log packing read from the coordinate chart.
struct PhaseSignature final {
    std::vector<double> composition;
    double log_eta = 0.0;
};

std::vector<PhaseSignature> signatures(const CompiledSystem& compiled, const std::vector<double>& y) {
    const auto phase_offsets = compiled.phase_offsets();
    const std::size_t components = compiled.model->parameters().component_count();
    std::vector<PhaseSignature> result;
    for (std::size_t phase = 0; phase < compiled.problem.phases.size(); ++phase) {
        const auto& spec = compiled.problem.phases[phase];
        PhaseSignature item;
        item.composition.assign(components, 0.0);
        double total = 0.0;
        for (std::size_t i = 0; i < spec.support.size(); ++i) {
            item.composition[spec.support[i]] = std::exp(y[phase_offsets[phase] + i]);
            total += item.composition[spec.support[i]];
        }
        for (double& value : item.composition) value /= total;
        std::size_t cursor = phase_offsets[phase] + spec.support.size();
        if (spec.amount.kind == AmountKind::Free) ++cursor;
        item.log_eta = y[cursor];
        result.push_back(std::move(item));
    }
    return result;
}

// Passage through a critical point between consecutive accepted states swaps a pair's
// identities: every nonzero component of the inter-phase difference vector (delta ln eta and
// each delta x_i) changes sign while both states pass the coincidence predicate. That is the
// pair's coalescence (owner decision 2026-09-16); zero components carry no sign information.
std::optional<std::pair<std::size_t, std::size_t>> swapped_pair(
    const CompiledSystem& compiled, const std::vector<double>& before, const std::vector<double>& after
) {
    const auto a = signatures(compiled, before);
    const auto b = signatures(compiled, after);
    for (std::size_t left = 0; left + 1 < a.size(); ++left) {
        for (std::size_t right = left + 1; right < a.size(); ++right) {
            std::vector<double> before_delta{a[left].log_eta - a[right].log_eta};
            std::vector<double> after_delta{b[left].log_eta - b[right].log_eta};
            for (std::size_t i = 0; i < a[left].composition.size(); ++i) {
                before_delta.push_back(a[left].composition[i] - a[right].composition[i]);
                after_delta.push_back(b[left].composition[i] - b[right].composition[i]);
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

}  // namespace

// Corrected-phase branch identity: the corrected density and the anchored locally unique
// root must share the anchor's window, with a positive corrected pressure slope. The physical
// state is recomputed from the coordinates, never taken from the supplied result's arrays.
BranchCheck check_density_branches(
    const HelmholtzModel& model, const CompiledSystem& compiled, const EquilibriumResult& result,
    const std::vector<double>& anchors, bool initial, core::Scope& scope
) {
    BranchCheck outcome;
    const PhasePoint prepared = compiled.prepare(result.coordinates);
    for (std::size_t phase = 0; phase < compiled.problem.phases.size(); ++phase) {
        const auto& spec = compiled.problem.phases[phase];
        std::vector<double> x = prepared.composition[phase];
        const double total = std::accumulate(x.begin(), x.end(), 0.0);
        for (double& value : x) value /= total;  // the driver's own normalisation of the log chart
        const double corrected = 1.0 / prepared.volume[phase];
        const auto branch = spec.kind == PhaseKind::Vapor ? DensityBranch::Vapor : DensityBranch::Liquid;
        const auto state = model.state(prepared.temperature, corrected, x);
        double anchor = anchors.empty() ? 0.0 : anchors[phase];
        DensityRoot root;
        try {
            if (initial) {
                root = solve_density(model, prepared.temperature, prepared.pressure, x, branch, std::nullopt);
                anchor = root.molar_density;
            }
            root = solve_density(model, prepared.temperature, prepared.pressure, x, branch, anchor);
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
        scope.probe("corrected_slope", state.pressure_density_derivative, static_cast<std::uint32_t>(phase));
        scope.probe("density_residual", root.convergence.residual, static_cast<std::uint32_t>(phase));
        const bool ok = std::isfinite(corrected) && corrected > 0.0
            && std::isfinite(state.pressure_density_derivative) && state.pressure_density_derivative > 0.0
            && root.convergence.residual <= policy::DENSITY_ROOT_RESIDUAL.value
            && displacement <= policy::DENSITY_ANCHOR_LOG_WINDOW.value
            && selected_displacement <= policy::DENSITY_ANCHOR_LOG_WINDOW.value;
        if (!ok) {
            outcome.terminal = true;
            outcome.message = "phase " + spec.name + " corrected density left the anchored branch window";
            return outcome;
        }
    }
    outcome.accepted = true;
    return outcome;
}

namespace {

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

Attempt admit_result(
    const HelmholtzModel& model, const CompiledSystem& compiled, EquilibriumResult result,
    const std::vector<double>& anchors, bool initial, core::Scope& scope
) {
    Attempt outcome;
    outcome.convergence = result.convergence;
    scope.probe("driver_result_present", 1.0);
    if (std::isfinite(result.convergence.residual)) scope.probe("row_residual", result.convergence.residual);
    if (!result.success || !result.has_state()) {
        outcome.stop = stop_of(result, outcome.terminal);
        outcome.message = result.message;
        outcome.result = std::move(result);
        return outcome;
    }
    const double recomputed = compiled.classified_residual_norm(result.coordinates);
    if (!(result.requested_tolerance_met.value_or(false)) || !(recomputed <= compiled.problem.row_tolerance)) {
        outcome.stop = ContinuationStop::stalled;
        outcome.message = "available state misses the caller row tolerance: " + result.message;
        outcome.result = std::move(result);
        return outcome;
    }
    const auto topology = check_topology(compiled, result.coordinates);
    if (topology.event != TopologyEvent::None) {
        result.topology_event = topology.event;
        result.topology_phase = topology.phase;
        result.topology_pair = topology.pair;
        outcome.stop = stop_of(result, outcome.terminal);
        outcome.message = "supplied state fails the candidate topology predicate";
        outcome.result = std::move(result);
        return outcome;
    }
    const auto branches = check_density_branches(model, compiled, result, anchors, initial, scope);
    outcome.accepted = branches.accepted;
    outcome.terminal = branches.terminal;
    outcome.message = branches.message;
    if (!branches.accepted) outcome.stop = ContinuationStop::branch_lost;
    outcome.result = std::move(result);
    return outcome;
}

ContinuationResult finish(ContinuationResult result, ContinuationStop stop, std::string message) {
    result.stop = stop;
    result.success = stop == ContinuationStop::reached_target;
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
        result.event_phase = attempt.result->topology_phase;
        result.event_pair = attempt.result->topology_pair;
    }
}

}  // namespace

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
    }
    return "unknown";
}

ContinuationResult continue_equilibrium(const HelmholtzModel& model, const ContinuationRequest& request) {
    core::Scope scope("continue_equilibrium");
    ContinuationResult result;
    const auto& options = request.options;
    if (!std::isfinite(options.initial_step) || options.initial_step <= 0.0 || options.initial_step > 1.0
        || !std::isfinite(options.minimum_step) || options.minimum_step <= 0.0
        || options.minimum_step > options.initial_step || options.max_steps == 0) {
        throw std::invalid_argument("invalid continuation options");
    }
    const std::string declared = declaration_difference(request.start_problem, request.target_problem);
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
        return std::make_unique<HelmholtzModel>(std::move(parameters), model.active_parameters(), model.thermochemistry());
    };
    auto target_model = parameter_model(1.0);
    CompiledSystem target = compile_problem(target_model ? *target_model : model, request.target_problem);
    const std::string layout = layout_difference(start, target);
    if (!layout.empty()) return finish(std::move(result), ContinuationStop::invalid_path, layout);
    if (!start.feasible || !target.feasible)
        return finish(std::move(result), ContinuationStop::invalid_path, "endpoint is not full rank at its generic point");
    Path path;
    path.base = request.start_problem;
    path.feed_start = start.n_feed;
    path.feed_target = target.n_feed;
    path.t0 = request.start_problem.temperature.value;
    path.t1 = request.target_problem.temperature.value;
    path.p0 = request.start_problem.pressure.value;
    path.p1 = request.target_problem.pressure.value;
    const bool moves = path.feed_start != path.feed_target
        || (!path.base.temperature.free && path.t0 != path.t1)
        || (!path.base.pressure.free && path.p0 != path.p1);
    if (parameter_path && moves)
        return finish(std::move(result), ContinuationStop::invalid_path, "parameter path must hold intensives and resolved feed fixed");
    if (!moves && (!parameter_path || parameter_start == request.parameter_targets))
        return finish(std::move(result), ContinuationStop::invalid_path,
            parameter_path ? "no active parameter differs" : "no fixed intensive or resolved feed differs");
    if (parameter_path) { path.parameter_start = parameter_start; path.parameter_target = request.parameter_targets; }

    // Initial admission at q = 0.
    std::vector<double> anchors;
    std::vector<double> y_current;
    CompiledSystem current = start;
    std::unique_ptr<HelmholtzModel> current_model;
    {
        core::Scope initial("equilibrium_continuation_initial");
        initial.probe("path_value", 0.0);
        Attempt admitted;
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
                ? solve_equilibrium_from_coordinates(model, request.start_problem, supplied.coordinates) : supplied, anchors, true, initial);
        } else {
            admitted = admit_result(model, start,
                solve_equilibrium_from_coordinates(model, request.start_problem, boundary_start(start)),
                anchors, true, initial);
        }
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
    while (true) {
        if (q >= 1.0) return finish(std::move(result), ContinuationStop::reached_target, "");
        if (result.accepted_steps >= options.max_steps)
            return finish(std::move(result), ContinuationStop::budget_exhausted, "max_steps reached before the target");
        const double q_try = std::min(q + dq, 1.0);
        ++result.attempts;
        core::Scope attempt_scope("equilibrium_continuation_attempt");
        attempt_scope.probe("path_value", q_try);
        Attempt attempt;
        const EquilibriumProblem problem = path.at(q_try);
        auto trial_model = parameter_model(q_try);
        const auto& trial = trial_model ? *trial_model : model;
        CompiledSystem compiled = compile_problem(trial, problem);
        const std::string generated = layout_difference(start, compiled);
        if (!generated.empty()) {
            attempt.stop = ContinuationStop::invalid_path;
            attempt.terminal = true;
            attempt.message = generated;
        } else if (!compiled.feasible) {
            // Same layout, generic-point rank loss at the generated problem.
            attempt.stop = ContinuationStop::singular_jacobian;
            attempt.terminal = true;
            attempt.message = "generated problem loses generic rank (nullity "
                + std::to_string(compiled.generic_nullity) + ")";
        } else {
            const auto jacobian = current.jacobian(y_current);
            const auto tangent = sensitivity(current, y_current, path.lever(2 + model.parameters().component_count() + model.active_parameters().size()));
            if (tangent.empty()) {
                // A state-local tangent cannot change by halving the step: singular is terminal
                // (owner decision 2026-09-16), the branch end at that accepted state.
                const bool numeric = std::all_of(jacobian.begin(), jacobian.end(), [](double v) { return std::isfinite(v); });
                const auto info = inspect_sensitivity(current, y_current, {});
                attempt.stop = numeric ? ContinuationStop::singular_jacobian : ContinuationStop::nonfinite;
                attempt.terminal = true;
                attempt.message = numeric
                    ? "implicit tangent is unavailable at the accepted state (rank "
                        + std::to_string(info.rank) + " of " + std::to_string(y_current.size())
                        + ", condition estimate " + std::to_string(info.condition_estimate) + ")"
                    : "nonfinite operands in the implicit tangent";
            } else {
                std::vector<double> predicted = y_current;
                for (std::size_t i = 0; i < predicted.size(); ++i) predicted[i] += (q_try - q) * tangent[i];
                std::size_t projected = 0;
                const bool inside = inside_domain(compiled, predicted, projected);
                attempt_scope.probe("predictor_projected", static_cast<double>(projected));
                if (!inside) {
                    attempt.stop = ContinuationStop::left_domain;
                    attempt.message = "predicted coordinates leave the admitted domain";
                } else {
                    attempt = admit_result(trial, compiled, solve_equilibrium_from_coordinates(trial, problem, predicted),
                                           anchors, false, attempt_scope);
                    if (attempt.accepted) {
                        if (const auto pair = swapped_pair(compiled, y_current, attempt.result->coordinates)) {
                            attempt.accepted = false;
                            attempt.terminal = true;
                            attempt.stop = ContinuationStop::phases_coalesced;
                            attempt.result->topology_event = TopologyEvent::PhasesCoalesced;
                            attempt.result->topology_pair = pair;
                            attempt.message = "phases " + compiled.problem.phases[pair->first].name + " and "
                                + compiled.problem.phases[pair->second].name
                                + " swapped identity between accepted states: passage through a critical point";
                        }
                    }
                }
            }
        }
        if (!attempt.result.has_value()) attempt_scope.probe("driver_result_present", 0.0);
        attempt_scope.probe("accepted", attempt.accepted ? 1.0 : 0.0);
        attempt_scope.probe(continuation_stop_name(attempt.accepted ? ContinuationStop::reached_target : attempt.stop), 1.0);
        if (attempt.convergence.has_value()) attempt_scope.finish(*attempt.convergence);
        if (attempt.accepted) {
            q = q_try;
            ++result.accepted_steps;
            local_retries = 0;
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

}  // namespace epcsaft::equilibrium
