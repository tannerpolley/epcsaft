#include "equilibrium/gibbs_minimizer.hpp"

#include "core/numerical_policy.hpp"
#include "core/trace.hpp"
#include "equilibrium/stoichiometry.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
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

template <class Value>
bool same_optional(const std::optional<Value>& left, const std::optional<Value>& right) {
    return left.has_value() == right.has_value() && (!left.has_value() || *left == *right);
}

bool same_correlation(const ReactionLogPolynomial& a, const ReactionLogPolynomial& b) {
    return a.a == b.a && a.b == b.b && a.c == b.c && a.d == b.d
        && a.reference_temperature == b.reference_temperature
        && a.temperature_min == b.temperature_min && a.temperature_max == b.temperature_max
        && a.standard_state_id == b.standard_state_id;
}

bool same_intensive(const FreeValue& a, const FreeValue& b) {
    if (a.free != b.free) return false;
    return !a.free || (a.absent == b.absent && a.value == b.value);
}

// Every field that is not a path lever must be identical at both endpoints.
std::string declaration_difference(const EquilibriumProblem& a, const EquilibriumProblem& b) {
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
        const bool same_law = x.name == y.name && x.stoichiometry == y.stoichiometry
            && same_optional(x.equilibrium_constant, y.equilibrium_constant)
            && same_optional(x.affinity_j_per_mol, y.affinity_j_per_mol)
            && x.correlation.has_value() == y.correlation.has_value()
            && (!x.correlation.has_value() || same_correlation(*x.correlation, *y.correlation));
        if (!same_law) return "reaction " + reaction_name(x, r) + " differs";
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
    const auto phase_offsets = offsets(compiled);
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

// Corrected-phase branch identity: the corrected density and the anchored locally unique
// root must share the anchor's window, with a positive corrected pressure slope. The physical
// state is recomputed from the coordinates, never taken from the supplied result's arrays.
Attempt admit_branches(
    const HelmholtzModel& model, const CompiledSystem& compiled, const EquilibriumResult& result,
    const std::vector<double>& anchors, bool initial, core::Scope& scope
) {
    Attempt outcome;
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
            outcome.stop = ContinuationStop::branch_lost;
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
            outcome.stop = ContinuationStop::branch_lost;
            outcome.terminal = true;
            outcome.message = "phase " + spec.name + " corrected density left the anchored branch window";
            return outcome;
        }
    }
    outcome.accepted = true;
    return outcome;
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
    outcome = admit_branches(model, compiled, result, anchors, initial, scope);
    outcome.convergence = result.convergence;
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
    CompiledSystem target = compile_problem(model, request.target_problem);
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
    if (!moves) return finish(std::move(result), ContinuationStop::invalid_path, "no fixed intensive or resolved feed differs");

    // Initial admission at q = 0.
    std::vector<double> anchors;
    std::vector<double> y_current;
    CompiledSystem current = start;
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
            admitted = admit_result(model, start, supplied, anchors, true, initial);
        } else {
            admitted = admit_result(model, start, solve_equilibrium(model, request.start_problem), anchors, true, initial);
        }
        initial.probe("accepted", admitted.accepted ? 1.0 : 0.0);
        initial.probe(continuation_stop_name(admitted.accepted ? ContinuationStop::reached_target : admitted.stop), 1.0);
        if (admitted.convergence.has_value()) initial.finish(*admitted.convergence);
        if (!admitted.accepted) {
            reject(result, admitted, 0.0);
            return finish(std::move(result), ContinuationStop::initial_failure, admitted.message);
        }
        result.path_values.push_back(0.0);
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
        CompiledSystem compiled = compile_problem(model, problem);
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
            const auto tangent = sensitivity(current, y_current, path.lever(current.lever_jacobian(y_current).size() / current.rows.size()));
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
                    attempt = admit_result(model, compiled, solve_equilibrium_from_coordinates(model, problem, predicted),
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
