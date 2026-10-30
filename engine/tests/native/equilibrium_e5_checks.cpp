// E5: solved-state actions against stencils of fresh nonlinear solves.
#include "equilibrium_checks_shared.hpp"
#ifdef _WIN32
// Exercise the same x87 control word with the selected Microsoft-targeting Clang.
using fpu_control_t = unsigned short;
#define _FPU_GETCW(word) __asm__ volatile("fnstcw %0" : "=m" (word))
#define _FPU_SETCW(word) __asm__ volatile("fldcw %0" : : "m" (word))
#else
#include <fpu_control.h>
#endif

namespace equilibrium_checks {

// ------------------------------------------------------------------ E5 solved-state actions

namespace e5 {

using epcsaft::PropertyObservable;
using epcsaft::equilibrium::AmountSpec;
using epcsaft::equilibrium::CompiledSystem;
using epcsaft::equilibrium::EquilibriumProblem;
using epcsaft::equilibrium::EquilibriumResult;
using epcsaft::equilibrium::FreeValue;
using epcsaft::equilibrium::PhaseKind;
using epcsaft::equilibrium::SolvedStateActionBatch;
using epcsaft::equilibrium::SolvedStateActionDirection;
using epcsaft::equilibrium::SolvedStateActionRequest;
using epcsaft::equilibrium::SolvedStateActionStatus;
using epcsaft::equilibrium::SolvedStateObservable;
using Kind = epcsaft::equilibrium::SolvedStateObservableKind;

constexpr double LADDER[3][3] = {{1.0e-2, 3.0e-3, 1.0e-3}, {3.0e-2, 1.0e-2, 3.0e-3}, {1.0e-1, 3.0e-2, 1.0e-2}};
constexpr double CRITERION[3] = {5.0e-6, 5.0e-5, 5.0e-4};
constexpr double STRUCTURAL_ZERO = 1.0e-12;

SolvedStateObservable obs(Kind kind, std::size_t phase = 0, std::size_t component = 0,
                          PropertyObservable property = PropertyObservable::IdealEnthalpy) {
    SolvedStateObservable result;
    result.kind = kind; result.phase = phase; result.component = component; result.property = property;
    return result;
}

SolvedStateActionDirection direction(double T, double P, std::vector<double> feed, std::vector<double> active = {}) {
    SolvedStateActionDirection result;
    result.temperature_k = T; result.pressure_pa = P;
    result.feed_amounts_mol = std::move(feed); result.active_parameters = std::move(active);
    return result;
}

// One qualification fixture: a model recipe, a problem, three directions and observables.
struct Fixture final {
    std::string label;
    epcsaft::Parameters parameters;
    std::vector<epcsaft::ActiveParameter> active;
    std::optional<epcsaft::ThermochemistryRecord> record;
    EquilibriumProblem problem;
    std::vector<SolvedStateActionDirection> directions;
    std::vector<SolvedStateObservable> observables;
    bool scaled = true;
    bool constant_ideal_cp = false; // declared constant-cp thermal reference, h is affine in T
    // Corner branch admission for a pure declared phase (case 6) and an analytic corner
    // reference for selected observables (the A->2B pure-state relation).
    std::optional<epcsaft::DensityBranch> corner_branch;
    std::function<std::optional<double>(const SolvedStateObservable&, double)> analytic;

    // A fresh model with every active parameter shifted (reciprocal pairs move together).
    epcsaft::HelmholtzModel model(const std::vector<double>& shift = {}) const {
        auto p = parameters;
        for (std::size_t a = 0; a < active.size() && a < shift.size(); ++a) {
            std::vector<double>& family = active[a].family == epcsaft::ParameterFamily::SegmentCount ? p.segment_count
                : active[a].family == epcsaft::ParameterFamily::SigmaConstant ? p.sigma_constant
                : active[a].family == epcsaft::ParameterFamily::DispersionEnergyOverK ? p.dispersion_energy_over_k
                : active[a].family == epcsaft::ParameterFamily::Kij ? p.k_ij
                : active[a].family == epcsaft::ParameterFamily::Lij ? p.l_ij
                : active[a].family == epcsaft::ParameterFamily::AssociationEnergyOverK ? p.association_energy_over_k
                : p.association_volume;
            for (const std::size_t index : active[a].indices) family[index] += shift[a];
        }
        return epcsaft::HelmholtzModel(p, active, record);
    }
};

// Every geometry, phase-pressure and residual-potential identity on the supported components, plus T and P.
std::vector<SolvedStateObservable> geometry_observables(const EquilibriumProblem& problem, std::size_t components) {
    std::vector<SolvedStateObservable> result = {obs(Kind::Temperature), obs(Kind::Pressure)};
    for (std::size_t phase = 0; phase < problem.phases.size(); ++phase) {
        for (const Kind kind : {Kind::PhaseAmount, Kind::PhaseVolume, Kind::PhaseMolarDensity, Kind::PhasePackingFraction, Kind::PhasePressure})
            result.push_back(obs(kind, phase));
        for (std::size_t i = 0; i < components; ++i) {
            result.push_back(obs(Kind::PhaseComponentAmount, phase, i));
            result.push_back(obs(Kind::PhaseMoleFraction, phase, i));
            const auto& support = problem.phases[phase].support;
            if (support.empty() || std::find(support.begin(), support.end(), i) != support.end()) {
                result.push_back(obs(Kind::PhaseResidualChemicalPotentialOverRt, phase, i));
                result.push_back(obs(Kind::PhaseLogActivity, phase, i));
            }
        }
    }
    return result;
}

// The physical observable map of an accepted result through the public owners; never E5.
std::optional<double> mapped(const SolvedStateObservable& item, const epcsaft::HelmholtzModel& model,
                             const EquilibriumResult& result) {
    const std::size_t components = model.parameters().component_count();
    const auto x = [&](std::size_t phase) {
        return std::vector<double>(result.mole_fractions.begin() + static_cast<std::ptrdiff_t>(phase * components),
                                   result.mole_fractions.begin() + static_cast<std::ptrdiff_t>((phase + 1) * components));
    };
    const auto state = [&](std::size_t phase) {
        return model.state(*result.temperature, result.molar_densities[phase], x(phase));
    };
    switch (item.kind) {
    case Kind::Temperature: return *result.temperature;
    case Kind::Pressure: return *result.pressure;
    case Kind::PhaseAmount: return result.amounts[item.phase];
    case Kind::PhaseComponentAmount: return result.amounts[item.phase] * x(item.phase)[item.component];
    case Kind::PhaseMoleFraction: return x(item.phase)[item.component];
    case Kind::PhaseVolume: return result.amounts[item.phase] / result.molar_densities[item.phase];
    case Kind::PhaseMolarDensity: return result.molar_densities[item.phase];
    case Kind::PhasePackingFraction: return result.packing_fractions[item.phase];
    case Kind::PhasePressure: return state(item.phase).pressure;
    case Kind::PhaseResidualChemicalPotentialOverRt: return state(item.phase).residual_chemical_potential_over_rt[item.component].value();
    case Kind::PhaseLogActivity:  // ln(x rho / (1 mol/m^3)) + mu_res/RT through the public State
        return std::log(x(item.phase)[item.component] * result.molar_densities[item.phase])
            + state(item.phase).residual_chemical_potential_over_rt[item.component].value();
    case Kind::PhaseProperty: {
        const auto s = state(item.phase);
        switch (item.property) {
        case PropertyObservable::IdealEnthalpy: return s.ideal_enthalpy;
        case PropertyObservable::IdealInternalEnergy: return s.ideal_internal_energy;
        case PropertyObservable::IdealEntropy: return s.ideal_entropy;
        case PropertyObservable::IdealIsobaricHeatCapacity: return s.ideal_isobaric_heat_capacity;
        case PropertyObservable::IdealIsochoricHeatCapacity: return s.ideal_isochoric_heat_capacity;
        case PropertyObservable::ResidualEnthalpy: return s.residual_enthalpy;
        case PropertyObservable::ResidualInternalEnergy: return s.residual_internal_energy;
        case PropertyObservable::ResidualEntropy: return s.residual_entropy;
        case PropertyObservable::ResidualIsobaricHeatCapacity: return s.residual_isobaric_heat_capacity;
        case PropertyObservable::ResidualIsochoricHeatCapacity: return s.residual_isochoric_heat_capacity;
        case PropertyObservable::TotalEnthalpy: return s.total_enthalpy;
        case PropertyObservable::TotalInternalEnergy: return s.total_internal_energy;
        case PropertyObservable::TotalEntropy: return s.total_entropy;
        case PropertyObservable::TotalIsobaricHeatCapacity: return s.total_isobaric_heat_capacity;
        case PropertyObservable::TotalIsochoricHeatCapacity: return s.total_isochoric_heat_capacity;
        }
        return std::nullopt;
    }
    }
    return std::nullopt;
}

// The sole test-only stencil scaling: max(|q_i|, one physical unit) per external coordinate.
double magnitude(const SolvedStateActionDirection& d, const CompiledSystem& compiled, const std::vector<double>& active_values) {
    double result = std::max(std::abs(d.temperature_k) / std::max(compiled.problem.temperature.value, 1.0),
                             std::abs(d.pressure_pa) / std::max(compiled.problem.pressure.value, 1.0));
    for (std::size_t i = 0; i < d.feed_amounts_mol.size(); ++i)
        result = std::max(result, std::abs(d.feed_amounts_mol[i]) / std::max(std::abs(compiled.n_feed[i]), 1.0));
    for (std::size_t a = 0; a < d.active_parameters.size(); ++a)
        result = std::max(result, std::abs(d.active_parameters[a]) / std::max(std::abs(active_values[a]), 1.0));
    return result;
}

std::vector<double> active_values(const Fixture& fixture) {
    std::vector<double> result;
    for (const auto& item : fixture.active) {
        const auto& p = fixture.parameters;
        const std::vector<double>& family = item.family == epcsaft::ParameterFamily::SegmentCount ? p.segment_count
            : item.family == epcsaft::ParameterFamily::SigmaConstant ? p.sigma_constant
            : item.family == epcsaft::ParameterFamily::DispersionEnergyOverK ? p.dispersion_energy_over_k
            : item.family == epcsaft::ParameterFamily::Kij ? p.k_ij : item.family == epcsaft::ParameterFamily::Lij ? p.l_ij
            : item.family == epcsaft::ParameterFamily::AssociationEnergyOverK ? p.association_energy_over_k
            : p.association_volume;
        result.push_back(family[item.indices.front()]);
    }
    return result;
}

// Independent reference: a fresh nonlinear solve at every stencil corner, mapped physically.
struct Reference final {
    std::vector<std::optional<double>> value;   // per observable, absent when a corner was rejected
    std::vector<double> rounding;              // rounding propagated through the stencil
    std::vector<std::string> rejected;          // corner failure names
};

Reference reference(const Fixture& fixture, const CompiledSystem& compiled, const EquilibriumResult& central,
                    std::size_t order, double multiplier) {
    const std::size_t components = compiled.model->parameters().component_count();
    const auto values = active_values(fixture);
    std::vector<double> scale(order, 1.0);
    for (std::size_t j = 0; j < order; ++j)
        if (fixture.scaled) scale[j] = magnitude(fixture.directions[j], compiled, values);
    // Keep the derivative reference independent of the multiplier-transfer method under test.
    auto reference_start = central;
    reference_start.optimizer_state.reset();
    Reference result;
    result.value.assign(fixture.observables.size(), 0.0);
    result.rounding.assign(fixture.observables.size(), 0.0);
    double product = 1.0;
    for (const double m : scale) product *= m;
    if (product == 0.0) return result;
    const double limit = std::min(fixture.problem.solver_options.value(epcsaft::policy::row), epcsaft::policy::EQUILIBRIUM_ROW_RESIDUAL.value);
    // Higher orders: one centered derivative of the preceding action on fresh solves.
    // First order below still has an independent physical value stencil.
    if (order > 1) {
        const auto& d = fixture.directions[order - 1];
        const double h = multiplier / scale.back();
        for (const double sign : {-1.0, 1.0}) {
            auto shifted = fixture;
            if (!shifted.problem.temperature.free) shifted.problem.temperature.value += sign * h * d.temperature_k;
            if (!shifted.problem.pressure.free) shifted.problem.pressure.value += sign * h * d.pressure_pa;
            for (std::size_t i = 0; i < components; ++i) shifted.problem.feed.values[i] += sign * h * d.feed_amounts_mol[i];
            std::vector<double> delta(fixture.active.size(), 0.0);
            for (std::size_t a = 0; a < delta.size(); ++a) delta[a] = sign * h * d.active_parameters[a];
            const auto model = shifted.model(delta);
            const auto solved = epcsaft::equilibrium::solve_equilibrium(model, shifted.problem, reference_start);
            if (!solved.success || !solved.requested_tolerance_met.value_or(false) || solved.convergence.residual > limit) {
                result.rejected.push_back("preceding-action corner: " + solved.message);
                for (auto& value : result.value) value.reset();
                return result;
            }
            if (fixture.corner_branch.has_value()) {
                const auto root = epcsaft::solve_density(model, shifted.problem.temperature.value,
                    shifted.problem.pressure.value, {1.0}, *fixture.corner_branch, std::nullopt);
                CHECK(std::abs(std::log(solved.molar_densities[0] / root.molar_density)) <= 1.0e-8);
            }
            const auto shifted_compiled = epcsaft::equilibrium::compile_problem(model, shifted.problem);
            SolvedStateActionRequest request{fixture.observables, {}};
            request.directions.assign(fixture.directions.begin(), fixture.directions.begin() + static_cast<std::ptrdiff_t>(order - 1));
            const auto batch = epcsaft::equilibrium::solved_state_actions(shifted_compiled, solved, request);
            for (std::size_t i = 0; i < result.value.size(); ++i) {
                // Preserve the stronger independent manufactured-value oracle
                // where it exists; only its caloric siblings use preceding actions.
                if (fixture.analytic && fixture.analytic(fixture.observables[i], fixture.problem.temperature.value)) continue;
                const auto& action = batch.results[i];
                if (action.status == SolvedStateActionStatus::Available && action.action && result.value[i]) {
                    *result.value[i] += sign * *action.action / (2.0 * h);
                    result.rounding[i] += 32.0 * std::numeric_limits<double>::epsilon() * std::abs(*action.action) / (2.0 * h);
                } else {
                    result.value[i].reset();
                    if (action.status != SolvedStateActionStatus::Available)
                        result.rejected.push_back(std::string("preceding action: ") + epcsaft::equilibrium::solved_state_status_name(action.status)
                            + " " + action.diagnostic_message.value_or("no diagnostic"));
                }
            }
        }
        if (!fixture.analytic) return result;
    }
    for (std::size_t corner = 0; corner < (std::size_t(1) << order); ++corner) {
        auto problem = fixture.problem;
        std::vector<double> shift(fixture.active.size(), 0.0);
        double weight = 1.0;
        for (std::size_t j = 0; j < order; ++j) {
            const double sign = corner & (std::size_t(1) << j) ? -1.0 : 1.0;
            weight *= sign;
            const double h = sign * multiplier / scale[j];
            const auto& d = fixture.directions[j];
            if (!problem.temperature.free) problem.temperature.value += h * d.temperature_k;
            if (!problem.pressure.free) problem.pressure.value += h * d.pressure_pa;
            for (std::size_t i = 0; i < components; ++i) problem.feed.values[i] += h * d.feed_amounts_mol[i];
            for (std::size_t a = 0; a < shift.size(); ++a) shift[a] += h * d.active_parameters[a];
        }
        const auto model = fixture.model(shift);
        const auto solved = epcsaft::equilibrium::solve_equilibrium(model, problem, reference_start);
        std::string failure;
        if (!solved.success) failure = "corner solve: " + solved.message;
        else if (solved.phase_count != central.phase_count) failure = "corner phase count";
        else if (!solved.requested_tolerance_met.value_or(false)) failure = "corner row tolerance";
        else if (!(solved.convergence.residual <= limit)) failure = "corner classified rows";
        else if (fixture.corner_branch.has_value()) {
            const auto root = epcsaft::solve_density(model, problem.temperature.value, problem.pressure.value, {1.0}, *fixture.corner_branch, std::nullopt);
            if (!(std::abs(std::log(solved.molar_densities[0] / root.molar_density)) <= 1.0e-8)) failure = "corner branch";
        }
        if (!failure.empty()) {
            result.rejected.push_back(failure);
            for (auto& value : result.value) value.reset();
            return result;
        }
        for (std::size_t i = 0; i < fixture.observables.size(); ++i) {
            const auto analytic = fixture.analytic ? fixture.analytic(fixture.observables[i], problem.temperature.value) : std::nullopt;
            if (order > 1 && !analytic.has_value()) continue;
            const auto value = analytic.has_value() ? analytic : mapped(fixture.observables[i], model, solved);
            if (value.has_value() && result.value[i].has_value()) {
                *result.value[i] += weight * *value;
                result.rounding[i] += 32.0 * std::numeric_limits<double>::epsilon() * std::abs(*value)
                    * product / std::pow(2.0 * multiplier, static_cast<double>(order));
            }
            else result.value[i].reset();
        }
    }
    for (std::size_t i = 0; i < result.value.size(); ++i)
        if (result.value[i] && (order == 1 || (fixture.analytic
            && fixture.analytic(fixture.observables[i], fixture.problem.temperature.value))))
            *result.value[i] *= product / std::pow(2.0 * multiplier, static_cast<double>(order));
    return result;
}

std::string status_field(const SolvedStateActionBatch& batch, std::size_t index) {
    return epcsaft::equilibrium::solved_state_status_name(batch.results[index].status);
}

// Sum the fixed chart-coordinate sub-actions; a cancelling sum cannot set its own scale.
std::vector<double> action_term_scale(const Fixture& fixture, const CompiledSystem& compiled,
                                     const EquilibriumResult& central, std::size_t order) {
    std::vector<double> result(fixture.observables.size(), 0.0);
    SolvedStateActionRequest request{fixture.observables, {}};
    const auto expand = [&](const auto& self, std::size_t j) -> void {
        if (j == order) {
            const auto batch = epcsaft::equilibrium::solved_state_actions(compiled, central, request);
            for (std::size_t i = 0; i < result.size(); ++i) {
                CHECK(batch.results[i].status == SolvedStateActionStatus::Available && batch.results[i].action);
                result[i] += std::abs(*batch.results[i].action);
            }
            return;
        }
        const auto& full = fixture.directions[j];
        for (int part = 0; part < 4; ++part) {
            auto d = full;
            if (part != 0) d.temperature_k = 0.0;
            if (part != 1) d.pressure_pa = 0.0;
            if (part != 2) std::fill(d.feed_amounts_mol.begin(), d.feed_amounts_mol.end(), 0.0);
            if (part != 3) {
                std::fill(d.active_parameters.begin(), d.active_parameters.end(), 0.0);
                for (auto& coefficients : d.reaction_coefficients) coefficients.fill(0.0);
            }
            bool nonzero = d.temperature_k != 0.0 || d.pressure_pa != 0.0;
            for (const auto* values : {&d.feed_amounts_mol, &d.active_parameters})
                for (const auto value : *values) nonzero = nonzero || value != 0.0;
            for (const auto& coefficients : d.reaction_coefficients)
                for (double value : coefficients) nonzero |= value != 0.0;
            if (nonzero) { request.directions.push_back(d); self(self, j + 1); request.directions.pop_back(); }
        }
    };
    expand(expand, 0);
    return result;
}

// Qualify one fixture through orders 1 to min(orders, 3); every observable/order/step is retained.
// analyses/2026-reactive-caloric-actions/scripts/main-2026-10-06/surface-equilibrium.cpp passes orders = 1.
bool qualify(const Fixture& fixture, const char* group, std::size_t orders = 3) {
    const auto model = fixture.model();
    const auto compiled = epcsaft::equilibrium::compile_problem(model, fixture.problem);
    const auto central = epcsaft::equilibrium::solve_equilibrium(model, fixture.problem);
    CHECK(central.success);
    for (std::size_t order = 1; order <= std::min<std::size_t>(orders, 3); ++order) {
        SolvedStateActionRequest request;
        request.observables = fixture.observables;
        request.directions.assign(fixture.directions.begin(), fixture.directions.begin() + static_cast<std::ptrdiff_t>(order));
        const auto batch = epcsaft::equilibrium::solved_state_actions(compiled, central, request);
        const auto scales = action_term_scale(fixture, compiled, central, order);
        // The ionic cp A2 stencil needs larger nodes to resolve its ~1e-7 A3
        // above fresh-solve rounding, without changing the action criterion.
        const std::array<double, 3> ionic_third_steps{2.0, 1.0, 0.5};
        // Normalized pinned/reactive directions need half-sized A3 nodes:
        // their coarse Richardson truncation exceeds c despite fine agreement.
        const std::array<double, 3> scaled_third_steps{0.05, 0.015, 0.005};
        // A->2B physical-value A3 needs larger nodes: the original propagated rounding exceeds c.
        const std::array<double, 3> doubling_third_steps{0.2, 0.06, 0.02};
        const double* steps = order == 3 && fixture.label == "a-to-2b" ? doubling_third_steps.data()
            : order == 3 && fixture.label.rfind("electrolyte", 0) == 0 ? ionic_third_steps.data()
            : order == 3 && fixture.scaled ? scaled_third_steps.data() : LADDER[order - 1];
        std::vector<Reference> ladder;
        for (std::size_t s = 0; s < 3; ++s) ladder.push_back(reference(fixture, compiled, central, order, steps[s]));
        for (std::size_t i = 0; i < fixture.observables.size(); ++i) {
            const auto& item = batch.results[i];
            // These are defining affine relations, not small nonzero actions:
            // p_phase=P, conserved nonreactive phase amounts with at most one
            // free amount, and the declared constant-cp ideal enthalpy.
            const bool affine_pressure = order > 1 && item.observable.kind == Kind::PhasePressure
                && !fixture.problem.pressure.free;
            const auto free_amounts = std::count_if(fixture.problem.phases.begin(), fixture.problem.phases.end(),
                [](const auto& phase) { return phase.amount.kind == epcsaft::equilibrium::AmountKind::Free; });
            const auto occupied_phases = std::count_if(fixture.problem.phases.begin(), fixture.problem.phases.end(),
                [](const auto& phase) { return phase.amount.kind == epcsaft::equilibrium::AmountKind::Free
                    || phase.amount.fraction > 0.0; });
            const bool affine_amount = order > 1 && fixture.problem.reactions.empty()
                && ((item.observable.kind == Kind::PhaseAmount && (free_amounts <= 1
                        || fixture.problem.phases[item.observable.phase].amount.kind == epcsaft::equilibrium::AmountKind::Pinned))
                    || (item.observable.kind == Kind::PhaseComponentAmount && occupied_phases == 1));
            const bool affine_ideal = order > 1 && fixture.constant_ideal_cp && item.observable.kind == Kind::PhaseProperty
                && item.observable.property == PropertyObservable::IdealEnthalpy;
            const bool constant_amount = order == 1 && fixture.problem.reactions.empty()
                && item.observable.kind == Kind::PhaseAmount
                && (free_amounts <= 1 || fixture.problem.phases[item.observable.phase].amount.kind == epcsaft::equilibrium::AmountKind::Pinned)
                && std::accumulate(fixture.directions[0].feed_amounts_mol.begin(), fixture.directions[0].feed_amounts_mol.end(), 0.0) == 0.0;
            const bool affine = affine_pressure || affine_amount || affine_ideal || constant_amount;
            if (affine) for (auto& rung : ladder) if (rung.value[i]) { rung.value[i] = 0.0; rung.rounding[i] = 0.0; }
            const auto id = std::string(group) + "/" + fixture.label + "/" + epcsaft::equilibrium::solved_state_observable_id(item.observable)
                + "/A" + std::to_string(order);
            const auto inputs = recording_support::fields({
                recording_support::field("case", fixture.label), recording_support::field("order", order),
                recording_support::field("status", status_field(batch, i)), recording_support::field("unit", item.unit),
                recording_support::field("value", epcsaft::core::optional_number(item.value)),
                recording_support::field("scaled", fixture.scaled ? 1U : 0U),
                recording_support::field("ladder", std::vector<double>(steps, steps + 3)),
                recording_support::field("rank", batch.diagnostics.rank.has_value() ? std::to_string(*batch.diagnostics.rank) : "absent"),
                recording_support::field("condition", epcsaft::core::optional_number(batch.diagnostics.condition_estimate)),
                recording_support::field("linear_residual", epcsaft::core::optional_number(batch.diagnostics.maximum_linear_residual)),
                recording_support::field("central_residual", epcsaft::core::optional_number(batch.diagnostics.central_residual)),
                recording_support::field("unscaled_residual", central.residual_norm()),
                recording_support::field("energy_order", batch.diagnostics.maximum_energy_order),
                recording_support::field("rejected_corners", ladder.back().rejected.empty() ? "none" : ladder.back().rejected.front())});
            const double action = item.action.value_or(std::numeric_limits<double>::quiet_NaN());
            const double nan = std::numeric_limits<double>::quiet_NaN();
            for (std::size_t step = 0; step + 1 < ladder.size(); ++step) {
                recording_support::observe((id + "/step" + std::to_string(step)).c_str(), "action", item.unit.c_str(), action,
                    ladder[step].value[i], inputs, "Relationship", "centered-solve-ladder",
                    "coarser reference step retained for the rate record", group, id.c_str());
            }
            // Richardson reference with rate certification (owner decision 2026-09-16): the
            // centered references scale as h^2, so the last two steps extrapolate to the
            // reference; the three steps certify the rate, or already agree within 2c.
            const auto& coarse = ladder[0].value[i];
            const auto& mid = ladder[1].value[i];
            const auto& fine = ladder[2].value[i];
            const double h1 = steps[0], h2 = steps[1], h3 = steps[2];
            const bool structural = fine.has_value() && mid.has_value() && *fine == 0.0 && *mid == 0.0;
            const bool complete = fine.has_value() && mid.has_value() && coarse.has_value();
            // A finest step whose difference grows again is rounding-bounded (the corner
            // solves agree to ~1e-15 relative); its ladder extrapolates from the two coarser
            // steps and certifies nothing beyond their agreement (Execution review, 2026-09-16).
            const bool noise_bounded = complete && std::abs(*mid - *fine) > std::abs(*coarse - *mid);
            std::optional<double> extrapolated;
            if (complete && noise_bounded) extrapolated = *mid + (*mid - *coarse) / (h1 * h1 / (h2 * h2) - 1.0);
            else if (fine.has_value() && mid.has_value()) extrapolated = *fine + (*fine - *mid) / (h2 * h2 / (h3 * h3) - 1.0);
            // Pressure is checked in reduced units P/P*, amount on the feed
            // basis; retain the existing dimensionless structural-zero criterion.
            const double affine_scale = affine_pressure ? fixture.problem.pressure.value
                : affine_amount || constant_amount ? compiled.feed_total : 1.0;
            const double criterion = structural ? STRUCTURAL_ZERO * affine_scale : CRITERION[order - 1] * scales[i];
            const double ratio_coarse = h1 * h1 / (h2 * h2), ratio_fine = h2 * h2 / (h3 * h3);
            double reference_uncertainty = nan;
            if (complete) {
                const double coarser = *mid + (*mid - *coarse) / (ratio_coarse - 1.0);
                const double finer = *fine + (*fine - *mid) / (ratio_fine - 1.0);
                reference_uncertainty = noise_bounded
                    ? std::abs(*coarse - *mid) / (ratio_coarse - 1.0)
                        + (ratio_coarse * ladder[1].rounding[i] + ladder[0].rounding[i]) / (ratio_coarse - 1.0)
                    : std::abs(finer - coarser)
                        + (ratio_fine * ladder[2].rounding[i] + ladder[1].rounding[i]) / (ratio_fine - 1.0);
                if (structural && *coarse == 0.0) reference_uncertainty = 0.0;
            }
            const bool resolved = reference_uncertainty <= criterion;
            const bool cp_third = fixture.label.rfind("electrolyte", 0) == 0 && order == 3 && item.observable.kind == Kind::PhaseProperty
                && (item.observable.property == PropertyObservable::ResidualIsobaricHeatCapacity
                    || item.observable.property == PropertyObservable::TotalIsobaricHeatCapacity);
            const bool nonzero_signal = !cp_third || std::abs(extrapolated.value_or(0.0)) > criterion;
            const double expected_rate = (h1 * h1 - h2 * h2) / (h2 * h2 - h3 * h3);
            bool refined = complete && (noise_bounded ? std::abs(*coarse - *mid) : std::abs(*fine - *mid)) <= 2.0 * criterion;
            double rate = nan;
            if (complete && !refined && !noise_bounded) {
                rate = (*coarse - *mid) / (*mid - *fine);
                refined = rate >= 0.5 * expected_rate && rate <= 2.0 * expected_rate;
            }
            const bool verdict = item.status == SolvedStateActionStatus::Available && item.action.has_value() && extrapolated.has_value()
                && std::abs(action - *extrapolated) <= criterion && resolved && nonzero_signal;
            // The best-step relative agreement is the falsifying strength independent of the
            // absolute floor: a row with |R| <= c would pass a zero action on the criterion alone.
            double best_step = nan;
            for (const auto* step : {&coarse, &mid, &fine})
                if (step->has_value() && std::isfinite(action))
                    best_step = std::min(std::isfinite(best_step) ? best_step : std::numeric_limits<double>::infinity(),
                                         std::abs(action - **step) / std::max(std::abs(**step), std::numeric_limits<double>::min()));
            const auto ladder_inputs = recording_support::fields({inputs, recording_support::field("reference_steps",
                std::vector<double>{coarse.value_or(nan), mid.value_or(nan), fine.value_or(nan)}),
                recording_support::field("rate", rate), recording_support::field("expected_rate", expected_rate),
                recording_support::field("noise_bounded", noise_bounded ? 1U : 0U),
                recording_support::field("absolute_floor", !structural && std::abs(extrapolated.value_or(0.0)) <= CRITERION[order - 1] ? 1U : 0U),
                recording_support::field("best_step_relative_agreement", best_step),
                recording_support::field("term_scale", scales[i]),
                recording_support::field("allowance", criterion),
                recording_support::field("reference_uncertainty", reference_uncertainty),
                recording_support::field("signal_over_allowance", std::abs(extrapolated.value_or(0.0)) / criterion)});
            recording_support::compare_custom(id.c_str(), "action", item.unit.c_str(), action, extrapolated.value_or(nan),
                ladder_inputs, "Relationship", affine ? "affine-constraint-identity" : "centered-solve-ladder",
                affine ? "the defining affine relation has zero higher action" : "solved-state action agrees with the h^2-extrapolated fresh-solve reference",
                "|action-reference|<=allowance && reference_uncertainty<=allowance && resolved_signal",
                verdict, group, id.c_str(), "D^k g", "case,observable,order", "|D^k g-R|\\leq c", "observable,order");
            // A ladder certified by step agreement records that agreement, not an uncomputed rate.
            const bool by_rate = std::isfinite(rate);
            const double spread = complete ? std::abs(noise_bounded ? *coarse - *mid : *fine - *mid) : nan;
            recording_support::compare_custom((id + "/rate").c_str(), by_rate ? "difference_ratio" : "step_agreement",
                by_rate ? "1" : item.unit.c_str(),
                by_rate ? rate : spread, by_rate ? expected_rate : 2.0 * criterion, ladder_inputs,
                "Relationship", "centered-solve-ladder",
                "the reference ladder certifies second-order convergence, or its last two informative steps agree within twice the criterion",
                "|fine-mid| <= 2c || expected/2 <= rate <= 2 expected || (noise_bounded && |coarse-mid| <= 2c)", refined,
                group, id.c_str(), by_rate ? "(\\Delta_1-\\Delta_2)/(\\Delta_2-\\Delta_3)" : "|R_i-R_{i+1}|", "case,observable,order",
                by_rate ? "\\approx(h_1^2-h_2^2)/(h_2^2-h_3^2)" : "\\leq 2c", "observable,order");
        }
    }
    return true;  // every verdict above throws when false; the surface-equilibrium replay reads this
}

bool record_bool(const char* id, const char* group, const char* claim, bool verdict, const std::string& inputs) {
    recording_support::compare(id, "verdict", "bool", verdict ? 1.0 : 0.0, 1.0, inputs, "Relationship", "solved-state-action",
        claim, "absolute", 0.0, 0.0, verdict, group, id, "verdict", "inputs", claim, "case");
    return verdict;
}

SolvedStateActionBatch actions(const Fixture& fixture, const epcsaft::HelmholtzModel& model, const EquilibriumResult& central,
                               std::vector<SolvedStateObservable> observables, std::size_t order) {
    SolvedStateActionRequest request;
    request.observables = std::move(observables);
    request.directions.assign(fixture.directions.begin(), fixture.directions.begin() + static_cast<std::ptrdiff_t>(order));
    return epcsaft::equilibrium::solved_state_actions(epcsaft::equilibrium::compile_problem(model, fixture.problem), central, request);
}

Fixture neutral_case() {
    Fixture f;
    f.label = "neutral-ternary";
    f.parameters = properties_support::neutral_ternary();
    f.active = {{epcsaft::ParameterFamily::Kij, {1, 3}}, {epcsaft::ParameterFamily::SegmentCount, {0}}};
    f.problem = equilibrium_support::flash(3);
    f.problem.phases = {equilibrium_support::phase("V", PhaseKind::Vapor), equilibrium_support::phase("L", PhaseKind::Liquid)};
    f.problem.temperature.value = 180.0;
    f.directions = {direction(0.2, 2000.0, {-0.01, 0.006, 0.004}, {0.001, 0.01}),
                    direction(-0.1, -1000.0, {0.003, -0.004, 0.001}, {-0.0005, 0.004}),
                    direction(0.05, 500.0, {0.002, 0.001, -0.003}, {0.0002, -0.002})};
    f.observables = geometry_observables(f.problem, 3);
    return f;
}

Fixture pinned_case(double fraction) {
    Fixture f;
    f.label = fraction > 0.0 ? "positive-pin" : "zero-pin";
    f.parameters = properties_support::neutral_binary();
    f.problem = equilibrium_support::boundary();
    f.problem.phases[1].amount = AmountSpec::pinned(fraction);
    f.problem.temperature = FreeValue::fixed(180.0);
    f.problem.feed.values = {0.5, 0.5};
    f.directions = {direction(0.2, 0.0, {0.01, 0.02}), direction(-0.1, 0.0, {0.005, -0.003}), direction(0.05, 0.0, {-0.002, 0.006})};
    f.observables = geometry_observables(f.problem, 2);
    f.observables.push_back(obs(Kind::PhaseProperty, 1, 0, PropertyObservable::ResidualEnthalpy));
    f.observables.push_back(obs(Kind::PhaseProperty, 1, 0, PropertyObservable::ResidualIsobaricHeatCapacity));
    return f;
}

}  // namespace e5

void check_e5_neutral() {
    using namespace e5;
    const auto fixture = neutral_case();
    qualify(fixture, "e5-1");
    // Permutations: values and actions follow physical identity, not vector position.
    const auto model = fixture.model();
    const auto central = epcsaft::equilibrium::solve_equilibrium(model, fixture.problem);
    const auto compiled = epcsaft::equilibrium::compile_problem(model, fixture.problem);
    auto independent_directions = fixture.directions;
    independent_directions.push_back(direction(0.0, 0.0, {0.0, 0.0, 0.0}, {0.0, 0.0}));
    const auto independent = epcsaft::equilibrium::solved_state_first_actions(
        compiled, central, fixture.observables, independent_directions);
    bool same_a1 = independent.size() == independent_directions.size();
    double worst_a1 = 0.0;
    std::size_t factorizations = 0;
    for (std::size_t j = 0; j < independent.size() && same_a1; ++j) {
        const auto separate = epcsaft::equilibrium::solved_state_actions(
            compiled, central, {fixture.observables, {independent_directions[j]}});
        factorizations += independent[j].diagnostics.state_factorizations;
        same_a1 = independent[j].results.size() == separate.results.size();
        for (std::size_t i = 0; i < separate.results.size() && same_a1; ++i) {
            const auto& actual = independent[j].results[i];
            const auto& expected = separate.results[i];
            same_a1 = actual.order == 1 && actual.unit == expected.unit && actual.status == expected.status
                && actual.diagnostic_message == expected.diagnostic_message
                && actual.value.has_value() == expected.value.has_value()
                && actual.action.has_value() == expected.action.has_value();
            if (actual.value && expected.value)
                worst_a1 = std::max(worst_a1, std::abs(*actual.value - *expected.value) / std::max(1.0, std::abs(*expected.value)));
            if (actual.action && expected.action)
                worst_a1 = std::max(worst_a1, std::abs(*actual.action - *expected.action) / std::max(1.0, std::abs(*expected.action)));
        }
    }
    record_bool("e5-1/independent-A1", "e5-1", "four independent A1 responses, including zero, equal separate actions with one center factorization",
        same_a1 && worst_a1 <= 2.0e-12 && factorizations == 1,
        recording_support::fields({recording_support::field("maximum_scaled_change", worst_a1),
            recording_support::field("shared_factorizations", factorizations),
            recording_support::field("directions", independent_directions.size())}));
    const auto reference_batch = actions(fixture, model, central, fixture.observables, 2);
    Fixture swapped = fixture;
    swapped.label = "phase-swap";
    std::swap(swapped.problem.phases[0], swapped.problem.phases[1]);
    std::swap(swapped.active[0], swapped.active[1]);
    for (auto& d : swapped.directions) std::swap(d.active_parameters[0], d.active_parameters[1]);
    for (auto& item : swapped.observables) if (item.kind != Kind::Temperature && item.kind != Kind::Pressure) item.phase = 1 - item.phase;
    const auto swapped_model = swapped.model();
    const auto swapped_batch = actions(swapped, swapped_model, epcsaft::equilibrium::solve_equilibrium(swapped_model, swapped.problem), swapped.observables, 2);
    Fixture permuted = fixture;
    permuted.label = "component-permutation";
    const std::vector<std::size_t> order = {2, 0, 1};
    permuted.parameters = equilibrium_support::permute_parameters(fixture.parameters, order);
    permuted.problem.feed.values = equilibrium_support::permute_vector(fixture.problem.feed.values, order);
    permuted.active = {{epcsaft::ParameterFamily::Kij, {5, 7}}, {epcsaft::ParameterFamily::SegmentCount, {1}}};
    for (auto& d : permuted.directions) d.feed_amounts_mol = equilibrium_support::permute_vector(d.feed_amounts_mol, order);
    for (auto& item : permuted.observables)
        if (item.kind == Kind::PhaseComponentAmount || item.kind == Kind::PhaseMoleFraction
            || item.kind == Kind::PhaseResidualChemicalPotentialOverRt || item.kind == Kind::PhaseLogActivity)
            item.component = static_cast<std::size_t>(std::find(order.begin(), order.end(), item.component) - order.begin());
    const auto permuted_model = permuted.model();
    const auto permuted_batch = actions(permuted, permuted_model, epcsaft::equilibrium::solve_equilibrium(permuted_model, permuted.problem), permuted.observables, 2);
    for (const auto* pair : {&swapped_batch, &permuted_batch}) {
        double worst = 0.0;
        bool statuses = true;
        for (std::size_t i = 0; i < reference_batch.results.size(); ++i) {
            const auto& a = reference_batch.results[i];
            const auto& b = pair->results[i];
            statuses = statuses && a.status == b.status && a.action.has_value() && b.action.has_value() && a.value.has_value() && b.value.has_value();
            if (statuses) worst = std::max(worst, std::max(std::abs(*a.action - *b.action) / std::max(1.0, std::abs(*a.action)),
                                                              std::abs(*a.value - *b.value) / std::max(1.0, std::abs(*a.value))));
        }
        const std::string label = pair == &swapped_batch ? "phase-swap" : "component-permutation";
        recording_support::compare(("e5-9/" + label).c_str(), "identity_difference", "1", worst, 0.0,
            recording_support::fields({recording_support::field("case", label), recording_support::field("results", reference_batch.results.size())}),
            "Relationship", "permutation", "A2 values and actions follow phase/component/active identity under permutation",
            "reference_or_one", 1.0e-8, 0.0, statuses && worst <= 1.0e-8, "e5-9", label.c_str(), "\\Delta", "identity", "\\Delta\\leq10^{-8}", "case");
    }
    // Warm reuse: a repeated identical request on the recorded tape reproduces the batch exactly.
    const auto again = actions(fixture, model, central, fixture.observables, 2);
    bool reproduced = again.results.size() == reference_batch.results.size();
    for (std::size_t i = 0; reproduced && i < again.results.size(); ++i)
        reproduced = again.results[i].status == reference_batch.results[i].status
            && again.results[i].value == reference_batch.results[i].value
            && again.results[i].action == reference_batch.results[i].action;
    record_bool("e5-9/warm-reuse", "e5-9", "a repeated request reproduces every status, value and action of the batch exactly",
        reproduced, recording_support::fields({recording_support::field("results", again.results.size())}));
}

void check_e5_pins() {
    using namespace e5;
    const auto positive = pinned_case(0.01);
    qualify(positive, "e5-2");
    const auto model = positive.model();
    const auto central = epcsaft::equilibrium::solve_equilibrium(model, positive.problem);
    const std::vector<SolvedStateObservable> requested = {obs(Kind::PhaseAmount, 1), obs(Kind::PhaseVolume, 1), obs(Kind::PhaseMolarDensity, 1)};
    std::vector<SolvedStateActionBatch> batches;
    for (std::size_t order = 1; order <= 3; ++order) batches.push_back(actions(positive, model, central, requested, order));
    // Pinned amount: N_V = 0.01 N_feed, so A1 = 0.01 sum(dfeed) and A2/A3 vanish exactly.
    const double expected = 0.01 * (positive.directions[0].feed_amounts_mol[0] + positive.directions[0].feed_amounts_mol[1]);
    const bool amount_ok = std::abs(batches[0].results[0].action.value_or(1.0) - expected) <= STRUCTURAL_ZERO
        && std::abs(batches[1].results[0].action.value_or(1.0)) <= STRUCTURAL_ZERO
        && std::abs(batches[2].results[0].action.value_or(1.0)) <= STRUCTURAL_ZERO;
    recording_support::compare("e5-2/positive-pin/pinned-amount-A1", "action", "mol", batches[0].results[0].action.value_or(std::numeric_limits<double>::quiet_NaN()),
        expected, recording_support::fields({recording_support::field("A2", epcsaft::core::optional_number(batches[1].results[0].action)),
            recording_support::field("A3", epcsaft::core::optional_number(batches[2].results[0].action))}),
        "Relationship", "pinned-amount", "pinned amount A1 is the fraction times the feed direction sum; A2 and A3 vanish",
        "strict_absolute", STRUCTURAL_ZERO, 0.0, amount_ok, "e5-2", "pinned-amount", "D N_V", "dfeed", "D N_V=0.01\\sum dn", "case");
    // Extensive volume follows the product rule V = N / rho from the returned N and rho actions.
    const auto& a1 = batches[0].results;
    const auto& a2 = batches[1].results;
    const double N = *a1[0].value, rho = *a1[2].value, Na = *a1[0].action, ra = *a1[2].action;
    const auto b1 = actions(positive, model, central, requested, 1);
    (void)b1;
    const double v1 = Na / rho - N * ra / (rho * rho);
    SolvedStateActionRequest second;
    second.observables = requested;
    second.directions = {positive.directions[1]};
    const auto other = epcsaft::equilibrium::solved_state_actions(epcsaft::equilibrium::compile_problem(model, positive.problem), central, second);
    const double Nb = *other.results[0].action, rb = *other.results[2].action, Nab = *a2[0].action, rab = *a2[2].action;
    const double v2 = Nab / rho - (Na * rb + Nb * ra) / (rho * rho) - N * rab / (rho * rho) + 2.0 * N * ra * rb / (rho * rho * rho);
    const bool volume_ok = std::abs(*a1[1].action - v1) <= 1.0e-10 * std::max(1.0, std::abs(v1))
        && std::abs(*a2[1].action - v2) <= 1.0e-10 * std::max(1.0, std::abs(v2));
    recording_support::compare("e5-2/positive-pin/volume-product-rule", "action", "m^3", *a2[1].action, v2,
        recording_support::fields({recording_support::field("A1_returned", *a1[1].action), recording_support::field("A1_product_rule", v1)}),
        "Relationship", "product-rule", "extensive volume actions equal the product rule of the returned amount and density actions",
        "reference_or_one", 0.0, 1.0e-10, volume_ok, "e5-2", "volume-product-rule", "D^2 V", "N,\\rho", "V=N/\\rho", "case");
    // Exactly zero pin: extensive rows are structural zero, intensive rows follow the ladders.
    auto zero = pinned_case(0.0);
    std::vector<SolvedStateObservable> intensive;
    for (const auto& item : zero.observables) {
        const bool extensive = item.phase == 1 && (item.kind == Kind::PhaseAmount || item.kind == Kind::PhaseComponentAmount || item.kind == Kind::PhaseVolume);
        if (!extensive) intensive.push_back(item);
    }
    const auto zero_model = zero.model();
    const auto zero_central = epcsaft::equilibrium::solve_equilibrium(zero_model, zero.problem);
    bool structural = zero_central.success;
    for (std::size_t order = 1; order <= 3 && structural; ++order) {
        const auto batch = actions(zero, zero_model, zero_central, {obs(Kind::PhaseAmount, 1), obs(Kind::PhaseComponentAmount, 1, 0), obs(Kind::PhaseVolume, 1)}, order);
        for (const auto& item : batch.results)
            structural = structural && item.status == SolvedStateActionStatus::Available && std::abs(item.value.value_or(1.0)) <= STRUCTURAL_ZERO
                && std::abs(item.action.value_or(1.0)) <= STRUCTURAL_ZERO;
    }
    record_bool("e5-2/zero-pin/structural-zero", "e5-2", "a zero pinned amount, its component amounts and its extensive volume are structural zeros through A3",
        structural, recording_support::fields({recording_support::field("success", zero_central.success ? 1U : 0U)}));
    zero.observables = intensive;
    qualify(zero, "e5-2");
}

namespace e5 {

Fixture reactive_case() {
    Fixture f;
    f.label = "reactive-two-liquid";
    f.parameters = equilibrium_support::synthetic_reactive_ternary();
    f.problem = equilibrium_support::synthetic_reactive_problem();
    f.directions = {direction(0.2, 1000.0, {0.001, -0.001, 0.0}), direction(-0.1, -500.0, {0.0, 0.001, -0.001}),
                    direction(0.05, 250.0, {-0.001, 0.0, 0.001})};
    f.observables = geometry_observables(f.problem, 3);
    return f;
}

Fixture isomerisation_case(double total, double a) {
    Fixture f;
    f.label = total == 1.0 ? "isomerisation" : "isomerisation-two-mol";
    f.parameters = e2::identical_ethane_pair();
    f.problem = equilibrium_support::flash(2);
    f.problem.phases = {equilibrium_support::phase("V", PhaseKind::Vapor)};
    f.problem.feed.values = {0.5 * total, 0.5 * total};
    f.problem.temperature = FreeValue::fixed(350.0);
    f.problem.pressure = FreeValue::fixed(2.0e5);
    f.problem.reactions = {{"A=B", {-1.0, 1.0}, std::nullopt, std::nullopt, legacy_isomerisation(1.0, a)}};
    f.directions = {direction(1.0, 0.0, {0.0, 0.0}), direction(1.0, 0.0, {0.0, 0.0}), direction(1.0, 0.0, {0.0, 0.0})};
    f.observables = {obs(Kind::PhaseMoleFraction, 0, 1), obs(Kind::PhaseComponentAmount, 0, 1)};
    f.scaled = false;
    return f;
}

}  // namespace e5

void check_e5_reactive() {
    using namespace e5;
    const auto fixture = reactive_case();
    qualify(fixture, "e5-3");
    {
        // At least one nonzero mixed action and the retained central affinities of both phases.
        const auto model = fixture.model();
        const auto central = epcsaft::equilibrium::solve_equilibrium(model, fixture.problem);
        const auto batch = actions(fixture, model, central, fixture.observables, 2);
        double largest = 0.0;
        for (const auto& item : batch.results) if (item.action.has_value()) largest = std::max(largest, std::abs(*item.action));
        std::vector<double> affinity(2, 0.0);
        for (std::size_t row = 0; row < central.rows.size(); ++row)
            if (central.rows[row].kind == epcsaft::equilibrium::RowKind::Chemical)
                affinity[central.rows[row].phase] += fixture.problem.reactions[0].stoichiometry[central.rows[row].index] * central.residuals[row];
        record_bool("e5-3/reactive-two-liquid/mixed-nonzero", "e5-3", "the reactive two-liquid state carries a nonzero A2 mixed action with zero reaction affinity in both phases",
            largest > 1.0e-6 && std::abs(affinity[0]) <= 1.0e-9 && std::abs(affinity[1]) <= 1.0e-9,
            recording_support::fields({recording_support::field("largest_mixed_action", largest), recording_support::field("affinity", affinity)}));
    }
    // Exact nonconstant-correlation isomerisation: p = K/(1+K) with analytic thermal derivatives.
    const double a = moderate_isomerisation_a(), T = 350.0;
    const double L1 = -1200.0 / (T * T) + 0.4 / T + 0.002, L2 = 2.0 * 1200.0 / (T * T * T) - 0.4 / (T * T), L3 = -6.0 * 1200.0 / (T * T * T * T) + 2.0 * 0.4 / (T * T * T);
    const double ln_k = a + 1200.0 / T + 0.4 * std::log(T) + 0.002 * T, p = std::exp(ln_k) / (1.0 + std::exp(ln_k));
    const auto logistic = [&](double p) { return std::array<double, 3>{p * (1.0 - p) * L1, p * (1.0 - p) * (L2 + (1.0 - 2.0 * p) * L1 * L1),
        p * (1.0 - p) * (L3 + 3.0 * (1.0 - 2.0 * p) * L1 * L2 + (1.0 - 6.0 * p + 6.0 * p * p) * L1 * L1 * L1)}; };
    const auto analytic = logistic(p);
    const double expected[3] = {-1.4784580498866213e-3, 8.4348805281749936e-6, -2.7304974653104138e-9};
    auto added_terms = isomerisation_case(1.0, a);
    added_terms.problem.reactions[0].correlation->e = 90000.0;
    added_terms.problem.reactions[0].correlation->f = 1e-6;
    qualify(added_terms, "e5-3-added-terms");
    for (const double total : {1.0, 2.0}) {
        auto fixture_t = isomerisation_case(total, a);
        const auto model = fixture_t.model();
        const auto central = epcsaft::equilibrium::solve_equilibrium(model, fixture_t.problem);
        CHECK(central.success);
        for (std::size_t order = 1; order <= 3; ++order) {
            for (const bool rescaled : {false, true}) {
                auto directions = fixture_t;
                if (rescaled) {
                    directions.directions = order == 2 ? std::vector<SolvedStateActionDirection>{direction(1.0e-6, 0.0, {0.0, 0.0}), direction(1.0e6, 0.0, {0.0, 0.0})}
                        : std::vector<SolvedStateActionDirection>{direction(1.0e-6, 0.0, {0.0, 0.0}), direction(1.0, 0.0, {0.0, 0.0}), direction(1.0e6, 0.0, {0.0, 0.0})};
                    if (order == 1) continue;
                }
                const auto batch = actions(directions, model, central, fixture_t.observables, order);
                const double bar = 2.0e-10 * std::abs(analytic[order - 1]) + 1.0e-14;
                const bool fraction_ok = batch.results[0].action.has_value() && std::abs(*batch.results[0].action - analytic[order - 1]) <= bar
                    && std::abs(analytic[order - 1] - expected[order - 1]) <= bar;
                const bool amount_ok = batch.results[1].action.has_value() && std::abs(*batch.results[1].action - total * analytic[order - 1]) <= total * bar;
                const auto id = std::string("e5-3/") + fixture_t.label + (rescaled ? "/rescaled" : "") + "/A" + std::to_string(order);
                recording_support::compare(id.c_str(), "action", "1/K^k", batch.results[0].action.value_or(std::numeric_limits<double>::quiet_NaN()), analytic[order - 1],
                    recording_support::fields({recording_support::field("total", total), recording_support::field("order", order), recording_support::field("rescaled", rescaled ? 1U : 0U),
                        recording_support::field("amount_action", epcsaft::core::optional_number(batch.results[1].action)), recording_support::field("bar", bar),
                        recording_support::field("status", status_field(batch, 0))}),
                    "Relationship", "analytic-correlation", "x_B thermal actions equal the analytic logistic derivatives of ln K(T); the B amount scales with the total",
                    "strict_absolute", bar, 0.0, fraction_ok && amount_ok, "e5-3", id.c_str(), "d^k x_B/dT^k", "T,K(T)", "p^{(k)}=p(1-p)[\\ldots]", "case");
            }
        }
        if (total == 2.0) {
            // Mixed amount/temperature chain: n_B = N p(T) is degree one in the feed, x_B degree zero.
            auto mixed = fixture_t;
            mixed.directions = {direction(0.0, 0.0, {0.05, 0.05}), direction(1.0, 0.0, {0.0, 0.0})};
            const auto batch = actions(mixed, model, central, fixture_t.observables, 2);
            const bool ok = batch.results[0].action.has_value() && std::abs(*batch.results[0].action) <= 1.0e-12
                && batch.results[1].action.has_value() && std::abs(*batch.results[1].action - 0.1 * analytic[0]) <= 2.0e-10 * std::abs(analytic[0]) + 1.0e-14;
            recording_support::compare("e5-3/isomerisation-two-mol/amount-temperature", "action", "mol/K", batch.results[1].action.value_or(std::numeric_limits<double>::quiet_NaN()), 0.1 * analytic[0],
                recording_support::fields({recording_support::field("fraction_mixed", epcsaft::core::optional_number(batch.results[0].action))}),
                "Relationship", "analytic-correlation", "the amount/temperature mixed action is sum(dfeed) p'(T) while the composition mixed action vanishes",
                "strict_absolute", 2.0e-10 * std::abs(analytic[0]) + 1.0e-14, 0.0, ok, "e5-3", "amount-temperature", "D^2 n_B", "dn,dT", "D^2 n_B=\\sum dn\\,p'", "case");
        }
    }
    {
        // Zero directions give exact zero actions after admission.
        auto zero = isomerisation_case(1.0, a);
        zero.directions = {direction(0.0, 0.0, {0.0, 0.0}), direction(1.0, 0.0, {0.0, 0.0}), direction(0.0, 0.0, {0.0, 0.0})};
        const auto model = zero.model();
        const auto central = epcsaft::equilibrium::solve_equilibrium(model, zero.problem);
        const auto batch = actions(zero, model, central, zero.observables, 3);
        record_bool("e5-3/isomerisation/zero-direction", "e5-3", "a zero direction yields an available exact zero action with the central value retained",
            batch.results[0].status == SolvedStateActionStatus::Available && batch.results[0].action == 0.0 && batch.results[0].value.has_value(),
            recording_support::fields({recording_support::field("status", status_field(batch, 0)), recording_support::field("value", epcsaft::core::optional_number(batch.results[0].value))}));
    }
    {
        // Fixed-T scalar K = 2: a conserved constant Gibbs shift leaves actions unchanged, but a
        // point-only supplied vector cannot supply the correlated thermal reference.
        auto scalar = isomerisation_case(1.0, a);
        scalar.label = "scalar-k";
        scalar.problem.reactions = {{"A=B", {-1.0, 1.0}, 2.0, std::nullopt}};
        scalar.directions = {direction(0.0, 1.0e4, {0.01, -0.01}), direction(0.0, -5.0e3, {0.0, 0.005}), direction(0.0, 2.0e3, {0.003, 0.0})};
        scalar.observables = geometry_observables(scalar.problem, 2);
        auto shifted = scalar;
        shifted.label = "scalar-k-shifted";
        // nu.g = -RT ln K with a conserved constant offset of 1500 J/mol on both species.
        shifted.problem.standard_gibbs_j_per_mol = {1500.0 + epcsaft::GAS_CONSTANT_J_PER_MOL_K * 350.0 * std::log(2.0) / 2.0,
                                                    1500.0 - epcsaft::GAS_CONSTANT_J_PER_MOL_K * 350.0 * std::log(2.0) / 2.0};
        const auto model = scalar.model();
        const auto central = epcsaft::equilibrium::solve_equilibrium(model, scalar.problem);
        const auto shifted_central = epcsaft::equilibrium::solve_equilibrium(model, shifted.problem);
        const auto batch = actions(scalar, model, central, scalar.observables, 2);
        const auto shifted_batch = actions(shifted, model, shifted_central, shifted.observables, 2);
        double worst = 0.0;
        bool available = true;
        for (std::size_t i = 0; i < batch.results.size(); ++i) {
            const auto& x = batch.results[i]; const auto& y = shifted_batch.results[i];
            available = available && x.action.has_value() && y.action.has_value();
            if (available) worst = std::max(worst, std::abs(*x.action - *y.action) / std::max(1.0, std::abs(*x.action)));
        }
        recording_support::compare("e5-3/scalar-k/gauge", "action_difference", "1", worst, 0.0,
            recording_support::fields({recording_support::field("results", batch.results.size())}),
            "Relationship", "gauge", "a conserved constant standard-state shift leaves fixed-T A2 actions unchanged",
            "reference_or_one", 1.0e-8, 0.0, available && worst <= 1.0e-8, "e5-3", "gauge", "\\Delta", "g", "\\Delta\\leq10^{-8}", "case");
        auto thermal = shifted;
        thermal.directions = {direction(1.0, 0.0, {0.0, 0.0})};
        const auto thermal_batch = actions(thermal, model, shifted_central, thermal.observables, 1);
        bool typed = true;
        for (const auto& item : thermal_batch.results) {
            const bool fixed_input = item.observable.kind == Kind::Temperature || item.observable.kind == Kind::Pressure;
            typed = typed && item.value.has_value()
                && (fixed_input ? item.status == SolvedStateActionStatus::Available && item.action.has_value()
                                : item.status == SolvedStateActionStatus::MissingThermochemistry && !item.action.has_value());
        }
        record_bool("e5-3/scalar-k/point-only-thermal", "e5-3", "a point-only supplied Gibbs vector returns missing_thermochemistry for thermal state actions and keeps fixed-input identities",
            typed, recording_support::fields({recording_support::field("status", status_field(thermal_batch, 2))}));
        const auto compiled = epcsaft::equilibrium::compile_problem(model, shifted.problem);
        const auto nonthermal = epcsaft::equilibrium::solved_state_actions(
            compiled, shifted_central, {shifted.observables, {scalar.directions[0]}});
        const auto independent = epcsaft::equilibrium::solved_state_first_actions(
            compiled, shifted_central, shifted.observables, {thermal.directions[0], scalar.directions[0]});
        bool isolated = independent.size() == 2;
        double worst_a1 = 0.0;
        for (std::size_t j = 0; j < independent.size() && isolated; ++j) {
            const auto& expected = j == 0 ? thermal_batch : nonthermal;
            isolated = independent[j].results.size() == expected.results.size();
            for (std::size_t i = 0; i < expected.results.size() && isolated; ++i) {
                const auto& actual = independent[j].results[i];
                const auto& single = expected.results[i];
                isolated = actual.status == single.status && actual.value.has_value() == single.value.has_value()
                    && actual.action.has_value() == single.action.has_value();
                if (actual.action && single.action)
                    worst_a1 = std::max(worst_a1, std::abs(*actual.action - *single.action) / std::max(1.0, std::abs(*single.action)));
            }
        }
        record_bool("e5-3/scalar-k/independent-thermal", "e5-3", "a refused thermal A1 does not suppress a following available feed A1 at the same center",
            isolated && worst_a1 <= 2.0e-12 && independent[0].results[2].status == SolvedStateActionStatus::MissingThermochemistry
                && independent[1].results[2].status == SolvedStateActionStatus::Available,
            recording_support::fields({recording_support::field("maximum_scaled_change", worst_a1)}));
    }
    {
        // A -> 2B with matching species thermochemistry: thermal A1-A3 of composition, amounts and calorics.
        const double R = epcsaft::GAS_CONSTANT_J_PER_MOL_K, tref = 300.0;
        const epcsaft::ReactionLogPolynomial density_law{-10.0 / R - std::log(R * tref / 1.0e5), 3000.0 / R, 10.0 / R - 1.0, 0.0,
            tref, 250.0, 450.0, epcsaft::equilibrium::EOS_STANDARD_STATE_ID};
        Fixture f;
        f.label = "a-to-2b";
        f.parameters = e2::identical_ethane_pair();
        f.record = e2::species_record(30.0, 20.0);
        f.problem = equilibrium_support::flash(2);
        f.problem.phases = {equilibrium_support::phase("V", PhaseKind::Vapor)};
        f.problem.feed.values = {0.5, 0.5};
        f.problem.temperature = FreeValue::fixed(350.0);
        f.problem.pressure = FreeValue::fixed(2.0e5);
        f.problem.reactions = {{"A=2B", {-1.0, 2.0}, std::nullopt, std::nullopt, density_law}};
        f.directions = {direction(20.0, 0.0, {0.0, 0.0}), direction(-10.0, 0.0, {0.0, 0.0}), direction(5.0, 0.0, {0.0, 0.0})};
        f.observables = {obs(Kind::PhaseMoleFraction, 0, 1), obs(Kind::PhaseComponentAmount, 0, 0), obs(Kind::PhaseComponentAmount, 0, 1), obs(Kind::PhaseAmount, 0),
            obs(Kind::PhaseProperty, 0, 0, PropertyObservable::TotalEnthalpy), obs(Kind::PhaseProperty, 0, 0, PropertyObservable::TotalIsobaricHeatCapacity)};
        f.scaled = false;
        // E2's independent pure-state relation y^2/(1-y) = K_P P0 / (phi P) at each corner
        // temperature supplies the composition and amount references; calorics use fresh solves.
        const epcsaft::ReactionLogPolynomial pressure_law{-10.0 / R, 3000.0 / R, 10.0 / R, 0.0, tref, 250.0, 450.0, "pressure-basis-P0"};
        const auto ethane = std::make_shared<epcsaft::HelmholtzModel>(properties_support::ethane());
        f.analytic = [ethane, pressure_law, R](const SolvedStateObservable& item, double T) -> std::optional<double> {
            if (item.kind != Kind::PhaseMoleFraction && item.kind != Kind::PhaseComponentAmount && item.kind != Kind::PhaseAmount) return std::nullopt;
            const auto pure = ethane->state_at_pressure(T, 2.0e5, {1.0}, epcsaft::DensityBranch::Vapor);
            const double compressibility = 2.0e5 / (pure.molar_density * R * T);
            const double phi = std::exp(pure.residual_chemical_potential_over_rt[0].value()) / compressibility;
            const double a = std::exp(epcsaft::evaluate_reaction_properties(pressure_law, T).ln_k) * 1.0e5 / (phi * 2.0e5);
            const double y = 2.0 * a / (a + std::sqrt(a * a + 4.0 * a));
            const double total = 1.5 / (2.0 - y);
            if (item.kind == Kind::PhaseMoleFraction) return item.component == 1 ? y : 1.0 - y;
            if (item.kind == Kind::PhaseAmount) return total;
            return item.component == 1 ? total * y : total * (1.0 - y);
        };
        qualify(f, "e5-3");
    }
    {
        // #322: x_B = K/(1+K) below EQUILIBRIUM_BOX_FLOOR (1e-20). Main refused the actions and the
        // declared seeds there, stopped continuation at x_B = 1.0e-20 and left the square bubble cell infeasible.
        const auto model = isomerisation_case(1.0, a).model();
        const auto isomer = [&](double k, double t, double b, double c, double d) {  // K(350 K) = k
            auto problem = isomerisation_case(1.0, a).problem;
            problem.temperature = FreeValue::fixed(t);
            problem.solver_options.values[epcsaft::policy::row] = 1.0e-10;
            problem.reactions[0].correlation = epcsaft::ReactionLogPolynomial{std::log(k) - b / T - c * std::log(T / 300.0) - d * T,
                b, c, d, 300.0, 250.0, 450.0, epcsaft::equilibrium::EOS_STANDARD_STATE_ID};
            return problem;
        };
        const auto row = [&](const std::string& id, double model_value, double reference, const char* claim, const std::string& inputs) {
            recording_support::compare(id.c_str(), "value", "1", model_value, reference, inputs, "Relationship",
                "analytic-correlation", claim, "reference_scale", 0.0, 1.0e-9, std::abs(model_value - reference) <= 1.0e-9 * std::abs(reference), "e5-3", id.c_str());
        };
        const double k = 1.0e-24;
        auto trace = isomerisation_case(1.0, a);
        trace.problem = isomer(k, T, 1200.0, 0.4, 0.002);
        const auto central = epcsaft::equilibrium::solve_equilibrium(model, trace.problem);
        const auto reference = logistic(k / (1.0 + k));
        for (std::size_t order = 1; order <= 3; ++order) {
            const auto batch = actions(trace, model, central, trace.observables, order);
            row("e5-3/sub-cutoff/A" + std::to_string(order), batch.results[0].action.value_or(std::numeric_limits<double>::quiet_NaN()),
                reference[order - 1], "x_B thermal actions below the cutoff equal the analytic logistic derivatives",
                recording_support::fields({recording_support::field("status", status_field(batch, 0))}));
        }
        for (const char* form : {"previous", "coordinates", "composition"}) {
            epcsaft::equilibrium::EquilibriumStart start;
            if (form == std::string("previous")) start.previous = central;
            else if (form == std::string("coordinates")) start.coordinates = central.coordinates;
            else start.compositions = {{1.0 - k, k}};
            const auto seeded = epcsaft::equilibrium::solve_equilibrium_starts(model, trace.problem, {start}, {});
            row(std::string("e5-3/sub-cutoff/seed/") + form, seeded.selected ? seeded.results[0]->mole_fractions[1] : std::numeric_limits<double>::quiet_NaN(),
                k / (1.0 + k), "a declared seed below the cutoff converges to the analytic root", recording_support::fields({recording_support::field("cause", seeded.causes[0])}));
        }
        // ln K = ln(1e-16) - 40000 (1/T - 1/450 K): the path 450 -> 350 K crosses the cutoff near 408 K.
        epcsaft::equilibrium::ContinuationRequest path;
        const double k_cross = 1.0e-16 * std::exp(-40000.0 * (1.0 / T - 1.0 / 450.0));
        path.start_problem = isomer(k_cross, 450.0, -40000.0, 0.0, 0.0);
        path.target_problem = isomer(k_cross, T, -40000.0, 0.0, 0.0);
        path.options = {0.25, 1.0e-3, 50, 8};
        const auto traced = epcsaft::equilibrium::continue_equilibrium(model, path);
        double worst = traced.success ? 0.0 : std::numeric_limits<double>::quiet_NaN();
        for (const auto& state : traced.states) {
            const double kt = 1.0e-16 * std::exp(-40000.0 * (1.0 / *state.temperature - 1.0 / 450.0));
            worst = std::max(worst, std::abs(state.mole_fractions[1] / (kt / (1.0 + kt)) - 1.0));
        }
        row("e5-3/sub-cutoff/continuation", 1.0 + worst, 1.0, "every continuation state through the cutoff is the analytic root",
            recording_support::fields({recording_support::field("stop", epcsaft::equilibrium::continuation_stop_name(traced.stop)),
                recording_support::field("accepted", traced.accepted_steps), recording_support::field("states", traced.states.size())}));
        // Square bubble cell at 260 K: the trace limit is linear in K, so K = 1e-24 is 1e-12 times K = 1e-12.
        std::array<double, 2> y{}, slope{};
        for (std::size_t i = 0; i < 2; ++i) {
            auto bubble = equilibrium_support::boundary();
            bubble.temperature = FreeValue::fixed(260.0);
            bubble.pressure = FreeValue::variable();
            bubble.solver_options.values[epcsaft::policy::row] = 1.0e-10;
            bubble.feed.values = {0.5, 0.5};
            bubble.reactions = isomer(i == 0 ? 1.0e-12 : k, T, 1200.0, 0.4, 0.002).reactions;
            const auto solved = epcsaft::equilibrium::solve_equilibrium(model, bubble);
            y[i] = slope[i] = std::numeric_limits<double>::quiet_NaN();
            if (!solved.success) continue;
            const auto batch = epcsaft::equilibrium::solved_state_actions(epcsaft::equilibrium::compile_problem(model, bubble), solved,
                {{obs(Kind::PhaseMoleFraction, 1, 1)}, {direction(1.0, 0.0, {0.0, 0.0})}});
            y[i] = solved.mole_fractions[3];
            slope[i] = batch.results[0].action.value_or(std::numeric_limits<double>::quiet_NaN()) / y[i];
        }
        row("e5-3/sub-cutoff/bubble/y_B", y[1], 1.0e-12 * y[0], "the square bubble cell reaches the trace root below the cutoff", recording_support::fields({}));
        row("e5-3/sub-cutoff/bubble/dlny_B/dT", slope[1], slope[0], "d ln y_B/dT below the cutoff equals its K = 1e-12 value", recording_support::fields({}));
    }
}

namespace e5 {

Fixture associating_case() {
    Fixture f;
    f.label = "associating";
    f.parameters = properties_support::associating_binary();
    f.active = {{epcsaft::ParameterFamily::AssociationEnergyOverK, {1, 4}}, {epcsaft::ParameterFamily::AssociationVolume, {1, 4}}};
    f.record = properties_support::p3_binary_record();
    f.problem = equilibrium_support::flash(2);
    f.problem.feed.values = {0.25, 0.75};
    f.problem.temperature = FreeValue::fixed(313.15);
    f.problem.pressure = FreeValue::fixed(epcsaft::HelmholtzModel(f.parameters).state(313.15, 42441.09427799721, {0.25, 0.75}).pressure);
    f.directions = {direction(0.2, 15000.0, {-0.01, 0.01}, {0.03, 0.001}), direction(-0.1, -10000.0, {0.005, -0.005}, {-0.02, -0.0005}),
                    direction(0.05, 6000.0, {0.002, -0.002}, {0.01, 0.0002})};
    f.observables = {obs(Kind::PhaseProperty, 0, 0, PropertyObservable::ResidualEnthalpy), obs(Kind::PhaseProperty, 0, 0, PropertyObservable::TotalEnthalpy),
                     obs(Kind::PhaseProperty, 0, 0, PropertyObservable::TotalIsobaricHeatCapacity), obs(Kind::PhaseAmount, 0), obs(Kind::PhaseMolarDensity, 0),
                     obs(Kind::PhaseLogActivity, 0, 0), obs(Kind::PhaseLogActivity, 0, 1)};
    f.scaled = false;  // owner decision 2026-09-16: the scaled A3 corners leave the feed domain
    return f;
}

}  // namespace e5

void check_e5_associating() {
    using namespace e5;
    const auto fixture = associating_case();
    qualify(fixture, "e5-4");
    const auto model = fixture.model();
    const auto central = epcsaft::equilibrium::solve_equilibrium(model, fixture.problem);
    // Selected results equal the corresponding full results.
    const auto full = actions(fixture, model, central, fixture.observables, 3);
    double worst = 0.0;
    bool ok = true;
    for (std::size_t i = 0; i < 3; ++i) {
        const auto selected = actions(fixture, model, central, {fixture.observables[i]}, 3);
        ok = ok && selected.results[0].action.has_value() && full.results[i].action.has_value();
        if (ok) worst = std::max(worst, std::abs(*selected.results[0].action - *full.results[i].action) / std::max(1.0, std::abs(*full.results[i].action)));
    }
    recording_support::compare("e5-4/associating/selected-vs-full", "action_difference", "1", worst, 0.0,
        recording_support::fields({recording_support::field("site_rank", full.results[0].property_diagnostics.has_value() ? full.results[0].property_diagnostics->site_rank : 0),
            recording_support::field("energy_order", full.diagnostics.maximum_energy_order)}),
        "Relationship", "selection", "selected A3 property results equal their full-batch values", "reference_or_one", 2.0e-12, 0.0,
        ok && worst <= 2.0e-12, "e5-4", "selected-vs-full", "\\Delta", "selection", "\\Delta\\leq2\\times10^{-12}", "case");
    // Mixed-action symmetry: the direction order changes neither A2 nor A3.
    auto rotated = fixture;
    rotated.directions = {fixture.directions[1], fixture.directions[2], fixture.directions[0]};
    auto swapped_pair = fixture;
    swapped_pair.directions = {fixture.directions[1], fixture.directions[0]};
    const auto natural_two = actions(fixture, model, central, fixture.observables, 2);
    const auto swapped_two = actions(swapped_pair, model, central, fixture.observables, 2);
    const auto rotated_three = actions(rotated, model, central, fixture.observables, 3);
    double asymmetry = 0.0;
    bool symmetric = true;
    for (std::size_t i = 0; i < fixture.observables.size(); ++i) {
        symmetric = symmetric && swapped_two.results[i].action.has_value() && natural_two.results[i].action.has_value()
            && rotated_three.results[i].action.has_value() && full.results[i].action.has_value();
        if (!symmetric) break;
        asymmetry = std::max({asymmetry, std::abs(*swapped_two.results[i].action - *natural_two.results[i].action) / std::max(1.0, std::abs(*natural_two.results[i].action)),
            std::abs(*rotated_three.results[i].action - *full.results[i].action) / std::max(1.0, std::abs(*full.results[i].action))});
    }
    recording_support::compare("e5-4/associating/mixed-symmetry", "action_difference", "1", asymmetry, 0.0,
        recording_support::fields({recording_support::field("results", fixture.observables.size())}),
        "Relationship", "symmetry", "A2 and A3 mixed actions are invariant under direction-order permutation", "reference_or_one", 1.0e-10, 0.0,
        symmetric && asymmetry <= 1.0e-10, "e5-4", "mixed-symmetry", "\\Delta", "order", "\\Delta\\leq10^{-10}", "case");
    // Radial identities: molar h has Dh[r] = 0 and D2h[v, r] = -Dh[v]; extensive H = N h has DH[r] = H, D2H[v, r] = 0.
    auto radial = fixture;
    radial.directions = {direction(0.0, 0.0, {0.02, -0.01}, {0.0, 0.0}), direction(0.0, 0.0, {0.25, 0.75}, {0.0, 0.0})};
    const std::vector<SolvedStateObservable> requested = {obs(Kind::PhaseProperty, 0, 0, PropertyObservable::TotalEnthalpy), obs(Kind::PhaseAmount, 0)};
    const auto dv = actions(radial, model, central, requested, 1);
    auto only_r = radial;
    only_r.directions = {radial.directions[1]};
    const auto dr = actions(only_r, model, central, requested, 1);
    const auto dvr = actions(radial, model, central, requested, 2);
    const double h = *dv.results[0].value, N = *dv.results[1].value;
    const double hv = *dv.results[0].action, hr = *dr.results[0].action, hvr = *dvr.results[0].action;
    const double Nv = *dv.results[1].action, Nr = *dr.results[1].action, Nvr = *dvr.results[1].action;
    const double tolerance = 1.0e-8 * std::max(1.0, std::abs(hv));
    const bool molar_ok = std::abs(hr) <= tolerance && std::abs(hvr + hv) <= tolerance;
    const double H_r = Nr * h + N * hr, H_vr = Nvr * h + Nv * hr + Nr * hv + N * hvr;
    const bool extensive_ok = std::abs(H_r - N * h) <= 1.0e-8 * std::max(1.0, std::abs(N * h)) && std::abs(H_vr) <= 1.0e-8 * std::max(1.0, std::abs(N * hv));
    recording_support::compare("e5-4/associating/radial-molar", "action", "J/mol", hvr, -hv,
        recording_support::fields({recording_support::field("Dh_r", hr), recording_support::field("Dh_v", hv), recording_support::field("H_r", H_r),
            recording_support::field("H", N * h), recording_support::field("H_vr", H_vr)}),
        "Relationship", "homogeneity", "molar enthalpy obeys Dh[r]=0 and D2h[v,r]=-Dh[v]; the extensive product H=N h obeys DH[r]=H and D2H[v,r]=0",
        "reference_or_one", 0.0, 1.0e-8, molar_ok && extensive_ok, "e5-4", "radial", "D^2 h[v,r]", "v,r", "D^2h[v,r]=-Dh[v]", "case");
}

namespace e5 {

Fixture electrolyte_case(bool with_record) {
    Fixture f;
    f.label = with_record ? "electrolyte" : "electrolyte-residual-only";
    f.parameters = properties_support::aqueous_electrolyte();
    if (with_record) {
        epcsaft::ThermochemistryRecord record;
        record.standard_pressure = 1.0e5;
        for (std::size_t i = 0; i < 3; ++i) record.component_intervals.push_back(properties_support::p3_record(properties_support::p3_forms()[0].second).component_intervals.front());
        f.record = record;
    }
    f.problem = equilibrium_support::flash(3);
    f.problem.feed.values = {0.98, 0.01, 0.01};
    f.problem.temperature = FreeValue::fixed(298.15);
    f.directions = {direction(0.1, 1000.0, {-0.002, 0.001, 0.001}), direction(-0.05, -500.0, {0.001, -0.0005, -0.0005}),
                    direction(0.02, 250.0, {-0.0004, 0.0002, 0.0002})};
    f.observables = {obs(Kind::PhaseProperty, 0, 0, PropertyObservable::ResidualEnthalpy), obs(Kind::PhaseProperty, 0, 0, PropertyObservable::ResidualIsobaricHeatCapacity)};
    if (with_record) {
        f.observables.push_back(obs(Kind::PhaseProperty, 0, 0, PropertyObservable::TotalEnthalpy));
        f.observables.push_back(obs(Kind::PhaseProperty, 0, 0, PropertyObservable::TotalIsobaricHeatCapacity));
    }
    f.scaled = false;
    return f;
}

Fixture two_root_case(PhaseKind kind) {
    Fixture f;
    f.label = kind == PhaseKind::Vapor ? "two-root-vapor" : "two-root-liquid";
    f.parameters = properties_support::ethane();
    f.problem = equilibrium_support::flash(1);
    f.problem.phases = {equilibrium_support::phase("P", kind)};
    f.problem.temperature = FreeValue::fixed(200.0);
    f.problem.pressure = FreeValue::fixed(217800.69494692178);
    f.directions = {direction(0.1, 2000.0, {0.0}), direction(-0.05, -1000.0, {0.0}), direction(0.03, 500.0, {0.0})};
    f.observables = {obs(Kind::PhasePressure, 0), obs(Kind::PhaseProperty, 0, 0, PropertyObservable::ResidualEnthalpy),
                     obs(Kind::PhaseProperty, 0, 0, PropertyObservable::ResidualIsobaricHeatCapacity), obs(Kind::PhaseMolarDensity, 0)};
    f.corner_branch = kind == PhaseKind::Vapor ? epcsaft::DensityBranch::Vapor : epcsaft::DensityBranch::Liquid;
    return f;
}

Fixture thermal_case(double temperature, bool with_record) {
    Fixture f;
    f.label = "thermal-" + recording_support::number(temperature) + (with_record ? "" : "-no-record");
    f.constant_ideal_cp = true; // p3_forms()[1] is Shomate cp=30 J/(mol K)
    f.parameters = properties_support::ethane();
    if (with_record) f.record = properties_support::p3_record(properties_support::p3_forms()[1].second);
    f.problem = equilibrium_support::flash(1);
    f.problem.phases = {equilibrium_support::phase("V", PhaseKind::Vapor)};
    f.problem.temperature = FreeValue::fixed(temperature);
    f.directions = {direction(20.0, 0.0, {0.0}), direction(-10.0, 0.0, {0.0}), direction(5.0, 0.0, {0.0})};
    f.observables = {obs(Kind::PhaseProperty, 0, 0, PropertyObservable::IdealEnthalpy), obs(Kind::PhaseProperty, 0, 0, PropertyObservable::TotalEnthalpy),
                     obs(Kind::PhaseProperty, 0, 0, PropertyObservable::IdealIsobaricHeatCapacity), obs(Kind::PhaseProperty, 0, 0, PropertyObservable::TotalIsobaricHeatCapacity)};
    f.scaled = false;
    return f;
}

}  // namespace e5

void check_e5_electrolyte() {
    using namespace e5;
    auto with_record = electrolyte_case(true);
    for (std::size_t i = 0; i < 3; ++i) with_record.observables.push_back(obs(Kind::PhaseLogActivity, 0, i));  // solvent and both ions
    qualify(with_record, "e5-5");
    // Residual-only without a record: residual siblings stay available, total rows are typed.
    auto residual_only = electrolyte_case(false);
    residual_only.observables.push_back(obs(Kind::PhaseProperty, 0, 0, PropertyObservable::TotalEnthalpy));
    const auto model = residual_only.model();
    const auto central = epcsaft::equilibrium::solve_equilibrium(model, residual_only.problem);
    const auto batch = actions(residual_only, model, central, residual_only.observables, 2);
    const bool ok = batch.results[0].status == SolvedStateActionStatus::Available && batch.results[0].action.has_value()
        && batch.results[1].status == SolvedStateActionStatus::Available && batch.results[1].action.has_value()
        && batch.results[2].status == SolvedStateActionStatus::MissingThermochemistry && !batch.results[2].action.has_value()
        && !batch.results[2].value.has_value() && batch.results[2].diagnostic_message.has_value();
    record_bool("e5-5/electrolyte-residual-only/sibling", "e5-5", "a missing thermal record types the total enthalpy row while residual siblings stay available",
        ok, recording_support::fields({recording_support::field("status_total", status_field(batch, 2)), recording_support::field("status_residual", status_field(batch, 0)),
            recording_support::field("message", batch.results[2].diagnostic_message.value_or("absent"))}));
    residual_only.observables.pop_back();
    qualify(residual_only, "e5-5");
}

void check_e5_roots() {
    using namespace e5;
    for (const PhaseKind kind : {PhaseKind::Vapor, PhaseKind::Liquid}) {
        const auto fixture = two_root_case(kind);
        const auto model = fixture.model();
        const auto central = epcsaft::equilibrium::solve_equilibrium(model, fixture.problem);
        const auto root = epcsaft::solve_density(model, 200.0, 217800.69494692178, {1.0},
            kind == PhaseKind::Vapor ? epcsaft::DensityBranch::Vapor : epcsaft::DensityBranch::Liquid, std::nullopt);
        const bool branch_ok = central.success && std::abs(std::log(central.molar_densities[0] / root.molar_density)) <= 1.0e-8
            && root.stable_root_count.value_or(0) == 2;
        record_bool((std::string("e5-6/") + fixture.label + "/central-root").c_str(), "e5-6",
            "the declared branch selects its own root of the two stable ethane roots", branch_ok,
            recording_support::fields({recording_support::field("rho", central.success ? central.molar_densities[0] : 0.0),
                recording_support::field("rho_branch", root.molar_density), recording_support::field("stable_roots", root.stable_root_count.value_or(0))}));
        qualify(fixture, "e5-6");
    }
}

void check_e5_thermal() {
    using namespace e5;
    qualify(thermal_case(320.0, true), "e5-7");
    const auto endpoint = thermal_case(250.0, true);
    const auto model = endpoint.model();
    const auto central = epcsaft::equilibrium::solve_equilibrium(model, endpoint.problem);
    const auto batch = actions(endpoint, model, central, endpoint.observables, 1);
    bool boundary = central.success;
    for (const auto& item : batch.results)
        boundary = boundary && item.value.has_value() && item.status == SolvedStateActionStatus::ThermalBoundary && !item.action.has_value();
    record_bool("e5-7/thermal-250/boundary", "e5-7", "at the closed interval endpoint central values are retained and thermal actions are typed thermal_boundary",
        boundary, recording_support::fields({recording_support::field("status", status_field(batch, 0)), recording_support::field("value", epcsaft::core::optional_number(batch.results[0].value))}));
    const auto outside = thermal_case(249.0, true);
    const auto outside_model = outside.model();
    const auto outside_central = epcsaft::equilibrium::solve_equilibrium(outside_model, outside.problem);
    const auto outside_batch = actions(outside, outside_model, outside_central, outside.observables, 1);
    bool interval = outside_central.success;
    for (const auto& item : outside_batch.results) interval = interval && item.status == SolvedStateActionStatus::OutsideInterval && !item.action.has_value();
    record_bool("e5-7/thermal-249/outside-interval", "e5-7", "below the closed interval dependent ideal and total rows are typed outside_interval",
        interval, recording_support::fields({recording_support::field("status", status_field(outside_batch, 0))}));
    const auto missing = thermal_case(320.0, false);
    const auto missing_model = missing.model();
    const auto missing_central = epcsaft::equilibrium::solve_equilibrium(missing_model, missing.problem);
    const auto missing_batch = actions(missing, missing_model, missing_central, missing.observables, 1);
    bool typed = missing_central.success;
    for (const auto& item : missing_batch.results) typed = typed && item.status == SolvedStateActionStatus::MissingThermochemistry;
    record_bool("e5-7/thermal-no-record/missing", "e5-7", "without a thermal record dependent rows are typed missing_thermochemistry",
        typed, recording_support::fields({recording_support::field("status", status_field(missing_batch, 0))}));
}

void check_e5_refusals() {
    using namespace e5;
    const auto fixture = neutral_case();
    const auto model = fixture.model();
    const auto compiled = epcsaft::equilibrium::compile_problem(model, fixture.problem);
    const auto central = epcsaft::equilibrium::solve_equilibrium(model, fixture.problem);
    CHECK(central.success);
    const auto invalid = [&](const char* id, const char* phrase, const SolvedStateActionRequest& request, const CompiledSystem& system, const EquilibriumResult& result) {
        std::string message;
        try { (void)epcsaft::equilibrium::solved_state_actions(system, result, request); }
        catch (const std::invalid_argument& error) { message = error.what(); }
        record_bool((std::string("e5-8/invalid/") + id).c_str(), "e5-8", phrase, message.find(phrase) != std::string::npos,
            recording_support::fields({recording_support::field("message", message)}));
    };
    SolvedStateActionRequest base;
    base.observables = {obs(Kind::PhaseAmount, 0)};
    base.directions = {fixture.directions[0]};
    auto empty = base; empty.observables.clear();
    invalid("no-observables", "no solved-state observables", empty, compiled, central);
    auto four = base; four.directions.assign(4, fixture.directions[0]);
    invalid("order-four", "solved-state action order", four, compiled, central);
    auto short_feed = base; short_feed.directions[0].feed_amounts_mol = {0.0, 0.0};
    invalid("feed-size", "solved-state action direction", short_feed, compiled, central);
    auto duplicate = base; duplicate.observables.push_back(obs(Kind::PhaseAmount, 0));
    invalid("duplicate", "duplicate solved-state observable", duplicate, compiled, central);
    auto ordinal = base; ordinal.observables = {obs(Kind::PhaseAmount, 5)};
    invalid("phase-ordinal", "solved-state observable phase ordinal", ordinal, compiled, central);
    auto unsupported = base; unsupported.observables = {obs(Kind::PhaseResidualChemicalPotentialOverRt, 0, 2)};
    auto restricted = fixture.problem; restricted.phases[0].support = {0, 1};
    const auto restricted_compiled = epcsaft::equilibrium::compile_problem(model, restricted);
    invalid("unsupported-potential", "residual chemical potential outside the phase support", unsupported, restricted_compiled, epcsaft::equilibrium::solve_equilibrium(model, restricted));
    unsupported.observables = {obs(Kind::PhaseLogActivity, 0, 2)};
    invalid("unsupported-log-activity", "log activity outside the phase support", unsupported, restricted_compiled, epcsaft::equilibrium::solve_equilibrium(model, restricted));
    auto free_t = base;
    const auto pinned = pinned_case(0.01);
    const auto pinned_model = pinned.model();
    const auto pinned_compiled = epcsaft::equilibrium::compile_problem(pinned_model, pinned.problem);
    const auto pinned_central = epcsaft::equilibrium::solve_equilibrium(pinned_model, pinned.problem);
    free_t.directions = {direction(0.0, 1.0, {0.0, 0.0})};
    invalid("free-pressure-control", "free pressure is a solved output", free_t, pinned_compiled, pinned_central);
    auto foreign = free_t; foreign.directions = {direction(0.0, 0.0, {0.0, 0.0})};
    invalid("mismatched-central", "central result does not match the compiled system", foreign, pinned_compiled, central);
    // Typed refusals with retained central values: a same-layout result of another problem
    // does not satisfy this problem's defining rows.
    auto other = fixture.problem; other.temperature.value = 181.0;
    const auto other_compiled = epcsaft::equilibrium::compile_problem(model, other);
    const auto other_batch = epcsaft::equilibrium::solved_state_actions(other_compiled, central, base);
    const double other_rows = epcsaft::equilibrium::maximum_abs(other_compiled.evaluate(central.coordinates).values);
    record_bool("e5-8/typed/defining-equations", "e5-8", "an accepted result of another same-layout problem is typed defining_equations_not_satisfied with its mapped value retained",
        other_rows > 1.0e-10 && other_batch.results[0].status == SolvedStateActionStatus::DefiningEquationsNotSatisfied
            && other_batch.results[0].value.has_value() && !other_batch.results[0].action.has_value(),
        recording_support::fields({recording_support::field("central_rows", other_rows), recording_support::field("status", status_field(other_batch, 0))}));
    auto exact = equilibrium_support::flash(2);
    exact.phases = {equilibrium_support::phase("L", PhaseKind::Liquid, AmountSpec::pinned(0.5)), equilibrium_support::phase("V", PhaseKind::Vapor, AmountSpec::pinned(0.5))};
    exact.temperature = FreeValue::variable(300.0); exact.pressure = FreeValue::variable(1.0e5);
    const auto& binary_model = pinned_model;
    const auto exact_compiled = epcsaft::equilibrium::compile_problem(binary_model, exact);
    const auto exact_central = epcsaft::equilibrium::solve_equilibrium(binary_model, exact);
    SolvedStateActionRequest exact_request;
    exact_request.observables = {obs(Kind::PhaseAmount, 0)};
    exact_request.directions = {direction(0.0, 0.0, {0.01, 0.0})};
    const auto not_run = epcsaft::equilibrium::solved_state_actions(exact_compiled, exact_central, exact_request);
    record_bool("e5-8/typed/central-unavailable", "e5-8", "the pre-solve NOT_RUN record is typed central_solve_unavailable",
        !exact_central.success && not_run.results[0].status == SolvedStateActionStatus::CentralSolveUnavailable && !not_run.results[0].value.has_value(),
        recording_support::fields({recording_support::field("solver_status", exact_central.solver_status)}));
    // Exact inventory faces admit central states; insertion directions change that support,
    // while a fixed intensive control retains its direct action.
    auto sparse = equilibrium_support::flash(3); sparse.feed.values = {0.5, 0.5, 0.0}; sparse.temperature.value = 180.0; sparse.pressure.value = 2.0e6;
    const auto sparse_compiled = epcsaft::equilibrium::compile_problem(model, sparse);
    const auto sparse_central = epcsaft::equilibrium::solve_equilibrium(model, sparse);
    SolvedStateActionRequest insertion;
    insertion.observables = {obs(Kind::PhaseAmount, 0), obs(Kind::Temperature)};
    insertion.directions = {direction(0.0, 0.0, {0.0, 0.0, 0.01}, {0.0, 0.0})};
    const auto inserted = epcsaft::equilibrium::solved_state_actions(sparse_compiled, sparse_central, insertion);
    record_bool("e5-8/typed/removed-species", "e5-8", "the admitted neutral face preserves the amount value while refusing insertion as support_change; fixed temperature remains available with zero action",
        sparse_central.success && !sparse_compiled.removed_species.empty() && sparse_central.mole_fractions[2] == 0
            && inserted.results[0].status == SolvedStateActionStatus::SupportChange
            && inserted.results[0].value.has_value() && !inserted.results[0].action.has_value()
            && inserted.results[1].status == SolvedStateActionStatus::Available
            && inserted.results[1].value == 180.0 && inserted.results[1].action == 0.0,
        recording_support::fields({recording_support::field("status", status_field(inserted, 0)), recording_support::field("message", sparse_central.message),
            recording_support::field("classified_residual_1", sparse_central.convergence.residual), recording_support::field("iterations", sparse_central.iterations),
            recording_support::field("mole_fractions_1", sparse_central.mole_fractions),
            recording_support::field("phase_amount_mol", epcsaft::core::optional_number(inserted.results[0].value)),
            recording_support::field("temperature_K", epcsaft::core::optional_number(inserted.results[1].value)),
            recording_support::field("temperature_action_K", epcsaft::core::optional_number(inserted.results[1].action)),
            recording_support::field("removed_species", sparse_compiled.removed_species.size())}));
    epcsaft::HelmholtzModel salt(properties_support::aqueous_electrolyte());
    auto water = equilibrium_support::flash(3); water.feed.values = {1.0, 0.0, 0.0}; water.temperature = FreeValue::fixed(298.15);
    const auto water_compiled = epcsaft::equilibrium::compile_problem(salt, water);
    const auto water_central = epcsaft::equilibrium::solve_equilibrium(salt, water);
    SolvedStateActionRequest ionic_request;
    ionic_request.observables = {obs(Kind::PhaseMolarDensity, 0)};
    ionic_request.directions = {direction(0.0, 0.0, {-0.02, 0.01, 0.01})};
    const auto ionic_batch = epcsaft::equilibrium::solved_state_actions(water_compiled, water_central, ionic_request);
    record_bool("e5-8/typed/removed-ions", "e5-8", "the admitted water face preserves the density value while refusing a charged insertion action as ionic_insertion",
        water_central.success && water_central.mole_fractions[1] == 0 && water_central.mole_fractions[2] == 0
            && ionic_batch.results[0].status == SolvedStateActionStatus::IonicInsertion
            && ionic_batch.results[0].value.has_value() && !ionic_batch.results[0].action.has_value(),
        recording_support::fields({recording_support::field("status", status_field(ionic_batch, 0)), recording_support::field("message", water_central.message),
            recording_support::field("classified_residual_1", water_central.convergence.residual), recording_support::field("iterations", water_central.iterations),
            recording_support::field("mole_fractions_1", water_central.mole_fractions),
            recording_support::field("density_mol_m3", epcsaft::core::optional_number(ionic_batch.results[0].value))}));
    // The matrix owner on a rank-deficient Jacobian: the exact all-pinned system loses one
    // rank at its generic point; nonfinite and excessive-residual paths remain unexecuted.
    const auto inspected = epcsaft::equilibrium::inspect_sensitivity(exact_compiled, exact_compiled.generic_point());
    record_bool("e5-8/matrix/rank-deficient", "e5-8", "the shared matrix owner reports the rank-deficient all-pinned Jacobian as singular with its rank below the dimension",
        !inspected.nonsingular && inspected.rank + 1 == exact_compiled.generic_point().size() && inspected.rank == exact_compiled.jacobian_rank,
        recording_support::fields({recording_support::field("rank", inspected.rank), recording_support::field("dimension", exact_compiled.generic_point().size()),
            recording_support::field("condition", inspected.condition_estimate)}));
    // Central preservation: an unavailable property action keeps the central value and its available siblings.
    auto partial = electrolyte_case(false);
    partial.observables = {obs(Kind::PhaseProperty, 0, 0, PropertyObservable::TotalEnthalpy), obs(Kind::PhaseMolarDensity, 0), obs(Kind::PhaseProperty, 0, 0, PropertyObservable::ResidualEnthalpy)};
    const auto partial_model = partial.model();
    const auto partial_central = epcsaft::equilibrium::solve_equilibrium(partial_model, partial.problem);
    const auto partial_batch = actions(partial, partial_model, partial_central, partial.observables, 1);
    record_bool("e5-8/typed/sibling-preservation", "e5-8", "an unavailable ideal dependency leaves the density and residual siblings available with finite values",
        partial_batch.results[0].status == SolvedStateActionStatus::MissingThermochemistry && partial_batch.results[1].status == SolvedStateActionStatus::Available
            && partial_batch.results[1].action.has_value() && partial_batch.results[2].status == SolvedStateActionStatus::Available && partial_batch.results[2].action.has_value(),
        recording_support::fields({recording_support::field("statuses", status_field(partial_batch, 0) + "," + status_field(partial_batch, 1) + "," + status_field(partial_batch, 2))}));
    // Forward-error bound 2 cond_inf eta on the equilibrated system (#130). An identical m = 1 pair at
    // its bubble point (zero-pinned vapor) has the exact feed response dP = drho = 0, dx_V = dx_L = dz,
    // dN_L = sum dn. At fixed T/(eps/k) every energy scale E is the same reduced state; only the raw-Pa
    // pressure column grows, so the raw cond_inf(J) grows with E while the equilibrated one must not.
    // The error is measured in the certified norm |C^-1 dy|/|C^-1 y| over the observable pinned-chart
    // coordinates [ln x_L, N_L/N_feed (feed total 1), ln v_L, ln x_V, ln v_V, P]; C holds the
    // power-of-two inverse column maxima of J. The mu/RT multipliers are not observables. The computed
    // J carries its own evaluation error dJ; the exact response z* has backward error eta_J in the computed
    // system, so the first-order bound is 2 cond_inf (eta + eta_J) (#272).
    const double limit = epcsaft::policy::EQUILIBRIUM_ACTION_FORWARD_ERROR.value;
    const std::size_t column[] = {7, 0, 1, 4, 5, 2, 3, 6};  // chart column of each observable below
    for (int rung = 0; rung < 16; ++rung) {
        const double energy = 1.0e4 * std::pow(2.0, 0.25 * rung);
        auto bubble = pinned_case(0.0);
        bubble.parameters = properties_support::parameters(2);
        bubble.parameters.segment_count = {1.0, 1.0};
        bubble.parameters.sigma_constant = {2.0, 2.0};
        bubble.parameters.dispersion_energy_over_k = {energy, energy};
        bubble.parameters.molar_mass_kg_per_mol = {0.03, 0.03};
        bubble.problem.temperature = FreeValue::fixed(0.7 * energy);
        bubble.directions = {direction(0.0, 0.0, {0.013, -0.007})};
        bubble.observables = {obs(Kind::Pressure), obs(Kind::PhaseMoleFraction, 0, 0), obs(Kind::PhaseMoleFraction, 0, 1), obs(Kind::PhaseMoleFraction, 1, 0),
            obs(Kind::PhaseMoleFraction, 1, 1), obs(Kind::PhaseAmount, 0), obs(Kind::PhaseMolarDensity, 0), obs(Kind::PhaseMolarDensity, 1)};
        const double exact[] = {0.0, 0.010, -0.010, 0.010, -0.010, 0.006, 0.0, 0.0};
        const auto bubble_model = bubble.model();
        const auto bubble_central = epcsaft::equilibrium::solve_equilibrium(bubble_model, bubble.problem);
        const auto batch = actions(bubble, bubble_model, bubble_central, bubble.observables, 1);
        const auto& diagnostics = batch.diagnostics;
        // An ASL host (#173): AMPL's fpinit sets the x87 precision control to 53 bits (control word
        // 0x27f), as pynumero's PyomoNLP does in a live campaign. The action must not depend on it.
        fpu_control_t clean_word, asl_word = 0x27f;
        _FPU_GETCW(clean_word);
        _FPU_SETCW(asl_word);
        const auto host = actions(bubble, bubble_model, bubble_central, bubble.observables, 1);
        _FPU_SETCW(clean_word);
        const auto n = static_cast<Eigen::Index>(bubble_central.coordinates.size());
        const auto values = epcsaft::equilibrium::compile_problem(bubble_model, bubble.problem).jacobian(bubble_central.coordinates);
        const Eigen::MatrixXd jacobian = Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(values.data(), n, n);
        const double raw_condition = jacobian.cwiseAbs().rowwise().sum().maxCoeff() * jacobian.inverse().cwiseAbs().rowwise().sum().maxCoeff();
        const bool certified = diagnostics.condition_estimate.has_value() && diagnostics.maximum_linear_residual.has_value();
        const double bound = certified ? 2.0 * *diagnostics.condition_estimate * *diagnostics.maximum_linear_residual : 0.0;
        double error = 0.0, norm = 0.0, shift = 0.0;
        Eigen::VectorXd z = Eigen::VectorXd::Zero(n);  // z*, chart coordinates
        bool identical = host.results.size() == std::size(exact);
        bool typed = bubble_central.success && certified && bound <= limit;
        for (std::size_t i = 0; i < std::size(exact); ++i) {
            const auto& item = batch.results[i];
            typed = typed && item.status == SolvedStateActionStatus::Available && item.value.has_value() && item.action.has_value();
            // Chart coordinate of the observable action (log coordinates for x and v), times C^-1.
            const double chart = i == 0 || i == 5 ? 1.0 : (i >= 6 ? -1.0 : 1.0) / item.value.value_or(1.0);
            const double weight = chart * std::exp2(std::round(std::log2(jacobian.col(static_cast<Eigen::Index>(column[i])).cwiseAbs().maxCoeff())));
            error = std::max(error, std::abs(weight * (item.action.value_or(0.0) - exact[i])));
            z(static_cast<Eigen::Index>(column[i])) = chart * exact[i];
            shift = std::max(shift, std::abs(weight * (host.results[i].action.value_or(0.0) - item.action.value_or(0.0))));
            identical = identical && item.action.has_value() && host.results[i].action == item.action;
            norm = std::max(norm, std::abs(weight * exact[i]));
        }
        // The first-order bound excludes the final rounding of each returned double observable.
        const double floor = std::numeric_limits<double>::epsilon();
        // dJ: the phase rows' composition entries are exactly +-1/2 for identical species. The multiplier
        // columns of dJ are zero, so z* needs no multiplier entries; eta_J is normed like the error.
        // eta_J excuses only rounding-level dJ; a gross composition-derivative defect belongs to the exact-derivative checks.
        Eigen::MatrixXd deviation = Eigen::MatrixXd::Zero(n, n);
        for (Eigen::Index r = 0; r < 6; ++r) for (const Eigen::Index c : {0, 1, 4, 5})
            if (jacobian(r, c) != 0.0) deviation(r, c) = jacobian(r, c) - std::copysign(0.5, jacobian(r, c));
        const auto power = [](double v) { return v > 0.0 ? std::exp2(-std::round(std::log2(v))) : 1.0; };
        const Eigen::VectorXd col = jacobian.cwiseAbs().colwise().maxCoeff().transpose().unaryExpr(power);
        const Eigen::VectorXd row = (jacobian * col.asDiagonal()).cwiseAbs().rowwise().maxCoeff().unaryExpr(power);
        const double eta_jacobian = (row.asDiagonal() * deviation * z).lpNorm<Eigen::Infinity>()
            / ((row.asDiagonal() * jacobian * col.asDiagonal()).cwiseAbs().rowwise().sum().maxCoeff() * norm
                + (row.asDiagonal() * jacobian * z).lpNorm<Eigen::Infinity>());
        const double allowed = bound + 2.0 * diagnostics.condition_estimate.value_or(0.0) * eta_jacobian;
        const bool verdict = typed && error <= (allowed + floor) * norm;
        const auto id = "e5-8/forward-bound/rung" + std::to_string(rung);
        recording_support::compare_custom(id.c_str(), "forward_error", "1", error / norm, allowed,
            recording_support::fields({recording_support::field("energy_over_k", energy), recording_support::field("temperature", 0.7 * energy),
                recording_support::field("pressure", epcsaft::core::optional_number(batch.results[0].value)), recording_support::field("status", status_field(batch, 0)),
                recording_support::field("raw_condition", raw_condition),
                recording_support::field("condition", epcsaft::core::optional_number(diagnostics.condition_estimate)),
                recording_support::field("backward_error", epcsaft::core::optional_number(diagnostics.maximum_linear_residual)),
                recording_support::field("jacobian_backward_error", eta_jacobian),
                recording_support::field("bound", allowed), recording_support::field("limit", limit), recording_support::field("reference_norm", norm)}),
            "Relationship", "analytic-bubble", "a badly scaled but regular state is admitted and its action agrees with the exact bubble response within 2 cond (eta + eta_J) in the equilibrated norm",
            "available, 2 kappa eta <= limit, |C^-1 dy|/|C^-1 y| <= 2 kappa (eta + eta_J) + eps", verdict, "e5-8", id.c_str(), "|C^{-1}dy|/|C^{-1}y|", "E",
            "|C^{-1}dy|/|C^{-1}y|\\leq2\\kappa_\\infty(\\eta+\\eta_J)+\\epsilon", "case");
        const auto host_id = "e5-8/asl-host/rung" + std::to_string(rung);
        record_bool(host_id.c_str(), "e5-8", "the action under the ASL x87 control word 0x27f is the clean-process action exactly",
            identical, recording_support::fields({recording_support::field("energy_over_k", energy), recording_support::field("control_word", 0x27f),
                recording_support::field("shift", shift / norm), recording_support::field("clean_error", error / norm)}));
    }
}

// S3 (#53): exact identities of solved-state actions on manufactured systems. Nonzero references
// use the retained CRITERION, zeros the retained STRUCTURAL_ZERO; no stencil and no new threshold.
void check_e5_identities() {
    using namespace e5;
    const auto within = [&](const std::string& id, const char* claim, double value, double reference, double tolerance,
                            const SolvedStateActionBatch& batch, const std::string& stated) {
        const bool verdict = std::isfinite(value) && std::abs(value - reference) <= tolerance;
        recording_support::compare(id.c_str(), "action_identity", "1", value, reference,
            recording_support::fields({stated, recording_support::field("tolerance", tolerance),
                recording_support::field("condition_estimate", batch.diagnostics.condition_estimate.value_or(0.0)),
                recording_support::field("linear_residual", batch.diagnostics.maximum_linear_residual.value_or(0.0))}),
            "Relationship", "exact-identity", claim, "absolute", tolerance, 0.0, verdict, "e5-10");
    };
    const auto action = [](const SolvedStateActionBatch& batch, std::size_t i) { return batch.results[i].action.value_or(std::nan("")); };

    // H3: Euler homogeneity in the feed amounts, radial and mixed with a non-radial direction.
    auto neutral = neutral_case();
    const auto neutral_model = neutral.model();
    const auto neutral_central = epcsaft::equilibrium::solve_equilibrium(neutral_model, neutral.problem);
    neutral.directions = {direction(0.0, 0.0, neutral.problem.feed.values, {0.0, 0.0}), direction(0.0, 0.0, {-0.01, 0.006, 0.004}, {0.0, 0.0})};
    std::vector<SolvedStateObservable> list;
    std::vector<bool> extensive;
    for (std::size_t phase = 0; phase < 2; ++phase) {
        list.push_back(obs(Kind::PhaseAmount, phase)); extensive.push_back(true);
        list.push_back(obs(Kind::PhaseMolarDensity, phase)); extensive.push_back(false);
        for (std::size_t i = 0; i < 3; ++i) {
            list.push_back(obs(Kind::PhaseComponentAmount, phase, i)); extensive.push_back(true);
            list.push_back(obs(Kind::PhaseMoleFraction, phase, i)); extensive.push_back(false);
            list.push_back(obs(Kind::PhaseLogActivity, phase, i)); extensive.push_back(false);
        }
    }
    const auto euler = recording_support::fields({recording_support::field("n_feed", neutral.directions[0].feed_amounts_mol),
        recording_support::field("v", neutral.directions[1].feed_amounts_mol)});
    const auto radial = actions(neutral, neutral_model, neutral_central, list, 1);
    const auto mixed = actions(neutral, neutral_model, neutral_central, list, 2);
    std::swap(neutral.directions[0], neutral.directions[1]);
    const auto skew = actions(neutral, neutral_model, neutral_central, list, 1);
    double wrong_sign = 0.0;
    for (std::size_t i = 0; i < list.size(); ++i) {
        const std::string id = "e5-10/euler/" + epcsaft::equilibrium::solved_state_observable_id(list[i]);
        const double value = radial.results[i].value.value_or(std::nan("")), first = action(skew, i);
        if (extensive[i]) {
            within(id + "/radial", "degree one: D y[n] = y", action(radial, i), value, CRITERION[0] * std::abs(value), radial, euler);
            within(id + "/mixed", "degree one: D2 y[n, v] = 0", action(mixed, i), 0.0, STRUCTURAL_ZERO * std::max(1.0, std::abs(first)), mixed, euler);
        } else {
            within(id + "/radial", "degree zero: D y[n] = 0", action(radial, i), 0.0, STRUCTURAL_ZERO * std::max(1.0, std::abs(value)), radial, euler);
            within(id + "/mixed", "degree zero: D2 y[n, v] = -D y[v]", action(mixed, i), -first, CRITERION[1] * std::abs(first), mixed, euler);
            wrong_sign = std::max(wrong_sign, std::abs(action(mixed, i) - first) / std::abs(first));
        }
    }
    record_bool("e5-10/euler/wrong-sign", "e5-10", "the mixed identity with the wrong sign misses by a relative two",
        wrong_sign > 1.0, recording_support::fields({recording_support::field("relative_defect", wrong_sign)}));

    // H4: in each liquid of A + B = C, sum_i nu_i ln a_i = ln K(T) along every admitted control path.
    const std::vector<double> nu = {-1.0, -1.0, 1.0};
    const double b = 1500.0, t0 = 298.15;  // labelled manufactured EOS-basis ln K = b (1/T - 1/T0), K(T0) = 1
    auto reactive = reactive_case();
    reactive.problem.standard_gibbs_j_per_mol.clear();
    reactive.problem.reactions[0].equilibrium_constant.reset();
    reactive.problem.reactions[0].correlation = epcsaft::ReactionLogPolynomial{-b / t0, b, 0.0, 0.0, t0, 250.0, 450.0,
        epcsaft::equilibrium::EOS_STANDARD_STATE_ID};
    const auto reactive_model = reactive.model();
    const auto reactive_central = epcsaft::equilibrium::solve_equilibrium(reactive_model, reactive.problem);
    CHECK(reactive_central.success);
    list.clear();
    for (std::size_t phase = 0; phase < 2; ++phase) for (std::size_t i = 0; i < 3; ++i) list.push_back(obs(Kind::PhaseLogActivity, phase, i));
    const auto manifold = [&](const char* label, std::vector<SolvedStateActionDirection> directions, double reference, bool values) {
        reactive.directions = std::move(directions);
        const std::size_t order = reactive.directions.size();
        const auto batch = actions(reactive, reactive_model, reactive_central, list, order);
        std::vector<double> d_t, d_p;
        for (const auto& item : reactive.directions) { d_t.push_back(item.temperature_k); d_p.push_back(item.pressure_pa); }
        const auto stated = recording_support::fields({recording_support::field("b", b), recording_support::field("T0", t0),
            recording_support::field("dT", d_t), recording_support::field("dP", d_p), recording_support::field("feed_direction", reactive.directions[0].feed_amounts_mol)});
        for (std::size_t phase = 0; phase < 2; ++phase) {
            double sum = 0.0, scale = 0.0, perturbed = 0.0;
            for (std::size_t i = 0; i < 3; ++i) {
                const double term = values ? batch.results[3 * phase + i].value.value_or(std::nan("")) : action(batch, 3 * phase + i);
                sum += nu[i] * term; scale += std::abs(nu[i] * term); perturbed += (i == 2 ? 2.0 : nu[i]) * term;
            }
            const std::string id = std::string("e5-10/manifold/") + label + "/phase" + std::to_string(phase);
            // The central identity holds to the retained affinity bound of e5-3; actions to the retained criterion.
            const double tolerance = values ? 1.0e-9 : CRITERION[order - 1] * scale;
            within(id, "sum_i nu_i ln a_i follows ln K(T) in each phase", sum, reference, tolerance, batch, stated);
            record_bool((id + "/perturbed-nu").c_str(), "e5-10", "a wrong stoichiometric coefficient fails the identity's own criterion",
                std::abs(perturbed - reference) > tolerance, recording_support::fields({stated, recording_support::field("defect", perturbed - reference),
                    recording_support::field("tolerance", tolerance), recording_support::field("defect_over_tolerance", std::abs(perturbed - reference) / tolerance),
                    recording_support::field("term_scale", scale)}));
        }
        return batch;
    };
    // Two phases, three species and one reaction leave no intensive freedom at fixed T and P, so a feed
    // direction moves only the phase amounts: the manifold identity would read 0 = 0 there. Pressure and
    // temperature directions do move the intensive state; K does not depend on pressure.
    const auto v = direction(0.0, 0.0, {0.001, -0.001, 0.0}), dP = direction(0.0, 1.0e5, {0.0, 0.0, 0.0}), dT = direction(1.0, 0.0, {0.0, 0.0, 0.0});
    const auto invariant = manifold("value", {v}, 0.0, true);
    manifold("pressure", {dP}, 0.0, false);
    manifold("pressure-temperature", {dP, dT}, 0.0, false);
    manifold("temperature", {dT}, -b / (t0 * t0), false);
    manifold("temperature-temperature", {dT, dT}, 2.0 * b / (t0 * t0 * t0), false);
    // Third order is where the weight j of the series logarithm first matters.
    manifold("temperature-temperature-temperature", {dT, dT, dT}, -6.0 * b / (t0 * t0 * t0 * t0), false);
    manifold("pressure-temperature-temperature", {dP, dT, dT}, 0.0, false);
    for (std::size_t i = 0; i < list.size(); ++i) {
        within("e5-10/phase-rule/" + epcsaft::equilibrium::solved_state_observable_id(list[i]), "a feed direction leaves the invariant intensive state unchanged",
               action(invariant, i), 0.0, STRUCTURAL_ZERO * std::max(1.0, std::abs(invariant.results[i].value.value_or(std::nan("")))), invariant,
               recording_support::fields({recording_support::field("feed_direction", v.feed_amounts_mol)}));
    }
}

}  // namespace equilibrium_checks
