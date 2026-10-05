#include "regression/observations.hpp"

#include "core/numerical_policy.hpp"
#include "properties/thermochemistry.hpp"
#include <algorithm>
#include <cmath>
#include <set>
#include <numeric>
#include <stdexcept>
#include <string>
#include <tuple>

namespace epcsaft::regression {
namespace {

using equilibrium::SolvedStateObservable;
using Kind = equilibrium::SolvedStateObservableKind;

void require(bool condition, const char* what) {
    if (!condition) throw std::invalid_argument(std::string("observation ") + what);
}

SolvedStateObservable observable(Kind kind, std::size_t phase, std::size_t component = 0) {
    SolvedStateObservable result;
    result.kind = kind;
    result.phase = phase;
    result.component = component;
    return result;
}

template<class X, class Fields>
bool optional_equal(const std::optional<X>& x, const std::optional<X>& y, const Fields& fields) {
    return x.has_value() == y.has_value() && (!x || fields(*x) == fields(*y));
}

double mean_molar_mass(const Parameters& parameters, const double* composition, bool& missing) {
    double result = 0.0;
    for (std::size_t i = 0; i < parameters.component_count(); ++i) {
        const double mass = parameters.molar_mass_kg_per_mol[i];
        if (composition[i] > 0.0 && !(mass > 0.0)) missing = true;
        result += composition[i] * mass;
    }
    return result;
}

double component_amount(const equilibrium::EquilibriumResult& state, std::size_t phase,
                        const std::vector<std::size_t>& species) {
    double value = 0.0;
    for (const auto i : species) value += state.amounts[phase] * state.mole_fractions[phase * (state.mole_fractions.size() / state.phase_count) + i];
    return value;
}

std::vector<SolvedStateObservable> enthalpy_observables(const equilibrium::EquilibriumProblem& problem) {
    std::vector<SolvedStateObservable> result;
    for (std::size_t p = 0; p < problem.phases.size(); ++p) {
        const auto& amount = problem.phases[p].amount;
        if (amount.kind == equilibrium::AmountKind::Pinned && amount.fraction == 0.0) continue;
        result.push_back(observable(Kind::PhaseAmount, p));
        result.push_back({Kind::PhaseProperty, p, 0, PropertyObservable::TotalEnthalpy});
    }
    return result;
}

void finish(Prediction& prediction, const Observation& observation, double value) {
    if (!std::isfinite(value)) {
        prediction.status = ObservationStatus::Nonfinite;
        prediction.message = "the predicted value is not finite";
        return;
    }
    prediction.value = value;  // an independently valid value survives an unavailable residual
    prediction.status = ObservationStatus::Available;
    prediction.residual = residual(observation, value, prediction.status);
    if (prediction.status == ObservationStatus::Nonfinite) prediction.message = "the residual is not finite";
}

}  // namespace

const char* observation_status_name(ObservationStatus status) noexcept {
    switch (status) {
    case ObservationStatus::NotEvaluated: return "not_evaluated";
    case ObservationStatus::ReferenceUnavailable: return "reference_unavailable";
    case ObservationStatus::Available: return "available";
    case ObservationStatus::ThermochemistryUnavailable: return "thermochemistry_unavailable";
    case ObservationStatus::DensityUnavailable: return "density_unavailable";
    case ObservationStatus::CentralSolveUnavailable: return "central_solve_unavailable";
    case ObservationStatus::TopologyEvent: return "topology_event";
    case ObservationStatus::MissingMolarMass: return "missing_molar_mass";
    case ObservationStatus::NonpositiveLogArgument: return "nonpositive_log_argument";
    case ObservationStatus::ActionUnavailable: return "action_unavailable";
    case ObservationStatus::Nonfinite: return "nonfinite";
    case ObservationStatus::InadmissibleTrial: return "inadmissible_trial";
    }
    return "unknown";
}

std::optional<double> residual(const Observation& observation, double value, ObservationStatus& status) {
    double result = 0.0;
    if (observation.form == ResidualForm::Difference) {
        result = (value - observation.observed) / observation.scale;
    } else if (!(value > 0.0) || !(observation.observed > 0.0)) {
        status = ObservationStatus::NonpositiveLogArgument;
        return std::nullopt;
    } else {
        // log1p keeps full accuracy near unity; the difference of logarithms covers the dynamic
        // range where the ratio itself would overflow or underflow.
        const double gap = value - observation.observed;
        result = std::abs(gap) < policy::LOG_RATIO_NEAR_UNITY.value * observation.observed
            ? std::log1p(gap / observation.observed) : std::log(value) - std::log(observation.observed);
    }
    if (std::isfinite(result)) return result;
    status = ObservationStatus::Nonfinite;
    return std::nullopt;
}

double residual_slope(const Observation& observation, double value) {
    return observation.form == ResidualForm::Difference ? 1.0 / observation.scale : 1.0 / value;
}

bool same_state(const Observation& left, const Observation& right) {
    if (electrolyte_observation(left.kind) || electrolyte_observation(right.kind))
        return electrolyte_observation(left.kind) && electrolyte_observation(right.kind)
            && (left.kind == ObservationKind::BulkRelativePermittivity) == (right.kind == ObservationKind::BulkRelativePermittivity)
            && left.temperature == right.temperature && left.pressure == right.pressure && left.branch == right.branch
            && left.composition == right.composition && left.ion_indices == right.ion_indices
            && left.stoichiometry == right.stoichiometry && left.molality == right.molality
            && left.reference_composition == right.reference_composition;
    // Pressure closures stay singletons; equilibrium projection and residual metadata do not select physics.
    if (pressure_closed(left.kind) || pressure_closed(right.kind)
        || left.kind == ObservationKind::LogEquilibriumConstant || right.kind == ObservationKind::LogEquilibriumConstant
        || left.kind == ObservationKind::AdditionEnthalpy || right.kind == ObservationKind::AdditionEnthalpy) return false;
    const auto& a = left.problem;
    const auto& b = right.problem;
    const auto intensive = [](const auto& v) { return std::tie(v.free, v.value, v.absent); };
    const auto feed = [](const auto& v) { return std::tie(v.kind, v.values, v.solvent, v.solvent_mass_kg); };
    const auto phase = [](const auto& v) {
        return std::tie(v.name, v.kind, v.support, v.support_declared, v.amount.kind, v.amount.fraction,
                        v.composition_guess, v.packing_guess);
    };
    const auto source = [](const auto& v) {
        return std::tie(v.solvent_index, v.solute_molality_ratios, v.standard_molality_mol_per_kg);
    };
    if (intensive(a.temperature) != intensive(b.temperature) || intensive(a.pressure) != intensive(b.pressure)
        || feed(a.feed) != feed(b.feed) || a.standard_gibbs_j_per_mol != b.standard_gibbs_j_per_mol
        || a.topology_declared != b.topology_declared || a.row_tolerance != b.row_tolerance
        || a.phases.size() != b.phases.size() || a.reactions.size() != b.reactions.size()
        || !optional_equal(a.source_reference, b.source_reference, source)) return false;
    for (std::size_t i = 0; i < a.phases.size(); ++i)
        if (phase(a.phases[i]) != phase(b.phases[i])) return false;
    for (std::size_t i = 0; i < a.reactions.size(); ++i)
        if (!same_reaction(a.reactions[i], b.reactions[i])) return false;
    return true;
}

void validate(const Observation& observation, const HelmholtzModel& model,
              const equilibrium::CompiledSystem* prepared) {
    const auto components = model.parameters().component_count();
    require(std::isfinite(observation.observed), "observed value");
    require(observation.form == ResidualForm::LogRatio
            || (std::isfinite(observation.scale) && observation.scale > 0.0), "scale");
    require(!observation.mass_basis || observation.kind == ObservationKind::Density
            || observation.kind == ObservationKind::PhaseDensity, "mass basis is defined only for a density");
    if (electrolyte_observation(observation.kind)) {
        require(std::isfinite(observation.temperature) && observation.temperature > 0.0
            && std::isfinite(observation.pressure) && observation.pressure > 0.0, "positive finite T/P");
        require(observation.branch == DensityBranch::Liquid && observation.composition.size() == components, "liquid composition");
        require(std::isfinite(observation.molality) && observation.molality >= 0.0, "molality");
        if (observation.kind == ObservationKind::BulkRelativePermittivity) {
            const auto admit = [&](const auto& x) {
                require(x.size() == components && std::all_of(x.begin(), x.end(), [](double v) { return std::isfinite(v) && v >= 0.0; })
                    && std::abs(std::accumulate(x.begin(), x.end(), 0.0) - 1.0) <= policy::COMPOSITION_SUM_TOLERANCE.value, "permittivity composition");
            };
            admit(observation.composition);
            if (!observation.reference_composition.empty()) admit(observation.reference_composition);
        } else {
            (void)formula_composition(observation.composition, model.parameters(), observation.ion_indices, observation.stoichiometry, 1.0);
            if (observation.kind == ObservationKind::OsmoticCoefficient)
                require(observation.molality >= policy::OSMOTIC_MOLALITY_FLOOR.value
                    && std::count_if(observation.composition.begin(), observation.composition.end(), [](double x) { return x > 0.0; }) == 1, "osmotic domain");
            if (observation.kind == ObservationKind::SolvationGibbsEnergy || observation.kind == ObservationKind::TransferGibbsEnergy)
                require(observation.species.size() == 1 && (observation.species[0] == observation.ion_indices[0]
                    || observation.species[0] == observation.ion_indices[1]), "selected ion");
            if (observation.kind == ObservationKind::TransferGibbsEnergy)
                (void)formula_composition(observation.reference_composition, model.parameters(), observation.ion_indices, observation.stoichiometry, 1.0);
        }
        require(observation.form != ResidualForm::LogRatio || observation.observed > 0.0, "positive LogRatio observation");
        return;
    }
    if (observation.kind == ObservationKind::LogEquilibriumConstant) {
        require(!observation.problem.temperature.free && !observation.problem.temperature.absent && std::isfinite(observation.problem.temperature.value)
            && observation.problem.temperature.value > 0.0 && observation.form == ResidualForm::Difference, "ln K needs fixed positive T and Difference");
        return validate(CoordinateSet{{{{}, 0.0, 1.0, -HUGE_VAL, HUGE_VAL, ReactionCoordinate{observation.reaction}}}}, model.parameters(), {observation}, model.thermochemistry());
    }
    if (pressure_closed(observation.kind) || observation.kind == ObservationKind::AdditionEnthalpy) {
        require(std::isfinite(observation.temperature) && observation.temperature > 0.0, "temperature");
        require(std::isfinite(observation.pressure) && observation.pressure > 0.0, "pressure");
        require(observation.composition.size() == components, "composition");
        if (observation.kind != ObservationKind::AdditionEnthalpy) return;
    }
    const auto owned = prepared ? std::optional<equilibrium::CompiledSystem>()
                                : equilibrium::compile_problem(model, observation.problem);
    const auto& compiled = prepared ? *prepared : *owned;
    if (observation.kind == ObservationKind::AdditionEnthalpy) {
        require(std::isfinite(observation.feed_amount) && observation.feed_amount > 0.0, "feed amount");
        const auto final = equilibrium::compile_problem(model, observation.final_problem);
        require(compiled.problem.reactions.size() == final.problem.reactions.size(), "endpoint reaction count");
        for (std::size_t r = 0; r < compiled.problem.reactions.size(); ++r)
            require(same_reaction(compiled.problem.reactions[r], final.problem.reactions[r]), "endpoint reaction declaration");
        const auto source = [](const auto& v) { return std::tie(v.solvent_index, v.solute_molality_ratios, v.standard_molality_mol_per_kg); };
        require(optional_equal(compiled.problem.source_reference, final.problem.source_reference, source), "endpoint source reference");
        for (std::size_t i = 0; i < components; ++i)
            require(std::abs(final.n_feed[i] - compiled.n_feed[i] - observation.feed_amount * observation.composition[i])
                <= policy::ADDITION_MATERIAL_BALANCE.value * std::max({1.0, std::abs(final.n_feed[i]), std::abs(compiled.n_feed[i])}),
                ("material balance component " + std::to_string(i)).c_str());
        return;
    }
    require(observation.phase < compiled.problem.phases.size(), "phase index");
    if (observation.kind != ObservationKind::PhaseMoleFraction && observation.kind != ObservationKind::PartialPressure
        && observation.kind != ObservationKind::PhaseAmountRatio) return;
    for (const auto& [phase, species_set] : {std::pair{observation.phase, observation.species},
             std::pair{observation.denominator_phase, observation.kind == ObservationKind::PhaseAmountRatio
                 ? observation.denominator_species : observation.species}}) {
        require(phase < compiled.problem.phases.size(), "denominator phase index");
        require(!species_set.empty(), "species set is empty");
        const auto& support = compiled.problem.phases[phase].support;
        std::set<std::size_t> seen;
        for (const auto species : species_set) {
            require(seen.insert(species).second, "species repeated");
            require(std::find(support.begin(), support.end(), species) != support.end(), "species outside the phase support");
        }
        if (observation.kind != ObservationKind::PhaseAmountRatio) break;
    }
}

Prediction predict(
    const HelmholtzModel& model, const Observation& observation,
    const std::optional<equilibrium::EquilibriumResult>& start,
    std::optional<equilibrium::CompiledSystem> prepared, ReferenceCache* references
) {
    if (electrolyte_observation(observation.kind)) {
        Prediction result;
        try {
            if (observation.kind == ObservationKind::BulkRelativePermittivity) {
                result.state = model.state_at_pressure(observation.temperature, observation.pressure, observation.composition, observation.branch);
                if (!observation.reference_composition.empty()) result.reference_state = model.state_at_pressure(
                    observation.temperature, observation.pressure, observation.reference_composition, observation.branch);
            } else {
                ReferenceCache local;
                auto& cache = references ? *references : local;
                const auto reference = [&](const auto& x) {
                    auto key = observation;
                    key.composition = x; key.molality = 0.0; key.reference_composition.clear();
                    for (const auto& entry : cache) if (same_state(key, entry.first)) return entry.second;
                    auto value = infinite_dilution_reference(model, key.temperature, key.pressure, x, key.ion_indices, key.stoichiometry);
                    cache.emplace_back(key, value);
                    return value;
                };
                result.reference = reference(observation.composition);
                if (observation.kind == ObservationKind::TransferGibbsEnergy) result.second_reference = reference(observation.reference_composition);
                if (observation.kind == ObservationKind::MeanIonicActivityCoefficient || observation.kind == ObservationKind::OsmoticCoefficient) {
                    if (observation.molality > 0.0) result.activity = electrolyte_activity(model, observation.temperature, observation.pressure,
                        observation.composition, observation.ion_indices, observation.stoichiometry, observation.molality, *result.reference);
                }
            }
            project_prediction(model, observation, result);
        } catch (const std::exception& error) {
            result.status = observation.kind == ObservationKind::BulkRelativePermittivity
                ? ObservationStatus::DensityUnavailable : ObservationStatus::ReferenceUnavailable;
            result.message = error.what();
        }
        return result;
    }
    if (observation.kind == ObservationKind::AdditionEnthalpy) {
        auto item = observation;
        item.kind = ObservationKind::SystemEnthalpy;
        item.form = ResidualForm::Difference; item.scale = 1.0; item.observed = 0.0;
        const auto failed = [](Prediction& p, const char* source) {
            if (p.status == ObservationStatus::Available) return false;
            p.message = std::string(source) + ":" + observation_status_name(p.status) + ": " + p.message;
            return true;
        };
        auto initial = predict(model, item, start, std::move(prepared));
        if (failed(initial, "initial")) return initial;
        item.problem = observation.final_problem;
        auto compiled = equilibrium::compile_problem(model, item.problem);
        const auto& a = initial.compiled->problem.phases;
        const auto& b = compiled.problem.phases;
        const bool compatible = a.size() == b.size() && std::equal(a.begin(), a.end(), b.begin(),
            [](const auto& x, const auto& y) { return x.kind == y.kind && x.support == y.support; });
        auto final = predict(model, item, compatible ? initial.solved : std::nullopt, std::move(compiled));
        if (failed(final, "final")) return final;
        item.kind = ObservationKind::FeedEnthalpy;
        auto feed = predict(model, item);
        if (failed(feed, "feed")) return feed;
        finish(initial, observation, (*final.value - *initial.value) / observation.feed_amount - *feed.value);
        initial.unit = "J/mol";
        initial.final_solved = std::move(final.solved); initial.final_compiled = std::move(final.compiled);
        return initial;
    }
    Prediction result;
    result.compiled = std::move(prepared);
    const auto& parameters = model.parameters();
    bool missing_mass = false;
    if (pressure_closed(observation.kind)) {
        result.unit = observation.kind == ObservationKind::FeedEnthalpy ? "J/mol"
            : observation.kind == ObservationKind::IsobaricHeatCapacity ? "J/(mol K)"
            : observation.mass_basis ? "kg/m^3" : "mol/m^3";
        State state;
        try {
            state = model.state_at_pressure(observation.temperature, observation.pressure,
                                            observation.composition, observation.branch);
        } catch (const core::ConvergenceFailure& failure) {
            result.status = ObservationStatus::DensityUnavailable;
            result.convergence = failure.record;
            result.message = failure.what();
            return result;
        } catch (const std::runtime_error& failure) {
            result.status = ObservationStatus::DensityUnavailable;
            result.message = failure.what();
            return result;
        }
        result.convergence = state.density_convergence;
        result.branch = state.density_branch;
        result.stable_root_count = state.stable_root_count;
        std::optional<double> caloric = state.total_isobaric_heat_capacity;
        auto thermal_status = state.thermochemistry_status;
        if (observation.kind == ObservationKind::FeedEnthalpy) {
            ThermochemistryRecord record;
            std::vector<double> support;
            if (model.thermochemistry()) {
                record.standard_pressure = model.thermochemistry()->standard_pressure;
                for (std::size_t i = 0; i < observation.composition.size(); ++i) if (observation.composition[i] > 0.0) {
                    support.push_back(observation.composition[i]);
                    record.component_intervals.push_back(model.thermochemistry()->component_intervals[i]);
                }
                const auto ideal = evaluate_ideal_properties(record, observation.temperature, state.molar_density, support);
                caloric = ideal.value ? std::optional<double>(ideal.value->enthalpy + state.residual_enthalpy) : std::nullopt;
                thermal_status = ideal.status;
            } else { caloric.reset(); thermal_status = ThermochemistryStatus::MissingRecord; }
        }
        if (observation.kind != ObservationKind::Density && !caloric) {
            result.status = ObservationStatus::ThermochemistryUnavailable;
            result.message = thermal_status == ThermochemistryStatus::MissingRecord ? "missing_record"
                : thermal_status == ThermochemistryStatus::OutsideInterval ? "outside_interval" : "nonfinite";
            return result;
        }
        double value = observation.kind == ObservationKind::Density ? state.molar_density : *caloric;
        if (observation.mass_basis)
            value *= mean_molar_mass(parameters, observation.composition.data(), missing_mass);
        if (missing_mass) {
            result.status = ObservationStatus::MissingMolarMass;
            return result;
        }
        finish(result, observation, value);
        return result;
    }
    if (observation.kind == ObservationKind::LogEquilibriumConstant) {
        result.unit = "1";
        for (const auto& r : observation.problem.reactions) if (r.name == observation.reaction)
            finish(result, observation, r.correlation ? evaluate_reaction_properties(*r.correlation, observation.problem.temperature.value).ln_k : std::log(*r.equilibrium_constant));
        return result;
    }
    equilibrium::EquilibriumResult solved;
    try {
        if (!result.compiled) result.compiled = equilibrium::compile_problem(model, observation.problem);
        solved = equilibrium::solve_compiled(*result.compiled, start);
    } catch (const core::ConvergenceFailure& failure) {
        result.status = ObservationStatus::CentralSolveUnavailable;
        result.convergence = failure.record;
        result.message = failure.what();
        return result;
    } catch (const std::runtime_error& failure) {
        result.status = ObservationStatus::CentralSolveUnavailable;
        result.message = failure.what();
        return result;
    }
    result.convergence = solved.convergence;
    result.requested_tolerance_met = solved.requested_tolerance_met;
    result.message = solved.solver_status + ": " + solved.message;
    // The driver sets the topology fields before its acceptance gate, so a vanished or
    // coalesced phase names the cause whatever `success` says.
    if (solved.topology_event != equilibrium::TopologyEvent::None) {
        result.status = ObservationStatus::TopologyEvent;
        return result;
    }
    if (!solved.success || !solved.has_state()) {
        result.status = ObservationStatus::CentralSolveUnavailable;
        return result;
    }
    result.solved = std::move(solved);
    project_prediction(model, observation, result);
    return result;
}

void project_prediction(const HelmholtzModel& model, const Observation& observation, Prediction& result) {
    if (electrolyte_observation(observation.kind)) {
        double value = 0.0;
        result.unit = "1";
        switch (observation.kind) {
        case ObservationKind::BulkRelativePermittivity:
            value = result.state->bulk_relative_permittivity;
            if (result.reference_state) value /= result.reference_state->bulk_relative_permittivity;
            break;
        case ObservationKind::MeanIonicActivityCoefficient: value = result.activity ? result.activity->mean_ionic_molality_scale : 1.0; break;
        case ObservationKind::OsmoticCoefficient:
            if (!result.activity->osmotic_coefficient || !std::holds_alternative<double>(*result.activity->osmotic_coefficient)) {
                result.status = ObservationStatus::ActionUnavailable; result.message = "osmotic precision floor or missing pure solvent"; return;
            }
            value = std::get<double>(*result.activity->osmotic_coefficient); break;
        default:
            result.unit = "J/mol";
            const auto ion = observation.species[0] == observation.ion_indices[0] ? 0 : 1;
            value = GAS_CONSTANT_J_PER_MOL_K * observation.temperature * (result.reference->log_fugacity[ion]
                - (result.second_reference ? result.second_reference->log_fugacity[ion] : 0.0));
        }
        result.residual.reset();
        finish(result, observation, value);
        return;
    }
    const auto& parameters = model.parameters();
    const auto components = parameters.component_count();
    if (pressure_closed(observation.kind) || observation.kind == ObservationKind::LogEquilibriumConstant
        || observation.kind == ObservationKind::AdditionEnthalpy) return;
    const auto& solved = *result.solved;
    bool missing_mass = false;
    result.value.reset();
    result.residual.reset();
    const auto phase = observation.phase;
    const double* fractions = solved.mole_fractions.data() + phase * components;
    double value = 0.0;
    switch (observation.kind) {
    case ObservationKind::SystemEnthalpy: {
        equilibrium::SolvedStateActionDirection zero{0.0, 0.0, std::vector<double>(components), std::vector<double>(model.active_parameters().size())};
        const auto entries = equilibrium::solved_state_actions(*result.compiled, solved,
            {enthalpy_observables(result.compiled->problem), {zero}}).results;
        for (std::size_t i = 0; i < entries.size(); i += 2) {
            for (const auto* e : {&entries[i], &entries[i + 1]}) if (!e->value) {
                result.status = e->status == equilibrium::SolvedStateActionStatus::MissingThermochemistry
                    || e->status == equilibrium::SolvedStateActionStatus::OutsideInterval
                    || e->status == equilibrium::SolvedStateActionStatus::ReferenceUnavailable
                    ? ObservationStatus::ThermochemistryUnavailable : ObservationStatus::ActionUnavailable;
                result.message = std::string(equilibrium::solved_state_status_name(e->status)) + ": " + e->diagnostic_message.value_or("");
                return;
            }
            value += *entries[i].value * *entries[i + 1].value;
        }
        result.unit = "J";
        break;
    }
    case ObservationKind::EquilibriumPressure:
        result.unit = "Pa";
        value = solved.pressure.front();
        break;
    case ObservationKind::EquilibriumTemperature:
        result.unit = "K";
        value = solved.temperature.front();
        break;
    case ObservationKind::PhaseMoleFraction:
    case ObservationKind::PartialPressure:
        for (const auto species : observation.species) value += fractions[species];
        result.unit = observation.kind == ObservationKind::PartialPressure ? "Pa" : "1";
        if (observation.kind == ObservationKind::PartialPressure) value *= solved.pressure.front();
        break;
    case ObservationKind::PhaseAmountRatio:
        result.unit = "1";
        value = component_amount(solved, phase, observation.species)
              / component_amount(solved, observation.denominator_phase, observation.denominator_species);
        break;
    case ObservationKind::PhaseFraction: {
        result.unit = "1";
        double total = 0.0;
        for (const auto amount : solved.amounts) total += amount;
        value = solved.amounts[phase] / total;
        break;
    }
    case ObservationKind::PhaseDensity:
        result.unit = observation.mass_basis ? "kg/m^3" : "mol/m^3";
        value = solved.molar_densities[phase];
        if (observation.mass_basis) value *= mean_molar_mass(parameters, fractions, missing_mass);
        break;
    case ObservationKind::AdditionEnthalpy:
    case ObservationKind::FeedEnthalpy:
    case ObservationKind::IsobaricHeatCapacity:
    case ObservationKind::LogEquilibriumConstant:
    case ObservationKind::Density:
        break;
    }
    if (missing_mass) {
        result.status = ObservationStatus::MissingMolarMass;
        return;
    }
    finish(result, observation, value);
}

PredictionGradient differentiate(
    const HelmholtzModel& model, const Observation& observation, const Prediction& prediction
) {
    if (prediction.status != ObservationStatus::Available)
        return {{}, prediction.status, prediction.message};
    if (electrolyte_observation(observation.kind)) {
        PredictionGradient gradient;
        const auto width = model.active_parameters().size();
        for (std::size_t a = 0; a < width; ++a) {
            double value = 0.0;
            if (observation.kind == ObservationKind::BulkRelativePermittivity) {
                value = prediction.state->permittivity_parameter_derivative.at(a);
                if (prediction.reference_state) {
                    const auto& ref = *prediction.reference_state;
                    value = (value - *prediction.value * ref.permittivity_parameter_derivative.at(a)) / ref.bulk_relative_permittivity;
                }
            } else if (observation.kind == ObservationKind::MeanIonicActivityCoefficient)
                value = prediction.activity ? prediction.activity->mean_ionic_parameter_derivative.at(a) : 0.0;
            else if (observation.kind == ObservationKind::OsmoticCoefficient) value = prediction.activity->osmotic_parameter_derivative.at(a);
            else {
                const auto ion = observation.species[0] == observation.ion_indices[0] ? 0 : 1;
                const auto column = prediction.reference->parameter_derivative.at(a)[ion];
                const auto ref = prediction.second_reference ? prediction.second_reference->parameter_derivative.at(a)[ion] : std::optional<double>(0.0);
                value = column && ref ? GAS_CONSTANT_J_PER_MOL_K * observation.temperature * (*column - *ref) : NAN;
            }
            if (!std::isfinite(value)) return {{}, ObservationStatus::ActionUnavailable, "electrolyte column: reference admission or pressure slope"};
            gradient.value.push_back(value);
        }
        gradient.status = ObservationStatus::Available;
        return gradient;
    }
    if (!pressure_closed(observation.kind))
        return differentiate_group(model, {&observation}, prediction).front();
    PredictionGradient gradient;
    const auto& active = model.active_parameters();
    const auto& parameters = model.parameters();
    const auto components = parameters.component_count();
    const std::size_t width = active.size();
    gradient.value.assign(width, 0.0);
    const auto unavailable = [&](const char* what) {
        gradient.value.clear();
        gradient.status = ObservationStatus::ActionUnavailable;
        gradient.message = what;
        return gradient;
    };
    bool missing_mass = false;
    const double factor = observation.mass_basis
        ? mean_molar_mass(parameters, observation.composition.data(), missing_mass) : 1.0;
    const bool caloric = observation.kind != ObservationKind::Density;
    const std::size_t count = caloric ? 1 : 3;
    const std::vector<PropertyObservable> observables = caloric
        ? std::vector<PropertyObservable>{observation.kind == ObservationKind::FeedEnthalpy
            ? PropertyObservable::ResidualEnthalpy : PropertyObservable::TotalIsobaricHeatCapacity} : std::vector<PropertyObservable>{};
    for (std::size_t first = 0; first < width; first += count) {
        std::vector<PressureActionDirection> directions;
        for (std::size_t k = first; k < width && k < first + count; ++k)
            directions.push_back({0.0, 0.0, std::vector<double>(components, 0.0),
                                  {{active[k].family, active[k].indices, 1.0}}});
        const auto batch = model.pressure_actions(observation.temperature, observation.pressure,
            observation.composition, observation.branch, observables, directions);
        if (batch.density_status != PropertyActionStatus::Available)
            return unavailable(equilibrium::solved_state_status_name(equilibrium::solved_state_status(batch.density_status)));
        if (caloric && batch.results[0].status != PropertyActionStatus::Available)
            return unavailable(equilibrium::solved_state_status_name(equilibrium::solved_state_status(batch.results[0].status)));
        for (std::size_t k = 0; k < directions.size(); ++k)
            gradient.value[first + k] = (caloric ? *batch.results[0].action : batch.density_action[k]) * factor;
    }
    gradient.status = ObservationStatus::Available;
    return gradient;
}

std::vector<PredictionGradient> differentiate_group(
    const HelmholtzModel& model, const std::vector<const Observation*>& observations, const Prediction& prediction, const CoordinateSet& coordinates
) {
    if (electrolyte_observation(observations.front()->kind)) {
        std::vector<PredictionGradient> result;
        for (const auto* row : observations) {
            auto projected = prediction;
            project_prediction(model, *row, projected);
            result.push_back(differentiate(model, *row, projected));
        }
        return result;
    }
    const auto& parameters = model.parameters();
    const auto components = parameters.component_count();
    const std::size_t active_width = model.active_parameters().size();
    std::vector<ReactionCoordinate> reactions;
    for (const auto& c : coordinates.coordinates) if (c.reaction) reactions.push_back(*c.reaction);
    const std::size_t width = active_width + reactions.size();
    const auto& problem = prediction.compiled ? prediction.compiled->problem : observations.front()->problem;
    if (observations.front()->kind == ObservationKind::AdditionEnthalpy) {
        auto item = *observations.front();
        item.kind = ObservationKind::SystemEnthalpy;
        auto initial = differentiate_group(model, {&item}, prediction, coordinates).front();
        if (initial.status != ObservationStatus::Available) { initial.message = "initial:" + initial.message; return {initial}; }
        auto center = prediction;
        std::swap(center.solved, center.final_solved); std::swap(center.compiled, center.final_compiled);
        item.problem = item.final_problem;
        auto final = differentiate_group(model, {&item}, center, coordinates).front();
        if (final.status != ObservationStatus::Available) { final.message = "final:" + final.message; return {final}; }
        item.kind = ObservationKind::FeedEnthalpy;
        const auto feed = differentiate(model, item, prediction);
        if (feed.status != ObservationStatus::Available) return {{{}, feed.status, "feed:" + feed.message}};
        for (std::size_t k = 0; k < width; ++k)
            initial.value[k] = (final.value[k] - initial.value[k]) / item.feed_amount - (k < feed.value.size() ? feed.value[k] : 0.0);
        return {initial};
    }
    // One unit direction per active parameter at one shared center (one factorization, one reference ladder).
    std::vector<equilibrium::SolvedStateActionDirection> directions(width);
    for (std::size_t k = 0; k < width; ++k) {
        directions[k].feed_amounts_mol.assign(components, 0.0);
        directions[k].active_parameters.assign(active_width, 0.0);
        if (k < active_width) directions[k].active_parameters[k] = 1.0;
        else {
            directions[k].reaction_log_constants.assign(problem.reactions.size(), 0.0);
            directions[k].reaction_log_constant_temperature_slopes.assign(problem.reactions.size(), 0.0);
            const auto& target = reactions[k - active_width];
            for (std::size_t j = 0; j < problem.reactions.size(); ++j) {
                const auto& reaction = problem.reactions[j];
                if (reaction.name != target.name) continue;
                const double T = prediction.solved ? prediction.solved->temperature.front() : problem.temperature.value;
                const double basis[] = {1.0, 1.0 / T,
                    reaction.correlation ? std::log(T / reaction.correlation->reference_temperature) : 0.0, T};
                directions[k].reaction_log_constants[j] = basis[static_cast<std::size_t>(target.coefficient)];
                const double slope[] = {0.0, -1.0 / (T * T), 1.0 / T, 1.0};
                directions[k].reaction_log_constant_temperature_slopes[j] = slope[static_cast<std::size_t>(target.coefficient)];
            }
        }
    }
    const auto* solved = prediction.solved ? &*prediction.solved : nullptr;
    std::vector<std::vector<std::size_t>> slots(observations.size());
    std::vector<SolvedStateObservable> observables;
    for (std::size_t row = 0; row < observations.size(); ++row) {
        const auto& observation = *observations[row];
        const auto phase = observation.phase;
        const auto add = [&](const SolvedStateObservable& item) {
            const auto found = std::find_if(observables.begin(), observables.end(), [&](const auto& prior) {
                return std::tie(prior.kind, prior.phase, prior.component, prior.property)
                    == std::tie(item.kind, item.phase, item.component, item.property);
            });
            slots[row].push_back(static_cast<std::size_t>(found - observables.begin()));
            if (found == observables.end()) observables.push_back(item);
        };
        switch (observation.kind) {
        case ObservationKind::SystemEnthalpy:
            for (const auto& entry : enthalpy_observables(problem)) add(entry);
            break;
        case ObservationKind::EquilibriumPressure: add(observable(Kind::Pressure, 0)); break;
        case ObservationKind::EquilibriumTemperature: add(observable(Kind::Temperature, 0)); break;
        case ObservationKind::PhaseMoleFraction:
        case ObservationKind::PartialPressure:
            for (const auto species : observation.species)
                add(observable(Kind::PhaseMoleFraction, phase, species));
            if (observation.kind == ObservationKind::PartialPressure) add(observable(Kind::Pressure, 0));
            break;
        case ObservationKind::PhaseAmountRatio:
            for (const auto i : observation.species) add(observable(Kind::PhaseComponentAmount, phase, i));
            for (const auto i : observation.denominator_species)
                add(observable(Kind::PhaseComponentAmount, observation.denominator_phase, i));
            break;
        case ObservationKind::PhaseFraction:
            for (std::size_t q = 0; q < solved->phase_count; ++q) add(observable(Kind::PhaseAmount, q));
            break;
        case ObservationKind::PhaseDensity:
            add(observable(Kind::PhaseMolarDensity, phase));
            if (observation.mass_basis)
                for (std::size_t i = 0; i < components; ++i)
                    add(observable(Kind::PhaseMoleFraction, phase, i));
            break;
        case ObservationKind::AdditionEnthalpy:
        case ObservationKind::FeedEnthalpy:
        case ObservationKind::IsobaricHeatCapacity:
        case ObservationKind::LogEquilibriumConstant:
        case ObservationKind::Density: break;
        }
    }
    const auto batches = !solved || width == 0 ? std::vector<equilibrium::SolvedStateActionBatch>()
        : equilibrium::solved_state_first_actions(*prediction.compiled, *solved, observables, directions);
    std::vector<PredictionGradient> gradients(observations.size());
    for (std::size_t row = 0; row < observations.size(); ++row) {
        const auto& observation = *observations[row];
        const auto phase = observation.phase;
        const double* fractions = solved ? solved->mole_fractions.data() + phase * components : nullptr;
        bool missing_mass = false;
        const double mass = observation.mass_basis ? mean_molar_mass(parameters, fractions, missing_mass) : 1.0;
        auto& gradient = gradients[row];
        gradient.value.assign(width, 0.0);
        gradient.status = ObservationStatus::Available;
        for (std::size_t k = 0; k < width; ++k) {
            if (observation.kind == ObservationKind::LogEquilibriumConstant) {
                for (std::size_t j = 0; k >= active_width && j < problem.reactions.size(); ++j)
                    if (problem.reactions[j].name == observation.reaction) gradient.value[k] = directions[k].reaction_log_constants[j];
                continue;
            }
            std::vector<double> action;
            for (const auto j : slots[row]) {
                const auto& item = batches[k].results[j];
                if (item.status != equilibrium::SolvedStateActionStatus::Available || !item.action.has_value()) {
                    gradient.value.clear();
                    gradient.status = ObservationStatus::ActionUnavailable;
                    gradient.message = std::string(equilibrium::solved_state_status_name(item.status)) + ": " + item.diagnostic_message.value_or("");
                    break;
                }
                action.push_back(*item.action);
            }
            if (gradient.status != ObservationStatus::Available) break;
            double value = 0.0;
            switch (observation.kind) {
            case ObservationKind::SystemEnthalpy:
                for (std::size_t i = 0; i < action.size(); i += 2)
                    value += action[i] * *batches[k].results[slots[row][i + 1]].value
                        + *batches[k].results[slots[row][i]].value * action[i + 1];
                break;
            case ObservationKind::PhaseAmountRatio: {
                const double A = component_amount(*solved, phase, observation.species);
                const double B = component_amount(*solved, observation.denominator_phase, observation.denominator_species);
                double dA = 0.0, dB = 0.0;
                for (std::size_t i = 0; i < action.size(); ++i)
                    (i < observation.species.size() ? dA : dB) += action[i];
                value = (dA * B - A * dB) / (B * B);
                break;
            }
            case ObservationKind::PhaseFraction: {
                double total = 0.0, total_action = 0.0;
                for (std::size_t q = 0; q < solved->phase_count; ++q) {
                    total += solved->amounts[q];
                    total_action += action[q];
                }
                value = (action[phase] * total - solved->amounts[phase] * total_action) / (total * total);
                break;
            }
            case ObservationKind::PhaseDensity:
                value = action[0] * mass;
                for (std::size_t i = 0; observation.mass_basis && i < components; ++i)
                    value += solved->molar_densities[phase] * parameters.molar_mass_kg_per_mol[i] * action[1 + i];
                break;
            case ObservationKind::PartialPressure:  // d(P sum y) = dP sum y + P sum dy
                for (const auto species : observation.species) value += fractions[species] * action.back();
                for (std::size_t i = 0; i + 1 < action.size(); ++i) value += solved->pressure.front() * action[i];
                break;
            default:
                for (const auto item : action) value += item;
            }
            gradient.value[k] = value;
        }
    }
    return gradients;
}

}  // namespace epcsaft::regression
