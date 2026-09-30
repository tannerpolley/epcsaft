#include "regression/observations.hpp"

#include "core/numerical_policy.hpp"
#include <algorithm>
#include <cmath>
#include <set>
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
    case ObservationStatus::Available: return "available";
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

bool same_reaction(const equilibrium::ReactionSpec& x, const equilibrium::ReactionSpec& y) {
    const auto correlation = [](const auto& v) {
        return std::tie(v.a, v.b, v.c, v.d, v.reference_temperature, v.temperature_min,
                        v.temperature_max, v.standard_state_id);
    };
    const auto reference = [](const auto& v) { return std::tie(v.reference_pressure_pa, v.source_basis); };
    return std::tie(x.name, x.stoichiometry, x.equilibrium_constant, x.affinity_j_per_mol)
        == std::tie(y.name, y.stoichiometry, y.equilibrium_constant, y.affinity_j_per_mol)
        && optional_equal(x.correlation, y.correlation, correlation)
        && optional_equal(x.reference, y.reference, reference);
}

bool same_state(const Observation& left, const Observation& right) {
    // Density stays a singleton; equilibrium projection and residual metadata do not select physics.
    if (left.kind == ObservationKind::Density || right.kind == ObservationKind::Density) return false;
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
    if (observation.kind == ObservationKind::Density) {
        require(std::isfinite(observation.temperature) && observation.temperature > 0.0, "temperature");
        require(std::isfinite(observation.pressure) && observation.pressure > 0.0, "pressure");
        require(observation.composition.size() == components, "composition");
        return;
    }
    const auto owned = prepared ? std::optional<equilibrium::CompiledSystem>()
                                : equilibrium::compile_problem(model, observation.problem);
    const auto& compiled = prepared ? *prepared : *owned;
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
    std::optional<equilibrium::CompiledSystem> prepared
) {
    Prediction result;
    result.compiled = std::move(prepared);
    const auto& parameters = model.parameters();
    bool missing_mass = false;
    if (observation.kind == ObservationKind::Density) {
        result.unit = observation.mass_basis ? "kg/m^3" : "mol/m^3";
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
        double value = state.molar_density;
        if (observation.mass_basis)
            value *= mean_molar_mass(parameters, observation.composition.data(), missing_mass);
        if (missing_mass) {
            result.status = ObservationStatus::MissingMolarMass;
            return result;
        }
        finish(result, observation, value);
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
    const auto& parameters = model.parameters();
    const auto components = parameters.component_count();
    const auto& solved = *result.solved;
    bool missing_mass = false;
    result.value.reset();
    result.residual.reset();
    const auto phase = observation.phase;
    const double* fractions = solved.mole_fractions.data() + phase * components;
    double value = 0.0;
    switch (observation.kind) {
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
    if (observation.kind != ObservationKind::Density)
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
    for (std::size_t first = 0; first < width; first += 3) {
        std::vector<PressureActionDirection> directions;
        for (std::size_t k = first; k < width && k < first + 3; ++k)
            directions.push_back({0.0, 0.0, std::vector<double>(components, 0.0),
                                  {{active[k].family, active[k].indices, 1.0}}});
        const auto batch = model.pressure_actions(observation.temperature, observation.pressure,
            observation.composition, observation.branch, {}, directions);
        if (batch.density_status != PropertyActionStatus::Available) return unavailable("density action");
        for (std::size_t k = 0; k < directions.size(); ++k)
            gradient.value[first + k] = batch.density_action[k] * factor;
    }
    gradient.status = ObservationStatus::Available;
    return gradient;
}

std::vector<PredictionGradient> differentiate_group(
    const HelmholtzModel& model, const std::vector<const Observation*>& observations, const Prediction& prediction, const CoordinateSet& coordinates
) {
    const auto& parameters = model.parameters();
    const auto components = parameters.component_count();
    const std::size_t active_width = model.active_parameters().size();
    std::vector<ReactionCoordinate> reactions;
    for (const auto& c : coordinates.coordinates) if (c.reaction) reactions.push_back(*c.reaction);
    const std::size_t width = active_width + reactions.size();
    const auto& solved = *prediction.solved;
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
            for (std::size_t q = 0; q < solved.phase_count; ++q) add(observable(Kind::PhaseAmount, q));
            break;
        case ObservationKind::PhaseDensity:
            add(observable(Kind::PhaseMolarDensity, phase));
            if (observation.mass_basis)
                for (std::size_t i = 0; i < components; ++i)
                    add(observable(Kind::PhaseMoleFraction, phase, i));
            break;
        case ObservationKind::Density: break;
        }
    }
    // One unit direction per active parameter at one shared center (one factorization, one reference ladder).
    std::vector<equilibrium::SolvedStateActionDirection> directions(width);
    for (std::size_t k = 0; k < width; ++k) {
        directions[k].feed_amounts_mol.assign(components, 0.0);
        directions[k].active_parameters.assign(active_width, 0.0);
        if (k < active_width) directions[k].active_parameters[k] = 1.0;
        else {
            directions[k].reaction_log_constants.assign(prediction.compiled->problem.reactions.size(), 0.0);
            const auto& target = reactions[k - active_width];
            for (std::size_t j = 0; j < prediction.compiled->problem.reactions.size(); ++j) {
                const auto& reaction = prediction.compiled->problem.reactions[j];
                if (reaction.name != target.name) continue;
                const double T = solved.temperature.front();
                const double basis[] = {1.0, 1.0 / T,
                    reaction.correlation ? std::log(T / reaction.correlation->reference_temperature) : 0.0, T};
                directions[k].reaction_log_constants[j] = basis[static_cast<std::size_t>(target.coefficient)];
            }
        }
    }
    const auto batches = width == 0 ? std::vector<equilibrium::SolvedStateActionBatch>()
        : equilibrium::solved_state_first_actions(*prediction.compiled, solved, observables, directions);
    std::vector<PredictionGradient> gradients(observations.size());
    for (std::size_t row = 0; row < observations.size(); ++row) {
        const auto& observation = *observations[row];
        const auto phase = observation.phase;
        const double* fractions = solved.mole_fractions.data() + phase * components;
        bool missing_mass = false;
        const double mass = observation.mass_basis ? mean_molar_mass(parameters, fractions, missing_mass) : 1.0;
        auto& gradient = gradients[row];
        gradient.value.assign(width, 0.0);
        gradient.status = ObservationStatus::Available;
        for (std::size_t k = 0; k < width; ++k) {
            std::vector<double> action;
            for (const auto j : slots[row]) {
                const auto& item = batches[k].results[j];
                if (item.status != equilibrium::SolvedStateActionStatus::Available || !item.action.has_value()) {
                    gradient.value.clear();
                    gradient.status = ObservationStatus::ActionUnavailable;
                    gradient.message = equilibrium::solved_state_status_name(item.status);
                    break;
                }
                action.push_back(*item.action);
            }
            if (gradient.status != ObservationStatus::Available) break;
            double value = 0.0;
            switch (observation.kind) {
            case ObservationKind::PhaseAmountRatio: {
                const double A = component_amount(solved, phase, observation.species);
                const double B = component_amount(solved, observation.denominator_phase, observation.denominator_species);
                double dA = 0.0, dB = 0.0;
                for (std::size_t i = 0; i < action.size(); ++i)
                    (i < observation.species.size() ? dA : dB) += action[i];
                value = (dA * B - A * dB) / (B * B);
                break;
            }
            case ObservationKind::PhaseFraction: {
                double total = 0.0, total_action = 0.0;
                for (std::size_t q = 0; q < solved.phase_count; ++q) {
                    total += solved.amounts[q];
                    total_action += action[q];
                }
                value = (action[phase] * total - solved.amounts[phase] * total_action) / (total * total);
                break;
            }
            case ObservationKind::PhaseDensity:
                value = action[0] * mass;
                for (std::size_t i = 0; observation.mass_basis && i < components; ++i)
                    value += solved.molar_densities[phase] * parameters.molar_mass_kg_per_mol[i] * action[1 + i];
                break;
            case ObservationKind::PartialPressure:  // d(P sum y) = dP sum y + P sum dy
                for (const auto species : observation.species) value += fractions[species] * action.back();
                for (std::size_t i = 0; i + 1 < action.size(); ++i) value += solved.pressure.front() * action[i];
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
