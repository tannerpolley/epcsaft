#include "regression/observations.hpp"

#include "equilibrium/gibbs_minimizer.hpp"

#include <cmath>
#include <set>
#include <stdexcept>
#include <string>

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

double mean_molar_mass(const Parameters& parameters, const double* composition, bool& missing) {
    double result = 0.0;
    for (std::size_t i = 0; i < parameters.component_count(); ++i) {
        const double mass = parameters.molar_mass_kg_per_mol[i];
        if (composition[i] > 0.0 && !(mass > 0.0)) missing = true;
        result += composition[i] * mass;
    }
    return result;
}

void finish(Prediction& prediction, const Observation& observation, double value) {
    if (!std::isfinite(value)) {
        prediction.status = ObservationStatus::Nonfinite;
        return;
    }
    prediction.value = value;
    prediction.status = ObservationStatus::Available;
    prediction.residual = residual(observation, value, prediction.status);
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
    if (observation.form == ResidualForm::Difference) return (value - observation.observed) / observation.scale;
    if (!(value > 0.0) || !(observation.observed > 0.0)) {
        status = ObservationStatus::NonpositiveLogArgument;
        return std::nullopt;
    }
    return std::log(value / observation.observed);
}

double residual_slope(const Observation& observation, double value) {
    return observation.form == ResidualForm::Difference ? 1.0 / observation.scale : 1.0 / value;
}

void validate(const Observation& observation, const HelmholtzModel& model) {
    const auto components = model.parameters().component_count();
    require(std::isfinite(observation.observed), "observed value");
    require(observation.form == ResidualForm::LogRatio
            || (std::isfinite(observation.scale) && observation.scale > 0.0), "scale");
    if (observation.kind == ObservationKind::Density) {
        require(std::isfinite(observation.temperature) && observation.temperature > 0.0, "temperature");
        require(std::isfinite(observation.pressure) && observation.pressure > 0.0, "pressure");
        require(observation.composition.size() == components, "composition");
        return;
    }
    const auto compiled = equilibrium::compile_problem(model, observation.problem);
    require(observation.phase < compiled.problem.phases.size(), "phase index");
    if (observation.kind != ObservationKind::PhaseMoleFraction) return;
    require(!observation.species.empty(), "species set is empty");
    const auto& support = compiled.problem.phases[observation.phase].support;
    std::set<std::size_t> seen;
    for (const auto species : observation.species) {
        require(seen.insert(species).second, "species repeated");
        require(std::find(support.begin(), support.end(), species) != support.end(),
                "species outside the phase support");
    }
}

Prediction predict(
    const HelmholtzModel& model, const Observation& observation,
    const std::optional<equilibrium::EquilibriumResult>& start
) {
    Prediction result;
    const auto& parameters = model.parameters();
    const auto components = parameters.component_count();
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
        solved = equilibrium::solve_equilibrium(model, observation.problem, start);
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
        result.unit = "1";
        for (const auto species : observation.species) value += fractions[species];
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
    result.solved = std::move(solved);
    if (missing_mass) {
        result.status = ObservationStatus::MissingMolarMass;
        return result;
    }
    finish(result, observation, value);
    return result;
}

PredictionGradient differentiate(
    const HelmholtzModel& model, const Observation& observation, const Prediction& prediction
) {
    PredictionGradient gradient;
    if (prediction.status != ObservationStatus::Available) {
        gradient.status = prediction.status;
        gradient.message = prediction.message;
        return gradient;
    }
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
    if (observation.kind == ObservationKind::Density) {
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
    const auto& solved = *prediction.solved;
    const auto compiled = equilibrium::compile_problem(model, observation.problem);
    const auto phase = observation.phase;
    const double* fractions = solved.mole_fractions.data() + phase * components;
    std::vector<SolvedStateObservable> observables;
    switch (observation.kind) {
    case ObservationKind::EquilibriumPressure: observables.push_back(observable(Kind::Pressure, 0)); break;
    case ObservationKind::EquilibriumTemperature: observables.push_back(observable(Kind::Temperature, 0)); break;
    case ObservationKind::PhaseMoleFraction:
        for (const auto species : observation.species)
            observables.push_back(observable(Kind::PhaseMoleFraction, phase, species));
        break;
    case ObservationKind::PhaseFraction:
        for (std::size_t q = 0; q < solved.phase_count; ++q) observables.push_back(observable(Kind::PhaseAmount, q));
        break;
    case ObservationKind::PhaseDensity:
        observables.push_back(observable(Kind::PhaseMolarDensity, phase));
        if (observation.mass_basis)
            for (std::size_t i = 0; i < components; ++i)
                observables.push_back(observable(Kind::PhaseMoleFraction, phase, i));
        break;
    case ObservationKind::Density: break;
    }
    const double mass = observation.mass_basis ? mean_molar_mass(parameters, fractions, missing_mass) : 1.0;
    for (std::size_t k = 0; k < width; ++k) {
        equilibrium::SolvedStateActionDirection direction;
        direction.feed_amounts_mol.assign(components, 0.0);
        direction.active_parameters.assign(width, 0.0);
        direction.active_parameters[k] = 1.0;
        const auto batch = equilibrium::solved_state_actions(compiled, solved, {observables, {direction}});
        std::vector<double> action;
        for (const auto& item : batch.results) {
            if (item.status != equilibrium::SolvedStateActionStatus::Available || !item.action.has_value())
                return unavailable(equilibrium::solved_state_status_name(item.status));
            action.push_back(*item.action);
        }
        double value = 0.0;
        switch (observation.kind) {
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
        default:
            for (const auto item : action) value += item;
        }
        gradient.value[k] = value;
    }
    gradient.status = ObservationStatus::Available;
    return gradient;
}

}  // namespace epcsaft::regression
