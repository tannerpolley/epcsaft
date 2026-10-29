#include "regression/observations.hpp"

#include "core/numerical_policy.hpp"
#include "core/trace.hpp"
#include "properties/thermochemistry.hpp"
#include <algorithm>
#include <cmath>
#include <set>
#include <numeric>
#include <new>
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
            && left.reference_composition == right.reference_composition
            && left.finite_solute_molalities == right.finite_solute_molalities;
    // Pressure closures stay singletons; equilibrium projection and residual metadata do not select physics.
    if (pressure_closed(left.kind) || pressure_closed(right.kind)
        || left.kind == ObservationKind::LogEquilibriumConstant || right.kind == ObservationKind::LogEquilibriumConstant
        || addition(left.kind) || addition(right.kind)) return false;
    const auto& a = left.problem;
    const auto& b = right.problem;
    const auto intensive = [](const auto& v) { return std::tie(v.free, v.value, v.absent); };
    const auto feed = [](const auto& v) { return std::tie(v.kind, v.values, v.solvent, v.solvent_mass_kg); };
    const auto phase = [](const auto& v) {
        return std::tie(v.name, v.kind, v.support, v.support_declared, v.amount.kind, v.amount.fraction,
                        v.composition_guess, v.packing_guess);
    };
    if (intensive(a.temperature) != intensive(b.temperature) || intensive(a.pressure) != intensive(b.pressure)
        || feed(a.feed) != feed(b.feed) || a.standard_gibbs_j_per_mol != b.standard_gibbs_j_per_mol
        || a.topology_declared != b.topology_declared || a.solver_options.value(policy::row) != b.solver_options.value(policy::row)
        || a.phases.size() != b.phases.size() || a.reactions.size() != b.reactions.size()
        || !(a.source_reference == b.source_reference)) return false;
    for (std::size_t i = 0; i < a.phases.size(); ++i)
        if (phase(a.phases[i]) != phase(b.phases[i])) return false;
    for (std::size_t i = 0; i < a.reactions.size(); ++i)
        if (!same_reaction(a.reactions[i], b.reactions[i])) return false;
    return true;
}

void validate(const Observation& observation, const HelmholtzModel& model,
              const equilibrium::CompiledSystem* prepared) {
    const auto components = model.parameters().component_count();
    require(!observation.problem.discovery && !observation.final_problem.discovery, "declared phases");
    require(observation.problem.solids.empty() && observation.final_problem.solids.empty(), "without solids");
    require(std::isfinite(observation.observed), "observed value");
    require(observation.form == ResidualForm::LogRatio
            || (std::isfinite(observation.scale) && observation.scale > 0.0), "scale");
    require(!observation.mass_basis || observation.kind == ObservationKind::Density
            || observation.kind == ObservationKind::PhaseDensity, "mass basis is defined only for a density");
    require(!observation.finite_solute_molalities || observation.kind == ObservationKind::MeanIonicActivityCoefficient,
            "full solute molalities require mean ionic activity");
    if (observation.finite_solute_molalities)
        (void)formula_composition(observation.composition, model.parameters(), observation.ion_indices,
            observation.stoichiometry, observation.molality, observation.finite_solute_molalities);
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
        } else if (observation.kind == ObservationKind::MixtureOsmoticCoefficient) {
            model.validate_electroneutrality(observation.composition);
            const auto& parameters = model.parameters();
            std::size_t solvents = 0, solvent = 0;
            double solute = 0.0;
            for (std::size_t i = 0; i < components; ++i) if (observation.composition[i] > 0.0) {
                if (parameters.charge_number[i] == 0) { ++solvents; solvent = i; } else solute += observation.composition[i];
            }
            require(solvents == 1 && solute / (observation.composition[solvent] * parameters.molar_mass_kg_per_mol[solvent])
                >= policy::OSMOTIC_MOLALITY_FLOOR.value, "osmotic domain: one neutral solvent, solute molality above the floor");
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
    if (observation.kind == ObservationKind::SaturationCappedAdditionEnthalpy) {
        require(observation.form == ResidualForm::Difference, "capped addition uses Difference");
        require(std::isfinite(observation.basis_amount) && observation.basis_amount > 0
            && std::isfinite(observation.loading_cap) && observation.loading_cap > 0, "positive basis and cap");
        require(observation.composition.size() == components
            && std::count(observation.composition.begin(), observation.composition.end(), 1.0) == 1
            && std::count(observation.composition.begin(), observation.composition.end(), 0.0) == components - 1, "pure added feed");
        const auto seed = equilibrium::compile_problem(model, observation.final_problem);
        const auto initial = equilibrium::compile_problem(model, observation.problem);
        require(!observation.problem.temperature.free && !observation.problem.pressure.free
            && observation.problem.temperature.value == observation.temperature && observation.problem.pressure.value == observation.pressure
            && !observation.final_problem.temperature.free && observation.final_problem.temperature.value == observation.temperature
            && observation.final_problem.pressure.free && seed.pinned_count == 1, "isobaric initial and bubble endpoint");
        auto fixed = observation; fixed.kind = ObservationKind::AdditionEnthalpy;
        const auto species = std::find(observation.composition.begin(), observation.composition.end(), 1.0) - observation.composition.begin();
        require(model.parameters().charge_number[species] == 0, "neutral added feed");
        fixed.feed_amount = seed.n_feed[species] - initial.n_feed[species];
        validate(fixed, model, &initial);
        return;
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
        require(compiled.problem.source_reference == final.problem.source_reference, "endpoint source reference");
        for (std::size_t i = 0; i < components; ++i)
            require(std::abs(final.n_feed[i] - compiled.n_feed[i] - observation.feed_amount * observation.composition[i])
                <= policy::ADDITION_MATERIAL_BALANCE.value * std::max({1.0, std::abs(final.n_feed[i]), std::abs(compiled.n_feed[i])}),
                ("material balance component " + std::to_string(i)).c_str());
        return;
    }
    if (observation.kind == ObservationKind::EquilibriumIsobaricHeatCapacity) {
        require(!compiled.problem.temperature.free && !compiled.problem.pressure.free
                && !compiled.problem.temperature.absent && !compiled.problem.pressure.absent,
                "apparent Cp requires driven T/P");
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
    std::optional<equilibrium::CompiledSystem> prepared, ReferenceCache* references,
    const Prediction* previous, const HelmholtzModel* previous_model, bool same_parameters,
    std::chrono::steady_clock::time_point deadline
) {
    core::Scope density_scope("fit_density_anchor");
    const auto prior = [&](const State* old, const auto& x, const auto& old_x, DensityBranch branch) -> std::optional<double> {
        if (!old || !previous_model || !(model.solver_options() == previous_model->solver_options()) || x != old_x || branch == DensityBranch::Unspecified
            || !properties::regular_density(old->molar_density, old->pressure_density_derivative, old->density_convergence)) return std::nullopt;
        if (!same_parameters) {
            const auto witness = old->negative_slope_density;
            const double w = policy::DENSITY_ANCHOR_LOG_WINDOW.value;
            const double lower = observation.pressure / (GAS_CONSTANT_J_PER_MOL_K * observation.temperature)
                * policy::DENSITY_LOWER_FRACTION.value;
            const double upper = model.molar_density_upper_bound(observation.temperature, x) * policy::DENSITY_UPPER_FRACTION.value;
            if (!witness || !std::isfinite(*witness) || !(*witness > lower && *witness < upper)
                || !(branch == DensityBranch::Liquid ? old->molar_density * std::exp(-w) > *witness
                                                     : old->molar_density * std::exp(w) < *witness)) return std::nullopt;
            try {
                const double slope = model.state(observation.temperature, *witness, x).pressure_density_derivative;
                density_scope.probe("witness", {*witness, slope});
                if (!std::isfinite(slope) || !(slope < 0.0)) return std::nullopt;
            } catch (const std::runtime_error& error) {
                density_scope.probe("witness_ineligible", {*witness}, 0, error.what());
                return std::nullopt;
            }
        }
        density_scope.probe("eligible", old->molar_density);
        return old->molar_density;
    };
    const auto accepted = [&](State& state, const State* old, std::optional<double> anchor) {
        if (!anchor || state.stable_root_count) return;  // An external anchor was retried cold.
        state.negative_slope_density = old->negative_slope_density;
        density_scope.probe("accepted", state.molar_density);
    };
    const auto close = [&](const auto& x, const State* old) {
        const auto anchor = prior(old, x, x, observation.branch);
        auto state = properties::close_density(model, observation.temperature, observation.pressure, x, observation.branch, anchor);
        accepted(state, old, anchor);
        return state;
    };
    if (electrolyte_observation(observation.kind)) {
        Prediction result;
        result.solver_options = observation.problem.solver_options.resolved(model.solver_options());
        try {
            if (observation.kind == ObservationKind::BulkRelativePermittivity) {
                result.state = close(observation.composition, previous && previous->state ? &*previous->state : nullptr);
                if (!observation.reference_composition.empty()) result.reference_state = close(observation.reference_composition,
                    previous && previous->reference_state ? &*previous->reference_state : nullptr);
            } else if (observation.kind == ObservationKind::MixtureOsmoticCoefficient) {
                std::vector<double> water(observation.composition.size(), 0.0);
                for (std::size_t i = 0; i < water.size(); ++i)
                    if (observation.composition[i] > 0.0 && model.parameters().charge_number[i] == 0) water[i] = 1.0;
                const State* pure = previous && previous->reference && previous->reference->pure_solvent_state ? &*previous->reference->pure_solvent_state : nullptr;
                const auto anchor = prior(pure, water, water, DensityBranch::Liquid);
                result.reference = solvent_reference(model, observation.temperature, observation.pressure, water, anchor);
                if (result.reference->pure_solvent_state) accepted(*result.reference->pure_solvent_state, pure, anchor);
                auto activity = solvent_activity(model, observation.temperature, observation.pressure, observation.composition, *result.reference,
                    close(observation.composition, previous && previous->activity ? &previous->activity->state : nullptr));
                if (!activity.solvent_activity) throw std::runtime_error(activity.cause);
                result.activity.emplace();
                result.activity->state = std::move(activity.state);
                result.activity->osmotic_coefficient = activity.osmotic_coefficient;
                for (const auto& column : activity.osmotic_parameter_derivative) result.activity->osmotic_parameter_derivative.push_back(column.value_or(NAN));
            } else {
                ReferenceCache local;
                auto& cache = references ? *references : local;
                const auto reference = [&](const auto& x, const InfiniteDilutionReference* old) {
                    auto key = observation;
                    key.finite_solute_molalities.reset();
                    key.composition = x; key.molality = 0.0; key.reference_composition.clear();
                    for (const auto& entry : cache) if (same_state(key, entry.first)) return entry.second;
                    const State* solvent = old && old->solvent_state ? &*old->solvent_state : nullptr;
                    const auto anchor = prior(solvent, x, x, DensityBranch::Liquid);
                    auto value = infinite_dilution_reference(model, key.temperature, key.pressure, x,
                        key.ion_indices, key.stoichiometry, anchor);
                    accepted(*value.solvent_state, solvent, anchor);
                    cache.emplace_back(key, value);
                    return value;
                };
                result.reference = reference(observation.composition, previous && previous->reference ? &*previous->reference : nullptr);
                if (observation.kind == ObservationKind::TransferGibbsEnergy) result.second_reference = reference(observation.reference_composition,
                    previous && previous->second_reference ? &*previous->second_reference : nullptr);
                if (observation.kind == ObservationKind::MeanIonicActivityCoefficient || observation.kind == ObservationKind::OsmoticCoefficient) {
                    if (observation.molality > 0.0) {
                        const auto x = formula_composition(observation.composition, model.parameters(), observation.ion_indices, observation.stoichiometry, observation.molality, observation.finite_solute_molalities);
                        const auto old_x = previous_model ? formula_composition(observation.composition, previous_model->parameters(),
                            observation.ion_indices, observation.stoichiometry, observation.molality, observation.finite_solute_molalities) : x;
                        const State* old = previous && previous->activity ? &previous->activity->state : nullptr;
                        const auto anchor = prior(old, x, old_x, DensityBranch::Liquid);
                        result.activity = electrolyte_activity(model, observation.temperature, observation.pressure,
                            observation.composition, observation.ion_indices, observation.stoichiometry, observation.molality, *result.reference, anchor, observation.finite_solute_molalities);
                        accepted(result.activity->state, old, anchor);
                    }
                }
            }
            project_prediction(model, observation, result);
        } catch (const std::bad_alloc&) {
            throw;
        } catch (const std::exception& error) {
            result.status = observation.kind == ObservationKind::BulkRelativePermittivity
                ? ObservationStatus::DensityUnavailable : ObservationStatus::ReferenceUnavailable;
            retain_density_failure(result.first_numerical_failure, error, "prediction");
            result.message = error.what();
        }
        return result;
    }
    if (observation.kind == ObservationKind::SaturationCappedAdditionEnthalpy) {
        validate(observation, model, prepared ? &*prepared : nullptr);
        const auto initial = equilibrium::compile_problem(model, observation.problem);
        equilibrium::AffineLoading affine{initial.n_feed, observation.composition,
            observation.loading_bounds[0], observation.loading_bounds[1], 0, observation.pressure};
        for (auto& value : affine.direction) value *= observation.basis_amount;
        auto search = std::make_shared<equilibrium::SaturationSearch>(equilibrium::search_loading(
            model, observation.final_problem, affine, observation.loading_grid, nullptr, deadline));
        Prediction seed; seed.saturation = search; seed.first_numerical_failure = search->first_numerical_failure;
        if (!search->accepted) {
            seed.status = ObservationStatus::CentralSolveUnavailable;
            for (const auto& cause : search->causes) seed.message += "saturation:" + cause + ";";
            return seed;
        }
        const auto& state = search->candidates[*search->accepted].state;
        const auto compiled = equilibrium::compile_problem(model, *state.declaration);
        seed.saturation_loading = state.coordinates[compiled.loading_column()];
        seed.endpoint_loading = std::min(observation.loading_cap, *seed.saturation_loading);
        seed.added_amount = observation.basis_amount * *seed.endpoint_loading;
        seed.cap_branch = observation.loading_cap == *seed.saturation_loading ? "switch"
            : observation.loading_cap < *seed.saturation_loading ? "source_cap" : "saturation_cap";
        auto fixed = observation; fixed.kind = ObservationKind::AdditionEnthalpy; fixed.feed_amount = *seed.added_amount;
        fixed.final_problem = observation.problem; fixed.final_problem.feed.kind = equilibrium::FeedKind::Amounts;
        fixed.final_problem.feed.values = initial.n_feed;
        for (std::size_t i = 0; i < observation.composition.size(); ++i) fixed.final_problem.feed.values[i] += fixed.feed_amount * fixed.composition[i];
        auto result = predict(model, fixed, start, prepared, references, &seed);
        result.saturation = search; result.saturation_loading = seed.saturation_loading;
        result.endpoint_loading = seed.endpoint_loading; result.added_amount = seed.added_amount; result.cap_branch = seed.cap_branch;
        core::retain_first(result.first_numerical_failure, seed.first_numerical_failure, true);
        return result;
    }
    if (observation.kind == ObservationKind::AdditionEnthalpy) {
        auto item = observation;
        item.kind = ObservationKind::SystemEnthalpy;
        item.form = ResidualForm::Difference; item.scale = 1.0; item.observed = 0.0;
        auto initial = predict(model, item, start, std::move(prepared));
        const auto failed = [&](Prediction& p, const char* source) {
            core::retain_first(p.first_numerical_failure, initial.first_numerical_failure, true);
            core::retain_first(initial.first_numerical_failure, p.first_numerical_failure);
            if (p.status == ObservationStatus::Available) return false;
            p.message = std::string(source) + ":" + observation_status_name(p.status) + ": " + p.message;
            return true;
        };
        if (failed(initial, "initial")) return initial;
        item.problem = observation.final_problem;
        auto compiled = equilibrium::compile_problem(model, item.problem);
        const auto& a = initial.compiled->problem.phases;
        const auto& b = compiled.problem.phases;
        const bool compatible = a.size() == b.size() && std::equal(a.begin(), a.end(), b.begin(),
            [](const auto& x, const auto& y) { return x.kind == y.kind && x.support == y.support; });
        Prediction final;
        if (previous && previous->saturation && previous->cap_branch != "source_cap") {
            final.solved = previous->saturation->candidates[*previous->saturation->accepted].state;
            final.compiled = equilibrium::compile_problem(model, *final.solved->declaration);
            project_prediction(model, item, final);
        } else final = predict(model, item, compatible ? initial.solved : std::nullopt, std::move(compiled));
        if (failed(final, "final")) return final;
        item.kind = ObservationKind::FeedEnthalpy;
        auto feed = predict(model, item);
        if (failed(feed, "feed")) return feed;
        initial.enthalpy_difference = *final.value - *initial.value;
        finish(initial, observation, *initial.enthalpy_difference / observation.feed_amount - *feed.value);
        initial.unit = "J/mol";
        initial.final_solved = std::move(final.solved); initial.final_compiled = std::move(final.compiled);
        return initial;
    }
    Prediction result;
    result.solver_options = observation.problem.solver_options.resolved(model.solver_options());
    result.compiled = std::move(prepared);
    const auto& parameters = model.parameters();
    bool missing_mass = false;
    if (pressure_closed(observation.kind)) {
        result.unit = observation.kind == ObservationKind::FeedEnthalpy ? "J/mol"
            : observation.kind == ObservationKind::IsobaricHeatCapacity ? "J/(mol K)"
            : observation.mass_basis ? "kg/m^3" : "mol/m^3";
        State state;
        try {
            state = close(observation.composition, previous && previous->state ? &*previous->state : nullptr);
        } catch (const std::runtime_error& failure) {
            result.status = ObservationStatus::DensityUnavailable;
            if (const auto* nonconverged = dynamic_cast<const core::ConvergenceFailure*>(&failure)) result.convergence = nonconverged->record;
            retain_density_failure(result.first_numerical_failure, failure, "prediction density");
            result.message = failure.what();
            return result;
        }
        result.convergence = state.density_convergence;
        core::retain_first(result.first_numerical_failure, state.first_numerical_failure);
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
        result.state = std::move(state);
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
    } catch (const std::runtime_error& failure) {
        result.status = ObservationStatus::CentralSolveUnavailable;
        if (const auto* nonconverged = dynamic_cast<const core::ConvergenceFailure*>(&failure)) result.convergence = nonconverged->record;
        retain_density_failure(result.first_numerical_failure, failure, "prediction solve");
        result.message = failure.what();
        return result;
    }
    result.convergence = solved.convergence;
    core::retain_first(result.first_numerical_failure, solved.first_numerical_failure);
    if (!solved.success) core::retain_first(result.first_numerical_failure, solved.convergence);
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
        case ObservationKind::MixtureOsmoticCoefficient:
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
        || addition(observation.kind)) return;
    const auto& solved = *result.solved;
    bool missing_mass = false;
    result.value.reset();
    result.residual.reset();
    const auto phase = observation.phase;
    const double* fractions = solved.mole_fractions.data() + phase * components;
    double value = 0.0;
    switch (observation.kind) {
    case ObservationKind::EquilibriumIsobaricHeatCapacity:
    case ObservationKind::SystemEnthalpy: {
        equilibrium::SolvedStateActionDirection zero{0.0, 0.0, std::vector<double>(components), std::vector<double>(model.active_parameters().size())};
        const auto kind = observation.kind == ObservationKind::SystemEnthalpy
            ? Kind::SystemEnthalpy : Kind::EquilibriumIsobaricHeatCapacity;
        const auto entry = equilibrium::solved_state_actions(*result.compiled, solved,
            {{observable(kind, 0)}, {zero}}).results.front();
        if (!entry.value) {
            result.status = entry.status == equilibrium::SolvedStateActionStatus::MissingThermochemistry
                || entry.status == equilibrium::SolvedStateActionStatus::OutsideInterval
                || entry.status == equilibrium::SolvedStateActionStatus::ReferenceUnavailable
                ? ObservationStatus::ThermochemistryUnavailable : ObservationStatus::ActionUnavailable;
            result.message = std::string(equilibrium::solved_state_status_name(entry.status)) + ": " + entry.diagnostic_message.value_or("");
            return;
        }
        value = *entry.value;
        result.unit = entry.unit;
        break;
    }
    case ObservationKind::EquilibriumPressure:
        result.unit = "Pa";
        value = *solved.pressure;
        break;
    case ObservationKind::EquilibriumTemperature:
        result.unit = "K";
        value = *solved.temperature;
        break;
    case ObservationKind::PhaseMoleFraction:
    case ObservationKind::PartialPressure:
        for (const auto species : observation.species) value += fractions[species];
        result.unit = observation.kind == ObservationKind::PartialPressure ? "Pa" : "1";
        if (observation.kind == ObservationKind::PartialPressure) value *= *solved.pressure;
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
    case ObservationKind::SaturationCappedAdditionEnthalpy:
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
        return {{}, prediction.status, prediction.message, {}, {}, {}, prediction.first_numerical_failure};
    if (electrolyte_observation(observation.kind)) {
        PredictionGradient gradient;
        const auto width = model.active_parameters().size();
        for (std::size_t a = 0; a < width; ++a) {
            double value = 0.0;
            ResidualEnergyAction action;
            if (observation.kind == ObservationKind::BulkRelativePermittivity) {
                value = prediction.state->permittivity_parameter_derivative.at(a);
                if (prediction.reference_state) {
                    const auto& ref = *prediction.reference_state;
                    value = (value - *prediction.value * ref.permittivity_parameter_derivative.at(a)) / ref.bulk_relative_permittivity;
                }
            } else if (observation.kind == ObservationKind::MeanIonicActivityCoefficient) {
                if (prediction.activity) action = prediction.activity->mean_ionic_parameter_derivative.at(a);
                value = prediction.activity ? action.action.value_or(NAN) : 0.0;
            } else if (observation.kind == ObservationKind::OsmoticCoefficient || observation.kind == ObservationKind::MixtureOsmoticCoefficient)
                value = prediction.activity->osmotic_parameter_derivative.at(a);
            else {
                const auto ion = observation.species[0] == observation.ion_indices[0] ? 0 : 1;
                const auto column = prediction.reference->parameter_derivative.at(a)[ion];
                const auto ref = prediction.second_reference ? prediction.second_reference->parameter_derivative.at(a)[ion]
                    : ResidualEnergyAction{0.0, PropertyActionStatus::Available, {}};
                action = !column.action || !prediction.second_reference ? column : ref;
                value = column.action && ref.action ? GAS_CONSTANT_J_PER_MOL_K * observation.temperature * (*column.action - *ref.action) : NAN;
            }
            if (std::isfinite(value)) { action.action = value; action.status = PropertyActionStatus::Available; }
            else if (action.action || action.status == PropertyActionStatus::Nonfinite) {
                action.action.reset(); action.status = PropertyActionStatus::Nonfinite;
            }
            core::retain_first(gradient.first_numerical_failure, action.diagnostics.first_numerical_failure);
            gradient.parameter_actions.push_back(action);
            gradient.value.push_back(value);
        }
        gradient.status = ObservationStatus::Available;
        for (const auto& action : gradient.parameter_actions) if (!action.action) {
            gradient.status = ObservationStatus::ActionUnavailable;
            gradient.message = "electrolyte column: native action cause " + std::to_string(int(action.status));
            break;
        }
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
    const auto unavailable = [&](const std::string& what) {
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
        core::Scope scope("fit_density_column");
        const auto& center = prediction.state;
        const auto anchor = center && observation.branch != DensityBranch::Unspecified
            && properties::regular_density(center->molar_density, center->pressure_density_derivative, center->density_convergence)
            ? std::optional<double>(center->molar_density) : std::nullopt;
        auto batch = model.pressure_actions(observation.temperature, observation.pressure,
            observation.composition, observation.branch, observables, directions, anchor);
        const auto retain = [&] {
            core::retain_first(gradient.first_numerical_failure, batch.diagnostics.first_numerical_failure);
            if (batch.density_status != PropertyActionStatus::Available)
                core::retain_first(gradient.first_numerical_failure, batch.density_convergence);
        };
        retain();
        const auto cause = [&] {
            std::string name = equilibrium::solved_state_status_name(equilibrium::solved_state_status(batch.density_status));
            return batch.density_message.empty() ? name : name + ": " + batch.density_message;
        };
        if (batch.density_structure_refused) return unavailable(cause());
        if (anchor && (batch.density_status != PropertyActionStatus::Available || !batch.molar_density
            || !properties::regular_density(*batch.molar_density, batch.diagnostics.pressure_density_derivative, batch.density_convergence))) {
            if (batch.density_convergence) scope.finish(*batch.density_convergence);
            scope.probe("anchor_failure", {*anchor}, 0, batch.density_status == PropertyActionStatus::Available
                ? "target density regularity unavailable or exceeds w/100" : !batch.density_message.empty() ? batch.density_message.c_str()
                : equilibrium::solved_state_status_name(equilibrium::solved_state_status(batch.density_status)));
            batch = model.pressure_actions(observation.temperature, observation.pressure,
                observation.composition, observation.branch, observables, directions);
            retain();
        }
        if (batch.density_status != PropertyActionStatus::Available) return unavailable(cause());
        if (caloric && batch.results[0].status != PropertyActionStatus::Available)
            return unavailable(equilibrium::solved_state_status_name(equilibrium::solved_state_status(batch.results[0].status)));
        for (std::size_t k = 0; k < directions.size(); ++k)
            gradient.value[first + k] = (caloric ? *batch.results[0].action : batch.density_action[k]) * factor;
    }
    gradient.status = ObservationStatus::Available;
    return gradient;
}

std::vector<PredictionGradient> differentiate_group(
    const HelmholtzModel& model, const std::vector<const Observation*>& observations, const Prediction& prediction, const CoordinateSet& coordinates,
    const std::vector<ObservationInputDirection>& inputs
) {
    if (electrolyte_observation(observations.front()->kind) || pressure_closed(observations.front()->kind)) {
        std::vector<PredictionGradient> result;
        for (const auto* row : observations) {
            if (row->finite_solute_molalities && !inputs.empty()) {
                PredictionGradient unavailable;
                unavailable.message = "full solute input actions are not declared";
                result.push_back(std::move(unavailable)); continue;
            }
            auto projected = prediction;
            if (!pressure_closed(row->kind)) project_prediction(model, *row, projected);
            auto gradient = differentiate(model, *row, projected);
            if (gradient.status != ObservationStatus::Available) { result.push_back(std::move(gradient)); continue; }
            const auto parameter_width = model.active_parameters().size()
                + std::count_if(coordinates.coordinates.begin(), coordinates.coordinates.end(), [](const auto& c) { return c.reaction.has_value(); });
            gradient.value.resize(parameter_width, 0.0);
            for (std::size_t j = 0; j < inputs.size(); ++j) {
                const auto& d = inputs[j].state;
                double value = NAN;
                if (electrolyte_observation(row->kind)) {
                    auto action = electrolyte_input_action(model, row->temperature, row->pressure, row->composition,
                        row->ion_indices, row->stoichiometry, row->molality, *prediction.reference, *prediction.activity,
                        {d.temperature_k, d.pressure_pa, inputs[j].molality});
                    value = row->kind == ObservationKind::MeanIonicActivityCoefficient ? action.mean_ionic : action.osmotic;
                    if (action.status != PropertyActionStatus::Available) gradient.message = action.message;
                    gradient.input_actions.push_back(std::move(action));
                } else {
                    const bool density = row->kind == ObservationKind::Density;
                    const auto batch = model.pressure_actions(row->temperature, row->pressure, row->composition, row->branch,
                        density ? std::vector<PropertyObservable>{} : std::vector<PropertyObservable>{PropertyObservable::TotalIsobaricHeatCapacity},
                        {{d.temperature_k, d.pressure_pa, d.feed_amounts_mol, {}}}, prediction.state->molar_density);
                    core::retain_first(gradient.first_numerical_failure, batch.diagnostics.first_numerical_failure);
                    if (batch.density_status == PropertyActionStatus::Available && (density || batch.results[0].status == PropertyActionStatus::Available)) {
                        value = density ? batch.density_action[0] : *batch.results[0].action;
                        if (row->mass_basis) {
                            bool missing = false;
                            value *= mean_molar_mass(model.parameters(), row->composition.data(), missing);
                            value += *batch.molar_density * std::inner_product(d.feed_amounts_mol.begin(), d.feed_amounts_mol.end(),
                                model.parameters().molar_mass_kg_per_mol.begin(), 0.0);
                        }
                    } else gradient.message = "pressure input action: " + batch.density_message;
                }
                if (!std::isfinite(value)) { gradient.status = ObservationStatus::ActionUnavailable; gradient.value.clear(); break; }
                gradient.value.push_back(value);
            }
            result.push_back(std::move(gradient));
        }
        return result;
    }
    const auto& parameters = model.parameters();
    const auto components = parameters.component_count();
    const std::size_t active_width = model.active_parameters().size();
    std::vector<ReactionCoordinate> reactions;
    for (const auto& c : coordinates.coordinates) if (c.reaction) reactions.push_back(*c.reaction);
    const std::size_t parameter_width = active_width + reactions.size(), width = parameter_width + inputs.size();
    const auto& problem = prediction.compiled ? prediction.compiled->problem : observations.front()->problem;
    // One unit direction per active parameter at one shared center (one factorization, one reference ladder).
    std::vector<equilibrium::SolvedStateActionDirection> directions(parameter_width);
    for (std::size_t k = 0; k < parameter_width; ++k) {
        directions[k].feed_amounts_mol.assign(components, 0.0);
        directions[k].active_parameters.assign(active_width, 0.0);
        if (k < active_width) directions[k].active_parameters[k] = 1.0;
        else {
            directions[k].reaction_coefficients.assign(problem.reactions.size(), {});
            const auto& target = reactions[k - active_width];
            for (std::size_t j = 0; j < problem.reactions.size(); ++j)
                if (problem.reactions[j].name == target.name)
                    directions[k].reaction_coefficients[j] = reaction_direction(target);
        }
    }
    for (const auto& input : inputs) directions.push_back(input.state);
    if (addition(observations.front()->kind)) {
        auto item = *observations.front();
        const bool capped = item.kind == ObservationKind::SaturationCappedAdditionEnthalpy;
        std::vector<double> loading_gradient(width, 0);
        if (capped) {
            const auto& search = *prediction.saturation;
            const auto& root = search.candidates[*search.accepted].state;
            const auto c = equilibrium::compile_problem(model, *root.declaration);
            const auto e = equilibrium::equilibrium_numerical_estimates(c, root);
            const auto columns = directions.empty() ? std::vector<equilibrium::SolvedStateActionBatch>()
                : equilibrium::solved_state_first_actions(c, root, {{equilibrium::SolvedStateObservableKind::Loading}}, directions);
            for (std::size_t k = 0; k < width; ++k) {
                if (!columns[k].results[0].action)
                    return {{{}, ObservationStatus::ActionUnavailable, std::string("saturation:")
                        + equilibrium::solved_state_status_name(columns[k].results[0].status) + ":"
                        + columns[k].results[0].diagnostic_message.value_or("")}};
                loading_gradient[k] = *columns[k].results[0].action;
            }
            if (prediction.cap_branch == "switch" || !e.loading.resolved
                || std::abs(item.loading_cap - *prediction.saturation_loading) <= e.loading.absolute) {
                PredictionGradient unavailable; unavailable.message = "saturation:switch or unresolved branch";
                unavailable.saturation_loading_gradient = loading_gradient; return {unavailable};
            }
            item.feed_amount = *prediction.added_amount;
        }
        item.kind = ObservationKind::SystemEnthalpy;
        auto initial = differentiate_group(model, {&item}, prediction, coordinates).front();
        if (initial.status != ObservationStatus::Available) { initial.message = "initial:" + initial.message; return {initial}; }
        auto center = prediction;
        std::swap(center.solved, center.final_solved); std::swap(center.compiled, center.final_compiled);
        item.problem = item.final_problem;
        auto final = differentiate_group(model, {&item}, center, coordinates).front();
        core::retain_first(final.first_numerical_failure, initial.first_numerical_failure, true);
        if (final.status != ObservationStatus::Available) { final.message = "final:" + final.message; return {final}; }
        item.kind = ObservationKind::FeedEnthalpy;
        auto feed = differentiate(model, item, prediction);
        core::retain_first(feed.first_numerical_failure, final.first_numerical_failure, true);
        if (feed.status != ObservationStatus::Available) { feed.message = "feed:" + feed.message; return {feed}; }
        core::retain_first(initial.first_numerical_failure, feed.first_numerical_failure);
        for (std::size_t k = 0; k < width; ++k)
            initial.value[k] = (final.value[k] - initial.value[k]) / item.feed_amount - (k < feed.value.size() ? feed.value[k] : 0.0);
        if (capped) {
            initial.saturation_loading_gradient = loading_gradient;
            initial.endpoint_loading_gradient.assign(width, 0);
            if (prediction.cap_branch == "saturation_cap") for (std::size_t k = 0; k < width; ++k) {
                initial.endpoint_loading_gradient[k] = loading_gradient[k];
                initial.value[k] -= *prediction.enthalpy_difference / item.feed_amount / *prediction.endpoint_loading * loading_gradient[k];
            }
        }
        return {initial};
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
        case ObservationKind::EquilibriumIsobaricHeatCapacity:
            add(observable(Kind::EquilibriumIsobaricHeatCapacity, 0));
            break;
        case ObservationKind::SystemEnthalpy:
            add(observable(Kind::SystemEnthalpy, 0));
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
            // NOLINTNEXTLINE(clang-analyzer-core.NullDereference): only a LogEquilibriumConstant row lacks a solved state, and same_state never groups it.
            for (std::size_t q = 0; q < solved->phase_count; ++q) add(observable(Kind::PhaseAmount, q));
            break;
        case ObservationKind::PhaseDensity:
            add(observable(Kind::PhaseMolarDensity, phase));
            if (observation.mass_basis)
                for (std::size_t i = 0; i < components; ++i)
                    add(observable(Kind::PhaseMoleFraction, phase, i));
            break;
        case ObservationKind::SaturationCappedAdditionEnthalpy:
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
        const double mass = solved && observation.mass_basis ? mean_molar_mass(parameters, fractions, missing_mass) : 1.0;
        auto& gradient = gradients[row];
        gradient.value.assign(width, 0.0);
        gradient.status = ObservationStatus::Available;
        for (const auto& batch : batches) gradient.state_diagnostics.push_back(batch.diagnostics);
        for (std::size_t k = 0; k < width; ++k) {
            if (observation.kind == ObservationKind::LogEquilibriumConstant) {
                for (std::size_t j = 0; j < problem.reactions.size(); ++j)
                    if (problem.reactions[j].name == observation.reaction) {
                        const auto& law = problem.reactions[j].correlation;
                        if (k >= active_width && k < parameter_width) {
                            const auto& d = directions[k].reaction_coefficients[j];
                            gradient.value[k] = reaction_log_value(ReactionLogPolynomial{d[0], d[1], d[2], d[3],
                                law ? law->reference_temperature : 1.0, 0, 0, {}}, problem.temperature.value);
                        } else if (law && k >= parameter_width) {
                            gradient.value[k] = directions[k].temperature_k * reaction_log_value(*law, problem.temperature.value, 1);
                        }
                    }
                continue;
            }
            if (!solved) {
                gradient.value.clear(); gradient.status = ObservationStatus::ActionUnavailable;
                gradient.message = "central solved state unavailable"; break;
            }
            std::vector<double> action;
            for (const auto j : slots[row]) {
                const auto& item = batches[k].results[j];
                if (item.property_diagnostics) core::retain_first(gradient.first_numerical_failure, item.property_diagnostics->first_numerical_failure);
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
                for (std::size_t i = 0; i + 1 < action.size(); ++i) value += *solved->pressure * action[i];
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
