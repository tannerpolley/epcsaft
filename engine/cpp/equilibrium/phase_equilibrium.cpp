#include "phase_equilibrium.hpp"

#include "grepe.hpp"
#include "phase_stability_coordinates.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

namespace epcsaft_engine::equilibrium {
namespace {

std::vector<double> normalize(const std::vector<double>& amounts) {
    if (amounts.empty()
        || std::any_of(amounts.begin(), amounts.end(), [](double value) {
            return !std::isfinite(value) || value <= 0.0;
        })) {
        throw std::invalid_argument("phase feed must be positive and finite");
    }
    const double total = std::accumulate(amounts.begin(), amounts.end(), 0.0);
    std::vector<double> result = amounts;
    for (double& value : result) {
        value /= total;
    }
    return result;
}

PhaseEquilibriumState all_component_state(
    const ProviderContext& provider,
    PhaseRole role,
    double temperature_k,
    const std::vector<double>& composition,
    double molar_volume_m3_per_mol
) {
    const auto charges = phase_provider_component_charges(provider);
    const auto phase = evaluate_nonreactive_phase(
        provider,
        charges,
        temperature_k,
        composition,
        molar_volume_m3_per_mol
    );
    return make_phase_equilibrium_state(
        provider,
        role,
        PhaseSupport::AllComponents,
        temperature_k,
        composition,
        molar_volume_m3_per_mol,
        phase
    );
}

}  // namespace

const char* phase_role_name(PhaseRole role) noexcept {
    return role == PhaseRole::Liquid ? "liquid" : "vapor";
}

const char* phase_support_name(PhaseSupport support) noexcept {
    if (support == PhaseSupport::AllComponents) {
        return "all_components";
    }
    return support == PhaseSupport::NeutralComponents
        ? "neutral_components"
        : "declared_neutral_components";
}

PhaseEquilibriumState audit_phase_state(
    const ProviderContext& provider,
    PhaseRole role,
    double temperature_k,
    const std::vector<double>& mole_fractions,
    double molar_volume_m3_per_mol
) {
    return all_component_state(
        provider,
        role,
        temperature_k,
        mole_fractions,
        molar_volume_m3_per_mol
    );
}

FixedTpPhaseResult solve_declared_fixed_tp(
    const ProviderContext& provider,
    const FixedTpPhaseProblem& problem
) {
    FixedTpPhaseResult result;
    result.parameter_fingerprint = provider.fingerprint();
    result.diagnostics.route = "grepe_declared_topology";
    result.diagnostics.certificate_scope = "fixed_tp_declared_topology";
    result.diagnostics.phase_selection_claim =
        "locally_certified_declared_topology";
    result.diagnostics.search_scope = "declared_topology_grepe";
    result.diagnostics.search_budget = "grepe_finite_phase_starts";
    result.diagnostics.search_exclusions = {
        "undeclared_phase_counts",
        "undeclared_component_supports",
    };
    result.diagnostics.completeness_criterion =
        "all_returned_phases_pass_grepe_local_equilibrium_certificates";
    result.diagnostics.hint_count = static_cast<int>(problem.hints.size());
    result.diagnostics.initial_guess_count =
        static_cast<int>(problem.initial_guesses.size());
    try {
        const auto component_ids = phase_provider_component_ids(provider);
        if (!std::isfinite(problem.temperature_k)
            || problem.temperature_k <= 0.0
            || !std::isfinite(problem.pressure_pa)
            || problem.pressure_pa <= 0.0
            || problem.feed_amounts_mol.size() != component_ids.size()
            || problem.phases.empty()
            || problem.phases.size() > 3
            || std::count_if(
                problem.phases.begin(),
                problem.phases.end(),
                [](const auto& phase) {
                    return phase.role == PhaseRole::Vapor;
                }
            ) > 1
            || std::any_of(
                problem.phases.begin(),
                problem.phases.end(),
                [](const auto& phase) {
                    return phase.support != PhaseSupport::AllComponents;
                }
            )) {
            throw std::invalid_argument(
                "fixed-T,P topology requires one to three all-component phases and at most one vapor"
            );
        }
        if (!problem.hints.empty() && !problem.initial_guesses.empty()) {
            throw std::invalid_argument(
                "phase hints and declared phase initial guesses are mutually exclusive"
            );
        }
        if (!problem.initial_guesses.empty()
            && problem.initial_guesses.size() != problem.phases.size()) {
            throw std::invalid_argument(
                "declared phase initial guesses must cover every topology phase"
            );
        }
        const double total = std::accumulate(
            problem.feed_amounts_mol.begin(),
            problem.feed_amounts_mol.end(),
            0.0
        );
        static_cast<void>(normalize(problem.feed_amounts_mol));
        GrepeRequest request;
        request.identity = "nonreactive_declared_phase_equilibrium";
        request.temperature.value = problem.temperature_k;
        request.pressure.value = problem.pressure_pa;
        request.reaction = make_nonreactive_reaction_system(
            provider,
            problem.temperature_k,
            problem.pressure_pa,
            problem.feed_amounts_mol
        );
        request.trace_floor = 0.0;
        for (const auto& hint : problem.hints) {
            if (hint.phase_index >= problem.phases.size()
                || hint.mole_fractions.size() != provider.component_count()
                || !std::isfinite(hint.molar_density_mol_m3)
                || hint.molar_density_mol_m3 <= 0.0) {
                throw std::invalid_argument("declared phase hint is invalid");
            }
            request.start_order_hints.push_back({
                hint.phase_index,
                normalize(hint.mole_fractions),
                1.0 / hint.molar_density_mol_m3,
            });
        }
        std::vector<std::size_t> support(provider.component_count());
        std::iota(support.begin(), support.end(), 0);
        for (std::size_t phase = 0; phase < problem.phases.size(); ++phase) {
            GrepePhase declared;
            declared.identity = "phase_" + std::to_string(phase);
            declared.role = problem.phases[phase].role;
            declared.model_reference_id = provider.fingerprint();
            declared.support_name = "all_components";
            declared.support_indices = support;
            declared.packing_fraction_interval = {
                kPhaseStabilityPackingFractionMinimum,
                kPhaseStabilityPackingFractionMaximum,
            };
            request.phases.push_back(std::move(declared));
        }
        if (!problem.initial_guesses.empty()) {
            request.declared_starts_only = true;
            std::vector<bool> assigned(problem.phases.size(), false);
            double fraction_sum = 0.0;
            for (const auto& guess : problem.initial_guesses) {
                if (guess.phase_index >= problem.phases.size()
                    || assigned[guess.phase_index]
                    || guess.mole_fractions.size() != provider.component_count()
                    || !std::isfinite(guess.molar_density_mol_m3)
                    || guess.molar_density_mol_m3 <= 0.0
                    || !std::isfinite(guess.phase_fraction)
                    || guess.phase_fraction <= 0.0) {
                    throw std::invalid_argument(
                        "declared phase initial guess is invalid"
                    );
                }
                const auto composition = normalize(guess.mole_fractions);
                ReactiveFinitePhaseStart start;
                for (double fraction : composition) {
                    start.amounts_mol.push_back(
                        total * guess.phase_fraction * fraction
                    );
                }
                start.molar_volume_m3_per_mol =
                    1.0 / guess.molar_density_mol_m3;
                request.phases[guess.phase_index].start = std::move(start);
                assigned[guess.phase_index] = true;
                fraction_sum += guess.phase_fraction;
            }
            if (std::abs(fraction_sum - 1.0) > 1.0e-12) {
                throw std::invalid_argument(
                    "declared phase initial-guess fractions must sum to one"
                );
            }
            auto& closure = request.phases.back().start->amounts_mol;
            for (std::size_t component = 0;
                 component < provider.component_count(); ++component) {
                double retained = problem.feed_amounts_mol[component];
                for (std::size_t phase = 0;
                     phase + 1 < request.phases.size(); ++phase) {
                    retained -= request.phases[phase].start->amounts_mol[component];
                }
                if (!(retained > 0.0)) {
                    throw std::invalid_argument(
                        "declared phase initial guesses cannot be balanced with positive amounts"
                    );
                }
                closure[component] = retained;
            }
        }
        GrepeResult solved = solve_grepe(
            provider.model(), provider.model(), request
        );
        const auto& equilibrium = solved.equilibrium;
        result.diagnostics.solver_status = equilibrium.solver_status;
        result.diagnostics.numerical_status = equilibrium.numerical_status;
        result.diagnostics.physical_status = equilibrium.physical_status;
        result.diagnostics.provider_evaluation_count = static_cast<int>(
            equilibrium.provider_evaluation_count
        );
        result.diagnostics.optimizer_iteration_count = static_cast<int>(
            equilibrium.optimizer_iterations
        );
        result.diagnostics.material_balance_inf = equilibrium.balance_inf_norm;
        result.diagnostics.pressure_residual_inf =
            equilibrium.pressure_relative_inf_norm;
        result.diagnostics.chemical_potential_residual_inf =
            equilibrium.transfer_chemical_potential_inf_norm;
        result.diagnostics.kkt_stationarity_inf =
            equilibrium.kkt_stationarity_inf_norm;
        result.diagnostics.phase_distinctness = equilibrium.phase_distinctness;
        if (!equilibrium.accepted
            || equilibrium.phase_amounts_mol.size() != problem.phases.size()
            || equilibrium.phase_volumes_m3.size() != problem.phases.size()) {
            result.diagnostics.failure_stage = "grepe";
            result.diagnostics.failure_reason = equilibrium.failure_reason.empty()
                ? solved.failure_code : equilibrium.failure_reason;
            return result;
        }
        std::vector<SolvedPhaseState> sensitivity_phases;
        for (std::size_t phase = 0; phase < problem.phases.size(); ++phase) {
            const auto& amounts = equilibrium.phase_amounts_mol[phase];
            const double amount = std::accumulate(
                amounts.begin(), amounts.end(), 0.0
            );
            if (!(amount > 0.0)) {
                throw std::runtime_error("GREPE returned an empty finite phase");
            }
            std::vector<double> composition = amounts;
            for (double& value : composition) {
                value /= amount;
            }
            const double molar_volume =
                equilibrium.phase_volumes_m3[phase] / amount;
            result.phases.push_back(all_component_state(
                provider,
                problem.phases[phase].role,
                problem.temperature_k,
                composition,
                molar_volume
            ));
            result.phase_amounts_mol.push_back(amount);
            result.phase_fractions.push_back(amount / total);
            sensitivity_phases.push_back({
                amounts,
                equilibrium.phase_volumes_m3[phase],
            });
        }
        result.total_free_energy_over_rt = equilibrium.objective_over_rt;
        result.diagnostics.numerical_status = "certified";
        result.diagnostics.physical_status = "accepted";
        result.diagnostics.accepted_start_count = 1;
        result.diagnostics.distinct_root_count = 1;
        result.diagnostics.hint_reordered_start_count =
            static_cast<int>(problem.hints.size());
        if (!provider.active_parameters().identities.empty()) {
            result.sensitivity = evaluate_declared_state_sensitivity(
                provider,
                problem.temperature_k,
                problem.pressure_pa,
                component_ids,
                problem.feed_amounts_mol,
                sensitivity_phases
            );
        }
    } catch (const std::exception& error) {
        result.phases.clear();
        result.phase_amounts_mol.clear();
        result.phase_fractions.clear();
        result.diagnostics.solver_status = "input_rejected";
        result.diagnostics.numerical_status = "not_evaluated";
        result.diagnostics.physical_status = "rejected";
        result.diagnostics.failure_stage = "problem_resolution";
        result.diagnostics.failure_reason = error.what();
    }
    return result;
}

}  // namespace epcsaft_engine::equilibrium
