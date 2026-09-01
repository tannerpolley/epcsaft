#include "phase_stability.hpp"

#include "grepe.hpp"
#include "phase_stability_tolerances.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <ctime>
#include <numeric>
#include <optional>
#include <stdexcept>
#include <utility>

namespace epcsaft_engine::equilibrium {
namespace {

constexpr double kGasConstant = 8.31446261815324;
constexpr double kLiquidCompressibilityMaximum = 0.25;
constexpr double kVaporCompressibilityMinimum = 0.45;
constexpr std::size_t kCandidatePoolLimit = 32;

template <typename Callable>
auto run_stage(int stage, Callable&& callable) {
    const auto wall_start = std::chrono::steady_clock::now();
    const std::clock_t cpu_start = std::clock();
    auto result = callable();
    result.timing.step = stage;
    result.timing.wall_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - wall_start
    ).count();
    result.timing.cpu_seconds = static_cast<double>(
        std::clock() - cpu_start
    ) / static_cast<double>(CLOCKS_PER_SEC);
    return result;
}

std::vector<double> independent_coordinates(
    const PhaseStabilityCoordinates& coordinates,
    const std::vector<double>& physical
) {
    const auto modified = phase_stability_transform_physical_fractions(
        coordinates, physical
    );
    std::vector<double> result;
    result.reserve(coordinates.independent_indices.size());
    for (std::size_t index : coordinates.independent_indices) {
        const auto found = std::find(
            coordinates.retained_indices.begin(),
            coordinates.retained_indices.end(),
            index
        );
        result.push_back(modified[static_cast<std::size_t>(
            found - coordinates.retained_indices.begin()
        )]);
    }
    return result;
}

std::optional<PhaseRole> mechanical_role(
    double pressure_pa,
    double temperature_k,
    double molar_volume_m3_per_mol
) {
    const double z = pressure_pa * molar_volume_m3_per_mol
        / (kGasConstant * temperature_k);
    if (z <= kLiquidCompressibilityMaximum) {
        return PhaseRole::Liquid;
    }
    if (z >= kVaporCompressibilityMinimum) {
        return PhaseRole::Vapor;
    }
    return std::nullopt;
}

std::optional<std::array<double, 2>> balanced_pair_seed_fractions(
    const std::vector<double>& feed,
    const std::vector<double>& first,
    const std::vector<double>& second
) {
    if (feed.size() != first.size() || feed.size() != second.size()) {
        return std::nullopt;
    }
    std::size_t pivot = feed.size();
    double largest_difference = 0.0;
    for (std::size_t component = 0; component < feed.size(); ++component) {
        const double difference = std::abs(first[component] - second[component]);
        if (difference > largest_difference) {
            largest_difference = difference;
            pivot = component;
        }
    }
    if (pivot == feed.size()) {
        return std::nullopt;
    }
    const double first_fraction = (feed[pivot] - second[pivot])
        / (first[pivot] - second[pivot]);
    if (!(first_fraction > 0.0 && first_fraction < 1.0)) {
        return std::nullopt;
    }
    for (std::size_t component = 0; component < feed.size(); ++component) {
        const double reconstructed = first_fraction * first[component]
            + (1.0 - first_fraction) * second[component];
        if (std::abs(reconstructed - feed[component]) > 1.0e-8) {
            return std::nullopt;
        }
    }
    return std::array<double, 2>{first_fraction, 1.0 - first_fraction};
}

struct CandidateSolve {
    GrepeResult solved;
    std::vector<std::uint64_t> stable_ids;
};

std::optional<CandidateSolve> solve_candidate(
    const PhaseStabilityThermodynamicAccess& thermodynamics,
    const FlashInput& input,
    const PhaseStabilityCoordinates& coordinates,
    const std::vector<PhaseStabilityMPoint>& candidates,
    const std::vector<std::size_t>& selected,
    const std::vector<double>& fractions
) {
    GrepeRequest request;
    request.identity = "nonreactive_phase_stability_candidate";
    request.temperature.value = input.temperature_k;
    request.pressure.value = input.pressure_pa;
    request.reaction = make_nonreactive_reaction_system(
        *thermodynamics.provider,
        input.temperature_k,
        input.pressure_pa,
        input.overall_mole_fractions
    );
    request.trace_floor = 0.0;
    request.declared_starts_only = true;
    std::vector<std::size_t> support(input.overall_mole_fractions.size());
    std::iota(support.begin(), support.end(), 0);
    CandidateSolve result;
    std::size_t vapor_count = 0;
    for (std::size_t phase = 0; phase < selected.size(); ++phase) {
        const PhaseStabilityMPoint& candidate = candidates[selected[phase]];
        const auto role = mechanical_role(
            input.pressure_pa, input.temperature_k, candidate.volume
        );
        if (!role.has_value()) {
            return std::nullopt;
        }
        vapor_count += *role == PhaseRole::Vapor ? 1u : 0u;
        if (vapor_count > 1) {
            return std::nullopt;
        }
        const auto physical = phase_stability_lift_independent_fractions(
            coordinates, candidate.independent_modified_fractions
        );
        GrepePhase declared;
        declared.identity = "phase_stability_card_"
            + std::to_string(selected[phase]);
        declared.role = *role;
        declared.model_reference_id = thermodynamics.provider_fingerprint;
        declared.support_name = "all_components";
        declared.support_indices = support;
        declared.packing_fraction_interval = {
            kPhaseStabilityPackingFractionMinimum,
            kPhaseStabilityPackingFractionMaximum,
        };
        ReactiveFinitePhaseStart start;
        start.molar_volume_m3_per_mol = candidate.volume;
        for (double component : physical) {
            start.amounts_mol.push_back(fractions[phase] * component);
        }
        declared.start = std::move(start);
        request.phases.push_back(std::move(declared));
        result.stable_ids.push_back(candidate.insertion_id);
    }
    try {
        result.solved = solve_grepe(
            thermodynamics.provider->model(),
            thermodynamics.provider->model(),
            request
        );
    } catch (const std::exception& error) {
        throw std::runtime_error(
            std::string("candidate_grepe_failure: ") + error.what()
        );
    }
    if (!result.solved.equilibrium.accepted) {
        const auto& equilibrium = result.solved.equilibrium;
        throw std::runtime_error(
            "candidate_grepe_rejected: "
            + equilibrium.failure_reason
            + "; solver=" + equilibrium.solver_status
            + "; iterations=" + std::to_string(equilibrium.optimizer_iterations)
            + "; seed_fraction_0=" + std::to_string(fractions.front())
            + "; balance=" + std::to_string(equilibrium.balance_inf_norm)
            + "; pressure=" + std::to_string(equilibrium.pressure_relative_inf_norm)
            + "; transfer=" + std::to_string(
                equilibrium.transfer_chemical_potential_inf_norm
            )
            + "; kkt=" + std::to_string(equilibrium.kkt_stationarity_inf_norm)
            + "; distinct=" + std::to_string(equilibrium.phase_distinctness)
            + (equilibrium.callback_error.empty()
                ? std::string{}
                : "; callback=" + equilibrium.callback_error)
        );
    }
    return result;
}

std::vector<PhaseStabilityPhase> map_phases(
    const PhaseStabilityThermodynamicAccess& thermodynamics,
    const FlashInput& input,
    const PhaseStabilityCoordinates& coordinates,
    const CandidateSolve& candidate
) {
    const auto& equilibrium = candidate.solved.equilibrium;
    std::vector<PhaseStabilityPhase> phases;
    for (std::size_t index = 0; index < equilibrium.phase_amounts_mol.size(); ++index) {
        const auto& amounts = equilibrium.phase_amounts_mol[index];
        const double amount = std::accumulate(amounts.begin(), amounts.end(), 0.0);
        std::vector<double> physical = amounts;
        for (double& value : physical) {
            value /= amount;
        }
        const double volume = equilibrium.phase_volumes_m3[index] / amount;
        const auto independent = independent_coordinates(coordinates, physical);
        const auto state = thermodynamics.evaluate(independent, std::log(volume));
        phases.push_back({
            candidate.stable_ids[index],
            amount,
            independent,
            physical,
            volume,
            equilibrium.phase_packing_fractions[index],
            state.helmholtz_over_rt_reference_amount,
            equilibrium.phase_pressures_pa[index],
            equilibrium.phase_chemical_potentials_over_rt[index],
            state.objective,
            state.gradient,
        });
    }
    if (std::any_of(phases.begin(), phases.end(), [&](const auto& phase) {
            return !mechanical_role(
                input.pressure_pa, input.temperature_k, phase.volume
            ).has_value();
        })) {
        throw std::runtime_error("role_classification_ambiguous");
    }
    std::stable_sort(phases.begin(), phases.end(), [&](const auto& left, const auto& right) {
        const auto left_role = mechanical_role(
            input.pressure_pa, input.temperature_k, left.volume
        );
        const auto right_role = mechanical_role(
            input.pressure_pa, input.temperature_k, right.volume
        );
        return std::pair(*left_role == PhaseRole::Vapor, left.volume)
            < std::pair(*right_role == PhaseRole::Vapor, right.volume);
    });
    return phases;
}

void retain_grepe_evidence(
    PhaseStabilityDiscoveryResult& result,
    const CandidateSolve& candidate,
    const std::vector<double>& charges
) {
    const auto& equilibrium = candidate.solved.equilibrium;
    result.material_balance_inf = equilibrium.balance_inf_norm;
    result.pressure_residual_inf = equilibrium.pressure_relative_inf_norm;
    result.chemical_potential_residual_inf =
        equilibrium.transfer_chemical_potential_inf_norm;
    result.kkt_stationarity_inf = equilibrium.kkt_stationarity_inf_norm;
    result.phase_distinctness = equilibrium.phase_distinctness;
    double charge_inf = 0.0;
    for (const auto& amounts : equilibrium.phase_amounts_mol) {
        const double total = std::accumulate(
            amounts.begin(), amounts.end(), 0.0
        );
        double charge = 0.0;
        for (std::size_t component = 0; component < amounts.size(); ++component) {
            charge += charges[component] * amounts[component] / total;
        }
        charge_inf = std::max(charge_inf, std::abs(charge));
    }
    result.electroneutrality_inf = charge_inf;
}

}  // namespace

PhaseStabilityDiscoveryResult run_phase_stability(
    const PhaseStabilityThermodynamicAccess& thermodynamics,
    const FlashInput& input,
    const PhaseStabilityResourceProfile& resources
) {
    PhaseStabilityDiscoveryResult result;
    const auto fail = [&](const char* stage, const std::string& reason) {
        result.failure_stage = stage;
        result.failure_reason = reason;
        return result;
    };
    if (thermodynamics.provider == nullptr
        || thermodynamics.provider_fingerprint.empty()
        || !thermodynamics.evaluate
        || !thermodynamics.volume_bounds_physical
        || !thermodynamics.packing_fraction) {
        return fail("stage1", "invalid_thermodynamic_access");
    }
    result.step1 = run_stage(1, [&] {
        return run_phase_stability_stage1(
            thermodynamics.component_ids,
            thermodynamics.charges,
            input.temperature_k,
            input.pressure_pa,
            input.overall_mole_fractions,
            thermodynamics.volume_bounds_physical,
            thermodynamics.provider_fingerprint,
            thermodynamics.total_ion_mole_fraction_max
        );
    });
    result.step_timings.push_back(result.step1.timing);
    if (result.step1.status != "complete") {
        return fail("stage1", result.step1.reason);
    }
    result.step2 = run_stage(2, [&] {
        return run_phase_stability_stage2(
            result.step1,
            thermodynamics.evaluate,
            resources.step2_search_work_budget
        );
    });
    result.step_timings.push_back(result.step2->timing);
    if (result.step2->outcome == PhaseStabilityStage2Outcome::Indeterminate) {
        return fail("stage2", result.step2->reason);
    }
    const auto& coordinates = *result.step1.coordinates;
    if (result.step2->outcome
        == PhaseStabilityStage2Outcome::NoNegativeWitnessDetected) {
        for (int stage = 3; stage <= 7; ++stage) {
            PhaseStabilityStageTiming skipped;
            skipped.step = stage;
            skipped.terminal_status = "not_applicable";
            skipped.terminal_reason = "skipped_after_no_witness";
            result.step_timings.push_back(std::move(skipped));
        }
        PhaseStabilityMPoint reference{
            0,
            *result.step1.independent_feed,
            result.step2->reference->volume,
            thermodynamics.packing_fraction(
                *result.step1.independent_feed,
                result.step2->reference->volume
            ),
            result.step2->reference->objective,
            result.step2->reference->gradient,
            "feed_reference",
        };
        std::optional<CandidateSolve> solved;
        try {
            solved = solve_candidate(
                thermodynamics, input, coordinates, {reference}, {0}, {1.0}
            );
        } catch (const std::exception& error) {
            return fail("grepe", error.what());
        }
        if (!solved.has_value()) {
            return fail("grepe", "feed_reference_grepe_rejected");
        }
        retain_grepe_evidence(result, *solved, thermodynamics.charges);
        try {
            result.phases = map_phases(
                thermodynamics, input, coordinates, *solved
            );
        } catch (const std::exception& error) {
            return fail("grepe", error.what());
        }
        result.total_free_energy_over_rt =
            solved->solved.equilibrium.objective_over_rt;
        result.outcome = "one_phase_no_negative_witness_detected";
        return result;
    }
    const PhaseStabilityPressureRootEvaluator pressure_roots = [&](const auto& composition) {
        return evaluate_phase_stability_pressure_envelope(
            composition,
            (*result.step1.volume_bounds)(composition),
            thermodynamics.evaluate,
            64,
            8
        );
    };
    result.step3 = run_stage(3, [&] {
        return run_phase_stability_stage3(
            result.step1, *result.step2, pressure_roots
        );
    });
    result.step_timings.push_back(result.step3->timing);
    if (result.step3->status != "complete" || !result.step3->state) {
        return fail("stage3", result.step3->reason);
    }
    PhaseStabilityPersistentState state = *result.step3->state;
    while (true) {
        result.step4_history.push_back(run_stage(4, [&] {
            return run_phase_stability_stage4(state);
        }));
        const auto& stage4 = result.step4_history.back();
        result.step_timings.push_back(stage4.timing);
        if (stage4.status != "complete") {
            return fail("stage4", stage4.reason);
        }
        result.step5_history.push_back(run_stage(5, [&] {
            return run_phase_stability_stage5(
                result.step1,
                stage4,
                state,
                thermodynamics.evaluate,
                resources
            );
        }));
        const auto& stage5 = result.step5_history.back();
        result.step_timings.push_back(stage5.timing);
        if (stage5.status != "complete") {
            return fail("stage5", stage5.reason);
        }
        result.step6_history.push_back(run_stage(6, [&] {
            return run_phase_stability_stage6(
                result.step1,
                stage4,
                state,
                thermodynamics.packing_fraction
            );
        }));
        const auto& stage6 = result.step6_history.back();
        result.step_timings.push_back(stage6.timing);
        if (stage6.status != "complete") {
            return fail("stage6", stage6.reason);
        }
        std::vector<PhaseStabilityMPoint> qualified_candidates{{
            0,
            *result.step1.independent_feed,
            result.step2->reference->volume,
            thermodynamics.packing_fraction(
                *result.step1.independent_feed,
                result.step2->reference->volume
            ),
            result.step2->reference->objective,
            result.step2->reference->gradient,
            "feed_reference",
        }};
        for (const auto& candidate : stage6.candidates) {
            try {
                const auto envelope = pressure_roots(
                    candidate.independent_modified_fractions
                );
                const double candidate_log_volume = std::log(candidate.volume);
                const std::size_t matches = static_cast<std::size_t>(std::count_if(
                    envelope.roots.begin(),
                    envelope.roots.end(),
                    [&](const auto& root) {
                        return std::abs(root.log_volume - candidate_log_volume)
                            <= 1.0e-8;
                    }
                ));
                const bool duplicate = std::any_of(
                    qualified_candidates.begin(),
                    qualified_candidates.end(),
                    [&](const auto& retained) {
                        if (std::abs(
                                candidate.packing_fraction
                                - retained.packing_fraction
                            ) >= kPhaseStabilityPaperStep6PackingDistinct.atol) {
                            return false;
                        }
                        for (std::size_t coordinate = 0;
                             coordinate
                                < candidate.independent_modified_fractions.size();
                             ++coordinate) {
                            if (std::abs(
                                    candidate.independent_modified_fractions[coordinate]
                                    - retained.independent_modified_fractions[coordinate]
                                ) >= kPhaseStabilityPaperStep6CompositionDistinct.atol) {
                                return false;
                            }
                        }
                        return true;
                    }
                );
                if (matches == 1 && !duplicate) {
                    qualified_candidates.push_back(candidate);
                }
            } catch (const std::exception&) {
                continue;
            }
        }
        std::optional<CandidateSolve> selected;
        const std::size_t limit = qualified_candidates.size();
        const std::size_t topology_count = limit < 2
            ? 0
            : limit * (limit - 1) / 2;
        if (topology_count > kCandidatePoolLimit) {
            return fail(
                "stage7",
                "candidate_subset_budget_exhausted: candidates="
                    + std::to_string(limit)
                    + ", topologies=" + std::to_string(topology_count)
            );
        }
        std::vector<std::string> candidate_failures;
        for (std::size_t left = 0;
             left < limit;
             ++left) {
            const auto left_physical = phase_stability_lift_independent_fractions(
                coordinates,
                qualified_candidates[left].independent_modified_fractions
            );
            for (std::size_t right = left + 1; right < limit; ++right) {
                const auto right_physical = phase_stability_lift_independent_fractions(
                    coordinates,
                    qualified_candidates[right].independent_modified_fractions
                );
                const auto fractions = balanced_pair_seed_fractions(
                    input.overall_mole_fractions,
                    left_physical,
                    right_physical
                );
                if (!fractions.has_value()) {
                    continue;
                }
                std::optional<CandidateSolve> solved;
                try {
                    solved = solve_candidate(
                        thermodynamics,
                        input,
                        coordinates,
                        qualified_candidates,
                        {left, right},
                        {(*fractions)[0], (*fractions)[1]}
                    );
                } catch (const std::exception& error) {
                    candidate_failures.push_back(error.what());
                    continue;
                }
                if (solved.has_value()
                    && (!selected.has_value()
                        || solved->solved.equilibrium.objective_over_rt
                            < selected->solved.equilibrium.objective_over_rt)) {
                    selected = std::move(solved);
                }
            }
        }
        if (selected.has_value()) {
            PhaseStabilityStageTiming handoff;
            handoff.step = 7;
            handoff.invocation_count = 1;
            handoff.terminal_status = "complete";
            handoff.terminal_reason = "handoff_to_grepe";
            result.step_timings.push_back(std::move(handoff));
            retain_grepe_evidence(result, *selected, thermodynamics.charges);
            try {
                result.phases = map_phases(
                    thermodynamics, input, coordinates, *selected
                );
            } catch (const std::exception& error) {
                return fail("grepe", error.what());
            }
            result.total_free_energy_over_rt =
                selected->solved.equilibrium.objective_over_rt;
            result.outcome = "physical_equilibrium_accepted";
            return result;
        }
        result.step7_history.push_back(run_stage(7, [&] {
            return run_phase_stability_stage7(
                state, stage5, stage6, resources
            );
        }));
        result.step_timings.push_back(result.step7_history.back().timing);
        if (result.step7_history.back().status != "complete") {
            if (!candidate_failures.empty()) {
                return fail("grepe", candidate_failures.back());
            }
            return fail("stage7", result.step7_history.back().reason);
        }
    }
}

}  // namespace epcsaft_engine::equilibrium
