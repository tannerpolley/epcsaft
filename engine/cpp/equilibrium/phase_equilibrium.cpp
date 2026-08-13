#include "phase_equilibrium.hpp"

#include "flash.hpp"
#include "held2.hpp"
#include "held2_step8.hpp"
#include "held2_tolerances.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

namespace epcsaft_engine::equilibrium {
namespace {

constexpr double kGasConstant = 8.31446261815324;
constexpr double kRoleZSplit = 0.35;
constexpr double kChemicalPotentialResidualMax = 1.0e-6;
constexpr double kDistinctTolerance = 1.0e-5;

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

std::vector<double> independent(
    const Held2Coordinates& coordinates,
    const std::vector<double>& physical
) {
    const auto modified = held2_transform_physical_fractions(coordinates, physical);
    std::vector<double> result;
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

class HeldLocalAccess {
public:
    HeldLocalAccess(
        const ProviderContext& provider,
        double temperature,
        double pressure,
        Held2CompositionDomain domain = Held2CompositionDomain::FiniteSearch
    )
        : provider_(provider),
          temperature_(temperature),
          pressure_(pressure),
          domain_(domain),
          charges_(phase_provider_component_charges(provider)),
          coordinates_(make_held2_coordinates(
              charges_, phase_provider_component_ids(provider)
          )) {}

    [[nodiscard]] const Held2Coordinates& coordinates() const {
        return coordinates_;
    }

    [[nodiscard]] Held2StateEvaluation evaluate(
        const std::vector<double>& composition,
        double log_volume
    ) const {
        const auto physical = held2_lift_independent_fractions(
            coordinates_, composition, domain_
        );
        const auto phase = evaluate_nonreactive_phase(
            provider_, charges_, temperature_, physical, std::exp(log_volume)
        );
        return evaluate_held2_phase_block(
            coordinates_,
            composition,
            log_volume,
            pressure_ / (kGasConstant * temperature_),
            pressure_,
            {phase.value, phase.gradient, phase.hessian, phase.pressure_pa},
            domain_
        );
    }

    [[nodiscard]] std::array<double, 2> volume_bounds(
        const std::vector<double>& composition
    ) const {
        return provider_.evaluate_molar_volume_bounds(
            temperature_,
            held2_lift_independent_fractions(coordinates_, composition, domain_),
            kHeld2PackingFractionMinimum,
            kHeld2PackingFractionMaximum
        );
    }

private:
    const ProviderContext& provider_;
    double temperature_;
    double pressure_;
    Held2CompositionDomain domain_;
    std::vector<double> charges_;
    Held2Coordinates coordinates_;
};

struct DeclaredStart {
    Held2StageIICandidate candidate;
    double phase_fraction = 0.0;
};

std::vector<DeclaredStart> discovery_starts(
    const ProviderContext& provider,
    const FixedTpPhaseProblem& problem,
    const Held2Coordinates& coordinates,
    const std::vector<double>& feed
) {
    const auto component_charges = phase_provider_component_charges(provider);
    if (!phase_topology_is_charged(component_charges)
        && provider.transport().component_count != 2) {
        return {};
    }
    FlashResult discovery;
    try {
        discovery = solve_tp_flash(
            provider,
            {problem.temperature_k, problem.pressure_pa, feed}
        );
    } catch (const std::exception&) {
        return {};
    }
    std::vector<DeclaredStart> result;
    if (const auto* held2 = std::get_if<Held2AlgorithmResult>(
            &discovery.solve
        ); held2 != nullptr && held2->phases.size() == problem.phases.size()) {
        result.reserve(held2->phases.size());
        for (const auto& phase : held2->phases) {
            result.push_back({
                {
                    independent(
                        coordinates,
                        phase.physical_fractions_provider_order
                    ),
                    phase.volume,
                    std::log(phase.volume),
                },
                phase.phase_fraction,
            });
        }
        return result;
    }
    const auto* held = std::get_if<HeldResult>(&discovery.solve);
    if (held == nullptr || problem.phases.size() != 2
        || held->outcome != "accepted" || held->stage_iii_attempts.empty()) {
        return result;
    }
    const auto& local = held->stage_iii_attempts.back().result.local;
    if (!local.accepted) {
        return result;
    }
    for (const auto* phase : {
             &local.evaluation.liquid,
             &local.evaluation.vapor,
         }) {
        const double amount = std::accumulate(
            phase->amounts_mol.begin(), phase->amounts_mol.end(), 0.0
        );
        if (!(amount > 0.0)) {
            return {};
        }
        const std::vector<double> composition{
            phase->amounts_mol[0] / amount,
            phase->amounts_mol[1] / amount,
        };
        const double volume = phase->volume_m3 / amount;
        result.push_back({
            {
                independent(coordinates, composition),
                volume,
                std::log(volume),
            },
            amount,
        });
    }
    return result;
}

void align_starts_with_declared_roles(
    std::vector<DeclaredStart>& starts,
    const std::vector<DeclaredPhaseSpecification>& phases
) {
    const auto vapor = std::find_if(
        phases.begin(), phases.end(), [](const auto& phase) {
            return phase.role == PhaseRole::Vapor;
        }
    );
    if (vapor == phases.end() || starts.size() != phases.size()) {
        return;
    }
    const auto vapor_start = std::max_element(
        starts.begin(), starts.end(), [](const auto& left, const auto& right) {
            return left.candidate.volume < right.candidate.volume;
        }
    );
    const std::size_t vapor_start_index = static_cast<std::size_t>(
        vapor_start - starts.begin()
    );
    const std::size_t vapor_index = static_cast<std::size_t>(
        vapor - phases.begin()
    );
    std::vector<DeclaredStart> aligned;
    aligned.reserve(starts.size());
    std::size_t next_liquid = 0;
    for (std::size_t phase = 0; phase < phases.size(); ++phase) {
        if (phase == vapor_index) {
            aligned.push_back(*vapor_start);
            continue;
        }
        while (next_liquid == vapor_start_index) {
            ++next_liquid;
        }
        aligned.push_back(starts[next_liquid]);
        ++next_liquid;
    }
    starts = std::move(aligned);
}

int reorder_starts_from_hints(
    std::vector<DeclaredStart>& starts,
    const std::vector<PhaseHint>& hints,
    const Held2Coordinates& coordinates
) {
    if (hints.empty() || starts.empty()) {
        return 0;
    }
    std::vector<DeclaredStart> reordered(starts.size());
    std::vector<bool> used(starts.size(), false);
    std::vector<bool> filled(starts.size(), false);
    int reordered_count = 0;
    for (const auto& hint : hints) {
        if (hint.phase_index >= starts.size()
            || hint.mole_fractions.size() != coordinates.charges.size()) {
            throw std::invalid_argument(
                "phase hint identity or composition is invalid"
            );
        }
        std::size_t best = starts.size();
        double best_distance = std::numeric_limits<double>::infinity();
        for (std::size_t candidate = 0; candidate < starts.size(); ++candidate) {
            if (used[candidate]) {
                continue;
            }
            const auto physical = held2_lift_independent_fractions(
                coordinates,
                starts[candidate].candidate.independent_modified_fractions
            );
            double distance = std::abs(std::log(
                (1.0 / starts[candidate].candidate.volume)
                / hint.molar_density_mol_m3
            ));
            for (std::size_t component = 0;
                 component < physical.size(); ++component) {
                distance += std::abs(
                    physical[component] - hint.mole_fractions[component]
                );
            }
            if (distance < best_distance) {
                best = candidate;
                best_distance = distance;
            }
        }
        if (best == starts.size()) {
            throw std::invalid_argument(
                "phase hints cannot be assigned uniquely"
            );
        }
        reordered[hint.phase_index] = starts[best];
        used[best] = true;
        filled[hint.phase_index] = true;
        if (best != hint.phase_index) {
            ++reordered_count;
        }
    }
    std::size_t next = 0;
    for (std::size_t position = 0; position < starts.size(); ++position) {
        if (filled[position]) {
            continue;
        }
        while (used[next]) {
            ++next;
        }
        reordered[position] = starts[next];
        used[next] = true;
    }
    starts = std::move(reordered);
    return reordered_count;
}

std::vector<Held2Problem67Phase> order_phases_from_starts(
    const std::vector<Held2Problem67Phase>& phases,
    const std::vector<DeclaredStart>& starts,
    const Held2Coordinates& coordinates
) {
    if (phases.size() != starts.size()) {
        return phases;
    }
    std::vector<Held2Problem67Phase> ordered(phases.size());
    std::vector<bool> used(phases.size(), false);
    for (std::size_t target = 0; target < starts.size(); ++target) {
        const auto target_composition = held2_lift_independent_fractions(
            coordinates,
            starts[target].candidate.independent_modified_fractions
        );
        std::size_t best = phases.size();
        double best_distance = std::numeric_limits<double>::infinity();
        for (std::size_t candidate = 0;
             candidate < phases.size(); ++candidate) {
            if (used[candidate]) {
                continue;
            }
            double distance = std::abs(std::log(
                phases[candidate].volume
                / starts[target].candidate.volume
            ));
            for (std::size_t component = 0;
                 component < target_composition.size(); ++component) {
                distance += std::abs(
                    phases[candidate].physical_fractions[component]
                    - target_composition[component]
                );
            }
            if (distance < best_distance) {
                best = candidate;
                best_distance = distance;
            }
        }
        ordered[target] = phases[best];
        used[best] = true;
    }
    return ordered;
}

PhaseEquilibriumState all_component_state(
    const ProviderContext& provider,
    PhaseRole role,
    double temperature,
    const std::vector<double>& composition,
    double volume
) {
    const auto component_charges = phase_provider_component_charges(provider);
    const auto phase = evaluate_nonreactive_phase(
        provider, component_charges, temperature, composition, volume
    );
    return make_phase_equilibrium_state(
        provider,
        role,
        PhaseSupport::AllComponents,
        temperature,
        composition,
        volume,
        phase
    );
}

double volume_seed(
    PhaseRole role,
    const std::array<double, 2>& bounds
) {
    const double lower = std::log(bounds[0]);
    const double upper = std::log(bounds[1]);
    if (role == PhaseRole::Liquid) {
        return std::min(bounds[0] * 1.35, std::nextafter(bounds[1], bounds[0]));
    }
    return std::exp(lower + 0.82 * (upper - lower));
}

std::vector<double> perturbed(
    const Held2Coordinates& coordinates,
    const std::vector<double>& base,
    int direction
) {
    if (direction == 0 || base.empty()) {
        return base;
    }
    double delta = 0.08 * static_cast<double>(direction);
    for (int attempt = 0; attempt < 24; ++attempt) {
        std::vector<double> trial = base;
        trial.front() += delta;
        try {
            static_cast<void>(held2_lift_independent_fractions(
                coordinates, trial
            ));
            return trial;
        } catch (const std::exception&) {
            delta *= 0.5;
        }
    }
    return base;
}

bool roles_certified(
    const std::vector<PhaseEquilibriumState>& phases,
    double& minimum_distance
) {
    minimum_distance = std::numeric_limits<double>::infinity();
    for (std::size_t left = 0; left < phases.size(); ++left) {
        if (phases[left].mechanical_class != "strict_stable") {
            return false;
        }
        for (std::size_t right = left + 1; right < phases.size(); ++right) {
            double distance = std::abs(std::log(
                phases[left].molar_volume_m3_per_mol
                / phases[right].molar_volume_m3_per_mol
            ));
            for (std::size_t component = 0;
                 component < phases[left].mole_fractions.size();
                 ++component) {
                distance = std::max(
                    distance,
                    std::abs(
                        phases[left].mole_fractions[component]
                        - phases[right].mole_fractions[component]
                    )
                );
            }
            minimum_distance = std::min(minimum_distance, distance);
            if (phases[left].role == PhaseRole::Liquid
                && phases[right].role == PhaseRole::Vapor
                && phases[left].molar_volume_m3_per_mol
                    >= phases[right].molar_volume_m3_per_mol) {
                return false;
            }
            if (phases[left].role == PhaseRole::Vapor
                && phases[right].role == PhaseRole::Liquid
                && phases[left].molar_volume_m3_per_mol
                    <= phases[right].molar_volume_m3_per_mol) {
                return false;
            }
        }
    }
    return phases.size() == 1 || minimum_distance > kDistinctTolerance;
}

FixedTpPhaseResult one_phase(
    const ProviderContext& provider,
    const FixedTpPhaseProblem& problem,
    const std::vector<double>& feed
) {
    FixedTpPhaseResult result;
    result.parameter_fingerprint = provider.fingerprint();
    result.diagnostics.route = "held_pressure_envelope_declared_topology";
    result.diagnostics.certificate_scope = "fixed_tp_declared_topology";
    result.diagnostics.phase_selection_claim =
        "mechanically_certified_declared_single_phase";
    result.diagnostics.search_scope =
        "fixed_composition_pressure_envelope_48_intervals";
    result.diagnostics.search_budget = "pressure_scan_intervals=48";
    result.diagnostics.search_exclusions = {
        "other_compositions",
        "other_phase_counts",
        "root_completeness_proof",
    };
    result.diagnostics.completeness_criterion =
        "detected_roots_only;pressure_envelope_root_completeness_not_proven";
    result.diagnostics.role_evidence = {
        "native_pressure_derivative_strict_stability",
        "detected_strict_stable_volume_extreme_for_declared_role",
    };
    result.diagnostics.hint_count = static_cast<int>(problem.hints.size());
    const HeldLocalAccess access(
        provider,
        problem.temperature_k,
        problem.pressure_pa,
        Held2CompositionDomain::TraceRefinement
    );
    const auto chart = independent(access.coordinates(), feed);
    const auto envelope = evaluate_held2_pressure_envelope(
        chart,
        access.volume_bounds(chart),
        [&access](const auto& composition, double log_volume) {
            return access.evaluate(composition, log_volume);
        },
        48
    );
    result.diagnostics.provider_evaluation_count = provider.evaluation_count();
    if (envelope.outcome != "selected"
        || envelope.selected_root_index < 0
        || static_cast<std::size_t>(envelope.selected_root_index)
            >= envelope.roots.size()) {
        result.diagnostics.solver_status = "search_exhausted";
        result.diagnostics.numerical_status = "uncertified";
        result.diagnostics.physical_status = "rejected";
        result.diagnostics.failure_stage = "phase_role_certificate";
        result.diagnostics.failure_reason = envelope.failure_reason.empty()
            ? "pressure_envelope_not_certified"
            : envelope.failure_reason;
        return result;
    }
    std::vector<const Held2PressureRoot*> stable;
    for (const auto& root : envelope.roots) {
        if (root.mechanical_class == "strict_stable") {
            stable.push_back(&root);
        }
    }
    std::sort(stable.begin(), stable.end(), [](const auto* left, const auto* right) {
        return left->volume < right->volume;
    });
    if (stable.empty()) {
        result.diagnostics.solver_status = "search_exhausted";
        result.diagnostics.numerical_status = "uncertified";
        result.diagnostics.physical_status = "rejected";
        result.diagnostics.failure_stage = "phase_role_certificate";
        result.diagnostics.failure_reason = "no_detected_strict_stable_root";
        return result;
    }
    const PhaseRole role = problem.phases.front().role;
    const Held2PressureRoot* selected =
        role == PhaseRole::Liquid ? stable.front() : stable.back();
    result.diagnostics.role_evidence.push_back(
        role == PhaseRole::Liquid
            ? "minimum_detected_strict_stable_volume_selected"
            : "maximum_detected_strict_stable_volume_selected"
    );
    if (stable.size() == 1u) {
        result.diagnostics.role_evidence.push_back(
            "single_detected_root_compressibility_role_split_z_0.35"
        );
        const double z = problem.pressure_pa * selected->volume
            / (kGasConstant * problem.temperature_k);
        if ((role == PhaseRole::Liquid && z >= kRoleZSplit)
            || (role == PhaseRole::Vapor && z <= kRoleZSplit)) {
            result.diagnostics.solver_status = "converged";
            result.diagnostics.numerical_status = "certified";
            result.diagnostics.physical_status = "rejected";
            result.diagnostics.failure_stage = "phase_role_certificate";
            result.diagnostics.failure_reason =
                "single_detected_stable_root_role_inconclusive";
            return result;
        }
    }
    result.phases.push_back(all_component_state(
        provider, role, problem.temperature_k, feed, selected->volume
    ));
    const double total = std::accumulate(
        problem.feed_amounts_mol.begin(),
        problem.feed_amounts_mol.end(),
        0.0
    );
    result.phase_amounts_mol = {total};
    result.phase_fractions = {1.0};
    result.total_free_energy_over_rt = total * selected->state.objective;
    result.diagnostics.solver_status = "converged";
    result.diagnostics.numerical_status = "certified";
    result.diagnostics.physical_status = "accepted";
    result.diagnostics.pressure_residual_inf = std::abs(
        result.phases.front().pressure_pa - problem.pressure_pa
    ) / problem.pressure_pa;
    result.diagnostics.phase_distinctness =
        std::numeric_limits<double>::infinity();
    result.diagnostics.provider_evaluation_count = provider.evaluation_count();
    return result;
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
    result.diagnostics.route = "held_step8_declared_topology";
    result.diagnostics.certificate_scope = "fixed_tp_declared_topology";
    result.diagnostics.phase_selection_claim =
        "locally_certified_declared_topology";
    result.diagnostics.search_scope =
        "declared_topology_held_step8_local_equations";
    result.diagnostics.search_budget =
        "problem67_major_iterations=32";
    result.diagnostics.search_exclusions = {
        "undeclared_phase_counts",
        "undeclared_component_supports",
    };
    result.diagnostics.completeness_criterion =
        "all_returned_phases_pass_local_equilibrium_and_role_certificates";
    result.diagnostics.role_evidence = {
        "native_pressure_derivative_strict_stability",
        "ordered_molar_volume_branch",
    };
    result.diagnostics.hint_count = static_cast<int>(problem.hints.size());
    try {
        const auto component_ids = phase_provider_component_ids(provider);
        const auto component_charges = phase_provider_component_charges(provider);
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
        const auto feed = normalize(problem.feed_amounts_mol);
        const auto coordinates = make_held2_coordinates(
            component_charges, component_ids
        );
        static_cast<void>(held2_transform_physical_fractions(
            coordinates, feed
        ));
        if (problem.phases.size() == 1) {
            result = one_phase(provider, problem, feed);
        } else {
            const HeldLocalAccess access(
                provider, problem.temperature_k, problem.pressure_pa
            );
            const auto base = independent(coordinates, feed);
            const std::size_t phase_count = problem.phases.size();
            const std::size_t block_size =
                coordinates.independent_indices.size() + 2;
            std::vector<DeclaredStart> starts = discovery_starts(
                provider, problem, coordinates, feed
            );
            const bool local_candidate_starts =
                starts.size() == phase_count;
            if (!local_candidate_starts) {
                std::vector<std::size_t> neutral_indices;
                for (std::size_t component = 0;
                     component < component_charges.size(); ++component) {
                    if (component_charges[component] == 0.0) {
                        neutral_indices.push_back(component);
                    }
                }
                starts.clear();
                const bool declares_vapor = std::any_of(
                    problem.phases.begin(),
                    problem.phases.end(),
                    [](const auto& phase) {
                        return phase.role == PhaseRole::Vapor;
                    }
                );
                for (std::size_t phase = 0; phase < phase_count; ++phase) {
                    std::vector<double> start;
                    if (phase_count <= neutral_indices.size()) {
                        std::vector<double> enriched = feed;
                        for (double& value : enriched) {
                            value *= 0.15;
                        }
                        enriched[neutral_indices[phase]] += 0.85;
                        start = independent(coordinates, enriched);
                    } else {
                        const int direction = phase_count == 2
                            ? (phase == 0 ? -1 : 1)
                            : static_cast<int>(phase) - 1;
                        start = perturbed(
                            coordinates, base, direction
                        );
                    }
                    const PhaseRole seed_role =
                        declares_vapor && phase == 0
                        ? PhaseRole::Vapor
                        : PhaseRole::Liquid;
                    const double volume = volume_seed(
                        seed_role,
                        access.volume_bounds(start)
                    );
                    starts.push_back({
                        {start, volume, std::log(volume)},
                        1.0 / static_cast<double>(phase_count),
                    });
                }
            }
            align_starts_with_declared_roles(starts, problem.phases);
            const auto canonical_starts = starts;
            result.diagnostics.hint_reordered_start_count =
                reorder_starts_from_hints(
                    starts, problem.hints, coordinates
                );
            std::vector<Held2StageIICandidate> candidates;
            std::vector<std::array<double, 2>> bounds;
            std::vector<double> initial(phase_count * block_size);
            for (std::size_t phase = 0; phase < phase_count; ++phase) {
                const auto& start = starts[phase];
                const auto volume_bounds = access.volume_bounds(
                    start.candidate.independent_modified_fractions
                );
                candidates.push_back(start.candidate);
                bounds.push_back({
                    std::log(volume_bounds[0]),
                    std::log(volume_bounds[1]),
                });
                const std::size_t offset = phase * block_size;
                initial[offset] = start.phase_fraction;
                std::copy(
                    start.candidate.independent_modified_fractions.begin(),
                    start.candidate.independent_modified_fractions.end(),
                    initial.begin() + static_cast<std::ptrdiff_t>(offset + 1)
                );
                initial[offset + block_size - 1] =
                    start.candidate.phase_coordinate;
            }
            const Held2Problem67Result solved = solve_held2_problem67(
                coordinates,
                feed,
                candidates,
                [&access](const auto& composition, double log_volume) {
                    return access.evaluate(composition, log_volume);
                },
                bounds,
                initial,
                32,
                false,
                [&access](const auto& composition) {
                    return access.volume_bounds(composition);
                },
                local_candidate_starts
                    ? kHeld2Problem67ExpandedRadius
                    : 1.0
            );
            result.diagnostics.provider_evaluation_count =
                provider.evaluation_count();
            result.diagnostics.optimizer_solve_count =
                solved.stage_iii_solve_count;
            result.diagnostics.optimizer_iteration_count =
                solved.optimizer_iteration_count;
            result.diagnostics.material_balance_inf =
                solved.ordinary_balance_inf_norm;
            result.diagnostics.electroneutrality_inf =
                solved.phase_charge_inf_norm;
            result.diagnostics.pressure_residual_inf =
                solved.pressure_stationarity_inf_norm;
            result.diagnostics.kkt_stationarity_inf =
                solved.kkt_stationarity_inf_norm;
            result.diagnostics.dual_sign_violation_inf =
                solved.dual_sign_violation_inf_norm;
            result.diagnostics.complementarity_inf =
                solved.bound_complementarity_inf_norm;
            if (solved.numerical_status != "converged"
                || solved.phases.size() != phase_count
                || !audit_held2_tolerance(
                    kHeld2Stage3Charge,
                    solved.phase_charge_inf_norm,
                    solved.phase_charge_scale
                ).passed
                || !audit_held2_tolerance(
                    kHeld2Stage3Stationarity,
                    solved.kkt_stationarity_inf_norm
                ).passed
                || !audit_held2_tolerance(
                    kHeld2Stage3DualSign,
                    solved.dual_sign_violation_inf_norm
                ).passed
                || !audit_held2_tolerance(
                    kHeld2Stage3Complementarity,
                    solved.bound_complementarity_inf_norm
                ).passed) {
                result.diagnostics.solver_status = solved.solver_status;
                result.diagnostics.numerical_status = "not_converged";
                result.diagnostics.physical_status = "rejected";
                result.diagnostics.failure_stage = "local_equilibrium";
                result.diagnostics.failure_reason =
                    solved.numerical_status != "converged"
                    || solved.phases.size() != phase_count
                    ? solved.failure_reason
                    : !audit_held2_tolerance(
                        kHeld2Stage3Charge,
                        solved.phase_charge_inf_norm,
                        solved.phase_charge_scale
                    ).passed
                        ? "stage_iii_electroneutrality_failed"
                        : "stage_iii_kkt_certificate_failed";
                return result;
            }
            std::vector<Held2Problem67Phase> selected =
                order_phases_from_starts(
                    solved.phases, canonical_starts, coordinates
                );
            const double total = std::accumulate(
                problem.feed_amounts_mol.begin(),
                problem.feed_amounts_mol.end(),
                0.0
            );
            for (std::size_t phase = 0; phase < phase_count; ++phase) {
                result.phases.push_back(all_component_state(
                    provider,
                    problem.phases[phase].role,
                    problem.temperature_k,
                    selected[phase].physical_fractions,
                    selected[phase].volume
                ));
                result.phase_fractions.push_back(selected[phase].phase_fraction);
                result.phase_amounts_mol.push_back(
                    total * selected[phase].phase_fraction
                );
            }
            if (!roles_certified(
                    result.phases,
                    result.diagnostics.phase_distinctness
                )) {
                result.diagnostics.provider_evaluation_count =
                    provider.evaluation_count();
                result.phases.clear();
                result.phase_fractions.clear();
                result.phase_amounts_mol.clear();
                result.diagnostics.solver_status = "converged";
                result.diagnostics.numerical_status = "certified";
                result.diagnostics.physical_status = "rejected";
                result.diagnostics.failure_stage = "phase_role_certificate";
                result.diagnostics.failure_reason =
                    "declared_roles_not_thermodynamically_certified";
                return result;
            }
            std::vector<Held2StateEvaluation> evaluations;
            for (const auto& phase : selected) {
                evaluations.push_back(access.evaluate(
                    independent(coordinates, phase.physical_fractions),
                    std::log(phase.volume)
                ));
            }
            bool potentials_certified = true;
            for (std::size_t phase = 1; phase < evaluations.size(); ++phase) {
                for (std::size_t potential = 0;
                     potential < evaluations.front().modified_potentials.size();
                     ++potential) {
                    const double numerator =
                        evaluations[phase].modified_potentials[potential]
                        - evaluations.front().modified_potentials[potential];
                    potentials_certified = potentials_certified
                        && std::abs(numerator)
                            <= kChemicalPotentialResidualMax;
                    result.diagnostics.chemical_potential_residual_inf =
                        std::max(
                            result.diagnostics.chemical_potential_residual_inf,
                            std::abs(numerator)
                        );
                }
            }
            if (!potentials_certified) {
                result.diagnostics.provider_evaluation_count =
                    provider.evaluation_count();
                result.phases.clear();
                result.phase_fractions.clear();
                result.phase_amounts_mol.clear();
                result.diagnostics.solver_status = "converged";
                result.diagnostics.numerical_status = "uncertified";
                result.diagnostics.physical_status = "rejected";
                result.diagnostics.failure_stage = "physical_certificate";
                result.diagnostics.failure_reason =
                    "modified_potential_equality_not_certified";
                return result;
            }
            result.total_free_energy_over_rt = total * solved.objective;
            result.diagnostics.solver_status = solved.solver_status;
            result.diagnostics.numerical_status = "certified";
            result.diagnostics.physical_status = "accepted";
        }
        if (!provider.active_parameters().identities.empty()
            && !result.phases.empty()) {
            std::vector<SolvedPhaseState> solved_phases;
            const double total = std::accumulate(
                problem.feed_amounts_mol.begin(),
                problem.feed_amounts_mol.end(),
                0.0
            );
            for (std::size_t phase = 0; phase < result.phases.size(); ++phase) {
                std::vector<double> amounts;
                for (double fraction : result.phases[phase].mole_fractions) {
                    amounts.push_back(
                        total * result.phase_fractions[phase] * fraction
                    );
                }
                solved_phases.push_back({
                    std::move(amounts),
                    total * result.phase_fractions[phase]
                        * result.phases[phase].molar_volume_m3_per_mol,
                });
            }
            result.sensitivity = evaluate_declared_state_sensitivity(
                provider,
                problem.temperature_k,
                problem.pressure_pa,
                component_ids,
                problem.feed_amounts_mol,
                solved_phases
            );
        }
        result.diagnostics.provider_evaluation_count = provider.evaluation_count();
    } catch (const std::exception& error) {
        result.diagnostics.provider_evaluation_count =
            provider.evaluation_count();
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
