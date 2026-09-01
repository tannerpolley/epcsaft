#include "flash.hpp"

#include "held2_algorithm.hpp"
#include "held2_tolerances.hpp"
#include "phase_equilibrium.hpp"
#include "physical_constants.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <memory>
#include <stdexcept>
#include <utility>

namespace epcsaft_engine::equilibrium {
namespace {

using epcsaft::detail::GAS_CONSTANT;
void require_mixture_eos(const epcsaft::NativeModel& eos) {
    if (!eos.supports_mixture_transport()
        && !eos.supports_fixed_topology_associating_mixture()) {
        throw std::invalid_argument("eos model is missing its mixture phase evaluator");
    }
}

void require_held2_eos(const epcsaft::NativeModel& eos) {
    if (eos.component_count() < 3
        || eos.component_ids().size() != eos.component_count()
        || eos.charge_numbers().size() != eos.component_count()
        || !eos.supports_electrolyte_transport()) {
        throw std::invalid_argument("eos model is missing the electrolyte phase requirements");
    }
}

void validate_input(const epcsaft::NativeModel& eos, const FlashInput& input) {
    if (!std::isfinite(input.temperature_k) || input.temperature_k <= 0.0
        || !std::isfinite(input.pressure_pa) || input.pressure_pa <= 0.0
        || input.overall_mole_fractions.size() != eos.component_count()
        || !std::all_of(
            input.overall_mole_fractions.begin(),
            input.overall_mole_fractions.end(),
            [](double value) { return std::isfinite(value) && value > 0.0; }
        )) {
        throw std::invalid_argument(
            "HELD discovery inputs must be positive, finite, and match the EOS model"
        );
    }
    double total = 0.0;
    for (double value : input.overall_mole_fractions) {
        total += value;
    }
    if (std::abs(total - 1.0) > 1.0e-12) {
        throw std::invalid_argument("HELD discovery composition must be normalized");
    }
}

class InstalledHeld2Problem {
public:
    InstalledHeld2Problem(
        const epcsaft::NativeModel& eos,
        const FlashInput& input
    )
        : eos_(eos),
          input_(input),
          coordinates_(
              make_held2_coordinates(phase_eos_component_charges(eos))
          ),
          total_ion_mole_fraction_max_(
              eos.eos_metadata().ion_mole_fraction_max
          ),
          pressure_over_rt_(
              input.pressure_pa / (GAS_CONSTANT * input.temperature_k)
          ) {
        static_cast<void>(held2_map_unit_cube_to_independent_fractions(
            coordinates_,
            std::vector<double>(coordinates_.independent_indices.size(), 0.0),
            total_ion_mole_fraction_max_
        ));
    }

    [[nodiscard]] Held2StateEvaluation evaluate(
        const std::vector<double>& independent,
        double log_volume
    ) const {
        return evaluate_physical(
            independent,
            held2_lift_independent_fractions(coordinates_, independent),
            log_volume,
            Held2CompositionDomain::FiniteSearch
        );
    }

    [[nodiscard]] Held2StateEvaluation evaluate_trace(
        const std::vector<double>& independent,
        double log_volume
    ) const {
        return evaluate_physical(
            independent,
            held2_lift_independent_fractions(
                coordinates_,
                independent,
                Held2CompositionDomain::TraceRefinement
            ),
            log_volume,
            Held2CompositionDomain::TraceRefinement
        );
    }

    [[nodiscard]] std::array<double, 2> volume_bounds(
        const std::vector<double>& independent
    ) const {
        return eos_.evaluate_molar_volume_bounds(
            input_.temperature_k,
            held2_lift_independent_fractions(coordinates_, independent),
            kHeld2PackingFractionMinimum,
            kHeld2PackingFractionMaximum
        ).bounds_m3_per_mol;
    }

    [[nodiscard]] double packing_fraction(
        const std::vector<double>& independent,
        double volume
    ) const {
        return eos_.evaluate_state_packing_fraction(
            input_.temperature_k,
            held2_lift_independent_fractions(coordinates_, independent),
            volume
        ).packing_fraction;
    }

    [[nodiscard]] std::vector<double> independent(
        const std::vector<double>& physical
    ) const {
        const std::vector<double> modified =
            held2_transform_physical_fractions(coordinates_, physical);
        std::vector<double> result;
        for (std::size_t eos : coordinates_.independent_indices) {
            const auto retained = std::find(
                coordinates_.retained_indices.begin(),
                coordinates_.retained_indices.end(),
                eos
            );
            result.push_back(modified[static_cast<std::size_t>(
                retained - coordinates_.retained_indices.begin()
            )]);
        }
        return result;
    }

    [[nodiscard]] double total_ion_mole_fraction_max() const {
        return total_ion_mole_fraction_max_;
    }

private:
    [[nodiscard]] Held2StateEvaluation evaluate_physical(
        const std::vector<double>& independent,
        const std::vector<double>& amounts,
        double log_volume,
        Held2CompositionDomain domain
    ) const {
        epcsaft::NativePhaseState eos_phase =
            eos_.evaluate_phase_state(
                input_.temperature_k, amounts, std::exp(log_volume)
            );
        Held2PhysicalPhaseBlock block;
        block.helmholtz_over_rt = eos_phase.value;
        block.gradient = std::move(eos_phase.gradient);
        block.hessian = std::move(eos_phase.hessian);
        block.pressure_pa = eos_phase.pressure_pa;
        return evaluate_held2_phase_block(
            coordinates_,
            independent,
            log_volume,
            pressure_over_rt_,
            input_.pressure_pa,
            block,
            domain
        );
    }

    const epcsaft::NativeModel& eos_;
    FlashInput input_;
    Held2Coordinates coordinates_;
    double total_ion_mole_fraction_max_;
    double pressure_over_rt_;
};

std::vector<std::string> component_ids(const epcsaft::NativeModel& eos) {
    return eos.component_ids();
}

Held2ThermodynamicAccess make_installed_held2_access(
    const epcsaft::NativeModel& eos,
    const FlashInput& input
) {
    require_held2_eos(eos);
    validate_input(eos, input);
    const auto problem = std::make_shared<InstalledHeld2Problem>(
        eos, input
    );
    Held2ThermodynamicAccess result;
    result.eos_fingerprint = eos.parameter_fingerprint();
    result.component_ids = component_ids(eos);
    result.charges = phase_eos_component_charges(eos);
    result.evaluate = [problem](const auto& composition, double log_volume) {
        return problem->evaluate(composition, log_volume);
    };
    result.evaluate_trace = [problem](
        const auto& composition, double log_volume
    ) {
        return problem->evaluate_trace(composition, log_volume);
    };
    result.volume_bounds_physical = [problem](const auto& physical) {
        return problem->volume_bounds(problem->independent(physical));
    };
    result.packing_fraction = [problem](
        const auto& composition, double volume
    ) {
        return problem->packing_fraction(composition, volume);
    };
    result.total_ion_mole_fraction_max =
        problem->total_ion_mole_fraction_max();
    return result;
}

void attach_neutral_sensitivity(
    const epcsaft::NativeModel& eos,
    FlashResult& result,
    const HeldResult& solve
) {
    if (eos.active_coordinates().empty()) {
        return;
    }
    std::vector<SolvedPhaseState> phases;
    if (solve.outcome == "one_phase" && solve.stage_i.has_reference) {
        phases.push_back({
            std::vector<double>(
                solve.stage_i.reference.amounts_mol.begin(),
                solve.stage_i.reference.amounts_mol.end()
            ),
            solve.stage_i.reference.volume_m3,
        });
    } else if (solve.outcome == "accepted" && !solve.stage_iii_attempts.empty()) {
        const auto& local = solve.stage_iii_attempts.back().result.local;
        if (local.accepted) {
            for (const auto* phase : {
                     &local.evaluation.liquid,
                     &local.evaluation.vapor,
                 }) {
                phases.push_back({
                    std::vector<double>(
                        phase->amounts_mol.begin(), phase->amounts_mol.end()
                    ),
                    phase->volume_m3,
                });
            }
        }
    }
    if (!phases.empty()) {
        result.sensitivity = evaluate_solved_state_sensitivity(
            eos,
            result.input.temperature_k,
            result.input.pressure_pa,
            component_ids(eos),
            result.input.overall_mole_fractions,
            phases
        );
        result.eos_evaluation_count +=
            result.sensitivity->eos_evaluation_count;
    }
}

void attach_electrolyte_sensitivity(
    const epcsaft::NativeModel& eos,
    FlashResult& result,
    const Held2AlgorithmResult& solve
) {
    const bool accepted =
        solve.outcome == "one_phase_no_negative_witness_detected"
        || solve.outcome == "physical_equilibrium_accepted";
    if (!accepted || eos.active_coordinates().empty()
        || solve.phases.empty()) {
        return;
    }
    std::vector<SolvedPhaseState> phases;
    phases.reserve(solve.phases.size());
    for (const Held2Phase& phase : solve.phases) {
        std::vector<double> amounts;
        amounts.reserve(phase.physical_fractions_eos_order.size());
        for (double fraction : phase.physical_fractions_eos_order) {
            amounts.push_back(phase.phase_fraction * fraction);
        }
        phases.push_back({std::move(amounts), phase.phase_fraction * phase.volume});
    }
    result.sensitivity = evaluate_solved_state_sensitivity(
        eos,
        result.input.temperature_k,
        result.input.pressure_pa,
        component_ids(eos),
        result.input.overall_mole_fractions,
        phases
    );
    result.eos_evaluation_count += result.sensitivity->eos_evaluation_count;
}

}  // namespace

FlashResult solve_tp_flash(
    const epcsaft::NativeModel& eos,
    const FlashInput& input
) {
    validate_input(eos, input);
    FlashResult result;
    result.input = input;
    result.parameter_fingerprint = eos.parameter_fingerprint();
    const bool electrolyte = std::any_of(
        eos.charge_numbers().begin(),
        eos.charge_numbers().end(),
        [](int charge) { return charge != 0; }
    );
    if (electrolyte) {
        require_held2_eos(eos);
        Held2AlgorithmResult solve = run_held2_algorithm(
            make_installed_held2_access(eos, input),
            input,
            {}
        );
        result.eos_evaluation_count = std::accumulate(
            solve.step_timings.begin(),
            solve.step_timings.end(),
            std::size_t{0},
            [](std::size_t count, const Held2StepTiming& timing) {
                return count + timing.eos_evaluations;
            }
        );
        attach_electrolyte_sensitivity(eos, result, solve);
        result.solve = std::move(solve);
        return result;
    }
    require_mixture_eos(eos);
    if (input.overall_mole_fractions.size() != 2
        || eos.component_count() != 2) {
        throw std::invalid_argument(
            "neutral HELD discovery requires a two-component model"
        );
    }
    HeldResult solve = solve_held(
        eos,
        input.temperature_k,
        input.pressure_pa,
        input.overall_mole_fractions.front()
    );
    result.eos_evaluation_count = solve.eos_evaluation_count;
    attach_neutral_sensitivity(eos, result, solve);
    result.solve = std::move(solve);
    return result;
}

}  // namespace epcsaft_engine::equilibrium
