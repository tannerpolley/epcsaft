#include "flash.hpp"

#include "phase_equilibrium.hpp"
#include "phase_stability_tolerances.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <utility>

namespace epcsaft_engine::equilibrium {
namespace {

constexpr double kGasConstantJPerMolK = 8.31446261815324;

void validate_input(const ProviderContext& provider, const FlashInput& input) {
    if (!provider.supports_mixture_phase()
        || !std::isfinite(input.temperature_k) || input.temperature_k <= 0.0
        || !std::isfinite(input.pressure_pa) || input.pressure_pa <= 0.0
        || input.overall_mole_fractions.size() != provider.component_count()
        || !std::all_of(
            input.overall_mole_fractions.begin(),
            input.overall_mole_fractions.end(),
            [](double value) { return std::isfinite(value) && value > 0.0; }
        )) {
        throw std::invalid_argument(
            "phase-stability inputs must be positive, finite, and match the provider"
        );
    }
    const double total = std::accumulate(
        input.overall_mole_fractions.begin(),
        input.overall_mole_fractions.end(),
        0.0
    );
    if (std::abs(total - 1.0) > 1.0e-12) {
        throw std::invalid_argument(
            "phase-stability discovery composition must be normalized"
        );
    }
}

class InstalledPhaseStabilityProblem {
public:
    InstalledPhaseStabilityProblem(
        const ProviderContext& provider,
        const FlashInput& input
    )
        : provider_(provider),
          input_(input),
          charges_(charges_from_provider(provider)),
          coordinates_(make_phase_stability_coordinates(
              charges_, provider.component_ids()
          )),
          total_ion_mole_fraction_max_(
              std::any_of(charges_.begin(), charges_.end(), [](double value) {
                  return value != 0.0;
              })
              ? provider.metadata().ion_mole_fraction_max
              : 1.0
          ),
          pressure_over_rt_(
              input.pressure_pa / (kGasConstantJPerMolK * input.temperature_k)
          ) {
        static_cast<void>(phase_stability_map_unit_cube_to_independent_fractions(
            coordinates_,
            std::vector<double>(coordinates_.independent_indices.size(), 0.0),
            total_ion_mole_fraction_max_
        ));
    }

    [[nodiscard]] PhaseStabilityStateEvaluation evaluate(
        const std::vector<double>& independent,
        double log_volume
    ) const {
        const auto physical = phase_stability_lift_independent_fractions(
            coordinates_, independent
        );
        MixturePhaseEvaluation phase = evaluate_nonreactive_phase(
            provider_,
            charges_,
            input_.temperature_k,
            physical,
            std::exp(log_volume)
        );
        return evaluate_phase_stability_phase_block(
            coordinates_,
            independent,
            log_volume,
            pressure_over_rt_,
            input_.pressure_pa,
            {
                phase.value,
                std::move(phase.gradient),
                std::move(phase.hessian),
                phase.pressure_pa,
            }
        );
    }

    [[nodiscard]] std::array<double, 2> volume_bounds(
        const std::vector<double>& independent
    ) const {
        return provider_.evaluate_molar_volume_bounds(
            input_.temperature_k,
            phase_stability_lift_independent_fractions(coordinates_, independent),
            kPhaseStabilityPackingFractionMinimum,
            kPhaseStabilityPackingFractionMaximum
        );
    }

    [[nodiscard]] double packing_fraction(
        const std::vector<double>& independent,
        double volume
    ) const {
        return provider_.evaluate_packing_fraction(
            input_.temperature_k,
            phase_stability_lift_independent_fractions(coordinates_, independent),
            volume
        ).value;
    }

    [[nodiscard]] std::vector<double> independent(
        const std::vector<double>& physical
    ) const {
        const auto modified = phase_stability_transform_physical_fractions(
            coordinates_, physical
        );
        std::vector<double> result;
        for (std::size_t provider_index : coordinates_.independent_indices) {
            const auto retained = std::find(
                coordinates_.retained_indices.begin(),
                coordinates_.retained_indices.end(),
                provider_index
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
    static std::vector<double> charges_from_provider(
        const ProviderContext& provider
    ) {
        std::vector<double> charges;
        for (int charge : provider.component_charges()) {
            charges.push_back(static_cast<double>(charge));
        }
        return charges;
    }

    const ProviderContext& provider_;
    FlashInput input_;
    std::vector<double> charges_;
    PhaseStabilityCoordinates coordinates_;
    double total_ion_mole_fraction_max_;
    double pressure_over_rt_;
};

PhaseStabilityThermodynamicAccess make_installed_access(
    const ProviderContext& provider,
    const FlashInput& input
) {
    validate_input(provider, input);
    const auto problem = std::make_shared<InstalledPhaseStabilityProblem>(
        provider, input
    );
    PhaseStabilityThermodynamicAccess result;
    result.provider = &provider;
    result.provider_fingerprint = provider.fingerprint();
    result.component_ids = provider.component_ids();
    for (int charge : provider.component_charges()) {
        result.charges.push_back(static_cast<double>(charge));
    }
    result.evaluate = [problem](const auto& composition, double log_volume) {
        return problem->evaluate(composition, log_volume);
    };
    result.volume_bounds_physical = [problem](const auto& physical) {
        return problem->volume_bounds(problem->independent(physical));
    };
    result.packing_fraction = [problem](const auto& composition, double volume) {
        return problem->packing_fraction(composition, volume);
    };
    result.total_ion_mole_fraction_max =
        problem->total_ion_mole_fraction_max();
    return result;
}

void attach_sensitivity(
    const ProviderContext& provider,
    FlashResult& result
) {
    if (provider.active_parameters().identities.empty()
        || result.solve.phases.empty()
        || (result.solve.outcome != "one_phase_no_negative_witness_detected"
            && result.solve.outcome != "physical_equilibrium_accepted")) {
        return;
    }
    std::vector<SolvedPhaseState> phases;
    for (const PhaseStabilityPhase& phase : result.solve.phases) {
        std::vector<double> amounts;
        for (double fraction : phase.physical_fractions_provider_order) {
            amounts.push_back(phase.phase_fraction * fraction);
        }
        phases.push_back({
            std::move(amounts),
            phase.phase_fraction * phase.volume,
        });
    }
    result.sensitivity = evaluate_solved_state_sensitivity(
        provider,
        result.input.temperature_k,
        result.input.pressure_pa,
        provider.component_ids(),
        result.input.overall_mole_fractions,
        phases
    );
}

}  // namespace

FlashResult solve_tp_flash(
    const ProviderContext& provider,
    const FlashInput& input
) {
    validate_input(provider, input);
    FlashResult result;
    result.input = input;
    result.parameter_fingerprint = provider.fingerprint();
    result.solve = run_phase_stability(
        make_installed_access(provider, input), input, {}
    );
    attach_sensitivity(provider, result);
    return result;
}

}  // namespace epcsaft_engine::equilibrium
