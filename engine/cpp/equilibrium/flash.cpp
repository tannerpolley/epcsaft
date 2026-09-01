#include "flash.hpp"

#include "held2_algorithm.hpp"
#include "held2_tolerances.hpp"

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

constexpr double kGasConstantJPerMolK = 8.31446261815324;
void require_mixture_provider(const ProviderContext& provider) {
    if (!provider.supports_mixture_phase()) {
        throw std::invalid_argument("provider model is missing its mixture phase evaluator");
    }
}

void require_held2_provider(const ProviderContext& provider) {
    if (provider.component_count() < 3
        || provider.component_ids().size() != provider.component_count()
        || provider.component_charges().size() != provider.component_count()
        || !provider.supports_electrolyte_phase()) {
        throw std::invalid_argument("provider model is missing the electrolyte phase requirements");
    }
}

void validate_input(const ProviderContext& provider, const FlashInput& input) {
    if (!std::isfinite(input.temperature_k) || input.temperature_k <= 0.0
        || !std::isfinite(input.pressure_pa) || input.pressure_pa <= 0.0
        || input.overall_mole_fractions.size() != provider.component_count()
        || !std::all_of(
            input.overall_mole_fractions.begin(),
            input.overall_mole_fractions.end(),
            [](double value) { return std::isfinite(value) && value > 0.0; }
        )) {
        throw std::invalid_argument(
            "HELD discovery inputs must be positive, finite, and match the Provider model"
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
        const ProviderContext& provider,
        const FlashInput& input
    )
        : provider_(provider),
          input_(input),
          coordinates_(
              make_held2_coordinates(charges_from_provider(provider))
          ),
          total_ion_mole_fraction_max_(
              provider.metadata().ion_mole_fraction_max
          ),
          pressure_over_rt_(
              input.pressure_pa / (kGasConstantJPerMolK * input.temperature_k)
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
        return provider_.evaluate_molar_volume_bounds(
            input_.temperature_k,
            held2_lift_independent_fractions(coordinates_, independent),
            kHeld2PackingFractionMinimum,
            kHeld2PackingFractionMaximum
        );
    }

    [[nodiscard]] double packing_fraction(
        const std::vector<double>& independent,
        double volume
    ) const {
        return provider_.evaluate_packing_fraction(
            input_.temperature_k,
            held2_lift_independent_fractions(coordinates_, independent),
            volume
        ).value;
    }

    [[nodiscard]] std::vector<double> independent(
        const std::vector<double>& physical
    ) const {
        const std::vector<double> modified =
            held2_transform_physical_fractions(coordinates_, physical);
        std::vector<double> result;
        for (std::size_t provider : coordinates_.independent_indices) {
            const auto retained = std::find(
                coordinates_.retained_indices.begin(),
                coordinates_.retained_indices.end(),
                provider
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
        MixturePhaseEvaluation provider_phase =
            provider_.evaluate_electrolyte(
                input_.temperature_k, amounts, std::exp(log_volume)
            );
        Held2PhysicalPhaseBlock block;
        block.helmholtz_over_rt = provider_phase.value;
        block.gradient = std::move(provider_phase.gradient);
        block.hessian = std::move(provider_phase.hessian);
        block.pressure_pa = provider_phase.pressure_pa;
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

    static std::vector<double> charges_from_provider(
        const ProviderContext& provider
    ) {
        const auto& provider_charges = provider.component_charges();
        std::vector<double> charges;
        charges.reserve(provider_charges.size());
        for (int charge : provider_charges) {
            charges.push_back(static_cast<double>(charge));
        }
        return charges;
    }

    const ProviderContext& provider_;
    FlashInput input_;
    Held2Coordinates coordinates_;
    double total_ion_mole_fraction_max_;
    double pressure_over_rt_;
};

std::vector<std::string> component_ids(const ProviderContext& provider) {
    return provider.component_ids();
}

Held2ThermodynamicAccess make_installed_held2_access(
    const ProviderContext& provider,
    const FlashInput& input
) {
    require_held2_provider(provider);
    validate_input(provider, input);
    const auto problem = std::make_shared<InstalledHeld2Problem>(
        provider, input
    );
    Held2ThermodynamicAccess result;
    result.provider_fingerprint = provider.fingerprint();
    result.component_ids = component_ids(provider);
    result.charges.reserve(provider.component_count());
    for (int charge : provider.component_charges()) {
        result.charges.push_back(static_cast<double>(charge));
    }
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
    const ProviderContext& provider,
    FlashResult& result,
    const HeldResult& solve
) {
    if (provider.active_parameters().identities.empty()) {
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
            provider,
            result.input.temperature_k,
            result.input.pressure_pa,
            component_ids(provider),
            result.input.overall_mole_fractions,
            phases
        );
    }
}

void attach_electrolyte_sensitivity(
    const ProviderContext& provider,
    FlashResult& result,
    const Held2AlgorithmResult& solve
) {
    const bool accepted =
        solve.outcome == "one_phase_no_negative_witness_detected"
        || solve.outcome == "physical_equilibrium_accepted";
    if (!accepted || provider.active_parameters().identities.empty()
        || solve.phases.empty()) {
        return;
    }
    std::vector<SolvedPhaseState> phases;
    phases.reserve(solve.phases.size());
    for (const Held2Phase& phase : solve.phases) {
        std::vector<double> amounts;
        amounts.reserve(phase.physical_fractions_provider_order.size());
        for (double fraction : phase.physical_fractions_provider_order) {
            amounts.push_back(phase.phase_fraction * fraction);
        }
        phases.push_back({std::move(amounts), phase.phase_fraction * phase.volume});
    }
    result.sensitivity = evaluate_solved_state_sensitivity(
        provider,
        result.input.temperature_k,
        result.input.pressure_pa,
        component_ids(provider),
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
    const bool electrolyte = std::any_of(
        provider.component_charges().begin(),
        provider.component_charges().end(),
        [](int charge) { return charge != 0; }
    );
    if (electrolyte) {
        require_held2_provider(provider);
        Held2AlgorithmResult solve = run_held2_algorithm(
            make_installed_held2_access(provider, input),
            input,
            {}
        );
        attach_electrolyte_sensitivity(provider, result, solve);
        result.solve = std::move(solve);
        return result;
    }
    require_mixture_provider(provider);
    if (input.overall_mole_fractions.size() != 2
        || provider.component_count() != 2) {
        throw std::invalid_argument(
            "neutral HELD discovery requires a two-component model"
        );
    }
    HeldResult solve = solve_held(
        provider,
        input.temperature_k,
        input.pressure_pa,
        input.overall_mole_fractions.front()
    );
    attach_neutral_sensitivity(provider, result, solve);
    result.solve = std::move(solve);
    return result;
}

}  // namespace epcsaft_engine::equilibrium
