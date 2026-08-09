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
void require_mixture_transport(const epcsaft_native_transport& transport) {
    if (transport.evaluate_mixture_phase == nullptr) {
        throw std::invalid_argument("provider capsule is missing its mixture phase evaluator");
    }
}

void require_held2_transport(const epcsaft_native_transport& transport) {
    if (transport.component_count < 3
        || transport.component_ids == nullptr
        || transport.component_charges == nullptr
        || transport.evaluate_electrolyte_phase == nullptr) {
        throw std::invalid_argument("provider capsule is missing the electrolyte phase contract");
    }
    if (transport.evaluate_molar_volume_bounds == nullptr) {
        throw std::invalid_argument(
            "provider capsule is missing the molar-volume domain contract"
        );
    }
    if (transport.evaluate_packing_fraction == nullptr) {
        throw std::invalid_argument(
            "provider capsule is missing the packing-fraction contract"
        );
    }
}

void validate_input(const ProviderContext& provider, const FlashInput& input) {
    if (!std::isfinite(input.temperature_k) || input.temperature_k <= 0.0
        || !std::isfinite(input.pressure_pa) || input.pressure_pa <= 0.0
        || input.overall_mole_fractions.size() != provider.transport().component_count
        || !std::all_of(
            input.overall_mole_fractions.begin(),
            input.overall_mole_fractions.end(),
            [](double value) { return std::isfinite(value) && value > 0.0; }
        )) {
        throw std::invalid_argument(
            "tp_flash inputs must be positive, finite, and match the Provider model"
        );
    }
    double total = 0.0;
    for (double value : input.overall_mole_fractions) {
        total += value;
    }
    if (std::abs(total - 1.0) > 1.0e-12) {
        throw std::invalid_argument("tp_flash composition must be normalized");
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
              make_held2_coordinates(charges_from_transport(provider.transport()))
          ),
          total_ion_mole_fraction_max_(
              provider.transport().total_ion_mole_fraction_max
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

    static std::vector<double> charges_from_transport(
        const epcsaft_native_transport& transport
    ) {
        if (transport.component_charges == nullptr) {
            throw std::invalid_argument(
                "Provider transport is missing electrolyte component charges"
            );
        }
        std::vector<double> charges;
        charges.reserve(transport.component_count);
        for (std::size_t component = 0;
             component < transport.component_count;
             ++component) {
            charges.push_back(
                static_cast<double>(transport.component_charges[component])
            );
        }
        return charges;
    }

    const ProviderContext& provider_;
    FlashInput input_;
    Held2Coordinates coordinates_;
    double total_ion_mole_fraction_max_;
    double pressure_over_rt_;
};

Held2ThermodynamicAccess make_installed_held2_access(
    const ProviderContext& provider,
    const FlashInput& input
) {
    require_held2_transport(provider.transport());
    validate_input(provider, input);
    const auto problem = std::make_shared<InstalledHeld2Problem>(
        provider, input
    );
    Held2ThermodynamicAccess result;
    result.provider_fingerprint = provider.fingerprint();
    result.component_ids.reserve(provider.transport().component_count);
    result.charges.reserve(provider.transport().component_count);
    for (std::size_t component = 0;
         component < provider.transport().component_count;
         ++component) {
        const char* id = provider.transport().component_ids[component];
        if (id == nullptr || id[0] == '\0') {
            throw std::invalid_argument(
                "Provider electrolyte component ID must not be empty"
            );
        }
        result.component_ids.emplace_back(id);
        result.charges.push_back(static_cast<double>(
            provider.transport().component_charges[component]
        ));
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

}  // namespace

FlashResult solve_tp_flash(
    const ProviderContext& provider,
    const FlashInput& input,
    Held2ProgressObserver* observer
) {
    validate_input(provider, input);
    FlashResult result;
    result.input = input;
    result.parameter_fingerprint = provider.fingerprint();
    const bool electrolyte = provider.transport().component_charges != nullptr
        && std::any_of(
            provider.transport().component_charges,
            provider.transport().component_charges
                + provider.transport().component_count,
            [](int32_t charge) { return charge != 0; }
        );
    if (electrolyte) {
        require_held2_transport(provider.transport());
        result.solve = run_held2_algorithm(
            make_installed_held2_access(provider, input),
            input,
            {},
            observer
        );
        return result;
    }
    require_mixture_transport(provider.transport());
    if (input.overall_mole_fractions.size() != 2
        || provider.transport().component_count != 2) {
        throw std::invalid_argument(
            "neutral tp_flash requires a two-component model"
        );
    }
    result.solve = solve_held(
        provider,
        input.temperature_k,
        input.pressure_pa,
        input.overall_mole_fractions.front()
    );
    return result;
}

}  // namespace epcsaft_engine::equilibrium
