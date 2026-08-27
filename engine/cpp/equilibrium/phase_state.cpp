#include "phase_equilibrium.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

namespace epcsaft_engine::equilibrium {
namespace {

constexpr double kGasConstant = 8.31446261815324;
constexpr double kMechanicalMargin = 1.0e-6;

}  // namespace

std::vector<std::string> phase_provider_component_ids(
    const ProviderContext& provider
) {
    const auto& component_ids = provider.component_ids();
    if (component_ids.empty()) {
        throw std::invalid_argument("phase Provider component identities are incomplete");
    }
    return component_ids;
}

std::vector<double> phase_provider_component_charges(
    const ProviderContext& provider
) {
    const auto& charges = provider.component_charges();
    std::vector<double> result(charges.size(), 0.0);
    for (std::size_t component = 0; component < charges.size(); ++component) {
        result[component] = static_cast<double>(charges[component]);
    }
    return result;
}

bool phase_topology_is_charged(
    const std::vector<double>& charges
) noexcept {
    return std::any_of(charges.begin(), charges.end(), [](double value) {
        return value != 0.0;
    });
}

MixturePhaseEvaluation evaluate_nonreactive_phase(
    const ProviderContext& provider,
    const std::vector<double>& charges,
    double temperature_k,
    const std::vector<double>& mole_fractions,
    double molar_volume_m3_per_mol
) {
    if (mole_fractions.size() == 1 && provider.supports_pure_phase()) {
        const auto pure = provider.evaluate(
            temperature_k, 1.0, molar_volume_m3_per_mol
        );
        MixturePhaseEvaluation result;
        result.value = pure.provider.helmholtz_over_rt_reference_amount;
        result.gradient = {
            pure.provider.gradient[0],
            pure.provider.gradient[1],
        };
        result.hessian.assign(
            pure.provider.hessian.begin(), pure.provider.hessian.end()
        );
        result.pressure_pa = pure.provider.pressure_pa;
        result.parameter_fingerprint = pure.parameter_fingerprint;
        return result;
    }
    return phase_topology_is_charged(charges)
        ? provider.evaluate_electrolyte(
            temperature_k, mole_fractions, molar_volume_m3_per_mol
        )
        : provider.evaluate_mixture(
            temperature_k, mole_fractions, molar_volume_m3_per_mol
        );
}

std::string phase_mechanical_class(
    const MixturePhaseEvaluation& phase,
    double temperature_k,
    double molar_volume_m3_per_mol
) {
    const std::size_t component_count = phase.gradient.size() - 1;
    const double derivative = -kGasConstant * temperature_k
        * molar_volume_m3_per_mol * phase.hessian[
            component_count * (component_count + 1) + component_count
        ];
    const double scale = std::max(1.0, std::abs(phase.pressure_pa));
    if (derivative < -kMechanicalMargin * scale) {
        return "strict_stable";
    }
    if (derivative > kMechanicalMargin * scale) {
        return "unstable";
    }
    return "marginal";
}

PhaseEquilibriumState make_phase_equilibrium_state(
    const ProviderContext& provider,
    PhaseRole role,
    PhaseSupport support,
    double temperature_k,
    const std::vector<double>& mole_fractions,
    double molar_volume_m3_per_mol,
    const MixturePhaseEvaluation& phase
) {
    PhaseEquilibriumState result;
    result.role = role;
    result.support = support;
    result.supported_component_ids = phase_provider_component_ids(provider);
    result.mole_fractions = mole_fractions;
    result.molar_volume_m3_per_mol = molar_volume_m3_per_mol;
    result.molar_density_mol_m3 = 1.0 / molar_volume_m3_per_mol;
    result.packing_fraction = provider.evaluate_packing_fraction(
        temperature_k, mole_fractions, molar_volume_m3_per_mol
    ).value;
    result.pressure_pa = phase.pressure_pa;
    result.chemical_potential_component_ids = result.supported_component_ids;
    result.chemical_potential_over_rt.assign(
        phase.gradient.begin(), phase.gradient.end() - 1
    );
    result.mechanical_class = phase_mechanical_class(
        phase, temperature_k, molar_volume_m3_per_mol
    );
    return result;
}

}  // namespace epcsaft_engine::equilibrium
