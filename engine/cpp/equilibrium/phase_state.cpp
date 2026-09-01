#include "phase_equilibrium.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

namespace epcsaft_engine::equilibrium {

std::vector<std::string> phase_eos_component_ids(
    const epcsaft::NativeModel& eos
) {
    const auto& component_ids = eos.component_ids();
    if (component_ids.empty()) {
        throw std::invalid_argument("phase EOS component identities are incomplete");
    }
    return component_ids;
}

std::vector<double> phase_eos_component_charges(
    const epcsaft::NativeModel& eos
) {
    const auto& charges = eos.charge_numbers();
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

PhaseEquilibriumState make_phase_equilibrium_state(
    const epcsaft::NativeModel& eos,
    PhaseRole role,
    PhaseSupport support,
    const std::vector<double>& mole_fractions,
    double molar_volume_m3_per_mol,
    const epcsaft::NativePhaseState& phase
) {
    PhaseEquilibriumState result;
    result.role = role;
    result.support = support;
    const auto component_ids = phase_eos_component_ids(eos);
    result.supported_component_ids.reserve(
        phase.supported_component_indices.size()
    );
    result.mole_fractions.reserve(phase.supported_component_indices.size());
    for (std::size_t component : phase.supported_component_indices) {
        result.supported_component_ids.push_back(component_ids[component]);
        result.mole_fractions.push_back(mole_fractions[component]);
    }
    result.molar_volume_m3_per_mol = molar_volume_m3_per_mol;
    result.molar_density_mol_m3 = 1.0 / molar_volume_m3_per_mol;
    result.packing_fraction = phase.packing_fraction;
    result.pressure_pa = phase.pressure_pa;
    result.chemical_potential_component_ids = result.supported_component_ids;
    result.chemical_potential_over_rt.assign(
        phase.gradient.begin(), phase.gradient.end() - 1
    );
    result.mechanical_class = phase.mechanical_class;
    return result;
}

}  // namespace epcsaft_engine::equilibrium
