#pragma once

#include "eos.hpp"

#include <pybind11/pybind11.h>

#include <cstdint>
#include <string>
#include <vector>

namespace epcsaft {

pybind11::capsule make_native_transport(
    const NativeModel& model,
    const std::string& parameter_fingerprint,
    double temperature_min_k,
    double temperature_max_k,
    const std::vector<std::string>& component_ids,
    double ion_mole_fraction_max,
    double pressure_min_pa,
    double pressure_max_pa,
    const std::string& topology_fingerprint,
    const std::string& component_order_fingerprint,
    const std::string& domain_fingerprint,
    const std::string& artifact_fingerprint,
    const std::vector<NativeActiveCoordinate>& active_parameters
);

pybind11::dict evaluate_source_reference_transfer(
    const pybind11::capsule& capsule,
    double temperature_k,
    double pressure_pa,
    const std::vector<std::string>& component_ids,
    const std::vector<double>& source_solvent_mole_fractions,
    const std::vector<std::array<std::size_t, 2>>& ion_pair_indices,
    const std::vector<std::array<int, 2>>& ion_pair_stoichiometry,
    const std::vector<std::size_t>& pure_component_indices,
    const std::vector<std::string>& pure_component_phases,
    const std::string& phase,
    const std::string& reference_convention,
    const std::string& reference_state_id,
    const std::string& activity_convention_id,
    const std::string& reference_state_fingerprint,
    double standard_molality_mol_per_kg,
    double source_reference_pressure_pa,
    std::uint32_t required_derivatives
);

}  // namespace epcsaft
