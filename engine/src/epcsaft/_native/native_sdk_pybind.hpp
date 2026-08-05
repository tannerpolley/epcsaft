#pragma once

#include "eos.hpp"

#include <pybind11/pybind11.h>

#include <cstdint>
#include <string>
#include <vector>

namespace epcsaft {

pybind11::capsule make_native_sdk_v1(
    const NativeModel& model,
    const std::string& parameter_fingerprint,
    double temperature_min_k,
    double temperature_max_k,
    const std::vector<std::string>& component_ids,
    double ion_mole_fraction_max,
    bool phase_callbacks_admitted,
    bool active_born_admitted,
    bool aqueous_kij_admitted,
    bool ion_fraction_suppression_admitted,
    double pressure_min_pa,
    double pressure_max_pa,
    bool neutral_reference_admitted,
    bool unknown_applicability_admitted,
    bool aqueous_solvation_factor_admitted,
    bool ion_solvation_kij_admitted,
    const std::string& topology_fingerprint,
    const std::string& component_order_fingerprint,
    const std::string& domain_fingerprint,
    const std::vector<std::string>& site_ids,
    const std::string& artifact_fingerprint
);

pybind11::dict evaluate_source_reference_transfer_v1(
    const pybind11::capsule& capsule,
    const std::string& expected_parameter_fingerprint,
    const std::string& expected_topology_fingerprint,
    double temperature_k,
    double pressure_pa,
    const std::vector<std::string>& component_ids,
    const std::vector<double>& source_solvent_mole_fractions,
    const std::string& reference_state_id,
    const std::string& activity_convention_id,
    const std::string& reference_state_fingerprint,
    double standard_molality_mol_per_kg,
    double source_reference_pressure_pa,
    std::uint32_t required_derivatives
);

}  // namespace epcsaft
