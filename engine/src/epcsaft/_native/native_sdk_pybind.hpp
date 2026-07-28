#pragma once

#include "eos.hpp"

#include <pybind11/pybind11.h>

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
    bool staged_trial_bundle_admitted,
    const std::string& topology_fingerprint
);

}  // namespace epcsaft
