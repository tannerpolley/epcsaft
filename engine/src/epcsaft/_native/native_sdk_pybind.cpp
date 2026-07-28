#include "native_sdk_pybind.hpp"

#include "native_sdk.hpp"

#include <epcsaft/native_sdk_v1.h>

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
    double pressure_min_pa,
    double pressure_max_pa,
    bool neutral_reference_admitted,
    bool unknown_applicability_admitted,
    bool aqueous_solvation_factor_admitted,
    bool staged_trial_bundle_admitted,
    const std::string& topology_fingerprint
) {
    const epcsaft_native_sdk_v1* sdk = create_native_sdk_v1(
        model,
        parameter_fingerprint,
        temperature_min_k,
        temperature_max_k,
        component_ids,
        ion_mole_fraction_max,
        phase_callbacks_admitted,
        active_born_admitted,
        aqueous_kij_admitted,
        pressure_min_pa,
        pressure_max_pa,
        neutral_reference_admitted,
        unknown_applicability_admitted,
        aqueous_solvation_factor_admitted,
        staged_trial_bundle_admitted,
        topology_fingerprint
    );
    return pybind11::capsule(
        const_cast<epcsaft_native_sdk_v1*>(sdk),
        EPCSAFT_NATIVE_SDK_V1_CAPSULE_NAME,
        [](PyObject* capsule) {
            void* pointer = PyCapsule_GetPointer(
                capsule, EPCSAFT_NATIVE_SDK_V1_CAPSULE_NAME
            );
            if (pointer == nullptr) {
                PyErr_Clear();
                return;
            }
            destroy_native_sdk_v1(
                static_cast<const epcsaft_native_sdk_v1*>(pointer)
            );
        }
    );
}

}  // namespace epcsaft
