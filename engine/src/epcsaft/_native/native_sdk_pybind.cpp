#include "native_sdk_pybind.hpp"

#include "native_sdk.hpp"

#include <epcsaft/native_sdk_v1.h>

#include <pybind11/stl.h>

#include <cstddef>

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
        ion_fraction_suppression_admitted,
        pressure_min_pa,
        pressure_max_pa,
        neutral_reference_admitted,
        unknown_applicability_admitted,
        aqueous_solvation_factor_admitted,
        ion_solvation_kij_admitted,
        topology_fingerprint,
        component_order_fingerprint,
        domain_fingerprint,
        site_ids,
        artifact_fingerprint
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
) {
    void* pointer = PyCapsule_GetPointer(
        capsule.ptr(), EPCSAFT_NATIVE_SDK_V1_CAPSULE_NAME
    );
    if (pointer == nullptr) {
        throw pybind11::value_error("source-reference Provider capsule is invalid");
    }
    const auto& sdk = *static_cast<const epcsaft_native_sdk_v1*>(pointer);
    const std::size_t required_table_size =
        offsetof(epcsaft_native_sdk_v1, evaluate_source_reference_transfer)
        + sizeof(sdk.evaluate_source_reference_transfer);
    if (sdk.table_size < required_table_size
        || sdk.source_reference_transfer_result_size
            != sizeof(epcsaft_source_reference_transfer_result_v1)
        || sdk.evaluate_source_reference_transfer == nullptr) {
        throw pybind11::value_error(
            "unsupported-source-reference-transfer"
        );
    }
    std::vector<const char*> component_id_pointers;
    component_id_pointers.reserve(component_ids.size());
    for (const auto& component_id : component_ids) {
        component_id_pointers.push_back(component_id.c_str());
    }
    const std::size_t row_count = sdk.neutral_reference_basis_row_count;
    std::vector<double> basis(row_count * component_ids.size());
    std::vector<double> fugacity(row_count);
    std::vector<double> activity(row_count);
    std::vector<double> transfer(row_count);
    std::vector<double> composition(component_ids.size());
    epcsaft_source_reference_request_v1 request{
        sizeof(epcsaft_source_reference_request_v1),
        component_ids.size(),
        component_id_pointers.data(),
        source_solvent_mole_fractions.size(),
        source_solvent_mole_fractions.data(),
        reference_state_id.c_str(),
        activity_convention_id.c_str(),
        reference_state_fingerprint.c_str(),
        standard_molality_mol_per_kg,
        source_reference_pressure_pa,
        required_derivatives,
    };
    epcsaft_source_reference_transfer_result_v1 result{};
    result.struct_size =
        sizeof(epcsaft_source_reference_transfer_result_v1);
    result.component_count = component_ids.size();
    result.neutral_basis_row_count = row_count;
    result.neutral_basis_capacity = basis.size();
    result.contraction_capacity = row_count;
    result.source_composition_capacity = composition.size();
    result.pressure_derivative_capacity = 0;
    result.neutral_basis = basis.data();
    result.log_fugacity_contractions = fugacity.data();
    result.activity_scale_log_contractions = activity.data();
    result.transfer_log_contractions = transfer.data();
    result.source_composition = composition.data();
    result.pressure_derivatives_per_pa = nullptr;
    const int status = sdk.evaluate_source_reference_transfer(
        sdk.model_context,
        expected_parameter_fingerprint.c_str(),
        expected_topology_fingerprint.c_str(),
        temperature_k,
        pressure_pa,
        &request,
        &result
    );
    pybind11::dict value;
    value["status"] = status;
    value["error"] = std::string(result.error);
    if (status != EPCSAFT_SOURCE_REFERENCE_OK_V1) {
        return value;
    }
    value["component_count"] = result.component_count;
    value["neutral_basis_row_count"] = result.neutral_basis_row_count;
    value["neutral_basis"] = basis;
    value["log_fugacity_contractions"] = fugacity;
    value["activity_scale_log_contractions"] = activity;
    value["transfer_log_contractions"] = transfer;
    value["source_composition"] = composition;
    value["derivative_availability"] = result.derivative_availability;
    value["temperature_k"] = result.temperature_k;
    value["pressure_pa"] = result.pressure_pa;
    value["source_reference_pressure_pa"] =
        result.source_reference_pressure_pa;
    value["standard_molality_mol_per_kg"] =
        result.standard_molality_mol_per_kg;
    value["solvent_molar_mass_kg_per_mol"] =
        result.solvent_molar_mass_kg_per_mol;
    value["reference_limit_molality_mol_per_kg"] =
        result.reference_limit_molality_mol_per_kg;
    value["reference_convergence_error"] =
        result.reference_convergence_error;
    value["source_temperature_min_k"] = result.source_temperature_min_k;
    value["source_temperature_max_k"] = result.source_temperature_max_k;
    value["source_pressure_min_pa"] = result.source_pressure_min_pa;
    value["source_pressure_max_pa"] = result.source_pressure_max_pa;
    value["parameter_fingerprint"] =
        std::string(result.parameter_fingerprint);
    value["topology_fingerprint"] =
        std::string(result.topology_fingerprint);
    value["component_order_fingerprint"] =
        std::string(result.component_order_fingerprint);
    value["reference_state_fingerprint"] =
        std::string(result.reference_state_fingerprint);
    value["domain_fingerprint"] =
        std::string(result.domain_fingerprint);
    value["artifact_fingerprint"] =
        std::string(result.artifact_fingerprint);
    value["helmholtz_basis_id"] =
        std::string(result.helmholtz_basis_id);
    return value;
}

}  // namespace epcsaft
