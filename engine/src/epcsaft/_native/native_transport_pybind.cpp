#include "native_transport_pybind.hpp"
#include "native_transport.hpp"

#include <pybind11/stl.h>

#include <cstddef>

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
) {
    const epcsaft_native_transport* transport = create_native_transport(
        model,
        parameter_fingerprint,
        temperature_min_k,
        temperature_max_k,
        component_ids,
        ion_mole_fraction_max,
        pressure_min_pa,
        pressure_max_pa,
        topology_fingerprint,
        component_order_fingerprint,
        domain_fingerprint,
        artifact_fingerprint,
        active_parameters
    );
    return pybind11::capsule(
        const_cast<epcsaft_native_transport*>(transport),
        EPCSAFT_NATIVE_TRANSPORT_CAPSULE_NAME,
        [](PyObject* capsule) {
            void* pointer = PyCapsule_GetPointer(
                capsule, EPCSAFT_NATIVE_TRANSPORT_CAPSULE_NAME
            );
            if (pointer == nullptr) {
                PyErr_Clear();
                return;
            }
            destroy_native_transport(
                static_cast<const epcsaft_native_transport*>(pointer)
            );
        }
    );
}

pybind11::dict evaluate_source_reference_transfer(
    const pybind11::capsule& capsule,
    double temperature_k,
    double pressure_pa,
    const std::vector<std::string>& component_ids,
    const std::vector<double>& source_solvent_mole_fractions,
    const std::vector<std::array<std::size_t, 2>>& ion_pair_indices,
    const std::vector<std::array<int, 2>>& ion_pair_stoichiometry,
    const std::string& phase,
    const std::string& reference_convention,
    const std::string& reference_state_id,
    const std::string& activity_convention_id,
    const std::string& reference_state_fingerprint,
    double standard_molality_mol_per_kg,
    double source_reference_pressure_pa,
    std::uint32_t required_derivatives
) {
    void* pointer = PyCapsule_GetPointer(
        capsule.ptr(), EPCSAFT_NATIVE_TRANSPORT_CAPSULE_NAME
    );
    if (pointer == nullptr) {
        throw pybind11::value_error("source-reference Provider capsule is invalid");
    }
    const auto& transport = *static_cast<const epcsaft_native_transport*>(pointer);
    if (transport.evaluate_source_reference_transfer == nullptr) {
        throw pybind11::value_error(
            "unsupported-source-reference-transfer"
        );
    }
    std::vector<const char*> component_id_pointers;
    component_id_pointers.reserve(component_ids.size());
    for (const auto& component_id : component_ids) {
        component_id_pointers.push_back(component_id.c_str());
    }
    const std::size_t row_count = component_ids.size() - 1u;
    std::vector<double> basis(row_count * component_ids.size());
    std::vector<double> fugacity(row_count);
    std::vector<double> activity(row_count);
    std::vector<double> transfer(row_count);
    std::vector<double> composition(component_ids.size());
    const bool pressure_derivatives_requested = (required_derivatives & 2u) != 0u;
    const bool parameter_derivatives_requested = (required_derivatives & 8u) != 0u;
    std::vector<double> pressure_derivatives(
        pressure_derivatives_requested ? row_count : 0u
    );
    std::vector<double> parameter_derivatives(
        parameter_derivatives_requested
            ? row_count * transport.active_parameter_count
            : 0u
    );
    if (ion_pair_indices.size() != ion_pair_stoichiometry.size()) {
        throw pybind11::value_error(
            "source-reference ion pair and stoichiometry counts differ"
        );
    }
    std::vector<std::size_t> cations;
    std::vector<std::size_t> anions;
    std::vector<std::int32_t> cation_stoichiometry;
    std::vector<std::int32_t> anion_stoichiometry;
    for (std::size_t index = 0; index < ion_pair_indices.size(); ++index) {
        cations.push_back(ion_pair_indices[index][0]);
        anions.push_back(ion_pair_indices[index][1]);
        cation_stoichiometry.push_back(ion_pair_stoichiometry[index][0]);
        anion_stoichiometry.push_back(ion_pair_stoichiometry[index][1]);
    }
    epcsaft_source_reference_request request{
        sizeof(epcsaft_source_reference_request),
        component_ids.size(),
        component_id_pointers.data(),
        source_solvent_mole_fractions.size(),
        source_solvent_mole_fractions.data(),
        ion_pair_indices.size(),
        cations.data(),
        anions.data(),
        cation_stoichiometry.data(),
        anion_stoichiometry.data(),
        phase.c_str(),
        reference_convention.c_str(),
        reference_state_id.c_str(),
        activity_convention_id.c_str(),
        reference_state_fingerprint.c_str(),
        standard_molality_mol_per_kg,
        source_reference_pressure_pa,
        required_derivatives,
    };
    epcsaft_source_reference_transfer_result result{};
    result.struct_size =
        sizeof(epcsaft_source_reference_transfer_result);
    result.component_count = component_ids.size();
    result.neutral_basis_row_count = row_count;
    result.neutral_basis_capacity = basis.size();
    result.contraction_capacity = row_count;
    result.source_composition_capacity = composition.size();
    result.pressure_derivative_capacity = pressure_derivatives.size();
    result.parameter_derivative_capacity = parameter_derivatives.size();
    result.neutral_basis = basis.data();
    result.log_fugacity_contractions = fugacity.data();
    result.activity_scale_log_contractions = activity.data();
    result.transfer_log_contractions = transfer.data();
    result.source_composition = composition.data();
    result.pressure_derivatives_per_pa = pressure_derivatives.data();
    result.parameter_derivatives = parameter_derivatives.data();
    const int status = transport.evaluate_source_reference_transfer(
        transport.model_context,
        temperature_k,
        pressure_pa,
        &request,
        &result
    );
    pybind11::dict value;
    value["status"] = status;
    value["error"] = std::string(result.error);
    if (status != EPCSAFT_SOURCE_REFERENCE_OK) {
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
    value["pressure_derivatives_per_pa"] = pressure_derivatives;
    value["parameter_derivatives"] = parameter_derivatives;
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
