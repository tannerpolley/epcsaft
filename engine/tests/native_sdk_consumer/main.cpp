#include <epcsaft/native_model_v1.h>
#include <epcsaft/native_sdk_v1.h>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <vector>

int main(int argc, char** argv) {
    if (argc != 2) {
        return 2;
    }
    char error[EPCSAFT_NATIVE_SDK_V1_ERROR_SIZE] = {};
    epcsaft_native_model_handle_v1* model =
        epcsaft_native_model_load_v1(argv[1], error, sizeof(error));
    if (model == nullptr) {
        std::cerr << error << '\n';
        return 3;
    }
    const epcsaft_native_sdk_v1* sdk = epcsaft_native_model_sdk_v1(model);
    const char* fingerprint = epcsaft_native_model_fingerprint_v1(model);
    bool valid = sdk != nullptr
        && fingerprint != nullptr
        && std::strncmp(fingerprint, "sha256:", 7) == 0
        && sdk->abi_version == EPCSAFT_NATIVE_SDK_V1_ABI_VERSION
        && sdk->component_count == 3
        && sdk->table_size
            >= offsetof(
                   epcsaft_native_sdk_v1,
                   evaluate_reacting_phase_parameters
               )
                + sizeof(sdk->evaluate_reacting_phase_parameters)
        && sdk->neutral_reference_derivative_result_size
            == sizeof(epcsaft_neutral_reference_derivative_result_v1)
        && sdk->reacting_phase_parameter_result_size
            == sizeof(epcsaft_reacting_phase_parameter_result_v1)
        && sdk->evaluate_neutral_reference_derivatives != nullptr
        && sdk->evaluate_reacting_phase_parameters != nullptr;
    const epcsaft_native_capability_descriptor_v1* diameter = nullptr;
    const epcsaft_native_capability_descriptor_v1* dispersion = nullptr;
    if (valid) {
        for (std::size_t index = 0; index < sdk->capability_count; ++index) {
            const auto& descriptor = sdk->capabilities[index];
            if (descriptor.coordinate_count != 1
                || descriptor.coordinates[0].component_index != 1) {
                continue;
            }
            if (descriptor.capability
                == EPCSAFT_NATIVE_CAPABILITY_REACTING_COMPONENT_DIAMETER_HELMHOLTZ_V1) {
                diameter = &descriptor;
            } else if (
                descriptor.capability
                == EPCSAFT_NATIVE_CAPABILITY_REACTING_COMPONENT_DISPERSION_HELMHOLTZ_V1
            ) {
                dispersion = &descriptor;
            }
        }
    }
    valid = valid && diameter != nullptr && dispersion != nullptr
        && diameter->parameter_family
            == EPCSAFT_NATIVE_PARAMETER_FAMILY_SEGMENT_DIAMETER_V1
        && dispersion->parameter_family
            == EPCSAFT_NATIVE_PARAMETER_FAMILY_DISPERSION_ENERGY_OVER_K_V1
        && std::strcmp(diameter->coordinates[0].unit, "angstrom") == 0
        && std::strcmp(dispersion->coordinates[0].unit, "kelvin") == 0
        && std::strcmp(
            diameter->topology_fingerprint,
            dispersion->topology_fingerprint
        ) == 0;

    epcsaft_active_parameter_request_v1 request{
        sizeof(request),
        EPCSAFT_NATIVE_PARAMETER_FAMILY_SEGMENT_DIAMETER_V1,
        EPCSAFT_NATIVE_PARAMETER_IDENTITY_COMPONENT_V1,
        1,
        -1,
        -1,
        2.8232,
        "angstrom",
    };
    const double amounts[] = {0.98, 0.01, 0.01};
    std::vector<double> state_gradient(4);
    std::vector<double> state_hessian(16);
    std::vector<double> parameter_gradient(1);
    std::vector<double> mixed(4);
    std::vector<double> packing_parameter(1);
    std::vector<double> volume_bound_parameter(2);
    std::vector<double> pressure_parameter(1);
    std::vector<double> chemical_parameter(3);
    epcsaft_reacting_phase_parameter_result_v1 phase{
        sizeof(phase),
        EPCSAFT_NATIVE_STATUS_OK_V1,
        4,
        1,
        state_gradient.size(),
        state_hessian.size(),
        parameter_gradient.size(),
        mixed.size(),
        packing_parameter.size(),
        volume_bound_parameter.size(),
        pressure_parameter.size(),
        chemical_parameter.size(),
        state_gradient.data(),
        state_hessian.data(),
        parameter_gradient.data(),
        mixed.data(),
        packing_parameter.data(),
        volume_bound_parameter.data(),
        pressure_parameter.data(),
        chemical_parameter.data(),
    };
    const int phase_status = valid
        ? sdk->evaluate_reacting_phase_parameters(
            sdk->model_context,
            fingerprint,
            diameter->topology_fingerprint,
            298.15,
            amounts,
            3,
            2.5e-5,
            0.01,
            0.5,
            &request,
            1,
            &phase
        )
        : EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT_V1;

    std::vector<double> basis(6);
    std::vector<double> contractions(2);
    std::vector<double> composition(3);
    std::vector<double> pressure_derivatives(2);
    std::vector<double> parameter_derivatives(2);
    epcsaft_neutral_reference_derivative_result_v1 reference{};
    reference.struct_size = sizeof(reference);
    reference.value.struct_size = sizeof(reference.value);
    reference.value.component_count = 3;
    reference.value.neutral_basis_row_count = 2;
    reference.value.neutral_basis_capacity = basis.size();
    reference.value.contraction_capacity = contractions.size();
    reference.value.reference_composition_capacity = composition.size();
    reference.value.neutral_basis = basis.data();
    reference.value.log_fugacity_contractions = contractions.data();
    reference.value.reference_composition = composition.data();
    reference.active_parameter_count = 1;
    reference.pressure_derivative_capacity = pressure_derivatives.size();
    reference.parameter_derivative_capacity = parameter_derivatives.size();
    reference.pressure_derivatives_per_pa = pressure_derivatives.data();
    reference.parameter_derivatives = parameter_derivatives.data();
    const int reference_status = valid
        ? sdk->evaluate_neutral_reference_derivatives(
            sdk->model_context,
            fingerprint,
            diameter->topology_fingerprint,
            298.15,
            100000.0,
            &request,
            1,
            &reference
        )
        : EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT_V1;
    valid = valid
        && phase_status == EPCSAFT_NATIVE_STATUS_OK_V1
        && reference_status == EPCSAFT_NATIVE_STATUS_OK_V1
        && phase.status == EPCSAFT_NATIVE_STATUS_OK_V1
        && reference.status == EPCSAFT_NATIVE_STATUS_OK_V1
        && reference.derivative_availability
            == (EPCSAFT_NEUTRAL_REFERENCE_DERIVATIVE_PRESSURE_V1
                | EPCSAFT_NEUTRAL_REFERENCE_DERIVATIVE_PARAMETERS_V1)
        && reference.value.derivative_availability
            == EPCSAFT_NEUTRAL_REFERENCE_DERIVATIVE_NONE_V1
        && std::strcmp(reference.reference_branch, "liquid") == 0
        && std::isfinite(phase.pressure_pa)
        && std::isfinite(mixed[0])
        && std::isfinite(pressure_derivatives[0])
        && std::isfinite(parameter_derivatives[0]);
    std::cout << "component_count=" << (sdk == nullptr ? 0 : sdk->component_count)
              << '\n'
              << "phase_pressure_pa=" << phase.pressure_pa << '\n'
              << "reference_branch=" << reference.reference_branch << '\n';
    epcsaft_native_model_destroy_v1(model);
    return valid ? 0 : 4;
}
