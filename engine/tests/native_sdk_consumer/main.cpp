#include <epcsaft/native_model_v1.h>
#include <epcsaft/native_sdk_v1.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>

bool neutral_inverse_packing_is_available(const char* path) {
    char error[EPCSAFT_NATIVE_SDK_V1_ERROR_SIZE] = {};
    epcsaft_native_model_handle_v1* model =
        epcsaft_native_model_load_v1(path, error, sizeof(error));
    if (model == nullptr) {
        std::cerr << error << '\n';
        return false;
    }
    const epcsaft_native_sdk_v1* sdk = epcsaft_native_model_sdk_v1(model);
    const char* fingerprint = epcsaft_native_model_fingerprint_v1(model);
    const char* topology = nullptr;
    if (sdk != nullptr) {
        for (std::size_t index = 0; index < sdk->capability_count; ++index) {
            if (sdk->capabilities[index].topology_fingerprint[0] != '\0') {
                topology = sdk->capabilities[index].topology_fingerprint;
                break;
            }
        }
    }
    const double amounts[] = {0.4, 0.6};
    std::vector<double> gradient(3);
    std::vector<double> hessian(9);
    epcsaft_inverse_packing_geometry_result_v1 result{};
    result.struct_size = sizeof(result);
    result.coordinate_count = 3;
    result.gradient_capacity = gradient.size();
    result.hessian_capacity = hessian.size();
    result.gradient = gradient.data();
    result.hessian = hessian.data();
    const int status =
        sdk != nullptr && fingerprint != nullptr && topology != nullptr
        ? sdk->evaluate_inverse_packing_geometry(
            sdk->model_context,
            fingerprint,
            topology,
            298.15,
            amounts,
            2,
            std::log(0.2),
            nullptr,
            0,
            &result
        )
        : EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT_V1;
    double packing_fraction = 0.0;
    std::vector<double> packing_gradient(3);
    std::vector<double> packing_hessian(9);
    const int packing_status = status == EPCSAFT_NATIVE_STATUS_OK_V1
        ? sdk->evaluate_packing_fraction(
            sdk->model_context,
            fingerprint,
            298.15,
            amounts,
            2,
            result.volume_m3,
            &packing_fraction,
            packing_gradient.data(),
            packing_gradient.size(),
            packing_hessian.data(),
            packing_hessian.size()
        )
        : EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT_V1;
    const bool valid = status == EPCSAFT_NATIVE_STATUS_OK_V1
        && packing_status == EPCSAFT_NATIVE_STATUS_OK_V1
        && result.active_parameter_count == 0
        && std::isfinite(result.volume_m3)
        && result.volume_m3 > 0.0
        && std::abs(packing_fraction - 0.2) <= std::numeric_limits<double>::epsilon()
        && std::abs(gradient[2] + result.volume_m3)
            <= 1.0e-14 * result.volume_m3
        && std::abs(hessian[8] - result.volume_m3)
            <= 1.0e-14 * result.volume_m3
        && std::abs(amounts[0] * gradient[0] + amounts[1] * gradient[1]
                    - result.volume_m3)
            <= 3.0e-14 * result.volume_m3;
    epcsaft_native_model_destroy_v1(model);
    return valid;
}

bool fixed_topology_association_is_available(
    const char* path,
    const char* expected_artifact_fingerprint
) {
    char error[EPCSAFT_NATIVE_SDK_V1_ERROR_SIZE] = {};
    epcsaft_native_model_handle_v1* model =
        epcsaft_native_model_load_v1(path, error, sizeof(error));
    if (model == nullptr) {
        std::cerr << error << '\n';
        return false;
    }
    const epcsaft_native_sdk_v1* sdk = epcsaft_native_model_sdk_v1(model);
    const char* fingerprint = epcsaft_native_model_fingerprint_v1(model);
    const auto* descriptor =
        sdk == nullptr ? nullptr : sdk->fixed_topology_association;
    bool valid = sdk != nullptr && fingerprint != nullptr
        && sdk->table_size == sizeof(epcsaft_native_sdk_v1)
        && sdk->fixed_topology_association_result_size
            == sizeof(epcsaft_fixed_topology_phase_result_v1)
        && descriptor != nullptr
        && sdk->evaluate_fixed_topology_association != nullptr
        && descriptor->schema_version
            == EPCSAFT_FIXED_TOPOLOGY_ASSOCIATION_SCHEMA_VERSION_V1
        && descriptor->site_count == 3 && descriptor->slot_count == 7
        && std::strcmp(
            descriptor->artifact_fingerprint,
            expected_artifact_fingerprint
        ) == 0;
    const double values[] = {
        1.4, 3.2, 220.0, 2400.0, 0.03, 1800.0, 0.025,
    };
    std::vector<epcsaft_fixed_topology_parameter_request_v1> requests;
    if (valid) {
        requests.reserve(descriptor->slot_count);
        for (std::size_t index = 0; index < descriptor->slot_count; ++index) {
            const auto& slot = descriptor->slots[index];
            requests.push_back({
                sizeof(epcsaft_fixed_topology_parameter_request_v1),
                slot.parameter_family,
                slot.component_index,
                slot.site_index_a,
                slot.site_index_b,
                slot.unit,
                values[index],
            });
        }
    }
    const std::size_t coordinate_count = 2 + requests.size();
    std::vector<double> gradient(coordinate_count);
    std::vector<double> hessian(coordinate_count * coordinate_count);
    epcsaft_fixed_topology_phase_result_v1 result{};
    result.struct_size = sizeof(result);
    result.coordinate_count = coordinate_count;
    result.active_parameter_count = requests.size();
    result.gradient_capacity = gradient.size();
    result.hessian_capacity = hessian.size();
    result.gradient = gradient.data();
    result.hessian = hessian.data();
    const int status = valid
        ? sdk->evaluate_fixed_topology_association(
            sdk->model_context,
            fingerprint,
            descriptor->topology_fingerprint,
            expected_artifact_fingerprint,
            300.0,
            1.0,
            1.0 / 6000.0,
            requests.data(),
            requests.size(),
            &result
        )
        : EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT_V1;
    valid = valid && status == EPCSAFT_NATIVE_STATUS_OK_V1
        && result.status == EPCSAFT_NATIVE_STATUS_OK_V1
        && std::strcmp(
            result.artifact_fingerprint,
            descriptor->artifact_fingerprint
        ) == 0
        && std::all_of(
            gradient.begin(), gradient.end(),
            [](double value) { return std::isfinite(value); }
        )
        && std::all_of(
            hessian.begin(), hessian.end(),
            [](double value) { return std::isfinite(value); }
        );
    if (valid) {
        std::cout << std::setprecision(17)
                  << "association_value="
                  << result.helmholtz_over_rt_reference_amount << '\n';
        const auto print_vector = [](const char* label, const auto& values) {
            std::cout << label << '=';
            for (std::size_t index = 0; index < values.size(); ++index) {
                if (index != 0) {
                    std::cout << ',';
                }
                std::cout << values[index];
            }
            std::cout << '\n';
        };
        print_vector("association_gradient", gradient);
        print_vector("association_hessian", hessian);
    }
    epcsaft_native_model_destroy_v1(model);
    return valid;
}

int main(int argc, char** argv) {
    if (argc != 6) {
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
                   evaluate_inverse_packing_geometry
               )
                + sizeof(sdk->evaluate_inverse_packing_geometry)
        && sdk->neutral_reference_derivative_result_size
            == sizeof(epcsaft_neutral_reference_derivative_result_v1)
        && sdk->reacting_phase_parameter_result_size
            == sizeof(epcsaft_reacting_phase_parameter_result_v1)
        && sdk->inverse_packing_geometry_result_size
            == sizeof(epcsaft_inverse_packing_geometry_result_v1)
        && sdk->evaluate_neutral_reference_derivatives != nullptr
        && sdk->evaluate_reacting_phase_parameters != nullptr
        && sdk->evaluate_inverse_packing_geometry != nullptr;
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
    std::vector<double> packing_state_gradient(4);
    std::vector<double> packing_state_hessian(16);
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
    phase.packing_state_gradient_capacity = packing_state_gradient.size();
    phase.packing_state_hessian_capacity = packing_state_hessian.size();
    phase.packing_state_gradient = packing_state_gradient.data();
    phase.packing_state_hessian = packing_state_hessian.data();
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
    std::vector<double> inverse_gradient(5);
    std::vector<double> inverse_hessian(25);
    epcsaft_inverse_packing_geometry_result_v1 inverse{};
    inverse.struct_size = sizeof(inverse);
    inverse.coordinate_count = 5;
    inverse.active_parameter_count = 1;
    inverse.gradient_capacity = inverse_gradient.size();
    inverse.hessian_capacity = inverse_hessian.size();
    inverse.gradient = inverse_gradient.data();
    inverse.hessian = inverse_hessian.data();
    const int inverse_status = valid
        ? sdk->evaluate_inverse_packing_geometry(
            sdk->model_context,
            fingerprint,
            diameter->topology_fingerprint,
            298.15,
            amounts,
            3,
            std::log(0.25),
            &request,
            1,
            &inverse
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
        && inverse_status == EPCSAFT_NATIVE_STATUS_OK_V1
        && reference_status == EPCSAFT_NATIVE_STATUS_OK_V1
        && phase.status == EPCSAFT_NATIVE_STATUS_OK_V1
        && inverse.status == EPCSAFT_NATIVE_STATUS_OK_V1
        && reference.status == EPCSAFT_NATIVE_STATUS_OK_V1
        && reference.derivative_availability
            == (EPCSAFT_NEUTRAL_REFERENCE_DERIVATIVE_PRESSURE_V1
                | EPCSAFT_NEUTRAL_REFERENCE_DERIVATIVE_PARAMETERS_V1)
        && reference.value.derivative_availability
            == EPCSAFT_NEUTRAL_REFERENCE_DERIVATIVE_NONE_V1
        && std::strcmp(reference.reference_branch, "liquid") == 0
        && std::isfinite(phase.pressure_pa)
        && std::isfinite(inverse.volume_m3)
        && inverse.volume_m3 > 0.0
        && std::abs(inverse_gradient[3] + inverse.volume_m3)
            <= 1.0e-14 * inverse.volume_m3
        && std::isfinite(mixed[0])
        && std::all_of(
            packing_state_gradient.begin(),
            packing_state_gradient.end(),
            [](double value) { return std::isfinite(value); }
        )
        && std::all_of(
            packing_state_hessian.begin(),
            packing_state_hessian.end(),
            [](double value) { return std::isfinite(value); }
        )
        && std::isfinite(pressure_derivatives[0])
        && std::isfinite(parameter_derivatives[0])
        && neutral_inverse_packing_is_available(argv[2])
        && fixed_topology_association_is_available(argv[3], argv[5])
        && !fixed_topology_association_is_available(argv[4], argv[5]);
    std::cout << "component_count=" << (sdk == nullptr ? 0 : sdk->component_count)
              << '\n'
              << "phase_pressure_pa=" << phase.pressure_pa << '\n'
              << "reference_branch=" << reference.reference_branch << '\n';
    epcsaft_native_model_destroy_v1(model);
    return valid ? 0 : 4;
}
