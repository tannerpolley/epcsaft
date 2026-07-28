#include "native_sdk.hpp"

#include "model.hpp"

#include <epcsaft/native_model_v1.h>

#include <nlohmann/json.hpp>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

struct epcsaft_native_model_handle_v1 {
    const epcsaft_native_sdk_v1* sdk;
    std::string fingerprint;
};

namespace {

using Json = nlohmann::json;

double optional_limit(const Json& object, const char* name) {
    const auto& value = object.at(name);
    return value.is_null()
        ? std::numeric_limits<double>::quiet_NaN()
        : value.get<double>();
}

void copy_error(const std::string& message, char* error, std::size_t capacity) {
    if (error == nullptr || capacity == 0) {
        return;
    }
    const std::size_t count = std::min(message.size(), capacity - 1);
    std::memcpy(error, message.data(), count);
    error[count] = '\0';
}

epcsaft::NativeParameters parameters(const Json& value) {
    return epcsaft::NativeParameters{
        value.at("segment_count").get<std::vector<double>>(),
        value.at("sigma_constant_angstrom").get<std::vector<double>>(),
        value.at("sigma_term_component").get<std::vector<int>>(),
        value.at("sigma_term_amplitude_angstrom").get<std::vector<double>>(),
        value.at("sigma_term_exponent_per_k").get<std::vector<double>>(),
        value.at("dispersion_energy_over_k_k").get<std::vector<double>>(),
        value.at("k_ij").get<std::vector<double>>(),
        value.at("l_ij").get<std::vector<double>>(),
        value.at("site_component").get<std::vector<int>>(),
        value.at("site_multiplicity").get<std::vector<int>>(),
        value.at("association_active").get<std::vector<int>>(),
        value.at("association_energy_over_k_k").get<std::vector<double>>(),
        value.at("association_volume").get<std::vector<double>>(),
        value.at("molar_mass_kg_per_mol").get<std::vector<double>>(),
        value.at("charge_number").get<std::vector<int>>(),
        value.at("relative_permittivity").get<std::vector<double>>(),
        value.at("dielectric_ion_suppression_coefficient").get<double>(),
        value.at("born_diameter_angstrom").get<std::vector<double>>(),
        value.at("solvation_factor").get<std::vector<double>>(),
        value.at("ionic_region_relative_permittivity").get<double>(),
        value.at("electrolyte_model").get<int>(),
        value.at("schreckenberg_dielectric_volume_m3_per_mol")
            .get<std::vector<double>>(),
        value.at("schreckenberg_dielectric_temperature_k")
            .get<std::vector<double>>(),
        value.at("zuber_ion_suppression_coefficient").get<std::vector<double>>(),
        value.at("rueben_dipole_scaling_debye2").get<std::vector<double>>(),
        value.at("rueben_polarizability_scaling_angstrom3")
            .get<std::vector<double>>(),
        value.at("rueben_correlation_integral_parameter")
            .get<std::vector<double>>(),
        value.at("permittivity_model").get<int>(),
    };
}

}  // namespace

extern "C" epcsaft_native_model_handle_v1* epcsaft_native_model_load_v1(
    const char* path,
    char* error,
    std::size_t error_capacity
) {
    try {
        if (path == nullptr || path[0] == '\0') {
            throw std::invalid_argument("resolved model path is required");
        }
        std::ifstream stream(path);
        if (!stream) {
            throw std::invalid_argument("resolved model file cannot be opened");
        }
        Json document;
        stream >> document;
        if (document.at("schema") != "epcsaft-resolved-native-model-v1") {
            throw std::invalid_argument("unsupported resolved model schema");
        }
        const Json& sdk = document.at("sdk");
        const epcsaft::NativeModel model(
            parameters(document.at("native_parameters"))
        );
        const std::string fingerprint =
            sdk.at("parameter_fingerprint").get<std::string>();
        const epcsaft_native_sdk_v1* table = epcsaft::create_native_sdk_v1(
            model,
            fingerprint,
            optional_limit(sdk, "temperature_min_k"),
            optional_limit(sdk, "temperature_max_k"),
            sdk.at("component_ids").get<std::vector<std::string>>(),
            optional_limit(sdk, "ion_mole_fraction_max"),
            sdk.at("phase_callbacks_admitted").get<bool>(),
            sdk.at("active_born_admitted").get<bool>(),
            sdk.at("aqueous_kij_admitted").get<bool>(),
            optional_limit(sdk, "pressure_min_pa"),
            optional_limit(sdk, "pressure_max_pa"),
            sdk.at("neutral_reference_admitted").get<bool>(),
            sdk.at("unknown_applicability_admitted").get<bool>(),
            sdk.at("aqueous_solvation_factor_admitted").get<bool>(),
            sdk.at("staged_trial_bundle_admitted").get<bool>(),
            sdk.value("topology_fingerprint", std::string{})
        );
        return new epcsaft_native_model_handle_v1{table, fingerprint};
    } catch (const std::exception& exception) {
        copy_error(exception.what(), error, error_capacity);
        return nullptr;
    } catch (...) {
        copy_error("unknown native model load failure", error, error_capacity);
        return nullptr;
    }
}

extern "C" const epcsaft_native_sdk_v1* epcsaft_native_model_sdk_v1(
    const epcsaft_native_model_handle_v1* model
) {
    return model == nullptr ? nullptr : model->sdk;
}

extern "C" const char* epcsaft_native_model_fingerprint_v1(
    const epcsaft_native_model_handle_v1* model
) {
    return model == nullptr ? nullptr : model->fingerprint.c_str();
}

extern "C" void epcsaft_native_model_destroy_v1(
    epcsaft_native_model_handle_v1* model
) {
    if (model != nullptr) {
        epcsaft::destroy_native_sdk_v1(model->sdk);
        delete model;
    }
}
