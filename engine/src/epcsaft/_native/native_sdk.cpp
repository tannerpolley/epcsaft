#include "native_sdk.hpp"
#include "phase_block.hpp"

#include <epcsaft/native_sdk_v1.h>

#include <pybind11/pybind11.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace epcsaft {
namespace {

struct NativeSdkEndpoint final {
    NativeModel model;
    std::string parameter_fingerprint;
    double temperature_min_k;
    double temperature_max_k;
    double ion_mole_fraction_max;
    std::vector<std::string> component_ids;
    std::vector<const char*> component_id_pointers;
    std::vector<int32_t> component_charges;
    epcsaft_native_sdk_v1 table;

    NativeSdkEndpoint(
        const NativeModel& source_model,
        std::string fingerprint,
        double minimum_temperature,
        double maximum_temperature,
        std::vector<std::string> identifiers,
        double maximum_ion_fraction
    )
        : model(source_model),
          parameter_fingerprint(std::move(fingerprint)),
          temperature_min_k(minimum_temperature),
          temperature_max_k(maximum_temperature),
          ion_mole_fraction_max(maximum_ion_fraction),
          component_ids(std::move(identifiers)),
          component_id_pointers([&] {
              std::vector<const char*> pointers;
              pointers.reserve(component_ids.size());
              for (const auto& identifier : component_ids) {
                  pointers.push_back(identifier.c_str());
              }
              return pointers;
          }()),
          component_charges(
              source_model.charge_numbers().begin(),
              source_model.charge_numbers().end()
          ),
          table{
              EPCSAFT_NATIVE_SDK_V1_ABI_VERSION,
              sizeof(epcsaft_native_sdk_v1),
              sizeof(epcsaft_phase_block_result_v1),
              this,
              &NativeSdkEndpoint::evaluate,
              sizeof(epcsaft_parameterized_phase_block_result_v1),
              &NativeSdkEndpoint::evaluate_parameters,
              source_model.component_count(),
              sizeof(epcsaft_mixture_phase_block_result_v1),
              &NativeSdkEndpoint::evaluate_mixture,
              &NativeSdkEndpoint::evaluate_mixture_kij,
              component_id_pointers.data(),
              component_charges.data(),
              &NativeSdkEndpoint::evaluate_electrolyte,
          } {}

    static int evaluate(
        const void* context,
        double temperature_k,
        double amount_mol,
        double volume_m3,
        epcsaft_phase_block_result_v1* result
    ) noexcept {
        if (result == nullptr) {
            return EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT_V1;
        }
        if (result->struct_size != sizeof(epcsaft_phase_block_result_v1)) {
            return EPCSAFT_NATIVE_STATUS_ABI_MISMATCH_V1;
        }
        const uint32_t result_size = result->struct_size;
        std::memset(result, 0, sizeof(*result));
        result->struct_size = result_size;

        auto fail = [&](int status, const char* message) {
            result->status = status;
            std::snprintf(result->error, sizeof(result->error), "%s", message);
            return status;
        };
        if (context == nullptr) {
            return fail(
                EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT_V1,
                "model context is null"
            );
        }
        if (!std::isfinite(temperature_k) || !std::isfinite(amount_mol)
            || !std::isfinite(volume_m3) || temperature_k <= 0.0
            || amount_mol <= 0.0 || volume_m3 <= 0.0) {
            return fail(
                EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT_V1,
                "temperature, amount, and volume must be positive and finite"
            );
        }

        const auto& endpoint = *static_cast<const NativeSdkEndpoint*>(context);
        if ((std::isfinite(endpoint.temperature_min_k)
             && temperature_k < endpoint.temperature_min_k)
            || (std::isfinite(endpoint.temperature_max_k)
                && temperature_k > endpoint.temperature_max_k)) {
            return fail(
                EPCSAFT_NATIVE_STATUS_DOMAIN_ERROR_V1,
                "temperature is outside the parameter source domain"
            );
        }

        try {
            const auto block = endpoint.model.evaluate_pure_phase(
                temperature_k, amount_mol, volume_m3
            );
            result->helmholtz_over_rt_reference_amount =
                block.helmholtz_over_rt_reference_amount;
            std::copy(
                block.gradient.begin(), block.gradient.end(), result->gradient
            );
            std::copy(
                block.hessian.begin(), block.hessian.end(), result->hessian
            );
            std::copy(block.third.begin(), block.third.end(), result->third);
            result->pressure_pa = block.pressure_pa;
            result->chemical_potential_over_rt = block.chemical_potential_over_rt;
            std::snprintf(
                result->parameter_fingerprint,
                sizeof(result->parameter_fingerprint),
                "%s",
                endpoint.parameter_fingerprint.c_str()
            );
            result->status = EPCSAFT_NATIVE_STATUS_OK_V1;
            return result->status;
        } catch (const std::invalid_argument& error) {
            return fail(EPCSAFT_NATIVE_STATUS_DOMAIN_ERROR_V1, error.what());
        } catch (const std::logic_error& error) {
            return fail(EPCSAFT_NATIVE_STATUS_UNSUPPORTED_MODEL_V1, error.what());
        } catch (const std::exception& error) {
            return fail(EPCSAFT_NATIVE_STATUS_EVALUATION_ERROR_V1, error.what());
        } catch (...) {
            return fail(
                EPCSAFT_NATIVE_STATUS_EVALUATION_ERROR_V1,
                "unknown provider evaluation failure"
            );
        }
    }

    static int evaluate_parameters(
        const void* context,
        double temperature_k,
        double amount_mol,
        double volume_m3,
        double segment_count,
        double segment_diameter_angstrom,
        double dispersion_energy_over_k,
        epcsaft_parameterized_phase_block_result_v1* result
    ) noexcept {
        if (result == nullptr) {
            return EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT_V1;
        }
        if (result->struct_size
            != sizeof(epcsaft_parameterized_phase_block_result_v1)) {
            return EPCSAFT_NATIVE_STATUS_ABI_MISMATCH_V1;
        }
        const uint32_t result_size = result->struct_size;
        std::memset(result, 0, sizeof(*result));
        result->struct_size = result_size;

        auto fail = [&](int status, const char* message) {
            result->status = status;
            std::snprintf(result->error, sizeof(result->error), "%s", message);
            return status;
        };
        if (context == nullptr) {
            return fail(
                EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT_V1,
                "model context is null"
            );
        }
        const auto& endpoint = *static_cast<const NativeSdkEndpoint*>(context);
        if ((std::isfinite(endpoint.temperature_min_k)
             && temperature_k < endpoint.temperature_min_k)
            || (std::isfinite(endpoint.temperature_max_k)
                && temperature_k > endpoint.temperature_max_k)) {
            return fail(
                EPCSAFT_NATIVE_STATUS_DOMAIN_ERROR_V1,
                "temperature is outside the parameter source domain"
            );
        }

        try {
            const auto block = endpoint.model.evaluate_pure_phase_parameters(
                temperature_k,
                amount_mol,
                volume_m3,
                segment_count,
                segment_diameter_angstrom,
                dispersion_energy_over_k
            );
            result->helmholtz_over_rt_reference_amount =
                block.helmholtz_over_rt_reference_amount;
            std::copy(block.gradient.begin(), block.gradient.end(), result->gradient);
            std::copy(block.hessian.begin(), block.hessian.end(), result->hessian);
            result->pressure_pa = block.pressure_pa;
            result->chemical_potential_over_rt = block.chemical_potential_over_rt;
            std::snprintf(
                result->parameter_fingerprint,
                sizeof(result->parameter_fingerprint),
                "%s",
                endpoint.parameter_fingerprint.c_str()
            );
            result->status = EPCSAFT_NATIVE_STATUS_OK_V1;
            return result->status;
        } catch (const std::invalid_argument& error) {
            return fail(EPCSAFT_NATIVE_STATUS_DOMAIN_ERROR_V1, error.what());
        } catch (const std::logic_error& error) {
            return fail(EPCSAFT_NATIVE_STATUS_UNSUPPORTED_MODEL_V1, error.what());
        } catch (const std::exception& error) {
            return fail(EPCSAFT_NATIVE_STATUS_EVALUATION_ERROR_V1, error.what());
        } catch (...) {
            return fail(
                EPCSAFT_NATIVE_STATUS_EVALUATION_ERROR_V1,
                "unknown provider evaluation failure"
            );
        }
    }

    template <typename Evaluator>
    static int evaluate_mixture_common(
        const void* context,
        double temperature_k,
        const double* amounts_mol,
        size_t amount_count,
        double volume_m3,
        std::size_t additional_coordinate_count,
        bool additional_inputs_valid,
        const char* invalid_state_message,
        epcsaft_mixture_phase_block_result_v1* result,
        Evaluator evaluator
    ) noexcept {
        if (result == nullptr) {
            return EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT_V1;
        }
        if (result->struct_size != sizeof(epcsaft_mixture_phase_block_result_v1)) {
            return EPCSAFT_NATIVE_STATUS_ABI_MISMATCH_V1;
        }
        if (context == nullptr) {
            result->status = EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT_V1;
            std::snprintf(
                result->error,
                sizeof(result->error),
                "%s",
                "model context is null"
            );
            return result->status;
        }

        const auto& endpoint = *static_cast<const NativeSdkEndpoint*>(context);
        const std::size_t coordinate_count =
            endpoint.table.component_count + additional_coordinate_count;
        if (result->coordinate_count != coordinate_count
            || result->gradient_capacity < coordinate_count
            || result->hessian_capacity < coordinate_count * coordinate_count
            || result->gradient == nullptr || result->hessian == nullptr) {
            result->status = EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT_V1;
            std::snprintf(
                result->error,
                sizeof(result->error),
                "%s",
                "mixture result buffers do not match the model coordinate count"
            );
            return result->status;
        }

        result->status = EPCSAFT_NATIVE_STATUS_OK_V1;
        result->helmholtz_over_rt_reference_amount = 0.0;
        result->pressure_pa = 0.0;
        result->parameter_fingerprint[0] = '\0';
        result->error[0] = '\0';
        auto fail = [&](int status, const char* message) {
            result->status = status;
            std::snprintf(result->error, sizeof(result->error), "%s", message);
            return status;
        };

        if (amounts_mol == nullptr || amount_count != endpoint.table.component_count) {
            return fail(
                EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT_V1,
                "amount buffer does not match the model component count"
            );
        }
        if (!std::isfinite(temperature_k) || !std::isfinite(volume_m3)
            || temperature_k <= 0.0 || volume_m3 <= 0.0
            || !additional_inputs_valid
            || !std::all_of(amounts_mol, amounts_mol + amount_count, [](double amount) {
                return std::isfinite(amount) && amount > 0.0;
            })) {
            return fail(
                EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT_V1,
                invalid_state_message
            );
        }
        if ((std::isfinite(endpoint.temperature_min_k)
             && temperature_k < endpoint.temperature_min_k)
            || (std::isfinite(endpoint.temperature_max_k)
                && temperature_k > endpoint.temperature_max_k)) {
            return fail(
                EPCSAFT_NATIVE_STATUS_DOMAIN_ERROR_V1,
                "temperature is outside the parameter source domain"
            );
        }

        try {
            const std::vector<double> amounts(amounts_mol, amounts_mol + amount_count);
            const auto block = evaluator(endpoint.model, amounts);
            std::copy(block.gradient.begin(), block.gradient.end(), result->gradient);
            std::copy(block.hessian.begin(), block.hessian.end(), result->hessian);
            result->helmholtz_over_rt_reference_amount =
                block.helmholtz_over_rt_reference_amount;
            result->pressure_pa = block.pressure_pa;
            std::snprintf(
                result->parameter_fingerprint,
                sizeof(result->parameter_fingerprint),
                "%s",
                endpoint.parameter_fingerprint.c_str()
            );
            result->status = EPCSAFT_NATIVE_STATUS_OK_V1;
            return result->status;
        } catch (const std::invalid_argument& error) {
            return fail(EPCSAFT_NATIVE_STATUS_DOMAIN_ERROR_V1, error.what());
        } catch (const std::logic_error& error) {
            return fail(EPCSAFT_NATIVE_STATUS_UNSUPPORTED_MODEL_V1, error.what());
        } catch (const std::exception& error) {
            return fail(EPCSAFT_NATIVE_STATUS_EVALUATION_ERROR_V1, error.what());
        } catch (...) {
            return fail(
                EPCSAFT_NATIVE_STATUS_EVALUATION_ERROR_V1,
                "unknown provider evaluation failure"
            );
        }
    }

    static int evaluate_mixture(
        const void* context,
        double temperature_k,
        const double* amounts_mol,
        size_t amount_count,
        double volume_m3,
        epcsaft_mixture_phase_block_result_v1* result
    ) noexcept {
        return evaluate_mixture_common(
            context,
            temperature_k,
            amounts_mol,
            amount_count,
            volume_m3,
            1,
            true,
            "temperature, amounts, and volume must be positive and finite",
            result,
            [temperature_k, volume_m3](
                const NativeModel& model,
                const std::vector<double>& amounts
            ) {
                return model.evaluate_mixture_phase(
                    temperature_k, amounts, volume_m3
                );
            }
        );
    }

    static int evaluate_mixture_kij(
        const void* context,
        double temperature_k,
        const double* amounts_mol,
        size_t amount_count,
        double volume_m3,
        double k_ij,
        epcsaft_mixture_phase_block_result_v1* result
    ) noexcept {
        return evaluate_mixture_common(
            context,
            temperature_k,
            amounts_mol,
            amount_count,
            volume_m3,
            2,
            std::isfinite(k_ij),
            "temperature, amounts, volume, and k_ij must be finite and state coordinates positive",
            result,
            [temperature_k, volume_m3, k_ij](
                const NativeModel& model,
                const std::vector<double>& amounts
            ) {
                return model.evaluate_mixture_phase_kij(
                    temperature_k, amounts, volume_m3, k_ij
                );
            }
        );
    }

    static int evaluate_electrolyte(
        const void* context,
        double temperature_k,
        const double* amounts_mol,
        size_t amount_count,
        double volume_m3,
        epcsaft_mixture_phase_block_result_v1* result
    ) noexcept {
        const auto* endpoint = static_cast<const NativeSdkEndpoint*>(context);
        return evaluate_mixture_common(
            context,
            temperature_k,
            amounts_mol,
            amount_count,
            volume_m3,
            1,
            true,
            "temperature, amounts, and volume must be positive and finite",
            result,
            [temperature_k, volume_m3, endpoint](
                const NativeModel& model,
                const std::vector<double>& amounts
            ) {
                return model.evaluate_electrolyte_phase(
                    temperature_k,
                    amounts,
                    volume_m3,
                    endpoint->ion_mole_fraction_max
                );
            }
        );
    }
};

}  // namespace

pybind11::capsule make_native_sdk_v1(
    const NativeModel& model,
    const std::string& parameter_fingerprint,
    double temperature_min_k,
    double temperature_max_k,
    const std::vector<std::string>& component_ids,
    double ion_mole_fraction_max
) {
    if (!model.supports_native_sdk_v1()
        && !model.supports_mixture_native_sdk_v1()
        && !model.supports_electrolyte_native_sdk_v1()) {
        throw std::invalid_argument(
            "native SDK v1 requires a supported phase-block model"
        );
    }
    if (parameter_fingerprint.empty()
        || parameter_fingerprint.size() >= EPCSAFT_NATIVE_SDK_V1_FINGERPRINT_SIZE) {
        throw std::invalid_argument("invalid parameter fingerprint");
    }
    if (component_ids.size() != model.component_count()
        || std::any_of(
            component_ids.begin(), component_ids.end(), [](const auto& identifier) {
                return identifier.empty();
            }
        )) {
        throw std::invalid_argument("invalid component identifiers");
    }
    for (std::size_t left = 0; left < component_ids.size(); ++left) {
        for (std::size_t right = left + 1; right < component_ids.size(); ++right) {
            if (component_ids[left] == component_ids[right]) {
                throw std::invalid_argument("duplicate component identifier");
            }
        }
    }
    auto* endpoint = new NativeSdkEndpoint(
        model,
        parameter_fingerprint,
        temperature_min_k,
        temperature_max_k,
        component_ids,
        ion_mole_fraction_max
    );
    return pybind11::capsule(
        &endpoint->table,
        EPCSAFT_NATIVE_SDK_V1_CAPSULE_NAME,
        [](PyObject* capsule) {
            void* pointer = PyCapsule_GetPointer(
                capsule, EPCSAFT_NATIVE_SDK_V1_CAPSULE_NAME
            );
            if (pointer == nullptr) {
                PyErr_Clear();
                return;
            }
            const auto* table = static_cast<const epcsaft_native_sdk_v1*>(pointer);
            delete static_cast<const NativeSdkEndpoint*>(table->model_context);
        }
    );
}

}  // namespace epcsaft
