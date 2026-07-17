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

namespace epcsaft {
namespace {

struct NativeSdkEndpoint final {
    NativeModel model;
    std::string parameter_fingerprint;
    double temperature_min_k;
    double temperature_max_k;
    epcsaft_native_sdk_v1 table;

    NativeSdkEndpoint(
        const NativeModel& source_model,
        std::string fingerprint,
        double minimum_temperature,
        double maximum_temperature
    )
        : model(source_model),
          parameter_fingerprint(std::move(fingerprint)),
          temperature_min_k(minimum_temperature),
          temperature_max_k(maximum_temperature),
          table{
              EPCSAFT_NATIVE_SDK_V1_ABI_VERSION,
              sizeof(epcsaft_native_sdk_v1),
              sizeof(epcsaft_phase_block_result_v1),
              this,
              &NativeSdkEndpoint::evaluate,
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
            std::copy(block.gradient.begin(), block.gradient.end(), result->gradient);
            std::copy(block.hessian.begin(), block.hessian.end(), result->hessian);
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
};

}  // namespace

pybind11::capsule make_native_sdk_v1(
    const NativeModel& model,
    const std::string& parameter_fingerprint,
    double temperature_min_k,
    double temperature_max_k
) {
    if (!model.supports_native_sdk_v1()) {
        throw std::invalid_argument(
            "native SDK v1 requires one pure nonassociating neutral model"
        );
    }
    if (parameter_fingerprint.empty()
        || parameter_fingerprint.size() >= EPCSAFT_NATIVE_SDK_V1_FINGERPRINT_SIZE) {
        throw std::invalid_argument("invalid parameter fingerprint");
    }
    auto* endpoint = new NativeSdkEndpoint(
        model,
        parameter_fingerprint,
        temperature_min_k,
        temperature_max_k
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
