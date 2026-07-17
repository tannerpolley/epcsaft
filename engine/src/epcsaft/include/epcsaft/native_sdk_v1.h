#ifndef EPCSAFT_NATIVE_SDK_V1_H
#define EPCSAFT_NATIVE_SDK_V1_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#define EPCSAFT_NATIVE_SDK_V1_CAPSULE_NAME "epcsaft.native_sdk.v1"
#define EPCSAFT_NATIVE_SDK_V1_ABI_VERSION 1u
#define EPCSAFT_NATIVE_SDK_V1_FINGERPRINT_SIZE 72u
#define EPCSAFT_NATIVE_SDK_V1_ERROR_SIZE 160u

enum epcsaft_native_status_v1 {
    EPCSAFT_NATIVE_STATUS_OK_V1 = 0,
    EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT_V1 = 1,
    EPCSAFT_NATIVE_STATUS_UNSUPPORTED_MODEL_V1 = 2,
    EPCSAFT_NATIVE_STATUS_DOMAIN_ERROR_V1 = 3,
    EPCSAFT_NATIVE_STATUS_EVALUATION_ERROR_V1 = 4,
    EPCSAFT_NATIVE_STATUS_ABI_MISMATCH_V1 = 5
};

typedef struct epcsaft_phase_block_result_v1 {
    uint32_t struct_size;
    int32_t status;
    double helmholtz_over_rt_reference_amount;
    double gradient[2];
    double hessian[4];
    double third[8];
    double pressure_pa;
    double chemical_potential_over_rt;
    char parameter_fingerprint[EPCSAFT_NATIVE_SDK_V1_FINGERPRINT_SIZE];
    char error[EPCSAFT_NATIVE_SDK_V1_ERROR_SIZE];
} epcsaft_phase_block_result_v1;

typedef int (*epcsaft_evaluate_pure_phase_v1)(
    const void* model_context,
    double temperature_k,
    double amount_mol,
    double volume_m3,
    epcsaft_phase_block_result_v1* result
);

typedef struct epcsaft_native_sdk_v1 {
    uint32_t abi_version;
    size_t table_size;
    size_t result_size;
    const void* model_context;
    epcsaft_evaluate_pure_phase_v1 evaluate_pure_phase;
} epcsaft_native_sdk_v1;

#ifdef __cplusplus
}
#endif

#endif
