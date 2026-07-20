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

enum epcsaft_parameterized_phase_coordinate_v1 {
    EPCSAFT_PHASE_AMOUNT_MOL_V1 = 0,
    EPCSAFT_PHASE_VOLUME_M3_V1 = 1,
    EPCSAFT_PHASE_SEGMENT_COUNT_V1 = 2,
    EPCSAFT_PHASE_SEGMENT_DIAMETER_ANGSTROM_V1 = 3,
    EPCSAFT_PHASE_DISPERSION_ENERGY_OVER_K_V1 = 4,
    EPCSAFT_PHASE_PARAMETERIZED_COORDINATE_COUNT_V1 = 5
};

typedef struct epcsaft_parameterized_phase_block_result_v1 {
    uint32_t struct_size;
    int32_t status;
    double helmholtz_over_rt_reference_amount;
    double gradient[EPCSAFT_PHASE_PARAMETERIZED_COORDINATE_COUNT_V1];
    double hessian[
        EPCSAFT_PHASE_PARAMETERIZED_COORDINATE_COUNT_V1 *
        EPCSAFT_PHASE_PARAMETERIZED_COORDINATE_COUNT_V1
    ];
    double pressure_pa;
    double chemical_potential_over_rt;
    char parameter_fingerprint[EPCSAFT_NATIVE_SDK_V1_FINGERPRINT_SIZE];
    char error[EPCSAFT_NATIVE_SDK_V1_ERROR_SIZE];
} epcsaft_parameterized_phase_block_result_v1;

typedef int (*epcsaft_evaluate_pure_phase_v1)(
    const void* model_context,
    double temperature_k,
    double amount_mol,
    double volume_m3,
    epcsaft_phase_block_result_v1* result
);

typedef int (*epcsaft_evaluate_pure_phase_parameters_v1)(
    const void* model_context,
    double temperature_k,
    double amount_mol,
    double volume_m3,
    double segment_count,
    double segment_diameter_angstrom,
    double dispersion_energy_over_k,
    epcsaft_parameterized_phase_block_result_v1* result
);

/*
 * A mixture result struct-size mismatch is return-code authoritative. The
 * callbacks return ABI_MISMATCH without writing an untrusted result layout.
 * After exact-size acceptance, failures also set result status and error.
 */
typedef struct epcsaft_mixture_phase_block_result_v1 {
    uint32_t struct_size;
    int32_t status;
    size_t coordinate_count;
    size_t gradient_capacity;
    size_t hessian_capacity;
    double* gradient;
    double* hessian;
    double helmholtz_over_rt_reference_amount;
    double pressure_pa;
    char parameter_fingerprint[EPCSAFT_NATIVE_SDK_V1_FINGERPRINT_SIZE];
    char error[EPCSAFT_NATIVE_SDK_V1_ERROR_SIZE];
} epcsaft_mixture_phase_block_result_v1;

typedef int (*epcsaft_evaluate_mixture_phase_v1)(
    const void* model_context,
    double temperature_k,
    const double* amounts_mol,
    size_t amount_count,
    double volume_m3,
    epcsaft_mixture_phase_block_result_v1* result
);

typedef int (*epcsaft_evaluate_mixture_phase_kij_v1)(
    const void* model_context,
    double temperature_k,
    const double* amounts_mol,
    size_t amount_count,
    double volume_m3,
    double k_ij,
    epcsaft_mixture_phase_block_result_v1* result
);

typedef struct epcsaft_native_sdk_v1 {
    uint32_t abi_version;
    size_t table_size;
    size_t result_size;
    const void* model_context;
    epcsaft_evaluate_pure_phase_v1 evaluate_pure_phase;
    size_t parameterized_result_size;
    epcsaft_evaluate_pure_phase_parameters_v1 evaluate_pure_phase_parameters;
    size_t component_count;
    size_t mixture_result_size;
    epcsaft_evaluate_mixture_phase_v1 evaluate_mixture_phase;
    epcsaft_evaluate_mixture_phase_kij_v1 evaluate_mixture_phase_kij;
    const char* const* component_ids;
    const int32_t* component_charges;
    epcsaft_evaluate_mixture_phase_v1 evaluate_electrolyte_phase;
} epcsaft_native_sdk_v1;

#ifdef __cplusplus
}
#endif

#endif
