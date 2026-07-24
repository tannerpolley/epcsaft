#ifndef EPCSAFT_NATIVE_MODEL_V1_H
#define EPCSAFT_NATIVE_MODEL_V1_H

#include <epcsaft/native_sdk_v1.h>

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct epcsaft_native_model_handle_v1 epcsaft_native_model_handle_v1;

/*
 * Load one immutable, Provider-resolved model artifact. The returned handle
 * owns the SDK table and model context; every table pointer remains valid until
 * epcsaft_native_model_destroy_v1 is called.
 */
epcsaft_native_model_handle_v1* epcsaft_native_model_load_v1(
    const char* path,
    char* error,
    size_t error_capacity
);

const epcsaft_native_sdk_v1* epcsaft_native_model_sdk_v1(
    const epcsaft_native_model_handle_v1* model
);

/* The returned string remains valid until the model handle is destroyed. */
const char* epcsaft_native_model_fingerprint_v1(
    const epcsaft_native_model_handle_v1* model
);

void epcsaft_native_model_destroy_v1(epcsaft_native_model_handle_v1* model);

#ifdef __cplusplus
}
#endif

#endif
