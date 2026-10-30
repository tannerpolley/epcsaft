#pragma once

#if defined(_WIN32)
#if defined(EPCSAFT_NATIVE_STATIC)
#define EPC_API
#elif defined(EPCSAFT_NATIVE_BUILD)
#define EPC_API __declspec(dllexport)
#else
#define EPC_API __declspec(dllimport)
#endif
#define EPC_EXCEPTION EPC_API
#else
#define EPC_API
#define EPC_EXCEPTION __attribute__((visibility("default")))
#endif
