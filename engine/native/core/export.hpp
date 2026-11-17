#pragma once

#if defined(_WIN32)
#if defined(EPCSAFT_NATIVE_STATIC)
#define EPC_API
#elif defined(EPCSAFT_NATIVE_BUILD)
#define EPC_API __declspec(dllexport)
#else
#define EPC_API __declspec(dllimport)
#endif
#else
// The library builds with hidden default visibility, so an unmarked symbol fails the Linux link as it fails on Windows.
#define EPC_API __attribute__((visibility("default")))
#endif
#define EPC_EXCEPTION EPC_API
