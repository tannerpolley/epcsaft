#pragma once
#include <chrono>
#ifdef _WIN32
#define NOMINMAX
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#else
#include <time.h>
#endif
struct QualificationClock {
    std::chrono::steady_clock::time_point elapsed = std::chrono::steady_clock::now();
    double cpu = process_cpu();
    static double process_cpu() {
#ifdef _WIN32
        FILETIME created, exited, kernel, user;
        if (!GetProcessTimes(GetCurrentProcess(), &created, &exited, &kernel, &user)) throw std::runtime_error("GetProcessTimes failed");
        ULARGE_INTEGER k{}, u{}; k.LowPart=kernel.dwLowDateTime; k.HighPart=kernel.dwHighDateTime;
        u.LowPart=user.dwLowDateTime; u.HighPart=user.dwHighDateTime;
        return (double(k.QuadPart)+double(u.QuadPart))*1e-7;
#else
        timespec value{}; if(clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&value)) throw std::runtime_error("clock_gettime failed");
        return value.tv_sec + value.tv_nsec*1e-9;
#endif
    }
    double seconds() const { return std::chrono::duration<double>(std::chrono::steady_clock::now()-elapsed).count(); }
    double cpu_seconds() const { return process_cpu()-cpu; }
};
