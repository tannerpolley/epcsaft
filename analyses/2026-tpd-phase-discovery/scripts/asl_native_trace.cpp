#include "core/trace.hpp"
#include <dlfcn.h>
#include <cstring>
#include <sstream>
struct AslTrace { epcsaft::core::Trace trace; epcsaft::core::Trace** slot; void* handle; };
extern "C" void* __tls_get_addr(void*);
extern "C" void* begin_asl_trace(const char* library, std::size_t offset) {
    void* handle = dlopen(library, RTLD_NOW | RTLD_LOCAL);
    std::size_t module = 0;
    if (!handle || dlinfo(handle, RTLD_DI_TLS_MODID, &module) || !module) return nullptr;
    std::size_t index[2] = {module, offset};
    auto slot = static_cast<epcsaft::core::Trace**>(__tls_get_addr(index));
    if (!slot || *slot) { dlclose(handle); return nullptr; }
    auto state = new AslTrace{{}, slot, handle};
    *slot = &state->trace;
    return state;
}
extern "C" const char* end_asl_trace(void* opaque) {
    auto state = static_cast<AslTrace*>(opaque);
    *state->slot = nullptr;
    static thread_local std::string output;
    std::ostringstream rows;
    std::size_t drivers = 0, iterations = 0, evaluations = 0, hessians = 0, trials = 0;
    for (const auto& span : state->trace.spans) trials += std::strcmp(span.name, "tangent_plane_trial") == 0;
    for (const auto& event : state->trace.events) if (std::strcmp(event.name, "hessian_calls") == 0) {
        const auto& span = state->trace.spans.at(event.span);
        if (!span.convergence) continue;
        if (drivers++) rows << ',';
        const auto& c = *span.convergence;
        iterations += c.iterations; evaluations += c.evaluations;
        const auto calls = static_cast<std::size_t>(event.values.at(0)); hessians += calls;
        rows << "{\"scope\":\"" << span.name << "\",\"iterations\":" << c.iterations
             << ",\"evaluations\":" << c.evaluations << ",\"hessian_calls\":" << calls << '}';
    }
    std::ostringstream summary;
    summary << "{\"spans\":" << state->trace.spans.size() << ",\"driver_finalizations\":" << drivers
            << ",\"iterations\":" << iterations << ",\"evaluations\":" << evaluations
            << ",\"hessian_calls\":" << hessians << ",\"tangent_plane_trials\":" << trials
            << ",\"drivers\":[" << rows.str() << "]}";
    output = summary.str();
    void* handle = state->handle; delete state; dlclose(handle);
    return output.c_str();
}
