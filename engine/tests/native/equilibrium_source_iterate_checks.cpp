#include <cstdio>

#if defined(__linux__)
#include "core/trace.hpp"
#include "equilibrium/gibbs_minimizer.hpp"
#include "surface/parameter_builder.hpp"
#include "recording_support.hpp"

#include <algorithm>
#include <cstdint>
#include <dlfcn.h>
#include <optional>
#include <stdexcept>
#include <string_view>

namespace {
using epcsaft::core::Trace;

bool temperature_source(const Trace& trace, const epcsaft::core::Span& span) {
    return span.parent && std::string_view(span.name) == "solve_equilibrium"
        && std::string_view(trace.spans[*span.parent].name) == "equilibrium_recovery";
}

std::optional<std::uint32_t> failed_source;
}

EPC_API bool Ipopt::TNLP::get_curr_iterate(
    const Ipopt::IpoptData* data, Ipopt::IpoptCalculatedQuantities* quantities,
    bool scaled, Ipopt::Index n, Ipopt::Number* x, Ipopt::Number* lower, Ipopt::Number* upper,
    Ipopt::Index m, Ipopt::Number* g, Ipopt::Number* lambda
) const {
    if (const auto* trace = epcsaft::core::current_trace()) {
        for (auto id = trace->current; id; id = trace->spans[*id].parent) {
            const auto& span = trace->spans[*id];
            if (std::string_view(span.name) != "solve_equilibrium") continue;
            if (temperature_source(*trace, span)) {
                if (!failed_source) failed_source = id;
                if (failed_source == id) return false;
            }
            break;
        }
    }
    // The Linux test executable interposes this symbol from Ipopt 3.14.20.
    using GetIterate = bool (*)(const TNLP*, const IpoptData*, IpoptCalculatedQuantities*,
        bool, Index, Number*, Number*, Number*, Index, Number*, Number*);
    static const auto real = reinterpret_cast<GetIterate>(dlsym(RTLD_NEXT,
        "_ZNK5Ipopt4TNLP16get_curr_iterateEPKNS_9IpoptDataEPNS_25IpoptCalculatedQuantitiesEbiPdS6_S6_iS6_S6_"));
    if (!real) throw std::runtime_error("could not resolve Ipopt::TNLP::get_curr_iterate");
    return real(this, data, quantities, scaled, n, x, lower, upper, m, g, lambda);
}
#endif

int main(int argc, char** argv) {
#if defined(__linux__)
    return recording_support::run(argc, argv, {
        {{"native.equilibrium.mdea-bubble-source-iterate", "mdea-bubble-source-iterate"}, [] {
            using namespace epcsaft;
            using namespace epcsaft::equilibrium;
#include "fixtures/mdea-411-bubble.inc"
            bubble.temperature = FreeValue::fixed(333.15);
            bubble.feed.values = {1.368, 8.615384615384617, 0, 0, 0, 0, 1, 0};
            const auto target = compile_problem(model, bubble);
            const auto limits = DriverLimits::recovery_allowance(std::chrono::steady_clock::time_point::max(),
                static_cast<std::size_t>(policy::EQUILIBRIUM_RECOVERY_WORK.value));
            core::Capture capture(true);
            bool propagated = false;
            try {
                const auto reached = recover_from_source(target, limits);
                std::printf("recovery returned: %s\n", reached.message.c_str());
            } catch (const std::logic_error& error) {
                std::printf("recovery threw: %s\n", error.what());
                propagated = std::string_view(error.what()) == "phase amount iterate or complementarity unavailable";
            }
            const auto& trace = capture.trace();
            const auto sources = std::count_if(trace.spans.begin(), trace.spans.end(),
                [&](const auto& span) { return temperature_source(trace, span); });
            std::printf("injected %s, temperature sources entered: %td\n", failed_source ? "yes" : "no", sources);
            if (!failed_source || !propagated || sources != 1)
                throw std::runtime_error("first source's Ipopt iterate failure must propagate before a second source");
        }},
    });
#else
    (void)argc;
    (void)argv;
    std::puts("Skipped: Ipopt iterate fault injection requires Linux symbol interposition.");
    return 77;
#endif
}
