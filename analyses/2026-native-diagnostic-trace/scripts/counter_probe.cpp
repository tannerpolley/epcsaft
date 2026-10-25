#include "case_inputs.hpp"
#include <cstdio>
int main() {
    for (int i = 0; i < 5; ++i) for (int mode = 0; mode < 2; ++mode) {
        auto [parameters, problem] = case_inputs(i);
        epcsaft::HelmholtzModel model(parameters);
        std::optional<epcsaft::core::Capture> capture;
        if (mode) capture.emplace();
        const auto result = epcsaft::equilibrium::solve_equilibrium(model, problem);
        const auto c = model.counters();
        std::printf("%d,%d,%d,%zu,%zu,%zu,%zu,%zu,%zu,%zu,%zu\n", i, mode, result.success, result.iterations, result.evaluations,
            c.forward, c.reverse, c.hessian_calls, c.recordings, c.directional_paths, c.maximum_energy_order);
    }
}
