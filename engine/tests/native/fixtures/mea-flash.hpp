#pragma once

#include "equilibrium/problem.hpp"
#include "surface/parameter_builder.hpp"

namespace mea_flash_fixture {

inline auto declaration(double temperature, double pressure, double loading) {
    using namespace epcsaft;
    using namespace epcsaft::equilibrium;
#include "mea-nominally-unloaded.inc"
    // Paper data archive v2, examples/mea-flash/mea_case.py and Listing 1.
    problem.temperature = FreeValue::fixed(temperature);
    problem.pressure = FreeValue::fixed(pressure);
    problem.feed.values = {loading, 1., 7.929776437398854, 0., 0., 0., 0., 0., 0.};
    problem.phases[0].composition_guess.clear();
    problem.phases.push_back(PhaseSpec{"V", PhaseKind::Vapor, {0, 1, 2}, true});
    for (auto& ratio : problem.source_reference->solute_molality_ratios) ratio *= 8;
    return std::pair{model.parameters(), std::move(problem)};
}

}
