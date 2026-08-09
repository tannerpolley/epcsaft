#pragma once

#include <string>

#include "phase_equilibrium.hpp"

namespace epcsaft_engine::equilibrium {

[[nodiscard]] std::string phase_equilibrium_result_to_json(
    const PhaseEquilibriumResult& result
);

}  // namespace epcsaft_engine::equilibrium
