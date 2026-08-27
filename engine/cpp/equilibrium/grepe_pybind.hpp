#pragma once

#include "grepe.hpp"

#include <pybind11/pybind11.h>

#include <string>

namespace epcsaft_engine::equilibrium {

[[nodiscard]] GrepeRequest parse_grepe_request(
    const pybind11::handle &problem, bool active_parameters_requested,
    const std::string &source_reference_fingerprint,
    const GrepeIntensive *temperature = nullptr,
    const GrepeIntensive *pressure = nullptr);

} // namespace epcsaft_engine::equilibrium
