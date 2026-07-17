#pragma once

#include "eos.hpp"

#include <pybind11/pybind11.h>

#include <string>

namespace epcsaft {

pybind11::capsule make_native_sdk_v1(
    const NativeModel& model,
    const std::string& parameter_fingerprint,
    double temperature_min_k,
    double temperature_max_k
);

}  // namespace epcsaft
