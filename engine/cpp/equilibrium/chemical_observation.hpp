#pragma once

#include <pybind11/pybind11.h>

namespace epcsaft_engine::equilibrium {

void bind_chemical_observation(pybind11::module_& module);

}  // namespace epcsaft_engine::equilibrium
