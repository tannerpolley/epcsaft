#pragma once

#include "core/parameters.hpp"

#include <cstddef>
#include <string>
#include <vector>

namespace epcsaft {

struct ParameterEntry final {
    std::string identity;
    std::string field;
    std::vector<std::size_t> index;
    double magnitude = 0.0;
    std::string unit;
};

struct ParameterDocument final {
    std::vector<std::string> component_ids;
    std::vector<int> site_component;
    std::vector<int> site_multiplicity;
    std::string electrolyte_choice;
    std::string ionic_dispersion_choice;
    std::string permittivity_choice;
    std::vector<ParameterEntry> entries;
};

Parameters build_parameters(const ParameterDocument& document);

}  // namespace epcsaft
