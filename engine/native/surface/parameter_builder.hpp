#pragma once

#include "core/parameters.hpp"
#include "properties/thermochemistry.hpp"

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
    std::string base_choice;
    std::string association_choice;
    std::string electrolyte_choice;
    std::string ionic_dispersion_choice;
    std::string permittivity_choice;
    std::string polar_choice;
    std::vector<ParameterEntry> entries;
};

EPC_API Parameters build_parameters(const ParameterDocument& document);
EPC_API ThermochemistryRecord prepare_thermochemistry(ThermochemistryRecord record, double source_gas_constant);
// The start document with each fitted slot's entry magnitude replaced; refuses unless it rebuilds `fitted` bitwise.
EPC_API ParameterDocument fitted_document(const ParameterDocument& start, const Parameters& fitted);

}  // namespace epcsaft
