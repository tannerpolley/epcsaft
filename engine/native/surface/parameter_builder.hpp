#pragma once

#include "core/parameters.hpp"
#include "properties/thermochemistry.hpp"

#include <cstddef>
#include <optional>
#include <string>
#include <vector>

namespace epcsaft {

enum class DomainQuantity { SoluteMolality, IonicStrength, IonMoleFraction, Temperature, Pressure };
enum class DomainLimit { Minimum, Maximum };
struct DomainBound final {
    std::string packet, domain_id, solute, unit;
    DomainQuantity quantity;
    DomainLimit limit;
    std::size_t solvent;
    std::vector<std::size_t> solutes;
    double bound;
};
struct OutsideDomain final {
    std::string packet, domain_id, variable, side, unit, solute;
    double value, bound;
    std::size_t row, phase;
};
struct DomainBounds final {
    std::vector<DomainBound> values;
};

EPC_API DomainBounds admit_domains(const Parameters&, std::vector<DomainBound>);
EPC_API void outside_domains(std::vector<OutsideDomain>&, const Parameters&, const DomainBounds&, const std::vector<double>&,
    std::optional<double> temperature, std::optional<double> pressure, std::size_t row = 0, std::size_t phase = 0,
    bool validate_composition = false);

struct ParameterEntry final {
    std::string identity;
    std::string field;
    std::vector<std::size_t> index;
    double magnitude = 0.0;
    std::string unit;
};

struct ParameterDocument final {
    std::vector<std::string> component_ids;
    std::vector<std::string> site_ids;  // "component_id/site_id" per site, or empty
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
