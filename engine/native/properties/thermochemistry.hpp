#pragma once

#include "contributions/ideal.hpp"

#include <cstddef>
#include <optional>
#include <string>
#include <vector>

namespace epcsaft {

enum class ThermochemistryStatus { MissingRecord, Available, OutsideInterval, Nonfinite };

struct ThermochemistryRecord final {
    double standard_pressure = 0.0;
    std::vector<std::vector<IdealInterval>> component_intervals;
};
struct IdealProperties final {
    double enthalpy = 0.0, internal_energy = 0.0, entropy = 0.0;
    double isobaric_heat_capacity = 0.0, isochoric_heat_capacity = 0.0;
};

struct IdealPropertiesEvaluation final {
    std::optional<IdealProperties> value;
    ThermochemistryStatus status = ThermochemistryStatus::Nonfinite;
};
void validate_thermochemistry(const ThermochemistryRecord&, std::size_t component_count);
IdealPropertiesEvaluation evaluate_ideal_properties(
    const ThermochemistryRecord&, double temperature, double molar_density,
    const std::vector<double>& mole_fractions
);
struct DensityStandardState final { double gibbs = 0.0, enthalpy = 0.0, entropy = 0.0; };
DensityStandardState pressure_to_density_standard(
    double gibbs, double enthalpy, double entropy, double temperature,
    double pressure, double density
);

struct ActivitySet final { std::vector<double> values; std::string standard_state_id; };
ActivitySet direct_activities(const std::vector<double>&, const std::string& standard_state_id);
ActivitySet mole_fraction_activities(
    const std::vector<double>&, const std::vector<double>&, const std::string& standard_state_id
);
ActivitySet molarity_activities(
    const std::vector<double>&, double molar_density, const std::vector<double>&,
    double standard_concentration_mol_per_m3, const std::string& standard_state_id
);
ActivitySet molality_activities(
    const std::vector<double>&, const std::vector<double>& amounts_mol, double solvent_mass_kg,
    double standard_molality_mol_per_kg, const std::string& standard_state_id
);

struct ReactionLogPolynomial final {
    double a = 0.0, b = 0.0, c = 0.0, d = 0.0;
    double reference_temperature = 0.0, temperature_min = 0.0, temperature_max = 0.0;
    std::string standard_state_id;
};
struct ReactionProperties final {
    double ln_k = 0.0, derivative = 0.0, enthalpy = 0.0;
    std::optional<double> k;
    std::string standard_state_id;
};

double reaction_ln_quotient(
    const std::vector<double>&, const ActivitySet&, const std::vector<int>&,
    const ReactionLogPolynomial& correlation
);
ReactionProperties evaluate_reaction_properties(const ReactionLogPolynomial&, double temperature);
ReactionLogPolynomial normalize_legacy_reaction_correlation(
    double a_old, double b, double c, double d, double reference_temperature,
    double temperature_min, double temperature_max, const std::string& standard_state_id
);
ReactionLogPolynomial normalize_log10_reciprocal_correlation(
    double a, double b, double c, double temperature_min, double temperature_max,
    const std::string& standard_state_id, double reference_temperature
);
ReactionProperties reaction_pressure_to_density(
    const ReactionProperties&, const std::vector<double>& stoichiometry,
    const std::string& density_standard_state_id, double temperature, double pressure, double density
);
ReactionProperties reaction_standard_state_conversion(
    const ReactionProperties&, const std::vector<double>& stoichiometry,
    const std::vector<double>& old_standard_gibbs_j_per_mol,
    const std::vector<double>& new_standard_gibbs_j_per_mol,
    const std::vector<double>& old_standard_enthalpy_j_per_mol,
    const std::vector<double>& new_standard_enthalpy_j_per_mol,
    const std::string& new_standard_state_id, double temperature
);

}  // namespace epcsaft
