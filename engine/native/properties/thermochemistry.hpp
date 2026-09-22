#pragma once

#include "contributions/ideal.hpp"

#include <cstddef>
#include <optional>
#include <stdexcept>
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

template <typename Scalar>
struct IdealMixtureValue final {
    Scalar enthalpy = Scalar(0.0);
    Scalar internal_energy = Scalar(0.0);
    Scalar entropy = Scalar(0.0);
    Scalar isobaric_heat_capacity = Scalar(0.0);
    Scalar isochoric_heat_capacity = Scalar(0.0);
};

template <typename Scalar>
IdealMixtureValue<Scalar> evaluate_ideal_mixture(
    const std::vector<IdealCorrelation>& correlations, double pressure,
    const Scalar& temperature, const Scalar& molar_density,
    const std::vector<Scalar>& mole_fractions,
    const std::vector<bool>& supported
) {
    if (correlations.size() != mole_fractions.size()
        || supported.size() != mole_fractions.size())
        throw std::invalid_argument("ideal mixture size");
    IdealMixtureValue<Scalar> result;
    const Scalar gas_constant = Scalar(GAS_CONSTANT_J_PER_MOL_K);
    for (std::size_t i = 0; i < mole_fractions.size(); ++i) {
        const auto thermal = evaluate_ideal_correlation(
            correlations[i], temperature, pressure);
        result.enthalpy += mole_fractions[i] * thermal.h;
        result.isobaric_heat_capacity += mole_fractions[i] * thermal.cp;
        if (supported[i]) result.entropy += mole_fractions[i] * (thermal.s
            - gas_constant * scalar_log(mole_fractions[i] * molar_density
                * gas_constant * temperature / Scalar(pressure)));
    }
    result.internal_energy = result.enthalpy - gas_constant * temperature;
    result.isochoric_heat_capacity = result.isobaric_heat_capacity - gas_constant;
    return result;
}

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
    // derivative, second_derivative, third_derivative are d^k ln K / dT^k.
    double ln_k = 0.0, derivative = 0.0, enthalpy = 0.0;
    double second_derivative = 0.0, third_derivative = 0.0;
    std::optional<double> k;
    std::string standard_state_id;
};
// Reaction Gibbs/enthalpy/heat capacity of supplied species thermochemistry on the
// EOS density basis (1 mol/m^3); absent when a participating interval excludes T.
struct ReactionThermochemistry final { double gibbs = 0.0, enthalpy = 0.0, heat_capacity = 0.0; };
std::optional<ReactionThermochemistry> reaction_thermochemistry_density(
    const ThermochemistryRecord&, const std::vector<double>& stoichiometry, double temperature
);

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
// Exact polynomial image of a standard-state change. Per species, old-minus-new standard
// enthalpy, entropy, heat capacity and heat-capacity slope at T0 give
// Delta g(T)/(RT) = C0/T + C1 + C2 ln(T/T0) + C3 T, which lands in the same ln K family;
// zero heat-capacity terms reproduce the constant-enthalpy shift exactly.
ReactionLogPolynomial reaction_standard_state_shift(
    const ReactionLogPolynomial&, const std::vector<double>& stoichiometry,
    const std::vector<double>& delta_enthalpy_j_per_mol,
    const std::vector<double>& delta_entropy_j_per_mol_k,
    const std::vector<double>& delta_heat_capacity_j_per_mol_k,
    const std::vector<double>& delta_heat_capacity_slope_j_per_mol_k2,
    double reference_temperature, const std::string& new_standard_state_id
);

}  // namespace epcsaft
