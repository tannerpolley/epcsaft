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
        // EqID: ideal_heat_capacity_cp
        result.isobaric_heat_capacity += mole_fractions[i] * thermal.cp;
        if (supported[i]) result.entropy += mole_fractions[i] * (thermal.s
            - gas_constant * scalar_log(mole_fractions[i] * molar_density
                * gas_constant * temperature / Scalar(pressure)));
    }
    result.internal_energy = result.enthalpy - gas_constant * temperature;
    // EqID: ideal_heat_capacity_cv
    result.isochoric_heat_capacity = result.isobaric_heat_capacity - gas_constant;
    return result;
}

void validate_thermochemistry(const ThermochemistryRecord&, std::size_t component_count);
IdealPropertiesEvaluation evaluate_ideal_properties(
    const ThermochemistryRecord&, double temperature, double molar_density,
    const std::vector<double>& mole_fractions
);
struct ReactionLogPolynomial final {
    double a = 0.0, b = 0.0, c = 0.0, d = 0.0;
    double reference_temperature = 0.0, temperature_min = 0.0, temperature_max = 0.0;
    std::string standard_state_id;
    double e = 0.0, f = 0.0;
};
// One scalar/AD expression owns the source law and its raw thermal derivatives.
template <typename Scalar>
// EqID: reaction_constant_ln_polynomial
Scalar reaction_log_value(const ReactionLogPolynomial& law, const Scalar& temperature,
                          std::size_t thermal_order = 0) {
    Scalar value;
    if (thermal_order == 0) {
        value = law.a + law.b / temperature + law.c * scalar_log(temperature / law.reference_temperature)
            + law.d * temperature;
        if (law.e == 0.0 && law.f == 0.0) return value;
        return value + law.e / (temperature * temperature) + law.f * temperature * temperature;
    }
    const Scalar t2 = temperature * temperature;
    if (thermal_order == 1) {
        value = -law.b / t2 + law.c / temperature + law.d;
        if (law.e == 0.0 && law.f == 0.0) return value;
        return value - 2.0 * law.e / (t2 * temperature) + 2.0 * law.f * temperature;
    }
    if (thermal_order == 2) {
        value = 2.0 * law.b / (t2 * temperature) - law.c / t2;
        if (law.e == 0.0 && law.f == 0.0) return value;
        return value + 6.0 * law.e / (t2 * t2) + 2.0 * law.f;
    }
    if (thermal_order == 3) {
        value = -6.0 * law.b / (t2 * t2) + 2.0 * law.c / (t2 * temperature);
        if (law.e == 0.0) return value;
        return value - 24.0 * law.e / (t2 * t2 * temperature);
    }
    throw std::invalid_argument("source thermal derivative order");
}
struct ReactionProperties final {
    // derivative, second_derivative, third_derivative are d^k ln K / dT^k.
    double ln_k = 0.0, derivative = 0.0, enthalpy = 0.0;
    double second_derivative = 0.0, third_derivative = 0.0;
    std::optional<double> k;
    std::string standard_state_id;
};
// Reaction Gibbs/enthalpy/heat capacity of supplied species thermochemistry on the
// EOS density basis (1 mol/m^3); absent when a participating interval excludes T.
struct ReactionThermochemistry final {
    double gibbs = NAN, enthalpy = NAN, heat_capacity = NAN, pressure_gibbs = NAN;
    double ln_k_pressure = NAN, ln_k_density = NAN;
    std::optional<double> k_pressure, k_density;
    std::array<double, 4> gibbs_derivatives = {NAN, NAN, NAN, NAN};
    ThermochemistryStatus status = ThermochemistryStatus::MissingRecord;
    std::string cause;
};
ReactionThermochemistry reaction_thermochemistry_density(
    const ThermochemistryRecord&, const std::vector<double>& stoichiometry, double temperature
);

ReactionProperties evaluate_reaction_properties(const ReactionLogPolynomial&, double temperature);
}  // namespace epcsaft
