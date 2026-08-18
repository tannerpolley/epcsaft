#pragma once

#include "model.hpp"
#include "physical_constants.hpp"

#include <cppad/cppad.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <type_traits>
#include <vector>

namespace epcsaft::detail {

inline constexpr double ELECTROLYTE_PI =
    3.141592653589793238462643383279502884;

inline bool dielectric_solvent_component(
    std::size_t component,
    const NativeParameters& parameters
) {
    if (parameters.charge_number[component] != 0) {
        return false;
    }
    if (parameters.relative_permittivity_formulation ==
            RELATIVE_PERMITTIVITY_FORMULATION_DENSITY_VOLUME_TEMPERATURE ||
        parameters.relative_permittivity_formulation ==
            RELATIVE_PERMITTIVITY_FORMULATION_DIPOLAR_DENSITY_PERTURBATION) {
        return true;
    }
    if (parameters.relative_permittivity[component] > 0.0) {
        return true;
    }
    return std::find(
        parameters.permittivity_term_component.begin(),
        parameters.permittivity_term_component.end(),
        static_cast<int>(component)
    ) != parameters.permittivity_term_component.end();
}

template <typename Scalar>
// EqID: component_permittivity_temperature_polynomial
std::vector<Scalar> component_relative_permittivities(
    const Scalar& temperature,
    const NativeParameters& parameters
) {
    std::vector<Scalar> values(
        parameters.relative_permittivity.begin(),
        parameters.relative_permittivity.end()
    );
    for (std::size_t term = 0;
         term < parameters.permittivity_term_component.size(); ++term) {
        const std::size_t component = static_cast<std::size_t>(
            parameters.permittivity_term_component[term]
        );
        const Scalar delta = temperature
            - Scalar(parameters.permittivity_reference_temperature_k[component]);
        Scalar temperature_power = Scalar(1.0);
        for (int exponent = 0;
             exponent < parameters.permittivity_term_power[term]; ++exponent) {
            temperature_power *= delta;
        }
        values[component] += Scalar(parameters.permittivity_term_coefficient[term])
            * temperature_power;
    }
    return values;
}

template <typename Scalar>
// EqExclude: scalar utility dispatch
Scalar electrolyte_sqrt(const Scalar& value) {
    if constexpr (std::is_same_v<Scalar, double>) {
        return std::sqrt(value);
    }
    return CppAD::sqrt(value);
}

template <typename Scalar>
// EqExclude: scalar utility dispatch
Scalar electrolyte_log(const Scalar& value) {
    if constexpr (std::is_same_v<Scalar, double>) {
        return std::log(value);
    }
    return CppAD::log(value);
}

template <typename Scalar>
// EqExclude: scalar utility dispatch
Scalar electrolyte_exp(const Scalar& value) {
    if constexpr (std::is_same_v<Scalar, double>) {
        return std::exp(value);
    }
    return CppAD::exp(value);
}

template <typename Scalar>
struct ElectrolyteResult final {
    Scalar bulk_relative_permittivity;
    Scalar debye_huckel;
};

template <typename Scalar>
Scalar ion_fraction_suppressed_relative_permittivity(
    const Scalar& salt_free_permittivity,
    const Scalar& ion_fraction,
    const Scalar& suppression_coefficient
) {
    // EqID: ion_fraction_suppressed_bulk_permittivity
    return salt_free_permittivity
        / (Scalar(1.0) + suppression_coefficient * ion_fraction);
}

template <typename Scalar>
Scalar bulk_relative_permittivity(
    const Scalar& temperature,
    const Scalar& molar_density,
    const std::vector<Scalar>& mole_fractions,
    const NativeParameters& parameters,
    const std::vector<Scalar>* relative_permittivity_override = nullptr,
    const Scalar* ion_fraction_suppression_override = nullptr
) {
    const std::vector<Scalar> correlated_relative_permittivity =
        relative_permittivity_override == nullptr
            ? component_relative_permittivities(temperature, parameters)
            : std::vector<Scalar>{};
    const std::vector<Scalar>* evaluated_relative_permittivity =
        relative_permittivity_override == nullptr
            ? &correlated_relative_permittivity
            : relative_permittivity_override;
    Scalar solvent_mass = Scalar(0.0);
    Scalar dielectric_mass = Scalar(0.0);
    Scalar solvent_fraction = Scalar(0.0);
    Scalar ion_fraction = Scalar(0.0);
    Scalar ion_dielectric = Scalar(0.0);
    Scalar dielectric_volume = Scalar(0.0);
    Scalar dielectric_temperature = Scalar(0.0);
    Scalar ion_specific_suppression = Scalar(0.0);
    Scalar effective_dipolar_density = Scalar(0.0);
    Scalar dipolar_correlation = Scalar(0.0);

    const Scalar number_density_angstrom3 =
        molar_density * PHYSICAL_AVOGADRO * 1.0e-30;
    const auto relative_permittivity = [&](std::size_t component) {
        return evaluated_relative_permittivity == nullptr
            ? Scalar(parameters.relative_permittivity[component])
            : (*evaluated_relative_permittivity)[component];
    };
    for (std::size_t component = 0; component < mole_fractions.size(); ++component) {
        const Scalar& fraction = mole_fractions[component];
        const int charge = parameters.charge_number[component];
        if (dielectric_solvent_component(component, parameters)) {
            solvent_fraction += fraction;
            const Scalar mass = fraction * parameters.molar_mass_kg_per_mol[component];
            solvent_mass += mass;
            dielectric_mass += mass * relative_permittivity(component);
            dielectric_volume += fraction
                * parameters.dielectric_volume_m3_per_mol[component];
            dielectric_temperature += fraction
                * parameters.dielectric_temperature_k[component];
        } else if (charge != 0) {
            ion_fraction += fraction;
            ion_dielectric += fraction * relative_permittivity(component);
            ion_specific_suppression += fraction
                * parameters.ion_specific_suppression_coefficient[component];
        }
        // EqID: effective_dipolar_density
        effective_dipolar_density += number_density_angstrom3 * fraction
            * (
                parameters.dipole_scaling_debye2[component]
                    * 1.0e-19 / (PHYSICAL_BOLTZMANN * temperature)
                + 3.0 * parameters.polarizability_scaling_angstrom3[component]
            );
        // EqID: dipolar_correlation_integral
        dipolar_correlation += fraction
            * parameters.correlation_integral_parameter[component];
    }

    if (parameters.relative_permittivity_formulation == RELATIVE_PERMITTIVITY_FORMULATION_DIPOLAR_DENSITY_PERTURBATION) {
        // EqID: effective_dipolar_density
        const Scalar y_effective = (4.0 * ELECTROLYTE_PI / 9.0)
            * effective_dipolar_density;
        // EqID: dipolar_correlation_integral
        const Scalar integral = Scalar(1.0) + dipolar_correlation
            * (electrolyte_exp(-y_effective) - Scalar(1.0));
        // EqID: dipolar_perturbation_bulk_permittivity
        return Scalar(1.0) + Scalar(3.0) * y_effective
            * (
                Scalar(1.0) + y_effective
                + ((17.0 / 16.0) * integral - Scalar(1.0))
                    * y_effective * y_effective
            );
    }
    if (parameters.relative_permittivity_formulation == RELATIVE_PERMITTIVITY_FORMULATION_DENSITY_VOLUME_TEMPERATURE) {
        // EqID: density_correlation_mixed_parameters
        const Scalar mixed_volume = dielectric_volume / solvent_fraction;
        const Scalar mixed_temperature = dielectric_temperature / solvent_fraction;
        // EqID: density_correlated_bulk_permittivity
        return Scalar(1.0) + molar_density * solvent_fraction * mixed_volume
            * (mixed_temperature / temperature - Scalar(1.0));
    }

    if (parameters.relative_permittivity_formulation == RELATIVE_PERMITTIVITY_FORMULATION_COMPONENT_MIXING) {
        const bool has_neutral_component = std::any_of(
            parameters.charge_number.begin(),
            parameters.charge_number.end(),
            [](int charge) { return charge == 0; }
        );
        if (!has_neutral_component) {
            // EqID: component_mixed_bulk_permittivity
            return ion_dielectric;
        }
    }

    // EqID: salt_free_solvent_permittivity
    const Scalar salt_free_permittivity = dielectric_mass / solvent_mass;
    if (parameters.relative_permittivity_formulation == RELATIVE_PERMITTIVITY_FORMULATION_COMPONENT_MIXING) {
        // EqID: component_mixed_bulk_permittivity
        return salt_free_permittivity * solvent_fraction + ion_dielectric;
    }
    if (parameters.relative_permittivity_formulation == RELATIVE_PERMITTIVITY_FORMULATION_SOLVENT_ONLY) {
        // EqID: salt_free_solvent_permittivity
        return salt_free_permittivity;
    }
    if (parameters.relative_permittivity_formulation == RELATIVE_PERMITTIVITY_FORMULATION_ION_SPECIFIC_SUPPRESSION) {
        // EqID: ion_specific_suppressed_bulk_permittivity
        return salt_free_permittivity / (Scalar(1.0) + ion_specific_suppression);
    }
    if (parameters.relative_permittivity_formulation == RELATIVE_PERMITTIVITY_FORMULATION_ION_FRACTION_SUPPRESSION) {
        const Scalar suppression_coefficient =
            ion_fraction_suppression_override == nullptr
                ? Scalar(parameters.ion_fraction_suppression_coefficient)
                : *ion_fraction_suppression_override;
        return ion_fraction_suppressed_relative_permittivity(
            salt_free_permittivity,
            ion_fraction,
            suppression_coefficient
        );
    }
    throw std::invalid_argument("invalid relative permittivity model");
}

template <typename Scalar>
// EqID: debye_finite_size_chi
Scalar finite_size_chi(const Scalar& reduced, bool use_series) {
    if (use_series) {
        const Scalar reduced_2 = reduced * reduced;
        const Scalar reduced_3 = reduced_2 * reduced;
        const Scalar reduced_4 = reduced_3 * reduced;
        const Scalar reduced_5 = reduced_4 * reduced;
        const Scalar reduced_6 = reduced_5 * reduced;
        return Scalar(1.0) - 0.75 * reduced + 0.6 * reduced_2
            - 0.5 * reduced_3 + (3.0 / 7.0) * reduced_4
            - 0.375 * reduced_5 + (1.0 / 3.0) * reduced_6;
    }
    return 3.0
        * (electrolyte_log(Scalar(1.0) + reduced) - reduced
           + 0.5 * reduced * reduced)
        / (reduced * reduced * reduced);
}

template <typename Scalar>
ElectrolyteResult<Scalar> electrolyte_contribution(
    const Scalar& temperature,
    const Scalar& molar_density,
    double reference_temperature,
    double reference_molar_density,
    const std::vector<Scalar>& mole_fractions,
    const std::vector<double>& reference_mole_fractions,
    const NativeParameters& parameters,
    const std::vector<Scalar>* relative_permittivity_override = nullptr,
    const Scalar* ion_fraction_suppression_override = nullptr
) {
    bool has_ion = false;
    for (int charge : parameters.charge_number) {
        has_ion = has_ion || charge != 0;
    }
    if (!has_ion) {
        return {0.0, Scalar(0.0)};
    }
    double reference_charge_square_sum = 0.0;
    double reference_solvent_mass = 0.0;
    for (std::size_t component = 0;
         component < reference_mole_fractions.size(); ++component) {
        const int charge = parameters.charge_number[component];
        if (dielectric_solvent_component(component, parameters)) {
            reference_solvent_mass += reference_mole_fractions[component]
                * parameters.molar_mass_kg_per_mol[component];
        } else if (charge != 0) {
            reference_charge_square_sum += reference_mole_fractions[component]
                * charge * charge;
        }
    }
    if (reference_charge_square_sum == 0.0) {
        if (reference_solvent_mass == 0.0) {
            return {0.0, Scalar(0.0)};
        }
        const Scalar neutral_permittivity = detail::bulk_relative_permittivity(
            temperature,
            molar_density,
            mole_fractions,
            parameters,
            relative_permittivity_override,
            ion_fraction_suppression_override
        );
        if constexpr (std::is_same_v<Scalar, double>) {
            if (!(neutral_permittivity > 0.0)
                || !std::isfinite(neutral_permittivity)) {
                throw std::invalid_argument(
                    "bulk relative permittivity is nonpositive"
                );
            }
        }
        return {neutral_permittivity, Scalar(0.0)};
    }

    Scalar solvent_mass = Scalar(0.0);
    Scalar charge_square_sum = Scalar(0.0);
    for (std::size_t component = 0; component < mole_fractions.size(); ++component) {
        const int charge = parameters.charge_number[component];
        if (dielectric_solvent_component(component, parameters)) {
            const Scalar mass = mole_fractions[component]
                * parameters.molar_mass_kg_per_mol[component];
            solvent_mass += mass;
        } else {
            // EqID: ionic_strength_charge_sum
            charge_square_sum += mole_fractions[component] * charge * charge;
        }
    }
    const Scalar bulk_permittivity = detail::bulk_relative_permittivity(
        temperature,
        molar_density,
        mole_fractions,
        parameters,
        relative_permittivity_override,
        ion_fraction_suppression_override
    );
    if constexpr (std::is_same_v<Scalar, double>) {
        if (parameters.relative_permittivity_formulation !=
                RELATIVE_PERMITTIVITY_FORMULATION_COMPONENT_MIXING &&
            (!(solvent_mass > 0.0) || !std::isfinite(solvent_mass))) {
            throw std::invalid_argument("electrolyte state has no solvent mass");
        }
        if (!(bulk_permittivity > 0.0) || !std::isfinite(bulk_permittivity)) {
            throw std::invalid_argument("bulk relative permittivity is nonpositive");
        }
        if (charge_square_sum == 0.0) {
            return {bulk_permittivity, Scalar(0.0)};
        }
    }

    // EqID: debye_screening_parameter
    const Scalar kappa = electrolyte_sqrt(
        molar_density * PHYSICAL_AVOGADRO
        * PHYSICAL_ELEMENTARY_CHARGE * PHYSICAL_ELEMENTARY_CHARGE * charge_square_sum
        / (PHYSICAL_VACUUM_PERMITTIVITY * bulk_permittivity
           * PHYSICAL_BOLTZMANN * temperature)
    );
    const double reference_permittivity_for_branch =
        detail::bulk_relative_permittivity(
            reference_temperature,
            reference_molar_density,
            reference_mole_fractions,
            parameters
        );
    const double reference_kappa_for_branch = std::sqrt(
        reference_molar_density * PHYSICAL_AVOGADRO
        * PHYSICAL_ELEMENTARY_CHARGE * PHYSICAL_ELEMENTARY_CHARGE
        * reference_charge_square_sum
        / (PHYSICAL_VACUUM_PERMITTIVITY * reference_permittivity_for_branch
            * PHYSICAL_BOLTZMANN * reference_temperature)
    );
    Scalar weighted_chi = Scalar(0.0);
    for (std::size_t component = 0; component < mole_fractions.size(); ++component) {
        const int charge = parameters.charge_number[component];
        if (charge == 0) {
            continue;
        }
        const double diameter_angstrom =
            parameters.debye_huckel_diameter_angstrom[component];
        const Scalar diameter_m = Scalar(1.0e-10 * diameter_angstrom);
        const Scalar reduced = kappa * diameter_m;
        const Scalar chi = finite_size_chi(
            reduced,
            reference_kappa_for_branch * 1.0e-10 * diameter_angstrom
                < 1.0e-3
        );
        weighted_chi += mole_fractions[component] * charge * charge * chi;
    }
    // EqID: debye_huckel_helmholtz
    const Scalar debye_huckel = -kappa * PHYSICAL_ELEMENTARY_CHARGE
        * PHYSICAL_ELEMENTARY_CHARGE * weighted_chi
        / (12.0 * ELECTROLYTE_PI * PHYSICAL_VACUUM_PERMITTIVITY
           * bulk_permittivity * PHYSICAL_BOLTZMANN * temperature);
    return {bulk_permittivity, debye_huckel};
}

}  // namespace epcsaft::detail
