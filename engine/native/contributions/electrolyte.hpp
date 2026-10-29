#pragma once

#include "core/numerical_policy.hpp"
#include "core/parameters.hpp"
#include "core/physical_constants.hpp"
#include "core/scalar.hpp"

#include <cppad/cppad.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <type_traits>
#include <vector>

namespace epcsaft {

template <class ParameterSet>
inline bool dielectric_solvent_component(
    std::size_t component,
    const ParameterSet& parameters
) {
    if (parameters.permittivity_form.is(PermittivityForm::SolventMoleWeighting))
        return parameters.permittivity_solvent[component] != 0;
    if (parameters.charge_number[component] != 0) {
        return false;
    }
    if (parameters.permittivity_form.is(PermittivityForm::DensityTemperature)
        || parameters.permittivity_form.is(PermittivityForm::DipolarPerturbation)) {
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

template <typename Scalar, class ParameterSet>
// EqID: component_permittivity_temperature_polynomial
std::vector<Scalar> component_relative_permittivities(
    const Scalar& temperature,
    const ParameterSet& parameters
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
struct ElectrolyteResult final {
    Scalar bulk_relative_permittivity;
    Scalar debye_huckel;
};

template <typename Scalar, class ParameterSet, typename Fixed>
Scalar bulk_relative_permittivity(
    const Fixed& temperature,
    const Scalar& molar_density,
    const std::vector<Fixed>& mole_fractions,
    const ParameterSet& parameters
) {
    const auto relative_permittivity = component_relative_permittivities(temperature, parameters);
    // EqID: solvent_mole_weighted_bulk_permittivity
    if (parameters.permittivity_form.is(PermittivityForm::SolventMoleWeighting)) {
        Fixed numerator = Fixed(0.0), denominator = Fixed(0.0);
        for (std::size_t i = 0; i < mole_fractions.size(); ++i) {
            if (!parameters.permittivity_solvent[i]) continue;
            numerator += mole_fractions[i] * relative_permittivity[i];
            denominator += mole_fractions[i];
        }
        if (parameters.permittivity_solvent_normalized) {
            if constexpr (std::is_same_v<Fixed, double>)
                if (!(denominator > 0.0)) throw std::invalid_argument("permittivity has no solvent denominator");
            return numerator / denominator;
        }
        return numerator;
    }
    Fixed solvent_mass = Fixed(0.0);
    Fixed dielectric_mass = Fixed(0.0);
    Fixed solvent_fraction = Fixed(0.0);
    Fixed ion_dielectric = Fixed(0.0);
    Fixed dielectric_volume = Fixed(0.0);
    Fixed dielectric_temperature = Fixed(0.0);
    Fixed ion_suppression = Fixed(0.0);
    for (std::size_t component = 0; component < mole_fractions.size(); ++component) {
        const Fixed& fraction = mole_fractions[component];
        const int charge = parameters.charge_number[component];
        if (dielectric_solvent_component(component, parameters)) {
            solvent_fraction += fraction;
            const Fixed mass = fraction * parameters.molar_mass_kg_per_mol[component];
            solvent_mass += mass;
            dielectric_mass += mass * relative_permittivity[component];
            dielectric_volume += fraction
                * parameters.dielectric_volume_m3_per_mol[component];
            dielectric_temperature += fraction
                * parameters.dielectric_temperature_k[component];
        } else if (charge != 0) {
            ion_dielectric += fraction * relative_permittivity[component];
            ion_suppression += fraction
                * parameters.ion_suppression_coefficient[component];
        }
    }

    if (parameters.permittivity_form.is(PermittivityForm::DipolarPerturbation)) {
        const Scalar number_density_angstrom3 = molar_density * NUMBER_DENSITY_PER_ANGSTROM3;
        Scalar effective_dipolar_density = Scalar(0.0);
        Fixed dipolar_correlation = Fixed(0.0);
        for (std::size_t component = 0; component < mole_fractions.size(); ++component) {
            const Fixed& fraction = mole_fractions[component];
            // EqID: effective_dipolar_density
            effective_dipolar_density += number_density_angstrom3 * fraction
                * (
                    Fixed(parameters.dipole_scaling_debye2[component])
                        * DEBYE_SQUARED_J_ANGSTROM3 / (BOLTZMANN_J_PER_K * temperature)
                    + Fixed(3.0) * parameters.polarizability_scaling_angstrom3[component]
                );
            // EqID: dipolar_correlation_integral
            dipolar_correlation += fraction
                * parameters.correlation_integral_parameter[component];
        }
        // EqID: effective_dipolar_density
        const Scalar y_effective = (4.0 * Scalar(PI) / 9.0)
            * effective_dipolar_density;
        // EqID: dipolar_correlation_integral
        const Scalar integral = Scalar(1.0) + dipolar_correlation
            * (scalar_exp(-y_effective) - Scalar(1.0));
        // EqID: dipolar_perturbation_bulk_permittivity
        return Scalar(1.0) + Scalar(3.0) * y_effective
            * (
                Scalar(1.0) + y_effective
                + ((17.0 / 16.0) * integral - Scalar(1.0))
                    * y_effective * y_effective
            );
    }
    if (parameters.permittivity_form.is(PermittivityForm::DensityTemperature)) {
        // EqID: density_correlation_mixed_parameters
        const Fixed mixed_volume = dielectric_volume / solvent_fraction;
        const Fixed mixed_temperature = dielectric_temperature / solvent_fraction;
        // EqID: density_correlated_bulk_permittivity
        return Scalar(1.0) + molar_density * solvent_fraction * mixed_volume
            * (mixed_temperature / temperature - Fixed(1.0));
    }

    if (parameters.permittivity_form.is(PermittivityForm::ComponentMixing)) {
        // EqID: component_mixed_bulk_permittivity
        if (std::none_of(parameters.charge_number.begin(), parameters.charge_number.end(),
                         [](int charge) { return charge == 0; }))
            return ion_dielectric;
        if (solvent_mass == Fixed(0.0)) {
            std::size_t solvent = 0, solvent_count = 0;
            for (std::size_t component = 0; component < mole_fractions.size(); ++component) {
                if (dielectric_solvent_component(component, parameters)) {
                    solvent = component;
                    ++solvent_count;
                }
            }
            if (solvent_count == 1)
                return mole_fractions[solvent] * relative_permittivity[solvent] + ion_dielectric;
        }
        // EqID: salt_free_solvent_permittivity
        return dielectric_mass / solvent_mass * solvent_fraction + ion_dielectric;
    }

    // EqID: salt_free_solvent_permittivity
    const Fixed salt_free_permittivity = dielectric_mass / solvent_mass;
    if (parameters.permittivity_form.is(PermittivityForm::IonSuppression)) {
        // EqID: ion_suppressed_bulk_permittivity
        // Covers ion-specific suppression directly, ion-fraction
        // suppression with alpha equal for every ion, and solvent-only
        // with alpha = 0.
        return salt_free_permittivity / (Fixed(1.0) + ion_suppression);
    }
    throw std::invalid_argument("no permittivity form selected");
}

template <typename Scalar>
// EqID: debye_finite_size_chi
Scalar finite_size_chi(const Scalar& reduced, bool use_series) {
    if (use_series) {
        // chi = sum_k 3 (-r)^k / (k + 3), Horner through k = 24 (see SCREENING_SERIES_SWITCH).
        Scalar sum = Scalar(0.0);
        for (int k = 24; k >= 0; --k) sum = sum * (-reduced) + 3.0 / (k + 3);
        return sum;
    }
    return 3.0
        * (scalar_log(Scalar(1.0) + reduced) - reduced
           + 0.5 * reduced * reduced)
        / (reduced * reduced * reduced);
}

template <typename Scalar, class ParameterSet, typename Fixed>
ElectrolyteResult<Scalar> electrolyte_contribution(
    const Fixed& temperature,
    const Scalar& molar_density,
    const std::vector<Fixed>& mole_fractions,
    const std::vector<double>& reference_mole_fractions,
    const ParameterSet& parameters
) {
    bool has_ion = false;
    for (int charge : parameters.charge_number) {
        has_ion = has_ion || charge != 0;
    }
    if (!has_ion) {
        return {Scalar(0.0), Scalar(0.0)};
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
            return {Scalar(0.0), Scalar(0.0)};
        }
        const Scalar neutral_permittivity = bulk_relative_permittivity(
            temperature,
            molar_density,
            mole_fractions,
            parameters
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

    Fixed solvent_mass = Fixed(0.0);
    Fixed charge_square_sum = Fixed(0.0);
    for (std::size_t component = 0; component < mole_fractions.size(); ++component) {
        const int charge = parameters.charge_number[component];
        if (dielectric_solvent_component(component, parameters)) {
            const Fixed mass = mole_fractions[component]
                * parameters.molar_mass_kg_per_mol[component];
            solvent_mass += mass;
        } else {
            // EqID: ionic_strength_charge_sum
            charge_square_sum += mole_fractions[component] * charge * charge;
        }
    }
    const Scalar bulk_permittivity = bulk_relative_permittivity(
        temperature,
        molar_density,
        mole_fractions,
        parameters
    );
    if constexpr (std::is_same_v<Scalar, double>) {
        if (!parameters.permittivity_form.is(PermittivityForm::ComponentMixing) &&
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
    const Scalar kappa = scalar_sqrt(
        molar_density * AVOGADRO_PER_MOL
        * ELEMENTARY_CHARGE_C * ELEMENTARY_CHARGE_C * charge_square_sum
        / (VACUUM_PERMITTIVITY_F_PER_M * bulk_permittivity
           * BOLTZMANN_J_PER_K * temperature)
    );
    Scalar weighted_chi = Scalar(0.0);
    for (std::size_t component = 0; component < mole_fractions.size(); ++component) {
        const int charge = parameters.charge_number[component];
        if (charge == 0) {
            continue;
        }
        const Fixed diameter_m = ANGSTROM_M * parameters.debye_huckel_diameter_angstrom[component];
        const Scalar reduced = kappa * diameter_m;
        const Scalar chi = scalar_select_less(reduced, Scalar(policy::SCREENING_SERIES_SWITCH.value),
            finite_size_chi(reduced, true), finite_size_chi(reduced, false));
        weighted_chi += mole_fractions[component] * charge * charge * chi;
    }
    // EqID: debye_huckel_helmholtz
    const Scalar debye_huckel = -kappa * ELEMENTARY_CHARGE_C
        * ELEMENTARY_CHARGE_C * weighted_chi
        / (12.0 * Scalar(PI) * VACUUM_PERMITTIVITY_F_PER_M
           * bulk_permittivity * BOLTZMANN_J_PER_K * temperature);
    return {bulk_permittivity, debye_huckel};
}

}  // namespace epcsaft
