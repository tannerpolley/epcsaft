#pragma once

#include "model.hpp"

#include <cppad/cppad.hpp>

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <type_traits>
#include <vector>

namespace epcsaft::detail {

inline constexpr double ELEMENTARY_CHARGE = 1.602176634e-19;
inline constexpr double VACUUM_PERMITTIVITY = 8.8541878128e-12;
inline constexpr double BOLTZMANN = 1.380649e-23;
inline constexpr double ELECTROLYTE_AVOGADRO = 6.02214076e23;
inline constexpr double ELECTROLYTE_PI =
    3.141592653589793238462643383279502884;

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
Scalar figiel_bulk_relative_permittivity(
    const Scalar& salt_free_permittivity,
    const Scalar& ion_fraction,
    const Scalar& suppression_coefficient
) {
    // EqID: figiel_bulk_permittivity
    return salt_free_permittivity
        / (Scalar(1.0) + suppression_coefficient * ion_fraction);
}

template <typename Scalar>
// EqID: figiel_bulk_permittivity
Scalar figiel_bulk_relative_permittivity_derivative(
    const Scalar& salt_free_permittivity,
    const Scalar& ion_fraction,
    const Scalar& suppression_coefficient
) {
    const Scalar denominator =
        Scalar(1.0) + suppression_coefficient * ion_fraction;
    return -salt_free_permittivity * ion_fraction
        / (denominator * denominator);
}

template <typename Scalar>
Scalar bulk_relative_permittivity(
    const Scalar& temperature,
    const Scalar& molar_density,
    const std::vector<Scalar>& mole_fractions,
    const NativeParameters& parameters,
    const std::vector<Scalar>* relative_permittivity_override = nullptr
) {
    Scalar solvent_mass = Scalar(0.0);
    Scalar dielectric_mass = Scalar(0.0);
    Scalar solvent_fraction = Scalar(0.0);
    Scalar ion_fraction = Scalar(0.0);
    Scalar ion_dielectric = Scalar(0.0);
    Scalar schreckenberg_volume = Scalar(0.0);
    Scalar schreckenberg_temperature = Scalar(0.0);
    Scalar zuber_suppression = Scalar(0.0);
    Scalar rueben_effective_density = Scalar(0.0);
    Scalar rueben_correlation = Scalar(0.0);

    const Scalar number_density_angstrom3 =
        molar_density * ELECTROLYTE_AVOGADRO * 1.0e-30;
    const auto relative_permittivity = [&](std::size_t component) {
        return relative_permittivity_override == nullptr
            ? Scalar(parameters.relative_permittivity[component])
            : (*relative_permittivity_override)[component];
    };
    for (std::size_t component = 0; component < mole_fractions.size(); ++component) {
        const Scalar fraction = mole_fractions[component];
        const int charge = parameters.charge_number[component];
        if (charge == 0) {
            solvent_fraction += fraction;
            const Scalar mass = fraction * parameters.molar_mass_kg_per_mol[component];
            solvent_mass += mass;
            dielectric_mass += mass * relative_permittivity(component);
            schreckenberg_volume += fraction
                * parameters.schreckenberg_dielectric_volume_m3_per_mol[component];
            schreckenberg_temperature += fraction
                * parameters.schreckenberg_dielectric_temperature_k[component];
        } else {
            ion_fraction += fraction;
            ion_dielectric += fraction * relative_permittivity(component);
            zuber_suppression += fraction
                * parameters.zuber_ion_suppression_coefficient[component];
        }
        // EqID: rueben_effective_dipolar_density
        rueben_effective_density += number_density_angstrom3 * fraction
            * (
                parameters.rueben_dipole_scaling_debye2[component]
                    * 1.0e-19 / (BOLTZMANN * temperature)
                + 3.0 * parameters.rueben_polarizability_scaling_angstrom3[component]
            );
        // EqID: rueben_correlation_integral
        rueben_correlation += fraction
            * parameters.rueben_correlation_integral_parameter[component];
    }

    if (parameters.permittivity_model == PERMITTIVITY_MODEL_RUEBEN_PERTURBATION) {
        // EqID: rueben_effective_dipolar_density
        const Scalar y_effective = (4.0 * ELECTROLYTE_PI / 9.0)
            * rueben_effective_density;
        // EqID: rueben_correlation_integral
        const Scalar integral = Scalar(1.0) + rueben_correlation
            * (electrolyte_exp(-y_effective) - Scalar(1.0));
        // EqID: rueben_bulk_permittivity
        return Scalar(1.0) + Scalar(3.0) * y_effective
            * (
                Scalar(1.0) + y_effective
                + ((17.0 / 16.0) * integral - Scalar(1.0))
                    * y_effective * y_effective
            );
    }
    if (parameters.permittivity_model == PERMITTIVITY_MODEL_SCHRECKENBERG_DENSITY) {
        // EqID: schreckenberg_mixed_parameters
        const Scalar mixed_volume = schreckenberg_volume / solvent_fraction;
        const Scalar mixed_temperature = schreckenberg_temperature / solvent_fraction;
        // EqID: schreckenberg_bulk_permittivity
        return Scalar(1.0) + molar_density * solvent_fraction * mixed_volume
            * (mixed_temperature / temperature - Scalar(1.0));
    }

    // EqID: salt_free_solvent_permittivity
    const Scalar salt_free_permittivity = dielectric_mass / solvent_mass;
    if (parameters.permittivity_model == PERMITTIVITY_MODEL_ASCANI_MOLE_FRACTION) {
        // EqID: ascani_bulk_permittivity
        return salt_free_permittivity * solvent_fraction + ion_dielectric;
    }
    if (parameters.permittivity_model == PERMITTIVITY_MODEL_HELD_WATER) {
        // EqID: salt_free_solvent_permittivity
        return salt_free_permittivity;
    }
    if (parameters.permittivity_model == PERMITTIVITY_MODEL_ZUBER_ION_SPECIFIC) {
        // EqID: zuber_bulk_permittivity
        return salt_free_permittivity / (Scalar(1.0) + zuber_suppression);
    }
    if (parameters.permittivity_model == PERMITTIVITY_MODEL_FIGIEL_COMPOSITION) {
        return figiel_bulk_relative_permittivity(
            salt_free_permittivity,
            ion_fraction,
            Scalar(parameters.dielectric_ion_suppression_coefficient)
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
    const std::vector<Scalar>& diameter,
    const std::vector<double>& reference_diameter,
    const std::vector<Scalar>* relative_permittivity_override = nullptr
) {
    bool has_ion = false;
    for (int charge : parameters.charge_number) {
        has_ion = has_ion || charge != 0;
    }
    if (!has_ion) {
        return {0.0, Scalar(0.0)};
    }
    Scalar solvent_mass = Scalar(0.0);
    Scalar charge_square_sum = Scalar(0.0);
    for (std::size_t component = 0; component < mole_fractions.size(); ++component) {
        const int charge = parameters.charge_number[component];
        if (charge == 0) {
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
        relative_permittivity_override
    );
    if constexpr (std::is_same_v<Scalar, double>) {
        if (!(solvent_mass > 0.0) || !std::isfinite(solvent_mass)) {
            throw std::invalid_argument("electrolyte state has no solvent mass");
        }
        if (!(bulk_permittivity > 0.0) || !std::isfinite(bulk_permittivity)) {
            throw std::invalid_argument("bulk relative permittivity is nonpositive");
        }
        if (charge_square_sum == 0.0) {
            return {bulk_permittivity, Scalar(0.0)};
        }
    }

    double reference_charge_square_sum = 0.0;
    for (std::size_t component = 0;
         component < reference_mole_fractions.size(); ++component) {
        const int charge = parameters.charge_number[component];
        if (charge != 0) {
            reference_charge_square_sum += reference_mole_fractions[component]
                * charge * charge;
        }
    }
    const std::vector<Scalar> reference_fractions(
        reference_mole_fractions.begin(), reference_mole_fractions.end()
    );
    const Scalar reference_permittivity = detail::bulk_relative_permittivity(
        Scalar(reference_temperature),
        Scalar(reference_molar_density),
        reference_fractions,
        parameters,
        relative_permittivity_override
    );
    if (reference_charge_square_sum == 0.0) {
        return {bulk_permittivity, Scalar(0.0)};
    }
    // EqID: debye_screening_parameter
    const Scalar kappa = electrolyte_sqrt(
        molar_density * ELECTROLYTE_AVOGADRO
        * ELEMENTARY_CHARGE * ELEMENTARY_CHARGE * charge_square_sum
        / (VACUUM_PERMITTIVITY * bulk_permittivity * BOLTZMANN * temperature)
    );
    const Scalar reference_kappa = electrolyte_sqrt(
        reference_molar_density * ELECTROLYTE_AVOGADRO
        * ELEMENTARY_CHARGE * ELEMENTARY_CHARGE * reference_charge_square_sum
        / (VACUUM_PERMITTIVITY * reference_permittivity * BOLTZMANN
            * reference_temperature)
    );
    const double reference_permittivity_for_branch =
        detail::bulk_relative_permittivity(
            reference_temperature,
            reference_molar_density,
            reference_mole_fractions,
            parameters
        );
    const double reference_kappa_for_branch = std::sqrt(
        reference_molar_density * ELECTROLYTE_AVOGADRO
        * ELEMENTARY_CHARGE * ELEMENTARY_CHARGE * reference_charge_square_sum
        / (VACUUM_PERMITTIVITY * reference_permittivity_for_branch
            * BOLTZMANN * reference_temperature)
    );
    Scalar weighted_chi = Scalar(0.0);
    for (std::size_t component = 0; component < mole_fractions.size(); ++component) {
        const int charge = parameters.charge_number[component];
        if (charge == 0) {
            continue;
        }
        const Scalar diameter_m = Scalar(1.0e-10) * diameter[component];
        const Scalar reduced = kappa * diameter_m;
        const Scalar chi = finite_size_chi(
            reduced,
            reference_kappa_for_branch
                    * 1.0e-10 * reference_diameter[component]
                < 1.0e-3
        );
        weighted_chi += mole_fractions[component] * charge * charge * chi;
    }
    // EqID: debye_huckel_helmholtz
    const Scalar debye_huckel = -kappa * ELEMENTARY_CHARGE * ELEMENTARY_CHARGE
        * weighted_chi
        / (12.0 * ELECTROLYTE_PI * VACUUM_PERMITTIVITY * bulk_permittivity
           * BOLTZMANN * temperature);
    return {bulk_permittivity, debye_huckel};
}

}  // namespace epcsaft::detail
