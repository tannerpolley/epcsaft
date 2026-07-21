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
Scalar electrolyte_sqrt(const Scalar& value) {
    if constexpr (std::is_same_v<Scalar, double>) {
        return std::sqrt(value);
    }
    return CppAD::sqrt(value);
}

template <typename Scalar>
Scalar electrolyte_log(const Scalar& value) {
    if constexpr (std::is_same_v<Scalar, double>) {
        return std::log(value);
    }
    return CppAD::log(value);
}

template <typename Scalar>
struct ElectrolyteResult final {
    Scalar bulk_relative_permittivity;
    Scalar debye_huckel;
};

template <typename Scalar>
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
    const std::vector<double>& reference_diameter
) {
    bool has_ion = false;
    for (int charge : parameters.charge_number) {
        has_ion = has_ion || charge != 0;
    }
    if (!has_ion) {
        return {0.0, Scalar(0.0)};
    }
    Scalar solvent_mass = Scalar(0.0);
    Scalar dielectric_mass = Scalar(0.0);
    Scalar solvent_fraction = Scalar(0.0);
    Scalar ion_fraction = Scalar(0.0);
    Scalar ion_dielectric = Scalar(0.0);
    Scalar charge_square_sum = Scalar(0.0);
    for (std::size_t component = 0; component < mole_fractions.size(); ++component) {
        const int charge = parameters.charge_number[component];
        if (charge == 0) {
            solvent_fraction += mole_fractions[component];
            const Scalar mass = mole_fractions[component]
                * parameters.molar_mass_kg_per_mol[component];
            solvent_mass += mass;
            dielectric_mass += mass * parameters.relative_permittivity[component];
        } else {
            ion_fraction += mole_fractions[component];
            ion_dielectric += mole_fractions[component]
                * parameters.relative_permittivity[component];
            charge_square_sum += mole_fractions[component] * charge * charge;
        }
    }
    const Scalar salt_free_permittivity = dielectric_mass / solvent_mass;
    const Scalar bulk_permittivity = parameters.electrolyte_model
            == ELECTROLYTE_MODEL_ASCANI_ORIGINAL_BORN
        ? salt_free_permittivity * solvent_fraction + ion_dielectric
        : salt_free_permittivity /
            (Scalar(1.0) + parameters.dielectric_ion_suppression_coefficient
                * ion_fraction);
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

    double reference_solvent_mass = 0.0;
    double reference_dielectric_mass = 0.0;
    double reference_solvent_fraction = 0.0;
    double reference_ion_fraction = 0.0;
    double reference_ion_dielectric = 0.0;
    double reference_charge_square_sum = 0.0;
    for (std::size_t component = 0;
         component < reference_mole_fractions.size(); ++component) {
        const int charge = parameters.charge_number[component];
        if (charge == 0) {
            reference_solvent_fraction += reference_mole_fractions[component];
            const double mass = reference_mole_fractions[component]
                * parameters.molar_mass_kg_per_mol[component];
            reference_solvent_mass += mass;
            reference_dielectric_mass += mass
                * parameters.relative_permittivity[component];
        } else {
            reference_ion_fraction += reference_mole_fractions[component];
            reference_ion_dielectric += reference_mole_fractions[component]
                * parameters.relative_permittivity[component];
            reference_charge_square_sum += reference_mole_fractions[component]
                * charge * charge;
        }
    }
    const double reference_salt_free_permittivity =
        reference_dielectric_mass / reference_solvent_mass;
    const double reference_permittivity = parameters.electrolyte_model
            == ELECTROLYTE_MODEL_ASCANI_ORIGINAL_BORN
        ? reference_salt_free_permittivity * reference_solvent_fraction
            + reference_ion_dielectric
        : reference_salt_free_permittivity
            / (1.0 + parameters.dielectric_ion_suppression_coefficient
                * reference_ion_fraction);
    if (reference_charge_square_sum == 0.0) {
        return {bulk_permittivity, Scalar(0.0)};
    }
    const Scalar kappa = electrolyte_sqrt(
        molar_density * ELECTROLYTE_AVOGADRO
        * ELEMENTARY_CHARGE * ELEMENTARY_CHARGE * charge_square_sum
        / (VACUUM_PERMITTIVITY * bulk_permittivity * BOLTZMANN * temperature)
    );
    const double reference_kappa = std::sqrt(
        reference_molar_density * ELECTROLYTE_AVOGADRO
        * ELEMENTARY_CHARGE * ELEMENTARY_CHARGE * reference_charge_square_sum
        / (VACUUM_PERMITTIVITY * reference_permittivity * BOLTZMANN
            * reference_temperature)
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
            reference_kappa * 1.0e-10 * reference_diameter[component] < 1.0e-3
        );
        weighted_chi += mole_fractions[component] * charge * charge * chi;
    }
    const Scalar debye_huckel = -kappa * ELEMENTARY_CHARGE * ELEMENTARY_CHARGE
        * weighted_chi
        / (12.0 * ELECTROLYTE_PI * VACUUM_PERMITTIVITY * bulk_permittivity
           * BOLTZMANN * temperature);
    return {bulk_permittivity, debye_huckel};
}

}  // namespace epcsaft::detail
