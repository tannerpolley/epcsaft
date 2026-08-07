#pragma once

#include "electrolyte.hpp"
#include "model.hpp"

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <type_traits>
#include <vector>

namespace epcsaft::detail {

template <
    typename Scalar,
    typename BornDiameter,
    typename SolvationFactor,
    typename IonicPermittivity
>
Scalar born_ssm_ds_contribution_impl(
    const Scalar& temperature,
    const std::vector<Scalar>& mole_fractions,
    const NativeParameters& parameters,
    const Scalar& bulk_relative_permittivity,
    const std::vector<Scalar>& sigma_angstrom,
    BornDiameter born_diameter_angstrom,
    SolvationFactor solvation_factor,
    IonicPermittivity ionic_region_relative_permittivity
) {
    bool has_ion = false;
    for (int charge : parameters.charge_number) {
        has_ion = has_ion || charge != 0;
    }
    if (!has_ion) {
        return Scalar(0.0);
    }

    Scalar fraction_sum = Scalar(0.0);
    Scalar weighted_solvation = Scalar(0.0);
    for (std::size_t component = 0; component < mole_fractions.size(); ++component) {
        fraction_sum += mole_fractions[component];
        weighted_solvation += mole_fractions[component]
            * solvation_factor(component);
    }
    if constexpr (std::is_same_v<Scalar, double>) {
        if (!(fraction_sum > 0.0) || !std::isfinite(weighted_solvation)) {
            throw std::invalid_argument("invalid solvation-factor mixture");
        }
    }
    // EqID: mixture_solvation_factor
    const Scalar mixture_solvation = weighted_solvation / fraction_sum;
    Scalar radial_work = Scalar(0.0);
    for (std::size_t component = 0; component < mole_fractions.size(); ++component) {
        const int charge = parameters.charge_number[component];
        if (charge == 0) {
            continue;
        }
        const Scalar born_diameter = (
            born_diameter_angstrom(component)
            + parameters.born_diameter_sigma_scale[component]
                * sigma_angstrom[component]
        ) * 1.0e-10;
        // EqID: born_shell_diameter
        const Scalar shell_diameter = born_diameter
            * (1.0 + (mixture_solvation - 1.0) / std::abs(charge));
        if constexpr (std::is_same_v<Scalar, double>) {
            if (!(shell_diameter > 0.0) || !std::isfinite(shell_diameter)) {
                throw std::invalid_argument("Born shell diameter is nonpositive");
            }
        }
        // EqID: ssm_ds_radial_work
        radial_work += mole_fractions[component] * charge * charge * (
            (1.0 - 1.0 / bulk_relative_permittivity) / shell_diameter
            + (1.0 - 1.0 / ionic_region_relative_permittivity)
                * (1.0 / born_diameter - 1.0 / shell_diameter)
        );
    }
    // EqID: ssm_ds_born_helmholtz
    return -PHYSICAL_ELEMENTARY_CHARGE * PHYSICAL_ELEMENTARY_CHARGE * radial_work
        / (4.0 * ELECTROLYTE_PI * PHYSICAL_VACUUM_PERMITTIVITY
            * PHYSICAL_BOLTZMANN
           * temperature);
}

template <typename Scalar>
// EqID: ssm_ds_born_helmholtz
Scalar born_ssm_ds_contribution(
    const Scalar& temperature,
    const std::vector<Scalar>& mole_fractions,
    const NativeParameters& parameters,
    const Scalar& bulk_relative_permittivity,
    const std::vector<Scalar>& sigma_angstrom
) {
    return born_ssm_ds_contribution_impl(
        temperature,
        mole_fractions,
        parameters,
        bulk_relative_permittivity,
        sigma_angstrom,
        [&parameters](std::size_t component) {
            return Scalar(parameters.born_diameter_angstrom[component]);
        },
        [&parameters](std::size_t component) {
            return Scalar(parameters.solvation_factor[component]);
        },
        Scalar(parameters.ionic_region_relative_permittivity)
    );
}

template <typename Scalar>
// EqID: ssm_ds_born_helmholtz
Scalar born_contribution(
    const Scalar& temperature,
    const std::vector<Scalar>& mole_fractions,
    const NativeParameters& parameters,
    const Scalar& bulk_relative_permittivity,
    const std::vector<Scalar>& sigma_angstrom
) {
    if (parameters.electrolyte_model
        == ELECTROLYTE_MODEL_FULLY_DISSOCIATED_DEBYE_HUCKEL) {
        return Scalar(0.0);
    }
    return born_ssm_ds_contribution(
        temperature,
        mole_fractions,
        parameters,
        bulk_relative_permittivity,
        sigma_angstrom
    );
}

template <typename Scalar>
// EqID: ssm_ds_born_helmholtz
Scalar born_ssm_ds_contribution_active(
    const Scalar& temperature,
    const std::vector<Scalar>& mole_fractions,
    const NativeParameters& parameters,
    const Scalar& bulk_relative_permittivity,
    const std::vector<Scalar>& sigma_angstrom,
    std::size_t active_component,
    const Scalar& active_born_diameter_angstrom
) {
    return born_ssm_ds_contribution_impl(
        temperature,
        mole_fractions,
        parameters,
        bulk_relative_permittivity,
        sigma_angstrom,
        [&parameters, active_component, &active_born_diameter_angstrom](
            std::size_t component
        ) {
            return component == active_component
                ? active_born_diameter_angstrom
                : Scalar(parameters.born_diameter_angstrom[component]);
        },
        [&parameters](std::size_t component) {
            return Scalar(parameters.solvation_factor[component]);
        },
        Scalar(parameters.ionic_region_relative_permittivity)
    );
}

template <typename Scalar>
// EqID: ssm_ds_born_helmholtz
Scalar born_ssm_ds_contribution_active_solvation_factor(
    const Scalar& temperature,
    const std::vector<Scalar>& mole_fractions,
    const NativeParameters& parameters,
    const Scalar& bulk_relative_permittivity,
    const std::vector<Scalar>& sigma_angstrom,
    std::size_t active_component,
    const Scalar& active_solvation_factor
) {
    return born_ssm_ds_contribution_impl(
        temperature,
        mole_fractions,
        parameters,
        bulk_relative_permittivity,
        sigma_angstrom,
        [&parameters](std::size_t component) {
            return Scalar(parameters.born_diameter_angstrom[component]);
        },
        [&parameters, active_component, &active_solvation_factor](
            std::size_t component
        ) {
            return component == active_component
                ? active_solvation_factor
                : Scalar(parameters.solvation_factor[component]);
        },
        Scalar(parameters.ionic_region_relative_permittivity)
    );
}

template <typename Scalar>
// EqID: ssm_ds_born_helmholtz
Scalar born_ssm_ds_contribution_active_ionic_permittivity(
    const Scalar& temperature,
    const std::vector<Scalar>& mole_fractions,
    const NativeParameters& parameters,
    const Scalar& bulk_relative_permittivity,
    const std::vector<Scalar>& sigma_angstrom,
    const Scalar& active_ionic_region_relative_permittivity
) {
    return born_ssm_ds_contribution_impl(
        temperature,
        mole_fractions,
        parameters,
        bulk_relative_permittivity,
        sigma_angstrom,
        [&parameters](std::size_t component) {
            return Scalar(parameters.born_diameter_angstrom[component]);
        },
        [&parameters](std::size_t component) {
            return Scalar(parameters.solvation_factor[component]);
        },
        active_ionic_region_relative_permittivity
    );
}

}  // namespace epcsaft::detail
