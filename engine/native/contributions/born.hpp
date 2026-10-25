#pragma once

#include "core/parameters.hpp"
#include "core/physical_constants.hpp"

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <type_traits>
#include <vector>

namespace epcsaft {

template <typename Scalar, class ParameterSet, typename Fixed>
Scalar born_radial_charging_work(std::size_t component, const ParameterSet& parameters,
    const Fixed& mixture_solvation, const Scalar& bulk_relative_permittivity) {
    const int charge = parameters.charge_number[component];
    if (charge == 0) return Scalar(0.0);
    const Fixed ionic_region_relative_permittivity = Fixed(parameters.ionic_region_relative_permittivity);
    const Fixed born_diameter = parameters.born_diameter_angstrom[component] * ANGSTROM_M;
    // EqID: born_outer_shell_diameter
    const Fixed outer_shell_diameter = born_diameter
        * (Fixed(1.0) + parameters.born_c_shell * (mixture_solvation - Fixed(1.0)) / std::abs(charge));
    if constexpr (std::is_same_v<Scalar, double>) {
        if (!(outer_shell_diameter > 0.0) || !std::isfinite(outer_shell_diameter))
            throw std::invalid_argument("Born outer-shell diameter is nonpositive");
    }
    // EqID: born_radial_charging_work
    return charge * charge
        * ((Scalar(1.0) - Scalar(1.0) / bulk_relative_permittivity) / outer_shell_diameter
           + parameters.born_c_dielectric
               * (Scalar(1.0) - Scalar(1.0) / ionic_region_relative_permittivity)
               * (Scalar(1.0) / born_diameter - Scalar(1.0) / outer_shell_diameter));
}

template <typename Scalar, class ParameterSet, typename Fixed>
Scalar born_contribution(
    const Fixed& temperature, const std::vector<Fixed>& mole_fractions,
    const ParameterSet& parameters, const Scalar& bulk_relative_permittivity
) {
    if (!parameters.electrolyte_form.is(ElectrolyteForm::DebyeHuckelBorn)) return Scalar(0.0);

    Fixed fraction_sum = Fixed(0.0);
    Fixed weighted_solvation = Fixed(0.0);
    for (std::size_t component = 0; component < mole_fractions.size(); ++component) {
        const Fixed solvation_factor = Fixed(parameters.solvation_factor[component]);
        fraction_sum += mole_fractions[component];
        weighted_solvation += mole_fractions[component] * solvation_factor;
    }
    if constexpr (std::is_same_v<Scalar, double>) {
        if (!(fraction_sum > 0.0) || !std::isfinite(weighted_solvation))
            throw std::invalid_argument("invalid solvation-factor mixture");
    }
    const Fixed mixture_solvation = weighted_solvation / fraction_sum;
    Scalar radial_work = Scalar(0.0);
    for (std::size_t component = 0; component < mole_fractions.size(); ++component) {
        const int charge = parameters.charge_number[component];
        if (charge == 0) continue;
        radial_work += mole_fractions[component] * born_radial_charging_work(
            component, parameters, mixture_solvation, bulk_relative_permittivity);
    }
    // EqID: born_residual_helmholtz
    return -ELEMENTARY_CHARGE_C * ELEMENTARY_CHARGE_C * radial_work
        / (4.0 * PI * VACUUM_PERMITTIVITY_F_PER_M * BOLTZMANN_J_PER_K * temperature);
}

}  // namespace epcsaft
