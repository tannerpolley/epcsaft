#pragma once

#include "electrolyte.hpp"
#include "model.hpp"

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace epcsaft::detail {

inline double born_ssm_ds_contribution(
    double temperature,
    const std::vector<double>& mole_fractions,
    const NativeParameters& parameters,
    double bulk_relative_permittivity
) {
    bool has_ion = false;
    for (int charge : parameters.charge_number) {
        has_ion = has_ion || charge != 0;
    }
    if (!has_ion) {
        return 0.0;
    }

    double fraction_sum = 0.0;
    double weighted_solvation = 0.0;
    for (std::size_t component = 0; component < mole_fractions.size(); ++component) {
        fraction_sum += mole_fractions[component];
        weighted_solvation +=
            mole_fractions[component] * parameters.solvation_factor[component];
    }
    if (!(fraction_sum > 0.0) || !std::isfinite(weighted_solvation)) {
        throw std::invalid_argument("invalid solvation-factor mixture");
    }
    const double mixture_solvation = weighted_solvation / fraction_sum;
    double radial_work = 0.0;
    for (std::size_t component = 0; component < mole_fractions.size(); ++component) {
        const int charge = parameters.charge_number[component];
        if (charge == 0 || mole_fractions[component] == 0.0) {
            continue;
        }
        const double born_diameter =
            parameters.born_diameter_angstrom[component] * 1.0e-10;
        const double shell_diameter = born_diameter
            * (1.0 + (mixture_solvation - 1.0) / std::abs(charge));
        if (!(shell_diameter > 0.0) || !std::isfinite(shell_diameter)) {
            throw std::invalid_argument("Born shell diameter is nonpositive");
        }
        radial_work += mole_fractions[component] * charge * charge * (
            (1.0 - 1.0 / bulk_relative_permittivity) / shell_diameter
            + (1.0 - 1.0 / parameters.ionic_region_relative_permittivity)
                * (1.0 / born_diameter - 1.0 / shell_diameter)
        );
    }
    return -ELEMENTARY_CHARGE * ELEMENTARY_CHARGE * radial_work
        / (4.0 * ELECTROLYTE_PI * VACUUM_PERMITTIVITY * BOLTZMANN
           * temperature);
}

}  // namespace epcsaft::detail
