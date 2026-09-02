#pragma once

#include "physics/parameters.hpp"
#include "physics/physical_constants.hpp"

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <type_traits>
#include <vector>

namespace epcsaft {

template <typename Scalar>
Scalar born_contribution(
    const Scalar &temperature, const std::vector<Scalar> &mole_fractions,
    const Parameters &parameters,
    const Scalar &bulk_relative_permittivity,
    const std::vector<Scalar> *born_diameter_override = nullptr,
    const std::vector<Scalar> *solvation_factor_override = nullptr,
    const Scalar *ionic_region_relative_permittivity_override = nullptr) {
  if (parameters.electrolyte_form != ElectrolyteForm::DebyeHuckelBorn) {
    return Scalar(0.0);
  }

  Scalar fraction_sum = Scalar(0.0);
  Scalar weighted_solvation = Scalar(0.0);
  for (std::size_t component = 0; component < mole_fractions.size(); ++component) {
    const Scalar solvation_factor =
        solvation_factor_override == nullptr
            ? Scalar(parameters.solvation_factor[component])
            : (*solvation_factor_override)[component];
    fraction_sum += mole_fractions[component];
    weighted_solvation += mole_fractions[component] * solvation_factor;
  }
  if constexpr (std::is_same_v<Scalar, double>) {
    if (!(fraction_sum > 0.0) || !std::isfinite(weighted_solvation)) {
      throw std::invalid_argument("invalid solvation-factor mixture");
    }
  }
  const Scalar mixture_solvation = weighted_solvation / fraction_sum;
  const Scalar ionic_region_relative_permittivity =
      ionic_region_relative_permittivity_override == nullptr
          ? Scalar(parameters.ionic_region_relative_permittivity)
          : *ionic_region_relative_permittivity_override;
  Scalar radial_work = Scalar(0.0);
  for (std::size_t component = 0; component < mole_fractions.size(); ++component) {
    const int charge = parameters.charge_number[component];
    if (charge == 0) {
      continue;
    }
    const Scalar born_diameter =
        (born_diameter_override == nullptr
             ? Scalar(parameters.born_diameter_angstrom[component])
             : (*born_diameter_override)[component]) *
        ANGSTROM_M;
    // EqID: born_outer_shell_diameter
    const Scalar outer_shell_diameter =
        born_diameter *
        (Scalar(1.0) + parameters.born_c_shell *
                           (mixture_solvation - Scalar(1.0)) /
                           std::abs(charge));
    if constexpr (std::is_same_v<Scalar, double>) {
      if (!(outer_shell_diameter > 0.0) ||
          !std::isfinite(outer_shell_diameter)) {
        throw std::invalid_argument("Born outer-shell diameter is nonpositive");
      }
    }
    // EqID: born_radial_charging_work
    radial_work += mole_fractions[component] * charge * charge *
                   ((Scalar(1.0) - Scalar(1.0) /
                                           bulk_relative_permittivity) /
                        outer_shell_diameter +
                    parameters.born_c_dielectric *
                        (Scalar(1.0) - Scalar(1.0) /
                                           ionic_region_relative_permittivity) *
                        (Scalar(1.0) / born_diameter -
                         Scalar(1.0) / outer_shell_diameter));
  }
  // EqID: born_residual_helmholtz
  return -ELEMENTARY_CHARGE_C * ELEMENTARY_CHARGE_C *
         radial_work /
         (4.0 * PI * VACUUM_PERMITTIVITY_F_PER_M *
          BOLTZMANN_J_PER_K * temperature);
}

}  // namespace epcsaft
