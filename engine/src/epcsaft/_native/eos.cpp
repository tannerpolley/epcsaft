#include "eos.hpp"

#include "association.hpp"
#include "born.hpp"
#include "density.hpp"
#include "derivatives.hpp"
#include "electrolyte.hpp"
#include "electrolyte_reference.hpp"
#include "helmholtz.hpp"
#include "phase_block.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace epcsaft {
namespace {

// EqExclude: input and domain validation
bool positive_finite(const std::vector<double> &values) {
  return std::all_of(values.begin(), values.end(), [](double value) {
    return std::isfinite(value) && value > 0.0;
  });
}

// EqExclude: input and domain validation
// Explicit zero dispersion is a valid source parameter for ions.
bool nonnegative_finite(const std::vector<double> &values) {
  return std::all_of(values.begin(), values.end(), [](double value) {
    return std::isfinite(value) && value >= 0.0;
  });
}

// EqExclude: input and domain validation
bool finite(const std::vector<double> &values) {
  return std::all_of(values.begin(), values.end(),
                     [](double value) { return std::isfinite(value); });
}

std::vector<std::string> split_identity(const std::string &identity) {
  std::vector<std::string> parts;
  std::size_t begin = 0;
  while (begin <= identity.size()) {
    const std::size_t end = identity.find('/', begin);
    parts.push_back(identity.substr(begin, end - begin));
    if (end == std::string::npos) {
      break;
    }
    begin = end + 1;
  }
  return parts;
}

std::size_t named_index(const std::vector<std::string> &names,
                        const std::string &name, const char *kind) {
  const auto found = std::find(names.begin(), names.end(), name);
  if (found == names.end()) {
    throw std::invalid_argument(std::string("unknown ") + kind +
                                " in active parameter identity: " + name);
  }
  return static_cast<std::size_t>(found - names.begin());
}

// EqExclude: input and domain validation
bool polar_composition_admitted(const NativeParameters &parameters,
                                const std::vector<double> &mole_fractions) {
  if (parameters.polar_formulation == POLAR_FORMULATION_NONE) {
    return true;
  }
  for (std::size_t component = 0; component < mole_fractions.size();
       ++component) {
    if (mole_fractions[component] < 0.0) {
      return false;
    }
    const bool carries_moment =
        parameters.dipole_moment_debye[component] != 0.0 ||
        parameters.quadrupole_moment_debye_angstrom[component] != 0.0;
    if (carries_moment && mole_fractions[component] <= 0.0) {
      return false;
    }
  }
  return true;
}

double log_fugacity_density(const NativeResult &state,
                            const detail::ThermodynamicDerivatives &derivatives,
                            const std::vector<double> &composition,
                            std::size_t component, std::size_t coordinates) {
  // EqID: log_fugacity_density_derivative
  auto hessian = [&derivatives, coordinates](std::size_t row,
                                             std::size_t column) {
    return derivatives.hessian[row * coordinates + column];
  };
  double weighted_composition_density = 0.0;
  for (std::size_t mixture_component = 0;
       mixture_component < composition.size(); ++mixture_component) {
    weighted_composition_density +=
        composition[mixture_component] * hessian(2 + mixture_component, 1);
  }
  const double compressibility_density =
      derivatives.density + state.molar_density * derivatives.density_density;
  const double chemical_potential_density =
      2.0 * derivatives.density +
      state.molar_density * derivatives.density_density +
      hessian(2 + component, 1) - weighted_composition_density;
  return chemical_potential_density -
         compressibility_density / state.compressibility_factor;
}

double fixed_pressure_log_fugacity(const NativeResult &state, double pressure,
                                   std::size_t component,
                                   double fugacity_density) {
  // EqID: fixed_pressure_fugacity_correction
  const double density_correction =
      -(state.pressure - pressure) / state.pressure_density_derivative;
  return (*state.log_fugacity_coefficient)[component] +
         fugacity_density * density_correction;
}

enum class ActiveSolvationParameter {
  BornDiameter,
  SolvationFactor,
  IonicRegionRelativePermittivity,
  SolventRelativePermittivity,
  ComponentParameter,
};

struct ResolvedSolvationParameter final {
  ActiveSolvationParameter family;
  std::size_t component;
  std::optional<NativeActiveParameter> component_parameter;
};

// EqID: fixed_pressure_total_derivative
detail::ReferenceObservation active_solvation_parameter_observation(
    const NativeModel &model, const NativeParameters &parameters,
    double temperature, double pressure, const std::vector<double> &composition,
    const std::array<std::size_t, 2> &ion_indices, const std::string &phase,
    const ResolvedSolvationParameter &active_parameter) {
  const auto state = model.evaluate_pressure(temperature, pressure, composition,
                                             phase);
  if (!state.log_fugacity_coefficient.has_value()) {
    throw std::runtime_error("reference fugacity is unavailable");
  }
  const auto sigma = detail::segment_sigma(temperature, parameters);
  const auto diameter =
      detail::segment_diameter(temperature, parameters, sigma);
  const auto site_fractions =
      detail::solve_site_fractions(temperature, state.molar_density,
                                   composition, parameters, sigma, diameter);
  const auto derivatives = [&] {
    if (active_parameter.family == ActiveSolvationParameter::BornDiameter) {
      return detail::thermodynamic_derivatives_active_born(
          temperature, state.molar_density, composition, parameters,
          site_fractions, active_parameter.component,
          parameters.born_diameter_angstrom[active_parameter.component]);
    }
    if (active_parameter.family == ActiveSolvationParameter::SolvationFactor) {
      return detail::thermodynamic_derivatives_active_solvation_factor(
          temperature, state.molar_density, composition, parameters,
          site_fractions, active_parameter.component,
          parameters.solvation_factor[active_parameter.component]);
    }
    if (active_parameter.family ==
        ActiveSolvationParameter::IonicRegionRelativePermittivity) {
      return detail::thermodynamic_derivatives_active_ionic_region_permittivity(
          temperature, state.molar_density, composition, parameters,
          site_fractions, parameters.ionic_region_relative_permittivity);
    }
    if (active_parameter.family ==
        ActiveSolvationParameter::SolventRelativePermittivity) {
      return detail::thermodynamic_derivatives_active_relative_permittivity(
          temperature, state.molar_density, composition, parameters,
          site_fractions, active_parameter.component,
          parameters.relative_permittivity[active_parameter.component]);
    }
    return detail::thermodynamic_derivatives_active_components(
        temperature, state.molar_density, composition, parameters,
        site_fractions, {*active_parameter.component_parameter});
  }();
  const std::optional<double> active_derivative =
      active_parameter.family == ActiveSolvationParameter::BornDiameter
          ? derivatives.active_born
      : active_parameter.family == ActiveSolvationParameter::SolvationFactor
          ? derivatives.active_solvation_factor
      : active_parameter.family ==
              ActiveSolvationParameter::IonicRegionRelativePermittivity
          ? derivatives.active_ionic_permittivity
      : active_parameter.family ==
              ActiveSolvationParameter::SolventRelativePermittivity
          ? derivatives.active_relative_permittivity
      : derivatives.active_component_parameters.size() == 1
          ? std::optional<double>(derivatives.active_component_parameters[0])
          : std::nullopt;
  if (!active_derivative.has_value()) {
    throw std::runtime_error("active solvation derivative is unavailable");
  }
  const std::size_t components = composition.size();
  const std::size_t coordinates = 3 + components;
  const std::size_t active_coordinate = 2 + components;
  const auto hessian = [&derivatives, coordinates](std::size_t row,
                                                   std::size_t column) {
    return derivatives.hessian[row * coordinates + column];
  };
  const double compressibility = state.compressibility_factor;
  const double pressure_diameter = state.molar_density * state.molar_density *
                                   detail::GAS_CONSTANT * temperature *
                                   hessian(1, active_coordinate);
  const double density_diameter =
      -pressure_diameter / state.pressure_density_derivative;
  std::vector<double> total_derivative(2);
  for (std::size_t ion = 0; ion < 2; ++ion) {
    const std::size_t component = ion_indices[ion];
    double weighted_composition_diameter = 0.0;
    double weighted_composition_density = 0.0;
    for (std::size_t mixture_component = 0; mixture_component < components;
         ++mixture_component) {
      weighted_composition_diameter +=
          composition[mixture_component] *
          hessian(2 + mixture_component, active_coordinate);
      weighted_composition_density +=
          composition[mixture_component] * hessian(2 + mixture_component, 1);
    }
    const double chemical_potential_diameter =
        *active_derivative +
        state.molar_density * hessian(1, active_coordinate) +
        hessian(2 + component, active_coordinate) -
        weighted_composition_diameter;
    const double log_fugacity_diameter =
        chemical_potential_diameter -
        state.molar_density * hessian(1, active_coordinate) / compressibility;
    const double compressibility_density =
        derivatives.density + state.molar_density * derivatives.density_density;
    const double chemical_potential_density =
        2.0 * derivatives.density +
        state.molar_density * derivatives.density_density +
        hessian(2 + component, 1) - weighted_composition_density;
    const double log_fugacity_density =
        chemical_potential_density - compressibility_density / compressibility;
    total_derivative[ion] =
        detail::GAS_CONSTANT * temperature *
        (log_fugacity_diameter + log_fugacity_density * density_diameter);
  }
  return {
      {
          (*state.log_fugacity_coefficient)[ion_indices[0]],
          (*state.log_fugacity_coefficient)[ion_indices[1]],
      },
      total_derivative,
  };
}

// EqID: infinite_dilution_solvation_energy
// EqID: fixed_pressure_total_derivative
NativeIonSolvationDerivativeBlock active_solvation_parameter_reference(
    NativeParameters active_parameters, double temperature, double pressure,
    const std::vector<double> &solvent_fractions,
    const std::array<std::size_t, 2> &ion_indices,
    const std::array<int, 2> &stoichiometry, const std::string &phase,
    const ResolvedSolvationParameter &active_parameter) {
  const NativeModel active_model(active_parameters);
  const auto result = detail::electrolyte_reference_sequence(
      solvent_fractions, active_parameters.molar_mass_kg_per_mol, ion_indices,
      stoichiometry,
      [&active_model, &active_parameters, temperature, pressure,
       &ion_indices, &phase,
       &active_parameter](const std::vector<double> &composition) {
        return active_solvation_parameter_observation(
            active_model, active_parameters, temperature, pressure, composition,
            ion_indices, phase, active_parameter);
      });
  if (!result.derivative.has_value() || result.derivative->size() != 2) {
    throw std::runtime_error("active solvation derivative is unavailable");
  }
  return {
      {
          detail::GAS_CONSTANT * temperature * result.log_fugacity[0],
          detail::GAS_CONSTANT * temperature * result.log_fugacity[1],
      },
      *result.derivative,
      result.molality,
      result.convergence_error,
  };
}

// EqID: fixed_pressure_total_derivative
detail::ReferenceObservation active_aqueous_kij_observation(
    const NativeModel &model, const NativeParameters &parameters,
    double temperature, double pressure, const std::vector<double> &composition,
    const std::array<double, 3> &k_ij) {
  const auto state =
      model.evaluate_pressure(temperature, pressure, composition, "liquid");
  if (!state.log_fugacity_coefficient.has_value()) {
    throw std::runtime_error("reference fugacity is unavailable");
  }
  const auto sigma = detail::segment_sigma(temperature, parameters);
  const auto diameter =
      detail::segment_diameter(temperature, parameters, sigma);
  const auto site_fractions =
      detail::solve_site_fractions(temperature, state.molar_density,
                                   composition, parameters, sigma, diameter);
  const auto derivatives = detail::thermodynamic_derivatives_active_aqueous_kij(
      temperature, state.molar_density, composition, parameters, site_fractions,
      k_ij);
  if (derivatives.active_aqueous_k_ij.size() != 3) {
    throw std::runtime_error("active aqueous k_ij derivative is unavailable");
  }
  const std::size_t components = composition.size();
  const std::size_t coordinates = 2 + components + 3;
  const std::size_t active_offset = 2 + components;
  const auto hessian = [&derivatives, coordinates](std::size_t row,
                                                   std::size_t column) {
    return derivatives.hessian[row * coordinates + column];
  };
  const double compressibility = state.compressibility_factor;
  std::array<double, 2> fugacity_density{};
  std::array<double, 2> corrected_log_fugacity{};
  for (std::size_t ion = 0; ion < 2; ++ion) {
    fugacity_density[ion] = log_fugacity_density(
        state, derivatives, composition, 1 + ion, coordinates);
    corrected_log_fugacity[ion] = fixed_pressure_log_fugacity(
        state, pressure, 1 + ion, fugacity_density[ion]);
  }
  std::vector<double> total_derivative(6);
  for (std::size_t active = 0; active < 3; ++active) {
    const std::size_t active_coordinate = active_offset + active;
    const double pressure_k = state.molar_density * state.molar_density *
                              detail::GAS_CONSTANT * temperature *
                              hessian(1, active_coordinate);
    const double density_k = -pressure_k / state.pressure_density_derivative;
    for (std::size_t ion = 0; ion < 2; ++ion) {
      const std::size_t component = 1 + ion;
      double weighted_composition_k = 0.0;
      for (std::size_t mixture_component = 0; mixture_component < components;
           ++mixture_component) {
        weighted_composition_k +=
            composition[mixture_component] *
            hessian(2 + mixture_component, active_coordinate);
      }
      const double chemical_potential_k =
          derivatives.active_aqueous_k_ij[active] +
          state.molar_density * hessian(1, active_coordinate) +
          hessian(2 + component, active_coordinate) - weighted_composition_k;
      const double log_fugacity_k =
          chemical_potential_k -
          state.molar_density * hessian(1, active_coordinate) / compressibility;
      total_derivative[ion * 3 + active] =
          log_fugacity_k + fugacity_density[ion] * density_k;
    }
  }
  return {
      corrected_log_fugacity,
      std::move(total_derivative),
  };
}

// EqID: fixed_pressure_total_derivative
detail::ReferenceObservation active_solvation_factor_observation(
    const NativeModel &model, const NativeParameters &parameters,
    double temperature, double pressure, const std::vector<double> &composition,
    double solvation_factor) {
  const auto state =
      model.evaluate_pressure(temperature, pressure, composition, "liquid");
  if (!state.log_fugacity_coefficient.has_value()) {
    throw std::runtime_error("reference fugacity is unavailable");
  }
  const auto sigma = detail::segment_sigma(temperature, parameters);
  const auto diameter =
      detail::segment_diameter(temperature, parameters, sigma);
  const auto site_fractions =
      detail::solve_site_fractions(temperature, state.molar_density,
                                   composition, parameters, sigma, diameter);
  const auto derivatives =
      detail::thermodynamic_derivatives_active_solvation_factor(
          temperature, state.molar_density, composition, parameters,
          site_fractions, 0, solvation_factor);
  if (!derivatives.active_solvation_factor.has_value()) {
    throw std::runtime_error(
        "active solvation-factor derivative is unavailable");
  }
  const std::size_t components = composition.size();
  const std::size_t coordinates = 2 + components + 1;
  const std::size_t active_coordinate = 2 + components;
  const auto hessian = [&derivatives, coordinates](std::size_t row,
                                                   std::size_t column) {
    return derivatives.hessian[row * coordinates + column];
  };
  const double pressure_factor = state.molar_density * state.molar_density *
                                 detail::GAS_CONSTANT * temperature *
                                 hessian(1, active_coordinate);
  const double density_factor =
      -pressure_factor / state.pressure_density_derivative;
  std::array<double, 2> corrected_log_fugacity{};
  std::vector<double> total_derivative(2);
  for (std::size_t ion = 0; ion < 2; ++ion) {
    const std::size_t component = 1 + ion;
    const double fugacity_density = log_fugacity_density(
        state, derivatives, composition, component, coordinates);
    corrected_log_fugacity[ion] = fixed_pressure_log_fugacity(
        state, pressure, component, fugacity_density);
    double weighted_composition_factor = 0.0;
    for (std::size_t mixture_component = 0; mixture_component < components;
         ++mixture_component) {
      weighted_composition_factor +=
          composition[mixture_component] *
          hessian(2 + mixture_component, active_coordinate);
    }
    const double chemical_potential_factor =
        *derivatives.active_solvation_factor +
        state.molar_density * hessian(1, active_coordinate) +
        hessian(2 + component, active_coordinate) - weighted_composition_factor;
    const double log_fugacity_factor =
        chemical_potential_factor - state.molar_density *
                                        hessian(1, active_coordinate) /
                                        state.compressibility_factor;
    total_derivative[ion] =
        log_fugacity_factor + fugacity_density * density_factor;
  }
  return {corrected_log_fugacity, std::move(total_derivative)};
}

// EqExclude: transport and buffer mechanics
std::optional<double> single_solvent_log_fugacity(
    const NativeModel &model, double temperature, double pressure,
    const std::vector<double> &solvent_fractions, const std::string &phase) {
  std::optional<std::size_t> solvent_index;
  for (std::size_t component = 0; component < solvent_fractions.size();
       ++component) {
    if (solvent_fractions[component] <= 0.0) {
      continue;
    }
    if (solvent_index.has_value()) {
      return std::nullopt;
    }
    solvent_index = component;
  }
  if (!solvent_index.has_value()) {
    return std::nullopt;
  }
  const auto state =
      model.evaluate_pressure(temperature, pressure, solvent_fractions, phase);
  if (!state.log_fugacity_coefficient.has_value()) {
    throw std::runtime_error("pure-solvent reference fugacity is unavailable");
  }
  return (*state.log_fugacity_coefficient)[*solvent_index];
}

} // namespace

// EqExclude: input and domain validation
NativeModel::NativeModel(NativeParameters parameters)
    : parameters_(std::move(parameters)) {
  const std::size_t components = parameters_.segment_count.size();
  const std::size_t terms = parameters_.sigma_term_component.size();
  const std::size_t sites = parameters_.site_component.size();
  if (components == 0 || parameters_.sigma_constant.size() != components ||
      parameters_.dispersion_energy_over_k.size() != components ||
      parameters_.k_ij.size() != components * components ||
      parameters_.l_ij.size() != components * components ||
      parameters_.sigma_term_amplitude.size() != terms ||
      parameters_.sigma_term_exponent.size() != terms ||
      parameters_.site_multiplicity.size() != sites ||
      parameters_.association_active.size() != sites * sites ||
      parameters_.association_energy_over_k.size() != sites * sites ||
      parameters_.association_volume.size() != sites * sites ||
      parameters_.molar_mass_kg_per_mol.size() != components ||
      parameters_.charge_number.size() != components ||
      parameters_.relative_permittivity.size() != components ||
      parameters_.born_diameter_angstrom.size() != components ||
      parameters_.born_diameter_sigma_scale.size() != components ||
      parameters_.solvation_factor.size() != components ||
      parameters_.dielectric_volume_m3_per_mol.size() != components ||
      parameters_.dielectric_temperature_k.size() != components ||
      parameters_.ion_specific_suppression_coefficient.size() != components ||
      parameters_.dipole_scaling_debye2.size() != components ||
      parameters_.polarizability_scaling_angstrom3.size() != components ||
      parameters_.correlation_integral_parameter.size() != components ||
      parameters_.dipole_moment_debye.size() != components ||
      parameters_.quadrupole_moment_debye_angstrom.size() != components) {
    throw std::invalid_argument("inconsistent native model dimensions");
  }
  if (!positive_finite(parameters_.segment_count) ||
      !positive_finite(parameters_.sigma_constant) ||
      !nonnegative_finite(parameters_.dispersion_energy_over_k) ||
      !finite(parameters_.k_ij) || !finite(parameters_.l_ij) ||
      !finite(parameters_.sigma_term_amplitude) ||
      !finite(parameters_.sigma_term_exponent) ||
      !finite(parameters_.association_energy_over_k) ||
      !finite(parameters_.association_volume) ||
      !finite(parameters_.molar_mass_kg_per_mol) ||
      !finite(parameters_.relative_permittivity) ||
      !std::isfinite(parameters_.ion_fraction_suppression_coefficient) ||
      !finite(parameters_.born_diameter_angstrom) ||
      !finite(parameters_.born_diameter_sigma_scale) ||
      !finite(parameters_.solvation_factor) ||
      !finite(parameters_.dielectric_volume_m3_per_mol) ||
      !finite(parameters_.dielectric_temperature_k) ||
      !finite(parameters_.ion_specific_suppression_coefficient) ||
      !finite(parameters_.dipole_scaling_debye2) ||
      !finite(parameters_.polarizability_scaling_angstrom3) ||
      !finite(parameters_.correlation_integral_parameter) ||
      !finite(parameters_.dipole_moment_debye) ||
      !finite(parameters_.quadrupole_moment_debye_angstrom) ||
      !std::isfinite(parameters_.ionic_region_relative_permittivity)) {
    throw std::invalid_argument("invalid native model parameters");
  }
  for (int component : parameters_.sigma_term_component) {
    if (component < 0 || static_cast<std::size_t>(component) >= components) {
      throw std::invalid_argument("invalid sigma correlation component");
    }
  }
  for (std::size_t site = 0; site < sites; ++site) {
    const int component = parameters_.site_component[site];
    if (component < 0 || static_cast<std::size_t>(component) >= components ||
        parameters_.site_multiplicity[site] <= 0) {
      throw std::invalid_argument("invalid association site");
    }
  }
  for (std::size_t pair = 0; pair < sites * sites; ++pair) {
    if (parameters_.association_active[pair] != 0 &&
        (parameters_.association_energy_over_k[pair] <= 0.0 ||
         parameters_.association_volume[pair] <= 0.0)) {
      throw std::invalid_argument("invalid association pair");
    }
  }
  const bool ionic = std::any_of(parameters_.charge_number.begin(),
                                 parameters_.charge_number.end(),
                                 [](int charge) { return charge != 0; });
  if (parameters_.polar_formulation != POLAR_FORMULATION_NONE &&
      parameters_.polar_formulation !=
          POLAR_FORMULATION_GROSS_VRABEC_POINT_MULTIPOLE) {
    throw std::invalid_argument("invalid polar formulation");
  }
  const bool has_multipole =
      std::any_of(parameters_.dipole_moment_debye.begin(),
                  parameters_.dipole_moment_debye.end(),
                  [](double value) { return value != 0.0; }) ||
      std::any_of(parameters_.quadrupole_moment_debye_angstrom.begin(),
                  parameters_.quadrupole_moment_debye_angstrom.end(),
                  [](double value) { return value != 0.0; });
  if (parameters_.polar_formulation == POLAR_FORMULATION_NONE &&
      has_multipole) {
    throw std::invalid_argument(
        "physical multipoles require a polar formulation");
  }
  for (std::size_t component = 0; component < components; ++component) {
    if (parameters_.dipole_moment_debye[component] < 0.0 ||
        parameters_.quadrupole_moment_debye_angstrom[component] < 0.0 ||
        (parameters_.charge_number[component] != 0 &&
         (parameters_.dipole_moment_debye[component] != 0.0 ||
          parameters_.quadrupole_moment_debye_angstrom[component] != 0.0))) {
      throw std::invalid_argument("invalid physical multipole parameters");
    }
  }
  if (ionic) {
    if (parameters_.electrolyte_model != ELECTROLYTE_MODEL_SSM_DS &&
        parameters_.electrolyte_model != ELECTROLYTE_MODEL_ORIGINAL_BORN &&
        parameters_.electrolyte_model !=
            ELECTROLYTE_MODEL_FULLY_DISSOCIATED_DEBYE_HUCKEL) {
      throw std::invalid_argument("invalid electrolyte model selection");
    }
    const bool ssm_ds =
        parameters_.electrolyte_model == ELECTROLYTE_MODEL_SSM_DS;
    const bool original_born =
        parameters_.electrolyte_model == ELECTROLYTE_MODEL_ORIGINAL_BORN;
    const bool fully_dissociated =
        parameters_.electrolyte_model ==
        ELECTROLYTE_MODEL_FULLY_DISSOCIATED_DEBYE_HUCKEL;
    const bool valid_ssm_permittivity =
        parameters_.relative_permittivity_formulation ==
            RELATIVE_PERMITTIVITY_FORMULATION_ION_FRACTION_SUPPRESSION ||
        parameters_.relative_permittivity_formulation ==
            RELATIVE_PERMITTIVITY_FORMULATION_DENSITY_VOLUME_TEMPERATURE ||
        parameters_.relative_permittivity_formulation ==
            RELATIVE_PERMITTIVITY_FORMULATION_ION_SPECIFIC_SUPPRESSION ||
        parameters_.relative_permittivity_formulation ==
            RELATIVE_PERMITTIVITY_FORMULATION_DIPOLAR_DENSITY_PERTURBATION;
    if ((ssm_ds && !valid_ssm_permittivity) ||
        (original_born &&
         parameters_.relative_permittivity_formulation !=
             RELATIVE_PERMITTIVITY_FORMULATION_COMPONENT_MIXING) ||
        (fully_dissociated &&
         parameters_.relative_permittivity_formulation !=
             RELATIVE_PERMITTIVITY_FORMULATION_SOLVENT_ONLY)) {
      throw std::invalid_argument(
          "incompatible electrolyte and permittivity models");
    }
    std::size_t solvent_count = 0;
    for (std::size_t component = 0; component < components; ++component) {
      const bool solvent = parameters_.charge_number[component] == 0;
      solvent_count += static_cast<std::size_t>(solvent);
      const bool needs_tabulated_permittivity =
          parameters_.relative_permittivity_formulation ==
              RELATIVE_PERMITTIVITY_FORMULATION_ION_FRACTION_SUPPRESSION ||
          parameters_.relative_permittivity_formulation ==
              RELATIVE_PERMITTIVITY_FORMULATION_COMPONENT_MIXING ||
          parameters_.relative_permittivity_formulation ==
              RELATIVE_PERMITTIVITY_FORMULATION_ION_SPECIFIC_SUPPRESSION ||
          parameters_.relative_permittivity_formulation ==
              RELATIVE_PERMITTIVITY_FORMULATION_SOLVENT_ONLY;
      if (solvent &&
          (!(parameters_.molar_mass_kg_per_mol[component] > 0.0) ||
           (needs_tabulated_permittivity &&
            !(parameters_.relative_permittivity[component] > 0.0)))) {
        throw std::invalid_argument("invalid solvent dielectric parameters");
      }
      if (ssm_ds) {
        if (!(parameters_.ionic_region_relative_permittivity > 0.0) ||
            !(parameters_.solvation_factor[component] > 0.0)) {
          throw std::invalid_argument("invalid SSM+DS electrolyte parameters");
        }
        if (!solvent &&
            !(parameters_.born_diameter_angstrom[component] > 0.0)) {
          throw std::invalid_argument("invalid fitted Born diameter");
        }
        if (parameters_.born_diameter_sigma_scale[component] != 0.0) {
          throw std::invalid_argument("invalid fitted Born diameter scale");
        }
      } else if (original_born) {
        const double expected_scale = solvent ? 0.0 : 1.0;
        if (parameters_.born_diameter_angstrom[component] != 0.0 ||
            parameters_.born_diameter_sigma_scale[component] !=
                expected_scale ||
            parameters_.solvation_factor[component] != 1.0 ||
            parameters_.ionic_region_relative_permittivity != 1.0 ||
            (!solvent && !(parameters_.relative_permittivity[component] > 0.0))) {
          throw std::invalid_argument("invalid original-Born limiting parameters");
        }
      } else if (parameters_.electrolyte_model ==
                     ELECTROLYTE_MODEL_FULLY_DISSOCIATED_DEBYE_HUCKEL &&
                 ((parameters_.charge_number[component] != 0 &&
                  parameters_.relative_permittivity[component] != 0.0) ||
                  parameters_.born_diameter_angstrom[component] != 0.0 ||
                  parameters_.born_diameter_sigma_scale[component] != 0.0 ||
                  parameters_.solvation_factor[component] != 0.0)) {
        throw std::invalid_argument(
            "invalid fully dissociated electrolyte parameters");
      }
      if (parameters_.relative_permittivity_formulation ==
              RELATIVE_PERMITTIVITY_FORMULATION_DENSITY_VOLUME_TEMPERATURE &&
          solvent &&
          (!(parameters_.dielectric_volume_m3_per_mol[component] > 0.0) ||
           !(parameters_.dielectric_temperature_k[component] > 0.0))) {
        throw std::invalid_argument(
            "invalid density formulation solvent parameters");
      }
      if (parameters_.relative_permittivity_formulation ==
              RELATIVE_PERMITTIVITY_FORMULATION_ION_SPECIFIC_SUPPRESSION &&
          !solvent &&
          !(parameters_.ion_specific_suppression_coefficient[component] >
            0.0)) {
        throw std::invalid_argument(
            "invalid ion-specific suppression parameter");
      }
      if (parameters_.relative_permittivity_formulation ==
          RELATIVE_PERMITTIVITY_FORMULATION_DIPOLAR_DENSITY_PERTURBATION) {
        if (parameters_.dipole_scaling_debye2[component] < 0.0 ||
            parameters_.polarizability_scaling_angstrom3[component] < 0.0 ||
            !(parameters_.correlation_integral_parameter[component] > 0.0)) {
          throw std::invalid_argument(
              "invalid dipolar perturbation parameters");
        }
      }
    }
    if (parameters_.relative_permittivity_formulation ==
            RELATIVE_PERMITTIVITY_FORMULATION_ION_FRACTION_SUPPRESSION &&
        !(parameters_.ion_fraction_suppression_coefficient > 0.0)) {
      throw std::invalid_argument("invalid ion-fraction suppression parameter");
    }
    if (parameters_.relative_permittivity_formulation ==
            RELATIVE_PERMITTIVITY_FORMULATION_ION_SPECIFIC_SUPPRESSION &&
        solvent_count != 1) {
      throw std::invalid_argument("mixed-solvent volume data are required");
    }
    if (fully_dissociated &&
        (parameters_.ion_fraction_suppression_coefficient != 0.0 ||
         parameters_.ionic_region_relative_permittivity != 0.0 ||
         std::any_of(parameters_.l_ij.begin(), parameters_.l_ij.end(),
                     [](double value) { return value != 0.0; }))) {
      throw std::invalid_argument(
          "invalid fully dissociated electrolyte parameters");
    }
  } else if (parameters_.electrolyte_model != ELECTROLYTE_MODEL_NONE ||
             parameters_.relative_permittivity_formulation !=
                 RELATIVE_PERMITTIVITY_FORMULATION_NONE) {
    throw std::invalid_argument("neutral model has electrolyte selection");
  }
}

NativeModel NativeModel::with_parameter_identities(
    const std::vector<std::string> &component_ids,
    const std::vector<std::string> &site_ids) const {
  if (component_ids.size() != parameters_.segment_count.size() ||
      site_ids.size() != parameters_.site_component.size() ||
      std::any_of(component_ids.begin(), component_ids.end(),
                  [](const std::string &value) { return value.empty(); }) ||
      std::any_of(site_ids.begin(), site_ids.end(),
                  [](const std::string &value) { return value.empty(); })) {
    throw std::invalid_argument("invalid native identity metadata");
  }
  for (std::size_t site = 0; site < site_ids.size(); ++site) {
    const auto parts = split_identity(site_ids[site]);
    const int component = parameters_.site_component[site];
    if (parts.size() != 2 || component < 0 ||
        static_cast<std::size_t>(component) >= component_ids.size() ||
        parts[0] != component_ids[static_cast<std::size_t>(component)] ||
        parts[1].empty()) {
      throw std::invalid_argument(
          "site identity metadata must use component/site identities");
    }
  }
  auto sorted_components = component_ids;
  auto sorted_sites = site_ids;
  std::sort(sorted_components.begin(), sorted_components.end());
  std::sort(sorted_sites.begin(), sorted_sites.end());
  if (std::adjacent_find(sorted_components.begin(), sorted_components.end()) !=
          sorted_components.end() ||
      std::adjacent_find(sorted_sites.begin(), sorted_sites.end()) !=
          sorted_sites.end()) {
    throw std::invalid_argument("duplicate native identity metadata");
  }
  NativeModel result(parameters_);
  result.component_ids_ = component_ids;
  result.site_ids_ = site_ids;
  return result;
}

// EqID: pure_active_parameter_mapping
NativeModel NativeModel::with_active_parameters(
    const std::vector<NativeParameterRecord> &active_parameters) const {
  if (active_parameters.empty()) {
    throw std::invalid_argument("active parameter set must be nonempty");
  }
  NativeParameters updated = parameters_;
  std::vector<std::pair<NativeFixedTopologyParameterKind, std::size_t>> seen;
  seen.reserve(active_parameters.size());
  for (const auto &record : active_parameters) {
    const auto parameter = resolve_fixed_topology_parameter(record);
    std::size_t identity = parameter.component;
    switch (parameter.kind) {
    case NativeFixedTopologyParameterKind::SegmentCount:
      updated.segment_count[parameter.component] = parameter.value;
      break;
    case NativeFixedTopologyParameterKind::SegmentDiameter:
    case NativeFixedTopologyParameterKind::SegmentDiameterCorrelationConstant:
      updated.sigma_constant[parameter.component] = parameter.value;
      break;
    case NativeFixedTopologyParameterKind::SegmentDiameterCorrelationAmplitude:
      identity = parameter.correlation_term;
      updated.sigma_term_amplitude[parameter.correlation_term] =
          parameter.value;
      break;
    case NativeFixedTopologyParameterKind::SegmentDiameterCorrelationExponent:
      identity = parameter.correlation_term;
      updated.sigma_term_exponent[parameter.correlation_term] = parameter.value;
      break;
    case NativeFixedTopologyParameterKind::DispersionEnergyOverK:
      updated.dispersion_energy_over_k[parameter.component] = parameter.value;
      break;
    case NativeFixedTopologyParameterKind::AssociationEnergyOverK:
    case NativeFixedTopologyParameterKind::AssociationVolume: {
      const std::size_t sites = updated.site_component.size();
      const std::size_t forward = parameter.site_a * sites + parameter.site_b;
      const std::size_t reverse = parameter.site_b * sites + parameter.site_a;
      identity = forward;
      auto &values =
          parameter.kind == NativeFixedTopologyParameterKind::AssociationEnergyOverK
              ? updated.association_energy_over_k
              : updated.association_volume;
      values[forward] = parameter.value;
      values[reverse] = parameter.value;
      break;
    }
    case NativeFixedTopologyParameterKind::PairDispersion:
    case NativeFixedTopologyParameterKind::PairSize:
    case NativeFixedTopologyParameterKind::CrossAssociationEnergyCorrection:
      throw std::logic_error(
          "pair coordinates require the resolved mixture-coordinate path");
    }
    const auto key = std::pair{parameter.kind, identity};
    if (std::find(seen.begin(), seen.end(), key) != seen.end()) {
      throw std::invalid_argument("active parameter is duplicated");
    }
    seen.push_back(key);
  }
  NativeModel result(std::move(updated));
  result.component_ids_ = component_ids_;
  result.site_ids_ = site_ids_;
  return result;
}

NativeFixedTopologyParameter NativeModel::resolve_fixed_topology_parameter(
    const NativeParameterRecord &parameter) const {
  if (!std::isfinite(parameter.value)) {
    throw std::invalid_argument("active parameter values must be finite");
  }
  if (component_ids_.size() != parameters_.segment_count.size() ||
      site_ids_.size() != parameters_.site_component.size()) {
    throw std::logic_error(
        "canonical active parameters require native identity metadata");
  }
  const auto parts = split_identity(parameter.identity);
  if (parts.size() == 3 && parts[0] == "component") {
    const std::size_t component =
        named_index(component_ids_, parts[1], "component");
    NativeFixedTopologyParameterKind kind;
    if (parts[2] == "segment_count") {
      kind = NativeFixedTopologyParameterKind::SegmentCount;
    } else if (parts[2] == "segment_diameter") {
      kind = NativeFixedTopologyParameterKind::SegmentDiameter;
    } else if (parts[2] == "dispersion_energy_over_k") {
      kind = NativeFixedTopologyParameterKind::DispersionEnergyOverK;
    } else {
      throw std::logic_error("unsupported active parameter identity: " +
                             parameter.identity);
    }
    if (parameter.value <= 0.0) {
      throw std::invalid_argument(
          "direct component active values must be positive");
    }
    return {kind, component, 0, 0, parameter.value};
  }
  if ((parts.size() == 5 || parts.size() == 6) && parts[0] == "component" &&
      parts[2] == "segment_diameter") {
    const std::size_t component =
        named_index(component_ids_, parts[1], "component");
    if (parts.size() == 5 && parts[4] == "constant") {
      if (parameter.value <= 0.0) {
        throw std::invalid_argument(
            "segment-diameter correlation constant must be positive");
      }
      return {
          NativeFixedTopologyParameterKind::SegmentDiameterCorrelationConstant,
          component,
          0,
          0,
          parameter.value,
      };
    }
    if (parts.size() == 6 && parts[4].compare(0, 5, "term-") == 0 &&
        (parts[5] == "amplitude" || parts[5] == "exponent_coefficient")) {
      const std::string ordinal_text = parts[4].substr(5);
      if (ordinal_text.empty() ||
          !std::all_of(
              ordinal_text.begin(), ordinal_text.end(),
              [](char value) { return value >= '0' && value <= '9'; })) {
        throw std::logic_error(
            "invalid segment-diameter correlation term identity");
      }
      std::size_t ordinal = 0;
      for (const char value : ordinal_text) {
        const std::size_t digit = static_cast<std::size_t>(value - '0');
        if (ordinal > (std::numeric_limits<std::size_t>::max() - digit) / 10) {
          throw std::logic_error(
              "segment-diameter correlation term index is too large");
        }
        ordinal = ordinal * 10 + digit;
      }
      std::size_t observed = 0;
      std::size_t term = parameters_.sigma_term_component.size();
      for (std::size_t index = 0;
           index < parameters_.sigma_term_component.size(); ++index) {
        if (parameters_.sigma_term_component[index] !=
            static_cast<int>(component)) {
          continue;
        }
        if (observed == ordinal) {
          term = index;
          break;
        }
        ++observed;
      }
      if (term == parameters_.sigma_term_component.size()) {
        throw std::logic_error(
            "segment-diameter correlation term is outside the model");
      }
      const auto kind = parts[5] == "amplitude"
                            ? NativeFixedTopologyParameterKind::
                                  SegmentDiameterCorrelationAmplitude
                            : NativeFixedTopologyParameterKind::
                                  SegmentDiameterCorrelationExponent;
      return {kind, component, 0, 0, parameter.value, term};
    }
  }
  if (parts.size() == 6 && parts[0] == "association") {
    const std::size_t site_a = named_index(site_ids_, parts[1] + "/" + parts[2],
                                           "association endpoint");
    const std::size_t site_b = named_index(site_ids_, parts[3] + "/" + parts[4],
                                           "association endpoint");
    NativeFixedTopologyParameterKind kind;
    if (parts[5] == "energy_over_k") {
      kind = NativeFixedTopologyParameterKind::AssociationEnergyOverK;
    } else if (parts[5] == "volume") {
      kind = NativeFixedTopologyParameterKind::AssociationVolume;
    } else {
      throw std::logic_error("unsupported active parameter identity: " +
                             parameter.identity);
    }
    if (parameter.value <= 0.0) {
      throw std::invalid_argument("association active values must be positive");
    }
    return {
        kind,
        0,
        std::min(site_a, site_b),
        std::max(site_a, site_b),
        parameter.value,
    };
  }
  throw std::logic_error("unsupported active parameter identity: " +
                         parameter.identity);
}

NativeActiveParameter NativeModel::resolve_component_parameter(
    const NativeParameterRecord &parameter) const {
  if (!std::isfinite(parameter.value) || parameter.value <= 0.0) {
    throw std::invalid_argument(
        "active parameter values must be positive and finite");
  }
  if (component_ids_.size() != parameters_.segment_count.size()) {
    throw std::logic_error(
        "canonical active parameters require native identity metadata");
  }
  const auto parts = split_identity(parameter.identity);
  if (parts.size() != 3 || parts[0] != "component") {
    throw std::logic_error("unsupported fugacity derivative identity: " +
                           parameter.identity);
  }
  const std::size_t component =
      named_index(component_ids_, parts[1], "component");
  if (parts[2] == "segment_diameter") {
    return {
        ActiveComponentParameter::SegmentDiameter,
        component,
        parameter.value,
    };
  }
  if (parts[2] == "dispersion_energy_over_k") {
    return {
        ActiveComponentParameter::DispersionEnergyOverK,
        component,
        parameter.value,
    };
  }
  throw std::logic_error("unsupported fugacity derivative identity: " +
                         parameter.identity);
}

// EqExclude: transport and buffer mechanics
bool NativeModel::has_temperature_independent_segment_diameter(
    std::size_t component) const {
  if (component >= parameters_.sigma_constant.size()) {
    throw std::out_of_range("segment-diameter component is out of range");
  }
  return std::none_of(parameters_.sigma_term_component.begin(),
                      parameters_.sigma_term_component.end(),
                      [component](int term_component) {
                        return term_component == static_cast<int>(component);
                      });
}

// EqExclude: input and domain validation
bool NativeModel::supports_pure_transport() const {
  return detail::supports_pure_phase_block(parameters_);
}

// EqExclude: input and domain validation
bool NativeModel::supports_fixed_topology_association() const {
  return detail::supports_fixed_topology_associating_phase_block(parameters_);
}

// EqExclude: input and domain validation
bool NativeModel::supports_fixed_topology_associating_mixture()
    const {
  return detail::supports_fixed_topology_associating_mixture_phase_block(
      parameters_);
}

// EqExclude: input and domain validation
bool NativeModel::supports_mixture_transport() const {
  return detail::supports_mixture_phase_block(parameters_);
}

// EqExclude: input and domain validation
bool NativeModel::supports_electrolyte_transport() const {
  return detail::supports_electrolyte_phase_block(parameters_);
}

// EqExclude: serialization
std::size_t NativeModel::component_count() const {
  return parameters_.segment_count.size();
}

// EqExclude: serialization
const std::vector<int> &NativeModel::charge_numbers() const {
  return parameters_.charge_number;
}

// EqExclude: serialization
const std::vector<int> &NativeModel::site_components() const {
  return parameters_.site_component;
}

// EqExclude: serialization
const std::vector<int> &NativeModel::site_multiplicities() const {
  return parameters_.site_multiplicity;
}

// EqExclude: serialization
const std::vector<int> &NativeModel::association_active_pairs() const {
  return parameters_.association_active;
}

// EqExclude: serialization
const std::vector<double> &NativeModel::molar_masses_kg_per_mol() const {
  return parameters_.molar_mass_kg_per_mol;
}

// EqExclude: transport and buffer mechanics
NativePurePhaseBlock NativeModel::evaluate_pure_phase(double temperature,
                                                      double amount,
                                                      double volume) const {
  return detail::pure_phase_block(temperature, amount, volume, parameters_);
}

// EqExclude: transport and buffer mechanics
NativeFixedTopologyPhaseBlock NativeModel::evaluate_fixed_topology_association(
    double temperature, double amount, double volume,
    const std::vector<NativeFixedTopologyParameter> &active_parameters) const {
  return detail::fixed_topology_parameterized_pure_phase_block(
      temperature, amount, volume, active_parameters, parameters_);
}

NativeFixedTopologyPhaseBlock
NativeModel::evaluate_fixed_topology_association_active_parameters(
    double temperature, double amount, double volume,
    const std::vector<NativeParameterRecord> &active_parameters) const {
  std::vector<NativeFixedTopologyParameter> resolved;
  resolved.reserve(active_parameters.size());
  for (const auto &parameter : active_parameters) {
    resolved.push_back(resolve_fixed_topology_parameter(parameter));
  }
  return evaluate_fixed_topology_association(temperature, amount, volume,
                                             resolved);
}

// EqExclude: transport and buffer mechanics
NativeMixturePhaseBlock
NativeModel::evaluate_fixed_topology_associating_mixture(
    double temperature, const std::vector<double> &amounts, double volume,
    const std::vector<NativeFixedTopologyParameter> &active_parameters) const {
  if (!detail::supports_fixed_topology_associating_mixture_phase_block(
          parameters_)) {
    throw std::logic_error(
        "fixed-topology associating-mixture parameters require at least "
        "two components and one active association edge");
  }
  return detail::fixed_topology_parameterized_mixture_phase_block(
      temperature, amounts, volume, active_parameters, parameters_);
}

NativeMixturePhaseBlock
NativeModel::evaluate_fixed_topology_associating_mixture_active_parameters(
    double temperature, const std::vector<double> &amounts, double volume,
    const std::vector<NativeParameterRecord> &active_parameters) const {
  std::vector<NativeFixedTopologyParameter> resolved;
  resolved.reserve(active_parameters.size());
  for (const auto &parameter : active_parameters) {
    resolved.push_back(resolve_fixed_topology_parameter(parameter));
  }
  return evaluate_fixed_topology_associating_mixture(temperature, amounts,
                                                     volume, resolved);
}

// EqExclude: transport and buffer mechanics
NativeMixturePhaseBlock
NativeModel::evaluate_mixture_phase(double temperature,
                                    const std::vector<double> &amounts,
                                    double volume) const {
  return detail::mixture_phase_block(temperature, amounts, volume, parameters_);
}

// EqExclude: transport and buffer mechanics
NativeMixturePhaseBlock
NativeModel::evaluate_mixture_phase_pair_coordinates(
    double temperature, const std::vector<double> &amounts, double volume,
    const std::vector<NativePairCoordinate> &coordinates) const {
  if (component_ids_.size() != parameters_.segment_count.size() ||
      coordinates.empty()) {
    throw std::invalid_argument(
        "resolved pair coordinates require native component identities");
  }
  std::vector<NativeFixedTopologyParameter> resolved;
  resolved.reserve(coordinates.size());
  for (const auto &coordinate : coordinates) {
    if (coordinate.component_a >= component_ids_.size() ||
        coordinate.component_b >= component_ids_.size() ||
        coordinate.component_a == coordinate.component_b ||
        !std::isfinite(coordinate.value)) {
      throw std::invalid_argument("invalid resolved pair-coordinate endpoints");
    }
    const auto endpoint_ids = std::minmax(
        component_ids_[coordinate.component_a],
        component_ids_[coordinate.component_b]);
    const auto component_indices =
        std::minmax(coordinate.component_a, coordinate.component_b);
    NativeFixedTopologyParameterKind kind;
    std::string expected_identity;
    if (coordinate.family == "k_ij" || coordinate.family == "l_ij") {
      if (!coordinate.association_role_a.empty() ||
          !coordinate.association_role_b.empty() ||
          !coordinate.association_edges.empty() ||
          !coordinate.association_energy_scales.empty()) {
        throw std::invalid_argument(
            "non-association pair coordinate contains association metadata");
      }
      kind = coordinate.family == "k_ij"
                 ? NativeFixedTopologyParameterKind::PairDispersion
                 : NativeFixedTopologyParameterKind::PairSize;
      expected_identity = "pair/" + endpoint_ids.first + "/" +
                          endpoint_ids.second + "/" + coordinate.family;
    } else if (coordinate.family == "k_hb_ij") {
      if (coordinate.association_role_a.empty() ||
          coordinate.association_role_a > coordinate.association_role_b) {
        throw std::invalid_argument(
            "cross-association coordinate role pair is not canonical");
      }
      kind = NativeFixedTopologyParameterKind::
          CrossAssociationEnergyCorrection;
      expected_identity = "pair/" + endpoint_ids.first + "/" +
                          endpoint_ids.second + "/association/" +
                          coordinate.association_role_a + "/" +
                          coordinate.association_role_b + "/k_hb_ij";
    } else {
      throw std::logic_error("unsupported resolved pair-coordinate family");
    }
    if (coordinate.identity != expected_identity) {
      throw std::invalid_argument(
          "resolved pair-coordinate identity does not match its application");
    }
    resolved.push_back({
        kind,
        component_indices.first,
        0,
        0,
        coordinate.value,
        0,
        component_indices.second,
        coordinate.association_edges,
        coordinate.association_energy_scales,
    });
  }
  return detail::fixed_topology_parameterized_mixture_phase_block(
      temperature, amounts, volume, resolved, parameters_);
}

// EqExclude: transport and buffer mechanics
NativeMixturePhaseBlock NativeModel::evaluate_electrolyte_phase(
    double temperature, const std::vector<double> &amounts, double volume,
    double ion_mole_fraction_max) const {
  return detail::electrolyte_phase_block(temperature, amounts, volume,
                                         ion_mole_fraction_max, parameters_);
}

// EqExclude: transport and buffer mechanics
NativeMixturePhaseBlock NativeModel::evaluate_reacting_phase_parameters(
    double temperature, const std::vector<double> &amounts, double volume,
    double ion_mole_fraction_max,
    const std::vector<NativeActiveParameter> &active_parameters) const {
  return detail::electrolyte_phase_block_active_parameters(
      temperature, amounts, volume, ion_mole_fraction_max, parameters_,
      active_parameters);
}

// EqExclude: transport and buffer mechanics
double NativeModel::evaluate_electrolyte_phase_value(
    double temperature, const std::vector<double> &amounts, double volume,
    double ion_mole_fraction_max) const {
  return evaluate_electrolyte_phase(temperature, amounts, volume,
                                    ion_mole_fraction_max)
      .helmholtz_over_rt_reference_amount;
}

// EqExclude: transport and buffer mechanics
NativePackingFractionBlock NativeModel::evaluate_packing_fraction(
    double temperature, const std::vector<double> &amounts, double volume,
    double ion_mole_fraction_max) const {
  return detail::packing_fraction_block(temperature, amounts, volume,
                                        ion_mole_fraction_max, parameters_);
}

// EqExclude: transport and buffer mechanics
NativePackingFractionBlock NativeModel::evaluate_packing_fraction_parameters(
    double temperature, const std::vector<double> &amounts, double volume,
    double ion_mole_fraction_max,
    const std::vector<NativeActiveParameter> &active_parameters) const {
  return detail::packing_fraction_block_active_parameters(
      temperature, amounts, volume, ion_mole_fraction_max, parameters_,
      active_parameters);
}

// EqExclude: transport and buffer mechanics
NativeInversePackingGeometryBlock
NativeModel::evaluate_inverse_packing_geometry(
    double temperature, const std::vector<double> &amounts,
    double log_packing_fraction,
    const std::vector<NativeActiveParameter> &active_parameters) const {
  return detail::inverse_packing_geometry_block(temperature, amounts,
                                                log_packing_fraction,
                                                parameters_, active_parameters);
}

// EqID: packing_volume_bounds
NativeMolarVolumeBoundsParameterBlock
NativeModel::molar_volume_bounds_parameters(
    double temperature, const std::vector<double> &mole_fractions,
    double packing_fraction_min, double packing_fraction_max,
    double ion_mole_fraction_max,
    const std::vector<NativeActiveParameter> &active_parameters) const {
  const NativeParameters active_values =
      detail::active_component_parameter_values(parameters_, active_parameters);
  const auto bounds = NativeModel(active_values)
                          .molar_volume_bounds(
                              temperature, mole_fractions, packing_fraction_min,
                              packing_fraction_max, ion_mole_fraction_max);
  const auto coefficient = detail::packing_fraction_block_active_parameters(
      temperature, mole_fractions, 1.0, ion_mole_fraction_max, parameters_,
      active_parameters);
  const std::size_t parameter_offset = mole_fractions.size() + 1;
  std::vector<double> derivatives(2 * active_parameters.size());
  for (std::size_t parameter = 0; parameter < active_parameters.size();
       ++parameter) {
    const double coefficient_derivative =
        coefficient.gradient[parameter_offset + parameter];
    derivatives[parameter] = coefficient_derivative / packing_fraction_max;
    derivatives[active_parameters.size() + parameter] =
        coefficient_derivative / packing_fraction_min;
  }
  return {bounds, std::move(derivatives)};
}

// EqID: fixed_pressure_total_derivative
NativeFugacityDerivativeBlock NativeModel::evaluate_fugacity_derivatives(
    double temperature, double pressure, const std::vector<double> &composition,
    const std::string &phase,
    const std::vector<NativeActiveParameter> &active_parameters) const {
  const NativeParameters active_values =
      detail::active_component_parameter_values(parameters_, active_parameters);
  const NativeModel active_model(active_values);
  const auto state =
      active_model.evaluate_pressure(temperature, pressure, composition, phase);
  if (!state.log_fugacity_coefficient.has_value() ||
      !state.density_diagnostics.has_value()) {
    throw std::runtime_error(
        "fixed-pressure fugacity diagnostics are unavailable");
  }
  const auto sigma = detail::segment_sigma(temperature, active_values);
  const auto diameter =
      detail::segment_diameter(temperature, active_values, sigma);
  const auto site_fractions =
      detail::solve_site_fractions(temperature, state.molar_density,
                                   composition, active_values, sigma, diameter);
  const auto derivatives = detail::thermodynamic_derivatives_active_components(
      temperature, state.molar_density, composition, active_values,
      site_fractions, active_parameters);
  const std::size_t components = composition.size();
  const std::size_t parameter_count = active_parameters.size();
  const std::size_t coordinates = 2 + components + parameter_count;
  const std::size_t parameter_offset = 2 + components;
  const auto hessian = [&derivatives, coordinates](std::size_t row,
                                                   std::size_t column) {
    return derivatives.hessian[row * coordinates + column];
  };
  std::vector<double> values(components);
  std::vector<double> pressure_derivatives(components);
  std::vector<double> parameter_derivatives(components * parameter_count);
  for (std::size_t component = 0; component < components; ++component) {
    const double fugacity_density = log_fugacity_density(
        state, derivatives, composition, component, coordinates);
    values[component] = fixed_pressure_log_fugacity(state, pressure, component,
                                                    fugacity_density);
    pressure_derivatives[component] =
        fugacity_density / state.pressure_density_derivative;
    for (std::size_t parameter = 0; parameter < parameter_count; ++parameter) {
      const std::size_t active = parameter_offset + parameter;
      double weighted_composition = 0.0;
      for (std::size_t mixture_component = 0; mixture_component < components;
           ++mixture_component) {
        weighted_composition += composition[mixture_component] *
                                hessian(2 + mixture_component, active);
      }
      const double pressure_parameter =
          state.molar_density * state.molar_density * detail::GAS_CONSTANT *
          temperature * hessian(1, active);
      const double density_parameter =
          -pressure_parameter / state.pressure_density_derivative;
      const double chemical_potential_parameter =
          derivatives.active_component_parameters[parameter] +
          state.molar_density * hessian(1, active) +
          hessian(2 + component, active) - weighted_composition;
      const double log_fugacity_parameter =
          chemical_potential_parameter - state.molar_density *
                                             hessian(1, active) /
                                             state.compressibility_factor;
      parameter_derivatives[component * parameter_count + parameter] =
          log_fugacity_parameter + fugacity_density * density_parameter;
    }
  }
  const auto &diagnostics = *state.density_diagnostics;
  const double condition =
      std::abs(pressure /
               (state.molar_density * diagnostics.pressure_density_derivative));
  if (!std::all_of(values.begin(), values.end(),
                   [](double item) { return std::isfinite(item); }) ||
      !std::all_of(pressure_derivatives.begin(), pressure_derivatives.end(),
                   [](double item) { return std::isfinite(item); }) ||
      !std::all_of(parameter_derivatives.begin(), parameter_derivatives.end(),
                   [](double item) { return std::isfinite(item); }) ||
      !std::isfinite(condition)) {
    throw std::runtime_error(
        "fixed-pressure fugacity derivative block is nonfinite");
  }
  return {
      std::move(values),
      std::move(pressure_derivatives),
      std::move(parameter_derivatives),
      state.molar_density,
      diagnostics.pressure_residual,
      diagnostics.pressure_density_derivative,
      condition,
      diagnostics.branch,
      diagnostics.bracket_lower,
      diagnostics.bracket_upper,
      diagnostics.stable_root_count,
      diagnostics.selected_stable_root_index,
  };
}

NativeFugacityDerivativeBlock NativeModel::evaluate_fugacity_derivatives_active(
    double temperature, double pressure, const std::vector<double> &composition,
    const std::string &phase,
    const std::vector<NativeParameterRecord> &active_parameters) const {
  if (active_parameters.empty()) {
    throw std::invalid_argument(
        "active fugacity derivative set must be nonempty");
  }
  std::vector<NativeActiveParameter> resolved;
  resolved.reserve(active_parameters.size());
  for (const auto &parameter : active_parameters) {
    resolved.push_back(resolve_component_parameter(parameter));
  }
  return evaluate_fugacity_derivatives(temperature, pressure, composition,
                                       phase, resolved);
}

// EqID: formula_unit_composition
std::vector<double> NativeModel::reference_composition(
    const std::vector<double> &solvent_fractions,
    const std::array<std::size_t, 2> &ion_indices,
    const std::array<int, 2> &stoichiometry, double molality) const {
  return detail::formula_composition(solvent_fractions,
                                     parameters_.molar_mass_kg_per_mol,
                                     ion_indices, stoichiometry, molality);
}

// EqID: reference_molality_sequence
NativeReferenceSequence NativeModel::electrolyte_reference_sequence(
    double temperature, double pressure,
    const std::vector<double> &solvent_fractions,
    const std::array<std::size_t, 2> &ion_indices,
    const std::array<int, 2> &stoichiometry, const std::string &phase,
    NativeEvaluationBudget *budget) const {
  auto checkpoint = [budget] {
    if (budget != nullptr) {
      budget->checkpoint();
    }
  };
  const auto result = detail::electrolyte_reference_sequence(
      solvent_fractions, parameters_.molar_mass_kg_per_mol, ion_indices,
      stoichiometry,
      [this, temperature, pressure, &ion_indices,
       &phase](const std::vector<double> &composition) {
        const auto state =
            evaluate_pressure(temperature, pressure, composition, phase);
        if (!state.log_fugacity_coefficient.has_value()) {
          throw std::runtime_error("reference fugacity is unavailable");
        }
        return detail::ReferenceObservation{
            {
                (*state.log_fugacity_coefficient)[ion_indices[0]],
                (*state.log_fugacity_coefficient)[ion_indices[1]],
            },
            std::nullopt,
        };
      },
      checkpoint);
  const auto pure_solvent_log_fugacity = single_solvent_log_fugacity(
      *this, temperature, pressure, solvent_fractions, phase);
  return {
      result.log_fugacity,
      pure_solvent_log_fugacity,
      result.molality,
      result.convergence_error,
  };
}

// EqID: infinite_dilution_solvation_energy
// EqID: fixed_pressure_total_derivative
NativeIonSolvationDerivativeBlock
NativeModel::evaluate_ion_solvation_derivatives(
    double temperature, double pressure,
    const std::vector<double> &solvent_fractions,
    const std::array<std::size_t, 2> &ion_indices,
    const std::array<int, 2> &stoichiometry, const std::string &phase,
    const std::vector<NativeParameterRecord> &active_parameters) const {
  const std::size_t components = parameters_.segment_count.size();
  if (!std::isfinite(temperature) || !std::isfinite(pressure) ||
      temperature <= 0.0 || pressure <= 0.0 || phase != "liquid" ||
      active_parameters.empty() || solvent_fractions.size() != components ||
      ion_indices[0] >= components || ion_indices[1] >= components ||
      ion_indices[0] == ion_indices[1] || stoichiometry[0] <= 0 ||
      stoichiometry[1] <= 0 ||
      parameters_.charge_number[ion_indices[0]] <= 0 ||
      parameters_.charge_number[ion_indices[1]] >= 0) {
    throw std::invalid_argument("invalid active ion-solvation reference inputs");
  }

  std::array<double, 2> values{};
  std::vector<double> derivatives(2 * active_parameters.size());
  double reference_molality = 0.0;
  double convergence_error = 0.0;
  for (std::size_t column = 0; column < active_parameters.size(); ++column) {
    const auto &record = active_parameters[column];
    if (!std::isfinite(record.value)) {
      throw std::invalid_argument("active parameter values must be finite");
    }
    const auto parts = split_identity(record.identity);
    NativeParameters updated = parameters_;
    ResolvedSolvationParameter resolved{
        ActiveSolvationParameter::BornDiameter, 0, std::nullopt};
    if (parts.size() == 3 && parts[0] == "component") {
      const std::size_t component =
          named_index(component_ids_, parts[1], "component");
      resolved.component = component;
      if (parts[2] == "born_diameter" && record.value > 0.0) {
        resolved.family = ActiveSolvationParameter::BornDiameter;
        updated.born_diameter_angstrom[component] = record.value;
      } else if (parts[2] == "solvation_factor" && record.value > 0.0) {
        resolved.family = ActiveSolvationParameter::SolvationFactor;
        updated.solvation_factor[component] = record.value;
      } else if ((parts[2] == "segment_diameter" ||
                  parts[2] == "dispersion_energy_over_k") &&
                 record.value > 0.0) {
        resolved.family = ActiveSolvationParameter::ComponentParameter;
        resolved.component_parameter = resolve_component_parameter(record);
        if (parts[2] == "segment_diameter") {
          updated.sigma_constant[component] = record.value;
        } else {
          updated.dispersion_energy_over_k[component] = record.value;
        }
      } else {
        throw std::logic_error("unsupported ion-solvation derivative identity: " +
                               record.identity);
      }
    } else if (parts.size() == 5 && parts[0] == "component" &&
               parts[2] == "relative_permittivity" &&
               parts[3] == "constant" && parts[4] == "constant" &&
               record.value > 1.0) {
      const std::size_t component =
          named_index(component_ids_, parts[1], "component");
      resolved.family = ActiveSolvationParameter::SolventRelativePermittivity;
      resolved.component = component;
      updated.relative_permittivity[component] = record.value;
    } else if (record.identity ==
                   "model/relative_permittivity/ionic_region_relative_permittivity" &&
               record.value > 0.0) {
      resolved.family =
          ActiveSolvationParameter::IonicRegionRelativePermittivity;
      updated.ionic_region_relative_permittivity = record.value;
    } else {
      throw std::logic_error("unsupported ion-solvation derivative identity: " +
                             record.identity);
    }

    const auto result = active_solvation_parameter_reference(
        std::move(updated), temperature, pressure, solvent_fractions,
        ion_indices, stoichiometry, phase, resolved);
    if (column == 0) {
      values = result.solvation_gibbs_j_per_mol;
      reference_molality = result.reference_molality_mol_per_kg;
      convergence_error = result.reference_convergence_error;
    }
    for (std::size_t ion = 0; ion < 2; ++ion) {
      derivatives[ion * active_parameters.size() + column] =
          result.parameter_derivatives_j_per_mol[ion];
    }
  }
  return {values, derivatives, reference_molality, convergence_error};
}

// EqID: infinite_dilution_solvation_energy
// EqID: fixed_pressure_total_derivative
NativeIonSolvationKij NativeModel::evaluate_ion_solvation_kij(
    double temperature, double pressure, std::size_t active_ion_component,
    const std::array<double, 3> &k_ij) const {
  if (!std::isfinite(temperature) || !std::isfinite(pressure) ||
      temperature <= 0.0 || pressure <= 0.0 ||
      (active_ion_component != 1 && active_ion_component != 2) ||
      !std::all_of(k_ij.begin(), k_ij.end(),
                   [](double value) { return std::isfinite(value); })) {
    throw std::invalid_argument("invalid active ion-solvation k_ij inputs");
  }
  if (parameters_.charge_number.size() != 3 ||
      parameters_.charge_number[0] != 0 || parameters_.charge_number[1] == 0 ||
      parameters_.charge_number[2] == 0 ||
      parameters_.charge_number[1] + parameters_.charge_number[2] != 0) {
    throw std::logic_error("unsupported active ion-solvation k_ij model");
  }
  NativeParameters active_parameters = parameters_;
  active_parameters.k_ij[1] = active_parameters.k_ij[3] = k_ij[0];
  active_parameters.k_ij[2] = active_parameters.k_ij[6] = k_ij[1];
  active_parameters.k_ij[5] = active_parameters.k_ij[7] = k_ij[2];
  const NativeModel active_model(active_parameters);
  const auto result = detail::electrolyte_reference_sequence(
      std::vector<double>{1.0, 0.0, 0.0},
      active_parameters.molar_mass_kg_per_mol, std::array<std::size_t, 2>{1, 2},
      std::array<int, 2>{1, 1},
      [&active_model, &active_parameters, temperature, pressure,
       &k_ij](const std::vector<double> &composition) {
        return active_aqueous_kij_observation(active_model, active_parameters,
                                              temperature, pressure,
                                              composition, k_ij);
      });
  if (!result.derivative.has_value() || result.derivative->size() != 6) {
    throw std::runtime_error(
        "active ion-solvation k_ij derivative is unavailable");
  }
  const std::size_t ion = active_ion_component - 1;
  std::array<double, 3> derivative{};
  for (std::size_t active = 0; active < 3; ++active) {
    derivative[active] = detail::GAS_CONSTANT * temperature *
                         (*result.derivative)[ion * 3 + active];
  }
  return {
      detail::GAS_CONSTANT * temperature * result.log_fugacity[ion],
      derivative,
      result.molality,
      result.convergence_error,
  };
}

// EqExclude: transport and buffer mechanics
NativeAqueousMiacKij NativeModel::evaluate_aqueous_miac_kij(
    double temperature, double pressure, double molality,
    const std::array<double, 3> &k_ij) const {
  return evaluate_aqueous_miac_kij_batch(temperature, pressure,
                                         std::vector<double>{molality}, k_ij)
      .front();
}

std::vector<NativeAqueousMiacKij>
// EqID: mean_ionic_activity_coefficient
// EqID: molality_activity_conversion
// EqID: fixed_pressure_total_derivative
NativeModel::evaluate_aqueous_miac_kij_batch(
    double temperature, double pressure, const std::vector<double> &molalities,
    const std::array<double, 3> &k_ij, NativeEvaluationBudget *budget) const {
  if (!std::isfinite(temperature) || !std::isfinite(pressure) ||
      temperature <= 0.0 || pressure <= 0.0 || molalities.empty() ||
      !positive_finite(molalities) ||
      !std::all_of(k_ij.begin(), k_ij.end(),
                   [](double value) { return std::isfinite(value); })) {
    throw std::invalid_argument("invalid active aqueous k_ij inputs");
  }
  if (parameters_.charge_number.size() != 3 ||
      parameters_.charge_number != std::vector<int>({0, 1, -1})) {
    throw std::logic_error("unsupported active aqueous k_ij model");
  }
  NativeParameters active_parameters = parameters_;
  active_parameters.k_ij[1] = active_parameters.k_ij[3] = k_ij[0];
  active_parameters.k_ij[2] = active_parameters.k_ij[6] = k_ij[1];
  active_parameters.k_ij[5] = active_parameters.k_ij[7] = k_ij[2];
  const NativeModel active_model(active_parameters);
  auto checkpoint = [budget] {
    if (budget != nullptr) {
      budget->checkpoint();
    }
  };
  const auto reference = detail::electrolyte_reference_sequence(
      std::vector<double>{1.0, 0.0, 0.0},
      active_parameters.molar_mass_kg_per_mol, std::array<std::size_t, 2>{1, 2},
      std::array<int, 2>{1, 1},
      [&active_model, &active_parameters, temperature, pressure,
       &k_ij](const std::vector<double> &reference_composition) {
        return active_aqueous_kij_observation(active_model, active_parameters,
                                              temperature, pressure,
                                              reference_composition, k_ij);
      },
      checkpoint);
  if (!reference.derivative.has_value() || reference.derivative->size() != 6 ||
      !reference.derivative_convergence_error.has_value()) {
    throw std::runtime_error("active aqueous k_ij derivative is unavailable");
  }
  std::vector<NativeAqueousMiacKij> results(molalities.size());
  detail::cppad_parallel_for(molalities.size(), [&](std::size_t index) {
    checkpoint();
    const double molality = molalities[index];
    const auto composition = detail::formula_composition(
        std::vector<double>{1.0, 0.0, 0.0},
        active_parameters.molar_mass_kg_per_mol,
        std::array<std::size_t, 2>{1, 2}, std::array<int, 2>{1, 1}, molality);
    const auto finite = active_aqueous_kij_observation(
        active_model, active_parameters, temperature, pressure, composition,
        k_ij);
    if (!finite.derivative.has_value() || finite.derivative->size() != 6) {
      throw std::runtime_error("active aqueous k_ij derivative is unavailable");
    }
    std::array<double, 3> total_derivative{};
    for (std::size_t active = 0; active < 3; ++active) {
      total_derivative[active] = 0.5 * ((*finite.derivative)[active] +
                                        (*finite.derivative)[3 + active] -
                                        (*reference.derivative)[active] -
                                        (*reference.derivative)[3 + active]);
    }
    const double conversion = std::log(
        1.0 + 2.0 * active_parameters.molar_mass_kg_per_mol[0] * molality);
    results[index] = {
        0.5 * (finite.log_fugacity[0] + finite.log_fugacity[1] -
               reference.log_fugacity[0] - reference.log_fugacity[1]) -
            conversion,
        total_derivative,
        reference.molality,
        reference.convergence_error,
        *reference.derivative_convergence_error,
    };
    const auto &result = results[index];
    if (!std::isfinite(result.log_mean_ionic_activity_coefficient_molality) ||
        !std::all_of(result.derivative.begin(), result.derivative.end(),
                     [](double value) { return std::isfinite(value); }) ||
        !std::isfinite(result.reference_molality_mol_per_kg) ||
        !std::isfinite(result.reference_convergence_error) ||
        !std::isfinite(result.reference_derivative_convergence_error)) {
      throw std::runtime_error(
          "active aqueous k_ij result contains nonfinite values");
    }
  });
  return results;
}

NativeAqueousMiacSolvationFactor
// EqExclude: transport and buffer mechanics
NativeModel::evaluate_aqueous_miac_solvation_factor(
    double temperature, double pressure, double molality,
    double solvation_factor) const {
  return evaluate_aqueous_miac_solvation_factor_batch(
             temperature, pressure, std::vector<double>{molality},
             solvation_factor)
      .front();
}

// EqExclude: transport and buffer mechanics
NativeIonFractionSuppression NativeModel::evaluate_ion_fraction_suppression(
    double total_ion_mole_fraction, double suppression_coefficient) const {
  if (!std::isfinite(total_ion_mole_fraction) ||
      !std::isfinite(suppression_coefficient) ||
      total_ion_mole_fraction <= 0.0 || total_ion_mole_fraction >= 1.0 ||
      suppression_coefficient <= 0.0) {
    throw std::invalid_argument(
        "invalid active ion-fraction suppression inputs");
  }
  if (parameters_.relative_permittivity_formulation !=
          RELATIVE_PERMITTIVITY_FORMULATION_ION_FRACTION_SUPPRESSION ||
      parameters_.charge_number.size() != 3 ||
      parameters_.charge_number[0] != 0 || parameters_.charge_number[1] <= 0 ||
      parameters_.charge_number[2] >= 0 ||
      !(parameters_.relative_permittivity[0] > 0.0)) {
    throw std::logic_error("unsupported active ion-fraction suppression model");
  }
  const double salt_free = parameters_.relative_permittivity[0];
  detail::CppAdSerialGuard cppad_guard;
  using ADScalar = CppAD::AD<double>;
  std::vector<ADScalar> active{suppression_coefficient};
  CppAD::Independent(active);
  const ADScalar permittivity =
      detail::ion_fraction_suppressed_relative_permittivity(
          ADScalar(salt_free), ADScalar(total_ion_mole_fraction), active[0]);
  CppAD::ADFun<double> function(active, {permittivity});
  const std::vector<double> point{suppression_coefficient};
  return {
      function.Forward(0, point)[0],
      function.Jacobian(point)[0],
      salt_free,
  };
}

std::vector<NativeAqueousMiacSolvationFactor>
// EqID: mean_ionic_activity_coefficient
// EqID: molality_activity_conversion
// EqID: fixed_pressure_total_derivative
NativeModel::evaluate_aqueous_miac_solvation_factor_batch(
    double temperature, double pressure, const std::vector<double> &molalities,
    double solvation_factor, NativeEvaluationBudget *budget) const {
  if (!std::isfinite(temperature) || !std::isfinite(pressure) ||
      !std::isfinite(solvation_factor) || temperature <= 0.0 ||
      pressure <= 0.0 || molalities.empty() || !positive_finite(molalities) ||
      solvation_factor <= 0.0) {
    throw std::invalid_argument("invalid active solvation-factor inputs");
  }
  if (parameters_.charge_number != std::vector<int>({0, 1, -1})) {
    throw std::logic_error("unsupported active solvation-factor model");
  }
  NativeParameters active_parameters = parameters_;
  active_parameters.solvation_factor[0] = solvation_factor;
  const NativeModel active_model(active_parameters);
  auto checkpoint = [budget] {
    if (budget != nullptr) {
      budget->checkpoint();
    }
  };
  const auto reference = detail::electrolyte_reference_sequence(
      std::vector<double>{1.0, 0.0, 0.0},
      active_parameters.molar_mass_kg_per_mol, std::array<std::size_t, 2>{1, 2},
      std::array<int, 2>{1, 1},
      [&active_model, &active_parameters, temperature, pressure,
       solvation_factor](const std::vector<double> &reference_composition) {
        return active_solvation_factor_observation(
            active_model, active_parameters, temperature, pressure,
            reference_composition, solvation_factor);
      },
      checkpoint);
  if (!reference.derivative.has_value() || reference.derivative->size() != 2 ||
      !reference.derivative_convergence_error.has_value()) {
    throw std::runtime_error(
        "active solvation-factor derivative is unavailable");
  }
  auto reference_log_fugacity = reference.log_fugacity;
  double reference_molality = reference.molality;
  double reference_convergence_error = reference.convergence_error;
  if (solvation_factor == parameters_.solvation_factor[0]) {
    const auto value_reference = active_model.electrolyte_reference_sequence(
        temperature, pressure, std::vector<double>{1.0, 0.0, 0.0},
        std::array<std::size_t, 2>{1, 2}, std::array<int, 2>{1, 1}, "liquid",
        budget);
    reference_log_fugacity = value_reference.log_fugacity;
    reference_molality = value_reference.reference_molality;
    reference_convergence_error = value_reference.convergence_error;
  }
  std::vector<NativeAqueousMiacSolvationFactor> results(molalities.size());
  detail::cppad_parallel_for(molalities.size(), [&](std::size_t index) {
    checkpoint();
    const double molality = molalities[index];
    const auto composition = detail::formula_composition(
        std::vector<double>{1.0, 0.0, 0.0},
        active_parameters.molar_mass_kg_per_mol,
        std::array<std::size_t, 2>{1, 2}, std::array<int, 2>{1, 1}, molality);
    const auto finite = active_solvation_factor_observation(
        active_model, active_parameters, temperature, pressure, composition,
        solvation_factor);
    if (!finite.derivative.has_value() || finite.derivative->size() != 2) {
      throw std::runtime_error(
          "active solvation-factor derivative is unavailable");
    }
    auto finite_log_fugacity = finite.log_fugacity;
    if (solvation_factor == parameters_.solvation_factor[0]) {
      const auto finite_state = active_model.evaluate_pressure(
          temperature, pressure, composition, "liquid");
      if (!finite_state.log_fugacity_coefficient.has_value()) {
        throw std::runtime_error("finite fugacity is unavailable");
      }
      finite_log_fugacity = {
          (*finite_state.log_fugacity_coefficient)[1],
          (*finite_state.log_fugacity_coefficient)[2],
      };
    }
    const double conversion = std::log(
        1.0 + 2.0 * active_parameters.molar_mass_kg_per_mol[0] * molality);
    results[index] = {
        0.5 * (finite_log_fugacity[0] + finite_log_fugacity[1] -
               reference_log_fugacity[0] - reference_log_fugacity[1]) -
            conversion,
        0.5 * ((*finite.derivative)[0] + (*finite.derivative)[1] -
               (*reference.derivative)[0] - (*reference.derivative)[1]),
        reference_molality,
        reference_convergence_error,
        *reference.derivative_convergence_error,
    };
    const auto &result = results[index];
    if (!std::isfinite(result.log_mean_ionic_activity_coefficient_molality) ||
        !std::isfinite(result.derivative) ||
        !std::isfinite(result.reference_molality_mol_per_kg) ||
        !std::isfinite(result.reference_convergence_error) ||
        !std::isfinite(result.reference_derivative_convergence_error)) {
      throw std::runtime_error(
          "active solvation-factor result contains nonfinite values");
    }
  });
  return results;
}

// EqID: packing_volume_bounds
std::array<double, 2> NativeModel::molar_volume_bounds(
    double temperature, const std::vector<double> &mole_fractions,
    double packing_fraction_min, double packing_fraction_max,
    double ion_mole_fraction_max) const {
  const std::size_t components = parameters_.segment_count.size();
  if (!std::isfinite(temperature) || temperature <= 0.0 ||
      mole_fractions.size() != components ||
      !std::isfinite(packing_fraction_min) ||
      !std::isfinite(packing_fraction_max) || packing_fraction_min <= 0.0 ||
      packing_fraction_min >= packing_fraction_max) {
    throw std::invalid_argument("invalid molar-volume bound inputs");
  }
  if (packing_fraction_max >= detail::close_packing_fraction()) {
    throw std::invalid_argument(
        "packing-fraction upper bound is outside the model domain");
  }

  double total_fraction = 0.0;
  double charge_fraction = 0.0;
  double ion_fraction = 0.0;
  for (std::size_t component = 0; component < components; ++component) {
    const double fraction = mole_fractions[component];
    if (!std::isfinite(fraction) || fraction < 0.0) {
      throw std::invalid_argument("invalid model composition");
    }
    total_fraction += fraction;
    const int charge = parameters_.charge_number[component];
    charge_fraction += fraction * static_cast<double>(charge);
    if (charge != 0) {
      ion_fraction += fraction;
    }
  }
  if (std::abs(total_fraction - 1.0) > 1.0e-12 ||
      std::abs(charge_fraction) > 1.0e-12 ||
      (std::isfinite(ion_mole_fraction_max) &&
       ion_fraction > ion_mole_fraction_max + 1.0e-12)) {
    throw std::invalid_argument("composition is outside the model domain");
  }

  const double coefficient = detail::packing_fraction_coefficient(
      temperature, mole_fractions, parameters_);
  if (!(coefficient > 0.0) || !std::isfinite(coefficient)) {
    throw std::invalid_argument("invalid model packing coefficient");
  }
  const std::array<double, 2> bounds = {
      coefficient / packing_fraction_max,
      coefficient / packing_fraction_min,
  };
  if (!(bounds[0] > 0.0) || !(bounds[0] < bounds[1]) ||
      !std::isfinite(bounds[1])) {
    throw std::invalid_argument("invalid model molar-volume bounds");
  }
  return bounds;
}

double NativeModel::molar_density_upper_bound(
    double temperature, const std::vector<double> &mole_fractions) const {
  const auto bounds =
      molar_volume_bounds(temperature, mole_fractions, 1.0e-12,
                          std::nextafter(detail::close_packing_fraction(), 0.0),
                          std::numeric_limits<double>::quiet_NaN());
  return 1.0 / bounds[0];
}

NativeResult
NativeModel::evaluate(double temperature, double molar_density,
                      const std::vector<double> &mole_fractions) const {
  const std::size_t components = parameters_.segment_count.size();
  if (!std::isfinite(temperature) || temperature <= 0.0 ||
      !std::isfinite(molar_density) || molar_density <= 0.0 ||
      mole_fractions.size() != components || !finite(mole_fractions)) {
    throw std::invalid_argument("invalid native state");
  }
  if (!polar_composition_admitted(parameters_, mole_fractions)) {
    throw std::invalid_argument("polar model mole fractions must be positive");
  }
  const auto sigma = detail::segment_sigma(temperature, parameters_);
  if (!positive_finite(sigma)) {
    throw std::invalid_argument("invalid segment diameter correlation");
  }
  const auto diameter =
      detail::segment_diameter(temperature, parameters_, sigma);
  const auto site_fraction = detail::solve_site_fractions(
      temperature, molar_density, mole_fractions, parameters_, sigma, diameter);
  const auto electrolyte = detail::electrolyte_contribution(
      temperature, molar_density, temperature, molar_density, mole_fractions,
      mole_fractions, parameters_, diameter, diameter);
  const double born_ssm_ds =
      detail::born_contribution(temperature, mole_fractions, parameters_,
                                electrolyte.bulk_relative_permittivity, sigma);
  const auto values = detail::residual_helmholtz(
      temperature, molar_density, mole_fractions, parameters_, sigma, diameter,
      site_fraction, electrolyte.debye_huckel, born_ssm_ds);

  const double residual = values.hard_chain + values.dispersion +
                          values.association + values.polar +
                          values.debye_huckel + values.born_ssm_ds;
  const auto derivatives = detail::thermodynamic_derivatives(
      temperature, molar_density, mole_fractions, parameters_, site_fraction);
  // EqID: compressibility_factor
  const double compressibility = 1.0 + molar_density * derivatives.density;
  double weighted_composition_derivative = 0.0;
  for (std::size_t component = 0; component < components; ++component) {
    weighted_composition_derivative +=
        mole_fractions[component] * derivatives.composition[component];
  }
  std::vector<double> chemical_potential(components);
  std::optional<std::vector<double>> log_fugacity;
  if (std::isfinite(compressibility) && compressibility > 0.0) {
    log_fugacity.emplace(components);
  }
  for (std::size_t component = 0; component < components; ++component) {
    // EqID: residual_chemical_potential
    chemical_potential[component] = residual + compressibility - 1.0 +
                                    derivatives.composition[component] -
                                    weighted_composition_derivative;
    if (log_fugacity.has_value()) {
      // EqID: log_fugacity_coefficient
      (*log_fugacity)[component] =
          chemical_potential[component] - std::log(compressibility);
    }
  }
  const double z_temperature = molar_density * derivatives.temperature_density;
  const double z_density =
      derivatives.density + molar_density * derivatives.density_density;
  const double pressure_slope_factor =
      compressibility + molar_density * z_density;
  const double cv_over_r =
      -2.0 * temperature * derivatives.temperature -
      temperature * temperature * derivatives.temperature_temperature;
  std::optional<double> residual_cp;
  if (std::isfinite(pressure_slope_factor) && pressure_slope_factor > 1.0e-10) {
    // EqID: residual_heat_capacity_cp
    const double cp_over_r =
        cv_over_r +
        (compressibility + temperature * z_temperature) *
            (compressibility + temperature * z_temperature) /
            pressure_slope_factor -
        1.0;
    if (std::isfinite(cp_over_r)) {
      residual_cp = detail::GAS_CONSTANT * cp_over_r;
    }
  }
  return {
      values.hard_chain,
      values.dispersion,
      values.association,
      values.polar,
      values.debye_huckel,
      values.born_ssm_ds,
      electrolyte.bulk_relative_permittivity,
      residual,
      compressibility,
      // EqID: pressure
      compressibility * molar_density * detail::GAS_CONSTANT * temperature,
      std::move(chemical_potential),
      std::move(log_fugacity),
      // EqID: residual_entropy
      detail::GAS_CONSTANT *
          (-residual - temperature * derivatives.temperature),
      // EqID: residual_internal_energy
      -detail::GAS_CONSTANT * temperature * temperature *
          derivatives.temperature,
      // EqID: residual_enthalpy
      detail::GAS_CONSTANT * temperature *
          (-temperature * derivatives.temperature + compressibility - 1.0),
      // EqID: residual_heat_capacity_cv
      detail::GAS_CONSTANT * cv_over_r,
      residual_cp,
      molar_density,
      // EqID: pressure_density_derivative
      detail::GAS_CONSTANT * temperature * pressure_slope_factor,
      std::nullopt,
  };
}

NativeResult
NativeModel::evaluate_pressure(double temperature, double pressure,
                               const std::vector<double> &mole_fractions,
                               const std::string &phase) const {
  const std::size_t components = parameters_.segment_count.size();
  if (!std::isfinite(temperature) || temperature <= 0.0 ||
      !std::isfinite(pressure) || pressure <= 0.0 ||
      mole_fractions.size() != components || !finite(mole_fractions)) {
    throw std::invalid_argument("invalid native pressure state");
  }
  if (!polar_composition_admitted(parameters_, mole_fractions)) {
    throw std::invalid_argument("polar model mole fractions must be positive");
  }
  // EqID: density_upper_bound
  const double packing_coefficient = detail::packing_fraction_coefficient(
      temperature, mole_fractions, parameters_);
  if (!(packing_coefficient > 0.0) || !std::isfinite(packing_coefficient)) {
    throw std::invalid_argument("invalid pressure-state packing coefficient");
  }
  const auto root = detail::solve_density(
      pressure,
      // EqID: ideal_gas_density
      pressure / (detail::GAS_CONSTANT * temperature),
      detail::close_packing_fraction() / packing_coefficient, phase,
      [this, temperature, &mole_fractions](double density) {
        return evaluate(temperature, density, mole_fractions);
      });
  auto result = evaluate(temperature, root.molar_density, mole_fractions);
  result.density_diagnostics = NativeDensityDiagnostics{
      root.bracket_lower,
      root.bracket_upper,
      root.pressure_residual,
      root.pressure_density_derivative,
      root.iterations,
      root.stable,
      root.branch,
      root.stable_root_count,
      root.selected_stable_root_index,
  };
  return result;
}

} // namespace epcsaft
