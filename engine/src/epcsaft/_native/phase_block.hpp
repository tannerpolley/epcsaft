#pragma once

#include "association.hpp"
#include "derivatives.hpp"
#include "helmholtz.hpp"
#include "model.hpp"

#include <cppad/cppad.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace epcsaft {

struct NativePurePhaseBlock final {
  double helmholtz_over_rt_reference_amount;
  std::array<double, 2> gradient;
  std::array<double, 4> hessian;
  std::array<double, 8> third;
  double pressure_pa;
  double chemical_potential_over_rt;
};

enum class NativeFixedTopologyParameterKind {
  SegmentCount,
  SegmentDiameter,
  SegmentDiameterCorrelationConstant,
  SegmentDiameterCorrelationAmplitude,
  SegmentDiameterCorrelationExponent,
  DispersionEnergyOverK,
  AssociationEnergyOverK,
  AssociationVolume,
};

struct NativeFixedTopologyParameter final {
  NativeFixedTopologyParameterKind kind;
  std::size_t component;
  std::size_t site_a;
  std::size_t site_b;
  double value;
  std::size_t correlation_term{0};
};

struct NativeFixedTopologyPhaseBlock final {
  double helmholtz_over_rt_reference_amount;
  std::vector<double> gradient;
  std::vector<double> hessian;
  double pressure_pa;
  double chemical_potential_over_rt;
};

struct NativeMixturePhaseBlock final {
  double helmholtz_over_rt_reference_amount;
  std::vector<double> gradient;
  std::vector<double> hessian;
  double pressure_pa;
};

struct NativePackingFractionBlock final {
  double packing_fraction;
  std::vector<double> gradient;
  std::vector<double> hessian;
};

struct NativeInversePackingGeometryBlock final {
  double volume_m3;
  std::vector<double> gradient;
  std::vector<double> hessian;
};

namespace detail {

// EqExclude: input and domain validation
inline bool supports_pure_phase_block(const NativeParameters &parameters) {
  return parameters.segment_count.size() == 1 &&
         parameters.site_component.empty() &&
         parameters.charge_number.size() == 1 &&
         parameters.charge_number[0] == 0;
}

// EqExclude: input and domain validation
inline bool supports_fixed_topology_associating_phase_block(
    const NativeParameters &parameters) {
  const std::size_t components = parameters.segment_count.size();
  const std::size_t sites = parameters.site_component.size();
  if (!(components > 0 && parameters.charge_number.size() == components &&
        sites > 0 && parameters.site_multiplicity.size() == sites &&
        parameters.association_active.size() == sites * sites &&
        parameters.association_energy_over_k.size() == sites * sites &&
        parameters.association_volume.size() == sites * sites &&
        std::all_of(parameters.site_component.begin(),
                    parameters.site_component.end(),
                    [components](int component) {
                      return component >= 0 &&
                             static_cast<std::size_t>(component) < components;
                    }) &&
        std::all_of(parameters.site_multiplicity.begin(),
                    parameters.site_multiplicity.end(),
                    [](int multiplicity) { return multiplicity > 0; }) &&
        std::any_of(parameters.association_active.begin(),
                    parameters.association_active.end(),
                    [](int active) { return active != 0; }))) {
    return false;
  }
  for (std::size_t left = 0; left < sites; ++left) {
    for (std::size_t right = left; right < sites; ++right) {
      const std::size_t forward = left * sites + right;
      const std::size_t reverse = right * sites + left;
      if (parameters.association_active[forward] !=
              parameters.association_active[reverse] ||
          parameters.association_energy_over_k[forward] !=
              parameters.association_energy_over_k[reverse] ||
          parameters.association_volume[forward] !=
              parameters.association_volume[reverse]) {
        return false;
      }
    }
  }
  return true;
}

// EqExclude: input and domain validation
inline bool supports_fixed_topology_associating_pure_phase_block(
    const NativeParameters &parameters) {
  return parameters.segment_count.size() == 1 &&
         supports_fixed_topology_associating_phase_block(parameters);
}

// EqExclude: input and domain validation
inline bool supports_fixed_topology_associating_mixture_phase_block(
    const NativeParameters &parameters) {
  return parameters.segment_count.size() >= 2 &&
         supports_fixed_topology_associating_phase_block(parameters);
}

// EqExclude: input and domain validation
inline bool supports_mixture_phase_block(const NativeParameters &parameters) {
  return parameters.segment_count.size() >= 2 &&
         parameters.site_component.empty() &&
         std::all_of(parameters.charge_number.begin(),
                     parameters.charge_number.end(),
                     [](int charge) { return charge == 0; });
}

// EqExclude: input and domain validation
inline bool
supports_binary_interaction_phase_block(const NativeParameters &parameters) {
  return parameters.segment_count.size() == 2 &&
         supports_mixture_phase_block(parameters);
}

// EqExclude: input and domain validation
inline bool
supports_electrolyte_phase_block(const NativeParameters &parameters) {
  return parameters.segment_count.size() >= 3 &&
         !parameters.site_component.empty() &&
         std::any_of(parameters.charge_number.begin(),
                     parameters.charge_number.end(),
                     [](int charge) { return charge != 0; });
}

// EqExclude: input and domain validation
inline bool supports_packing_geometry(const NativeParameters &parameters) {
  const std::size_t components = parameters.segment_count.size();
  return components != 0 && parameters.sigma_constant.size() == components &&
         parameters.dispersion_energy_over_k.size() == components &&
         parameters.charge_number.size() == components &&
         parameters.sigma_term_component.size() ==
             parameters.sigma_term_amplitude.size() &&
         parameters.sigma_term_component.size() ==
             parameters.sigma_term_exponent.size() &&
         std::all_of(parameters.sigma_term_component.begin(),
                     parameters.sigma_term_component.end(),
                     [components](int component) {
                       return component >= 0 &&
                              static_cast<std::size_t>(component) < components;
                     }) &&
         std::all_of(parameters.segment_count.begin(),
                     parameters.segment_count.end(),
                     [](double value) {
                       return std::isfinite(value) && value > 0.0;
                     }) &&
         std::all_of(parameters.sigma_constant.begin(),
                     parameters.sigma_constant.end(),
                     [](double value) {
                       return std::isfinite(value) && value > 0.0;
                     }) &&
         std::all_of(
             parameters.dispersion_energy_over_k.begin(),
             parameters.dispersion_energy_over_k.end(),
             [](double value) { return std::isfinite(value) && value > 0.0; });
}

template <typename Scalar> struct PhaseBlockParameters final {
  std::vector<Scalar> segment_count;
  std::vector<Scalar> sigma_constant;
  std::vector<int> sigma_term_component;
  std::vector<Scalar> sigma_term_amplitude;
  std::vector<Scalar> sigma_term_exponent;
  std::vector<Scalar> dispersion_energy_over_k;
  std::vector<Scalar> k_ij;
  std::vector<Scalar> l_ij;
  std::vector<int> site_component;
  std::vector<int> site_multiplicity;
  std::vector<int> association_active;
  std::vector<int> charge_number;
  int electrolyte_model;
  std::vector<double> dipole_moment_debye;
  std::vector<double> quadrupole_moment_debye_angstrom;
  int polar_formulation;
};

template <typename Scalar>
// EqExclude: scalar utility dispatch
std::vector<Scalar> phase_block_scalars(const std::vector<double> &values) {
  return std::vector<Scalar>(values.begin(), values.end());
}

template <typename Scalar>
// EqExclude: scalar utility dispatch
PhaseBlockParameters<Scalar>
phase_block_parameters(const NativeParameters &parameters) {
  return {
      phase_block_scalars<Scalar>(parameters.segment_count),
      phase_block_scalars<Scalar>(parameters.sigma_constant),
      parameters.sigma_term_component,
      phase_block_scalars<Scalar>(parameters.sigma_term_amplitude),
      phase_block_scalars<Scalar>(parameters.sigma_term_exponent),
      phase_block_scalars<Scalar>(parameters.dispersion_energy_over_k),
      phase_block_scalars<Scalar>(parameters.k_ij),
      phase_block_scalars<Scalar>(parameters.l_ij),
      parameters.site_component,
      parameters.site_multiplicity,
      parameters.association_active,
      parameters.charge_number,
      parameters.electrolyte_model,
      parameters.dipole_moment_debye,
      parameters.quadrupole_moment_debye_angstrom,
      parameters.polar_formulation,
  };
}

template <typename Scalar>
// EqID: active_binary_interaction_mapping
PhaseBlockParameters<Scalar>
binary_interaction_parameters(const NativeParameters &parameters,
                              const Scalar &k_ij) {
  auto values = phase_block_parameters<Scalar>(parameters);
  values.k_ij = {Scalar(parameters.k_ij[0]), k_ij, k_ij,
                 Scalar(parameters.k_ij[3])};
  return values;
}

template <typename Scalar>
// EqID: active_binary_size_interaction_mapping
PhaseBlockParameters<Scalar>
binary_size_interaction_parameters(const NativeParameters &parameters,
                                   const Scalar &l_ij) {
  auto values = phase_block_parameters<Scalar>(parameters);
  values.l_ij = {Scalar(parameters.l_ij[0]), l_ij, l_ij,
                 Scalar(parameters.l_ij[3])};
  return values;
}

template <typename Scalar, typename Parameters>
Scalar mixture_phase_helmholtz_over_rt_reference_amount(
    double temperature, const std::vector<Scalar> &amounts,
    const Scalar &volume, const Parameters &parameters) {
  // EqID: phase_total_amount
  Scalar total = Scalar(0.0);
  for (const Scalar &amount : amounts) {
    total += amount;
  }
  // EqID: phase_molar_density
  const Scalar density = total / volume;
  std::vector<Scalar> composition;
  composition.reserve(amounts.size());
  for (const Scalar &amount : amounts) {
    // EqID: phase_mole_fraction
    composition.push_back(amount / total);
  }
  const auto sigma = segment_sigma(Scalar(temperature), parameters);
  const auto diameter =
      segment_diameter(Scalar(temperature), parameters, sigma);
  const auto contributions =
      residual_helmholtz(Scalar(temperature), density, composition, parameters,
                         sigma, diameter, {}, Scalar(0.0), Scalar(0.0));
  Scalar ideal = Scalar(0.0);
  for (const Scalar &amount : amounts) {
    ideal += amount * (scalar_log(amount / volume) - Scalar(1.0));
  }
  // EqID: phase_dimensionless_helmholtz
  return ideal + total * (contributions.hard_chain + contributions.dispersion +
                          contributions.polar);
}

// EqExclude: input and domain validation
inline void validate_mixture_phase_state(double temperature,
                                         const std::vector<double> &amounts,
                                         double volume,
                                         const NativeParameters &parameters) {
  if (!std::isfinite(temperature) || temperature <= 0.0 ||
      !std::isfinite(volume) || volume <= 0.0 ||
      amounts.size() != parameters.segment_count.size() ||
      !std::all_of(amounts.begin(), amounts.end(), [](double amount) {
        return std::isfinite(amount) && amount > 0.0;
      })) {
    throw std::invalid_argument(
        "temperature, amounts, and volume must be positive and finite");
  }
}

// EqExclude: input and domain validation
inline void validate_electrolyte_phase_state(
    double temperature, const std::vector<double> &amounts, double volume,
    double ion_mole_fraction_max, const NativeParameters &parameters) {
  validate_mixture_phase_state(temperature, amounts, volume, parameters);
  double total_amount = 0.0;
  double ionic_amount = 0.0;
  double charge_amount = 0.0;
  for (std::size_t component = 0; component < amounts.size(); ++component) {
    total_amount += amounts[component];
    const int charge = parameters.charge_number[component];
    if (charge != 0) {
      ionic_amount += amounts[component];
    }
    charge_amount += amounts[component] * static_cast<double>(charge);
  }
  if (std::abs(charge_amount / total_amount) > 1.0e-12) {
    throw std::invalid_argument(
        "electrolyte phase amounts must be electroneutral");
  }
  if (std::isfinite(ion_mole_fraction_max) &&
      ionic_amount / total_amount > ion_mole_fraction_max + 1.0e-12) {
    throw std::invalid_argument(
        "ion mole fraction exceeds the parameter source domain");
  }
}

template <typename PhaseFunction>
// EqID: phase_pressure
// EqExclude: CppAD-generated derivative machinery
NativeMixturePhaseBlock differentiated_mixture_phase_block(
    double temperature, const std::vector<double> &point,
    std::size_t volume_coordinate, PhaseFunction phase_function) {
  CppAdSerialGuard cppad_guard;
  using ADScalar = CppAD::AD<double>;
  std::vector<ADScalar> independent(point.begin(), point.end());
  CppAD::Independent(independent);
  const ADScalar phase = phase_function(independent);
  CppAD::ADFun<double> function(independent, {phase});

  const auto gradient = function.Jacobian(point);
  const auto hessian = function.Hessian(point, 0);
  const double value = function.Forward(0, point)[0];
  if (!std::isfinite(value) ||
      !std::all_of(gradient.begin(), gradient.end(),
                   [](double item) { return std::isfinite(item); }) ||
      !std::all_of(hessian.begin(), hessian.end(),
                   [](double item) { return std::isfinite(item); })) {
    throw std::runtime_error("mixture phase block contains nonfinite values");
  }

  return {
      value,
      gradient,
      hessian,
      -GAS_CONSTANT * temperature * gradient[volume_coordinate],
  };
}

// EqExclude: transport and buffer mechanics
inline NativeMixturePhaseBlock
mixture_phase_block(double temperature, const std::vector<double> &amounts,
                    double volume, const NativeParameters &parameters) {
  if (!supports_mixture_phase_block(parameters)) {
    throw std::logic_error("mixture phase block requires at least two neutral "
                           "nonassociating components");
  }
  validate_mixture_phase_state(temperature, amounts, volume, parameters);

  mixture_phase_helmholtz_over_rt_reference_amount(temperature, amounts, volume,
                                                   parameters);

  std::vector<double> point(amounts.begin(), amounts.end());
  point.push_back(volume);
  const std::size_t volume_coordinate = amounts.size();
  return differentiated_mixture_phase_block(
      temperature, point, volume_coordinate, [&](const auto &independent) {
        using ADScalar = CppAD::AD<double>;
        const std::vector<ADScalar> active_amounts(
            independent.begin(), independent.begin() + volume_coordinate);
        return mixture_phase_helmholtz_over_rt_reference_amount(
            temperature, active_amounts, independent[volume_coordinate],
            parameters);
      });
}

// EqExclude: transport and buffer mechanics
inline NativeMixturePhaseBlock
mixture_phase_block_kij(double temperature, const std::vector<double> &amounts,
                        double volume, double k_ij,
                        const NativeParameters &parameters) {
  if (!supports_binary_interaction_phase_block(parameters)) {
    throw std::logic_error("active k_ij phase block requires two neutral "
                           "nonassociating components");
  }
  if (!std::isfinite(k_ij)) {
    throw std::invalid_argument("k_ij must be finite");
  }
  validate_mixture_phase_state(temperature, amounts, volume, parameters);
  const auto value_parameters = binary_interaction_parameters(parameters, k_ij);
  mixture_phase_helmholtz_over_rt_reference_amount(temperature, amounts, volume,
                                                   value_parameters);

  std::vector<double> point(amounts.begin(), amounts.end());
  point.push_back(volume);
  point.push_back(k_ij);
  const std::size_t volume_coordinate = amounts.size();
  return differentiated_mixture_phase_block(
      temperature, point, volume_coordinate, [&](const auto &independent) {
        using ADScalar = CppAD::AD<double>;
        const std::vector<ADScalar> active_amounts(
            independent.begin(), independent.begin() + volume_coordinate);
        const auto active_parameters = binary_interaction_parameters(
            parameters, independent[volume_coordinate + 1]);
        return mixture_phase_helmholtz_over_rt_reference_amount(
            temperature, active_amounts, independent[volume_coordinate],
            active_parameters);
      });
}

// EqExclude: transport and buffer mechanics
inline NativeMixturePhaseBlock
mixture_phase_block_lij(double temperature, const std::vector<double> &amounts,
                        double volume, double l_ij,
                        const NativeParameters &parameters) {
  if (!supports_binary_interaction_phase_block(parameters)) {
    throw std::logic_error("active l_ij phase block requires two neutral "
                           "nonassociating components");
  }
  if (!std::isfinite(l_ij)) {
    throw std::invalid_argument("l_ij must be finite");
  }
  validate_mixture_phase_state(temperature, amounts, volume, parameters);
  const auto value_parameters =
      binary_size_interaction_parameters(parameters, l_ij);
  mixture_phase_helmholtz_over_rt_reference_amount(temperature, amounts, volume,
                                                   value_parameters);

  std::vector<double> point(amounts.begin(), amounts.end());
  point.push_back(volume);
  point.push_back(l_ij);
  const std::size_t volume_coordinate = amounts.size();
  return differentiated_mixture_phase_block(
      temperature, point, volume_coordinate, [&](const auto &independent) {
        using ADScalar = CppAD::AD<double>;
        const std::vector<ADScalar> active_amounts(
            independent.begin(), independent.begin() + volume_coordinate);
        const auto active_parameters = binary_size_interaction_parameters(
            parameters, independent[volume_coordinate + 1]);
        return mixture_phase_helmholtz_over_rt_reference_amount(
            temperature, active_amounts, independent[volume_coordinate],
            active_parameters);
      });
}

// EqExclude: transport and buffer mechanics
inline NativeParameters active_component_parameter_values(
    const NativeParameters &parameters,
    const std::vector<NativeActiveParameter> &active_parameters) {
  NativeParameters values = parameters;
  for (const auto &active : active_parameters) {
    if (active.component >= values.segment_count.size() ||
        !std::isfinite(active.value) || active.value <= 0.0) {
      throw std::invalid_argument(
          "active component parameter is outside its domain");
    }
    if (active.family == ActiveComponentParameter::SegmentDiameter) {
      values.sigma_constant[active.component] = active.value;
    } else {
      values.dispersion_energy_over_k[active.component] = active.value;
    }
  }
  return values;
}

template <typename Scalar>
// EqExclude: scalar utility dispatch
std::vector<Scalar> active_packing_diameters(
    double temperature, const NativeParameters &parameters,
    const std::vector<NativeActiveParameter> &active_parameters,
    const std::vector<Scalar> &active_values) {
  auto sigma = segment_sigma(Scalar(temperature), parameters);
  std::vector<Scalar> dispersion(parameters.dispersion_energy_over_k.begin(),
                                 parameters.dispersion_energy_over_k.end());
  for (std::size_t coordinate = 0; coordinate < active_parameters.size();
       ++coordinate) {
    const auto &parameter = active_parameters[coordinate];
    if (parameter.family == ActiveComponentParameter::SegmentDiameter) {
      sigma[parameter.component] +=
          active_values[coordinate] -
          parameters.sigma_constant[parameter.component];
    } else {
      dispersion[parameter.component] = active_values[coordinate];
    }
  }
  return segment_diameter(Scalar(temperature), parameters, sigma, &dispersion);
}

// EqID: electrolyte_phase_coordinate_jacobian
// EqID: electrolyte_phase_hessian_chain_rule
inline NativeMixturePhaseBlock electrolyte_phase_block_impl(
    double temperature, const std::vector<double> &amounts, double volume,
    double ion_mole_fraction_max, const NativeParameters &parameters,
    const std::vector<NativeActiveParameter> &active_parameters) {
  if (!supports_electrolyte_phase_block(parameters)) {
    throw std::logic_error(
        "electrolyte phase block requires an ionic associating model");
  }
  const NativeParameters value_parameters =
      active_component_parameter_values(parameters, active_parameters);
  validate_electrolyte_phase_state(temperature, amounts, volume,
                                   ion_mole_fraction_max, value_parameters);
  // Keep the analytic phase objective available to constrained-optimizer
  // trial points. Packing admissibility is owned by the separate packing
  // and molar-volume-bound callbacks.

  const std::size_t components = amounts.size();
  const std::size_t parameter_count = active_parameters.size();
  const std::size_t coordinates = components + 1 + parameter_count;
  const std::size_t state_count = 1 + components + parameter_count;
  const std::size_t derivative_coordinates = 2 + components + parameter_count;
  // EqID: phase_total_amount
  double total_amount = 0.0;
  for (std::size_t component = 0; component < components; ++component) {
    total_amount += amounts[component];
  }

  // EqID: phase_molar_density
  const double density = total_amount / volume;
  std::vector<double> composition(components);
  for (std::size_t component = 0; component < components; ++component) {
    // EqID: phase_mole_fraction
    composition[component] = amounts[component] / total_amount;
  }
  const auto sigma = segment_sigma(temperature, value_parameters);
  const auto diameter = segment_diameter(temperature, value_parameters, sigma);
  const auto site_fractions = solve_site_fractions(
      temperature, density, composition, value_parameters, sigma, diameter);
  const auto derivatives =
      active_parameters.empty()
          ? thermodynamic_derivatives(temperature, density, composition,
                                      value_parameters, site_fractions)
          : thermodynamic_derivatives_active_components(
                temperature, density, composition, value_parameters,
                site_fractions, active_parameters);

  std::vector<double> state_gradient(state_count);
  state_gradient[0] = derivatives.density;
  std::copy(derivatives.composition.begin(), derivatives.composition.end(),
            state_gradient.begin() + 1);
  std::copy(derivatives.active_component_parameters.begin(),
            derivatives.active_component_parameters.end(),
            state_gradient.begin() + 1 +
                static_cast<std::ptrdiff_t>(components));
  std::vector<double> state_jacobian(state_count * coordinates, 0.0);
  for (std::size_t amount = 0; amount < components; ++amount) {
    state_jacobian[amount] = 1.0 / volume;
    for (std::size_t component = 0; component < components; ++component) {
      state_jacobian[(1 + component) * coordinates + amount] =
          ((component == amount ? 1.0 : 0.0) - composition[component]) /
          total_amount;
    }
  }
  state_jacobian[components] = -total_amount / (volume * volume);
  for (std::size_t parameter = 0; parameter < parameter_count; ++parameter) {
    state_jacobian[(1 + components + parameter) * coordinates + components + 1 +
                   parameter] = 1.0;
  }

  std::vector<double> gradient(coordinates, 0.0);
  std::vector<double> hessian(coordinates * coordinates, 0.0);
  double ideal = 0.0;
  for (std::size_t amount = 0; amount < components; ++amount) {
    ideal += amounts[amount] * (std::log(amounts[amount] / volume) - 1.0);
    gradient[amount] = std::log(amounts[amount] / volume);
    hessian[amount * coordinates + amount] = 1.0 / amounts[amount];
    hessian[amount * coordinates + components] = -1.0 / volume;
    hessian[components * coordinates + amount] = -1.0 / volume;
  }
  gradient[components] = -total_amount / volume;
  hessian[components * coordinates + components] =
      total_amount / (volume * volume);

  auto state_second = [&](std::size_t state, std::size_t left,
                          std::size_t right) {
    if (state > components) {
      return 0.0;
    }
    if (state == 0) {
      if (left == components && right == components) {
        return 2.0 * total_amount / (volume * volume * volume);
      }
      if ((left == components && right < components) ||
          (right == components && left < components)) {
        return -1.0 / (volume * volume);
      }
      return 0.0;
    }
    if (left >= components || right >= components) {
      return 0.0;
    }
    const std::size_t component = state - 1;
    return (2.0 * composition[component] - (component == left ? 1.0 : 0.0) -
            (component == right ? 1.0 : 0.0)) /
           (total_amount * total_amount);
  };

  for (std::size_t left = 0; left < coordinates; ++left) {
    const double total_left = left < components ? 1.0 : 0.0;
    double residual_gradient = total_left * derivatives.residual_helmholtz;
    for (std::size_t state = 0; state < state_count; ++state) {
      residual_gradient += total_amount * state_gradient[state] *
                           state_jacobian[state * coordinates + left];
    }
    gradient[left] += residual_gradient;

    for (std::size_t right = 0; right < coordinates; ++right) {
      const double total_right = right < components ? 1.0 : 0.0;
      double residual_hessian = 0.0;
      for (std::size_t state = 0; state < state_count; ++state) {
        const double left_derivative =
            state_jacobian[state * coordinates + left];
        const double right_derivative =
            state_jacobian[state * coordinates + right];
        residual_hessian +=
            state_gradient[state] *
            (total_left * right_derivative + total_right * left_derivative +
             total_amount * state_second(state, left, right));
        for (std::size_t other = 0; other < state_count; ++other) {
          residual_hessian +=
              total_amount *
              derivatives
                  .hessian[(1 + state) * derivative_coordinates + 1 + other] *
              left_derivative * state_jacobian[other * coordinates + right];
        }
      }
      hessian[left * coordinates + right] += residual_hessian;
    }
  }

  // EqID: phase_dimensionless_helmholtz
  const double value = ideal + total_amount * derivatives.residual_helmholtz;
  if (!std::isfinite(value) ||
      !std::all_of(gradient.begin(), gradient.end(),
                   [](double item) { return std::isfinite(item); }) ||
      !std::all_of(hessian.begin(), hessian.end(),
                   [](double item) { return std::isfinite(item); })) {
    throw std::runtime_error(
        "electrolyte phase block contains nonfinite values");
  }
  // EqID: phase_pressure
  const double pressure = -GAS_CONSTANT * temperature * gradient[components];
  return {
      value,
      std::move(gradient),
      std::move(hessian),
      pressure,
  };
}

// EqExclude: transport and buffer mechanics
inline NativeMixturePhaseBlock
electrolyte_phase_block(double temperature, const std::vector<double> &amounts,
                        double volume, double ion_mole_fraction_max,
                        const NativeParameters &parameters) {
  return electrolyte_phase_block_impl(temperature, amounts, volume,
                                      ion_mole_fraction_max, parameters, {});
}

// EqExclude: transport and buffer mechanics
inline NativeMixturePhaseBlock electrolyte_phase_block_active_parameters(
    double temperature, const std::vector<double> &amounts, double volume,
    double ion_mole_fraction_max, const NativeParameters &parameters,
    const std::vector<NativeActiveParameter> &active_parameters) {
  if (active_parameters.empty()) {
    throw std::invalid_argument(
        "reacting phase requires at least one active parameter");
  }
  return electrolyte_phase_block_impl(temperature, amounts, volume,
                                      ion_mole_fraction_max, parameters,
                                      active_parameters);
}

// EqID: phase_pressure
// EqID: association_mass_action
// EqExclude: CppAD-generated derivative machinery
// EqID: phase_packing_fraction
// EqExclude: CppAD-generated derivative machinery
inline NativePackingFractionBlock
packing_fraction_block(double temperature, const std::vector<double> &amounts,
                       double volume, double ion_mole_fraction_max,
                       const NativeParameters &parameters) {
  if (!supports_packing_geometry(parameters)) {
    throw std::logic_error(
        "packing fraction geometry requires compatible component parameters");
  }
  if (std::any_of(parameters.charge_number.begin(),
                  parameters.charge_number.end(),
                  [](int charge) { return charge != 0; })) {
    validate_electrolyte_phase_state(temperature, amounts, volume,
                                     ion_mole_fraction_max, parameters);
  } else {
    validate_mixture_phase_state(temperature, amounts, volume, parameters);
  }
  const double checked_value =
      phase_packing_fraction(temperature, amounts, volume, parameters);
  if (!std::isfinite(checked_value)) {
    throw std::invalid_argument("packing fraction is outside the model domain");
  }

  CppAdSerialGuard cppad_guard;
  using ADScalar = CppAD::AD<double>;
  std::vector<double> point(amounts.begin(), amounts.end());
  point.push_back(volume);
  std::vector<ADScalar> independent(point.begin(), point.end());
  CppAD::Independent(independent);
  const std::size_t volume_coordinate = amounts.size();
  const std::vector<ADScalar> active_amounts(
      independent.begin(), independent.begin() + volume_coordinate);
  const ADScalar active_value = phase_packing_fraction(
      temperature, active_amounts, independent[volume_coordinate], parameters);
  CppAD::ADFun<double> function(independent, {active_value});
  const double value = function.Forward(0, point)[0];
  const auto gradient = function.Jacobian(point);
  const auto hessian = function.Hessian(point, 0);
  if (!std::isfinite(value) ||
      !std::all_of(gradient.begin(), gradient.end(),
                   [](double item) { return std::isfinite(item); }) ||
      !std::all_of(hessian.begin(), hessian.end(),
                   [](double item) { return std::isfinite(item); })) {
    throw std::runtime_error(
        "packing fraction block contains nonfinite values");
  }
  return {value, gradient, hessian};
}

// EqID: phase_packing_fraction
// EqExclude: CppAD-generated derivative machinery
inline NativePackingFractionBlock packing_fraction_block_active_parameters(
    double temperature, const std::vector<double> &amounts, double volume,
    double ion_mole_fraction_max, const NativeParameters &parameters,
    const std::vector<NativeActiveParameter> &active_parameters) {
  if (active_parameters.empty()) {
    throw std::invalid_argument(
        "active packing fraction requires at least one parameter");
  }
  const NativeParameters value_parameters =
      active_component_parameter_values(parameters, active_parameters);
  validate_electrolyte_phase_state(temperature, amounts, volume,
                                   ion_mole_fraction_max, value_parameters);
  const double checked_value =
      phase_packing_fraction(temperature, amounts, volume, value_parameters);
  if (!std::isfinite(checked_value) ||
      checked_value >= close_packing_fraction()) {
    throw std::invalid_argument("packing fraction is outside the model domain");
  }

  CppAdSerialGuard cppad_guard;
  using ADScalar = CppAD::AD<double>;
  std::vector<double> point(amounts.begin(), amounts.end());
  point.push_back(volume);
  for (const auto &parameter : active_parameters) {
    point.push_back(parameter.value);
  }
  std::vector<ADScalar> independent(point.begin(), point.end());
  CppAD::Independent(independent);
  const std::size_t components = amounts.size();
  const std::size_t parameter_offset = components + 1;
  const std::vector<ADScalar> active_amounts(independent.begin(),
                                             independent.begin() + components);
  ADScalar total = 0.0;
  for (const auto &amount : active_amounts) {
    total += amount;
  }
  std::vector<ADScalar> composition;
  composition.reserve(components);
  for (const auto &amount : active_amounts) {
    composition.push_back(amount / total);
  }
  const std::vector<ADScalar> active_values(
      independent.begin() + parameter_offset, independent.end());
  const auto diameter = active_packing_diameters(
      temperature, parameters, active_parameters, active_values);
  const ADScalar active_value = packing_fraction_from_geometry(
      total / independent[components], composition, parameters, diameter);
  CppAD::ADFun<double> function(independent, {active_value});
  const double value = function.Forward(0, point)[0];
  const auto gradient = function.Jacobian(point);
  const auto hessian = function.Hessian(point, 0);
  if (!std::isfinite(value) ||
      !std::all_of(gradient.begin(), gradient.end(),
                   [](double item) { return std::isfinite(item); }) ||
      !std::all_of(hessian.begin(), hessian.end(),
                   [](double item) { return std::isfinite(item); })) {
    throw std::runtime_error(
        "active packing fraction block contains nonfinite values");
  }
  return {value, gradient, hessian};
}

// EqID: phase_packing_fraction
// EqExclude: CppAD-generated derivative machinery
inline NativeInversePackingGeometryBlock inverse_packing_geometry_block(
    double temperature, const std::vector<double> &amounts,
    double log_packing_fraction, const NativeParameters &parameters,
    const std::vector<NativeActiveParameter> &active_parameters) {
  if (!supports_packing_geometry(parameters)) {
    throw std::logic_error(
        "inverse packing geometry requires compatible component parameters");
  }
  const NativeParameters value_parameters =
      active_component_parameter_values(parameters, active_parameters);
  if (!supports_packing_geometry(value_parameters)) {
    throw std::invalid_argument(
        "active inverse packing parameters are outside the geometry domain");
  }
  if (!std::isfinite(temperature) || temperature <= 0.0 ||
      amounts.size() != parameters.segment_count.size() ||
      !std::all_of(amounts.begin(), amounts.end(),
                   [](double amount) {
                     return std::isfinite(amount) && amount > 0.0;
                   }) ||
      !std::isfinite(log_packing_fraction)) {
    throw std::invalid_argument(
        "inverse packing temperature, amounts, and log eta are invalid");
  }
  const double packing_fraction = std::exp(log_packing_fraction);
  if (!std::isfinite(packing_fraction) || packing_fraction <= 0.0 ||
      packing_fraction >= close_packing_fraction()) {
    throw std::invalid_argument(
        "inverse packing fraction is outside the model domain");
  }

  CppAdSerialGuard cppad_guard;
  using ADScalar = CppAD::AD<double>;
  std::vector<double> point(amounts.begin(), amounts.end());
  point.push_back(log_packing_fraction);
  for (const auto &parameter : active_parameters) {
    point.push_back(parameter.value);
  }
  std::vector<ADScalar> independent(point.begin(), point.end());
  CppAD::Independent(independent);
  const std::size_t components = amounts.size();
  const std::size_t parameter_offset = components + 1;
  const std::vector<ADScalar> active_amounts(independent.begin(),
                                             independent.begin() + components);
  ADScalar total = 0.0;
  for (const auto &amount : active_amounts) {
    total += amount;
  }
  std::vector<ADScalar> composition;
  composition.reserve(components);
  for (const auto &amount : active_amounts) {
    composition.push_back(amount / total);
  }
  const std::vector<ADScalar> active_values(
      independent.begin() + parameter_offset, independent.end());
  const auto diameter = active_packing_diameters(
      temperature, parameters, active_parameters, active_values);
  const ADScalar coefficient =
      packing_fraction_from_geometry(total, composition, parameters, diameter);
  const ADScalar volume = coefficient / scalar_exp(independent[components]);
  CppAD::ADFun<double> function(independent, {volume});
  const double value = function.Forward(0, point)[0];
  const auto gradient = function.Jacobian(point);
  const auto hessian = function.Hessian(point, 0);
  if (!std::isfinite(value) || value <= 0.0 ||
      !std::all_of(gradient.begin(), gradient.end(),
                   [](double item) { return std::isfinite(item); }) ||
      !std::all_of(hessian.begin(), hessian.end(),
                   [](double item) { return std::isfinite(item); })) {
    throw std::runtime_error(
        "inverse packing geometry contains nonfinite values");
  }
  return {value, gradient, hessian};
}

template <typename Scalar, typename Parameters>
Scalar pure_phase_helmholtz_over_rt_reference_amount(
    double temperature, const Scalar &amount, const Scalar &volume,
    const Parameters &parameters) {
  // EqID: phase_molar_density
  const Scalar density = amount / volume;
  const std::vector<Scalar> composition = {Scalar(1.0)};
  const auto sigma = segment_sigma(Scalar(temperature), parameters);
  const auto diameter =
      segment_diameter(Scalar(temperature), parameters, sigma);
  const auto contributions =
      residual_helmholtz(Scalar(temperature), density, composition, parameters,
                         sigma, diameter, {}, Scalar(0.0), Scalar(0.0));
  const Scalar residual =
      contributions.hard_chain + contributions.dispersion + contributions.polar;
  // EqID: phase_dimensionless_helmholtz
  return amount * (scalar_log(density) - Scalar(1.0) + residual);
}

// EqExclude: CppAD-generated derivative machinery
inline double
cubic_directional_derivative(CppAD::ADFun<double> &function,
                             const std::array<double, 2> &point,
                             const std::array<double, 2> &direction) {
  const std::vector<double> point_vector(point.begin(), point.end());
  const std::vector<double> direction_vector(direction.begin(),
                                             direction.end());
  const std::vector<double> zero(2, 0.0);
  function.Forward(0, point_vector);
  function.Forward(1, direction_vector);
  function.Forward(2, zero);
  return 6.0 * function.Forward(3, zero)[0];
}

// EqID: phase_pressure
// EqExclude: CppAD-generated derivative machinery
inline NativePurePhaseBlock
pure_phase_block(double temperature, double amount, double volume,
                 const NativeParameters &parameters) {
  if (!supports_pure_phase_block(parameters)) {
    throw std::logic_error(
        "native transport v1 requires one pure nonassociating neutral model");
  }
  if (!std::isfinite(temperature) || !std::isfinite(amount) ||
      !std::isfinite(volume) || temperature <= 0.0 || amount <= 0.0 ||
      volume <= 0.0) {
    throw std::invalid_argument(
        "temperature, amount, and volume must be positive and finite");
  }

  pure_phase_helmholtz_over_rt_reference_amount(temperature, amount, volume,
                                                parameters);

  CppAdSerialGuard cppad_guard;
  using ADScalar = CppAD::AD<double>;
  std::vector<ADScalar> independent = {amount, volume};
  CppAD::Independent(independent);
  const ADScalar phase = pure_phase_helmholtz_over_rt_reference_amount(
      temperature, independent[0], independent[1], parameters);
  CppAD::ADFun<double> function(independent, {phase});

  const std::array<double, 2> point = {amount, volume};
  const std::vector<double> point_vector(point.begin(), point.end());
  const auto jacobian = function.Jacobian(point_vector);
  const auto hessian_vector = function.Hessian(point_vector, 0);
  const std::array<double, 2> gradient = {jacobian[0], jacobian[1]};
  const std::array<double, 4> hessian = {
      hessian_vector[0],
      hessian_vector[1],
      hessian_vector[2],
      hessian_vector[3],
  };

  const double t_nnn =
      cubic_directional_derivative(function, point, {1.0, 0.0});
  const double t_vvv =
      cubic_directional_derivative(function, point, {0.0, 1.0});
  const double plus = cubic_directional_derivative(function, point, {1.0, 1.0});
  const double minus =
      cubic_directional_derivative(function, point, {1.0, -1.0});
  const double sum_mixed = (plus - t_nnn - t_vvv) / 3.0;
  const double mixed_difference = (minus - t_nnn + t_vvv) / 3.0;
  const double t_nnv = 0.5 * (sum_mixed - mixed_difference);
  const double t_nvv = 0.5 * (sum_mixed + mixed_difference);
  const std::array<double, 8> third = {
      t_nnn, t_nnv, t_nnv, t_nvv, t_nnv, t_nvv, t_nvv, t_vvv,
  };

  const double value = function.Forward(0, point_vector)[0];
  if (!std::isfinite(value) ||
      !std::all_of(gradient.begin(), gradient.end(),
                   [](double item) { return std::isfinite(item); }) ||
      !std::all_of(hessian.begin(), hessian.end(),
                   [](double item) { return std::isfinite(item); }) ||
      !std::all_of(third.begin(), third.end(),
                   [](double item) { return std::isfinite(item); })) {
    throw std::runtime_error("phase block contains nonfinite values");
  }

  return {
      value,
      gradient,
      hessian,
      third,
      -GAS_CONSTANT * temperature * gradient[1],
      gradient[0],
  };
}

// EqID: phase_pressure
// EqID: association_mass_action
// EqExclude: CppAD-generated derivative machinery
inline NativeMixturePhaseBlock fixed_topology_parameterized_mixture_phase_block(
    double temperature, const std::vector<double> &amounts, double volume,
    const std::vector<NativeFixedTopologyParameter> &active_parameters,
    const NativeParameters &base_parameters) {
  if (base_parameters.segment_count.empty()) {
    throw std::logic_error("active parameters require at least one component");
  }
  const bool electrolyte = std::any_of(base_parameters.charge_number.begin(),
                                       base_parameters.charge_number.end(),
                                       [](int charge) { return charge != 0; });
  if (electrolyte) {
    validate_electrolyte_phase_state(temperature, amounts, volume,
                                     std::numeric_limits<double>::infinity(),
                                     base_parameters);
  } else {
    validate_mixture_phase_state(temperature, amounts, volume, base_parameters);
  }
  const std::size_t components = base_parameters.segment_count.size();
  const std::size_t sites = base_parameters.site_component.size();
  NativeParameters value_parameters = base_parameters;
  std::vector<std::pair<NativeFixedTopologyParameterKind, std::size_t>>
      identities;
  identities.reserve(active_parameters.size());
  for (const auto &parameter : active_parameters) {
    const bool signed_correlation =
        parameter.kind == NativeFixedTopologyParameterKind::
                              SegmentDiameterCorrelationAmplitude ||
        parameter.kind == NativeFixedTopologyParameterKind::
                              SegmentDiameterCorrelationExponent;
    if (!std::isfinite(parameter.value) ||
        (!signed_correlation && parameter.value <= 0.0)) {
      throw std::invalid_argument(
          "fixed-topology active values must be positive and finite");
    }
    std::size_t identity = 0;
    switch (parameter.kind) {
    case NativeFixedTopologyParameterKind::SegmentCount:
      if (parameter.component >= components) {
        throw std::invalid_argument(
            "component parameter index is outside the fixed topology");
      }
      identity = parameter.component;
      value_parameters.segment_count[parameter.component] = parameter.value;
      break;
    case NativeFixedTopologyParameterKind::SegmentDiameter:
      if (parameter.component >= components) {
        throw std::invalid_argument(
            "component parameter index is outside the fixed topology");
      }
      identity = parameter.component;
      value_parameters.sigma_constant[parameter.component] = parameter.value;
      break;
    case NativeFixedTopologyParameterKind::SegmentDiameterCorrelationConstant:
      if (parameter.component >= components) {
        throw std::invalid_argument(
            "correlation component index is outside the model");
      }
      identity = parameter.component;
      value_parameters.sigma_constant[parameter.component] = parameter.value;
      break;
    case NativeFixedTopologyParameterKind::SegmentDiameterCorrelationAmplitude:
      if (parameter.correlation_term >=
          value_parameters.sigma_term_amplitude.size()) {
        throw std::invalid_argument(
            "correlation term index is outside the model");
      }
      identity = parameter.correlation_term;
      value_parameters.sigma_term_amplitude[parameter.correlation_term] =
          parameter.value;
      break;
    case NativeFixedTopologyParameterKind::SegmentDiameterCorrelationExponent:
      if (parameter.correlation_term >=
          value_parameters.sigma_term_exponent.size()) {
        throw std::invalid_argument(
            "correlation term index is outside the model");
      }
      identity = parameter.correlation_term;
      value_parameters.sigma_term_exponent[parameter.correlation_term] =
          parameter.value;
      break;
    case NativeFixedTopologyParameterKind::DispersionEnergyOverK:
      if (parameter.component >= components) {
        throw std::invalid_argument(
            "component parameter index is outside the fixed topology");
      }
      identity = parameter.component;
      value_parameters.dispersion_energy_over_k[parameter.component] =
          parameter.value;
      break;
    case NativeFixedTopologyParameterKind::AssociationEnergyOverK:
    case NativeFixedTopologyParameterKind::AssociationVolume: {
      if (parameter.site_a > parameter.site_b || parameter.site_b >= sites) {
        throw std::invalid_argument(
            "association endpoints must be canonical model site indices");
      }
      const std::size_t forward = parameter.site_a * sites + parameter.site_b;
      const std::size_t reverse = parameter.site_b * sites + parameter.site_a;
      if (base_parameters.association_active[forward] == 0 ||
          base_parameters.association_active[reverse] == 0) {
        throw std::invalid_argument(
            "association endpoint pair is not active in the fixed topology");
      }
      identity = forward;
      auto &values =
          parameter.kind ==
                  NativeFixedTopologyParameterKind::AssociationEnergyOverK
              ? value_parameters.association_energy_over_k
              : value_parameters.association_volume;
      values[forward] = parameter.value;
      values[reverse] = parameter.value;
      break;
    }
    }
    const auto key = std::pair{parameter.kind, identity};
    if (std::find(identities.begin(), identities.end(), key) !=
        identities.end()) {
      throw std::invalid_argument(
          "fixed-topology active parameter is duplicated");
    }
    identities.push_back(key);
  }

  const double total_amount =
      std::accumulate(amounts.begin(), amounts.end(), 0.0);
  const double density = total_amount / volume;
  std::vector<double> composition;
  composition.reserve(components);
  for (const double amount : amounts) {
    composition.push_back(amount / total_amount);
  }
  const auto value_sigma = segment_sigma(temperature, value_parameters);
  const auto value_diameter =
      segment_diameter(temperature, value_parameters, value_sigma);
  const auto site_fractions =
      solve_site_fractions(temperature, density, composition, value_parameters,
                           value_sigma, value_diameter);

  CppAdSerialGuard cppad_guard;
  using ADScalar = CppAD::AD<double>;
  const std::size_t volume_coordinate = components;
  const std::size_t parameter_coordinate = volume_coordinate + 1;
  const std::size_t coordinates =
      parameter_coordinate + active_parameters.size();
  const std::size_t width = coordinates + sites;
  std::vector<ADScalar> independent(width);
  std::copy(amounts.begin(), amounts.end(), independent.begin());
  independent[volume_coordinate] = volume;
  for (std::size_t index = 0; index < active_parameters.size(); ++index) {
    independent[parameter_coordinate + index] = active_parameters[index].value;
  }
  for (std::size_t site = 0; site < sites; ++site) {
    independent[coordinates + site] = site_fractions[site];
  }
  CppAD::Independent(independent);

  auto parameters = phase_block_parameters<ADScalar>(base_parameters);
  std::vector<ADScalar> association_energies(
      base_parameters.association_energy_over_k.begin(),
      base_parameters.association_energy_over_k.end());
  std::vector<ADScalar> association_volumes(
      base_parameters.association_volume.begin(),
      base_parameters.association_volume.end());
  for (std::size_t index = 0; index < active_parameters.size(); ++index) {
    const auto &parameter = active_parameters[index];
    const ADScalar &active = independent[parameter_coordinate + index];
    switch (parameter.kind) {
    case NativeFixedTopologyParameterKind::SegmentCount:
      parameters.segment_count[parameter.component] = active;
      break;
    case NativeFixedTopologyParameterKind::SegmentDiameter:
      parameters.sigma_constant[parameter.component] = active;
      break;
    case NativeFixedTopologyParameterKind::SegmentDiameterCorrelationConstant:
      parameters.sigma_constant[parameter.component] = active;
      break;
    case NativeFixedTopologyParameterKind::SegmentDiameterCorrelationAmplitude:
      parameters.sigma_term_amplitude[parameter.correlation_term] = active;
      break;
    case NativeFixedTopologyParameterKind::SegmentDiameterCorrelationExponent:
      parameters.sigma_term_exponent[parameter.correlation_term] = active;
      break;
    case NativeFixedTopologyParameterKind::DispersionEnergyOverK:
      parameters.dispersion_energy_over_k[parameter.component] = active;
      break;
    case NativeFixedTopologyParameterKind::AssociationEnergyOverK:
    case NativeFixedTopologyParameterKind::AssociationVolume: {
      const std::size_t forward = parameter.site_a * sites + parameter.site_b;
      const std::size_t reverse = parameter.site_b * sites + parameter.site_a;
      auto &values =
          parameter.kind ==
                  NativeFixedTopologyParameterKind::AssociationEnergyOverK
              ? association_energies
              : association_volumes;
      values[forward] = active;
      values[reverse] = active;
      break;
    }
    }
  }

  std::vector<ADScalar> active_amounts(independent.begin(),
                                       independent.begin() + volume_coordinate);
  ADScalar active_total = ADScalar(0.0);
  for (const ADScalar &amount : active_amounts) {
    active_total += amount;
  }
  const ADScalar active_density = active_total / independent[volume_coordinate];
  std::vector<ADScalar> active_composition;
  active_composition.reserve(components);
  for (const ADScalar &amount : active_amounts) {
    active_composition.push_back(amount / active_total);
  }
  const std::vector<ADScalar> active_site_fractions(
      independent.begin() + static_cast<std::ptrdiff_t>(coordinates),
      independent.end());
  const auto sigma = segment_sigma(ADScalar(temperature), parameters);
  const auto diameter =
      segment_diameter(ADScalar(temperature), parameters, sigma);
  const auto reference_sigma = segment_sigma(temperature, value_parameters);
  const auto reference_diameter =
      segment_diameter(temperature, value_parameters, reference_sigma);
  const std::vector<ADScalar> reference_diameter_ad(reference_diameter.begin(),
                                                    reference_diameter.end());
  const auto electrolyte_contributions = electrolyte_contribution(
      ADScalar(temperature), active_density, temperature, density,
      active_composition, composition, base_parameters, diameter,
      reference_diameter_ad,
      static_cast<const std::vector<ADScalar> *>(nullptr), &reference_diameter);
  const auto born = born_contribution(
      ADScalar(temperature), active_composition, base_parameters,
      electrolyte_contributions.bulk_relative_permittivity, sigma);
  const auto contributions = residual_helmholtz(
      ADScalar(temperature), active_density, active_composition, parameters,
      sigma, diameter, active_site_fractions,
      electrolyte_contributions.debye_huckel, born);
  ADScalar ideal = ADScalar(0.0);
  for (const ADScalar &amount : active_amounts) {
    ideal += amount * (scalar_log(amount / independent[volume_coordinate]) -
                       ADScalar(1.0));
  }
  const ADScalar phase =
      ideal +
      active_total * (contributions.hard_chain + contributions.dispersion +
                      contributions.association + contributions.polar +
                      contributions.debye_huckel + contributions.born_ssm_ds);
  const auto delta = association_delta_values(
      ADScalar(temperature), active_density, active_composition, parameters,
      sigma, diameter, association_energies, association_volumes);
  const auto mass_action =
      association_residuals(active_density, active_composition, parameters,
                            active_site_fractions, delta);
  std::vector<ADScalar> dependent{phase};
  dependent.insert(dependent.end(), mass_action.begin(), mass_action.end());
  CppAD::ADFun<double> function(independent, dependent);

  std::vector<double> point(amounts.begin(), amounts.end());
  point.push_back(volume);
  for (const auto &parameter : active_parameters) {
    point.push_back(parameter.value);
  }
  point.insert(point.end(), site_fractions.begin(), site_fractions.end());
  const auto jacobian = function.Jacobian(point);
  std::vector<double> association_jacobian(sites * sites);
  for (std::size_t row = 0; row < sites; ++row) {
    for (std::size_t column = 0; column < sites; ++column) {
      association_jacobian[row * sites + column] =
          jacobian[(1 + row) * width + coordinates + column];
    }
  }
  std::vector<std::vector<double>> site_sensitivity(
      coordinates, std::vector<double>(sites, 0.0));
  for (std::size_t coordinate = 0; coordinate < coordinates; ++coordinate) {
    std::vector<double> right_hand_side(sites);
    for (std::size_t row = 0; row < sites; ++row) {
      right_hand_side[row] = -jacobian[(1 + row) * width + coordinate];
    }
    site_sensitivity[coordinate] =
        solve_linear(association_jacobian, std::move(right_hand_side));
  }
  std::vector<std::vector<double>> directions(coordinates,
                                              std::vector<double>(width, 0.0));
  for (std::size_t coordinate = 0; coordinate < coordinates; ++coordinate) {
    directions[coordinate][coordinate] = 1.0;
    for (std::size_t site = 0; site < sites; ++site) {
      directions[coordinate][coordinates + site] =
          site_sensitivity[coordinate][site];
    }
  }
  const auto phase_hessian = function.Hessian(point, 0);
  std::vector<std::vector<double>> association_hessians(sites);
  for (std::size_t site = 0; site < sites; ++site) {
    association_hessians[site] = function.Hessian(point, 1 + site);
  }
  std::vector<double> gradient(coordinates);
  std::vector<double> hessian(coordinates * coordinates);
  for (std::size_t coordinate = 0; coordinate < coordinates; ++coordinate) {
    gradient[coordinate] = jacobian[coordinate];
    for (std::size_t site = 0; site < sites; ++site) {
      gradient[coordinate] +=
          jacobian[coordinates + site] * site_sensitivity[coordinate][site];
    }
  }
  for (std::size_t left = 0; left < coordinates; ++left) {
    for (std::size_t right = 0; right < coordinates; ++right) {
      std::vector<double> right_hand_side(sites);
      for (std::size_t site = 0; site < sites; ++site) {
        right_hand_side[site] = -directional_second(
            association_hessians[site], directions[left], directions[right]);
      }
      const auto site_second =
          solve_linear(association_jacobian, std::move(right_hand_side));
      double second = directional_second(phase_hessian, directions[left],
                                         directions[right]);
      for (std::size_t site = 0; site < sites; ++site) {
        second += jacobian[coordinates + site] * site_second[site];
      }
      hessian[left * coordinates + right] = second;
    }
  }
  const double value = function.Forward(0, point)[0];
  if (!std::isfinite(value) ||
      !std::all_of(gradient.begin(), gradient.end(),
                   [](double item) { return std::isfinite(item); }) ||
      !std::all_of(hessian.begin(), hessian.end(),
                   [](double item) { return std::isfinite(item); })) {
    throw std::runtime_error(
        "fixed-topology phase block contains nonfinite values");
  }
  return {
      value,
      gradient,
      hessian,
      -GAS_CONSTANT * temperature * gradient[volume_coordinate],
  };
}

// EqExclude: transport and buffer mechanics
inline NativeFixedTopologyPhaseBlock fixed_topology_parameterized_pure_phase_block(
    double temperature, double amount, double volume,
    const std::vector<NativeFixedTopologyParameter> &active_parameters,
    const NativeParameters &base_parameters) {
  if (base_parameters.segment_count.size() != 1) {
    throw std::logic_error("active pure parameters require one component");
  }
  const auto block = fixed_topology_parameterized_mixture_phase_block(
      temperature, {amount}, volume, active_parameters, base_parameters);
  return {
      block.helmholtz_over_rt_reference_amount,
      block.gradient,
      block.hessian,
      block.pressure_pa,
      block.gradient[0],
  };
}

} // namespace detail
} // namespace epcsaft
