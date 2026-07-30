#pragma once

#include "association.hpp"
#include "derivatives.hpp"
#include "helmholtz.hpp"
#include "model.hpp"

#include <cppad/cppad.hpp>

#include <algorithm>
#include <array>
#include <cmath>
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

struct NativeParameterizedPurePhaseBlock final {
    double helmholtz_over_rt_reference_amount;
    std::array<double, 5> gradient;
    std::array<double, 25> hessian;
    double pressure_pa;
    double chemical_potential_over_rt;
};

struct NativeAssociatingParameterizedPurePhaseBlock final {
    double helmholtz_over_rt_reference_amount;
    std::array<double, 7> gradient;
    std::array<double, 49> hessian;
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

namespace detail {

// EqExclude: input and domain validation
inline bool supports_pure_phase_block(const NativeParameters& parameters) {
    return parameters.segment_count.size() == 1
        && parameters.site_component.empty()
        && parameters.charge_number.size() == 1
        && parameters.charge_number[0] == 0;
}

// EqExclude: input and domain validation
inline bool supports_associating_pure_phase_block(
    const NativeParameters& parameters
) {
    if (parameters.segment_count.size() != 1
        || parameters.site_component.size() != 2
        || parameters.charge_number.size() != 1
        || parameters.charge_number[0] != 0
        || !std::all_of(
            parameters.site_component.begin(),
            parameters.site_component.end(),
            [](int component) { return component == 0; }
        )) {
        return false;
    }
    const auto& active = parameters.association_active;
    return active.size() == 4
        && parameters.site_multiplicity
            == std::vector<int>({1, 1})
        && active[0] == 0 && active[1] != 0
        && active[2] != 0 && active[3] == 0
        && parameters.association_energy_over_k[1]
            == parameters.association_energy_over_k[2]
        && parameters.association_volume[1]
            == parameters.association_volume[2];
}

// EqExclude: input and domain validation
inline bool supports_mixture_phase_block(const NativeParameters& parameters) {
    return parameters.segment_count.size() >= 2
        && parameters.site_component.empty()
        && std::all_of(
            parameters.charge_number.begin(),
            parameters.charge_number.end(),
            [](int charge) { return charge == 0; }
        );
}

// EqExclude: input and domain validation
inline bool supports_binary_interaction_phase_block(
    const NativeParameters& parameters
) {
    return parameters.segment_count.size() == 2
        && supports_mixture_phase_block(parameters);
}

// EqExclude: input and domain validation
inline bool supports_electrolyte_phase_block(const NativeParameters& parameters) {
    return parameters.segment_count.size() >= 3
        && !parameters.site_component.empty()
        && std::any_of(
            parameters.charge_number.begin(),
            parameters.charge_number.end(),
            [](int charge) { return charge != 0; }
        );
}

template <typename Scalar>
struct PhaseBlockParameters final {
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
};

template <typename Scalar>
// EqExclude: scalar utility dispatch
std::vector<Scalar> phase_block_scalars(const std::vector<double>& values) {
    return std::vector<Scalar>(values.begin(), values.end());
}

template <typename Scalar>
// EqID: active_binary_interaction_mapping
PhaseBlockParameters<Scalar> binary_interaction_parameters(
    const NativeParameters& parameters,
    const Scalar& k_ij
) {
    return {
        phase_block_scalars<Scalar>(parameters.segment_count),
        phase_block_scalars<Scalar>(parameters.sigma_constant),
        parameters.sigma_term_component,
        phase_block_scalars<Scalar>(parameters.sigma_term_amplitude),
        phase_block_scalars<Scalar>(parameters.sigma_term_exponent),
        phase_block_scalars<Scalar>(parameters.dispersion_energy_over_k),
        {Scalar(parameters.k_ij[0]), k_ij, k_ij, Scalar(parameters.k_ij[3])},
        phase_block_scalars<Scalar>(parameters.l_ij),
        parameters.site_component,
        parameters.site_multiplicity,
        parameters.association_active,
        parameters.charge_number,
        parameters.electrolyte_model,
    };
}

template <typename Scalar>
// EqID: active_binary_size_interaction_mapping
PhaseBlockParameters<Scalar> binary_size_interaction_parameters(
    const NativeParameters& parameters,
    const Scalar& l_ij
) {
    return {
        phase_block_scalars<Scalar>(parameters.segment_count),
        phase_block_scalars<Scalar>(parameters.sigma_constant),
        parameters.sigma_term_component,
        phase_block_scalars<Scalar>(parameters.sigma_term_amplitude),
        phase_block_scalars<Scalar>(parameters.sigma_term_exponent),
        phase_block_scalars<Scalar>(parameters.dispersion_energy_over_k),
        phase_block_scalars<Scalar>(parameters.k_ij),
        {Scalar(parameters.l_ij[0]), l_ij, l_ij, Scalar(parameters.l_ij[3])},
        parameters.site_component,
        parameters.site_multiplicity,
        parameters.association_active,
        parameters.charge_number,
        parameters.electrolyte_model,
    };
}

template <typename Scalar, typename Parameters>
Scalar mixture_phase_helmholtz_over_rt_reference_amount(
    double temperature,
    const std::vector<Scalar>& amounts,
    const Scalar& volume,
    const Parameters& parameters
) {
    // EqID: phase_total_amount
    Scalar total = Scalar(0.0);
    for (const Scalar& amount : amounts) {
        total += amount;
    }
    // EqID: phase_molar_density
    const Scalar density = total / volume;
    std::vector<Scalar> composition;
    composition.reserve(amounts.size());
    for (const Scalar& amount : amounts) {
        // EqID: phase_mole_fraction
        composition.push_back(amount / total);
    }
    const auto sigma = segment_sigma(Scalar(temperature), parameters);
    const auto diameter = segment_diameter(Scalar(temperature), parameters, sigma);
    const auto contributions = residual_helmholtz(
        Scalar(temperature),
        density,
        composition,
        parameters,
        sigma,
        diameter,
        {},
        Scalar(0.0),
        Scalar(0.0)
    );
    Scalar ideal = Scalar(0.0);
    for (const Scalar& amount : amounts) {
        ideal += amount * (scalar_log(amount / volume) - Scalar(1.0));
    }
    // EqID: phase_dimensionless_helmholtz
    return ideal + total * (contributions.hard_chain + contributions.dispersion);
}

// EqExclude: input and domain validation
inline void validate_mixture_phase_state(
    double temperature,
    const std::vector<double>& amounts,
    double volume,
    const NativeParameters& parameters
) {
    if (!std::isfinite(temperature) || temperature <= 0.0
        || !std::isfinite(volume) || volume <= 0.0
        || amounts.size() != parameters.segment_count.size()
        || !std::all_of(amounts.begin(), amounts.end(), [](double amount) {
            return std::isfinite(amount) && amount > 0.0;
        })) {
        throw std::invalid_argument(
            "temperature, amounts, and volume must be positive and finite"
        );
    }
}

// EqExclude: input and domain validation
inline void validate_electrolyte_phase_state(
    double temperature,
    const std::vector<double>& amounts,
    double volume,
    double ion_mole_fraction_max,
    const NativeParameters& parameters
) {
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
        throw std::invalid_argument("electrolyte phase amounts must be electroneutral");
    }
    if (std::isfinite(ion_mole_fraction_max)
        && ionic_amount / total_amount > ion_mole_fraction_max + 1.0e-12) {
        throw std::invalid_argument("ion mole fraction exceeds the parameter source domain");
    }
}

template <typename PhaseFunction>
// EqID: phase_pressure
// EqExclude: CppAD-generated derivative machinery
NativeMixturePhaseBlock differentiated_mixture_phase_block(
    double temperature,
    const std::vector<double>& point,
    std::size_t volume_coordinate,
    PhaseFunction phase_function
) {
    CppAdSerialGuard cppad_guard;
    using ADScalar = CppAD::AD<double>;
    std::vector<ADScalar> independent(point.begin(), point.end());
    CppAD::Independent(independent);
    const ADScalar phase = phase_function(independent);
    CppAD::ADFun<double> function(independent, {phase});

    const auto gradient = function.Jacobian(point);
    const auto hessian = function.Hessian(point, 0);
    const double value = function.Forward(0, point)[0];
    if (!std::isfinite(value)
        || !std::all_of(gradient.begin(), gradient.end(), [](double item) {
            return std::isfinite(item);
        })
        || !std::all_of(hessian.begin(), hessian.end(), [](double item) {
            return std::isfinite(item);
        })) {
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
inline NativeMixturePhaseBlock mixture_phase_block(
    double temperature,
    const std::vector<double>& amounts,
    double volume,
    const NativeParameters& parameters
) {
    if (!supports_mixture_phase_block(parameters)) {
        throw std::logic_error(
            "mixture phase block requires at least two neutral nonassociating components"
        );
    }
    validate_mixture_phase_state(temperature, amounts, volume, parameters);

    mixture_phase_helmholtz_over_rt_reference_amount(
        temperature, amounts, volume, parameters
    );

    std::vector<double> point(amounts.begin(), amounts.end());
    point.push_back(volume);
    const std::size_t volume_coordinate = amounts.size();
    return differentiated_mixture_phase_block(
        temperature,
        point,
        volume_coordinate,
        [&](const auto& independent) {
            using ADScalar = CppAD::AD<double>;
            const std::vector<ADScalar> active_amounts(
                independent.begin(), independent.begin() + volume_coordinate
            );
            return mixture_phase_helmholtz_over_rt_reference_amount(
                temperature,
                active_amounts,
                independent[volume_coordinate],
                parameters
            );
        }
    );
}

// EqExclude: transport and buffer mechanics
inline NativeMixturePhaseBlock mixture_phase_block_kij(
    double temperature,
    const std::vector<double>& amounts,
    double volume,
    double k_ij,
    const NativeParameters& parameters
) {
    if (!supports_binary_interaction_phase_block(parameters)) {
        throw std::logic_error(
            "active k_ij phase block requires two neutral nonassociating components"
        );
    }
    if (!std::isfinite(k_ij)) {
        throw std::invalid_argument("k_ij must be finite");
    }
    validate_mixture_phase_state(temperature, amounts, volume, parameters);
    const auto value_parameters = binary_interaction_parameters(parameters, k_ij);
    mixture_phase_helmholtz_over_rt_reference_amount(
        temperature, amounts, volume, value_parameters
    );

    std::vector<double> point(amounts.begin(), amounts.end());
    point.push_back(volume);
    point.push_back(k_ij);
    const std::size_t volume_coordinate = amounts.size();
    return differentiated_mixture_phase_block(
        temperature,
        point,
        volume_coordinate,
        [&](const auto& independent) {
            using ADScalar = CppAD::AD<double>;
            const std::vector<ADScalar> active_amounts(
                independent.begin(), independent.begin() + volume_coordinate
            );
            const auto active_parameters = binary_interaction_parameters(
                parameters, independent[volume_coordinate + 1]
            );
            return mixture_phase_helmholtz_over_rt_reference_amount(
                temperature,
                active_amounts,
                independent[volume_coordinate],
                active_parameters
            );
        }
    );
}

// EqExclude: transport and buffer mechanics
inline NativeMixturePhaseBlock mixture_phase_block_lij(
    double temperature,
    const std::vector<double>& amounts,
    double volume,
    double l_ij,
    const NativeParameters& parameters
) {
    if (!supports_binary_interaction_phase_block(parameters)) {
        throw std::logic_error(
            "active l_ij phase block requires two neutral nonassociating components"
        );
    }
    if (!std::isfinite(l_ij)) {
        throw std::invalid_argument("l_ij must be finite");
    }
    validate_mixture_phase_state(temperature, amounts, volume, parameters);
    const auto value_parameters =
        binary_size_interaction_parameters(parameters, l_ij);
    mixture_phase_helmholtz_over_rt_reference_amount(
        temperature, amounts, volume, value_parameters
    );

    std::vector<double> point(amounts.begin(), amounts.end());
    point.push_back(volume);
    point.push_back(l_ij);
    const std::size_t volume_coordinate = amounts.size();
    return differentiated_mixture_phase_block(
        temperature,
        point,
        volume_coordinate,
        [&](const auto& independent) {
            using ADScalar = CppAD::AD<double>;
            const std::vector<ADScalar> active_amounts(
                independent.begin(), independent.begin() + volume_coordinate
            );
            const auto active_parameters = binary_size_interaction_parameters(
                parameters, independent[volume_coordinate + 1]
            );
            return mixture_phase_helmholtz_over_rt_reference_amount(
                temperature,
                active_amounts,
                independent[volume_coordinate],
                active_parameters
            );
        }
    );
}

// EqID: electrolyte_phase_coordinate_jacobian
// EqID: electrolyte_phase_hessian_chain_rule
inline NativeMixturePhaseBlock electrolyte_phase_block(
    double temperature,
    const std::vector<double>& amounts,
    double volume,
    double ion_mole_fraction_max,
    const NativeParameters& parameters
) {
    if (!supports_electrolyte_phase_block(parameters)) {
        throw std::logic_error(
            "electrolyte phase block requires an ionic associating model"
        );
    }
    validate_electrolyte_phase_state(
        temperature, amounts, volume, ion_mole_fraction_max, parameters
    );

    const std::size_t components = amounts.size();
    const std::size_t coordinates = components + 1;
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
    const auto sigma = segment_sigma(temperature, parameters);
    const auto diameter = segment_diameter(temperature, parameters, sigma);
    const auto site_fractions = solve_site_fractions(
        temperature, density, composition, parameters, sigma, diameter
    );
    const auto derivatives = thermodynamic_derivatives(
        temperature, density, composition, parameters, site_fractions
    );

    std::vector<double> state_gradient(1 + components);
    state_gradient[0] = derivatives.density;
    std::copy(
        derivatives.composition.begin(),
        derivatives.composition.end(),
        state_gradient.begin() + 1
    );
    std::vector<double> state_jacobian((1 + components) * coordinates, 0.0);
    for (std::size_t amount = 0; amount < components; ++amount) {
        state_jacobian[amount] = 1.0 / volume;
        for (std::size_t component = 0; component < components; ++component) {
            state_jacobian[(1 + component) * coordinates + amount] =
                ((component == amount ? 1.0 : 0.0) - composition[component])
                / total_amount;
        }
    }
    state_jacobian[components] = -total_amount / (volume * volume);

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

    auto state_second = [&](std::size_t state, std::size_t left, std::size_t right) {
        if (state == 0) {
            if (left == components && right == components) {
                return 2.0 * total_amount / (volume * volume * volume);
            }
            if (left == components || right == components) {
                return -1.0 / (volume * volume);
            }
            return 0.0;
        }
        if (left == components || right == components) {
            return 0.0;
        }
        const std::size_t component = state - 1;
        return (
            2.0 * composition[component]
            - (component == left ? 1.0 : 0.0)
            - (component == right ? 1.0 : 0.0)
        ) / (total_amount * total_amount);
    };

    for (std::size_t left = 0; left < coordinates; ++left) {
        const double total_left = left < components ? 1.0 : 0.0;
        double residual_gradient = total_left * derivatives.residual_helmholtz;
        for (std::size_t state = 0; state < 1 + components; ++state) {
            residual_gradient += total_amount * state_gradient[state]
                * state_jacobian[state * coordinates + left];
        }
        gradient[left] += residual_gradient;

        for (std::size_t right = 0; right < coordinates; ++right) {
            const double total_right = right < components ? 1.0 : 0.0;
            double residual_hessian = 0.0;
            for (std::size_t state = 0; state < 1 + components; ++state) {
                const double left_derivative =
                    state_jacobian[state * coordinates + left];
                const double right_derivative =
                    state_jacobian[state * coordinates + right];
                residual_hessian += state_gradient[state]
                    * (total_left * right_derivative
                       + total_right * left_derivative
                       + total_amount * state_second(state, left, right));
                for (std::size_t other = 0; other < 1 + components; ++other) {
                    residual_hessian += total_amount
                        * derivatives.hessian[(1 + state) * (2 + components) + 1 + other]
                        * left_derivative
                        * state_jacobian[other * coordinates + right];
                }
            }
            hessian[left * coordinates + right] += residual_hessian;
        }
    }

    // EqID: phase_dimensionless_helmholtz
    const double value = ideal + total_amount * derivatives.residual_helmholtz;
    if (!std::isfinite(value)
        || !std::all_of(gradient.begin(), gradient.end(), [](double item) {
            return std::isfinite(item);
        })
        || !std::all_of(hessian.begin(), hessian.end(), [](double item) {
            return std::isfinite(item);
        })) {
        throw std::runtime_error("electrolyte phase block contains nonfinite values");
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

// EqID: phase_pressure
// EqID: association_mass_action
// EqExclude: CppAD-generated derivative machinery
inline NativeMixturePhaseBlock associating_pure_phase_parameter_block(
    double temperature,
    double amount,
    double volume,
    double parameter_value,
    bool active_energy,
    const NativeParameters& parameters
) {
    if (!supports_associating_pure_phase_block(parameters)) {
        throw std::logic_error(
            "active association parameters require one neutral associating component"
        );
    }
    if (!std::isfinite(temperature) || !std::isfinite(amount)
        || !std::isfinite(volume) || !std::isfinite(parameter_value)
        || temperature <= 0.0 || amount <= 0.0 || volume <= 0.0
        || parameter_value <= 0.0) {
        throw std::invalid_argument(
            "association phase inputs must be positive and finite"
        );
    }

    NativeParameters value_parameters = parameters;
    auto& active_values = active_energy
        ? value_parameters.association_energy_over_k
        : value_parameters.association_volume;
    for (std::size_t pair = 0; pair < value_parameters.association_active.size();
         ++pair) {
        if (value_parameters.association_active[pair] != 0) {
            active_values[pair] = parameter_value;
        }
    }

    const double density = amount / volume;
    const std::vector<double> composition{1.0};
    const auto sigma = segment_sigma(temperature, value_parameters);
    const auto diameter = segment_diameter(temperature, value_parameters, sigma);
    const auto site_fractions = solve_site_fractions(
        temperature,
        density,
        composition,
        value_parameters,
        sigma,
        diameter
    );
    const auto derivatives =
        thermodynamic_derivatives_active_pure_association_parameter(
            temperature,
            density,
            composition,
            parameters,
            site_fractions,
            parameter_value,
            active_energy
        );

    const std::size_t state_coordinates = 4;
    const std::size_t parameter_coordinate = 3;
    const double a = derivatives.residual_helmholtz;
    const double a_rho = derivatives.density;
    const double a_parameter = *derivatives.active_association_parameter;
    const double a_rho_rho =
        derivatives.hessian[1 * state_coordinates + 1];
    const double a_rho_parameter =
        derivatives.hessian[1 * state_coordinates + parameter_coordinate];
    const double a_parameter_parameter =
        derivatives.hessian[
            parameter_coordinate * state_coordinates + parameter_coordinate
        ];

    std::vector<double> gradient{
        std::log(density) + a + density * a_rho,
        -density - density * density * a_rho,
        amount * a_parameter,
    };
    std::vector<double> hessian(9);
    hessian[0] = 1.0 / amount
        + (2.0 * a_rho + density * a_rho_rho) / volume;
    hessian[1] = -1.0 / volume
        - 2.0 * density * a_rho / volume
        - density * density * a_rho_rho / volume;
    hessian[2] = a_parameter + density * a_rho_parameter;
    hessian[3] = hessian[1];
    hessian[4] = density / volume
        * (1.0 + 2.0 * density * a_rho
           + density * density * a_rho_rho);
    hessian[5] = -density * density * a_rho_parameter;
    hessian[6] = hessian[2];
    hessian[7] = hessian[5];
    hessian[8] = amount * a_parameter_parameter;

    const double value =
        amount * (std::log(density) - 1.0 + derivatives.residual_helmholtz);
    const double pressure = -GAS_CONSTANT * temperature * gradient[1];
    if (!std::isfinite(value) || !std::isfinite(pressure)
        || !std::all_of(gradient.begin(), gradient.end(), [](double item) {
            return std::isfinite(item);
        })
        || !std::all_of(hessian.begin(), hessian.end(), [](double item) {
            return std::isfinite(item);
        })) {
        throw std::runtime_error(
            "active association phase block contains nonfinite values"
        );
    }
    return {value, std::move(gradient), std::move(hessian), pressure};
}

// EqID: phase_packing_fraction
// EqExclude: CppAD-generated derivative machinery
inline NativePackingFractionBlock packing_fraction_block(
    double temperature,
    const std::vector<double>& amounts,
    double volume,
    double ion_mole_fraction_max,
    const NativeParameters& parameters
) {
    if (!supports_electrolyte_phase_block(parameters)) {
        throw std::logic_error(
            "packing fraction block requires an ionic associating model"
        );
    }
    validate_electrolyte_phase_state(
        temperature, amounts, volume, ion_mole_fraction_max, parameters
    );
    const double checked_value = phase_packing_fraction(
        temperature, amounts, volume, parameters
    );
    if (!std::isfinite(checked_value)
        || checked_value >= close_packing_fraction()) {
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
        independent.begin(), independent.begin() + volume_coordinate
    );
    const ADScalar active_value = phase_packing_fraction(
        temperature,
        active_amounts,
        independent[volume_coordinate],
        parameters
    );
    CppAD::ADFun<double> function(independent, {active_value});
    const double value = function.Forward(0, point)[0];
    const auto gradient = function.Jacobian(point);
    const auto hessian = function.Hessian(point, 0);
    if (!std::isfinite(value)
        || !std::all_of(gradient.begin(), gradient.end(), [](double item) {
            return std::isfinite(item);
        })
        || !std::all_of(hessian.begin(), hessian.end(), [](double item) {
            return std::isfinite(item);
        })) {
        throw std::runtime_error("packing fraction block contains nonfinite values");
    }
    return {value, gradient, hessian};
}

template <typename Scalar>
// EqID: pure_active_parameter_mapping
PhaseBlockParameters<Scalar> pure_phase_parameters(
    const Scalar& segment_count,
    const Scalar& segment_diameter_angstrom,
    const Scalar& dispersion_energy_over_k
) {
    return {
        {segment_count},
        {segment_diameter_angstrom},
        {},
        {},
        {},
        {dispersion_energy_over_k},
        {Scalar(0.0)},
        {Scalar(0.0)},
        {},
        {},
        {},
        {0},
        ELECTROLYTE_MODEL_NONE,
    };
}

template <typename Scalar, typename Parameters>
Scalar pure_phase_helmholtz_over_rt_reference_amount(
    double temperature,
    const Scalar& amount,
    const Scalar& volume,
    const Parameters& parameters
) {
    // EqID: phase_molar_density
    const Scalar density = amount / volume;
    const std::vector<Scalar> composition = {Scalar(1.0)};
    const auto sigma = segment_sigma(Scalar(temperature), parameters);
    const auto diameter = segment_diameter(
        Scalar(temperature), parameters, sigma
    );
    const auto contributions = residual_helmholtz(
        Scalar(temperature),
        density,
        composition,
        parameters,
        sigma,
        diameter,
        {},
        Scalar(0.0),
        Scalar(0.0)
    );
    const Scalar residual = contributions.hard_chain + contributions.dispersion;
    // EqID: phase_dimensionless_helmholtz
    return amount * (scalar_log(density) - Scalar(1.0) + residual);
}

// EqExclude: CppAD-generated derivative machinery
inline double cubic_directional_derivative(
    CppAD::ADFun<double>& function,
    const std::array<double, 2>& point,
    const std::array<double, 2>& direction
) {
    const std::vector<double> point_vector(point.begin(), point.end());
    const std::vector<double> direction_vector(direction.begin(), direction.end());
    const std::vector<double> zero(2, 0.0);
    function.Forward(0, point_vector);
    function.Forward(1, direction_vector);
    function.Forward(2, zero);
    return 6.0 * function.Forward(3, zero)[0];
}

// EqID: phase_pressure
// EqExclude: CppAD-generated derivative machinery
inline NativePurePhaseBlock pure_phase_block(
    double temperature,
    double amount,
    double volume,
    const NativeParameters& parameters
) {
    if (!supports_pure_phase_block(parameters)) {
        throw std::logic_error("native SDK v1 requires one pure nonassociating neutral model");
    }
    if (!std::isfinite(temperature) || !std::isfinite(amount)
        || !std::isfinite(volume) || temperature <= 0.0 || amount <= 0.0
        || volume <= 0.0) {
        throw std::invalid_argument("temperature, amount, and volume must be positive and finite");
    }

    pure_phase_helmholtz_over_rt_reference_amount(
        temperature, amount, volume, parameters
    );

    CppAdSerialGuard cppad_guard;
    using ADScalar = CppAD::AD<double>;
    std::vector<ADScalar> independent = {amount, volume};
    CppAD::Independent(independent);
    const ADScalar phase = pure_phase_helmholtz_over_rt_reference_amount(
        temperature, independent[0], independent[1], parameters
    );
    CppAD::ADFun<double> function(independent, {phase});

    const std::array<double, 2> point = {amount, volume};
    const std::vector<double> point_vector(point.begin(), point.end());
    const auto jacobian = function.Jacobian(point_vector);
    const auto hessian_vector = function.Hessian(point_vector, 0);
    const std::array<double, 2> gradient = {jacobian[0], jacobian[1]};
    const std::array<double, 4> hessian = {
        hessian_vector[0], hessian_vector[1],
        hessian_vector[2], hessian_vector[3],
    };

    const double t_nnn = cubic_directional_derivative(
        function, point, {1.0, 0.0}
    );
    const double t_vvv = cubic_directional_derivative(
        function, point, {0.0, 1.0}
    );
    const double plus = cubic_directional_derivative(
        function, point, {1.0, 1.0}
    );
    const double minus = cubic_directional_derivative(
        function, point, {1.0, -1.0}
    );
    const double sum_mixed = (plus - t_nnn - t_vvv) / 3.0;
    const double mixed_difference = (minus - t_nnn + t_vvv) / 3.0;
    const double t_nnv = 0.5 * (sum_mixed - mixed_difference);
    const double t_nvv = 0.5 * (sum_mixed + mixed_difference);
    const std::array<double, 8> third = {
        t_nnn, t_nnv, t_nnv, t_nvv,
        t_nnv, t_nvv, t_nvv, t_vvv,
    };

    const double value = function.Forward(0, point_vector)[0];
    if (!std::isfinite(value)
        || !std::all_of(gradient.begin(), gradient.end(), [](double item) {
            return std::isfinite(item);
        })
        || !std::all_of(hessian.begin(), hessian.end(), [](double item) {
            return std::isfinite(item);
        })
        || !std::all_of(third.begin(), third.end(), [](double item) {
            return std::isfinite(item);
        })) {
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
// EqExclude: CppAD-generated derivative machinery
inline NativeParameterizedPurePhaseBlock parameterized_pure_phase_block(
    double temperature,
    double amount,
    double volume,
    double segment_count,
    double segment_diameter_angstrom,
    double dispersion_energy_over_k,
    const NativeParameters& base_parameters
) {
    if (!supports_pure_phase_block(base_parameters)) {
        throw std::logic_error(
            "parameterized phase block requires one pure nonassociating neutral model"
        );
    }
    const std::array<double, 6> inputs = {
        temperature,
        amount,
        volume,
        segment_count,
        segment_diameter_angstrom,
        dispersion_energy_over_k,
    };
    if (!std::all_of(inputs.begin(), inputs.end(), [](double item) {
            return std::isfinite(item) && item > 0.0;
        })) {
        throw std::invalid_argument(
            "parameterized phase inputs must be positive and finite"
        );
    }

    const auto value_parameters = pure_phase_parameters(
        segment_count, segment_diameter_angstrom, dispersion_energy_over_k
    );
    pure_phase_helmholtz_over_rt_reference_amount(
        temperature, amount, volume, value_parameters
    );

    CppAdSerialGuard cppad_guard;
    using ADScalar = CppAD::AD<double>;
    std::vector<ADScalar> independent = {
        amount,
        volume,
        segment_count,
        segment_diameter_angstrom,
        dispersion_energy_over_k,
    };
    CppAD::Independent(independent);
    const auto active_parameters = pure_phase_parameters(
        independent[2], independent[3], independent[4]
    );
    const ADScalar phase = pure_phase_helmholtz_over_rt_reference_amount(
        temperature, independent[0], independent[1], active_parameters
    );
    CppAD::ADFun<double> function(independent, {phase});

    const std::array<double, 5> point = {
        amount,
        volume,
        segment_count,
        segment_diameter_angstrom,
        dispersion_energy_over_k,
    };
    const std::vector<double> point_vector(point.begin(), point.end());
    const auto jacobian = function.Jacobian(point_vector);
    const auto hessian_vector = function.Hessian(point_vector, 0);
    std::array<double, 5> gradient{};
    std::array<double, 25> hessian{};
    std::copy(jacobian.begin(), jacobian.end(), gradient.begin());
    std::copy(hessian_vector.begin(), hessian_vector.end(), hessian.begin());

    const double value = function.Forward(0, point_vector)[0];
    if (!std::isfinite(value)
        || !std::all_of(gradient.begin(), gradient.end(), [](double item) {
            return std::isfinite(item);
        })
        || !std::all_of(hessian.begin(), hessian.end(), [](double item) {
            return std::isfinite(item);
        })) {
        throw std::runtime_error(
            "parameterized phase block contains nonfinite values"
        );
    }

    return {
        value,
        gradient,
        hessian,
        -GAS_CONSTANT * temperature * gradient[1],
        gradient[0],
    };
}

// EqID: phase_pressure
// EqID: association_mass_action
// EqExclude: CppAD-generated derivative machinery
inline NativeAssociatingParameterizedPurePhaseBlock associating_parameterized_pure_phase_block(
    double temperature,
    double amount,
    double volume,
    double segment_count,
    double segment_diameter_angstrom,
    double dispersion_energy_over_k,
    double association_energy_over_k,
    double association_volume,
    bool temperature_dependent_size_basis,
    const NativeParameters& base_parameters
) {
    if (!supports_associating_pure_phase_block(base_parameters)
        || !base_parameters.sigma_term_component.empty()) {
        throw std::logic_error(
            "joint associating parameters require one neutral constant-diameter 2B model"
        );
    }
    const std::array<double, 8> inputs{
        temperature,
        amount,
        volume,
        segment_count,
        segment_diameter_angstrom,
        dispersion_energy_over_k,
        association_energy_over_k,
        association_volume,
    };
    if (!std::all_of(inputs.begin(), inputs.end(), [](double item) {
            return std::isfinite(item) && item > 0.0;
        })) {
        throw std::invalid_argument(
            "joint associating phase inputs must be positive and finite"
        );
    }

    NativeParameters value_parameters = base_parameters;
    value_parameters.segment_count[0] = segment_count;
    value_parameters.sigma_constant[0] = segment_diameter_angstrom;
    value_parameters.dispersion_energy_over_k[0] =
        dispersion_energy_over_k;
    for (std::size_t pair = 0;
         pair < value_parameters.association_active.size();
         ++pair) {
        if (value_parameters.association_active[pair] != 0) {
            value_parameters.association_energy_over_k[pair] =
                association_energy_over_k;
            value_parameters.association_volume[pair] = association_volume;
        }
    }
    const double density = amount / volume;
    const std::vector<double> composition{1.0};
    const auto value_sigma = segment_sigma(temperature, value_parameters);
    const auto value_diameter = segment_diameter(
        temperature, value_parameters, value_sigma
    );
    const auto site_fractions = solve_site_fractions(
        temperature,
        density,
        composition,
        value_parameters,
        value_sigma,
        value_diameter,
        temperature_dependent_size_basis
    );

    CppAdSerialGuard cppad_guard;
    using ADScalar = CppAD::AD<double>;
    constexpr std::size_t coordinates = 7;
    const std::size_t sites = site_fractions.size();
    const std::size_t width = coordinates + sites;
    std::vector<ADScalar> independent(width);
    independent[0] = amount;
    independent[1] = volume;
    independent[2] = segment_count;
    independent[3] = segment_diameter_angstrom;
    independent[4] = dispersion_energy_over_k;
    independent[5] = association_energy_over_k;
    independent[6] = association_volume;
    for (std::size_t site = 0; site < sites; ++site) {
        independent[coordinates + site] = site_fractions[site];
    }
    CppAD::Independent(independent);

    PhaseBlockParameters<ADScalar> parameters{
        {independent[2]},
        {independent[3]},
        {},
        {},
        {},
        {independent[4]},
        {ADScalar(0.0)},
        {ADScalar(0.0)},
        base_parameters.site_component,
        base_parameters.site_multiplicity,
        base_parameters.association_active,
        {0},
        ELECTROLYTE_MODEL_NONE,
    };
    const ADScalar active_density = independent[0] / independent[1];
    const std::vector<ADScalar> active_composition{ADScalar(1.0)};
    const std::vector<ADScalar> active_site_fractions(
        independent.begin() + static_cast<std::ptrdiff_t>(coordinates),
        independent.end()
    );
    const auto sigma = segment_sigma(ADScalar(temperature), parameters);
    const auto diameter = segment_diameter(
        ADScalar(temperature), parameters, sigma
    );
    const auto contributions = residual_helmholtz(
        ADScalar(temperature),
        active_density,
        active_composition,
        parameters,
        sigma,
        diameter,
        active_site_fractions,
        ADScalar(0.0),
        ADScalar(0.0)
    );
    const ADScalar phase = independent[0] * (
        CppAD::log(active_density) - ADScalar(1.0)
        + contributions.hard_chain
        + contributions.dispersion
        + contributions.association
    );
    std::vector<ADScalar> association_energies(
        base_parameters.association_energy_over_k.begin(),
        base_parameters.association_energy_over_k.end()
    );
    std::vector<ADScalar> association_volumes(
        base_parameters.association_volume.begin(),
        base_parameters.association_volume.end()
    );
    for (std::size_t pair = 0;
         pair < base_parameters.association_active.size();
         ++pair) {
        if (base_parameters.association_active[pair] != 0) {
            association_energies[pair] = independent[5];
            association_volumes[pair] = independent[6];
        }
    }
    const auto delta = association_delta_values(
        ADScalar(temperature),
        active_density,
        active_composition,
        parameters,
        sigma,
        diameter,
        association_energies,
        association_volumes,
        temperature_dependent_size_basis
    );
    const auto mass_action = association_residuals(
        active_density,
        active_composition,
        parameters,
        active_site_fractions,
        delta
    );
    std::vector<ADScalar> dependent{phase};
    dependent.insert(dependent.end(), mass_action.begin(), mass_action.end());
    CppAD::ADFun<double> function(independent, dependent);

    std::vector<double> point{
        amount,
        volume,
        segment_count,
        segment_diameter_angstrom,
        dispersion_energy_over_k,
        association_energy_over_k,
        association_volume,
    };
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
        coordinates, std::vector<double>(sites, 0.0)
    );
    for (std::size_t coordinate = 0; coordinate < coordinates; ++coordinate) {
        std::vector<double> right_hand_side(sites);
        for (std::size_t row = 0; row < sites; ++row) {
            right_hand_side[row] =
                -jacobian[(1 + row) * width + coordinate];
        }
        site_sensitivity[coordinate] = solve_linear(
            association_jacobian, std::move(right_hand_side)
        );
    }
    std::vector<std::vector<double>> directions(
        coordinates, std::vector<double>(width, 0.0)
    );
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
    std::array<double, coordinates> gradient{};
    std::array<double, coordinates * coordinates> hessian{};
    for (std::size_t coordinate = 0; coordinate < coordinates; ++coordinate) {
        gradient[coordinate] = jacobian[coordinate];
        for (std::size_t site = 0; site < sites; ++site) {
            gradient[coordinate] += jacobian[coordinates + site]
                * site_sensitivity[coordinate][site];
        }
    }
    for (std::size_t left = 0; left < coordinates; ++left) {
        for (std::size_t right = 0; right < coordinates; ++right) {
            std::vector<double> right_hand_side(sites);
            for (std::size_t site = 0; site < sites; ++site) {
                right_hand_side[site] = -directional_second(
                    association_hessians[site],
                    directions[left],
                    directions[right]
                );
            }
            const auto site_second = solve_linear(
                association_jacobian, std::move(right_hand_side)
            );
            double second = directional_second(
                phase_hessian, directions[left], directions[right]
            );
            for (std::size_t site = 0; site < sites; ++site) {
                second += jacobian[coordinates + site] * site_second[site];
            }
            hessian[left * coordinates + right] = second;
        }
    }
    const double value = function.Forward(0, point)[0];
    if (!std::isfinite(value)
        || !std::all_of(gradient.begin(), gradient.end(), [](double item) {
            return std::isfinite(item);
        })
        || !std::all_of(hessian.begin(), hessian.end(), [](double item) {
            return std::isfinite(item);
        })) {
        throw std::runtime_error(
            "joint associating phase block contains nonfinite values"
        );
    }
    return {
        value,
        gradient,
        hessian,
        -GAS_CONSTANT * temperature * gradient[1],
        gradient[0],
    };
}

}  // namespace detail
}  // namespace epcsaft
