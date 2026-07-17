#pragma once

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

namespace detail {

inline bool supports_pure_phase_block(const NativeParameters& parameters) {
    return parameters.segment_count.size() == 1
        && parameters.site_component.empty()
        && parameters.charge_number.size() == 1
        && parameters.charge_number[0] == 0;
}

template <typename Scalar>
Scalar pure_phase_helmholtz_over_rt_reference_amount(
    double temperature,
    const Scalar& amount,
    const Scalar& volume,
    const NativeParameters& parameters
) {
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
    return amount * (scalar_log(density) - Scalar(1.0) + residual);
}

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

}  // namespace detail
}  // namespace epcsaft
