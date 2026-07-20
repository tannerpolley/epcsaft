#pragma once

#include "association.hpp"
#include "born.hpp"
#include "electrolyte.hpp"
#include "helmholtz.hpp"
#include "model.hpp"

#include <cppad/cppad.hpp>

#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

namespace epcsaft::detail {

struct ThermodynamicDerivatives final {
    double temperature;
    double density;
    std::vector<double> composition;
    double temperature_temperature;
    double temperature_density;
    double density_density;
    double residual_helmholtz;
    std::vector<double> hessian;
};

inline double directional_second(
    const std::vector<double>& hessian,
    const std::vector<double>& left,
    const std::vector<double>& right
) {
    const std::size_t size = left.size();
    double value = 0.0;
    for (std::size_t row = 0; row < size; ++row) {
        for (std::size_t column = 0; column < size; ++column) {
            value += left[row] * hessian[row * size + column] * right[column];
        }
    }
    return value;
}

inline ThermodynamicDerivatives thermodynamic_derivatives(
    double temperature,
    double molar_density,
    const std::vector<double>& mole_fractions,
    const NativeParameters& parameters,
    const std::vector<double>& site_fractions
) {
    using ADScalar = CppAD::AD<double>;
    const std::size_t components = mole_fractions.size();
    const std::size_t sites = site_fractions.size();
    const std::size_t coordinates = 2 + components;
    const std::size_t width = coordinates + sites;

    std::vector<ADScalar> independent(width);
    independent[0] = temperature;
    independent[1] = molar_density;
    for (std::size_t component = 0; component < components; ++component) {
        independent[2 + component] = mole_fractions[component];
    }
    for (std::size_t site = 0; site < sites; ++site) {
        independent[coordinates + site] = site_fractions[site];
    }
    CppAD::Independent(independent);

    const std::vector<ADScalar> composition(
        independent.begin() + 2, independent.begin() + coordinates
    );
    const std::vector<ADScalar> fractions(
        independent.begin() + coordinates, independent.end()
    );
    const auto sigma = segment_sigma(independent[0], parameters);
    const auto diameter = segment_diameter(independent[0], parameters, sigma);
    const auto reference_sigma = segment_sigma(temperature, parameters);
    const auto electrolyte = electrolyte_contribution(
        independent[0],
        independent[1],
        temperature,
        molar_density,
        composition,
        mole_fractions,
        parameters,
        sigma,
        reference_sigma
    );
    const auto born = born_ssm_ds_contribution(
        independent[0], composition, parameters,
        electrolyte.bulk_relative_permittivity
    );
    const auto contributions = residual_helmholtz(
        independent[0],
        independent[1],
        composition,
        parameters,
        sigma,
        diameter,
        fractions,
        electrolyte.debye_huckel,
        born
    );
    const auto delta = association_delta(
        independent[0],
        independent[1],
        composition,
        parameters,
        sigma,
        diameter
    );
    const auto mass_action = association_residuals(
        independent[1], composition, parameters, fractions, delta
    );
    std::vector<ADScalar> dependent = {
        contributions.hard_chain + contributions.dispersion
            + contributions.association + contributions.debye_huckel
            + contributions.born_ssm_ds,
    };
    dependent.insert(dependent.end(), mass_action.begin(), mass_action.end());
    CppAD::ADFun<double> function(independent, dependent);

    std::vector<double> point = {temperature, molar_density};
    point.insert(point.end(), mole_fractions.begin(), mole_fractions.end());
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
        if (sites == 0) {
            break;
        }
        std::vector<double> right_hand_side(sites);
        for (std::size_t row = 0; row < sites; ++row) {
            right_hand_side[row] = -jacobian[(1 + row) * width + coordinate];
        }
        site_sensitivity[coordinate] = solve_linear(
            association_jacobian, std::move(right_hand_side)
        );
    }

    std::vector<double> gradient(coordinates);
    for (std::size_t coordinate = 0; coordinate < coordinates; ++coordinate) {
        gradient[coordinate] = jacobian[coordinate];
        for (std::size_t site = 0; site < sites; ++site) {
            gradient[coordinate] += jacobian[coordinates + site]
                * site_sensitivity[coordinate][site];
        }
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
    const auto helmholtz_hessian = function.Hessian(point, 0);

    std::vector<std::vector<double>> association_hessians(sites);
    for (std::size_t site = 0; site < sites; ++site) {
        association_hessians[site] = function.Hessian(point, 1 + site);
    }
    auto reduced_second = [&](std::size_t left, std::size_t right) {
        std::vector<double> site_second(sites, 0.0);
        if (sites != 0) {
            std::vector<double> right_hand_side(sites);
            for (std::size_t site = 0; site < sites; ++site) {
                right_hand_side[site] = -directional_second(
                    association_hessians[site], directions[left], directions[right]
                );
            }
            site_second = solve_linear(
                association_jacobian, std::move(right_hand_side)
            );
        }
        double value = directional_second(
            helmholtz_hessian, directions[left], directions[right]
        );
        for (std::size_t site = 0; site < sites; ++site) {
            value += jacobian[coordinates + site] * site_second[site];
        }
        return value;
    };

    std::vector<double> reduced_hessian(coordinates * coordinates);
    for (std::size_t row = 0; row < coordinates; ++row) {
        for (std::size_t column = 0; column < coordinates; ++column) {
            reduced_hessian[row * coordinates + column] =
                reduced_second(row, column);
        }
    }
    const double residual_helmholtz = function.Forward(0, point)[0];

    return {
        gradient[0],
        gradient[1],
        std::vector<double>(gradient.begin() + 2, gradient.end()),
        reduced_hessian[0],
        reduced_hessian[1],
        reduced_hessian[coordinates + 1],
        residual_helmholtz,
        std::move(reduced_hessian),
    };
}

}  // namespace epcsaft::detail
