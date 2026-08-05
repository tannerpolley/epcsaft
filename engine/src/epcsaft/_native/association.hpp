#pragma once

#include "model.hpp"

#include <cppad/cppad.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace epcsaft::detail {

inline constexpr double ASSOCIATION_NUMBER_DENSITY_FACTOR = 6.02214076e-7;
inline constexpr double ASSOCIATION_PI = 3.141592653589793238462643383279502884;

template <typename Scalar>
// EqExclude: scalar utility dispatch
Scalar association_exp(const Scalar& value) {
    if constexpr (std::is_same_v<Scalar, double>) {
        return std::exp(value);
    }
    return CppAD::exp(value);
}

template <typename Scalar, typename Parameters>
std::vector<Scalar> association_delta_values(
    const Scalar& temperature,
    const Scalar& molar_density,
    const std::vector<Scalar>& mole_fractions,
    const Parameters& parameters,
    const std::vector<Scalar>& sigma,
    const std::vector<Scalar>& diameter,
    const std::vector<Scalar>& association_energy_over_k,
    const std::vector<Scalar>& association_volume
) {
    const std::size_t components = parameters.segment_count.size();
    const std::size_t sites = parameters.site_component.size();
    const Scalar number_density =
        molar_density * ASSOCIATION_NUMBER_DENSITY_FACTOR;
    Scalar zeta_2 = Scalar(0.0);
    Scalar zeta_3 = Scalar(0.0);
    for (std::size_t component = 0; component < components; ++component) {
        const Scalar weight = mole_fractions[component]
            * parameters.segment_count[component];
        zeta_2 += ASSOCIATION_PI / 6.0 * number_density * weight
            * diameter[component] * diameter[component];
        zeta_3 += ASSOCIATION_PI / 6.0 * number_density * weight
            * diameter[component] * diameter[component] * diameter[component];
    }
    const Scalar one_minus_eta = Scalar(1.0) - zeta_3;

    std::vector<Scalar> delta(sites * sites, Scalar(0.0));
    for (std::size_t a = 0; a < sites; ++a) {
        const auto i = static_cast<std::size_t>(parameters.site_component[a]);
        for (std::size_t b = 0; b < sites; ++b) {
            const std::size_t pair = a * sites + b;
            if (parameters.association_active[pair] == 0) {
                continue;
            }
            const auto j = static_cast<std::size_t>(parameters.site_component[b]);
            const Scalar reduced_diameter =
                diameter[i] * diameter[j] / (diameter[i] + diameter[j]);
            // EqID: association_contact
            const Scalar contact =
                Scalar(1.0) / one_minus_eta
                + 3.0 * reduced_diameter * zeta_2
                    / (one_minus_eta * one_minus_eta)
                + 2.0 * reduced_diameter * reduced_diameter * zeta_2 * zeta_2
                    / (one_minus_eta * one_minus_eta * one_minus_eta);
            const Scalar size_ij = Scalar(0.5) * (sigma[i] + sigma[j])
                * (Scalar(1.0) - parameters.l_ij[i * components + j]);
            // EqID: association_strength
            delta[pair] = contact
                * (association_exp(
                       association_energy_over_k[pair] / temperature
                   ) - Scalar(1.0))
                * size_ij * size_ij * size_ij
                * association_volume[pair];
        }
    }
    return delta;
}

template <typename Scalar>
// EqExclude: scalar utility dispatch
std::vector<Scalar> association_delta(
    const Scalar& temperature,
    const Scalar& molar_density,
    const std::vector<Scalar>& mole_fractions,
    const NativeParameters& parameters,
    const std::vector<Scalar>& sigma,
    const std::vector<Scalar>& diameter
) {
    return association_delta_values(
        temperature,
        molar_density,
        mole_fractions,
        parameters,
        sigma,
        diameter,
        std::vector<Scalar>(
            parameters.association_energy_over_k.begin(),
            parameters.association_energy_over_k.end()
        ),
        std::vector<Scalar>(
            parameters.association_volume.begin(),
            parameters.association_volume.end()
        )
    );
}

template <typename Scalar, typename Parameters>
std::vector<Scalar> association_residuals(
    const Scalar& molar_density,
    const std::vector<Scalar>& mole_fractions,
    const Parameters& parameters,
    const std::vector<Scalar>& site_fractions,
    const std::vector<Scalar>& delta
) {
    const std::size_t sites = parameters.site_component.size();
    const Scalar number_density =
        molar_density * ASSOCIATION_NUMBER_DENSITY_FACTOR;
    std::vector<Scalar> residuals(sites, Scalar(0.0));
    for (std::size_t a = 0; a < sites; ++a) {
        Scalar bonded_density = Scalar(0.0);
        for (std::size_t b = 0; b < sites; ++b) {
            const auto component =
                static_cast<std::size_t>(parameters.site_component[b]);
            bonded_density += mole_fractions[component]
                * parameters.site_multiplicity[b]
                * site_fractions[b]
                * delta[a * sites + b];
        }
        // EqID: association_mass_action
        residuals[a] = site_fractions[a]
            * (Scalar(1.0) + number_density * bonded_density)
            - Scalar(1.0);
    }
    return residuals;
}

// EqExclude: numerical iteration mechanics
inline std::vector<double> solve_site_fractions(
    double temperature,
    double molar_density,
    const std::vector<double>& mole_fractions,
    const NativeParameters& parameters,
    const std::vector<double>& sigma,
    const std::vector<double>& diameter
) {
    const std::size_t sites = parameters.site_component.size();
    if (sites == 0) {
        return {};
    }
    const auto delta = association_delta_values(
        temperature,
        molar_density,
        mole_fractions,
        parameters,
        sigma,
        diameter,
        parameters.association_energy_over_k,
        parameters.association_volume
    );
    const double number_density =
        molar_density * ASSOCIATION_NUMBER_DENSITY_FACTOR;
    std::vector<double> fractions(sites, 1.0);
    constexpr std::size_t max_iterations = 1000;
    constexpr double tolerance = 2.0e-13;
    for (std::size_t iteration = 0; iteration < max_iterations; ++iteration) {
        std::vector<double> next(sites, 1.0);
        for (std::size_t a = 0; a < sites; ++a) {
            double bonded_density = 0.0;
            for (std::size_t b = 0; b < sites; ++b) {
                const auto component =
                    static_cast<std::size_t>(parameters.site_component[b]);
                bonded_density += mole_fractions[component]
                    * parameters.site_multiplicity[b]
                    * fractions[b]
                    * delta[a * sites + b];
            }
            const double fixed_point = 1.0 / (1.0 + number_density * bonded_density);
            next[a] = 0.5 * (fractions[a] + fixed_point);
        }
        const auto residuals = association_residuals(
            molar_density, mole_fractions, parameters, next, delta
        );
        double residual_norm = 0.0;
        for (std::size_t site = 0; site < sites; ++site) {
            if (!std::isfinite(next[site]) || next[site] <= 0.0
                || next[site] > 1.0) {
                throw std::runtime_error("association solve left physical domain");
            }
            residual_norm = std::max(residual_norm, std::abs(residuals[site]));
        }
        fractions = std::move(next);
        if (residual_norm <= tolerance) {
            return fractions;
        }
    }
    throw std::runtime_error("association solve did not converge");
}

// EqExclude: numerical iteration mechanics
inline std::vector<double> solve_linear(
    std::vector<double> matrix,
    std::vector<double> right_hand_side
) {
    const std::size_t size = right_hand_side.size();
    double largest_pivot = 0.0;
    double smallest_pivot = std::numeric_limits<double>::infinity();
    for (std::size_t column = 0; column < size; ++column) {
        std::size_t pivot_row = column;
        for (std::size_t row = column + 1; row < size; ++row) {
            if (std::abs(matrix[row * size + column])
                > std::abs(matrix[pivot_row * size + column])) {
                pivot_row = row;
            }
        }
        const double pivot = std::abs(matrix[pivot_row * size + column]);
        largest_pivot = std::max(largest_pivot, pivot);
        smallest_pivot = std::min(smallest_pivot, pivot);
        if (!std::isfinite(pivot) || pivot == 0.0) {
            throw std::runtime_error("association Jacobian is singular");
        }
        if (pivot_row != column) {
            for (std::size_t entry = column; entry < size; ++entry) {
                std::swap(
                    matrix[column * size + entry],
                    matrix[pivot_row * size + entry]
                );
            }
            std::swap(right_hand_side[column], right_hand_side[pivot_row]);
        }
        for (std::size_t row = column + 1; row < size; ++row) {
            const double factor = matrix[row * size + column]
                / matrix[column * size + column];
            for (std::size_t entry = column; entry < size; ++entry) {
                matrix[row * size + entry] -=
                    factor * matrix[column * size + entry];
            }
            right_hand_side[row] -= factor * right_hand_side[column];
        }
    }
    if (smallest_pivot <= 1.0e-13 * largest_pivot) {
        throw std::runtime_error("association Jacobian is ill-conditioned");
    }
    std::vector<double> solution(size, 0.0);
    for (std::size_t offset = 0; offset < size; ++offset) {
        const std::size_t row = size - 1 - offset;
        double value = right_hand_side[row];
        for (std::size_t column = row + 1; column < size; ++column) {
            value -= matrix[row * size + column] * solution[column];
        }
        solution[row] = value / matrix[row * size + row];
    }
    return solution;
}

}  // namespace epcsaft::detail
