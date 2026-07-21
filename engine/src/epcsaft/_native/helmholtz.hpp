#pragma once

#include "model.hpp"

#include <cppad/cppad.hpp>

#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <type_traits>
#include <vector>

namespace epcsaft::detail {

inline constexpr double AVOGADRO = 6.02214076e23;
inline constexpr double GAS_CONSTANT = 8.31446261815324;
inline constexpr double PI = 3.141592653589793238462643383279502884;

inline constexpr std::array<double, 7> A0 = {
    0.9105631445, 0.6361281449, 2.6861347891, -26.547362491,
    97.759208784, -159.59154087, 91.297774084,
};
inline constexpr std::array<double, 7> A1 = {
    -0.3084016918, 0.1860531159, -2.5030047259, 21.419793629,
    -65.255885330, 83.318680481, -33.746922930,
};
inline constexpr std::array<double, 7> A2 = {
    -0.0906148351, 0.4527842806, 0.5962700728, -1.7241829131,
    -4.1302112531, 13.776631870, -8.6728470368,
};
inline constexpr std::array<double, 7> B0 = {
    0.7240946941, 2.2382791861, -4.0025849485, -21.003576815,
    26.855641363, 206.55133841, -355.60235612,
};
inline constexpr std::array<double, 7> B1 = {
    -0.5755498075, 0.6995095521, 3.8925673390, -17.215471648,
    192.67226447, -161.82646165, -165.20769346,
};
inline constexpr std::array<double, 7> B2 = {
    0.0976883116, -0.2557574982, -9.1558561530, 20.642075974,
    -38.804430052, 93.626774077, -29.666905585,
};

template <typename Scalar>
struct HelmholtzResult final {
    Scalar hard_chain;
    Scalar dispersion;
    Scalar association;
    Scalar debye_huckel;
    Scalar born_ssm_ds;
};

template <typename Scalar>
Scalar scalar_log(const Scalar& value) {
    if constexpr (std::is_same_v<Scalar, double>) {
        return std::log(value);
    }
    return CppAD::log(value);
}

template <typename Scalar>
Scalar scalar_exp(const Scalar& value) {
    if constexpr (std::is_same_v<Scalar, double>) {
        return std::exp(value);
    }
    return CppAD::exp(value);
}

template <typename Scalar>
Scalar scalar_sqrt(const Scalar& value) {
    if constexpr (std::is_same_v<Scalar, double>) {
        return std::sqrt(value);
    }
    return CppAD::sqrt(value);
}

template <typename Scalar>
Scalar polynomial(const std::array<Scalar, 7>& coefficients, const Scalar& eta) {
    Scalar value = coefficients.back();
    for (std::size_t index = coefficients.size() - 1; index > 0; --index) {
        value = value * eta + coefficients[index - 1];
    }
    return value;
}

template <typename Scalar, typename Parameters>
std::vector<Scalar> segment_sigma(
    const Scalar& temperature,
    const Parameters& parameters
) {
    std::vector<Scalar> sigma(
        parameters.sigma_constant.begin(), parameters.sigma_constant.end()
    );
    for (std::size_t term = 0; term < parameters.sigma_term_component.size(); ++term) {
        const std::size_t component =
            static_cast<std::size_t>(parameters.sigma_term_component[term]);
        sigma[component] += parameters.sigma_term_amplitude[term] *
            scalar_exp(Scalar(parameters.sigma_term_exponent[term]) * temperature);
    }
    return sigma;
}

template <typename Scalar, typename Parameters>
std::vector<Scalar> segment_diameter(
    const Scalar& temperature,
    const Parameters& parameters,
    const std::vector<Scalar>& sigma
) {
    std::vector<Scalar> diameter(sigma.size());
    for (std::size_t component = 0; component < sigma.size(); ++component) {
        diameter[component] = parameters.charge_number[component] == 0
            ? sigma[component] *
                (Scalar(1.0) - Scalar(0.12) * scalar_exp(
                    Scalar(-3.0 * parameters.dispersion_energy_over_k[component])
                        / temperature
                ))
            : Scalar(0.88) * sigma[component];
    }
    return diameter;
}

inline double close_packing_fraction() {
    return PI / (3.0 * std::sqrt(2.0));
}

template <typename Scalar, typename Parameters>
Scalar packing_fraction_from_geometry(
    const Scalar& molar_density,
    const std::vector<Scalar>& mole_fractions,
    const Parameters& parameters,
    const std::vector<Scalar>& diameter
) {
    Scalar moment = Scalar(0.0);
    for (std::size_t component = 0; component < mole_fractions.size(); ++component) {
        moment += mole_fractions[component]
            * parameters.segment_count[component]
            * diameter[component] * diameter[component] * diameter[component];
    }
    return (PI / 6.0) * molar_density * (AVOGADRO / 1.0e30) * moment;
}

template <typename Parameters>
double packing_fraction_coefficient(
    double temperature,
    const std::vector<double>& mole_fractions,
    const Parameters& parameters
) {
    const auto sigma = segment_sigma(temperature, parameters);
    const auto diameter = segment_diameter(temperature, parameters, sigma);
    for (std::size_t component = 0; component < mole_fractions.size(); ++component) {
        if (!(sigma[component] > 0.0) || !std::isfinite(sigma[component])
            || !(diameter[component] > 0.0)
            || !std::isfinite(diameter[component])) {
            throw std::invalid_argument("invalid segment diameter correlation");
        }
    }
    return packing_fraction_from_geometry(
        1.0, mole_fractions, parameters, diameter
    );
}

template <typename Scalar, typename Parameters>
Scalar phase_packing_fraction(
    double temperature,
    const std::vector<Scalar>& amounts,
    const Scalar& volume,
    const Parameters& parameters
) {
    Scalar total = Scalar(0.0);
    for (const Scalar& amount : amounts) {
        total += amount;
    }
    std::vector<Scalar> composition;
    composition.reserve(amounts.size());
    for (const Scalar& amount : amounts) {
        composition.push_back(amount / total);
    }
    const auto sigma = segment_sigma(Scalar(temperature), parameters);
    const auto diameter = segment_diameter(Scalar(temperature), parameters, sigma);
    return packing_fraction_from_geometry(
        total / volume, composition, parameters, diameter
    );
}

template <typename Scalar, typename Parameters>
HelmholtzResult<Scalar> residual_helmholtz(
    const Scalar& temperature,
    const Scalar& molar_density,
    const std::vector<Scalar>& mole_fractions,
    const Parameters& parameters,
    const std::vector<Scalar>& sigma,
    const std::vector<Scalar>& diameter,
    const std::vector<Scalar>& site_fraction,
    const Scalar& debye_huckel,
    const Scalar& born_ssm_ds
) {
    const std::size_t size = parameters.segment_count.size();
    const Scalar number_density = molar_density * (AVOGADRO / 1.0e30);
    Scalar mean_segment_count = Scalar(0.0);
    for (std::size_t i = 0; i < size; ++i) {
        mean_segment_count += mole_fractions[i] * parameters.segment_count[i];
    }

    std::array<Scalar, 4> zeta = {
        Scalar(0.0), Scalar(0.0), Scalar(0.0), Scalar(0.0),
    };
    for (std::size_t order = 0; order < zeta.size() - 1; ++order) {
        Scalar moment = Scalar(0.0);
        for (std::size_t i = 0; i < size; ++i) {
            Scalar power = Scalar(1.0);
            for (std::size_t exponent = 0; exponent < order; ++exponent) {
                power *= diameter[i];
            }
            moment += mole_fractions[i] * parameters.segment_count[i] * power;
        }
        zeta[order] = (PI / 6.0) * number_density * moment;
    }
    zeta[3] = packing_fraction_from_geometry(
        molar_density, mole_fractions, parameters, diameter
    );
    if constexpr (std::is_same_v<Scalar, double>) {
        if (zeta[3] >= close_packing_fraction()) {
            throw std::invalid_argument("packing fraction is at or above close packing");
        }
    }

    const Scalar one_minus_eta = Scalar(1.0) - zeta[3];
    const Scalar hard_sphere =
        (3.0 * zeta[1] * zeta[2] / one_minus_eta +
         zeta[2] * zeta[2] * zeta[2] /
             (zeta[3] * one_minus_eta * one_minus_eta) +
         (zeta[2] * zeta[2] * zeta[2] / (zeta[3] * zeta[3]) - zeta[0]) *
             scalar_log(one_minus_eta)) /
        zeta[0];

    Scalar hard_chain = mean_segment_count * hard_sphere;
    for (std::size_t i = 0; i < size; ++i) {
        const Scalar contact_diameter = diameter[i] / Scalar(2.0);
        const Scalar contact =
            Scalar(1.0) / one_minus_eta +
            contact_diameter * 3.0 * zeta[2] /
                (one_minus_eta * one_minus_eta) +
            contact_diameter * contact_diameter * 2.0 * zeta[2] * zeta[2] /
                (one_minus_eta * one_minus_eta * one_minus_eta);
        hard_chain -= mole_fractions[i] * (parameters.segment_count[i] - 1.0) *
            scalar_log(contact);
    }

    const Scalar chain_factor_1 =
        (mean_segment_count - 1.0) / mean_segment_count;
    const Scalar chain_factor_2 = chain_factor_1 *
        (mean_segment_count - 2.0) / mean_segment_count;
    std::array<Scalar, 7> a{};
    std::array<Scalar, 7> b{};
    for (std::size_t i = 0; i < a.size(); ++i) {
        a[i] = A0[i] + chain_factor_1 * A1[i] + chain_factor_2 * A2[i];
        b[i] = B0[i] + chain_factor_1 * B1[i] + chain_factor_2 * B2[i];
    }
    const Scalar integral_1 = polynomial(a, zeta[3]);
    const Scalar integral_2 = polynomial(b, zeta[3]);
    const Scalar eta_2 = zeta[3] * zeta[3];
    const Scalar eta_3 = eta_2 * zeta[3];
    const Scalar eta_4 = eta_3 * zeta[3];
    const Scalar c1 = Scalar(1.0) /
        (Scalar(1.0) +
         mean_segment_count * (8.0 * zeta[3] - 2.0 * eta_2) /
             (one_minus_eta * one_minus_eta * one_minus_eta * one_minus_eta) +
         (1.0 - mean_segment_count) *
             (20.0 * zeta[3] - 27.0 * eta_2 + 12.0 * eta_3 - 2.0 * eta_4) /
             ((one_minus_eta * (Scalar(2.0) - zeta[3])) *
              (one_minus_eta * (Scalar(2.0) - zeta[3]))));

    Scalar first_moment = Scalar(0.0);
    Scalar second_moment = Scalar(0.0);
    for (std::size_t i = 0; i < size; ++i) {
        for (std::size_t j = 0; j < size; ++j) {
            const Scalar sigma_ij = Scalar(0.5) * (sigma[i] + sigma[j])
                * (Scalar(1.0) - parameters.l_ij[i * size + j]);
            const Scalar epsilon_ij =
                parameters.charge_number[i] * parameters.charge_number[j] > 0
                ? Scalar(0.0)
                : scalar_sqrt(
                    parameters.dispersion_energy_over_k[i] *
                    parameters.dispersion_energy_over_k[j]
                ) * (1.0 - parameters.k_ij[i * size + j]);
            const Scalar weight = mole_fractions[i] * mole_fractions[j] *
                parameters.segment_count[i] * parameters.segment_count[j] *
                sigma_ij * sigma_ij * sigma_ij;
            first_moment += weight * epsilon_ij / temperature;
            second_moment += weight * epsilon_ij * epsilon_ij /
                (temperature * temperature);
        }
    }
    const Scalar dispersion =
        -2.0 * PI * number_density * integral_1 * first_moment -
        PI * number_density * mean_segment_count * c1 * integral_2 * second_moment;

    Scalar association = 0.0;
    for (std::size_t site = 0; site < site_fraction.size(); ++site) {
        const std::size_t component =
            static_cast<std::size_t>(parameters.site_component[site]);
        association += mole_fractions[component] *
            parameters.site_multiplicity[site] *
            (scalar_log(site_fraction[site]) - site_fraction[site] / 2.0 + 0.5);
    }
    return {hard_chain, dispersion, association, debye_huckel, born_ssm_ds};
}

}  // namespace epcsaft::detail
