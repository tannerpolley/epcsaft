#pragma once

#include "physics/numerical_policy.hpp"
#include "physics/parameters.hpp"
#include "physics/physical_constants.hpp"
#include "physics/scalar.hpp"

#include <cppad/cppad.hpp>

#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <type_traits>
#include <vector>

namespace epcsaft {

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
    Scalar born;
};

template <typename Scalar>
// EqID: dispersion_polynomial
Scalar polynomial(const std::array<Scalar, 7>& coefficients, const Scalar& eta) {
    Scalar value = coefficients.back();
    for (std::size_t index = coefficients.size() - 1; index > 0; --index) {
        value = value * eta + coefficients[index - 1];
    }
    return value;
}

template <typename Scalar, typename ParameterSet>
// EqID: segment_sigma_correlation
std::vector<Scalar> segment_sigma(
    const Scalar& temperature,
    const ParameterSet& parameters
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

template <typename Scalar, typename ParameterSet>
// EqID: neutral_segment_diameter
// EqID: ionic_packing_diameter
std::vector<Scalar> segment_diameter(
    const Scalar& temperature,
    const ParameterSet& parameters,
    const std::vector<Scalar>& sigma,
    const std::vector<Scalar>* dispersion_energy_override = nullptr
) {
    std::vector<Scalar> diameter(sigma.size());
    for (std::size_t component = 0; component < sigma.size(); ++component) {
        diameter[component] = parameters.charge_number[component] == 0
            ? sigma[component] *
                (Scalar(1.0) - Scalar(0.12) * scalar_exp(
                    Scalar(-3.0)
                        * (dispersion_energy_override == nullptr
                               ? Scalar(parameters.dispersion_energy_over_k[component])
                               : (*dispersion_energy_override)[component])
                        / temperature
                ))
            : Scalar(parameters.packing_diameter_angstrom[component]);
    }
    return diameter;
}

template <typename Scalar, typename ParameterSet>
// EqID: packing_fraction
Scalar packing_fraction_from_geometry(
    const Scalar& molar_density,
    const std::vector<Scalar>& mole_fractions,
    const ParameterSet& parameters,
    const std::vector<Scalar>& diameter
) {
    Scalar moment = Scalar(0.0);
    for (std::size_t component = 0; component < mole_fractions.size(); ++component) {
        moment += mole_fractions[component]
            * parameters.segment_count[component]
            * diameter[component] * diameter[component] * diameter[component];
    }
    return (PI / 6.0) * molar_density * NUMBER_DENSITY_PER_ANGSTROM3 * moment;
}

template <typename Scalar, typename ParameterSet>
HelmholtzResult<Scalar> residual_helmholtz(
    const Scalar& temperature,
    const Scalar& molar_density,
    const std::vector<Scalar>& mole_fractions,
    const ParameterSet& parameters,
    const std::vector<Scalar>& sigma,
    const std::vector<Scalar>& diameter,
    const std::vector<Scalar>& site_fraction,
    const Scalar& debye_huckel,
    const Scalar& born,
    const std::vector<Scalar>* dispersion_energy_override = nullptr
) {
    const std::size_t size = parameters.segment_count.size();
    // EqID: angstrom_number_density
    const Scalar number_density = molar_density * NUMBER_DENSITY_PER_ANGSTROM3;
    // EqID: mean_segment_count
    Scalar mean_segment_count = Scalar(0.0);
    for (std::size_t i = 0; i < size; ++i) {
        mean_segment_count += mole_fractions[i] * parameters.segment_count[i];
    }

    // EqID: packing_moments
    std::array<Scalar, 4> zeta = {
        Scalar(0.0), Scalar(0.0), Scalar(0.0), Scalar(0.0),
    };
    const std::size_t moment_count =
        size == 1 ? zeta.size() : zeta.size() - 1;
    for (std::size_t order = 0; order < moment_count; ++order) {
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
    if (size != 1) {
        zeta[3] = packing_fraction_from_geometry(
            molar_density, mole_fractions, parameters, diameter
        );
    }
    if constexpr (std::is_same_v<Scalar, double>) {
        if (zeta[3] >= policy::CLOSE_PACKING_FRACTION.value) {
            throw std::invalid_argument("packing fraction is at or above close packing");
        }
    }

    const Scalar one_minus_eta = Scalar(1.0) - zeta[3];
    // EqID: hard_sphere_helmholtz
    const Scalar hard_sphere =
        (3.0 * zeta[1] * zeta[2] / one_minus_eta +
         zeta[2] * zeta[2] * zeta[2] /
             (zeta[3] * one_minus_eta * one_minus_eta) +
         (zeta[2] * zeta[2] * zeta[2] / (zeta[3] * zeta[3]) - zeta[0]) *
             scalar_log(one_minus_eta)) /
        zeta[0];

    // EqID: hard_chain_helmholtz
    Scalar hard_chain = mean_segment_count * hard_sphere;
    for (std::size_t i = 0; i < size; ++i) {
        const Scalar contact_diameter = diameter[i] / Scalar(2.0);
        // EqID: hard_sphere_contact
        const Scalar contact =
            Scalar(1.0) / one_minus_eta +
            contact_diameter * 3.0 * zeta[2] /
                (one_minus_eta * one_minus_eta) +
            contact_diameter * contact_diameter * 2.0 * zeta[2] * zeta[2] /
                (one_minus_eta * one_minus_eta * one_minus_eta);
        hard_chain -= mole_fractions[i] * (parameters.segment_count[i] - 1.0) *
            scalar_log(contact);
    }

    // EqID: dispersion_chain_coefficients
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
    // EqID: dispersion_integral_one
    const Scalar integral_1 = polynomial(a, zeta[3]);
    // EqID: dispersion_integral_two
    const Scalar integral_2 = polynomial(b, zeta[3]);
    const Scalar eta_2 = zeta[3] * zeta[3];
    const Scalar eta_3 = eta_2 * zeta[3];
    const Scalar eta_4 = eta_3 * zeta[3];
    // EqID: dispersion_compressibility_factor
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
            // EqID: cross_segment_diameter
            const Scalar sigma_ij = Scalar(0.5) * (sigma[i] + sigma[j])
                * (Scalar(1.0) - parameters.l_ij[i * size + j]);
            const std::size_t pair = i * size + j;
            const double reference_temperature =
                parameters.k_ij_reference_temperature_k[pair];
            Scalar k_ij = parameters.k_ij[pair];
            if (reference_temperature > 0.0) {
                k_ij += parameters.k_ij_reciprocal_temperature_slope[pair]
                    * (Scalar(1.0) / temperature
                        - Scalar(1.0 / reference_temperature));
            }
            // EqID: cross_dispersion_energy
            // EqID: same_sign_ionic_dispersion
            const Scalar epsilon_i = dispersion_energy_override == nullptr
                ? Scalar(parameters.dispersion_energy_over_k[i])
                : (*dispersion_energy_override)[i];
            const Scalar epsilon_j = dispersion_energy_override == nullptr
                ? Scalar(parameters.dispersion_energy_over_k[j])
                : (*dispersion_energy_override)[j];
            const Scalar epsilon_ij = scalar_sqrt(epsilon_i * epsilon_j)
                * (Scalar(1.0) - k_ij);
            const Scalar weight = mole_fractions[i] * mole_fractions[j] *
                parameters.segment_count[i] * parameters.segment_count[j] *
                sigma_ij * sigma_ij * sigma_ij;
            // EqID: dispersion_first_moment
            first_moment += weight * epsilon_ij / temperature;
            // EqID: dispersion_second_moment
            second_moment += weight * epsilon_ij * epsilon_ij /
                (temperature * temperature);
        }
    }
    // EqID: dispersion_helmholtz
    const Scalar dispersion =
        -2.0 * PI * number_density * integral_1 * first_moment -
        PI * number_density * mean_segment_count * c1 * integral_2 * second_moment;

    // EqID: association_helmholtz
    Scalar association = 0.0;
    for (std::size_t site = 0; site < site_fraction.size(); ++site) {
        const std::size_t component =
            static_cast<std::size_t>(parameters.site_component[site]);
        association += mole_fractions[component] *
            parameters.site_multiplicity[site] *
            (scalar_log(site_fraction[site]) - site_fraction[site] / 2.0 + 0.5);
    }
    // EqID: residual_helmholtz
    return {hard_chain, dispersion, association, debye_huckel, born};
}

}  // namespace epcsaft
