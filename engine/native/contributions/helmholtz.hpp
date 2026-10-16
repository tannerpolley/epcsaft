#pragma once

#include "core/numerical_policy.hpp"
#include "core/parameters.hpp"
#include "core/physical_constants.hpp"
#include "core/scalar.hpp"

#include <cppad/cppad.hpp>

#include <array>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace epcsaft {

// EqID: phase_helmholtz_amount_volume
// The mechanical ideal term is deliberately independent of the residual tape.
// Its caller may use an off-normalized composition chart; no normalization is
// implied here.
template <typename Scalar>
Scalar ideal_mechanical_helmholtz(
    const std::vector<Scalar>& amounts, const Scalar& volume, const Scalar& rho0,
    bool clamp_subnormal = true
) {
    if constexpr (std::is_same_v<Scalar, double>) {
        if (!(volume > 0.0) || !std::isfinite(volume)
            || !(rho0 > 0.0) || !std::isfinite(rho0))
            throw std::invalid_argument("invalid ideal mechanical scale");
    }
    Scalar result = Scalar(0.0);
    const Scalar log_scale = scalar_log(volume) + scalar_log(rho0);
    for (const Scalar& amount : amounts) {
        if constexpr (std::is_same_v<Scalar, double>) {
            if (!std::isfinite(amount) || amount < 0.0)
                throw std::invalid_argument("invalid ideal amount");
            if (amount > 0.0) {
                const double log_amount =
                    std::log(clamp_subnormal
                                 ? std::max(amount, std::numeric_limits<double>::min())
                                 : amount);
                result += amount * (log_amount - log_scale - 1.0);
            }
        } else {
            const Scalar log_amount = scalar_log(
                clamp_subnormal
                    ? scalar_select_less(
                        amount, Scalar(std::numeric_limits<double>::min()),
                        Scalar(std::numeric_limits<double>::min()), amount)
                    : amount);
            result += scalar_select_less(
                Scalar(0.0), amount,
                amount * (log_amount - log_scale - Scalar(1.0)),
                Scalar(0.0));
        }
    }
    return result;
}

// Gross and Sadowski (2001) Table 1; FeOs carries extra digits (#222).
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

// Gross and Vrabec (2006) Table 1, rows n = 0..4 of the dipolar J integrals.
inline constexpr std::array<std::array<double, 5>, 3> DIPOLE_A = {{
    {0.3043504, -0.1358588, 1.4493329, 0.3556977, -2.0653308},
    {0.9534641, -1.8396383, 2.0131180, -7.3724958, 8.2374135},
    {-1.1610080, 4.5258607, 0.9751222, -12.281038, 5.9397575},
}};
inline constexpr std::array<std::array<double, 5>, 3> DIPOLE_B = {{
    {0.2187939, -1.1896431, 1.1626889, 0.0, 0.0},
    {-0.5873164, 1.2489132, -0.5085280, 0.0, 0.0},
    {3.4869576, -14.915974, 15.372022, 0.0, 0.0},
}};
inline constexpr std::array<std::array<double, 5>, 3> DIPOLE_C = {{
    {-0.0646774, 0.1975882, -0.8087562, 0.6902849, 0.0},
    {-0.9520876, 2.9924258, -2.3802636, -0.2701261, 0.0},
    {-0.6260979, 1.2924686, 1.6542783, -3.4396744, 0.0},
}};

template <typename Scalar>
struct HelmholtzResult final {
    Scalar hard_chain;
    Scalar dispersion;
    Scalar association;
    Scalar debye_huckel;
    Scalar born;
    Scalar dipolar;
};

// Fixed: the scalar of quantities that do not depend on density (temperature and composition), so a
// density series carries only density-dependent arithmetic; the tape and double instantiations take Fixed = Scalar.
template <typename Scalar, typename Fixed>
// EqID: dispersion_polynomial
Scalar polynomial(const std::array<Fixed, 7>& coefficients, const Scalar& eta) {
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
    const std::vector<Scalar>& sigma
) {
    std::vector<Scalar> diameter(sigma.size());
    for (std::size_t component = 0; component < sigma.size(); ++component) {
        diameter[component] = parameters.charge_number[component] == 0
            ? sigma[component] *
                (Scalar(1.0) - Scalar(0.12) * scalar_exp(
                    Scalar(-3.0)
                        * Scalar(parameters.dispersion_energy_over_k[component])
                        / temperature
                ))
            : Scalar(parameters.packing_diameter_angstrom[component]);
    }
    return diameter;
}

template <typename Scalar, typename ParameterSet, typename Fixed>
// EqID: packing_fraction
Scalar packing_fraction_from_geometry(
    const Scalar& molar_density,
    const std::vector<Fixed>& mole_fractions,
    const ParameterSet& parameters,
    const std::vector<Fixed>& diameter
) {
    Fixed moment = Fixed(0.0);
    for (std::size_t component = 0; component < mole_fractions.size(); ++component) {
        moment += mole_fractions[component]
            * parameters.segment_count[component]
            * diameter[component] * diameter[component] * diameter[component];
    }
    return (PI / 6.0) * molar_density * NUMBER_DENSITY_PER_ANGSTROM3 * moment;
}

// EqID: packing_moments
// zeta_3 comes from packing_fraction_from_geometry, the owner of the density bound.
template <typename Scalar, typename ParameterSet, typename Fixed>
std::array<Scalar, 4> packing_moments(
    const Scalar& molar_density,
    const std::vector<Fixed>& mole_fractions,
    const ParameterSet& parameters,
    const std::vector<Fixed>& diameter
) {
    const Scalar number_density = molar_density * NUMBER_DENSITY_PER_ANGSTROM3;
    std::array<Scalar, 4> zeta;
    for (std::size_t order = 0; order < 3; ++order) {
        Fixed moment = Fixed(0.0);
        for (std::size_t i = 0; i < mole_fractions.size(); ++i) {
            Fixed power = Fixed(1.0);
            for (std::size_t exponent = 0; exponent < order; ++exponent) power *= diameter[i];
            moment += mole_fractions[i] * parameters.segment_count[i] * power;
        }
        zeta[order] = (PI / 6.0) * number_density * moment;
    }
    zeta[3] = packing_fraction_from_geometry(molar_density, mole_fractions, parameters, diameter);
    return zeta;
}

// EqID: hard_sphere_contact
// Contact value at reduced diameter d_i d_j/(d_i + d_j); d_i/2 for a like pair.
template <typename Scalar, typename Fixed>
Scalar hard_sphere_contact(const Fixed& reduced_diameter, const Scalar& zeta_2, const Scalar& one_minus_eta) {
    return Scalar(1.0) / one_minus_eta
        + 3.0 * reduced_diameter * zeta_2 / (one_minus_eta * one_minus_eta)
        + 2.0 * reduced_diameter * reduced_diameter * zeta_2 * zeta_2
            / (one_minus_eta * one_minus_eta * one_minus_eta);
}

template <typename Scalar, typename Fixed>
// EqID: dipolar_integrals
// Coefficient a_n(m) = t0 + (m-1)/m t1 + (m-1)(m-2)/m^2 t2 with m capped at 2.
Scalar dipolar_integral(
    const std::array<std::array<double, 5>, 3>& table,
    const std::array<std::array<double, 5>, 3>* energy_table,
    const Fixed& chain, const Fixed& reduced_energy, const Scalar& eta
) {
    const Fixed m = scalar_select_less(chain, Fixed(2.0), chain, Fixed(2.0));
    const Fixed first = (m - 1.0) / m;
    const Fixed second = first * (m - 2.0) / m;
    Scalar value = Scalar(0.0);
    for (std::size_t n = table[0].size(); n-- > 0;) {
        Fixed coefficient = table[0][n] + first * table[1][n] + second * table[2][n];
        if (energy_table != nullptr) {
            const auto& b = *energy_table;
            coefficient += (b[0][n] + first * b[1][n] + second * b[2][n]) * reduced_energy;
        }
        value = value * eta + coefficient;
    }
    return value;
}

template <typename Scalar, typename ParameterSet, typename Fixed>
// EqID: dipolar_helmholtz
// Gross and Vrabec (2006) Eqs. 7-9 with q_i = n_i mu_i^2 / (4 pi eps0 k T m_i sigma_i^3).
// The unlike energy in J2 is sqrt(eps_i eps_j) without k_ij (source p. 1196).
Scalar dipolar_helmholtz(
    const Fixed& temperature, const Scalar& number_density, const Scalar& eta,
    const std::vector<Fixed>& x, const ParameterSet& parameters,
    const std::vector<Fixed>& sigma, const std::vector<Fixed>& epsilon
) {
    const std::size_t size = x.size();
    constexpr double scale = DEBYE_SQUARED_J_ANGSTROM3 / BOLTZMANN_J_PER_K;
    std::vector<Fixed> q(size, Fixed(0.0));
    std::vector<std::size_t> polar;
    for (std::size_t i = 0; i < size; ++i) {
        const double mu = parameters.dipole_moment_debye[i];
        if (mu == 0.0 || parameters.dipole_count[i] == 0.0) continue;
        polar.push_back(i);
        q[i] = parameters.dipole_count[i] * mu * mu * scale
            / (temperature * parameters.segment_count[i] * sigma[i] * sigma[i] * sigma[i]);
    }
    Scalar second = Scalar(0.0);
    Scalar third = Scalar(0.0);
    for (const std::size_t i : polar) {
        const Fixed volume_i = sigma[i] * sigma[i] * sigma[i];
        for (const std::size_t j : polar) {
            const Fixed volume_j = sigma[j] * sigma[j] * sigma[j];
            const Fixed sigma_ij = 0.5 * (sigma[i] + sigma[j]);
            const Scalar j2 = dipolar_integral<Scalar, Fixed>(
                DIPOLE_A, &DIPOLE_B,
                scalar_sqrt(parameters.segment_count[i] * parameters.segment_count[j]),
                scalar_sqrt(epsilon[i] * epsilon[j]) / temperature, eta);
            second += x[i] * x[j] * volume_i * volume_j
                / (sigma_ij * sigma_ij * sigma_ij) * q[i] * q[j] * j2;
            for (const std::size_t k : polar) {
                const Scalar j3 = dipolar_integral<Scalar, Fixed>(
                    DIPOLE_C, nullptr,
                    scalar_exp(scalar_log(parameters.segment_count[i]
                        * parameters.segment_count[j] * parameters.segment_count[k]) / 3.0),
                    Fixed(0.0), eta);
                third += x[i] * x[j] * x[k] * volume_i * volume_j
                    * sigma[k] * sigma[k] * sigma[k]
                    / (sigma_ij * 0.5 * (sigma[i] + sigma[k]) * 0.5 * (sigma[j] + sigma[k]))
                    * q[i] * q[j] * q[k] * j3;
            }
        }
    }
    // EqID: dipolar_second_order
    const Scalar a2 = -PI * number_density * second;
    // EqID: dipolar_third_order
    const Scalar a3 = -4.0 * PI * PI / 3.0 * number_density * number_density * third;
    // a_dd = A2 / (1 - A3/A2). With every dipole absent A2 = A3 = 0 and the Pade factor
    // is taken as its limit 1; both selected operands avoid 0/0 on the tape.
    return a2 * scalar_select_less(Scalar(0.0), a2 * a2, a2, Scalar(1.0))
        / scalar_select_less(Scalar(0.0), a2 * a2, a2 - a3, Scalar(1.0));
}

// EqID: pair_parameter_temperature_law
// q(T) = a_0 + a_1 (T - T_r) + a_2 (T - T_r)^2 + a_3 (T - T_r)^3 + s (1/T - 1/T_r) of every declared law.
template <typename Fixed, typename ParameterSet>
std::vector<Fixed> temperature_law_values(const Fixed& temperature, const ParameterSet& parameters) {
    std::vector<Fixed> values(parameters.temperature_law_reference_k.size(), Fixed(0.0));
    for (std::size_t term = 0; term < parameters.temperature_law_term_law.size(); ++term) {
        const auto law = static_cast<std::size_t>(parameters.temperature_law_term_law[term]);
        const double reference = parameters.temperature_law_reference_k[law];
        const int power = parameters.temperature_law_term_power[term];
        Fixed basis = power < 0 ? Fixed(1.0) / temperature - Fixed(1.0 / reference) : Fixed(1.0);
        for (int p = 0; p < power; ++p) basis *= temperature - Fixed(reference);
        values[law] += parameters.temperature_law_coefficient[term] * basis;
    }
    return values;
}

inline std::string law_domain_failure(const char* quantity, std::size_t slot, std::size_t size, double temperature) {
    return std::string("invalid evaluated ") + quantity + " for pair " + std::to_string(slot / size) + "/"
        + std::to_string(slot % size) + " at T = " + std::to_string(temperature) + " K";
}

// The law value owning `slot` (row-major over `size`, smaller index first). An evaluated value outside its
// quantity's domain is refused with its target and temperature, never clamped: finite k_ij, l_ij < 1
// (positive unlike diameter), nonnegative association energy and volume.
template <typename Fixed, typename ParameterSet>
std::optional<Fixed> temperature_law_value(const Fixed& temperature, const ParameterSet& parameters,
    const std::vector<Fixed>& laws, TemperatureLaw target, std::size_t slot, std::size_t size) {
    for (const auto& [owner, owned, law] : parameters.temperature_law_slot) {
        if (owner != static_cast<int>(target) || static_cast<std::size_t>(owned) != slot) continue;
        const Fixed& value = laws[static_cast<std::size_t>(law)];
        if constexpr (std::is_same_v<Fixed, double>) {
            if (!std::isfinite(value) || (target == TemperatureLaw::Lij && !(value < 1.0))
                || (target >= TemperatureLaw::AssociationEnergyOverK && !(value >= 0.0)))
                throw std::invalid_argument(law_domain_failure(
                    target == TemperatureLaw::Kij ? "k_ij" : target == TemperatureLaw::Lij ? "l_ij"
                    : target == TemperatureLaw::AssociationVolume ? "association_volume" : "association_energy_over_k",
                    slot, size, temperature));
        }
        return value;
    }
    return std::nullopt;
}

// k_ij or l_ij of components (i, j) at this temperature: its law, else the stored constant with the
// legacy reciprocal-temperature k_ij term in its original operation order.
template <typename Fixed, typename ParameterSet>
Fixed pair_correction(const Fixed& temperature, const ParameterSet& parameters, const std::vector<Fixed>& laws,
    TemperatureLaw target, std::size_t i, std::size_t j) {
    const std::size_t size = parameters.segment_count.size(), pair = i * size + j;
    if (!laws.empty())
        if (auto value = temperature_law_value(temperature, parameters, laws, target, std::min(i, j) * size + std::max(i, j), size))
            return *value;
    if (target == TemperatureLaw::Lij) return parameters.l_ij[pair];
    const double reference_temperature = parameters.k_ij_reference_temperature_k[pair];
    Fixed k_ij = parameters.k_ij[pair];
    if (reference_temperature > 0.0) {
        k_ij += parameters.k_ij_reciprocal_temperature_slope[pair]
            * (Fixed(1.0) / temperature - Fixed(1.0 / reference_temperature));
    }
    return k_ij;
}

template <typename Scalar, typename ParameterSet, typename Fixed>
HelmholtzResult<Scalar> residual_helmholtz(
    const Fixed& temperature,
    const Scalar& molar_density,
    const std::vector<Fixed>& mole_fractions,
    const ParameterSet& parameters,
    const std::vector<Fixed>& sigma,
    const std::vector<Fixed>& diameter,
    const std::vector<Scalar>& site_fraction,
    const Scalar& debye_huckel,
    const Scalar& born
) {
    const std::size_t size = parameters.segment_count.size();
    // EqID: angstrom_number_density
    const Scalar number_density = molar_density * NUMBER_DENSITY_PER_ANGSTROM3;
    // EqID: mean_segment_count
    Fixed mean_segment_count = Fixed(0.0);
    for (std::size_t i = 0; i < size; ++i) {
        mean_segment_count += mole_fractions[i] * parameters.segment_count[i];
    }

    const auto zeta = packing_moments(molar_density, mole_fractions, parameters, diameter);
    if constexpr (std::is_same_v<Scalar, double>) {
        // EqID: close_packing_fraction
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
        hard_chain -= mole_fractions[i] * (parameters.segment_count[i] - 1.0) *
            scalar_log(hard_sphere_contact(diameter[i] / Fixed(2.0), zeta[2], one_minus_eta));
    }

    // EqID: dispersion_chain_coefficients
    const Fixed chain_factor_1 =
        (mean_segment_count - 1.0) / mean_segment_count;
    const Fixed chain_factor_2 = chain_factor_1 *
        (mean_segment_count - 2.0) / mean_segment_count;
    std::array<Fixed, 7> a{};
    std::array<Fixed, 7> b{};
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

    Fixed first_moment = Fixed(0.0);
    Fixed second_moment = Fixed(0.0);
    const auto laws = temperature_law_values(temperature, parameters);
    for (std::size_t i = 0; i < size; ++i) {
        for (std::size_t j = 0; j < size; ++j) {
            // EqID: cross_segment_diameter
            const Fixed sigma_ij = Fixed(0.5) * (sigma[i] + sigma[j])
                * (Fixed(1.0) - pair_correction(temperature, parameters, laws, TemperatureLaw::Lij, i, j));
            const Fixed k_ij = pair_correction(temperature, parameters, laws, TemperatureLaw::Kij, i, j);
            // EqID: cross_dispersion_energy
            // EqID: same_sign_ionic_dispersion
            const int charge_i = parameters.charge_number[i];
            const int charge_j = parameters.charge_number[j];
            const bool same_sign_excluded =
                parameters.ionic_dispersion.is(IonicDispersion::ExcludeSameSign)
                && charge_i * charge_j > 0;
            const Fixed epsilon_i = Fixed(parameters.dispersion_energy_over_k[i]);
            const Fixed epsilon_j = Fixed(parameters.dispersion_energy_over_k[j]);
            const Fixed epsilon_ij = same_sign_excluded
                ? Fixed(0.0)
                : scalar_sqrt(epsilon_i * epsilon_j) * (Fixed(1.0) - k_ij);
            const Fixed weight = mole_fractions[i] * mole_fractions[j] *
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

    // Optional: an undeclared polar family leaves the dipole arrays empty.
    Scalar dipolar = Scalar(0.0);
    if (!parameters.dipole_moment_debye.empty()) {
        const std::vector<Fixed> epsilon(
            parameters.dispersion_energy_over_k.begin(), parameters.dispersion_energy_over_k.end());
        dipolar = dipolar_helmholtz(
            temperature, number_density, zeta[3], mole_fractions, parameters, sigma, epsilon);
    }

    // EqID: association_helmholtz
    Scalar association = Scalar(0.0);
    for (std::size_t site = 0; site < site_fraction.size(); ++site) {
        const std::size_t component =
            static_cast<std::size_t>(parameters.site_component[site]);
        association += mole_fractions[component] *
            parameters.site_multiplicity[site] *
            (scalar_log(site_fraction[site]) - site_fraction[site] / 2.0 + 0.5);
    }
    // EqID: residual_helmholtz
    return {hard_chain, dispersion, association, debye_huckel, born, dipolar};
}

}  // namespace epcsaft
