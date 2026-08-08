#pragma once

#include "model.hpp"
#include "physical_constants.hpp"

#include <cppad/cppad.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <vector>

namespace epcsaft::detail {

inline constexpr double POLAR_PI = 3.141592653589793238462643383279502884;
inline constexpr double DQ_ALPHA = 1.19374;

class PolarDomainError final : public std::invalid_argument {
public:
    PolarDomainError() : std::invalid_argument("polar-domain-error") {}
};

inline constexpr std::array<double, 5> DD_A0 = {
    0.3043504, -0.1358588, 1.4493329, 0.3556977, -2.0653308,
};
inline constexpr std::array<double, 5> DD_A1 = {
    0.9534641, -1.8396383, 2.0131180, -7.3724958, 8.2374135,
};
inline constexpr std::array<double, 5> DD_A2 = {
    -1.1610080, 4.5258607, 0.9751222, -12.281038, 5.9397575,
};
inline constexpr std::array<double, 5> DD_B0 = {
    0.2187939, -1.1896431, 1.1626889, 0.0, 0.0,
};
inline constexpr std::array<double, 5> DD_B1 = {
    -0.5873164, 1.2489132, -0.5085280, 0.0, 0.0,
};
inline constexpr std::array<double, 5> DD_B2 = {
    3.4869576, -14.915974, 15.372022, 0.0, 0.0,
};
inline constexpr std::array<double, 5> DD_C0 = {
    -0.0646774, 0.1975882, -0.8087562, 0.6902849, 0.0,
};
inline constexpr std::array<double, 5> DD_C1 = {
    -0.9520876, 2.9924258, -2.3802636, -0.2701261, 0.0,
};
inline constexpr std::array<double, 5> DD_C2 = {
    -0.6260979, 1.2924686, 1.6542783, -3.4396744, 0.0,
};

inline constexpr std::array<double, 5> QQ_A0 = {
    1.2378308, 2.4355031, 1.6330905, -1.6118152, 6.9771185,
};
inline constexpr std::array<double, 5> QQ_A1 = {
    1.2854109, -11.465615, 22.086893, 7.4691383, -17.197772,
};
inline constexpr std::array<double, 5> QQ_A2 = {
    1.7942954, 0.7695103, 7.2647923, 94.486699, -77.148458,
};
inline constexpr std::array<double, 5> QQ_B0 = {
    0.4542718, -4.5016264, 3.5858868, 0.0, 0.0,
};
inline constexpr std::array<double, 5> QQ_B1 = {
    -0.8137340, 10.064030, -10.876631, 0.0, 0.0,
};
inline constexpr std::array<double, 5> QQ_B2 = {
    6.8682675, -5.1732238, -17.240207, 0.0, 0.0,
};
inline constexpr std::array<double, 5> QQ_C0 = {
    -0.5000437, 6.5318692, -16.014780, 14.425970, 0.0,
};
inline constexpr std::array<double, 5> QQ_C1 = {
    2.0002094, -6.7838658, 20.383246, -10.895984, 0.0,
};
inline constexpr std::array<double, 5> QQ_C2 = {
    3.1358271, 7.2475888, 3.0759478, 0.0, 0.0,
};

inline constexpr std::array<double, 5> DQ_A0 = {
    0.6970950, -0.6335541, 2.9455090, -1.4670273, 0.0,
};
inline constexpr std::array<double, 5> DQ_A1 = {
    -0.6734593, -1.4258991, 4.1944139, 1.0266216, 0.0,
};
inline constexpr std::array<double, 5> DQ_A2 = {
    0.6703408, -4.3384718, 7.2341684, 0.0, 0.0,
};
inline constexpr std::array<double, 5> DQ_B0 = {
    -0.4840383, 1.9704055, -2.1185727, 0.0, 0.0,
};
inline constexpr std::array<double, 5> DQ_B1 = {
    0.6765101, -3.0138675, 0.4674266, 0.0, 0.0,
};
inline constexpr std::array<double, 5> DQ_B2 = {
    -1.1675601, 2.1348843, 0.0, 0.0, 0.0,
};
inline constexpr std::array<double, 5> DQ_C0 = {
    7.846431, 33.42700, 4.689111, 0.0, 0.0,
};
inline constexpr std::array<double, 5> DQ_C1 = {
    -20.72202, -58.63904, -1.764887, 0.0, 0.0,
};

template <typename Scalar>
struct PolarTerms final {
    Scalar dd_second;
    Scalar dd_third;
    Scalar qq_second;
    Scalar qq_third;
    Scalar dq_second;
    Scalar dq_third;
};

template <typename Scalar>
// EqExclude: scalar utility dispatch
double polar_primal(const Scalar& value) {
    if constexpr (std::is_same_v<Scalar, double>) {
        return value;
    } else {
        return CppAD::Value(CppAD::Var2Par(value));
    }
}

template <typename Scalar>
// EqExclude: scalar utility dispatch
Scalar polar_sqrt(const Scalar& value) {
    if constexpr (std::is_same_v<Scalar, double>) {
        return std::sqrt(value);
    } else {
        return CppAD::sqrt(value);
    }
}

template <typename Scalar>
// EqExclude: scalar utility dispatch
Scalar polar_cbrt(const Scalar& value) {
    if constexpr (std::is_same_v<Scalar, double>) {
        return std::cbrt(value);
    } else {
        return CppAD::pow(value, 1.0 / 3.0);
    }
}

template <typename Scalar>
// EqExclude: scalar utility dispatch
Scalar polar_power(const Scalar& value, std::size_t exponent) {
    Scalar result = Scalar(1.0);
    for (std::size_t index = 0; index < exponent; ++index) {
        result *= value;
    }
    return result;
}

template <typename Scalar>
// EqID: polar_chain_coefficient
Scalar chain_coefficient(
    const Scalar& chain,
    double zero,
    double one,
    double two
) {
    const Scalar first = (chain - 1.0) / chain;
    return zero + first * one + first * (chain - 2.0) / chain * two;
}

template <typename Scalar>
// EqExclude: scalar utility dispatch
Scalar polar_polynomial(const std::array<Scalar, 5>& coefficients, const Scalar& eta) {
    Scalar result = coefficients.back();
    for (std::size_t index = coefficients.size() - 1; index > 0; --index) {
        result = result * eta + coefficients[index - 1];
    }
    return result;
}

template <typename Scalar>
// EqID: polar_pair_integral
Scalar polar_j2(
    const std::array<double, 5>& a0,
    const std::array<double, 5>& a1,
    const std::array<double, 5>& a2,
    const std::array<double, 5>& b0,
    const std::array<double, 5>& b1,
    const std::array<double, 5>& b2,
    const Scalar& chain,
    const Scalar& epsilon_over_temperature,
    const Scalar& eta
) {
    std::array<Scalar, 5> coefficients{};
    for (std::size_t n = 0; n < coefficients.size(); ++n) {
        coefficients[n] = chain_coefficient(chain, a0[n], a1[n], a2[n])
            + chain_coefficient(chain, b0[n], b1[n], b2[n])
                * epsilon_over_temperature;
    }
    return polar_polynomial(coefficients, eta);
}

template <typename Scalar>
// EqID: polar_triplet_integral
Scalar polar_j3(
    const std::array<double, 5>& c0,
    const std::array<double, 5>& c1,
    const std::array<double, 5>& c2,
    const Scalar& chain,
    const Scalar& eta
) {
    std::array<Scalar, 5> coefficients{};
    for (std::size_t n = 0; n < coefficients.size(); ++n) {
        coefficients[n] = chain_coefficient(chain, c0[n], c1[n], c2[n]);
    }
    return polar_polynomial(coefficients, eta);
}

template <typename Scalar>
// EqID: polar_triplet_integral
Scalar dq_j3(const Scalar& chain, const Scalar& eta) {
    std::array<Scalar, 5> coefficients{};
    for (std::size_t n = 0; n < coefficients.size(); ++n) {
        coefficients[n] = DQ_C0[n] + (chain - 1.0) / chain * DQ_C1[n];
    }
    return polar_polynomial(coefficients, eta);
}

template <typename Scalar>
// EqID: polar_reduced_dipole_moment
Scalar reduced_dipole_squared(
    double moment_debye,
    const Scalar& segments,
    const Scalar& epsilon_over_k,
    const Scalar& sigma_angstrom
) {
    const double moment = moment_debye * PHYSICAL_DEBYE_C_M;
    const Scalar sigma_m = sigma_angstrom * PHYSICAL_ANGSTROM_M;
    return moment * moment /
        (4.0 * POLAR_PI * PHYSICAL_VACUUM_PERMITTIVITY * segments
         * PHYSICAL_BOLTZMANN * epsilon_over_k * polar_power(sigma_m, 3));
}

template <typename Scalar>
// EqID: polar_reduced_quadrupole_moment
Scalar reduced_quadrupole_squared(
    double moment_debye_angstrom,
    const Scalar& segments,
    const Scalar& epsilon_over_k,
    const Scalar& sigma_angstrom
) {
    const double moment = moment_debye_angstrom * PHYSICAL_DEBYE_C_M
        * PHYSICAL_ANGSTROM_M;
    const Scalar sigma_m = sigma_angstrom * PHYSICAL_ANGSTROM_M;
    return moment * moment /
        (4.0 * POLAR_PI * PHYSICAL_VACUUM_PERMITTIVITY * segments
         * PHYSICAL_BOLTZMANN * epsilon_over_k * polar_power(sigma_m, 5));
}

template <typename Scalar, typename Parameters>
// EqExclude: scalar utility dispatch
Scalar polar_epsilon(
    std::size_t component,
    const Parameters& parameters,
    const std::vector<Scalar>* override
) {
    return override == nullptr
        ? Scalar(parameters.dispersion_energy_over_k[component])
        : (*override)[component];
}

template <typename Scalar, typename Parameters>
// EqID: polar_cross_energy
Scalar polar_cross_epsilon(
    std::size_t i,
    std::size_t j,
    const Parameters& parameters,
    const std::vector<Scalar>* energy_override
) {
    const Scalar k =
        Scalar(parameters.k_ij[i * parameters.segment_count.size() + j]);
    return (Scalar(1.0) - k)
        * polar_sqrt(
            polar_epsilon<Scalar>(i, parameters, energy_override)
            * polar_epsilon<Scalar>(j, parameters, energy_override)
        );
}

template <typename Scalar>
// EqExclude: input and domain validation
bool polar_present(const std::vector<Scalar>& fractions, const std::vector<double>& moments) {
    for (std::size_t i = 0; i < moments.size(); ++i) {
        if (moments[i] > 0.0 && polar_primal(fractions[i]) > 0.0) {
            return true;
        }
    }
    return false;
}

template <typename Scalar, typename Parameters>
// EqExclude: scalar utility dispatch
PolarTerms<Scalar> polar_terms(
    const Scalar& temperature,
    const Scalar& number_density,
    const std::vector<Scalar>& fractions,
    const Parameters& parameters,
    const std::vector<Scalar>& sigma,
    const Scalar& eta,
    const std::vector<Scalar>* energy_override = nullptr
) {
    PolarTerms<Scalar> result{
        Scalar(0.0), Scalar(0.0), Scalar(0.0),
        Scalar(0.0), Scalar(0.0), Scalar(0.0),
    };
    if (parameters.polar_formulation == POLAR_FORMULATION_NONE) {
        return result;
    }

    const std::size_t size = fractions.size();
    std::vector<Scalar> epsilon(size);
    std::vector<Scalar> mu2(size, Scalar(0.0));
    std::vector<Scalar> q2(size, Scalar(0.0));
    for (std::size_t i = 0; i < size; ++i) {
        epsilon[i] = polar_epsilon<Scalar>(i, parameters, energy_override);
        if (parameters.dipole_moment_debye[i] > 0.0) {
            mu2[i] = reduced_dipole_squared<Scalar>(
                parameters.dipole_moment_debye[i], Scalar(parameters.segment_count[i]),
                epsilon[i], sigma[i]
            );
        }
        if (parameters.quadrupole_moment_debye_angstrom[i] > 0.0) {
            q2[i] = reduced_quadrupole_squared<Scalar>(
                parameters.quadrupole_moment_debye_angstrom[i],
                Scalar(parameters.segment_count[i]), epsilon[i], sigma[i]
            );
        }
    }

    const bool dd_active = polar_present(fractions, parameters.dipole_moment_debye);
    const bool qq_active = polar_present(
        fractions, parameters.quadrupole_moment_debye_angstrom
    );

    if (dd_active) {
        Scalar pair_sum = Scalar(0.0);
        Scalar triplet_sum = Scalar(0.0);
        for (std::size_t i = 0; i < size; ++i) {
            if (parameters.dipole_moment_debye[i] == 0.0) continue;
            for (std::size_t j = 0; j < size; ++j) {
                if (parameters.dipole_moment_debye[j] == 0.0) continue;
                const Scalar sigma_ij = Scalar(0.5) * (sigma[i] + sigma[j]);
                const Scalar pair_chain = polar_sqrt(
                    parameters.segment_count[i] * parameters.segment_count[j]
                );
                const Scalar chain = polar_primal(pair_chain) > 2.0
                    ? Scalar(2.0) : pair_chain;
                const Scalar j2 = polar_j2(
                    DD_A0, DD_A1, DD_A2, DD_B0, DD_B1, DD_B2, chain,
                    polar_cross_epsilon(i, j, parameters, energy_override)
                        / temperature,
                    eta
                );
                pair_sum += fractions[i] * fractions[j]
                    * epsilon[i] * epsilon[j] / (temperature * temperature)
                    * polar_power(sigma[i], 3) * polar_power(sigma[j], 3)
                    / polar_power(sigma_ij, 3) * mu2[i] * mu2[j] * j2;
                for (std::size_t k = 0; k < size; ++k) {
                    if (parameters.dipole_moment_debye[k] == 0.0) continue;
                    const Scalar sigma_ik = Scalar(0.5) * (sigma[i] + sigma[k]);
                    const Scalar sigma_jk = Scalar(0.5) * (sigma[j] + sigma[k]);
                    const Scalar raw_triplet_chain = polar_cbrt(
                        parameters.segment_count[i] * parameters.segment_count[j]
                        * parameters.segment_count[k]
                    );
                    const Scalar triplet_chain = polar_primal(raw_triplet_chain) > 2.0
                        ? Scalar(2.0) : raw_triplet_chain;
                    const Scalar j3 = polar_j3(
                        DD_C0, DD_C1, DD_C2, triplet_chain, eta
                    );
                    triplet_sum += fractions[i] * fractions[j] * fractions[k]
                        * epsilon[i] * epsilon[j] * epsilon[k]
                        / polar_power(temperature, 3)
                        * polar_power(sigma[i], 3) * polar_power(sigma[j], 3)
                        * polar_power(sigma[k], 3) / (sigma_ij * sigma_ik * sigma_jk)
                        * mu2[i] * mu2[j] * mu2[k] * j3;
                }
            }
        }
        // EqID: polar_dd_second_order
        result.dd_second = -POLAR_PI * number_density * pair_sum;
        // EqID: polar_dd_third_order
        result.dd_third = -4.0 * POLAR_PI * POLAR_PI / 3.0
            * number_density * number_density * triplet_sum;
    }

    if (qq_active) {
        Scalar pair_sum = Scalar(0.0);
        Scalar triplet_sum = Scalar(0.0);
        for (std::size_t i = 0; i < size; ++i) {
            if (parameters.quadrupole_moment_debye_angstrom[i] == 0.0) continue;
            for (std::size_t j = 0; j < size; ++j) {
                if (parameters.quadrupole_moment_debye_angstrom[j] == 0.0) continue;
                const Scalar sigma_ij = Scalar(0.5) * (sigma[i] + sigma[j]);
                const Scalar chain = polar_sqrt(
                    parameters.segment_count[i] * parameters.segment_count[j]
                );
                const Scalar j2 = polar_j2(
                    QQ_A0, QQ_A1, QQ_A2, QQ_B0, QQ_B1, QQ_B2, chain,
                    polar_cross_epsilon(i, j, parameters, energy_override)
                        / temperature,
                    eta
                );
                pair_sum += fractions[i] * fractions[j]
                    * epsilon[i] * epsilon[j] / (temperature * temperature)
                    * polar_power(sigma[i], 5) * polar_power(sigma[j], 5)
                    / polar_power(sigma_ij, 7) * q2[i] * q2[j] * j2;
                for (std::size_t k = 0; k < size; ++k) {
                    if (parameters.quadrupole_moment_debye_angstrom[k] == 0.0) continue;
                    const Scalar sigma_ik = Scalar(0.5) * (sigma[i] + sigma[k]);
                    const Scalar sigma_jk = Scalar(0.5) * (sigma[j] + sigma[k]);
                    const Scalar triplet_chain = polar_cbrt(
                        parameters.segment_count[i] * parameters.segment_count[j]
                        * parameters.segment_count[k]
                    );
                    const Scalar j3 = polar_j3(
                        QQ_C0, QQ_C1, QQ_C2, triplet_chain, eta
                    );
                    triplet_sum += fractions[i] * fractions[j] * fractions[k]
                        * epsilon[i] * epsilon[j] * epsilon[k]
                        / polar_power(temperature, 3)
                        * polar_power(sigma[i], 5) * polar_power(sigma[j], 5)
                        * polar_power(sigma[k], 5)
                        / (polar_power(sigma_ij, 3) * polar_power(sigma_ik, 3)
                           * polar_power(sigma_jk, 3))
                        * q2[i] * q2[j] * q2[k] * j3;
                }
            }
        }
        // EqID: polar_qq_second_order
        result.qq_second = -9.0 * POLAR_PI / 16.0 * number_density * pair_sum;
        // EqID: polar_qq_third_order
        result.qq_third = 9.0 * POLAR_PI * POLAR_PI / 16.0
            * number_density * number_density * triplet_sum;
    }

    if (dd_active && qq_active) {
        Scalar pair_sum = Scalar(0.0);
        Scalar triplet_sum = Scalar(0.0);
        for (std::size_t i = 0; i < size; ++i) {
            for (std::size_t j = 0; j < size; ++j) {
                const Scalar sigma_ij = Scalar(0.5) * (sigma[i] + sigma[j]);
                if (parameters.dipole_moment_debye[i] > 0.0
                    && parameters.quadrupole_moment_debye_angstrom[j] > 0.0) {
                    const Scalar chain = polar_sqrt(
                        parameters.segment_count[i] * parameters.segment_count[j]
                    );
                    const Scalar j2 = polar_j2(
                        DQ_A0, DQ_A1, DQ_A2, DQ_B0, DQ_B1, DQ_B2, chain,
                        polar_cross_epsilon(i, j, parameters, energy_override)
                            / temperature,
                        eta
                    );
                    pair_sum += fractions[i] * fractions[j]
                        * epsilon[i] * epsilon[j] / (temperature * temperature)
                        * polar_power(sigma[i], 3) * polar_power(sigma[j], 5)
                        / polar_power(sigma_ij, 5) * mu2[i] * q2[j] * j2;
                }
                for (std::size_t k = 0; k < size; ++k) {
                    const Scalar moment_product = mu2[i] * mu2[j] * q2[k]
                        + DQ_ALPHA * mu2[i] * q2[j] * q2[k];
                    if (polar_primal(moment_product) == 0.0) continue;
                    const Scalar sigma_ik = Scalar(0.5) * (sigma[i] + sigma[k]);
                    const Scalar sigma_jk = Scalar(0.5) * (sigma[j] + sigma[k]);
                    const Scalar chain = polar_cbrt(
                        parameters.segment_count[i] * parameters.segment_count[j]
                        * parameters.segment_count[k]
                    );
                    triplet_sum += fractions[i] * fractions[j] * fractions[k]
                        * epsilon[i] * epsilon[j] * epsilon[k]
                        / polar_power(temperature, 3)
                        * polar_power(sigma[i], 4) * polar_power(sigma[j], 4)
                        * polar_power(sigma[k], 4)
                        / (polar_power(sigma_ij, 2) * polar_power(sigma_ik, 2)
                           * polar_power(sigma_jk, 2))
                        * moment_product * dq_j3(chain, eta);
                }
            }
        }
        // EqID: polar_dq_second_order
        result.dq_second = -9.0 * POLAR_PI / 4.0 * number_density * pair_sum;
        // EqID: polar_dq_third_order
        result.dq_third = -number_density * number_density * triplet_sum;
    }
    return result;
}

template <typename Scalar>
// EqID: polar_pade
Scalar polar_pade(const Scalar& second, const Scalar& third, bool active) {
    if (!active) {
        return Scalar(0.0);
    }
    const Scalar denominator = second - third;
    const double second_value = polar_primal(second);
    const double third_value = polar_primal(third);
    const double denominator_value = polar_primal(denominator);
    const double scale = std::max(std::abs(second_value), std::abs(third_value));
    if (!std::isfinite(second_value) || !std::isfinite(third_value)
        || !std::isfinite(denominator_value)
        || std::abs(denominator_value)
            <= 64.0 * std::numeric_limits<double>::epsilon() * scale) {
        throw PolarDomainError();
    }
    return second * second / denominator;
}

template <typename Scalar, typename Parameters>
// EqID: polar_helmholtz
Scalar polar_helmholtz(
    const PolarTerms<Scalar>& terms,
    const std::vector<Scalar>& fractions,
    const Parameters& parameters
) {
    const bool dd_active = polar_present(fractions, parameters.dipole_moment_debye);
    const bool qq_active = polar_present(
        fractions, parameters.quadrupole_moment_debye_angstrom
    );
    return polar_pade(terms.dd_second, terms.dd_third, dd_active)
        + polar_pade(terms.qq_second, terms.qq_third, qq_active)
        + polar_pade(terms.dq_second, terms.dq_third, dd_active && qq_active);
}

}  // namespace epcsaft::detail
