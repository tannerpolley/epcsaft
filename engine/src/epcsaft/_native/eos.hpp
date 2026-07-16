#pragma once

#include <cppad/cppad.hpp>

#include <array>
#include <cmath>
#include <cstddef>
#include <type_traits>
#include <vector>

namespace epcsaft {

struct NativeResult final {
    double hard_chain;
    double dispersion;
    double residual_helmholtz;
    double compressibility_factor;
    double pressure;
};

class NativeModel final {
public:
    NativeModel(
        std::vector<double> m,
        std::vector<double> sigma,
        std::vector<double> epsilon_over_k,
        std::vector<double> kij
    );
    NativeResult evaluate(
        double temperature,
        double molar_density,
        const std::vector<double>& mole_fractions
    ) const;

private:
    const std::vector<double> m_;
    const std::vector<double> sigma_;
    const std::vector<double> epsilon_over_k_;
    const std::vector<double> kij_;
};

namespace detail {

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
};

template <typename Scalar>
Scalar scalar_log(const Scalar& value) {
    if constexpr (std::is_same_v<Scalar, double>) {
        return std::log(value);
    } else {
        return CppAD::log(value);
    }
}

template <typename Scalar>
Scalar polynomial(const std::array<double, 7>& coefficients, const Scalar& eta) {
    Scalar value = coefficients.back();
    for (std::size_t index = coefficients.size() - 1; index > 0; --index) {
        value = value * eta + coefficients[index - 1];
    }
    return value;
}

template <typename Scalar>
HelmholtzResult<Scalar> residual_helmholtz(
    double temperature,
    const Scalar& molar_density,
    const std::vector<double>& mole_fractions,
    const std::vector<double>& m,
    const std::vector<double>& sigma,
    const std::vector<double>& epsilon_over_k,
    const std::vector<double>& kij
) {
    const std::size_t size = m.size();
    const Scalar number_density = molar_density * (AVOGADRO / 1.0e30);
    std::vector<double> diameter(size);
    double mean_segment_count = 0.0;
    for (std::size_t i = 0; i < size; ++i) {
        diameter[i] = sigma[i] *
            (1.0 - 0.12 * std::exp(-3.0 * epsilon_over_k[i] / temperature));
        mean_segment_count += mole_fractions[i] * m[i];
    }

    std::array<Scalar, 4> zeta = {
        Scalar(0.0), Scalar(0.0), Scalar(0.0), Scalar(0.0),
    };
    for (std::size_t order = 0; order < zeta.size(); ++order) {
        double moment = 0.0;
        for (std::size_t i = 0; i < size; ++i) {
            moment += mole_fractions[i] * m[i] *
                std::pow(diameter[i], static_cast<int>(order));
        }
        zeta[order] = (PI / 6.0) * number_density * moment;
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
        const double contact_diameter =
            diameter[i] * diameter[i] / (diameter[i] + diameter[i]);
        const Scalar contact =
            Scalar(1.0) / one_minus_eta +
            contact_diameter * 3.0 * zeta[2] /
                (one_minus_eta * one_minus_eta) +
            contact_diameter * contact_diameter * 2.0 * zeta[2] * zeta[2] /
                (one_minus_eta * one_minus_eta * one_minus_eta);
        hard_chain -= mole_fractions[i] * (m[i] - 1.0) * scalar_log(contact);
    }

    const double chain_factor_1 =
        (mean_segment_count - 1.0) / mean_segment_count;
    const double chain_factor_2 = chain_factor_1 *
        (mean_segment_count - 2.0) / mean_segment_count;
    std::array<double, 7> a{};
    std::array<double, 7> b{};
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

    double first_moment = 0.0;
    double second_moment = 0.0;
    for (std::size_t i = 0; i < size; ++i) {
        for (std::size_t j = 0; j < size; ++j) {
            const double sigma_ij = 0.5 * (sigma[i] + sigma[j]);
            const double epsilon_ij =
                std::sqrt(epsilon_over_k[i] * epsilon_over_k[j]) *
                (1.0 - kij[i * size + j]);
            const double weight = mole_fractions[i] * mole_fractions[j] *
                m[i] * m[j] * sigma_ij * sigma_ij * sigma_ij;
            first_moment += weight * epsilon_ij / temperature;
            second_moment += weight * epsilon_ij * epsilon_ij /
                (temperature * temperature);
        }
    }
    const Scalar dispersion =
        -2.0 * PI * number_density * integral_1 * first_moment -
        PI * number_density * mean_segment_count * c1 * integral_2 * second_moment;
    return {hard_chain, dispersion};
}

}  // namespace detail

}  // namespace epcsaft
