#pragma once

#include "contributions/msa.hpp"

namespace epcsaft {

class MsaScreeningAtomic final : public CppAD::atomic_four<double> {
    template <std::size_t N = 0, typename Operation>
    static bool at_order(std::size_t order, const Operation& run) {
        if (order == N) {
            run(std::integral_constant<std::size_t, N>{});
            return true;
        }
        if constexpr (N < 5) return at_order<N + 1>(order, run);
        return false;
    }

    template <std::size_t N, typename Base>
    static std::vector<Taylor<N, Base>> series(const CppAD::vector<Base>& x) {
        constexpr std::size_t width = N + 1;
        std::vector<Taylor<N, Base>> result(x.size() / width);
        for (std::size_t i = 0; i < result.size(); ++i)
            for (std::size_t k = 0; k < width; ++k) result[i].c[k] = x[i * width + k];
        return result;
    }

    template <typename Base>
    bool forward_coefficients(std::size_t low, std::size_t up, const CppAD::vector<Base>& x, CppAD::vector<Base>& y) {
        if (low == 0) {
            std::vector<Base> central;
            for (std::size_t i = 0; i < x.size(); i += up + 1) central.push_back(x[i]);
            if constexpr (std::is_same_v<Base, double>) y[0] = solve_msa_gamma(central);
            else {
                std::vector<Base> value(1);
                (*this)(central, value);
                y[0] = value[0];
            }
        }
        if (up == 0) return true;
        const std::size_t first = std::max(std::size_t{1}, low);
        return at_order(up, [&](auto order) {
            Taylor<order, Base> gamma;
            for (std::size_t k = 0; k < first; ++k) gamma.c[k] = y[k];
            msa_gamma_coefficients(series<order>(x), gamma, first, up);
            for (std::size_t k = first; k <= up; ++k) y[k] = gamma.c[k];
        });
    }

public:
    MsaScreeningAtomic() : atomic_four("msa_screening") {}

    bool for_type(std::size_t, const CppAD::vector<CppAD::ad_type_enum>& x,
                  CppAD::vector<CppAD::ad_type_enum>& y) override {
        y[0] = *std::max_element(x.begin(), x.end());
        return true;
    }

    bool rev_depend(std::size_t, const CppAD::vector<bool>&, CppAD::vector<bool>& x,
                    const CppAD::vector<bool>& y) override {
        for (std::size_t i = 0; i < x.size(); ++i) x[i] = y[0];
        return true;
    }

    bool forward(std::size_t, const CppAD::vector<bool>&, std::size_t low, std::size_t up,
                 const CppAD::vector<double>& x, CppAD::vector<double>& y) override {
        return forward_coefficients(low, up, x, y);
    }

    bool forward(std::size_t, const CppAD::vector<bool>&, std::size_t low, std::size_t up,
                 const CppAD::vector<CppAD::AD<double>>& x, CppAD::vector<CppAD::AD<double>>& y) override {
        return forward_coefficients(low, up, x, y);
    }

    bool reverse(std::size_t, const CppAD::vector<bool>& selected, std::size_t up, const CppAD::vector<double>& x,
                 const CppAD::vector<double>& y, CppAD::vector<double>& dx, const CppAD::vector<double>& dy) override {
        return at_order(up, [&](auto order) {
            using Series = Taylor<order>;
            const auto inputs = series<order>(x);
            const auto gamma = series<order>(y).front();
            const auto m = msa_moments<MsaEvaluation::Slope>(gamma, inputs);
            const Series rhs = scalar_sqrt(PI * inputs[0]) * scalar_sqrt(m.q);
            // Transpose the Taylor convolution of -R_q/R_gamma with the input variations.
            const Series factor = PI / (2.0 * rhs - PI * inputs[0] * m.q_gamma);
            const auto pullback = [&](std::size_t i, const Series& tangent) {
                for (std::size_t k = 0; k <= up; ++k) {
                    double coefficient = 0.0;
                    if (selected[i])
                        for (std::size_t j = k; j <= up; ++j) coefficient += dy[j] * tangent.c[j - k];
                    dx[i * (up + 1) + k] = coefficient;
                }
            };
            pullback(0, factor * m.q);
            const Series dielectric = factor * inputs[0];
            const Series coupling = 2.0 * m.h_moment / m.denominator;
            for (std::size_t i = 1; i < inputs.size(); i += 3) {
                const Series& n = inputs[i];
                const Series& d = inputs[i + 1];
                const Series& z = inputs[i + 2];
                const Series h = 1.0 + gamma * d;
                const Series e = (z - m.coupling * d * d) / h;
                const Series qn = e * e - coupling * d * (z / h - m.coupling * d * d * (1.0 / h - 1.0 / 3.0));
                const Series qd = n * (-2.0 * e * (gamma * e + 2.0 * m.coupling * d) / h
                    - coupling * (z / (h * h) - m.coupling * d * d * (3.0 / h - gamma * d / (h * h) - 1.0)));
                const Series qz = n * (2.0 * e - coupling * d) / h;
                pullback(i, dielectric * qn);
                pullback(i + 1, dielectric * qd);
                pullback(i + 2, dielectric * qz);
            }
        });
    }
};

inline MsaScreeningAtomic& msa_screening_atomic() {
    // CppAD keeps the registration for the lifetime of the parallel runtime.
    static auto* const atomic = new MsaScreeningAtomic;
    return *atomic;
}

inline CppAD::AD<double> msa_gamma(const std::vector<CppAD::AD<double>>& inputs) {
    std::vector<CppAD::AD<double>> result(1);
    msa_screening_atomic()(inputs, result);
    return result.front();
}

}  // namespace epcsaft
