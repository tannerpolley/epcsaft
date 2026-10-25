#pragma once

#include <cppad/cppad.hpp>

#include <array>
#include <cmath>
#include <cstddef>
#include <type_traits>

// One dispatch for the transcendentals, so a kernel reads the same whether it is
// evaluating a double, a Taylor series or recording a tape.

namespace epcsaft {

// Truncated Taylor series in one variable, c[k] = f^(k)(t0)/k!: forward-mode derivatives of the
// shared kernels along one coordinate without a tape (density closure). Value arithmetic is the
// double arithmetic of the same expression; the operators are forced inline for the closure's inner loop.
template <std::size_t N>
struct Taylor final {
    std::array<double, N + 1> c{};
    Taylor() = default;
    Taylor(double value) : c{value} {}  // NOLINT(google-explicit-constructor): constants enter the kernels as values
    [[gnu::always_inline]] Taylor& operator+=(const Taylor& o) { _Pragma("GCC unroll 8") for (std::size_t k = 0; k <= N; ++k) c[k] += o.c[k]; return *this; }
    [[gnu::always_inline]] Taylor& operator-=(const Taylor& o) { _Pragma("GCC unroll 8") for (std::size_t k = 0; k <= N; ++k) c[k] -= o.c[k]; return *this; }
    [[gnu::always_inline]] Taylor& operator*=(const Taylor& o) { return *this = *this * o; }
    [[gnu::always_inline]] Taylor& operator/=(const Taylor& o) { return *this = *this / o; }
    [[gnu::always_inline]] friend Taylor operator-(Taylor a) { _Pragma("GCC unroll 8") for (double& v : a.c) v = -v; return a; }
    [[gnu::always_inline]] friend Taylor operator+(Taylor a, const Taylor& b) { return a += b; }
    [[gnu::always_inline]] friend Taylor operator+(Taylor a, double b) { a.c[0] += b; return a; }
    [[gnu::always_inline]] friend Taylor operator+(double a, Taylor b) { b.c[0] = a + b.c[0]; return b; }
    [[gnu::always_inline]] friend Taylor operator-(Taylor a, const Taylor& b) { return a -= b; }
    [[gnu::always_inline]] friend Taylor operator-(Taylor a, double b) { a.c[0] -= b; return a; }
    [[gnu::always_inline]] friend Taylor operator-(double a, Taylor b) { b = -b; b.c[0] = a + b.c[0]; return b; }
    [[gnu::always_inline]] friend Taylor operator*(const Taylor& a, const Taylor& b) {
        Taylor r;
        _Pragma("GCC unroll 8")
        for (std::size_t k = 0; k <= N; ++k) {
            r.c[k] = a.c[0] * b.c[k];
            _Pragma("GCC unroll 8")
            for (std::size_t j = 1; j <= k; ++j) r.c[k] += a.c[j] * b.c[k - j];
        }
        return r;
    }
    [[gnu::always_inline]] friend Taylor operator*(Taylor a, double b) { _Pragma("GCC unroll 8") for (double& v : a.c) v *= b; return a; }
    [[gnu::always_inline]] friend Taylor operator*(double a, Taylor b) { _Pragma("GCC unroll 8") for (double& v : b.c) v = a * v; return b; }
    // Quotients keep the value's double division; the derivative coefficients multiply by one reciprocal.
    [[gnu::always_inline]] friend Taylor operator/(const Taylor& a, const Taylor& b) {
        Taylor q;
        q.c[0] = a.c[0] / b.c[0];
        const double inverse = 1.0 / b.c[0];
        _Pragma("GCC unroll 8")
        for (std::size_t k = 1; k <= N; ++k) {
            double s = a.c[k];
            _Pragma("GCC unroll 8")
            for (std::size_t j = 1; j <= k; ++j) s -= b.c[j] * q.c[k - j];
            q.c[k] = s * inverse;
        }
        return q;
    }
    [[gnu::always_inline]] friend Taylor operator/(Taylor a, double b) {
        const double inverse = 1.0 / b;
        a.c[0] /= b;
        _Pragma("GCC unroll 8")
        for (std::size_t k = 1; k <= N; ++k) a.c[k] *= inverse;
        return a;
    }
    [[gnu::always_inline]] friend Taylor operator/(double a, const Taylor& b) { return Taylor(a) / b; }
    // exp: k e_k = sum_j j a_j e_(k-j); log of (shift + a): (shift + a_0) l_k = a_k - sum_(j<k) j l_j a_(k-j) / k.
    [[gnu::always_inline]] friend Taylor exp_series(const Taylor& a, double e0, double shift) {
        Taylor e;
        e.c[0] = e0;
        _Pragma("GCC unroll 8")
        for (std::size_t k = 1; k <= N; ++k) {
            double s = 0.0;
            _Pragma("GCC unroll 8")
            for (std::size_t j = 1; j <= k; ++j)
                s += (static_cast<double>(j) / static_cast<double>(k)) * a.c[j] * (j == k ? e0 + shift : e.c[k - j]);
            e.c[k] = s;
        }
        return e;
    }
    [[gnu::always_inline]] friend Taylor log_series(const Taylor& a, double l0, double shift) {
        Taylor l;
        l.c[0] = l0;
        const double inverse = 1.0 / (shift + a.c[0]);
        _Pragma("GCC unroll 8")
        for (std::size_t k = 1; k <= N; ++k) {
            double s = a.c[k];
            _Pragma("GCC unroll 8")
            for (std::size_t j = 1; j < k; ++j) s -= (static_cast<double>(j) / static_cast<double>(k)) * l.c[j] * a.c[k - j];
            l.c[k] = s * inverse;
        }
        return l;
    }
    [[gnu::always_inline]] friend Taylor sqrt(const Taylor& a) {
        Taylor s;
        s.c[0] = std::sqrt(a.c[0]);
        const double inverse = 0.5 / s.c[0];
        _Pragma("GCC unroll 8")
        for (std::size_t k = 1; k <= N; ++k) {
            double v = a.c[k];
            _Pragma("GCC unroll 8")
            for (std::size_t j = 1; j < k; ++j) v -= s.c[j] * s.c[k - j];
            s.c[k] = v * inverse;
        }
        return s;
    }
};

template <typename Scalar>
inline constexpr bool is_taylor = false;
template <std::size_t N>
inline constexpr bool is_taylor<Taylor<N>> = true;

template <typename Scalar>
// EqExclude: scalar utility dispatch
Scalar scalar_log1p(const Scalar& value) {
    if constexpr (std::is_same_v<Scalar, double>) return std::log1p(value);
    else if constexpr (is_taylor<Scalar>) return log_series(value, std::log1p(value.c[0]), 1.0);
    else return CppAD::log1p(value);
}

template <typename Scalar>
// EqExclude: scalar utility dispatch
Scalar scalar_expm1(const Scalar& value) {
    if constexpr (std::is_same_v<Scalar, double>) return std::expm1(value);
    else if constexpr (is_taylor<Scalar>) return exp_series(value, std::expm1(value.c[0]), 1.0);
    else return CppAD::expm1(value);
}

template <typename Scalar>
// EqExclude: scalar utility dispatch; condition remains live on an AD tape
Scalar scalar_select_less(const Scalar& left, const Scalar& right,
                          const Scalar& yes, const Scalar& no) {
    if constexpr (std::is_same_v<Scalar, double>) return left < right ? yes : no;
    else if constexpr (is_taylor<Scalar>) return left.c[0] < right.c[0] ? yes : no;
    else return CppAD::CondExpLt(left, right, yes, no);
}

template <typename Scalar>
// EqExclude: scalar utility dispatch
Scalar scalar_log(const Scalar& value) {
    if constexpr (std::is_same_v<Scalar, double>) return std::log(value);
    else if constexpr (is_taylor<Scalar>) return log_series(value, std::log(value.c[0]), 0.0);
    else return CppAD::log(value);
}

template <typename Scalar>
// EqExclude: scalar utility dispatch
Scalar scalar_exp(const Scalar& value) {
    if constexpr (std::is_same_v<Scalar, double>) return std::exp(value);
    else if constexpr (is_taylor<Scalar>) return exp_series(value, std::exp(value.c[0]), 0.0);
    else return CppAD::exp(value);
}

template <typename Scalar>
// EqExclude: scalar utility dispatch
Scalar scalar_sqrt(const Scalar& value) {
    if constexpr (std::is_same_v<Scalar, double>) return std::sqrt(value);
    else if constexpr (is_taylor<Scalar>) return sqrt(value);
    else return CppAD::sqrt(value);
}

}  // namespace epcsaft
