#pragma once

#include <cppad/cppad.hpp>

#include <cmath>
#include <type_traits>

// One dispatch for the transcendentals, so a kernel reads the same whether it is
// evaluating a double or recording a tape. The tree this replaces defined these
// three times under three prefixes.

namespace epcsaft {

template <typename Scalar>
// EqExclude: scalar utility dispatch
Scalar scalar_log1p(const Scalar& value) {
    if constexpr (std::is_same_v<Scalar, double>) return std::log1p(value);
    return CppAD::log1p(value);
}

template <typename Scalar>
// EqExclude: scalar utility dispatch
Scalar scalar_expm1(const Scalar& value) {
    if constexpr (std::is_same_v<Scalar, double>) return std::expm1(value);
    return CppAD::expm1(value);
}

template <typename Scalar>
// EqExclude: scalar utility dispatch; condition remains live on an AD tape
Scalar scalar_select_less(const Scalar& left, const Scalar& right,
                          const Scalar& yes, const Scalar& no) {
    if constexpr (std::is_same_v<Scalar, double>) return left < right ? yes : no;
    return CppAD::CondExpLt(left, right, yes, no);
}

template <typename Scalar>
// EqExclude: scalar utility dispatch
Scalar scalar_log(const Scalar& value) {
    if constexpr (std::is_same_v<Scalar, double>) {
        return std::log(value);
    }
    return CppAD::log(value);
}

template <typename Scalar>
// EqExclude: scalar utility dispatch
Scalar scalar_exp(const Scalar& value) {
    if constexpr (std::is_same_v<Scalar, double>) {
        return std::exp(value);
    }
    return CppAD::exp(value);
}

template <typename Scalar>
// EqExclude: scalar utility dispatch
Scalar scalar_sqrt(const Scalar& value) {
    if constexpr (std::is_same_v<Scalar, double>) {
        return std::sqrt(value);
    }
    return CppAD::sqrt(value);
}

}  // namespace epcsaft
