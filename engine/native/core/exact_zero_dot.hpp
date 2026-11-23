#pragma once
#include "core/wider_scalar.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace epcsaft::core {

inline bool exact_zero_dot(const std::vector<double>& left, const std::vector<double>& right) {
    if (left.size() != right.size()) throw std::invalid_argument("exact dot width");
    using Limits = std::numeric_limits<double>;
    // At this product exponent an exact FMA residual fits in binary64, including subnormals.
    constexpr int minimum_product_exponent = Limits::min_exponent + Limits::digits - 1;
    int shift = std::numeric_limits<int>::min();
    for (std::size_t i = 0; i < left.size(); ++i) {
        if (!std::isfinite(left[i]) || !std::isfinite(right[i])) return false;
        if (left[i] == 0.0 || right[i] == 0.0) continue;
        int a_exponent, b_exponent;
        const double a = std::frexp(std::abs(left[i]), &a_exponent);
        const double b = std::frexp(std::abs(right[i]), &b_exponent);
        const int exponent = a_exponent + b_exponent - (std::fma(a, b, -0.5) < 0 ? 2 : 1);
        shift = std::max(shift, minimum_product_exponent - exponent);
    }
    // Wider exponent range also admits partial sums of finite products before cancellation.
    std::vector<WiderScalar> expansion;
    const auto add = [&](WiderScalar value) {
        std::vector<WiderScalar> next;
        for (const auto& term : expansion) {
            const WiderScalar sum = value + term;
            const WiderScalar recovered = sum - value;
            const WiderScalar error = (value - (sum - recovered)) + (term - recovered);
            if (error != 0.0) next.push_back(error);
            value = sum;
        }
        if (value != 0.0) next.push_back(value);
        expansion = std::move(next);
    };
    for (std::size_t i = 0; i < left.size(); ++i) {
        if (left[i] == 0.0 || right[i] == 0.0) continue;
        int a_exponent, b_exponent;
        const double a = std::frexp(left[i], &a_exponent);
        const double b = std::frexp(right[i], &b_exponent);
        const int exponent = a_exponent + b_exponent + shift;
        if (exponent > Limits::max_exponent + 1) return false;
        // Apply the common shift to the left exponent, then redistribute powers of two
        // between the operands to keep each scaling exact and each operand normal.
        const double scaled_left = std::ldexp(a, exponent / 2);
        const double scaled_right = std::ldexp(b, exponent - exponent / 2);
        const double product = scaled_left * scaled_right;
        if (!std::isfinite(product)) return false;
        add(std::fma(scaled_left, scaled_right, -product));
        add(product);
    }
    return expansion.empty();
}

}  // namespace epcsaft::core
