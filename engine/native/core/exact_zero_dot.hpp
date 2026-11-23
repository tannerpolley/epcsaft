#pragma once
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace epcsaft::core {

inline bool exact_zero_dot(const std::vector<double>& left, const std::vector<double>& right) {
    if (left.size() != right.size()) throw std::invalid_argument("exact dot width");
    std::vector<double> expansion;
    const auto add = [&](double value) {
        std::vector<double> next;
        for (double term : expansion) {
            const double sum = value + term, recovered = sum - value;
            const double error = (value - (sum - recovered)) + (term - recovered);
            if (error != 0.0) next.push_back(error);
            value = sum;
        }
        if (value != 0.0) next.push_back(value);
        expansion = std::move(next);
    };
    for (std::size_t i = 0; i < left.size(); ++i) {
        if (left[i] == 0.0 || right[i] == 0.0) continue;
        const double product = left[i] * right[i];
        // Decline certification if an error-free product residual could underflow.
        if (!std::isfinite(product) || std::abs(product) < std::ldexp(std::numeric_limits<double>::min(), 53)) return false;
        add(std::fma(left[i], right[i], -product));
        add(product);
    }
    return expansion.empty();
}

}  // namespace epcsaft::core
