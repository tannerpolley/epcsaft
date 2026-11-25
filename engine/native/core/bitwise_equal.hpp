#pragma once

#include <cstddef>
#include <cstring>

namespace epcsaft::core {

inline bool bitwise_equal(const double* a, std::size_t a_size, const double* b, std::size_t b_size) {
    return a_size == b_size && (a_size == 0 || std::memcmp(a, b, a_size * sizeof(double)) == 0);
}

inline bool bitwise_equal(double a, double b) { return bitwise_equal(&a, 1, &b, 1); }

template <class Doubles>
bool bitwise_equal(const Doubles& a, const Doubles& b) {
    return bitwise_equal(a.data(), a.size(), b.data(), b.size());
}

}  // namespace epcsaft::core
