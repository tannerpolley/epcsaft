#pragma once

#include <cstring>

namespace epcsaft::core {

template <class Doubles>
bool bitwise_equal(const Doubles& a, const Doubles& b) {
    return a.size() == b.size()
        && (a.size() == 0 || std::memcmp(a.data(), b.data(), a.size() * sizeof(double)) == 0);
}

}  // namespace epcsaft::core
