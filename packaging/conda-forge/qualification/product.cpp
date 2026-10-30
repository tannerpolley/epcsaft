#include "core/linear_error.hpp"
#include <iostream>
int main() {
    using namespace epcsaft;
    const WiderScalar delta = WiderScalar(std::ldexp(1.0, -40));
    const WiderScalar a = 1 + delta, b = 1 - delta;
    const WiderScalar expected = -delta * delta;
    const WiderScalar residual = core::product_residual(a, b, WiderScalar(1));
    std::cout << std::setprecision(25) << residual << " expected " << expected << '\n';
    return residual == expected ? 0 : 1;
}
