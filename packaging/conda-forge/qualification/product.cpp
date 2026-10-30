#include "core/linear_error.hpp"
#include <iomanip>
#include <iostream>
int main() {
    using namespace epcsaft;
    const WiderScalar delta = WiderScalar(std::ldexp(1.0, -40));
    const WiderScalar a = 1 + delta, b = 1 - delta;
    const WiderScalar expected = -delta * delta;
#ifdef EPCSAFT_SOFTWARE_WIDER
    const char* branch = "product_residual";
    const WiderScalar residual = core::product_residual(a, b, WiderScalar(1));
#else
    const char* branch = "fma";
    const WiderScalar residual = std::fma(a, b, WiderScalar(-1));
#endif
    std::cout << std::setprecision(25) << branch << ' ' << residual << " expected " << expected << '\n';
    return residual == expected ? 0 : 1;
}
