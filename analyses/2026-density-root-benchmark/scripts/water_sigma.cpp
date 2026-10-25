// Materialize the Held packet's diameter with the Engine's EqID segment_sigma_correlation.
#include "contributions/helmholtz.hpp"
#include <iomanip>
#include <iostream>

struct SigmaCoefficients {
    std::vector<double> sigma_constant, sigma_term_amplitude, sigma_term_exponent;
    std::vector<int> sigma_term_component;
};

int main() {
    double temperature, constant;
    std::size_t terms;
    while (std::cin >> temperature >> constant >> terms) {
        SigmaCoefficients coefficients;
        coefficients.sigma_constant = {constant};
        for (std::size_t i = 0; i < terms; ++i) {
            double amplitude, exponent;
            std::cin >> amplitude >> exponent;
            coefficients.sigma_term_component.push_back(0);
            coefficients.sigma_term_amplitude.push_back(amplitude);
            coefficients.sigma_term_exponent.push_back(exponent);
        }
        if (!std::cin) return 1;
        std::cout << std::setprecision(17) << epcsaft::segment_sigma(temperature, coefficients)[0] << '\n';
    }
}
