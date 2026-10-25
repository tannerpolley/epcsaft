#include "properties/thermochemistry.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>

int main() {
    std::string reaction;
    double temperature, a1, a2, a3, a4;
    std::cout << "reaction_id,T_K,engine_ln_K,engine_log10_K\n" << std::setprecision(17);
    while (std::cin >> reaction >> temperature >> a1 >> a2 >> a3 >> a4) {
        epcsaft::ReactionLogPolynomial law;
        law.reference_temperature = 298.15;
        law.temperature_min = 298.15;
        law.temperature_max = 423.15;
        law.standard_state_id = "PHREEQC-molality-activity";
        // PHREEQC manual SOLUTION_SPECIES: exact logarithm-base and Tref mapping.
        law.a = std::log(10.0) * a1 + a4 * std::log(law.reference_temperature);
        law.b = std::log(10.0) * a3;
        law.c = a4;
        law.d = std::log(10.0) * a2;
        const auto value = epcsaft::evaluate_reaction_properties(law, temperature);
        if (!std::isfinite(value.ln_k)) throw std::runtime_error("nonfinite calcium law");
        std::cout << reaction << ',' << temperature << ',' << value.ln_k << ','
                  << value.ln_k / std::log(10.0) << '\n';
    }
    if (!std::cin.eof()) throw std::runtime_error("invalid source-law input");
}
