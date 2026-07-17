#pragma once

#include <vector>

namespace epcsaft {

struct NativeParameters final {
    std::vector<double> segment_count;
    std::vector<double> sigma_constant;
    std::vector<int> sigma_term_component;
    std::vector<double> sigma_term_amplitude;
    std::vector<double> sigma_term_exponent;
    std::vector<double> dispersion_energy_over_k;
    std::vector<double> k_ij;
    std::vector<int> site_component;
    std::vector<int> site_multiplicity;
    std::vector<int> association_active;
    std::vector<double> association_energy_over_k;
    std::vector<double> association_volume;
};

}  // namespace epcsaft
