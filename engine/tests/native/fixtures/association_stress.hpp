#pragma once
#include "../properties_support.hpp"
#include "contributions/association.hpp"

namespace properties_support {
inline epcsaft::Parameters association_stress(std::size_t sites, double q) {
    constexpr double temperature = 300.0;
    constexpr double density = 1000.0;
    auto parameters = ethane();
    parameters.site_component.assign(sites, 0);
    parameters.site_multiplicity.assign(sites, 1);
    parameters.association_active.assign(sites * sites, 0);
    parameters.association_energy_over_k.assign(sites * sites, 0.0);
    parameters.association_volume.assign(sites * sites, 0.0);
    parameters.association_active[0] = 1;
    parameters.association_volume[0] = 0.1;
    parameters.association_energy_over_k[0] = temperature * std::log(2.0);
    const auto sigma = epcsaft::segment_sigma(temperature, parameters);
    const auto diameter = epcsaft::segment_diameter(temperature, parameters, sigma);
    const double unit_delta = epcsaft::association_delta(
        temperature, density, {1.0}, parameters, sigma, diameter)[0];
    const double number_density = density * epcsaft::NUMBER_DENSITY_PER_ANGSTROM3;
    parameters.association_energy_over_k[0] = temperature
        * std::log1p(q / (number_density * unit_delta));
    return parameters;
}

}  // namespace properties_support
