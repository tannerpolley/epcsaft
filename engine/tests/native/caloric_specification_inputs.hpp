#pragma once
// Manufactured DIPPR107 records adopted by issue 224; Gross (2001) residual parameters.
namespace caloric_specification_checks {
using namespace epcsaft;
inline ThermochemistryRecord records(bool propane = false) {
    ThermochemistryRecord r; r.standard_pressure = 100000.0;
    const std::vector<std::array<double, 5>> coefficients = propane
        ? std::vector<std::array<double, 5>>{{45,30,500,15,1000}}
        : std::vector<std::array<double, 5>>{{25,10,500,5,1000},{35,20,500,10,1000}};
    for (const auto& c : coefficients)
        r.component_intervals.push_back({{100,550,true,true,{IdealDIPPR107{c,298.15,0,0},100000}}});
    return r;
}
inline Parameters propane_parameters() {
    auto p = properties_support::parameters(1);
    p.segment_count = {2.0020}; p.sigma_constant = {3.6184}; p.dispersion_energy_over_k = {208.11};
    p.molar_mass_kg_per_mol = {0.044096}; return p;
}
}
