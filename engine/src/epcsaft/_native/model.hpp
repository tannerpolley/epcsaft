#pragma once

#include <vector>

namespace epcsaft {

inline constexpr int ELECTROLYTE_MODEL_NONE = 0;
inline constexpr int ELECTROLYTE_MODEL_FIGIEL_SSM_DS = 1;
inline constexpr int ELECTROLYTE_MODEL_ASCANI_ORIGINAL_BORN = 2;
inline constexpr int ELECTROLYTE_MODEL_HELD_2008_DEBYE_HUCKEL = 3;

struct NativeParameters final {
    std::vector<double> segment_count;
    std::vector<double> sigma_constant;
    std::vector<int> sigma_term_component;
    std::vector<double> sigma_term_amplitude;
    std::vector<double> sigma_term_exponent;
    std::vector<double> dispersion_energy_over_k;
    std::vector<double> k_ij;
    std::vector<double> l_ij;
    std::vector<int> site_component;
    std::vector<int> site_multiplicity;
    std::vector<int> association_active;
    std::vector<double> association_energy_over_k;
    std::vector<double> association_volume;
    std::vector<double> molar_mass_kg_per_mol;
    std::vector<int> charge_number;
    std::vector<double> relative_permittivity;
    double dielectric_ion_suppression_coefficient;
    std::vector<double> born_diameter_angstrom;
    std::vector<double> solvation_factor;
    double ionic_region_relative_permittivity;
    int electrolyte_model;
};

}  // namespace epcsaft
