#pragma once

#include <vector>

namespace epcsaft {

inline constexpr int ELECTROLYTE_MODEL_NONE = 0;
inline constexpr int ELECTROLYTE_MODEL_SSM_DS = 1;
inline constexpr int ELECTROLYTE_MODEL_ORIGINAL_BORN = 2;
inline constexpr int ELECTROLYTE_MODEL_FULLY_DISSOCIATED_DEBYE_HUCKEL = 3;
inline constexpr int RELATIVE_PERMITTIVITY_FORMULATION_NONE = 0;
inline constexpr int RELATIVE_PERMITTIVITY_FORMULATION_ION_FRACTION_SUPPRESSION = 1;
inline constexpr int RELATIVE_PERMITTIVITY_FORMULATION_COMPONENT_MIXING = 2;
inline constexpr int RELATIVE_PERMITTIVITY_FORMULATION_DENSITY_VOLUME_TEMPERATURE = 3;
inline constexpr int RELATIVE_PERMITTIVITY_FORMULATION_ION_SPECIFIC_SUPPRESSION = 4;
inline constexpr int RELATIVE_PERMITTIVITY_FORMULATION_DIPOLAR_DENSITY_PERTURBATION = 5;
inline constexpr int RELATIVE_PERMITTIVITY_FORMULATION_SOLVENT_ONLY = 6;
inline constexpr int POLAR_FORMULATION_NONE = 0;
inline constexpr int POLAR_FORMULATION_GROSS_VRABEC_POINT_MULTIPOLE = 1;

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
    double ion_fraction_suppression_coefficient;
    std::vector<double> born_diameter_angstrom;
    std::vector<double> solvation_factor;
    double ionic_region_relative_permittivity;
    int electrolyte_model;
    std::vector<double> dielectric_volume_m3_per_mol;
    std::vector<double> dielectric_temperature_k;
    std::vector<double> ion_specific_suppression_coefficient;
    std::vector<double> dipole_scaling_debye2;
    std::vector<double> polarizability_scaling_angstrom3;
    std::vector<double> correlation_integral_parameter;
    int relative_permittivity_formulation;
    std::vector<double> dipole_moment_debye;
    std::vector<double> quadrupole_moment_debye_angstrom;
    int polar_formulation;
};

}  // namespace epcsaft
