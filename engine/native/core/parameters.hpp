#pragma once

#include "core/declared.hpp"

#include <vector>

// The parameter arrays the equations read. Plain data, no validation: the model
// layer builds and checks one of these once, and the physics layer trusts it.
//
// Not in the target layout's six-file list. It is here because every physics
// header needs the arrays and physics may not include the model layer, so the
// data has to live at the bottom.

namespace epcsaft {

// Four published forms survive. The five the old tree carried collapsed to
// these because three of them were one equation: ion suppression
// eps_sf / (1 + sum_i alpha_i x_i) covers ion-specific suppression directly,
// ion-fraction suppression with alpha_i equal for every ion, and solvent-only
// with alpha_i = 0. The model layer fills ion_suppression_coefficient
// accordingly, so the equation here has one branch instead of three.
enum class PermittivityForm {
    IonSuppression = 1,      // eps_sf / (1 + sum alpha_i x_i)
    ComponentMixing = 2,     // linear mixing over components
    DensityTemperature = 3,  // density and temperature correlation
    DipolarPerturbation = 4, // dipolar density perturbation
};

enum class ElectrolyteForm {
    DebyeHuckel = 1,
    DebyeHuckelBorn = 2,
};

// The registered case split same_sign_ionic_dispersion. Declared by the parameter
// source, never defaulted: a set with a charged component and no declaration is
// rejected by the model constructor, and a neutral set may not declare one.
enum class IonicDispersion {
    AllIonPairs = 1,      // every pair takes sqrt(eps_i eps_j)(1 - k_ij)
    ExcludeSameSign = 2,  // pairs with z_i z_j > 0, the diagonal included, take 0
};

struct FixedParameters {
    // Chain and segment geometry
    std::vector<int> sigma_term_component;
    std::vector<double> sigma_term_amplitude;
    std::vector<double> sigma_term_exponent;
    std::vector<double> molar_mass_kg_per_mol;

    // Pair corrections, row-major over component pairs
    std::vector<double> k_ij_reference_temperature_k;

    // Association sites; pair arrays are row-major over sites
    std::vector<int> site_component;
    std::vector<int> site_multiplicity;
    std::vector<int> association_active;  // 0 inactive, 1 active, 2 Wolbach--Sandler rule edge

    // Electrolyte
    std::vector<int> charge_number;
    std::vector<double> packing_diameter_angstrom;
    std::vector<double> debye_huckel_diameter_angstrom;
    std::vector<double> born_diameter_angstrom;
    std::vector<double> solvation_factor;
    double ionic_region_relative_permittivity = 0.0;
    double born_c_shell = 0.0;
    double born_c_dielectric = 0.0;
    core::Declared<ElectrolyteForm> electrolyte_form;
    core::Declared<IonicDispersion> ionic_dispersion;

    // Permittivity
    core::Declared<PermittivityForm> permittivity_form;
    std::vector<double> relative_permittivity;
    std::vector<int> permittivity_term_component;
    std::vector<double> permittivity_term_coefficient;
    std::vector<int> permittivity_term_power;
    std::vector<double> permittivity_reference_temperature_k;
    std::vector<double> ion_suppression_coefficient;  // per component, alpha_i
    std::vector<double> dielectric_volume_m3_per_mol;
    std::vector<double> dielectric_temperature_k;
    std::vector<double> dipole_scaling_debye2;
    std::vector<double> polarizability_scaling_angstrom3;
    std::vector<double> correlation_integral_parameter;

    // Optional residual dipole-dipole term (ADR 0011); empty unless declared
    std::vector<double> dipole_moment_debye;
    std::vector<double> dipole_count;

};

// Precisely the eight active families; all electrolyte and topology inputs stay fixed.
// Scalar is supplied by the consuming layer, so Core has no AD dependency.
template <class Scalar>
struct ParameterRecord final : FixedParameters {
    std::vector<Scalar> segment_count;
    std::vector<Scalar> sigma_constant;
    std::vector<Scalar> dispersion_energy_over_k;
    std::vector<Scalar> k_ij;
    std::vector<Scalar> k_ij_reciprocal_temperature_slope;
    std::vector<Scalar> l_ij;
    std::vector<Scalar> association_energy_over_k;
    std::vector<Scalar> association_volume;

    std::size_t component_count() const { return segment_count.size(); }
    std::size_t site_count() const { return site_component.size(); }
};

using Parameters = ParameterRecord<double>;

}  // namespace epcsaft
