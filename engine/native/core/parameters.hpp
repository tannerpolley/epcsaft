#pragma once

#include "core/declared.hpp"

#include <array>
#include <string>
#include <vector>

// The parameter arrays the equations read. Plain data, no validation: the model
// layer builds and checks one of these once, and the physics layer trusts it.
// Shared parameter data lives in Core because physics may not include the model layer.

namespace epcsaft {

// Published forms. Ion suppression eps_sf / (1 + sum_i alpha_i x_i) covers
// ion-specific suppression directly, ion-fraction suppression with alpha_i equal
// for every ion, and solvent-only with alpha_i = 0; the model layer fills
// ion_suppression_coefficient accordingly.
enum class PermittivityForm {
    SolventMoleWeighting = 5, // declared solvent mole fractions
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
    ExcludeSameSign = 2,  // pairs with z_i z_j > 0, the diagonal included, take 0 unless both are charged molecules
};

// Quantities a temperature law may own (equation pair_parameter_temperature_law).
enum class TemperatureLaw {
    Kij = 1,                    // component pair
    Lij = 2,                    // component pair
    AssociationEnergyOverK = 3, // explicit site-pair edge, K
    AssociationVolume = 4,      // explicit site-pair edge
};

struct FixedParameters {
    bool ideal_gas = false;
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

    // Temperature laws q(T) = sum_p a_p (T - T_r)^p + s (1/T - 1/T_r): a reference per law; per coefficient
    // its law, power p (-1 for s) and source identity; per owned slot (target, row-major slot with the
    // smaller index first, law). Empty unless declared.
    std::vector<double> temperature_law_reference_k;
    std::vector<int> temperature_law_term_law;
    std::vector<int> temperature_law_term_power;
    std::vector<std::string> temperature_law_coefficient_identity;
    std::vector<std::array<int, 3>> temperature_law_slot;
    // Combining-rule edges evaluated at each temperature: edge and mirror, its two explicit self edges (site-pair
    // slots), and the law of its volume correction c, kappa_ij = sqrt(kappa_ii kappa_jj) (1 - c), or -1.
    std::vector<std::array<int, 5>> association_thermal_rule;

    // Electrolyte
    std::vector<int> charge_number;
    std::vector<int> molecular_ion;  // 1: an ion keeping the neutral short-range terms (#641); empty means none
    std::vector<double> solvation_factor;
    std::vector<std::array<int, 4>> cation_mixing;  // solvent, unordered cations, common anion
    double born_c_shell = 0.0;
    core::Declared<ElectrolyteForm> electrolyte_form;
    core::Declared<IonicDispersion> ionic_dispersion;

    // Permittivity
    core::Declared<PermittivityForm> permittivity_form;
    std::vector<int> permittivity_solvent;
    bool permittivity_solvent_normalized = false;
    std::vector<double> relative_permittivity;
    std::vector<int> permittivity_term_component;
    std::vector<double> permittivity_term_coefficient;
    std::vector<int> permittivity_term_power;
    std::vector<double> permittivity_reference_temperature_k;
    std::vector<double> dipole_scaling_debye2;
    std::vector<double> polarizability_scaling_angstrom3;
    std::vector<double> correlation_integral_parameter;

    // Optional residual dipole-dipole (ADR 0011) or quadrupole-quadrupole (ADR 0014) term; empty unless declared
    std::vector<double> dipole_moment_debye;
    std::vector<double> dipole_count;
    std::vector<double> quadrupole_moment_debye_angstrom;
    std::vector<double> quadrupole_count;

    std::vector<int> form_identity() const {
        std::vector<int> result = {static_cast<int>(ideal_gas), electrolyte_form.declared() ? static_cast<int>(electrolyte_form.require("electrolyte")) : 0,
                ionic_dispersion.declared() ? static_cast<int>(ionic_dispersion.require("dispersion")) : 0,
                permittivity_form.declared() ? static_cast<int>(permittivity_form.require("permittivity")) : 0};
        for (const auto& record : cation_mixing) result.insert(result.end(), record.begin(), record.end());
        if (permittivity_form.is(PermittivityForm::SolventMoleWeighting)) {
            result.push_back(permittivity_solvent_normalized);
            result.insert(result.end(), permittivity_solvent.begin(), permittivity_solvent.end());
        }
        return result;
    }

};

// Active parameter families; topology and source conventions remain fixed.
// Scalar is supplied by the consuming layer, so Core has no AD dependency.
template <class Scalar>
struct ParameterRecord final : FixedParameters {
    Scalar ionic_region_relative_permittivity = 0.0;
    Scalar born_c_dielectric = 0.0;
    std::vector<Scalar> packing_diameter_angstrom;
    std::vector<Scalar> debye_huckel_diameter_angstrom;
    std::vector<Scalar> born_diameter_angstrom;
    std::vector<Scalar> ion_suppression_coefficient;
    std::vector<Scalar> segment_count;
    std::vector<Scalar> sigma_constant;
    std::vector<Scalar> dispersion_energy_over_k;
    std::vector<Scalar> k_ij;
    std::vector<Scalar> k_ij_reciprocal_temperature_slope;
    std::vector<Scalar> l_ij;
    std::vector<Scalar> association_energy_over_k;
    std::vector<Scalar> association_volume;
    std::vector<Scalar> dielectric_volume_m3_per_mol;
    std::vector<Scalar> dielectric_temperature_k;
    std::vector<Scalar> temperature_law_coefficient;
    std::vector<Scalar> cation_pair_mixing;
    std::vector<Scalar> cation_triple_mixing;

    std::size_t component_count() const { return ideal_gas ? molar_mass_kg_per_mol.size() : segment_count.size(); }
    std::size_t site_count() const { return site_component.size(); }
};

using Parameters = ParameterRecord<double>;

}  // namespace epcsaft
