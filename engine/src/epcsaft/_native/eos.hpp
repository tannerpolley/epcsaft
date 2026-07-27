#pragma once

#include "model.hpp"

#include <array>
#include <cstddef>
#include <optional>
#include <string>
#include <vector>

namespace epcsaft {

struct NativeDensityDiagnostics final {
    double bracket_lower;
    double bracket_upper;
    double pressure_residual;
    double pressure_density_derivative;
    std::size_t iterations;
    bool stable;
    std::string branch;
};

struct NativeResult final {
    double hard_chain;
    double dispersion;
    double association;
    double debye_huckel;
    double born_ssm_ds;
    double bulk_relative_permittivity;
    double residual_helmholtz;
    double compressibility_factor;
    double pressure;
    std::vector<double> residual_chemical_potential_over_rt;
    std::optional<std::vector<double>> log_fugacity_coefficient;
    double residual_entropy;
    double residual_internal_energy;
    double residual_enthalpy;
    double residual_isochoric_heat_capacity;
    std::optional<double> residual_isobaric_heat_capacity;
    double molar_density;
    double pressure_density_derivative;
    std::optional<NativeDensityDiagnostics> density_diagnostics;
};

struct NativePurePhaseBlock;
struct NativeParameterizedPurePhaseBlock;
struct NativeMixturePhaseBlock;
struct NativePackingFractionBlock;

struct NativeReferenceSequence final {
    std::array<double, 2> log_fugacity;
    std::optional<double> pure_solvent_log_fugacity_coefficient;
    double reference_molality;
    double convergence_error;
};

struct NativeIonSolvationBorn final {
    double solvation_gibbs_j_per_mol;
    double derivative_j_per_mol_per_angstrom;
    double reference_molality_mol_per_kg;
    double reference_convergence_error;
};

struct NativeAqueousMiacKij final {
    double log_mean_ionic_activity_coefficient_molality;
    std::array<double, 3> derivative;
    double reference_molality_mol_per_kg;
    double reference_convergence_error;
    double reference_derivative_convergence_error;
};

struct NativeAqueousMiacSolvationFactor final {
    double log_mean_ionic_activity_coefficient_molality;
    double derivative;
    double reference_molality_mol_per_kg;
    double reference_convergence_error;
    double reference_derivative_convergence_error;
};

class NativeModel final {
public:
    explicit NativeModel(NativeParameters parameters);
    NativeResult evaluate(
        double temperature,
        double molar_density,
        const std::vector<double>& mole_fractions
    ) const;
    NativeResult evaluate_pressure(
        double temperature,
        double pressure,
        const std::vector<double>& mole_fractions,
        const std::string& phase
    ) const;
    bool supports_native_sdk_v1() const;
    bool supports_mixture_native_sdk_v1() const;
    bool supports_electrolyte_native_sdk_v1() const;
    std::size_t component_count() const;
    const std::vector<int>& charge_numbers() const;
    const std::vector<double>& molar_masses_kg_per_mol() const;
    NativePurePhaseBlock evaluate_pure_phase(
        double temperature,
        double amount,
        double volume
    ) const;
    NativeParameterizedPurePhaseBlock evaluate_pure_phase_parameters(
        double temperature,
        double amount,
        double volume,
        double segment_count,
        double segment_diameter_angstrom,
        double dispersion_energy_over_k
    ) const;
    NativeMixturePhaseBlock evaluate_mixture_phase(
        double temperature,
        const std::vector<double>& amounts,
        double volume
    ) const;
    NativeMixturePhaseBlock evaluate_mixture_phase_kij(
        double temperature,
        const std::vector<double>& amounts,
        double volume,
        double k_ij
    ) const;
    NativeMixturePhaseBlock evaluate_electrolyte_phase(
        double temperature,
        const std::vector<double>& amounts,
        double volume,
        double ion_mole_fraction_max
    ) const;
    std::array<double, 2> molar_volume_bounds(
        double temperature,
        const std::vector<double>& mole_fractions,
        double packing_fraction_min,
        double packing_fraction_max,
        double ion_mole_fraction_max
    ) const;
    NativePackingFractionBlock evaluate_packing_fraction(
        double temperature,
        const std::vector<double>& amounts,
        double volume,
        double ion_mole_fraction_max
    ) const;
    std::vector<double> reference_composition(
        const std::vector<double>& solvent_fractions,
        const std::array<std::size_t, 2>& ion_indices,
        const std::array<int, 2>& stoichiometry,
        double molality
    ) const;
    NativeReferenceSequence electrolyte_reference_sequence(
        double temperature,
        double pressure,
        const std::vector<double>& solvent_fractions,
        const std::array<std::size_t, 2>& ion_indices,
        const std::array<int, 2>& stoichiometry,
        const std::string& phase
    ) const;
    NativeIonSolvationBorn evaluate_ion_solvation_born(
        double temperature,
        double pressure,
        double born_diameter_angstrom
    ) const;
    NativeAqueousMiacKij evaluate_aqueous_miac_kij(
        double temperature,
        double pressure,
        double formula_unit_molality_mol_per_kg,
        const std::array<double, 3>& k_ij
    ) const;
    std::vector<NativeAqueousMiacKij> evaluate_aqueous_miac_kij_batch(
        double temperature,
        double pressure,
        const std::vector<double>& formula_unit_molalities_mol_per_kg,
        const std::array<double, 3>& k_ij
    ) const;
    NativeAqueousMiacSolvationFactor evaluate_aqueous_miac_solvation_factor(
        double temperature,
        double pressure,
        double formula_unit_molality_mol_per_kg,
        double solvation_factor
    ) const;
    std::vector<NativeAqueousMiacSolvationFactor>
    evaluate_aqueous_miac_solvation_factor_batch(
        double temperature,
        double pressure,
        const std::vector<double>& formula_unit_molalities_mol_per_kg,
        double solvation_factor
    ) const;

private:
    const NativeParameters parameters_;
};

}  // namespace epcsaft
