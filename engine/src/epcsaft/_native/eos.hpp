#pragma once

#include "model.hpp"

#include <array>
#include <atomic>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <stdexcept>
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
    std::size_t stable_root_count;
    std::size_t selected_stable_root_index;
};

struct NativeResult final {
    double hard_chain;
    double dispersion;
    double association;
    double ares_polar;
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
struct NativeFixedTopologyParameter;
struct NativeFixedTopologyPhaseBlock;
struct NativeMixturePhaseBlock;
struct NativePackingFractionBlock;
struct NativeInversePackingGeometryBlock;

struct NativeReferenceSequence final {
    std::array<double, 2> log_fugacity;
    std::optional<double> pure_solvent_log_fugacity_coefficient;
    double reference_molality;
    double convergence_error;
};

struct NativeIonSolvationDerivativeBlock final {
    std::array<double, 2> solvation_gibbs_j_per_mol;
    std::vector<double> parameter_derivatives_j_per_mol;
    double reference_molality_mol_per_kg;
    double reference_convergence_error;
};

struct NativeIonSolvationKij final {
    double solvation_gibbs_j_per_mol;
    std::array<double, 3> derivative;
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

struct NativeIonFractionSuppression final {
    double bulk_relative_permittivity;
    double derivative;
    double salt_free_relative_permittivity;
};

class NativeEvaluationCancelled final : public std::runtime_error {
public:
    NativeEvaluationCancelled()
        : std::runtime_error("aqueous MIAC evaluation work budget exhausted") {}
};

class NativeEvaluationDeadlineExceeded final : public std::runtime_error {
public:
    NativeEvaluationDeadlineExceeded()
        : std::runtime_error("aqueous MIAC evaluation deadline exceeded") {}
};

class NativeEvaluationBudget final {
public:
    NativeEvaluationBudget(std::uint64_t work_budget, std::uint64_t deadline_ns)
        : remaining_work_(work_budget), deadline_ns_(deadline_ns) {}

    void checkpoint() {
        if (deadline_ns_ != 0 && monotonic_now_ns() >= deadline_ns_) {
            throw NativeEvaluationDeadlineExceeded();
        }
        if (remaining_work_.load(std::memory_order_relaxed)
            == std::numeric_limits<std::uint64_t>::max()) {
            return;
        }
        std::uint64_t remaining = remaining_work_.load(
            std::memory_order_relaxed
        );
        while (remaining != 0) {
            if (remaining_work_.compare_exchange_weak(
                    remaining,
                    remaining - 1,
                    std::memory_order_relaxed,
                    std::memory_order_relaxed
                )) {
                return;
            }
        }
        throw NativeEvaluationCancelled();
    }

private:
    static std::uint64_t monotonic_now_ns() {
        const auto now = std::chrono::steady_clock::now().time_since_epoch();
        return static_cast<std::uint64_t>(
            std::chrono::duration_cast<std::chrono::nanoseconds>(now).count()
        );
    }

    std::atomic<std::uint64_t> remaining_work_;
    const std::uint64_t deadline_ns_;
};

enum class ActiveComponentParameter {
    SegmentDiameter,
    DispersionEnergyOverK,
};

struct NativeActiveParameter final {
    ActiveComponentParameter family;
    std::size_t component;
    double value;
};

struct NativeFugacityDerivativeBlock final {
    std::vector<double> log_fugacity_coefficient;
    std::vector<double> pressure_derivative_per_pa;
    std::vector<double> parameter_derivatives;
    double molar_density_mol_per_m3;
    double root_residual_pa;
    double pressure_density_derivative_pa_m3_per_mol;
    double density_condition_number;
    std::string branch;
    double bracket_lower_mol_per_m3;
    double bracket_upper_mol_per_m3;
    std::size_t stable_root_count;
    std::size_t selected_stable_root_index;
};

struct NativeMolarVolumeBoundsParameterBlock final {
    std::array<double, 2> bounds_m3_per_mol;
    std::vector<double> parameter_derivatives_m3_per_mol;
};

// Canonical EOS-facing request.  The native resolver translates identities to
// the existing typed derivative primitives; Python never supplies a native
// array slot or callback.
struct NativeParameterRecord final {
    std::string identity;
    double value;
};

struct NativePairCoordinate final {
    std::string identity;
    std::string family;
    std::size_t component_a;
    std::size_t component_b;
    std::string association_role_a;
    std::string association_role_b;
    std::vector<std::array<std::size_t, 2>> association_edges;
    std::vector<double> association_energy_scales;
    double value;
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
    bool supports_pure_transport() const;
    bool supports_fixed_topology_association() const;
    bool supports_fixed_topology_associating_mixture() const;
    bool supports_mixture_transport() const;
    bool supports_electrolyte_transport() const;
    std::size_t component_count() const;
    const std::vector<int>& charge_numbers() const;
    const std::vector<int>& site_components() const;
    const std::vector<int>& site_multiplicities() const;
    const std::vector<int>& association_active_pairs() const;
    bool has_temperature_independent_segment_diameter(
        std::size_t component
    ) const;
    const std::vector<double>& molar_masses_kg_per_mol() const;
    NativePurePhaseBlock evaluate_pure_phase(
        double temperature,
        double amount,
        double volume
    ) const;
    NativeModel with_parameter_identities(
        const std::vector<std::string>& component_ids,
        const std::vector<std::string>& site_ids
    ) const;
    NativeModel with_active_parameters(
        const std::vector<NativeParameterRecord>& active_parameters
    ) const;
    NativeFixedTopologyPhaseBlock evaluate_fixed_topology_association(
        double temperature,
        double amount,
        double volume,
        const std::vector<NativeFixedTopologyParameter>& active_parameters
    ) const;
    NativeFixedTopologyPhaseBlock
    evaluate_fixed_topology_association_active_parameters(
        double temperature,
        double amount,
        double volume,
        const std::vector<NativeParameterRecord>& active_parameters
    ) const;
    NativeMixturePhaseBlock evaluate_fixed_topology_associating_mixture(
        double temperature,
        const std::vector<double>& amounts,
        double volume,
        const std::vector<NativeFixedTopologyParameter>& active_parameters
    ) const;
    NativeMixturePhaseBlock
    evaluate_fixed_topology_associating_mixture_active_parameters(
        double temperature,
        const std::vector<double>& amounts,
        double volume,
        const std::vector<NativeParameterRecord>& active_parameters
    ) const;
    NativeMixturePhaseBlock evaluate_mixture_phase(
        double temperature,
        const std::vector<double>& amounts,
        double volume
    ) const;
    NativeMixturePhaseBlock evaluate_mixture_phase_pair_coordinates(
        double temperature,
        const std::vector<double>& amounts,
        double volume,
        const std::vector<NativePairCoordinate>& coordinates
    ) const;
    NativeMixturePhaseBlock evaluate_electrolyte_phase(
        double temperature,
        const std::vector<double>& amounts,
        double volume,
        double ion_mole_fraction_max
    ) const;
    NativeMixturePhaseBlock evaluate_reacting_phase_parameters(
        double temperature,
        const std::vector<double>& amounts,
        double volume,
        double ion_mole_fraction_max,
        const std::vector<NativeActiveParameter>& active_parameters
    ) const;
    double evaluate_electrolyte_phase_value(
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
    double molar_density_upper_bound(
        double temperature,
        const std::vector<double>& mole_fractions
    ) const;
    NativePackingFractionBlock evaluate_packing_fraction(
        double temperature,
        const std::vector<double>& amounts,
        double volume,
        double ion_mole_fraction_max
    ) const;
    NativePackingFractionBlock evaluate_packing_fraction_parameters(
        double temperature,
        const std::vector<double>& amounts,
        double volume,
        double ion_mole_fraction_max,
        const std::vector<NativeActiveParameter>& active_parameters
    ) const;
    NativeInversePackingGeometryBlock evaluate_inverse_packing_geometry(
        double temperature,
        const std::vector<double>& amounts,
        double log_packing_fraction,
        const std::vector<NativeActiveParameter>& active_parameters
    ) const;
    NativeMolarVolumeBoundsParameterBlock molar_volume_bounds_parameters(
        double temperature,
        const std::vector<double>& mole_fractions,
        double packing_fraction_min,
        double packing_fraction_max,
        double ion_mole_fraction_max,
        const std::vector<NativeActiveParameter>& active_parameters
    ) const;
    NativeFugacityDerivativeBlock evaluate_fugacity_derivatives(
        double temperature,
        double pressure,
        const std::vector<double>& composition,
        const std::string& phase,
        const std::vector<NativeActiveParameter>& active_parameters
    ) const;
    NativeFugacityDerivativeBlock evaluate_fugacity_derivatives_active(
        double temperature,
        double pressure,
        const std::vector<double>& composition,
        const std::string& phase,
        const std::vector<NativeParameterRecord>& active_parameters
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
        const std::string& phase,
        NativeEvaluationBudget* budget = nullptr
    ) const;
    NativeIonSolvationDerivativeBlock evaluate_ion_solvation_derivatives(
        double temperature,
        double pressure,
        const std::vector<double>& solvent_fractions,
        const std::array<std::size_t, 2>& ion_indices,
        const std::array<int, 2>& stoichiometry,
        const std::string& phase,
        const std::vector<NativeParameterRecord>& active_parameters
    ) const;

    NativeIonSolvationKij evaluate_ion_solvation_kij(
        double temperature,
        double pressure,
        std::size_t active_ion_component,
        const std::array<double, 3>& k_ij
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
        const std::array<double, 3>& k_ij,
        NativeEvaluationBudget* budget = nullptr
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
        double solvation_factor,
        NativeEvaluationBudget* budget = nullptr
    ) const;
    NativeIonFractionSuppression evaluate_ion_fraction_suppression(
        double total_ion_mole_fraction,
        double ion_fraction_suppression_coefficient
    ) const;

private:
    NativeFixedTopologyParameter resolve_fixed_topology_parameter(
        const NativeParameterRecord& parameter
    ) const;
    NativeActiveParameter resolve_component_parameter(
        const NativeParameterRecord& parameter
    ) const;
    const NativeParameters parameters_;
    std::vector<std::string> component_ids_;
    std::vector<std::string> site_ids_;
};

}  // namespace epcsaft
