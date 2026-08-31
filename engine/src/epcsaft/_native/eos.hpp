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
    std::size_t stable_root_count;
    std::size_t selected_stable_root_index;
};

struct NativeResult final {
    double hard_chain;
    double dispersion;
    double association;
    double debye_huckel;
    double born;
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

struct NativeMolarVolumeBoundsParameterBlock final {
    std::array<double, 2> bounds_m3_per_mol;
    std::vector<double> parameter_derivatives_m3_per_mol;
};

struct FixedPressureCompositionProjection final {
    std::string failure_code;
    std::vector<std::size_t> coordinate_component_indices;
    std::vector<std::size_t> dependent_component_indices;
    std::vector<double> composition_derivatives;
    std::vector<double> fixed_volume_chemical_potential_derivatives_over_rt;
    std::vector<double> chemical_potential_derivatives_over_rt;
    std::vector<double> molar_density_derivatives;
    std::vector<double> fixed_volume_pressure_derivatives_pa;
    double pressure_closure_condition_measure = 0.0;
};

FixedPressureCompositionProjection fixed_pressure_composition_projection(
    double temperature,
    const std::vector<double>& mole_fractions,
    double molar_volume,
    const std::vector<int>& charge_numbers,
    const std::vector<double>& hessian
);

struct NativeActiveCoordinate final {
    std::string identity;
    std::string family;
    std::string application;
    std::vector<std::size_t> indices;
    std::vector<std::array<std::size_t, 2>> association_edges;
    std::vector<double> association_energy_scales;
    double value;
    double reference_temperature_k;
};

struct NativeActiveParameterConstraint final {
    std::string identity;
    double temperature_k;
    double offset;
    std::vector<double> coefficients;
    double upper;
};

struct NativeProviderMetadata final {
    std::string parameter_fingerprint;
    std::string topology_fingerprint;
    std::string component_order_fingerprint;
    std::string domain_fingerprint;
    std::string artifact_fingerprint;
    std::string candidate_qualification;
    std::string candidate_source_fingerprint;
    double temperature_min_k;
    double temperature_max_k;
    double pressure_min_pa;
    double pressure_max_pa;
    double ion_mole_fraction_max;
    double source_temperature_min_k;
    double source_temperature_max_k;
    double source_pressure_min_pa;
    double source_pressure_max_pa;
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
    NativeResult evaluate_pressure_continuation(
        double temperature,
        double pressure,
        const std::vector<double>& mole_fractions,
        double reference_molar_density
    ) const;
    bool supports_pure_transport() const;
    bool supports_fixed_topology_associating_mixture() const;
    bool supports_mixture_transport() const;
    bool supports_electrolyte_transport() const;
    std::size_t component_count() const;
    const std::vector<int>& charge_numbers() const;
    const std::vector<int>& site_components() const;
    const std::vector<int>& site_multiplicities() const;
    const std::vector<int>& association_active_pairs() const;
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
    NativeModel with_provider_metadata(
        std::string parameter_fingerprint,
        double temperature_min_k,
        double temperature_max_k,
        double ion_mole_fraction_max,
        double pressure_min_pa,
        double pressure_max_pa,
        std::string topology_fingerprint,
        std::string component_order_fingerprint,
        std::string domain_fingerprint,
        std::string artifact_fingerprint,
        double source_temperature_min_k,
        double source_temperature_max_k,
        double source_pressure_min_pa,
        double source_pressure_max_pa,
        std::string candidate_qualification,
        std::string candidate_source_fingerprint
    ) const;
    NativeModel with_active_coordinates(
        const std::vector<NativeActiveCoordinate>& active_coordinates
    ) const;
    const NativeProviderMetadata& provider_metadata() const;
    const std::vector<std::string>& component_ids() const;
    const std::vector<NativeActiveCoordinate>& active_coordinates() const;
    std::vector<NativeActiveParameterConstraint>
    active_parameter_constraints() const;
    NativeFixedTopologyPhaseBlock evaluate_fixed_topology_association(
        double temperature,
        double amount,
        double volume,
        const std::vector<NativeFixedTopologyParameter>& active_parameters
    ) const;
    NativeMixturePhaseBlock evaluate_mixture_phase(
        double temperature,
        const std::vector<double>& amounts,
        double volume
    ) const;
    NativeMixturePhaseBlock evaluate_phase_derivatives(
        double temperature,
        const std::vector<double>& amounts,
        double volume,
        const std::vector<NativeActiveCoordinate>& coordinates,
        std::size_t derivative_order
    ) const;
    NativeMixturePhaseBlock evaluate_electrolyte_phase(
        double temperature,
        const std::vector<double>& amounts,
        double volume,
        double ion_mole_fraction_max
    ) const;
    NativeMixturePhaseBlock evaluate_parameterized_phase(
        double temperature,
        const std::vector<double>& amounts,
        double volume,
        double ion_mole_fraction_max,
        const std::vector<NativeActiveCoordinate>& active_parameters
    ) const;
    NativeMixturePhaseBlock evaluate_supported_phase_derivatives(
        double temperature,
        const std::vector<double>& amounts,
        double volume,
        double ion_mole_fraction_max,
        const std::vector<NativeActiveCoordinate>& active_parameters
    ) const;
    FixedPressureCompositionProjection
    evaluate_fixed_pressure_composition_derivatives(
        double temperature,
        const std::vector<double>& mole_fractions,
        double molar_volume,
        double ion_mole_fraction_max
    ) const;
    NativeMixturePhaseBlock evaluate_temperature_phase_derivatives(
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
        const std::vector<NativeActiveCoordinate>& active_parameters
    ) const;
    NativeInversePackingGeometryBlock evaluate_inverse_packing_geometry(
        double temperature,
        const std::vector<double>& amounts,
        double log_packing_fraction,
        const std::vector<NativeActiveCoordinate>& active_parameters
    ) const;
    NativeMolarVolumeBoundsParameterBlock molar_volume_bounds_parameters(
        double temperature,
        const std::vector<double>& mole_fractions,
        double packing_fraction_min,
        double packing_fraction_max,
        double ion_mole_fraction_max,
        const std::vector<NativeActiveCoordinate>& active_parameters
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
private:
    NativeFixedTopologyParameter resolve_fixed_topology_parameter(
        const std::string& identity,
        double value
    ) const;
    NativeFixedTopologyParameter resolve_active_coordinate(
        const NativeActiveCoordinate& coordinate
    ) const;
    const NativeParameters parameters_;
    std::vector<std::string> component_ids_;
    std::vector<std::string> site_ids_;
    std::optional<NativeProviderMetadata> provider_metadata_;
    std::vector<NativeActiveCoordinate> active_coordinates_;
};

}  // namespace epcsaft
