#pragma once

#include "model.hpp"
#include "phase_block.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
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
    std::string requested_phase;
    std::string certificate_kind;
    std::string local_search_outcome;
    std::size_t local_root_count;
    double reachability_log_density_limit;
    std::size_t reachability_scan_interval_count;
    std::size_t certified_stable_root_index;
    double certified_origin_molar_density;
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
    std::size_t eos_evaluation_count = 1;
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
    std::size_t eos_evaluation_count = 0;
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

[[nodiscard]] std::vector<std::string> active_parameter_identities(
    const std::vector<NativeActiveCoordinate>& coordinates
);

struct NativeActiveParameterConstraint final {
    std::string identity;
    double temperature_k;
    double offset;
    std::vector<double> coefficients;
    double upper;
};

struct NativeEosMetadata final {
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

inline constexpr const char* EOS_HELMHOLTZ_BASIS_ID =
    "A_over_RT_reference_amount:n_ref=1mol:rho_ref=1mol_per_m3";

class EosDomainError final : public std::domain_error {
public:
    using std::domain_error::domain_error;
};

struct NativeEosEvaluationCounts final {
    std::size_t state = 0;
    std::size_t volume_bound = 0;
    std::size_t packing = 0;

    [[nodiscard]] std::size_t total() const noexcept {
        return state + volume_bound + packing;
    }

    NativeEosEvaluationCounts& operator+=(
        const NativeEosEvaluationCounts& other
    ) noexcept {
        state += other.state;
        volume_bound += other.volume_bound;
        packing += other.packing;
        return *this;
    }
};

struct NativePhaseState final {
    double value = 0.0;
    std::vector<double> gradient;
    std::vector<double> hessian;
    double pressure_pa = 0.0;
    double packing_fraction = 0.0;
    std::string mechanical_class;
    std::vector<std::size_t> supported_component_indices;
    std::string topology_fingerprint;
    NativeEosEvaluationCounts evaluations;
};

struct NativeParameterizedPhaseState final {
    NativeMixturePhaseBlock phase;
    std::vector<std::size_t> supported_component_indices;
    std::string topology_fingerprint;
    std::vector<double> parameter_derivatives;
    std::vector<double> state_parameter_derivatives;
    std::vector<double> pressure_parameter_derivatives_pa;
    std::vector<double> chemical_potential_parameter_derivatives_over_rt;
    NativeEosEvaluationCounts evaluations;
};

struct NativeParameterizedPackingFraction final {
    NativePackingFractionBlock packing;
    std::vector<double> parameter_derivatives;
    NativeEosEvaluationCounts evaluations;
};

struct NativeTemperaturePhaseDerivatives final {
    std::vector<double> chemical_potential_temperature_derivatives_over_rt_per_k;
    double pressure_temperature_derivative_pa_per_k = 0.0;
    NativeEosEvaluationCounts evaluations;
};

struct NativeInversePackingGeometryResult final {
    double volume_m3 = 0.0;
    std::vector<double> gradient;
    std::vector<double> hessian;
    NativeEosEvaluationCounts evaluations;
};

enum NativeSourceReferenceStatus {
    SOURCE_REFERENCE_OK = 0,
    SOURCE_REFERENCE_UNSUPPORTED = 1,
    SOURCE_REFERENCE_INVALID_SOURCE = 2,
    SOURCE_REFERENCE_COMPONENT_MISMATCH = 3,
    SOURCE_REFERENCE_MISSING_MODEL = 4,
    SOURCE_REFERENCE_DOMAIN_REJECTION = 5,
    SOURCE_REFERENCE_LIMIT_FAILURE = 6,
    SOURCE_REFERENCE_DERIVATIVE_UNAVAILABLE = 7,
    SOURCE_REFERENCE_NONFINITE = 8,
    SOURCE_REFERENCE_IDENTITY_MISMATCH = 9,
    SOURCE_REFERENCE_EVALUATION_ERROR = 10,
};

struct NativeSourceReferenceRequest final {
    std::vector<std::string> component_ids;
    std::vector<double> source_solvent_mole_fractions;
    std::vector<std::array<std::size_t, 2>> ion_pair_indices;
    std::vector<std::array<int, 2>> ion_pair_stoichiometry;
    std::vector<std::size_t> pure_component_indices;
    std::vector<std::string> pure_component_phases;
    std::string phase;
    std::string reference_convention;
    std::string reference_state_id;
    std::string activity_convention_id;
    std::string reference_state_fingerprint;
    double standard_molality_mol_per_kg;
    double source_reference_pressure_pa;
    std::uint32_t required_derivatives;
};

struct NativeSourceReferenceResult final {
    int status = SOURCE_REFERENCE_OK;
    std::string error;
    std::size_t component_count = 0;
    std::size_t neutral_basis_row_count = 0;
    std::vector<double> neutral_basis;
    std::vector<double> log_fugacity_contractions;
    std::vector<double> activity_scale_log_contractions;
    std::vector<double> transfer_log_contractions;
    std::vector<double> source_composition;
    std::uint32_t derivative_availability = 0;
    std::vector<double> pressure_derivatives_per_pa;
    std::vector<double> parameter_derivatives;
    double temperature_k = 0.0;
    double pressure_pa = 0.0;
    double source_reference_pressure_pa = 0.0;
    double standard_molality_mol_per_kg = 0.0;
    double solvent_molar_mass_kg_per_mol = 0.0;
    double reference_limit_molality_mol_per_kg = 0.0;
    double reference_convergence_error = 0.0;
    double source_temperature_min_k = 0.0;
    double source_temperature_max_k = 0.0;
    double source_pressure_min_pa = 0.0;
    double source_pressure_max_pa = 0.0;
    double candidate_temperature_min_k = 0.0;
    double candidate_temperature_max_k = 0.0;
    double candidate_pressure_min_pa = 0.0;
    double candidate_pressure_max_pa = 0.0;
    std::string candidate_qualification;
    std::string candidate_source_fingerprint;
    std::string parameter_fingerprint;
    std::string topology_fingerprint;
    std::string component_order_fingerprint;
    std::string reference_state_fingerprint;
    std::string domain_fingerprint;
    std::string artifact_fingerprint;
    std::string helmholtz_basis_id;
    NativeEosEvaluationCounts evaluations;
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
        double reference_molar_density,
        const std::string& phase,
        std::size_t certified_stable_root_index,
        double certified_origin_molar_density
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
    NativeModel with_eos_metadata(
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
    const NativeEosMetadata& eos_metadata() const;
    [[nodiscard]] const std::string& parameter_fingerprint() const;
    [[nodiscard]] std::vector<std::string> active_parameter_identities() const;
    [[nodiscard]] std::string parameterized_topology_fingerprint() const;
    [[nodiscard]] NativePhaseState evaluate_phase_state(
        double temperature,
        const std::vector<double>& amounts,
        double volume
    ) const;
    [[nodiscard]] NativeMolarVolumeBoundsParameterBlock evaluate_molar_volume_bounds(
        double temperature,
        const std::vector<double>& mole_fractions,
        double packing_fraction_min,
        double packing_fraction_max
    ) const;
    [[nodiscard]] NativePackingFractionBlock evaluate_state_packing_fraction(
        double temperature,
        const std::vector<double>& amounts,
        double volume
    ) const;
    [[nodiscard]] NativeParameterizedPhaseState evaluate_parameterized_phase_state(
        double temperature,
        const std::vector<double>& amounts,
        double volume,
        const std::vector<NativeActiveCoordinate>& active_parameters
    ) const;
    [[nodiscard]] NativeParameterizedPackingFraction
    evaluate_parameterized_packing_fraction(
        double temperature,
        const std::vector<double>& amounts,
        double volume,
        const std::vector<NativeActiveCoordinate>& active_parameters
    ) const;
    [[nodiscard]] NativeTemperaturePhaseDerivatives evaluate_temperature_phase(
        double temperature,
        const std::vector<double>& amounts,
        double volume
    ) const;
    [[nodiscard]] NativeInversePackingGeometryResult evaluate_state_inverse_packing_geometry(
        double temperature,
        const std::vector<double>& amounts,
        double log_packing_fraction,
        const std::vector<NativeActiveCoordinate>& active_parameters
    ) const;
    [[nodiscard]] NativeSourceReferenceResult evaluate_source_reference_transfer(
        double temperature,
        double pressure,
        const NativeSourceReferenceRequest& request
    ) const;
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
    NativeMixturePhaseBlock evaluate_supported_phase_derivatives(
        double temperature,
        const std::vector<double>& amounts,
        double volume,
        double ion_mole_fraction_max,
        const std::vector<NativeActiveCoordinate>& active_parameters,
        const std::vector<std::size_t>& supported_component_indices
    ) const;
    FixedPressureCompositionProjection
    evaluate_fixed_pressure_composition_derivatives(
        double temperature,
        const std::vector<double>& mole_fractions,
        double molar_volume
    ) const;
    NativeMixturePhaseBlock evaluate_temperature_phase_derivatives(
        double temperature,
        const std::vector<double>& amounts,
        double volume,
        double ion_mole_fraction_max,
        const std::vector<std::size_t>& supported_component_indices
    ) const;
    std::array<double, 2> molar_volume_bounds(
        double temperature,
        const std::vector<double>& mole_fractions,
        double packing_fraction_min,
        double packing_fraction_max
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
        const std::vector<NativeActiveCoordinate>& active_parameters,
        const std::vector<std::size_t>& supported_component_indices
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
    void validate_temperature(double temperature) const;
    void validate_pressure(double pressure) const;
    void validate_composition(const std::vector<double>& mole_fractions) const;
    [[nodiscard]] std::vector<std::size_t> validate_amounts(
        const std::vector<double>& amounts,
        double volume
    ) const;
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
    std::optional<NativeEosMetadata> eos_metadata_;
    std::vector<NativeActiveCoordinate> active_coordinates_;
};

}  // namespace epcsaft
