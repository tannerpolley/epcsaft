#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

#include "eos.hpp"
#include "phase_block.hpp"

namespace epcsaft_engine::equilibrium {

inline constexpr const char* PROVIDER_HELMHOLTZ_BASIS_ID =
    "A_over_RT_reference_amount:n_ref=1mol:rho_ref=1mol_per_m3";

class ProviderDomainError final : public std::domain_error {
public:
    using std::domain_error::domain_error;
};

struct PhaseEvaluation {
    double amount_mol = 0.0;
    double volume_m3 = 0.0;
    epcsaft::NativePurePhaseBlock provider{};
    std::string parameter_fingerprint;
};

struct MixturePhaseEvaluation {
    double value = 0.0;
    std::vector<double> gradient;
    std::vector<double> hessian;
    double pressure_pa = 0.0;
    std::string parameter_fingerprint;
};

struct PackingFractionEvaluation {
    double value = 0.0;
    std::vector<double> gradient;
    std::vector<double> hessian;
};

struct ProviderActiveParameterSet {
    std::vector<std::string> identities;
};

struct ParameterizedPhaseEvaluation {
    MixturePhaseEvaluation phase;
    PackingFractionEvaluation packing;
    std::array<double, 2> molar_volume_bounds_m3_per_mol{};
    std::vector<double> parameter_derivatives;
    std::vector<double> state_parameter_derivatives;
    std::vector<double> packing_parameter_derivatives;
    std::vector<double> molar_volume_bound_parameter_derivatives_m3_per_mol;
    std::vector<double> pressure_parameter_derivatives_pa;
    std::vector<double> chemical_potential_parameter_derivatives_over_rt;
    std::vector<std::size_t> supported_component_indices;
    std::string topology_fingerprint;
};

struct TemperaturePhaseEvaluation {
    std::vector<double> chemical_potential_temperature_derivatives_over_rt_per_k;
    double pressure_temperature_derivative_pa_per_k = 0.0;
};

struct InversePackingGeometryEvaluation {
    double volume_m3 = 0.0;
    std::vector<double> gradient;
    std::vector<double> hessian;
    std::vector<std::size_t> supported_component_indices;
    std::string topology_fingerprint;
};

enum SourceReferenceStatus {
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

struct SourceReferenceRequest final {
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

struct SourceReferenceResult final {
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
    std::string parameter_fingerprint;
    std::string topology_fingerprint;
    std::string component_order_fingerprint;
    std::string reference_state_fingerprint;
    std::string domain_fingerprint;
    std::string artifact_fingerprint;
    std::string helmholtz_basis_id;
};

class ProviderContext {
public:
    explicit ProviderContext(const epcsaft::NativeModel& model);

    [[nodiscard]] std::string parameterized_topology_fingerprint() const;

    [[nodiscard]] PhaseEvaluation evaluate(
        double temperature_k,
        double amount_mol,
        double volume_m3
    ) const;

    [[nodiscard]] MixturePhaseEvaluation evaluate_mixture(
        double temperature_k,
        const std::vector<double>& amounts_mol,
        double volume_m3
    ) const;

    [[nodiscard]] MixturePhaseEvaluation evaluate_electrolyte(
        double temperature_k,
        const std::vector<double>& amounts_mol,
        double volume_m3
    ) const;

    [[nodiscard]] std::array<double, 2> evaluate_molar_volume_bounds(
        double temperature_k,
        const std::vector<double>& mole_fractions,
        double packing_fraction_min,
        double packing_fraction_max
    ) const;

    [[nodiscard]] PackingFractionEvaluation evaluate_packing_fraction(
        double temperature_k,
        const std::vector<double>& amounts_mol,
        double volume_m3
    ) const;

    [[nodiscard]] ProviderActiveParameterSet active_parameters() const;

    [[nodiscard]] ParameterizedPhaseEvaluation evaluate_parameterized_phase(
        double temperature_k,
        const std::vector<double>& amounts_mol,
        double volume_m3,
        double packing_fraction_min,
        double packing_fraction_max,
        const ProviderActiveParameterSet& active_parameters
    ) const;
    [[nodiscard]] TemperaturePhaseEvaluation evaluate_temperature_phase(
        double temperature_k,
        const std::vector<double>& amounts_mol,
        double volume_m3
    ) const;

    [[nodiscard]] InversePackingGeometryEvaluation
    evaluate_inverse_packing_geometry(
        double temperature_k,
        const std::vector<double>& amounts_mol,
        double log_packing_fraction,
        const ProviderActiveParameterSet& active_parameters
    ) const;

    [[nodiscard]] SourceReferenceResult evaluate_source_reference_transfer(
        double temperature_k,
        double pressure_pa,
        const SourceReferenceRequest& request
    ) const;

    [[nodiscard]] const std::string& fingerprint() const;
    [[nodiscard]] const epcsaft::NativeModel& model() const noexcept;
    [[nodiscard]] const epcsaft::NativeProviderMetadata& metadata() const;
    [[nodiscard]] std::size_t component_count() const noexcept;
    [[nodiscard]] const std::vector<std::string>& component_ids() const noexcept;
    [[nodiscard]] const std::vector<int>& component_charges() const noexcept;
    [[nodiscard]] bool supports_pure_phase() const noexcept;
    [[nodiscard]] bool supports_mixture_phase() const noexcept;
    [[nodiscard]] bool supports_electrolyte_phase() const noexcept;
    [[nodiscard]] int evaluation_count() const noexcept;

private:
    epcsaft::NativeModel model_;
    std::string fingerprint_;
    mutable int evaluation_count_ = 0;
};

}  // namespace epcsaft_engine::equilibrium
