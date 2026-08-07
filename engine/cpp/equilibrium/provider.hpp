#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include "native_transport.hpp"

namespace epcsaft_engine::equilibrium {

struct PhaseEvaluation {
    double amount_mol = 0.0;
    double volume_m3 = 0.0;
    epcsaft_phase_block_result provider{};
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

struct ProviderActiveParameterInput {
    std::string family;
    std::string identity;
    std::vector<std::string> component_ids;
    double value = 0.0;
    std::string unit;
};

struct ProviderActiveParameter {
    std::string family;
    std::string identity;
    std::vector<std::string> component_ids;
    double value = 0.0;
    std::string unit;
    std::uint32_t family_code = 0;
    std::uint32_t identity_code = 0;
    std::int32_t component_index = -1;
    std::int32_t pair_component_index_a = -1;
    std::int32_t pair_component_index_b = -1;
};

struct ProviderActiveParameterSet {
    std::vector<ProviderActiveParameter> parameters;
    std::string topology_fingerprint;
};

struct ParameterizedPhaseEvaluation {
    MixturePhaseEvaluation phase;
    PackingFractionEvaluation packing;
    std::array<double, 2> molar_volume_bounds_m3_per_mol{};
    std::vector<double> state_parameter_derivatives;
    std::vector<double> pressure_parameter_derivatives_pa;
    std::vector<double> chemical_potential_parameter_derivatives_over_rt;
    std::string topology_fingerprint;
};

struct InversePackingGeometryEvaluation {
    double volume_m3 = 0.0;
    std::vector<double> gradient;
    std::vector<double> hessian;
    std::string topology_fingerprint;
};

class ProviderContext {
public:
    ProviderContext(const epcsaft_native_transport& transport, std::string fingerprint);

    [[nodiscard]] std::string reacting_topology_fingerprint() const;

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

    [[nodiscard]] ProviderActiveParameterSet resolve_active_parameters(
        double temperature_k,
        const std::vector<ProviderActiveParameterInput>& requests
    ) const;

    [[nodiscard]] ParameterizedPhaseEvaluation evaluate_reacting_phase_parameters(
        double temperature_k,
        const std::vector<double>& amounts_mol,
        double volume_m3,
        double packing_fraction_min,
        double packing_fraction_max,
        const ProviderActiveParameterSet& active_parameters
    ) const;

    [[nodiscard]] InversePackingGeometryEvaluation
    evaluate_inverse_packing_geometry(
        double temperature_k,
        const std::vector<double>& amounts_mol,
        double log_packing_fraction,
        const ProviderActiveParameterSet& active_parameters
    ) const;

    [[nodiscard]] const std::string& fingerprint() const;
    [[nodiscard]] const epcsaft_native_transport& transport() const noexcept;

private:
    const epcsaft_native_transport& transport_;
    std::string fingerprint_;
};

}  // namespace epcsaft_engine::equilibrium
