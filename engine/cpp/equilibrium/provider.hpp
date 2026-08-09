#pragma once

#include <array>
#include <cstddef>
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
    explicit ProviderContext(const epcsaft_native_transport& transport);

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
