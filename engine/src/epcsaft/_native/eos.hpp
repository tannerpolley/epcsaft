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

private:
    const NativeParameters parameters_;
};

}  // namespace epcsaft
