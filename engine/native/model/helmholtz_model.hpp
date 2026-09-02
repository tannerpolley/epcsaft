#pragma once

#include "model/helmholtz_tape.hpp"

#include <array>
#include <cstddef>
#include <mutex>
#include <optional>
#include <vector>

namespace epcsaft {

using ParameterFamily = detail::ParameterFamily;
using ActiveParameter = detail::ActiveParameter;

enum class DerivativeOrder { Value, Gradient, PressureSlope, Hessian };

struct HelmholtzDerivatives final {
    double value = 0.0;
    std::vector<double> gradient;
    std::vector<double> hessian;
    double volume_second_derivative = 0.0;
    std::vector<double> site_fractions;
};

struct SweepCounters final {
    std::size_t forward = 0;
    std::size_t reverse = 0;
    std::size_t hessian_calls = 0;
    std::size_t recordings = 0;
};

struct State final {
    double hard_chain = 0.0;
    double dispersion = 0.0;
    double association = 0.0;
    double debye_huckel = 0.0;
    double born = 0.0;
    double bulk_relative_permittivity = 0.0;
    double residual_helmholtz = 0.0;
    double compressibility_factor = 0.0;
    double pressure = 0.0;
    std::vector<double> residual_chemical_potential_over_rt;
    std::optional<std::vector<double>> log_fugacity_coefficient;
    double residual_entropy = 0.0;
    double residual_internal_energy = 0.0;
    double residual_enthalpy = 0.0;
    double residual_isochoric_heat_capacity = 0.0;
    std::optional<double> residual_isobaric_heat_capacity;
    double molar_density = 0.0;
    double pressure_density_derivative = 0.0;
    std::vector<double> site_fractions;
};

enum class DensityBranch { Unspecified, Vapor, Liquid };

struct DensityRoot final {
    double molar_density = 0.0;
    double pressure_residual = 0.0;
    double pressure_density_derivative = 0.0;
    std::size_t iterations = 0;
    std::size_t stable_root_count = 0;
    DensityBranch branch = DensityBranch::Unspecified;
};

class HelmholtzModel final {
public:
    explicit HelmholtzModel(
        Parameters parameters,
        std::vector<ActiveParameter> active = {}
    );

    HelmholtzDerivatives evaluate(
        double temperature,
        double volume,
        const std::vector<double>& amounts,
        DerivativeOrder order
    ) const;
    State state(
        double temperature,
        double molar_density,
        const std::vector<double>& mole_fractions
    ) const;
    State state_at_pressure(
        double temperature,
        double pressure,
        const std::vector<double>& mole_fractions,
        DensityBranch branch,
        std::optional<double> anchor_density = std::nullopt
    ) const;
    double packing_fraction_coefficient(
        double temperature,
        const std::vector<double>& mole_fractions
    ) const;
    double phase_packing_fraction(
        double temperature,
        const std::vector<double>& amounts,
        double volume
    ) const;
    double molar_density_upper_bound(
        double temperature,
        const std::vector<double>& mole_fractions
    ) const;
    SweepCounters counters() const;
    void reset_counters() const;
    const Parameters& parameters() const { return parameters_; }

private:
    friend DensityRoot solve_density(
        const HelmholtzModel&, double, double, const std::vector<double>&,
        DensityBranch, std::optional<double>
    );

    void validate_point(
        double temperature,
        double volume,
        const std::vector<double>& amounts
    ) const;
    void validate_composition(const std::vector<double>& composition) const;
    detail::PredicateSignature predicate_signature(
        double temperature,
        double density,
        const std::vector<double>& composition
    ) const;

    Parameters parameters_;
    std::vector<ActiveParameter> active_;
    mutable std::unique_ptr<detail::HelmholtzTape> tape_;
    mutable SweepCounters counters_;
    mutable std::mutex mutex_;
};

DensityRoot solve_density(
    const HelmholtzModel& model,
    double temperature,
    double pressure,
    const std::vector<double>& mole_fractions,
    DensityBranch branch,
    std::optional<double> anchor_density
);

struct InfiniteDilutionReference final {
    std::array<double, 2> log_fugacity{};
    double molality = 0.0;
    double convergence_error = 0.0;
    std::optional<double> pure_solvent_log_fugacity;
};

struct ElectrolyteActivity final {
    std::array<double, 2> single_ion{};
    double mean_ionic_mole_fraction_scale = 0.0;
    double mean_ionic_log_molality_scale = 0.0;
    double mean_ionic_molality_scale = 0.0;
    std::optional<double> osmotic_coefficient;
};

std::vector<double> formula_composition(
    const std::vector<double>& solvent_fractions,
    const std::vector<double>& molar_masses,
    const std::array<std::size_t, 2>& ion_indices,
    const std::array<int, 2>& stoichiometry,
    double molality
);

InfiniteDilutionReference infinite_dilution_reference(
    const HelmholtzModel& model,
    double temperature,
    double pressure,
    const std::vector<double>& solvent_fractions,
    const std::array<std::size_t, 2>& ion_indices,
    const std::array<int, 2>& stoichiometry
);

ElectrolyteActivity electrolyte_activity(
    const HelmholtzModel& model,
    double temperature,
    double pressure,
    const std::vector<double>& solvent_fractions,
    const std::array<std::size_t, 2>& ion_indices,
    const std::array<int, 2>& stoichiometry,
    double molality,
    const InfiniteDilutionReference& reference
);

}  // namespace epcsaft
