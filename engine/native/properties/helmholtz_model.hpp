#pragma once

#include "core/convergence.hpp"
#include "properties/helmholtz_tape.hpp"
#include "properties/thermochemistry.hpp"

#include <array>
#include <cstddef>
#include <mutex>
#include <optional>
#include <vector>

namespace epcsaft {

using ParameterFamily = detail::ParameterFamily;
using ActiveParameter = detail::ActiveParameter;

enum class DerivativeOrder { Value, Gradient, PressureSlope, Hessian };

enum class NeutralProperty { ResidualHelmholtz, Pressure };

struct NeutralPropertyDerivatives final {
    double value = 0.0;
    std::vector<double> gradient;
    std::vector<double> hessian;
};

struct HelmholtzDerivatives final {
    double value = 0.0;
    std::vector<double> gradient;
    std::vector<double> hessian;
    double volume_second_derivative = 0.0;
    // Per-order association residual floor; zero when no association sites are present.
    double gradient_floor = 0.0;
    std::vector<double> site_fractions;
};

struct SweepCounters final {
    std::size_t forward = 0;
    std::size_t reverse = 0;
    std::size_t hessian_calls = 0;
    std::size_t recordings = 0;
    std::size_t directional_paths = 0;
    std::size_t maximum_energy_order = 0;
};

struct ParameterActionDirection final {
    ParameterFamily family;
    std::vector<std::size_t> indices;
    double value = 0.0;
};
enum class PropertyObservable {
    IdealEnthalpy, IdealInternalEnergy, IdealEntropy, IdealIsobaricHeatCapacity,
    IdealIsochoricHeatCapacity, ResidualEnthalpy, ResidualInternalEnergy,
    ResidualEntropy, ResidualIsobaricHeatCapacity, ResidualIsochoricHeatCapacity,
    TotalEnthalpy, TotalInternalEnergy, TotalEntropy, TotalIsobaricHeatCapacity,
    TotalIsochoricHeatCapacity
};
enum class PropertyActionStatus {
    Available, MissingThermochemistry, OutsideInterval, ThermalBoundary,
    SupportChange, IonicInsertion, AssociationSingular, AssociationIllConditioned,
    AssociationResidual, DensityUnavailable, PressureSingular,
    CentralPropertyUnavailable, Nonfinite
};
struct DensityActionDirection final {
    double temperature = 0.0, molar_density = 0.0;
    std::vector<double> composition;
    std::vector<ParameterActionDirection> parameter;
};
struct PressureActionDirection final {
    double temperature = 0.0, pressure = 0.0;
    std::vector<double> composition;
    std::vector<ParameterActionDirection> parameter;
};
struct PropertyActionResult final {
    PropertyObservable observable;
    std::optional<double> value, action;
    PropertyActionStatus status = PropertyActionStatus::Nonfinite;
    std::size_t energy_order = 0;
};

enum class DensityBranch { Unspecified, Vapor, Liquid };
struct PropertyActionDiagnostics final {
    std::size_t site_rank = 0, site_factorizations = 0, site_solves = 0;
    double site_pivot_ratio = 0.0, association_residual = 0.0;
    double pressure_density_derivative = 0.0, pressure_residual = 0.0;
    double stability_determinant = 0.0;
    std::size_t pressure_factorizations = 0, pressure_solves = 0;
    std::size_t thermal_evaluations = 0, recordings = 0, directional_paths = 0;
    std::size_t forward_sweeps = 0, dense_hessian_calls = 0;
    DensityBranch branch = DensityBranch::Unspecified;
    std::optional<std::size_t> stable_root_count;
};
struct PropertyActionBatch final {
    std::vector<PropertyActionResult> results;
    std::size_t maximum_energy_order = 0;
    PropertyActionDiagnostics diagnostics;
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
    std::optional<double> ideal_enthalpy;
    std::optional<double> ideal_internal_energy;
    std::optional<double> ideal_entropy;
    std::optional<double> ideal_isobaric_heat_capacity;
    std::optional<double> ideal_isochoric_heat_capacity;
    std::optional<double> total_enthalpy;
    std::optional<double> total_internal_energy;
    std::optional<double> total_entropy;
    std::optional<double> total_isobaric_heat_capacity;
    std::optional<double> total_isochoric_heat_capacity;
    ThermochemistryStatus thermochemistry_status = ThermochemistryStatus::MissingRecord;
    double molar_density = 0.0;
    double pressure_density_derivative = 0.0;
    std::vector<double> site_fractions;
    std::optional<core::Convergence> density_convergence;
    DensityBranch density_branch = DensityBranch::Unspecified;
    std::optional<std::size_t> stable_root_count;
};

struct DensityRoot final {
    double molar_density = 0.0;
    double pressure_residual = 0.0;
    double pressure_density_derivative = 0.0;
    std::size_t iterations = 0;
    std::optional<std::size_t> stable_root_count;
    DensityBranch branch = DensityBranch::Unspecified;
    core::Convergence convergence;
};

class HelmholtzModel final {
public:
    explicit HelmholtzModel(
        Parameters parameters,
        std::vector<ActiveParameter> active = {},
        std::optional<ThermochemistryRecord> thermochemistry = std::nullopt
    );

    HelmholtzDerivatives evaluate(
        double temperature,
        double volume,
        const std::vector<double>& amounts,
        DerivativeOrder order
    ) const;
    double third_directional(double, double, const std::vector<double>&, const std::array<std::vector<double>, 3>&) const;
    NeutralPropertyDerivatives neutral_property(
        double temperature,
        double molar_density,
        NeutralProperty property,
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
    PropertyActionBatch density_actions(
        double, double, const std::vector<double>&,
        const std::vector<PropertyObservable>&,
        const std::vector<DensityActionDirection>&
    ) const;
    PropertyActionBatch pressure_actions(
        double, double, const std::vector<double>&, DensityBranch,
        const std::vector<PropertyObservable>&,
        const std::vector<PressureActionDirection>&,
        std::optional<double> = std::nullopt
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
    const std::vector<ActiveParameter>& active_parameters() const { return active_; }

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
        const std::vector<double>& composition,
        bool residual = true, bool ideal = false
    ) const;
    PropertyActionBatch directional_actions(double, double, const std::vector<double>&, const std::vector<PropertyObservable>&, const std::vector<DensityActionDirection>&, const DensityRoot* = nullptr, const std::vector<double>& = {}) const;
    HelmholtzDerivatives evaluate_locked(
        double temperature,
        double volume,
        const std::vector<double>& amounts,
        DerivativeOrder order
    ) const;

    Parameters parameters_;
    std::vector<ActiveParameter> active_;
    std::optional<ThermochemistryRecord> thermochemistry_;
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
