#pragma once

#include "core/convergence.hpp"
#include "properties/helmholtz_tape.hpp"
#include "properties/thermochemistry.hpp"

#include <array>
#include <cstddef>
#include <mutex>
#include <limits>
#include <optional>
#include <string>
#include <variant>
#include <vector>

namespace epcsaft {

using ParameterFamily = detail::ParameterFamily;
using ActiveParameter = detail::ActiveParameter;

// PressureSlope returns the reduced gradient, a_VV and the association floor, not a Hessian.
// Phase Hessians differentiate only V, positive amounts and (unless isothermal) T; other
// gradient/Hessian entries and the floor's coordinates are uncomputed padding, not physical zeros.
enum class DerivativeOrder { Value, Gradient, PressureSlope, Hessian, PhaseHessian, IsothermalPhaseHessian };

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
    // #86 pressure chart: root, first-order density derivatives per direction, status, and the root search's failure message and record.
    std::optional<double> molar_density;
    std::vector<double> density_action;
    PropertyActionStatus density_status = PropertyActionStatus::CentralPropertyUnavailable;
    std::string density_message;
    std::optional<core::Convergence> density_convergence;
};
struct ResidualEnergyAction final {
    std::optional<double> action; PropertyActionStatus status = PropertyActionStatus::Nonfinite; PropertyActionDiagnostics diagnostics;
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
    std::vector<double> log_fugacity_parameter_derivative;  // component x active, fixed T/P/x
    std::vector<double> density_parameter_derivative, permittivity_parameter_derivative;
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

// Fixed-composition pressure crossings, slope zeros and finite-evaluation brackets (ALG-DENSITY).
struct DensityDiagnosticPoint final {
    double molar_density = 0.0, lower = 0.0, upper = 0.0;
    double packing_fraction = std::numeric_limits<double>::quiet_NaN();
    double pressure = std::numeric_limits<double>::quiet_NaN();
    double pressure_density_derivative = std::numeric_limits<double>::quiet_NaN();
    double residual = std::numeric_limits<double>::quiet_NaN();
    int mechanical_sign = 0;
    bool linked = false, declared = true;
    std::string cause, stop;
};
struct DensityDiagnosticBracket final {
    double lower = 0.0, upper = 0.0, lower_residual = 0.0, upper_residual = 0.0;
};
struct DensityDiagnostics final {
    std::vector<DensityDiagnosticPoint> samples, pressure_roots, spinodals;
    std::vector<DensityDiagnosticBracket> declared_brackets;
    std::vector<DensityDiagnosticPoint> unresolved;
    std::optional<DensityDiagnosticPoint> finite_domain;
    DerivativeOrder evaluation_order = DerivativeOrder::PressureSlope;
    double geometric_cap = 0.0;
    std::size_t declared_crossings = 0, evaluations = 0;
    std::string stop;
};

// A Stage-A neutral reference is a pure, neutral solvent plus one explicit
// full-width solute molality-ratio vector.  Mixed-solvent declarations are
// intentionally outside this first native consumer contract.
struct NeutralReferenceDeclaration final {
    std::optional<std::size_t> solvent_index;
    std::vector<double> solute_molality_ratios;
    double standard_molality_mol_per_kg = 0.0;
};

enum class ReferenceSourceBasis {
    CommonMolalityInfiniteDilution,
    RawMoleFractionInfiniteDilution,
};

struct ReactionReference final {
    std::optional<double> reference_pressure_pa;
    ReferenceSourceBasis source_basis =
        ReferenceSourceBasis::CommonMolalityInfiniteDilution;
};

enum class NeutralReferenceStatus {
    Available,
    InvalidDeclaration,
    NonNeutralDeclaration,
    UnsupportedSourceBasis,
    FreeTemperature,
    FreePressure,
    DensityUnavailable,
    BranchChange,
    Nonfinite,
    TerminalChangeExceeded,
};

struct NeutralReferenceReactionDiagnostic final {
    std::size_t reaction_index = 0;
    std::optional<double> terminal_contraction;
    std::optional<double> observed_terminal_change;
    // Pressure derivatives of the finite terminal contraction.  These are native
    // reference-path diagnostics; explicit fixed-pressure sources carry exact zero
    // system-pressure derivatives.
    std::optional<double> terminal_contraction_pressure_derivative;
    std::optional<double> terminal_contraction_pressure_second_derivative;
    // dC/dT at the reference pressure, set only when requested and when T times its
    // last ladder change is within the value's terminal criterion.
    std::optional<double> terminal_contraction_temperature_derivative;
    std::optional<double> observed_temperature_derivative_change;
    // d2C/dT2 at the reference pressure, set only for temperature order two and when T^2
    // times its last ladder change is within the same criterion.
    std::optional<double> terminal_contraction_temperature_second_derivative;
    std::optional<double> observed_temperature_second_derivative_change;
    // d2C/dTdP, only with temperature and pressure jets; admitted like the first slope at scales T, P.
    std::optional<double> terminal_contraction_temperature_pressure_derivative;
    // dC/dtheta_a at the reference pressure per active parameter, only when requested; nullopt
    // where max(|theta_a|, 1) times its last ladder change exceeds the value's criterion.
    std::vector<std::optional<double>> terminal_contraction_parameter_derivative;
    std::vector<std::optional<double>> terminal_contraction_temperature_parameter_derivative;  // d2C/dTdtheta
    std::vector<std::array<std::optional<double>, 3>> temperature_parameter_ladder;  // last three C_Ttheta values per active coordinate
    std::vector<std::optional<double>> observed_temperature_parameter_change;  // T max(1, |theta|) times ladder change
    double reference_pressure_pa = 0.0;
    double solvent_density = 0.0;  // liquid density of the trace's pure-solvent point
    NeutralReferenceStatus status = NeutralReferenceStatus::Nonfinite;
    std::string message;
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
    // Mixed contraction of the extensive residual Helmholtz scalar with three or four
    // directions of width 2 + C + A in [T, V, n, active parameters]: E3 or E4.
    ResidualEnergyAction residual_energy_action(double, double, const std::vector<double>&, const std::vector<std::vector<double>>&) const;
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
    DensityDiagnostics density_diagnostics(
        double temperature, double pressure, const std::vector<double>& composition,
        const std::vector<double>& nodes, double density_tolerance, double log_tolerance,
        std::size_t max_evaluations, double max_seconds, bool adaptive,
        std::optional<double> finite_anchor, double finite_upper, DerivativeOrder evaluation_order
    ) const;
    SweepCounters counters() const;
    void reset_counters() const;
    const Parameters& parameters() const { return parameters_; }
    const std::vector<ActiveParameter>& active_parameters() const { return active_; }
    std::vector<double> active_parameter_values() const;
    const std::optional<ThermochemistryRecord>& thermochemistry() const { return thermochemistry_; }

private:
    friend DensityRoot solve_density(
        const HelmholtzModel&, double, double, const std::vector<double>&,
        DensityBranch, std::optional<double>, DensityRoot*
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
    bool ensure_tape(const std::vector<double>&, const std::vector<double>&, const detail::PredicateSignature&, const std::vector<IdealCorrelation>&, double, bool) const;
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
    mutable std::unique_ptr<detail::HelmholtzTape> previous_tape_;
    mutable SweepCounters counters_;
    mutable std::mutex mutex_;
};

DensityRoot solve_density(
    const HelmholtzModel& model,
    double temperature,
    double pressure,
    const std::vector<double>& mole_fractions,
    DensityBranch branch,
    std::optional<double> anchor_density,
    DensityRoot* sibling_root = nullptr  // unanchored only: the opposite extreme stable root
);

double reaction_reference_log_scale(
    const HelmholtzModel& model,
    const NeutralReferenceDeclaration& declaration,
    const std::vector<double>& stoichiometry
);

NeutralReferenceReactionDiagnostic reaction_reference(
    const HelmholtzModel& model,
    double temperature,
    double system_pressure,
    const NeutralReferenceDeclaration& declaration,
    const std::vector<double>& stoichiometry,
    const ReactionReference& reference
);

// Evaluate multiple reaction references while reusing one validated ladder per
// resolved source pressure within the call.  This remains a native owner; the
// Python surface exposes only the established singular diagnostic/result path.
std::vector<NeutralReferenceReactionDiagnostic> reaction_references(
    const HelmholtzModel& model,
    double temperature,
    double system_pressure,
    const NeutralReferenceDeclaration& declaration,
    const std::vector<std::vector<double>>& stoichiometries,
    const std::vector<ReactionReference>& references,
    bool pressure_jets = true,
    std::optional<double> solvent_anchor = std::nullopt,  // tried first for the pure-solvent point
    std::size_t temperature_order = 0,  // 1: dC/dT; 2: also d2C/dT2
    bool parameter_jets = false  // dC/dtheta for every active parameter
);

struct InfiniteDilutionReference final {
    std::array<double, 2> log_fugacity{};
    double molality = 0.0;
    double convergence_error = 0.0;
    std::optional<double> pure_solvent_log_fugacity;
    std::optional<State> pure_solvent_state;
    std::vector<std::array<std::optional<double>, 2>> parameter_derivative;
    std::vector<double> parameter_convergence_error;
};

struct ElectrolyteActivity final {
    std::array<double, 2> single_ion{};
    double mean_ionic_mole_fraction_scale = 0.0;
    double mean_ionic_log_molality_scale = 0.0;
    double mean_ionic_molality_scale = 0.0;
    std::optional<std::variant<double, core::Stop>> osmotic_coefficient;
    std::vector<double> mean_ionic_parameter_derivative, osmotic_parameter_derivative;
};

std::vector<double> formula_composition(
    const std::vector<double>& solvent_fractions,
    const Parameters& parameters,
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
