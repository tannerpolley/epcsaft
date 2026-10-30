#pragma once
#include "core/wider_scalar.hpp"
#include "core/export.hpp"
#include "core/trace.hpp"

#include "core/convergence.hpp"
#include "properties/helmholtz_tape.hpp"
#include "properties/thermochemistry.hpp"

#include <array>
#include <cstddef>
#include <deque>
#include <mutex>
#include <limits>
#include <optional>
#include <string>
#include <stdexcept>
#include <tuple>
#include <variant>
#include <vector>

namespace epcsaft {

using ParameterFamily = detail::ParameterFamily;
using ActiveParameter = detail::ActiveParameter;
void validate_cation_mixing(const Parameters& parameters);

// PressureSlope returns the reduced gradient, a_VV and the association floor, not a Hessian.
// Phase Hessians differentiate only V, positive amounts and (unless isothermal) T; other
// gradient/Hessian entries and the floor's coordinates are uncomputed padding, not physical zeros.
// Dry cation-mixing transverse entries are defined zero storage, not physical insertion derivatives (full EOS 12).
// StateHessian keeps the full gradient and T, V and active-parameter Hessian columns.
enum class DerivativeOrder { Value, Gradient, PressureSlope, Hessian, PhaseHessian, IsothermalPhaseHessian, StateHessian };

enum class NeutralProperty { ResidualHelmholtz, Pressure };

struct NeutralPropertyDerivatives final {
    double value = 0.0;
    std::vector<double> gradient;
    std::vector<double> hessian;
};

template <typename Base>
struct HelmholtzDerivativeValues final {
    Base value = 0.0;
    std::vector<Base> gradient;
    std::vector<Base> hessian;
    Base volume_second_derivative = 0.0;
    // Per-order association residual floor; zero when no association sites are present.
    Base gradient_floor = 0.0;
    std::vector<Base> gradient_error;
    std::vector<Base> site_fractions;
    std::optional<core::Convergence> first_numerical_failure;
};

using HelmholtzDerivatives = HelmholtzDerivativeValues<double>;

// Fixed-temperature pure-fluid projection of the shared volume series (EqID pure_isothermal_projection).
template <typename Base>
struct PureIsothermalValues final {
    HelmholtzDerivativeValues<Base> derivative;
    Base pressure = 0.0, chemical_potential_over_rt = 0.0, pressure_density_derivative = 0.0;
};

using PureIsothermal = PureIsothermalValues<double>;

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
    IdealEnthalpy,
    IdealInternalEnergy,
    IdealEntropy,
    IdealIsobaricHeatCapacity,
    IdealIsochoricHeatCapacity,
    ResidualEnthalpy,
    ResidualInternalEnergy,
    ResidualEntropy,
    ResidualIsobaricHeatCapacity,
    ResidualIsochoricHeatCapacity,
    TotalEnthalpy,
    TotalInternalEnergy,
    TotalEntropy,
    TotalIsobaricHeatCapacity,
    TotalIsochoricHeatCapacity,
    MolarVolumeTemperatureDerivative,
    MolarVolumePressureDerivative,
    LogFugacityCoefficient  // pressure chart only; one result per component, tagged by `component`
};
// The first CALORIC_PROPERTY_COUNT members are three families (ideal 0, residual 1, total 2) of
// five kinds (enthalpy 0, internal energy 1, entropy 2, Cp 3, Cv 4), in that order.
inline constexpr std::size_t CALORIC_PROPERTY_COUNT = 15;
inline std::size_t property_family(PropertyObservable o) noexcept { return static_cast<std::size_t>(o) / 5; }
inline std::size_t property_kind(PropertyObservable o) noexcept { return static_cast<std::size_t>(o) % 5; }
inline PropertyObservable property_observable(std::size_t family, std::size_t kind) noexcept {
    return static_cast<PropertyObservable>(5 * family + kind);
}
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
    std::optional<std::size_t> component;
    std::optional<double> value, action;
    PropertyActionStatus status = PropertyActionStatus::Nonfinite;
    std::size_t energy_order = 0;
};

enum class DensityBranch { Unspecified, Vapor, Liquid };
struct EPC_API DensityStructureError final : std::runtime_error {
    using std::runtime_error::runtime_error;
    DensityStructureError(const std::string& message, const std::optional<core::Convergence>& first)
        : std::runtime_error(message) { core::retain_first(first_numerical_failure, first); }
    std::optional<core::Convergence> first_numerical_failure;
};
inline void retain_density_failure(std::optional<core::Convergence>& first, const std::exception& error, const char* operation) {
    if (const auto* structure = dynamic_cast<const DensityStructureError*>(&error))
        core::retain_first(first, structure->first_numerical_failure);
    core::retain_failure(first, error, operation);
}
struct PropertyActionDiagnostics final {
    std::size_t site_rank = 0, site_factorizations = 0, site_solves = 0;
    double site_pivot_ratio = 0.0, association_residual = 0.0;
    double pressure_density_derivative = 0.0, pressure_residual = 0.0;
    double stability_determinant = 0.0;
    std::size_t pressure_factorizations = 0, pressure_solves = 0;
    std::size_t thermal_evaluations = 0, recordings = 0, directional_paths = 0, forward_sweeps = 0;
    DensityBranch branch = DensityBranch::Unspecified;
    std::optional<std::size_t> stable_root_count;
    std::optional<core::Convergence> first_numerical_failure;
};
struct PropertyActionBatch final {
    policy::SolverOptions solver_options;
    std::vector<PropertyActionResult> results;
    std::size_t maximum_energy_order = 0;
    PropertyActionDiagnostics diagnostics;
    // #86 pressure chart: root, first-order density derivatives per direction, status, and the root search's failure message and record.
    std::optional<double> molar_density;
    std::vector<double> density_action;
    std::optional<double> density_mixed_action;  // d^k rho over all k >= 2 directions, when observables are requested
    PropertyActionStatus density_status = PropertyActionStatus::CentralPropertyUnavailable;
    std::string density_message;
    bool density_structure_refused = false;
    std::optional<core::Convergence> density_convergence;
};
struct ResidualEnergyAction final {
    std::optional<double> action; PropertyActionStatus status = PropertyActionStatus::Nonfinite; PropertyActionDiagnostics diagnostics;
};
// Ordinary Taylor coefficients along a fixed-composition pressure-closed path.
struct ReferenceTaylor final {
    std::vector<double> contraction, volume;
    PropertyActionStatus status = PropertyActionStatus::Nonfinite;
    PropertyActionDiagnostics diagnostics;
};
struct ReferenceTaylorError final {
    std::array<double, 2> fixed{}, pressure{}, pressure_floor{};
    double temperature_coefficient = 0.0;
};
inline double require_component_property(const std::optional<double>& value) {
    if (!value) throw std::invalid_argument("ionic_insertion: absent charged component property");
    return *value;
}
struct State final {
    policy::SolverOptions solver_options;
    double hard_chain = 0.0;
    double dispersion = 0.0;
    double association = 0.0;
    double debye_huckel = 0.0;
    double born = 0.0;
    double dipolar = 0.0;
    double cation_mixing = 0.0;
    double bulk_relative_permittivity = 0.0;
    double residual_helmholtz = 0.0;
    double compressibility_factor = 0.0;
    double pressure = 0.0;
    std::vector<std::optional<double>> residual_chemical_potential_over_rt;
    std::optional<std::vector<std::optional<double>>> log_fugacity_coefficient;
    std::vector<PropertyActionStatus> component_property_status;
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
    std::vector<std::optional<double>> log_fugacity_parameter_derivative;  // component x active, fixed T/P/x
    std::vector<double> density_parameter_derivative, permittivity_parameter_derivative;
    std::vector<double> site_fractions;
    std::optional<core::Convergence> density_convergence;
    DensityBranch density_branch = DensityBranch::Unspecified;
    DensityBranch requested_density_branch = DensityBranch::Unspecified;
    EPC_API DensityBranch returned_density_branch() const;
    EPC_API bool density_branch_mismatch() const;
    std::optional<std::size_t> stable_root_count;
    std::optional<double> negative_slope_density;  // mol/m3, sampled cold mechanical-loop witness
    // Set only by state(..., derived_properties = true); fixed composition, SI molar units.
    std::optional<double> isothermal_compressibility, isobaric_expansivity, heat_capacity_ratio,
        isentropic_compressibility, speed_of_sound, joule_thomson_coefficient;
    std::optional<std::vector<std::optional<double>>> partial_molar_volume, partial_molar_enthalpy,
        partial_molar_entropy;  // component order
    std::optional<core::Convergence> first_numerical_failure;
};

EPC_API std::optional<std::array<double, 2>> superancillary_densities(const class HelmholtzModel&, double temperature);

struct DensityRoot final {
    double molar_density = 0.0;
    double pressure_residual = 0.0;
    double pressure_density_derivative = 0.0;
    std::size_t iterations = 0;
    std::optional<std::size_t> stable_root_count;
    DensityBranch branch = DensityBranch::Unspecified;
    core::Convergence convergence;
    std::optional<double> negative_slope_density;  // no new EOS evaluation
    std::optional<core::Convergence> first_numerical_failure;
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
    policy::SolverOptions solver_options;
    std::vector<DensityDiagnosticPoint> samples, pressure_roots, spinodals;
    std::vector<DensityDiagnosticBracket> declared_brackets;
    std::vector<DensityDiagnosticPoint> unresolved;
    std::optional<DensityDiagnosticPoint> finite_domain;
    DerivativeOrder evaluation_order = DerivativeOrder::PressureSlope;
    double geometric_cap = 0.0;
    std::size_t declared_crossings = 0, evaluations = 0;
    std::string stop;
    std::optional<core::Convergence> first_numerical_failure;
};

// A Stage-A neutral reference is a pure, neutral solvent plus one explicit
// full-width solute molality-ratio vector.  Mixed-solvent declarations are
// intentionally outside this first native consumer contract.
enum class SpeciesReference { InfiniteDilution, PureLiquid, PureGas };

struct NeutralReferenceDeclaration final {
    std::optional<std::size_t> solvent_index;
    std::vector<double> solute_molality_ratios;
    double standard_molality_mol_per_kg = 0.0;
    std::vector<SpeciesReference> species_reference;  // empty retains the existing solvent reference
};
inline bool operator==(const NeutralReferenceDeclaration& a, const NeutralReferenceDeclaration& b) {
    return std::tie(a.solvent_index, a.solute_molality_ratios, a.standard_molality_mol_per_kg, a.species_reference)
        == std::tie(b.solvent_index, b.solute_molality_ratios, b.standard_molality_mol_per_kg, b.species_reference);
}

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
};

struct NeutralReferenceTrace;
struct NeutralReferenceReactionDiagnostic final {
    std::shared_ptr<const NeutralReferenceTrace> trace;
    std::array<double, 2> numerical_error{NAN, NAN};
    std::string numerical_error_cause;
    policy::SolverOptions solver_options;
    std::size_t reaction_index = 0;
    std::optional<double> terminal_contraction;
    PropertyActionDiagnostics diagnostics;
    std::optional<core::Convergence> density_convergence;
    std::optional<double> observed_terminal_change;
    // Pressure derivatives of the reference contraction, admitted like dC/dT at scales P and P^2;
    // explicit fixed-pressure sources carry exact zero system-pressure derivatives.
    std::optional<double> terminal_contraction_pressure_derivative;
    std::optional<double> terminal_contraction_pressure_second_derivative;
    // Exact solvent-limit partials, supplied only when requested and numerically regular.
    std::optional<double> terminal_contraction_temperature_derivative;
    std::optional<double> observed_temperature_derivative_change;
    // d2C/dT2 at the reference pressure, set only for temperature order two.
    std::optional<double> terminal_contraction_temperature_second_derivative;
    std::optional<double> observed_temperature_second_derivative_change;
    // d2C/dTdP, only with temperature and pressure jets.
    std::optional<double> terminal_contraction_temperature_pressure_derivative;
    // dC/dtheta_a at the reference pressure, nullopt for an unavailable requested column.
    std::vector<std::optional<double>> terminal_contraction_parameter_derivative;
    std::vector<std::optional<double>> terminal_contraction_temperature_parameter_derivative;  // d2C/dTdtheta
    double reference_pressure_pa = 0.0;
    double solvent_density = 0.0;  // liquid density of the trace's pure-solvent point
    NeutralReferenceStatus status = NeutralReferenceStatus::Nonfinite;
    std::string message;
    struct Derivative final {
        std::size_t temperature_order = 0, pressure_order = 0;
        std::optional<std::size_t> parameter;
        PropertyActionDiagnostics diagnostics;
        std::optional<double> value;
        std::string message;
        std::optional<PropertyActionStatus> action_status;
    };
    std::vector<Derivative> derivatives;  // raw SI partials; system-P zeros for fixed sources
    std::vector<Derivative> lower_derivatives;
    // d^(a+b+[p]) C / dT^a dP^b [dtheta_p], nullopt where it was not admitted or not requested.
    std::optional<double> derivative_value(std::size_t a, std::size_t b, std::optional<std::size_t> p) const {
        const auto parameter = [&p](const std::vector<std::optional<double>>& v) {
            return *p < v.size() ? v[*p] : std::nullopt; };
        if (p && a + b == 0) return parameter(terminal_contraction_parameter_derivative);
        if (p && a == 1 && b == 0) return parameter(terminal_contraction_temperature_parameter_derivative);
        if (!p && a + b == 0) return terminal_contraction;
        if (!p && a == 1 && b == 0) return terminal_contraction_temperature_derivative;
        if (!p && a == 2 && b == 0) return terminal_contraction_temperature_second_derivative;
        if (!p && a == 0 && b == 1) return terminal_contraction_pressure_derivative;
        if (!p && a == 0 && b == 2) return terminal_contraction_pressure_second_derivative;
        if (!p && a == 1 && b == 1) return terminal_contraction_temperature_pressure_derivative;
        for (const auto& column : derivatives)
            if (column.temperature_order == a && column.pressure_order == b && column.parameter == p) return column.value;
        return std::nullopt;
    }
    std::string derivative_cause(std::size_t a, std::size_t b, std::optional<std::size_t> p) const {
        for (const auto* columns : {&lower_derivatives, &derivatives})
            for (const auto& column : *columns)
                if (column.temperature_order == a && column.pressure_order == b && column.parameter == p && !column.value)
                    return column.message;
        return message;
    }
};

class HelmholtzModel final {
public:
    EPC_API explicit HelmholtzModel(
        Parameters parameters,
        std::vector<ActiveParameter> active = {},
        std::optional<ThermochemistryRecord> thermochemistry = std::nullopt,
        policy::SolverOptions solver_options = {}
    );

    EPC_API ~HelmholtzModel();

    EPC_API HelmholtzDerivatives evaluate(
        double temperature,
        double volume,
        const std::vector<double>& amounts,
        DerivativeOrder order
    ) const;
    // Mixed contraction of the extensive residual Helmholtz scalar with three to five
    // directions of width 2 + C + A in [T, V, n, active parameters]: E3, E4 or E5.
    // Errors of E's ordinary Taylor coefficients on a requested path, including the site recurrence.
    EPC_API std::vector<double> residual_coefficient_errors(const std::vector<std::vector<double>>& regular,
        const std::vector<double>& sites) const;
    EPC_API ResidualEnergyAction residual_energy_action(double, double, const std::vector<double>&, const std::vector<std::vector<double>>&) const;
    EPC_API ReferenceTaylor reference_taylor(double temperature, double volume,
        const std::vector<double>& composition, const std::vector<double>& stoichiometry,
        const PressureActionDirection&, std::size_t, bool infinite_dilution = false, ReferenceTaylorError* error = nullptr) const;
    EPC_API NeutralPropertyDerivatives neutral_property(
        double temperature,
        double molar_density,
        NeutralProperty property,
        DerivativeOrder order
    ) const;
    EPC_API State state(
        double temperature,
        double molar_density,
        const std::vector<double>& mole_fractions,
        bool derived_properties = false
    ) const;
    EPC_API State state_at_pressure(
        double temperature,
        double pressure,
        const std::vector<double>& mole_fractions,
        DensityBranch branch,
        std::optional<double> anchor_density = std::nullopt,
        bool derived_properties = false
    ) const;
    EPC_API PropertyActionBatch density_actions(
        double, double, const std::vector<double>&,
        const std::vector<PropertyObservable>&,
        const std::vector<DensityActionDirection>&
    ) const;
    EPC_API PropertyActionBatch pressure_actions(
        double, double, const std::vector<double>&, DensityBranch,
        const std::vector<PropertyObservable>&,
        const std::vector<PressureActionDirection>&,
        std::optional<double> = std::nullopt
    ) const;
    EPC_API double packing_fraction_coefficient(
        double temperature,
        const std::vector<double>& mole_fractions
    ) const;
    EPC_API double molar_density_upper_bound(
        double temperature,
        const std::vector<double>& mole_fractions
    ) const;
    EPC_API DensityDiagnostics density_diagnostics(
        double temperature, double pressure, const std::vector<double>& composition,
        const std::vector<double>& nodes, double density_tolerance, double log_tolerance,
        std::size_t max_evaluations, double max_seconds, bool adaptive,
        std::optional<double> finite_anchor, double finite_upper, DerivativeOrder evaluation_order
    ) const;
    // Physical-state admission (full-eos.md §12): |sum z_i x_i| <= COMPOSITION_SUM_TOLERANCE (absolute).
    // Public state requests and equilibrium feeds only; internal iterates use the unconstrained extension.
    EPC_API void validate_electroneutrality(const std::vector<double>& mole_fractions) const;
    EPC_API SweepCounters counters() const;
    // reaction_references results of this immutable model by exact inputs: a short FIFO (#380).
    // Inputs [first, last) only add columns: a stored result covers smaller ones there (#595).
    std::optional<std::vector<NeutralReferenceReactionDiagnostic>> stored_references(const std::vector<double>& inputs,
        std::size_t first, std::size_t last) const;
    void store_references(std::vector<double> inputs, const std::vector<NeutralReferenceReactionDiagnostic>& result) const;
    EPC_API void reset_counters() const;
    const policy::SolverOptions& solver_options() const { return solver_options_; }
    template <typename Base>
    PureIsothermalValues<Base> pure_isothermal(Base temperature, Base volume, Base amount = 1.0) const;
    HelmholtzDerivativeValues<WiderScalar> isothermal(WiderScalar temperature, WiderScalar volume,
        const std::vector<WiderScalar>& amounts) const;
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
        double temperature, const std::vector<double>& composition, bool ideal = false
    ) const;
    PropertyActionBatch directional_actions(double, double, const std::vector<double>&, const std::vector<PropertyObservable>&, const std::vector<DensityActionDirection>&, const DensityRoot* = nullptr, const std::vector<double>& = {}, double = 0.0) const;
    bool ensure_tape(const std::vector<double>&, const std::vector<double>&, const detail::PredicateSignature&, const std::vector<IdealCorrelation>&, double, bool) const;
    HelmholtzDerivatives evaluate_locked(
        double temperature,
        double volume,
        const std::vector<double>& amounts,
        DerivativeOrder order
    ) const;

    const policy::SolverOptions solver_options_;
    Parameters parameters_;
    std::vector<ActiveParameter> active_;
    std::optional<ThermochemistryRecord> thermochemistry_;
    mutable std::unique_ptr<detail::HelmholtzTape> tape_;
    mutable std::unique_ptr<detail::HelmholtzTape> previous_tape_;
    mutable std::unique_ptr<detail::HelmholtzTape> reference_tape_;
    mutable SweepCounters counters_;
    mutable std::deque<std::pair<std::vector<double>, std::vector<NeutralReferenceReactionDiagnostic>>> references_;
    // reference_taylor central Hessians of this immutable model by exact (T, V, n): a short FIFO (#449).
    mutable std::deque<std::pair<std::vector<double>, HelmholtzDerivatives>> centrals_;
    mutable std::mutex mutex_;
};
#if defined(_WIN32) && !defined(EPCSAFT_NATIVE_BUILD) && !defined(EPCSAFT_NATIVE_STATIC)
extern template EPC_API PureIsothermalValues<double> HelmholtzModel::pure_isothermal(double, double, double) const;
#endif


EPC_API DensityRoot solve_density(
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

EPC_API NeutralReferenceReactionDiagnostic reaction_reference(
    const HelmholtzModel& model,
    double temperature,
    double system_pressure,
    const NeutralReferenceDeclaration& declaration,
    const std::vector<double>& stoichiometry,
    const ReactionReference& reference
);

// Raw source-pressure partial; infinite_dilution selects the ion-free first-insertion limit.
void reference_errors(const HelmholtzModel&, double, const std::vector<double>&, NeutralReferenceReactionDiagnostic&);
EPC_API ResidualEnergyAction reference_partial_from_state(const HelmholtzModel&, double temperature,
    const std::vector<double>& composition, const std::vector<double>& stoichiometry,
    std::size_t temperature_order, std::size_t pressure_order,
    std::optional<std::size_t> parameter, const State& state, bool infinite_dilution = false);

struct SpeciesActivities final {
    std::optional<double> standard_molality_mol_per_kg;
    std::optional<std::size_t> solvent_index;
    std::vector<std::optional<double>> log_activities, activities;
    std::vector<std::string> causes;
};
// EqID reaction_activity_quotient; a null state selects ideal mole fractions.
EPC_API SpeciesActivities species_activities(const std::vector<double>&, const State*,
    const std::vector<NeutralReferenceReactionDiagnostic::Derivative>& = {},
    const Parameters* = nullptr, const NeutralReferenceDeclaration* = nullptr,
    ReferenceSourceBasis = ReferenceSourceBasis::RawMoleFractionInfiniteDilution);

// Evaluate multiple reaction references while reusing one solvent closure per
// resolved source pressure within the call.  This remains a native owner; the
// Python surface exposes only the established singular diagnostic/result path.
EPC_API std::vector<NeutralReferenceReactionDiagnostic> reaction_references(
    const HelmholtzModel& model,
    double temperature,
    double system_pressure,
    const NeutralReferenceDeclaration& declaration,
    const std::vector<std::vector<double>>& stoichiometries,
    const std::vector<ReactionReference>& references,
    bool pressure_jets = true,
    std::optional<double> solvent_anchor = std::nullopt,  // tried first for the pure-solvent point
    std::size_t temperature_order = 0,  // 1: dC/dT; 2: also d2C/dT2
    bool parameter_jets = false,  // dC/dtheta for every active parameter
    std::size_t reference_order = 0,  // consumed C3/C4 and one-theta columns; 0 preserves lower work
    std::vector<std::vector<NeutralReferenceReactionDiagnostic::Derivative>>* species_values = nullptr
);

// Fugacity-based Henry constant H_i^px = lim f_i/x_i = RT rho0 exp(L_i) of a neutral solute in one
// pure solvent of an all-neutral model (EqID: henry_constant_px), from one unit reaction_references row.
struct HenryConstant final {
    std::size_t component_index = 0;
    std::optional<double> raw_log_reference;  // admitted limit of ln(rho/rho0) + mu_i^res/RT, rho0 = 1 mol/m3
    std::optional<double> henry_constant_px;  // Pa
    std::optional<double> observed_terminal_change;
    NeutralReferenceStatus status = NeutralReferenceStatus::Nonfinite;
    std::string message;
};
EPC_API std::vector<HenryConstant> henry_constants(const HelmholtzModel& model, double temperature, double pressure,
    const std::vector<double>& solvent_fractions, const std::vector<std::size_t>& solute_indices);

struct EPC_EXCEPTION IonicReferenceBackground final : std::invalid_argument {
    using std::invalid_argument::invalid_argument;
};

struct InfiniteDilutionReference final {
    policy::SolverOptions solver_options;
    const HelmholtzModel* model = nullptr;
    double temperature = 0.0, pressure = 0.0;
    std::vector<double> solvent_fractions;
    PropertyActionStatus solvent_status = PropertyActionStatus::CentralPropertyUnavailable;
    std::string solvent_cause;
    std::optional<core::Convergence> solvent_convergence;
    std::optional<std::array<std::size_t, 2>> ion_indices;
    std::array<double, 2> log_fugacity{};  // analytic ion ln(phi) limit
    double molality = 0.0;  // the zero-solute limit, mol/kg
    std::optional<double> convergence_error;  // no terminal sequence is measured
    std::optional<double> pure_solvent_log_fugacity;
    std::optional<State> pure_solvent_state, solvent_state;
    std::array<PropertyActionDiagnostics, 2> diagnostics;
    std::vector<std::array<ResidualEnergyAction, 2>> parameter_derivative;
};

struct SolventActivity final {
    std::optional<double> solvent_activity;
    std::optional<std::variant<double, core::Stop>> osmotic_coefficient;
    std::vector<std::optional<double>> solvent_activity_parameter_derivative, osmotic_parameter_derivative;
    std::vector<PropertyActionStatus> parameter_status;
    std::vector<std::string> parameter_cause;
    PropertyActionStatus status = PropertyActionStatus::CentralPropertyUnavailable;
    std::string cause, osmotic_cause;
    State state;
};
EPC_API InfiniteDilutionReference solvent_reference(const HelmholtzModel&, double, double, const std::vector<double>&,
    std::optional<double> density = std::nullopt);
EPC_API SolventActivity solvent_activity(const HelmholtzModel&, double, double, const std::vector<double>&,
    const InfiniteDilutionReference&, std::optional<State> state = std::nullopt,
    std::optional<double> formula_solute_ratio = std::nullopt);

struct ElectrolyteActivity final {
    policy::SolverOptions solver_options;
    std::array<double, 2> single_ion{};
    double mean_ionic_mole_fraction_scale = 0.0;
    double mean_ionic_log_molality_scale = 0.0;
    double mean_ionic_molality_scale = 0.0;
    std::optional<std::variant<double, core::Stop>> osmotic_coefficient;
    std::vector<double> osmotic_parameter_derivative;
    std::vector<ResidualEnergyAction> mean_ionic_parameter_derivative;
    State state;
};

struct ElectrolyteInputDirection final {
    double temperature = 0.0, pressure = 0.0, molality = 0.0;
};
struct ElectrolyteInputAction final {
    policy::SolverOptions solver_options;
    double mean_ionic = NAN, osmotic = NAN;
    std::array<double, 2> reference_column{};  // analytic limit column
    std::vector<PropertyActionDiagnostics> diagnostics;
    PropertyActionStatus status = PropertyActionStatus::Nonfinite;
    std::string message;
};
ElectrolyteInputAction electrolyte_input_action(const HelmholtzModel&, double temperature, double pressure,
    const std::vector<double>& solvent, const std::array<std::size_t, 2>& ions, const std::array<int, 2>& stoichiometry,
    double molality, const InfiniteDilutionReference&, const ElectrolyteActivity&,
    const ElectrolyteInputDirection&);

EPC_API std::vector<double> formula_composition(
    const std::vector<double>& solvent_fractions,
    const Parameters& parameters,
    const std::array<std::size_t, 2>& ion_indices,
    const std::array<int, 2>& stoichiometry,
    double molality, const std::optional<std::vector<double>>& finite_solute_molalities = std::nullopt,
    double* declared_solute_ratio = nullptr
);

EPC_API InfiniteDilutionReference infinite_dilution_reference(
    const HelmholtzModel& model,
    double temperature,
    double pressure,
    const std::vector<double>& solvent_fractions,
    const std::array<std::size_t, 2>& ion_indices,
    const std::array<int, 2>& stoichiometry,
    std::optional<double> solvent_density = std::nullopt
);

EPC_API ElectrolyteActivity electrolyte_activity(
    const HelmholtzModel& model,
    double temperature,
    double pressure,
    const std::vector<double>& solvent_fractions,
    const std::array<std::size_t, 2>& ion_indices,
    const std::array<int, 2>& stoichiometry,
    double molality,
    const InfiniteDilutionReference& reference,
    std::optional<double> finite_density = std::nullopt,
    const std::optional<std::vector<double>>& finite_solute_molalities = std::nullopt
);

namespace properties {
// Unit-amount residual Helmholtz coefficients in V+t, used by density pressure closure.
EPC_API Taylor<3> density_volume_taylor(const HelmholtzModel& model, double temperature, double density,
    const std::vector<double>& composition, double& gradient_floor);

// The factorized mass-action Jacobian of the solved site fractions; a singular or ill-conditioned one
// throws, or with `status` is reported there.
Eigen::FullPivLU<Eigen::MatrixXd> site_factorization(const std::vector<double>& sites, const std::vector<double>& composition,
    const Parameters& parameters, double number_density, const std::vector<double>& delta, PropertyActionStatus* status = nullptr);
void site_factorization(Eigen::FullPivLU<Eigen::MatrixXd>& factorization, const std::vector<double>& sites,
    const std::vector<double>& composition, const Parameters& parameters, double number_density,
    const std::vector<double>& delta, PropertyActionStatus* status = nullptr);
bool regular_density(double density, double slope, const std::optional<core::Convergence>& record);
// ALG-DENSITY single loop: an Unspecified root below the stored negative-slope sample is vapor-side, above it
// liquid-side; explicit labels and roots without a sample are kept. vapor_side refuses a liquid reference (#423).
DensityBranch density_side(DensityBranch branch, double density, std::optional<double> sample);
std::optional<std::string> vapor_side(const State& state);
// Anchored closure, regularity-checked; a failed anchor is traced and retried once cold on the same branch.
State close_density(const HelmholtzModel& model, double temperature, double pressure,
                    const std::vector<double>& composition, DensityBranch branch, std::optional<double> prior_density);
}  // namespace properties

}  // namespace epcsaft
