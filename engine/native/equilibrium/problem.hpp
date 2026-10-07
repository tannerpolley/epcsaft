#pragma once

#include "core/convergence.hpp"
#include "properties/helmholtz_model.hpp"

#include <cstddef>
#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>
#include <map>
#include <memory>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace epcsaft::core { class Scope; }
namespace epcsaft::equilibrium {

enum class PhaseKind { Liquid, Vapor };
enum class AmountKind { Free, Pinned };
enum class FeedKind { MoleFractions, Amounts, MassFractions, Molalities };

inline const char* phase_kind_name(PhaseKind kind) noexcept {
    return kind == PhaseKind::Liquid ? "liquid" : "vapor";
}

inline const char* amount_kind_name(AmountKind kind) noexcept {
    return kind == AmountKind::Free ? "free" : "pinned";
}

inline const char* feed_kind_name(FeedKind kind) noexcept {
    switch (kind) {
    case FeedKind::MoleFractions: return "mole_fractions";
    case FeedKind::Amounts: return "amounts";
    case FeedKind::MassFractions: return "mass_fractions";
    case FeedKind::Molalities: return "molalities";
    }
    return "unknown";
}

struct AmountSpec final {
    AmountKind kind = AmountKind::Free;
    double fraction = 0.0;
    static AmountSpec free() { return {}; }
    static AmountSpec pinned(double value) { return {AmountKind::Pinned, value}; }
};

struct FreeValue final {
    bool free = false;
    double value = 0.0;
    bool absent = false;
    static FreeValue fixed(double value) { return {false, value, false}; }
    static FreeValue variable() { return {true, 0.0, true}; }
    static FreeValue variable(double guess) { return {true, guess, false}; }

    double seed_or(double fallback) const noexcept {
        return free && absent ? fallback : value;
    }

    bool is_positive() const noexcept {
        return !absent && std::isfinite(value) && value > 0.0;
    }
};

struct PhaseSpec final {
    std::string name;
    PhaseKind kind = PhaseKind::Liquid;
    std::vector<std::size_t> support;
    // An omitted support is completed; an explicitly declared support never gains an index,
    // and a declared empty support is rejected.
    bool support_declared = false;
    AmountSpec amount;
    std::vector<double> composition_guess;
    std::optional<double> packing_guess;

    bool is_pinned() const noexcept { return amount.kind == AmountKind::Pinned; }
    bool is_free() const noexcept { return amount.kind == AmountKind::Free; }
};

struct FeedSpec final {
    FeedKind kind = FeedKind::Amounts;
    std::vector<double> values;
    std::vector<std::size_t> solvent;
    double solvent_mass_kg = 1.0;
};

// The only admitted correlation identity: activity (rho x / (1 mol/m^3)) exp(mu_res / RT).
inline constexpr const char* EOS_STANDARD_STATE_ID = "eos-helmholtz-coordinate-basis";

// Exactly one target: scalar K, affinity in J/mol, or the P3 ln K(T) correlation on the
// EOS density basis. No target silently overrides another.
struct ReactionSpec final {
    std::string name;
    std::vector<double> stoichiometry;
    std::optional<double> equilibrium_constant;
    std::optional<double> affinity_j_per_mol;
    std::optional<ReactionLogPolynomial> correlation;
    std::optional<ReactionReference> reference;
};

// Exact declaration equality, shared by fitting admission/grouping and continuation.
inline bool same_reaction(const ReactionSpec& x, const ReactionSpec& y) {
    const auto correlation = [](const auto& v) {
        return std::tie(v.a, v.b, v.c, v.d, v.reference_temperature, v.temperature_min,
                        v.temperature_max, v.standard_state_id);
    };
    const auto reference = [](const auto& v) { return std::tie(v.reference_pressure_pa, v.source_basis); };
    const auto optional_equal = [](const auto& a, const auto& b, const auto& fields) {
        return a.has_value() == b.has_value() && (!a || fields(*a) == fields(*b));
    };
    return std::tie(x.name, x.stoichiometry, x.equilibrium_constant, x.affinity_j_per_mol)
        == std::tie(y.name, y.stoichiometry, y.equilibrium_constant, y.affinity_j_per_mol)
        && optional_equal(x.correlation, y.correlation, correlation)
        && optional_equal(x.reference, y.reference, reference);
}

// Explicit finite TPD search options; unset supports select the ALG-STABILITY defaults.
struct StabilityOptions final {
    bool evaluate_supplied_trials = false;
    std::optional<std::vector<PhaseSpec>> trials;
    std::optional<std::size_t> declaration_budget;
    std::optional<std::size_t> trial_budget;
    std::vector<double> standard_uncertainty_over_rt;
};
struct DiscoveryRecord;

struct EquilibriumProblem final {
    std::vector<PhaseSpec> phases;
    FreeValue temperature;
    FreeValue pressure;
    FeedSpec feed;
    std::vector<ReactionSpec> reactions;
    std::vector<double> standard_gibbs_j_per_mol;
    std::optional<NeutralReferenceDeclaration> source_reference;
    bool topology_declared = true;
    std::optional<StabilityOptions> discovery;
    std::optional<double> enthalpy_j;
    std::optional<double> entropy_j_per_k;
    // Caller row criterion for the solve; the named policy value is the only default.
    double row_tolerance = policy::EQUILIBRIUM_ROW_RESIDUAL.value;

    std::size_t component_count() const noexcept {
        return feed.values.size();
    }
};

enum class RowKind { Chemical, Pressure, Moiety, Electro, Normalization, Enthalpy, Entropy };

inline const char* row_kind_name(RowKind kind) noexcept {
    switch (kind) {
    case RowKind::Chemical: return "chemical_potential_stationarity";
    case RowKind::Pressure: return "pressure_equality";
    case RowKind::Moiety: return "material_balance";
    case RowKind::Electro: return "phase_electroneutrality";
    case RowKind::Normalization: return "phase_normalization";
    case RowKind::Enthalpy: return "extensive_enthalpy_specification";
    case RowKind::Entropy: return "extensive_entropy_specification";
    }
    return "unknown";
}

struct Row final {
    RowKind kind = RowKind::Chemical;
    std::size_t phase = 0;
    std::size_t index = 0;
    const char* eqid = "chemical_potential_stationarity";
};

struct PointEvaluation final {
    std::vector<double> values;
    std::vector<double> jacobian;
    double objective = 0.0;
    double floor = 0.0;

    double residual_norm() const noexcept {
        double result = 0.0;
        for (const double value : values) {
            const double magnitude = std::abs(value);
            if (!std::isfinite(magnitude)) return std::numeric_limits<double>::infinity();
            result = std::max(result, magnitude);
        }
        return result;
    }
};

struct Layout final {
    std::vector<std::size_t> x;
    std::optional<std::size_t> amount;
    std::size_t volume = 0;
};

struct PackingTerms final {
    double coefficient = 0.0;
    double temperature_derivative = 0.0;
    std::vector<double> component;
};

// derivative[k][i] is d^k (g_i / RT) / dT^k at the prepared temperature.  The
// pressure fields are the matching first and second derivatives at fixed T.
struct StandardTerm final {
    std::vector<std::vector<double>> derivative;
    std::vector<double> pressure_derivative;
    std::vector<double> pressure_second_derivative;
    std::vector<double> temperature_pressure_derivative;  // d2(g_i/RT)/dTdP; empty unless requested
    // parameter_derivative[a][i] is d(g_i/RT)/dtheta_a at fixed T and P; empty unless requested.
    std::vector<std::vector<double>> parameter_derivative, temperature_parameter_derivative,
        pressure_parameter_derivative;
    std::vector<NeutralReferenceReactionDiagnostic> reference_diagnostics;
};

// One latest-point phase evaluation shared by rows, Jacobian, gradient and curvature.
struct PhasePoint final {
    std::vector<double> point;
    std::vector<Layout> shape;
    double temperature = 0.0;
    double pressure = 0.0;
    std::vector<std::vector<double>> composition;
    std::vector<double> amount;
    std::vector<double> volume;
    std::vector<double> eta;
    std::vector<PackingTerms> packing;
    std::vector<HelmholtzDerivatives> derivative;
    StandardTerm standard;
    // A referenced problem's standard term with jets along every active parameter, built on the
    // first parameter-driven row series at this point and reused by every later direction.
    mutable std::optional<StandardTerm> parameter_standard;
};

// Taylor series of one phase along a path: point[k] is the k-th coefficient of the phase as
// the model sees it, (T, v, x, active parameters); gradient[k] the matching coefficient of
// the extensive residual Helmholtz gradient; amount[k] and packing[k] the phase amount and
// packing fraction.
struct PhaseSeries final {
    std::vector<std::vector<double>> point, gradient;
    std::vector<double> amount, packing;
};

// Defining rows and physical state along y(t) = sum_k y_k t^k with controls
// q(t) = q0 + t dq, dq = [T, P, n_feed (C), active parameters (A)]; a free T or P entry of
// dq is unused. `energy` is the worst status of the model's residual contractions.
struct RowSeries final {
    std::vector<std::vector<double>> rows;
    std::vector<double> temperature, pressure;
    std::vector<PhaseSeries> phase;
    std::array<std::vector<std::vector<double>>, 3> reaction_shift; // L, L_T, L_TT Taylor coefficients
    PropertyActionStatus energy = PropertyActionStatus::Available;
    PropertyActionDiagnostics energy_diagnostics;
    std::size_t maximum_energy_order = 0;
};

void restrict_composition_guess(PhaseSpec& phase, const std::vector<std::size_t>& support);

struct CompiledSystem final {
    const HelmholtzModel* model = nullptr;
    EquilibriumProblem problem;
    std::vector<double> n_feed;
    std::vector<double> feed_composition;
    std::vector<double> balance;
    std::vector<double> standard_gibbs;
    std::vector<std::vector<double>> B;
    std::vector<std::vector<double>> inventory_restrictions;
    std::vector<double> inventory_upper_bounds;
    std::vector<std::vector<std::size_t>> declared_supports;
    std::vector<Row> rows;
    std::size_t pinned_count = 0;
    std::size_t free_intensive_count = 0;
    std::size_t rank = 0;
    std::size_t charge_row_count = 0;
    std::size_t primal_size = 0;
    long rows_minus_unknowns = 0;
    std::size_t generic_nullity = 0;
    bool feasible = false;
    std::vector<double> F_y;
    std::vector<std::size_t> removed_species;
    std::vector<std::string> row_names;
    double feed_total = 0.0;
    bool has_ions = false;
    bool has_reaction = false;
    bool conserves_total_moles = true;
    std::size_t jacobian_rank = 0;
    std::vector<NeutralReferenceReactionDiagnostic> reference_diagnostics;

    // State-only preparation differentiates V, support amounts and T when free (phase Hessians);
    // fixed controls and active parameters are padding. Sensitivity and lever callers need full.
    PhasePoint prepare(const std::vector<double>& point, bool state_only = false) const;
    PointEvaluation evaluate(const std::vector<double>& point) const;
    // Residuals that classify a solve: for a free-amount (Gibbs) cell the chemical and
    // pressure rows are scaled by N_p x_i / N_feed and N_p / N_feed, the first-order
    // conditions of the constrained Gibbs minimum in the solved chart; a pinned
    // (feasibility) cell classifies its rows unscaled. Balance, charge and normalization
    // rows are unchanged in both routes.
    std::vector<double> classified_residuals(const PhasePoint& prepared, const std::vector<double>& rows) const;
    double classified_residual_norm(const std::vector<double>& point) const;
    PointEvaluation evaluate_at(const PhasePoint& prepared) const;
    std::vector<double> generic_point() const;
    // First solved-chart coordinate of each declared phase: its log compositions, then a free amount and log packing.
    std::vector<std::size_t> phase_offsets() const;
    std::vector<double> jacobian(const std::vector<double>& point) const;
    std::vector<double> jacobian_at(const PhasePoint& prepared) const;
    std::vector<double> lever_jacobian(const std::vector<double>& point) const;
    std::vector<double> objective_gradient(const std::vector<double>& point) const;
    std::vector<double> objective_gradient_at(const PhasePoint& prepared) const;
    std::vector<double> hessian(
        const std::vector<double>& point,
        const std::vector<double>& constraint_multipliers = {}
    ) const;
    std::vector<double> hessian_at(
        const PhasePoint& prepared, const std::vector<double>& constraint_multipliers
    ) const;
    // Multiplier-weighted constraint curvature only; the objective curvature is separate.
    std::vector<double> constraint_hessian_at(
        const PhasePoint& prepared, const std::vector<double>& constraint_multipliers
    ) const;
    // Row and physical Taylor coefficients through `order` (at most 3) along the state
    // coefficients `state` (state[0] is the prepared point) and the control direction.
    RowSeries row_series(
        const PhasePoint& prepared, const std::vector<std::vector<double>>& state,
        const std::vector<double>& control, std::size_t order
    ) const;

    std::size_t unknown_count() const noexcept {
        return primal_size + rank + charge_row_count;
    }

    std::size_t constraint_count() const noexcept { return rows.size(); }

    bool has_caloric_condition() const noexcept {
        return problem.enthalpy_j.has_value() || problem.entropy_j_per_k.has_value();
    }

    bool feasibility_cell() const noexcept {
        return pinned_count != 0 || has_caloric_condition();
    }

    bool is_square() const noexcept { return rows_minus_unknowns == 0; }
};

CompiledSystem compile_problem(
    const HelmholtzModel& model, const EquilibriumProblem& problem
);

// Conserved inventory only: maximum is mol in amount mode, dimensionless in fraction mode.
// witness is n/N_feed, or (y, tau=N_feed/N); finite_witness uses the same chart.
struct InventoryBound final {
    std::size_t species = 0;
    std::string status, termination;
    std::optional<double> maximum, feed_normalized_amount, log_fraction;
    double balance_residual = NAN, normalization_residual = NAN, lower_bound_residual = NAN;
    std::vector<double> witness, amounts_mol, finite_witness;
    bool limiting = false, supported = false;
    double floor = 0.0, lp_error = std::numeric_limits<double>::quiet_NaN();
    double finite_residual = std::numeric_limits<double>::quiet_NaN();
    std::string feasibility_termination;
    std::vector<double> conserved_functional;
};
std::vector<InventoryBound> inventory_bounds(const CompiledSystem& compiled,
    const std::vector<std::size_t>& species, bool fraction, double floor);

// Candidate-topology outcome of the existing producer loop: the first vanished free phase
// in declaration order, else the first coincident pair in lexicographic ordinal order.
enum class TopologyEvent { None, PhaseDisappeared, PhasesCoalesced };

struct EquilibriumResult final {
    const HelmholtzModel* model_identity = nullptr;
    std::shared_ptr<const EquilibriumProblem> declaration;
    std::shared_ptr<DiscoveryRecord> discovery;
    bool success = false;
    std::string message;
    std::string solver_status;
    std::vector<double> temperature;
    std::vector<double> pressure;
    std::vector<double> amounts;
    std::vector<double> mole_fractions;
    std::vector<double> molar_densities;
    std::vector<double> packing_fractions;
    std::vector<double> multipliers;
    std::vector<double> constraint_multipliers;
    std::vector<double> residuals;
    // The residuals the classification actually used (see CompiledSystem::classified_residuals).
    std::vector<double> classified_residuals;
    std::vector<double> sensitivity;
    std::vector<double> coordinates;
    std::size_t phase_count = 0;
    std::size_t iterations = 0;
    std::size_t evaluations = 0;
    std::size_t hessian_calls = 0;
    double objective = 0.0;
    double reduced_hessian_minimum = 0.0;
    core::Convergence convergence;
    std::vector<Row> rows;
    std::size_t unknown_count = 0;
    std::size_t jacobian_rank = 0;
    std::size_t generic_nullity = 0;
    // Absent without a finite final row residual; otherwise residual <= caller row_tolerance.
    std::optional<bool> requested_tolerance_met;
    TopologyEvent topology_event = TopologyEvent::None;
    std::optional<std::size_t> topology_phase;
    std::optional<std::pair<std::size_t, std::size_t>> topology_pair;
    std::vector<NeutralReferenceReactionDiagnostic> reference_diagnostics;

    bool has_state() const noexcept {
        return success && !temperature.empty() && !pressure.empty()
            && !amounts.empty() && !mole_fractions.empty();
    }

    double residual_norm() const noexcept {
        double result = 0.0;
        for (const double value : residuals) {
            const double magnitude = std::abs(value);
            if (!std::isfinite(magnitude)) return std::numeric_limits<double>::infinity();
            result = std::max(result, magnitude);
        }
        return result;
    }
};

struct StabilityTrial final {
    PhaseSpec phase;
    std::string seed, evidence, message;
    EquilibriumResult optimizer;
    std::vector<double> raw_coordinates, raw_residuals, composition, split_coordinates;
    std::vector<std::vector<double>> basis, row_transformation;
    // Quantities carry units in their names; errors are local estimates, never enclosures.
    std::map<std::string, double> quantities;
    std::vector<std::map<std::string, double>> split_checks;
    std::optional<std::size_t> donor;
    bool resolved_descent = false;
};
struct StabilityReport final {
    bool completed = false;
    std::string message;
    std::vector<StabilityTrial> trials;
    std::vector<double> candidate_coordinates, standard_over_rt, plane_over_rt;
    std::vector<std::map<std::string, double>> present_phases;
    std::map<std::string, double> work;
};
struct DiscoveryRecord final {
    EquilibriumProblem declaration;
    EquilibriumResult candidate;
    std::vector<StabilityReport> rounds;
    std::vector<EquilibriumResult> declarations;
    std::vector<std::string> declaration_roles;
    std::vector<std::vector<double>> declaration_attempt_coordinates;
    std::map<std::string, double> work;
};
CompiledSystem compile_stability_trial(const HelmholtzModel&, double, double, const PhaseSpec&,
    const std::vector<double>&, const std::vector<double>&, std::vector<std::vector<double>>&);
double present_phase_residual_distance(double, double);
bool descent_resolved(double, double, double);
std::vector<double> reconstruct_trial_composition(std::vector<double>, const std::vector<std::size_t>&, const std::vector<int>&);
std::pair<EquilibriumProblem, std::vector<double>> remove_discovery_phase(const CompiledSystem&,
    const std::vector<double>&, std::map<std::string, double>* = nullptr);
StabilityReport stability_test(const CompiledSystem&, const EquilibriumResult&, const StabilityOptions& = {});
EquilibriumResult discover_phases(const HelmholtzModel&, const EquilibriumProblem&);

EquilibriumResult solve_equilibrium(
    const HelmholtzModel& model,
    const EquilibriumProblem& problem,
    const std::optional<EquilibriumResult>& start = std::nullopt
);
// The same solve of a compiled problem (its model must outlive it), inside the caller's open span if given.
EquilibriumResult solve_compiled(const CompiledSystem& compiled, const std::optional<EquilibriumResult>& start = std::nullopt,
                                 core::Scope* outer = nullptr);

// ---- Solved-state derivative actions (E5): native records, no public promise.

enum class SolvedStateObservableKind {
    Temperature,
    Pressure,
    PhaseAmount,
    PhaseComponentAmount,
    PhaseMoleFraction,
    PhaseVolume,
    PhaseMolarDensity,
    PhasePackingFraction,
    PhasePressure,
    PhaseResidualChemicalPotentialOverRt,
    PhaseProperty,
    // ln a_i = ln(x_i rho / (1 mol/m^3)) + mu_res_i/RT on EOS_STANDARD_STATE_ID; it excludes any
    // standard-state term, and f_i = a_i R T (1 mol/m^3) in Pa.
    PhaseLogActivity,
    PhaseGibbsMolarVolume,
    SystemEnthalpy,
    EquilibriumIsobaricHeatCapacity
};

struct SolvedStateObservable final {
    SolvedStateObservableKind kind = SolvedStateObservableKind::Temperature;
    std::size_t phase = 0;
    std::size_t component = 0;
    PropertyObservable property = PropertyObservable::IdealEnthalpy;
};

// Physical-unit external controls: a free T or P component must be zero.
struct SolvedStateActionDirection final {
    double temperature_k = 0.0;
    double pressure_pa = 0.0;
    std::vector<double> feed_amounts_mol;
    std::vector<double> active_parameters;
    std::vector<std::array<double, 4>> reaction_coefficients; // increments of (a,b,c,d), evaluated at T(t)
    double enthalpy_j = 0.0;
    double entropy_j_per_k = 0.0;
};

struct SolvedStateActionRequest final {
    std::vector<SolvedStateObservable> observables;
    std::vector<SolvedStateActionDirection> directions;
};

enum class SolvedStateActionStatus {
    Available, CentralSolveUnavailable, DefiningEquationsNotSatisfied,
    PhaseBoundaryNotDifferentiable, NonuniqueCoalescedPhases, RankDeficient, NonfiniteJacobian,
    LinearSolveResidual, SupportChange, IonicInsertion, EnergyOrderUnavailable,
    DensityUnavailable, PressureSingular, AssociationSingular, AssociationIllConditioned,
    AssociationResidual, MissingThermochemistry, ReferenceUnavailable,
    OutsideInterval, ThermalBoundary,
    CentralPropertyUnavailable, Nonfinite, IllConditioned
};

// Per-species ideal-gas h (J/mol) and, on request, cp (J/(mol K)) at the standard pressure
// of a referenced problem (EqID: reference_species_calorics); a refusal clears both.
struct ReferenceCalorics final {
    struct Column final {
        std::vector<double> values;
        SolvedStateActionStatus status = SolvedStateActionStatus::Available;
        std::string message;
    };
    struct Derivative final {
        std::size_t temperature_order = 0, pressure_order = 0;
        std::optional<std::size_t> parameter;
        Column column; // raw derivatives of the record-anchored density-standard g, J/mol
    };
    std::vector<Derivative> gibbs;
    std::vector<std::vector<double>> reaction_map; // anchored inverse of nu_D by reaction
    std::vector<std::vector<double>> reaction_defect; // I - nu_D nu_D^+, admits combined source motions
    bool thermal_boundary = false;
    SolvedStateActionStatus status = SolvedStateActionStatus::Available;
    std::string message;
};
ReferenceCalorics reference_calorics(const CompiledSystem &compiled, double temperature, double pressure,
                                     bool heat_capacity, bool pressure_slope = false,
                                     bool parameter_slopes = false, std::size_t derivative_order = 2);

struct SolvedStateActionResult final {
    SolvedStateObservable observable;
    std::string unit;
    std::size_t order = 0;
    std::optional<double> value;
    std::optional<double> action;
    SolvedStateActionStatus status = SolvedStateActionStatus::CentralSolveUnavailable;
    std::optional<PropertyActionDiagnostics> property_diagnostics;
    std::optional<std::string> diagnostic_message;
};

struct SolvedStateActionDiagnostics final {
    std::size_t equation_count = 0;
    std::size_t state_dimension = 0;
    std::size_t maximum_energy_order = 0;
    std::size_t state_factorizations = 0;
    std::size_t linear_solves = 0;
    std::optional<std::size_t> rank;
    std::optional<double> condition_estimate;
    std::optional<double> maximum_linear_residual;
    std::optional<double> central_residual;
    std::optional<double> central_floor;
};

struct SolvedStateActionBatch final {
    std::vector<SolvedStateActionResult> results;
    SolvedStateActionDiagnostics diagnostics;
};

inline const char* solved_state_status_name(SolvedStateActionStatus status) noexcept {
    switch (status) {
    case SolvedStateActionStatus::Available: return "available";
    case SolvedStateActionStatus::CentralSolveUnavailable: return "central_solve_unavailable";
    case SolvedStateActionStatus::DefiningEquationsNotSatisfied: return "defining_equations_not_satisfied";
    case SolvedStateActionStatus::PhaseBoundaryNotDifferentiable: return "phase_boundary_not_differentiable";
    case SolvedStateActionStatus::NonuniqueCoalescedPhases: return "nonunique_coalesced_phases";
    case SolvedStateActionStatus::RankDeficient: return "rank_deficient";
    case SolvedStateActionStatus::NonfiniteJacobian: return "nonfinite_jacobian";
    case SolvedStateActionStatus::LinearSolveResidual: return "linear_solve_residual";
    case SolvedStateActionStatus::SupportChange: return "support_change";
    case SolvedStateActionStatus::IonicInsertion: return "ionic_insertion";
    case SolvedStateActionStatus::EnergyOrderUnavailable: return "energy_order_unavailable";
    case SolvedStateActionStatus::DensityUnavailable: return "density_unavailable";
    case SolvedStateActionStatus::PressureSingular: return "pressure_singular";
    case SolvedStateActionStatus::AssociationSingular: return "association_singular";
    case SolvedStateActionStatus::AssociationIllConditioned: return "association_ill_conditioned";
    case SolvedStateActionStatus::AssociationResidual: return "association_residual";
    case SolvedStateActionStatus::MissingThermochemistry: return "missing_thermochemistry";
    case SolvedStateActionStatus::ReferenceUnavailable: return "reference_unavailable";
    case SolvedStateActionStatus::OutsideInterval: return "outside_interval";
    case SolvedStateActionStatus::ThermalBoundary: return "thermal_boundary";
    case SolvedStateActionStatus::CentralPropertyUnavailable: return "central_property_unavailable";
    case SolvedStateActionStatus::Nonfinite: return "nonfinite";
    case SolvedStateActionStatus::IllConditioned: return "ill_conditioned";
    }
    return "unknown";
}

inline SolvedStateActionStatus solved_state_status(PropertyActionStatus status) noexcept {
    switch (status) {
    case PropertyActionStatus::Available: return SolvedStateActionStatus::Available;
    case PropertyActionStatus::MissingThermochemistry: return SolvedStateActionStatus::MissingThermochemistry;
    case PropertyActionStatus::OutsideInterval: return SolvedStateActionStatus::OutsideInterval;
    case PropertyActionStatus::ThermalBoundary: return SolvedStateActionStatus::ThermalBoundary;
    case PropertyActionStatus::SupportChange: return SolvedStateActionStatus::SupportChange;
    case PropertyActionStatus::IonicInsertion: return SolvedStateActionStatus::IonicInsertion;
    case PropertyActionStatus::AssociationSingular: return SolvedStateActionStatus::AssociationSingular;
    case PropertyActionStatus::AssociationIllConditioned: return SolvedStateActionStatus::AssociationIllConditioned;
    case PropertyActionStatus::AssociationResidual: return SolvedStateActionStatus::AssociationResidual;
    case PropertyActionStatus::DensityUnavailable: return SolvedStateActionStatus::DensityUnavailable;
    case PropertyActionStatus::PressureSingular: return SolvedStateActionStatus::PressureSingular;
    case PropertyActionStatus::CentralPropertyUnavailable: return SolvedStateActionStatus::CentralPropertyUnavailable;
    case PropertyActionStatus::Nonfinite: return SolvedStateActionStatus::Nonfinite;
    }
    return SolvedStateActionStatus::Nonfinite;
}

std::string solved_state_observable_id(const SolvedStateObservable& observable);

// Exact first-, second- or third-order mixed directional actions of physical solved-state
// observables at one accepted result; no nonlinear solve, no compilation, no mutation.
SolvedStateActionBatch solved_state_actions(
    const CompiledSystem& compiled,
    const EquilibriumResult& central,
    const SolvedStateActionRequest& request
);
// Independent A1 directions at one accepted center; each response retains its own status.
std::vector<SolvedStateActionBatch> solved_state_first_actions(
    const CompiledSystem& compiled, const EquilibriumResult& central,
    const std::vector<SolvedStateObservable>& observables,
    const std::vector<SolvedStateActionDirection>& directions
);

}  // namespace epcsaft::equilibrium
