#pragma once

#include "core/convergence.hpp"
#include "core/linear_error.hpp"
#include "properties/helmholtz_model.hpp"

#include <cstddef>
#include <chrono>
#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>
#include <numeric>
#include <map>
#include <memory>
#include <string>
#include <string_view>
#include <tuple>
#include <utility>
#include <vector>

namespace epcsaft::core { class Scope; }
namespace epcsaft::equilibrium {

// ALG-EQUILIBRIUM: scaled fresh-Newton solve admission, infinity if unresolved.
double polish_forward_bound(const std::vector<double>& jacobian, const std::vector<double>& rows,
    const std::vector<double>& delta, const std::vector<double>& scale);

enum class PhaseKind { Liquid, Vapor };
enum class AmountKind { Free, Pinned };
enum class FeedKind { MoleFractions, Amounts, MassFractions, Molalities };

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
    // Calculated and supplied starts retain their compositions when reused as guesses.
    bool declared_composition_guess = true;
    std::optional<double> packing_guess;

    bool is_pinned() const noexcept { return amount.kind == AmountKind::Pinned; }
    bool is_zero_pinned() const noexcept { return is_pinned() && amount.fraction == 0.0; }
    bool is_free() const noexcept { return amount.kind == AmountKind::Free; }
};

struct FeedSpec final {
    FeedKind kind = FeedKind::Amounts;
    std::vector<double> values;
    std::vector<std::size_t> solvent;
    double solvent_mass_kg = 1.0;
};

// EOS correlation identity: activity (rho x / (1 mol/m^3)) exp(mu_res / RT).
inline constexpr const char* EOS_STANDARD_STATE_ID = "eos-helmholtz-coordinate-basis";

inline constexpr const char* IDEAL_MOLE_FRACTION_STANDARD_STATE_ID = "ideal-mole-fraction";
inline constexpr const char* IDEAL_MOLARITY_STANDARD_STATE_ID = "ideal-molarity-1-mol-per-litre";
struct IdealMoleFraction final {
    std::vector<std::string> component_ids;
    std::vector<int> charge_number;
    policy::SolverOptions solver_options;
    std::string activity_basis;
    std::optional<double> solution_volume_m3;
    std::optional<std::string> unit_activity_solvent;
    std::optional<std::size_t> solvent_index;
    IdealMoleFraction(std::vector<std::string> ids, std::vector<int> charges,
                      policy::SolverOptions options = {}, std::optional<std::string> basis = {},
                      std::optional<double> volume = {}, std::optional<std::string> solvent = {});
    const char* standard_state_id() const {
        return activity_basis == "molarity" ? IDEAL_MOLARITY_STANDARD_STATE_ID : IDEAL_MOLE_FRACTION_STANDARD_STATE_ID;
    }
    // EqID ideal_molarity_reaction; log_amount is ln(n_i / mol), V is in m^3.
    template<class Scalar> Scalar molarity_log_activity(const Scalar& log_amount, std::size_t i) const {
        return solvent_index == i ? Scalar(0.0) : Scalar(log_amount - std::log(*solution_volume_m3) - std::log(1000.0));
    }
};

// Exactly one target: scalar K, affinity in J/mol, or the P3 ln K(T) correlation on the
// EOS density basis. No target silently overrides another.
struct ReactionSpec final {
    std::string name;
    std::vector<double> stoichiometry;
    std::optional<double> equilibrium_constant;
    std::optional<double> affinity_j_per_mol;
    std::optional<ReactionLogPolynomial> correlation;
    std::optional<ReactionReference> reference;
    bool uses_species_thermochemistry() const { return !equilibrium_constant && !affinity_j_per_mol && !correlation; }
};

// EqID solid_standard_gibbs: a pure stoichiometric solid S -> sum_i nu_i A_i (charge neutral) with exactly one
// source constant, its reference at the source pressure P0, and V_s in m^3/mol; Free or Pinned(0) amount.
struct SolidSpec final {
    std::string name;
    std::vector<double> stoichiometry;
    std::optional<double> equilibrium_constant;
    std::optional<ReactionLogPolynomial> correlation;
    std::optional<ReactionReference> reference;
    double molar_volume_m3_per_mol = 0.0;
    AmountSpec amount;
};

// Exact declaration equality, shared by fitting admission/grouping and continuation.
inline bool same_reaction(const ReactionSpec& x, const ReactionSpec& y) {
    const auto correlation = [](const auto& v) {
        return std::tie(v.a, v.b, v.c, v.d, v.reference_temperature, v.temperature_min,
                        v.temperature_max, v.standard_state_id, v.e, v.f);
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

// Critical-equation cell seeds (EqID critical_equations): molar density [mol/m^3] and an optional
// unnormalized amount-perturbation direction u; empty u selects the lowest mode of M at the seed.
struct CriticalDeclaration final {
    double density_guess_mol_per_m3 = 0.0;
    std::vector<double> mode_guess;
};

// EqID envelope_local_condition: ln(T/1 K) - s, ln(P/1 Pa) - s or ln x_k^inc - ln x_k^feed - s.
enum class EnvelopeCoordinate { LogTemperature, LogPressure, LogK };
struct EnvelopeCondition final {
    EnvelopeCoordinate kind = EnvelopeCoordinate::LogPressure;
    std::size_t component = 0;
    double target = 0.0;
};

// EqID finite_addition_enthalpy: n(alpha) = n0 + alpha * direction, at target bubble pressure.
struct AffineLoading final {
    std::vector<double> initial, direction;
    double lower = 0, upper = 0, seed = 0, pressure = 0;
};

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
    std::optional<CriticalDeclaration> critical;  // one free neutral phase, T and P free
    std::optional<EnvelopeCondition> envelope;    // feed plus a zero-pinned incipient phase
    std::optional<AffineLoading> loading;
    bool azeotrope = false;  // EqID azeotrope_condition: a binary bubble/dew cell whose feed composition is the seed
    // Caller row criterion for the solve; the named policy value is the only default.
    policy::SolverOptions solver_options;
    std::vector<SolidSpec> solids;
};

enum class RowKind {
    Chemical, Pressure, Moiety, Electro, Normalization, Enthalpy, Entropy, ReactionMassAction,
    CriticalMode, CriticalNormalization, CriticalCubic, EnvelopeCondition, AzeotropeCondition, SolidAffinity, ReactionMolarity, SaturationPressure
};

// Rows a Gibbs cell classifies at its minimum; they are never optimizer constraints there.
inline bool stationarity_row(RowKind kind) noexcept {
    return kind == RowKind::Chemical || kind == RowKind::Pressure || kind == RowKind::SolidAffinity;
}

inline const char* row_kind_name(RowKind kind) noexcept {
    switch (kind) {
    case RowKind::SaturationPressure: return "saturation_pressure";
    case RowKind::ReactionMassAction: return "ideal_mole_fraction_reaction";
    case RowKind::ReactionMolarity: return "ideal_molarity_reaction";
    case RowKind::Chemical: return "chemical_potential_stationarity";
    case RowKind::Pressure: return "pressure_equality";
    case RowKind::Moiety: return "material_balance";
    case RowKind::Electro: return "phase_electroneutrality";
    case RowKind::Normalization: return "phase_normalization";
    case RowKind::Enthalpy: return "extensive_enthalpy_specification";
    case RowKind::Entropy: return "extensive_entropy_specification";
    case RowKind::CriticalMode: case RowKind::CriticalNormalization: case RowKind::CriticalCubic:
        return "critical_equations";
    case RowKind::EnvelopeCondition: return "envelope_local_condition";
    case RowKind::AzeotropeCondition: return "azeotrope_condition";
    case RowKind::SolidAffinity: return "solid_saturation_index";
    }
    return "unknown";
}

struct Row final {
    RowKind kind = RowKind::Chemical;
    std::size_t phase = 0;
    std::size_t index = 0;
};

// max |v_i|, infinite if any entry is nonfinite.
template <typename Base>
Base maximum_abs(const std::vector<Base>& values) noexcept {
    Base result = 0.0;
    for (const Base value : values) {
        const Base magnitude = std::abs(value);
        if (!std::isfinite(magnitude)) return std::numeric_limits<Base>::infinity();
        result = std::max(result, magnitude);
    }
    return result;
}

template <typename Base>
struct RowUncertaintyValues final {
    std::vector<Base> terms, inner;
    std::string cause;
};
using RowUncertainty = RowUncertaintyValues<double>;

template <typename Base>
struct PointEvaluationValues final {
    std::vector<Base> values;
    Base objective = 0.0;
    Base floor = 0.0;
};
using PointEvaluation = PointEvaluationValues<double>;

struct Layout final {
    std::vector<std::size_t> x;
    std::optional<std::size_t> amount;
    std::size_t volume = 0;
};

// The fluid packing domain shared by seeds, prepared points and stability states.
inline bool admitted_packing_fraction(double eta) {
    return std::isfinite(eta) && eta >= policy::PACKING_FRACTION_FLOOR.value
        && eta <= policy::PACKING_FRACTION_CEILING.value;
}

template <typename Base>
struct PackingTermsValues final {
    Base coefficient = 0.0;
    Base temperature_derivative = 0.0;
    std::vector<Base> component;
};
using PackingTerms = PackingTermsValues<double>;

// derivative[k][i] is d^k (g_i / RT) / dT^k at the prepared temperature.  The
// pressure fields are the matching first and second derivatives at fixed T.
struct StandardTerm final {
    struct NumericalError final {
        std::array<std::vector<double>, 2> absolute;
        std::string cause;
    };
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
template <typename Base>
struct PhasePointValues final {
    std::vector<Base> point;
    std::vector<Layout> shape;
    Base temperature = 0.0;
    Base pressure = 0.0;
    std::vector<std::vector<Base>> composition;
    std::vector<Base> amount;
    std::vector<Base> volume;
    std::vector<Base> eta;
    std::vector<PackingTermsValues<Base>> packing;
    std::vector<HelmholtzDerivativeValues<Base>> derivative;
    std::vector<PureIsothermalValues<Base>> pure;
    StandardTerm standard;
    // A referenced problem's standard term with jets along every active parameter, built on the
    // first parameter-driven row series at this point and reused by every later direction.
    mutable std::optional<StandardTerm> parameter_standard;
    mutable std::map<std::string, ResidualEnergyAction> energy_actions;
};
using PhasePoint = PhasePointValues<double>;

// Taylor series of one phase along a path: point[k] is the k-th coefficient of the phase as
// the model sees it, (T, v, x, active parameters); gradient[k] the matching coefficient of
// the extensive residual Helmholtz gradient; amount[k] and packing[k] the phase amount and
// packing fraction.
template <typename Base>
struct PhaseSeriesValues final {
    std::vector<std::vector<Base>> point, gradient;
    std::vector<Base> amount, packing;
};
using PhaseSeries = PhaseSeriesValues<double>;

// Defining rows and physical state along y(t) = sum_k y_k t^k with controls
// q(t) = q0 + t dq, dq = [T, P, n_feed (C), active parameters (A)]; a free T or P entry of
// dq is unused. `energy` is the worst status of the model's residual contractions.
struct RowSeries final {
    std::vector<std::vector<double>> rows;
    std::vector<double> temperature, pressure, loading;
    std::vector<PhaseSeries> phase;
    std::array<std::vector<std::vector<double>>, 3> reaction_shift; // L, L_T, L_TT Taylor coefficients
    std::vector<std::vector<double>> solid_amount, solid_affinity;  // [solid][k]: n_s and A_s/RT
    PropertyActionStatus energy = PropertyActionStatus::Available;
    PropertyActionDiagnostics energy_diagnostics;
    std::size_t maximum_energy_order = 0;
};

void restrict_composition_guess(PhaseSpec& phase, const std::vector<std::size_t>& support);

// A compiled solid: its amount column when free, B nu_s, ln K^rho(T, P0) and its T slope where admitted. A
// present solid's row is solid_saturation_index; an absent one's is n_s = 0 (ALG-EQUILIBRIUM classification).
struct CompiledSolid final {
    std::optional<std::size_t> column;
    std::vector<double> balance;
    double log_k = 0.0;
    std::optional<double> log_k_slope;
    bool present = false;
    std::optional<NeutralReferenceReactionDiagnostic> reference;
};

struct CompiledSystem final {
    std::optional<core::Convergence> first_numerical_failure;
    std::vector<core::NumericalRank> rank_records;
    std::optional<core::NumericalRank> last_rank;
    std::vector<double> optimizer_identity() const;
    const HelmholtzModel* model = nullptr;
    const IdealMoleFraction* ideal = nullptr;
    // EqID ideal_reaction_amount_coordinates: column j is ln(n_support[j]/N_feed).
    std::vector<std::size_t> amount_coordinates;
    std::vector<double> ideal_log_k;
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
    std::vector<std::pair<std::size_t, std::size_t>> log_charge_rows;  // numerical row, phase
    std::size_t pinned_count = 0;
    std::size_t mode_count = 0;  // critical eigenvector coordinates u, just before the free intensives
    std::size_t free_intensive_count = 0;
    std::size_t rank = 0;
    std::size_t charge_row_count = 0;
    std::size_t primal_size = 0;
    std::size_t generic_nullity = 0;
    bool feasible = false;
    std::vector<std::size_t> removed_species;
    double feed_total = 0.0;
    bool has_ions = false;
    bool has_reaction = false;
    bool conserves_total_moles = true;
    std::size_t jacobian_rank = 0;
    std::vector<NeutralReferenceReactionDiagnostic> reference_diagnostics;
    std::vector<CompiledSolid> solids;

    // EqIDs solid_standard_gibbs and solid_saturation_index: g_s/RT and A_s/RT = lambda.B nu_s - g_s/RT.
    double solid_gibbs(const PhasePoint& prepared, std::size_t solid) const;
    double solid_affinity(const PhasePoint& prepared, std::size_t solid) const;

    // State-only preparation differentiates V, support amounts and T when free (phase Hessians);
    // fixed controls and active parameters are padding. Sensitivity and lever callers need full.
    PhasePoint prepare(const std::vector<double>& point, bool state_only = false,
        std::optional<core::Convergence>* first_numerical_failure = nullptr) const;
    template <typename Base>
    PhasePointValues<Base> prepare_values(const std::vector<Base>&, bool,
        std::optional<core::Convergence>* = nullptr) const;
    template <typename Base> PointEvaluationValues<Base> evaluate_at(const PhasePointValues<Base>&) const;
    template <typename Base> std::vector<Base> jacobian_at(const PhasePointValues<Base>&) const;
    template <typename Base> RowUncertaintyValues<Base> row_uncertainty_at(const PhasePointValues<Base>&,
        const std::vector<std::vector<Base>>&) const;
    PointEvaluation evaluate(const std::vector<double>& point) const;
    // Residuals that classify a solve: for a free-amount (Gibbs) cell the chemical and
    // pressure rows are scaled by N_p x_i / N_feed and N_p / N_feed, the first-order
    // conditions of the constrained Gibbs minimum in the solved chart; a pinned
    // (feasibility) cell classifies its rows unscaled. Trace moiety rows are relative to their
    // size in both routes (trace_moiety_relative); charged Vapor also retains the log-charge criterion.
    std::vector<double> classified_residuals(const PhasePoint& prepared, const std::vector<double>& rows) const;
    // Moiety rows (Jacobian rows for width n) divided by their size, the larger of the feed inventory and the
    // row's term magnitude per mol of feed, when it is at most EQUILIBRIUM_TRACE_MOIETY_FRACTION (#633).
    template <typename Base> std::vector<Base> trace_moiety_relative(const PhasePointValues<Base>&, std::vector<Base>, std::size_t = 1) const;
    std::vector<double> trace_moiety_relative(const PhasePoint& prepared, std::vector<double> values,
                                              std::size_t width = 1) const;
    std::vector<double> numerical_rows(const PhasePoint&, std::vector<double>, std::size_t width = 1) const;
    double classified_residual_norm(const std::vector<double>& point) const;
    RowUncertainty row_uncertainty(const PhasePoint&) const;
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
    // Row and physical Taylor coefficients through `order` (at most 3; 2 in a critical cell) along the state
    // coefficients `state` (state[0] is the prepared point) and the control direction.
    RowSeries row_series(
        const PhasePoint& prepared, const std::vector<std::vector<double>>& state,
        const std::vector<double>& control, std::size_t order
    ) const;

    std::size_t loading_column() const noexcept { return primal_size - free_intensive_count - 1; }
    template<class Scalar> std::vector<Scalar> feed_at(const std::vector<Scalar>& point) const {
        std::vector<Scalar> feed(n_feed.begin(), n_feed.end());
        if (problem.loading) for (std::size_t i = 0; i < feed.size(); ++i)
            feed[i] = problem.loading->initial[i] + point[loading_column()] * problem.loading->direction[i];
        return feed;
    }
    double feed_total_at(const PhasePoint& at) const {
        const auto feed = feed_at(at.point);
        return std::accumulate(feed.begin(), feed.end(), 0.0);
    }

    std::size_t unknown_count() const noexcept {
        return ideal ? primal_size : primal_size + rank + charge_row_count;
    }


    bool has_caloric_condition() const noexcept {
        return problem.enthalpy_j.has_value() || problem.entropy_j_per_k.has_value();
    }

    bool feasibility_cell() const noexcept {
        return ideal || pinned_count != 0 || has_caloric_condition() || problem.critical.has_value();
    }
    // One minimized phase: no optimizer ln eta bounds, prepare refuses eta outside the domain (#314).
    bool log_packing_unbounded() const noexcept { return problem.phases.size() == 1 && !feasibility_cell(); }

    // Scalar conditions e: a caloric row, two net critical equations or an envelope row.
    std::size_t condition_count() const noexcept {
        return has_caloric_condition() + 2 * problem.critical.has_value() + problem.envelope.has_value();
    }

    std::size_t mode_offset() const noexcept { return primal_size - free_intensive_count - mode_count; }

    // Control direction [T, P, n_feed (C), active parameters (A)]; the extended layout appends the
    // (a, b, c, d) increments of every reaction law and then the H and S targets.
    std::size_t control_width() const noexcept { return 2 + n_feed.size() + (model ? model->active_parameters().size() : 0); }
    std::size_t extended_control_width() const noexcept { return control_width() + 4 * problem.reactions.size() + 2; }

    // A free amount is a coordinate seeded at 1 mol; a pinned amount is its fraction of the feed.
    double phase_amount(const PhaseSpec& phase) const noexcept { return phase.is_free() ? 1.0 : phase.amount.fraction * feed_total; }

    // The declared value, else the default seed of an unseeded free intensive: 300 K and 1e5 Pa.
    double seed_temperature() const noexcept { return problem.temperature.is_positive() ? problem.temperature.value : 300.0; }
    double seed_pressure() const noexcept { return problem.pressure.is_positive() ? problem.pressure.value : 1.0e5; }

};

CompiledSystem compile_problem(
    const HelmholtzModel& model, const EquilibriumProblem& problem
);
CompiledSystem compile_problem(const IdealMoleFraction&, const EquilibriumProblem&);
// The classified declaration of a solved point: a free solid is present exactly when its amount is positive.
CompiledSystem classify_solids(CompiledSystem compiled, const std::vector<double>& point);

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

// Spectrum of M, the cubic form and the augmented-Jacobian diagnostics of a critical-equation
// candidate; it certifies no stability.
struct CriticalRecord final {
    std::vector<double> eigenvalues, mode;
    double cubic = NAN, zero_band = NAN, condition_estimate = NAN, backward_error = NAN;
    std::size_t rank = 0;
};

struct EquilibriumNumericalEstimates final {
    std::string basis = "canonical_defining_rows", cause = "not_evaluated";
    std::string source_reference_method;
    std::size_t rank = 0, dimension = 0, precision_bits = std::numeric_limits<double>::digits;
    double threshold = 0, condition = INFINITY, inverse_norm = INFINITY, backward_error = INFINITY, physical_radius = INFINITY;
    std::vector<double> residuals, row_terms, row_inner, row_floor, physical_scale;
    core::NumericalEstimate temperature, pressure, loading;
    std::vector<core::NumericalEstimate> densities, packing, composition, amounts;
};

template <typename Base>
struct PhysicalProjection final {
    using Vector = Eigen::Matrix<Base, Eigen::Dynamic, 1>;
    Eigen::Matrix<Base, Eigen::Dynamic, Eigen::Dynamic> gradient;
    Vector values, scale;
};
template <typename Base> PhysicalProjection<Base> physical_projection(const CompiledSystem&, const PhasePointValues<Base>&);
template <typename Base, typename Vector>
Base physical_radius(const CompiledSystem& c, const PhasePointValues<Base>& at, const Vector& error) {
    Base radius = error[1] / std::max(std::abs(at.pressure), Base(1));
    for (Eigen::Index i = 2; i < error.size(); ++i)
        radius = std::max(radius, error[i] / ((i - 2) % (3 + c.n_feed.size()) == 2 ? Base(c.feed_total) : Base(1)));
    return radius;
}
template <typename Base, typename Jacobian, typename Factor, typename Correction>
auto physical_error(const CompiledSystem& c, const PhasePointValues<Base>& at,
    const Jacobian& j, const Eigen::Matrix<Base, Eigen::Dynamic, 1>& f,
    const Factor& lu, const Correction& d, const Eigen::Matrix<Base, Eigen::Dynamic, 1>& column) {
    using Vector = Eigen::Matrix<Base, Eigen::Dynamic, 1>;
    const auto projection = physical_projection(c, at);
    std::vector<std::vector<Base>> inner;
    for (const auto& derivative : at.derivative) inner.push_back(derivative.gradient_error.empty()
        ? std::vector<Base>(derivative.gradient.size(), derivative.gradient_floor) : derivative.gradient_error);
    const auto u = c.row_uncertainty_at(at, inner);
    const Vector terms = Eigen::Map<const Vector>(u.terms.data(), f.size()), errors = Eigen::Map<const Vector>(u.inner.data(), f.size());
    const Vector eta = core::canonical_row_floor(j, projection.scale, terms, errors), row = Vector::Ones(f.size());
    return core::linear_error(j, f, projection.gradient, eta, row, column, &lu, &d);
}

struct OptimizerState;
struct LinearSolveDiagnostics final {
    std::string system = "one_component_two_density_reduction", reason;
    double residual_norm = std::numeric_limits<double>::infinity(), evaluation_floor = std::numeric_limits<double>::infinity();
    double jacobian_inverse_norm = std::numeric_limits<double>::infinity(), condition_estimate = std::numeric_limits<double>::infinity();
    double forward_error_estimate = std::numeric_limits<double>::infinity();
    std::size_t precision_bits = std::numeric_limits<double>::digits, refinement_iterations = 0;
};

struct IpoptIteration final {
    int mode = 0, iteration = 0, line_search_trials = 0;
    double objective = 0, inf_pr = 0, inf_du = 0, mu = 0, step_norm = 0, regularization = 0, alpha_du = 0, alpha_pr = 0;
};
struct IpoptPoint final {
    std::vector<double> tnlp_coordinates, coordinate_scale, tnlp_lower_multipliers, tnlp_upper_multipliers, tnlp_constraint_multipliers;
    std::shared_ptr<const EquilibriumProblem> declaration;
    std::vector<Row> rows;
    std::vector<std::size_t> constraint_rows;
    double objective_factor = 0;
};
struct IpoptDiagnostics final {
    std::optional<IpoptIteration> last;
    std::optional<int> first_restoration_iteration, initialization_return, application_return;
    std::size_t restoration_entries = 0, callbacks = 0;
    bool warm_initialized = false, predicted = false, physical_amounts = false;
    std::optional<IpoptPoint> failed_point;
};
struct EquilibriumResult final {
    std::shared_ptr<const EquilibriumResult> rejected_reduction;
    std::optional<core::Convergence> terminal_numerical_failure, last_full_evaluation;
    std::vector<double> phase_pressures;
    std::optional<core::Convergence> first_numerical_failure;
    std::vector<core::NumericalRank> rank_records;
    std::optional<IpoptDiagnostics> ipopt_diagnostics, recovered_ipopt_diagnostics;
    std::optional<core::NumericalRank> last_rank;
    std::map<std::string, double> solver_settings;
    std::string diagnostic_options_file, linear_solver;  // linear_solver: empty when no Ipopt solve ran
    std::shared_ptr<const OptimizerState> optimizer_state;
    const HelmholtzModel* model_identity = nullptr;
    const IdealMoleFraction* ideal_identity = nullptr;
    std::shared_ptr<const EquilibriumProblem> declaration;
    std::shared_ptr<DiscoveryRecord> discovery;
    bool success = false;
    std::string message;
    std::string solver_status;
    std::optional<double> temperature;
    std::optional<double> pressure;
    std::vector<double> amounts;
    std::vector<double> mole_fractions;
    std::vector<double> molar_densities;
    std::vector<double> packing_fractions;
    std::vector<double> multipliers;
    std::vector<double> constraint_multipliers;
    std::vector<double> residuals;
    // The residuals the classification actually used (see CompiledSystem::classified_residuals).
    std::vector<double> classified_residuals;
    std::vector<double> coordinates;
    std::size_t phase_count = 0;
    std::size_t iterations = 0;
    std::size_t evaluations = 0;
    std::size_t hessian_calls = 0;
    double objective = 0.0;
    double reduced_hessian_minimum = 0.0;
    core::Convergence convergence;
    std::optional<LinearSolveDiagnostics> linear_solve_diagnostics;
    std::optional<EquilibriumNumericalEstimates> numerical_estimates;
    std::shared_ptr<const PhasePoint> prepared_endpoint;
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
    std::optional<CriticalRecord> critical;
    std::vector<double> solid_amounts, solid_affinity_decades;  // per declared solid: mol, and SI_s in decades

    bool has_state() const noexcept {
        return success && temperature && pressure
            && !amounts.empty() && !mole_fractions.empty();
    }

    double residual_norm() const noexcept { return maximum_abs(residuals); }
};

EquilibriumNumericalEstimates equilibrium_numerical_estimates(const CompiledSystem&, const EquilibriumResult&, const PhasePoint* = nullptr);
struct DriverLimits;
EquilibriumNumericalEstimates equilibrium_numerical_estimates(const CompiledSystem&, const EquilibriumResult&, const std::vector<long double>&, const DriverLimits* = nullptr);
EquilibriumNumericalEstimates equilibrium_numerical_estimates(const HelmholtzModel&, const EquilibriumResult&);
EquilibriumNumericalEstimates equilibrium_numerical_estimates(const IdealMoleFraction&, const EquilibriumResult&);

// EqID reaction_activity_quotient; per declared phase/reaction, in component order.
struct ReactionActivities final {
    std::size_t phase_index = 0, reaction_index = 0;
    std::string basis, cause;
    std::optional<double> reference_pressure_pa, log_quotient;
    SpeciesActivities species;
};
std::vector<ReactionActivities> reaction_activities(const HelmholtzModel&, const EquilibriumResult&);
std::vector<ReactionActivities> reaction_activities(const IdealMoleFraction&, const EquilibriumResult&);

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
    std::optional<core::Convergence> first_numerical_failure;
    std::optional<core::NumericalRank> last_rank;
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
std::vector<double> reconstruct_trial_composition(std::vector<double>, const std::vector<std::size_t>&, const std::vector<int>&, PhaseKind = PhaseKind::Liquid);
// n dln(phi_i)/dn_j at fixed T, P on the support (row-major), from an IsothermalPhaseHessian evaluation at unit amount and v.
std::vector<double> tangent_plane_jacobian(const HelmholtzDerivatives&, double, const std::vector<std::size_t>&);
std::pair<EquilibriumProblem, std::vector<double>> remove_discovery_phase(const CompiledSystem&,
    const std::vector<double>&, std::map<std::string, double>* = nullptr);
StabilityReport stability_test(const CompiledSystem&, const EquilibriumResult&, const StabilityOptions& = {});
std::pair<EquilibriumResult, std::vector<StabilityReport>> cooling_onset(const CompiledSystem&);
EquilibriumResult discover_phases(const HelmholtzModel&, const EquilibriumProblem&);

EquilibriumResult solve_equilibrium(
    const HelmholtzModel& model,
    const EquilibriumProblem& problem,
    const std::optional<EquilibriumResult>& start = std::nullopt
);
EquilibriumResult solve_equilibrium(const IdealMoleFraction&, const EquilibriumProblem&,
    const std::optional<EquilibriumResult>& start = std::nullopt);
std::map<std::string, double> effective_water_amounts(const IdealMoleFraction&, const EquilibriumResult&,
    const std::string& water, const std::vector<std::string>& ions, const std::vector<std::string>& molecular_species);
// The same solve of a compiled problem (its model must outlive it), inside the caller's open span if given.
// `recover`: the public no-start solve and cooling onset only (ALG-EQUILIBRIUM neutral recovery, #565).
struct DriverLimits;
EquilibriumResult solve_compiled(const CompiledSystem& compiled, const std::optional<EquilibriumResult>& start = std::nullopt,
                                 core::Scope* outer = nullptr, const DriverLimits* limits = nullptr,
                                 std::shared_ptr<const OptimizerState> warm = {}, bool recover = false, EquilibriumResult* audit = nullptr, bool select_endpoint = false);

EquilibriumResult solve_saturation(const CompiledSystem&);
struct SolvedStateActionDiagnostics;
// All operands are generated from one exact trial declaration in declared phase order.
struct SaturationOperand final {
    core::Quantity quantity = core::Quantity::Linear;
    double value = NAN;
    core::NumericalEstimate error;
};
struct SaturationCandidate final {
    EquilibriumResult state;
    std::vector<SaturationOperand> operands;
    bool admitted = false;
};
struct SaturationComparison final {
    std::size_t left = 0, right = 0;
    std::vector<core::EquilibriumComparison> operands;
    std::string identity;
};
struct SaturationNode final {
    double loading = 0, mismatch = NAN, slope = NAN;
    EquilibriumResult state;
    EquilibriumNumericalEstimates estimates;
    std::string action_cause;
    double action_seconds = 0, estimate_seconds = 0;
    EquilibriumNumericalEstimates local_estimates;
    std::shared_ptr<SolvedStateActionDiagnostics> action_diagnostics;
};
struct SaturationSearch final {
    std::vector<double> grid;
    std::vector<SaturationNode> nodes;
    std::vector<std::array<std::size_t, 2>> events;
    std::vector<SaturationCandidate> candidates;
    std::vector<SaturationComparison> comparisons;
    std::vector<std::string> causes;
    std::vector<std::vector<std::size_t>> groups;
    std::optional<std::size_t> accepted;
    std::optional<core::Convergence> first_numerical_failure;
    std::size_t brackets = 0, iterations = 0, evaluations = 0, hessian_calls = 0, attempts = 0;
    bool complete = false;
    AffineLoading declaration;
    std::array<std::size_t, 4> work_at_entry{};
    std::string uniqueness_scope = "unique within the declared search";
    std::string detection_limits = "unsampled tangencies, paired crossings or multiple roots in a cell, turns between nodes and undiscovered branches can escape detection";
};
void compare_loading_candidates(SaturationSearch&, core::ComparisonEligibility);
void admit_loading_search(SaturationSearch&);
SaturationSearch search_loading(const HelmholtzModel&, EquilibriumProblem, const AffineLoading&,
    const std::vector<double>&, const DriverLimits* = nullptr,
    std::chrono::steady_clock::time_point = std::chrono::steady_clock::time_point::max());

// ---- Solved-state derivative actions: native records, no public promise.

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
    EquilibriumIsobaricHeatCapacity,
    SolidAmount,  // `phase` is the solid ordinal
    SolidSaturationIndex,
    Loading
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
    std::string numerical_error_cause;
    struct Column final {
        std::vector<double> values, numerical_error;
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
                                     bool parameter_slopes = false, std::size_t derivative_order = 2, bool numerical_errors = false);

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
// Independent first-order directions at one accepted center; each response retains its own status.
std::vector<SolvedStateActionBatch> solved_state_first_actions(
    const CompiledSystem& compiled, const EquilibriumResult& central,
    const std::vector<SolvedStateObservable>& observables,
    const std::vector<SolvedStateActionDirection>& directions
);

}  // namespace epcsaft::equilibrium
