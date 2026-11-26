#pragma once
#include "core/export.hpp"

#include "core/convergence.hpp"
#include "core/linear_error.hpp"
#include "core/retained.hpp"
#include "properties/helmholtz_model.hpp"

#include <cstddef>
#include <exception>
#include <chrono>
#include <algorithm>
#include <cstdint>
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
#include <variant>
#include <vector>

namespace epcsaft::core { class Scope; }
namespace epcsaft::equilibrium {

enum class PhaseKind { Liquid, Vapor };
inline DensityBranch branch_of(PhaseKind kind) { return kind == PhaseKind::Vapor ? DensityBranch::Vapor : DensityBranch::Liquid; }
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
private:
    const std::shared_ptr<const char> identity_ = std::make_shared<const char>();
public:
    const auto& identity() const { return identity_; }
    IdealMoleFraction(const IdealMoleFraction&) = delete;
    const std::vector<std::string> component_ids;
    const std::vector<int> charge_number;
    const policy::SolverOptions solver_options;
    const std::string activity_basis;
    const std::optional<double> solution_volume_m3;
    const std::optional<std::string> unit_activity_solvent;
    const std::optional<std::size_t> solvent_index;
    EPC_API IdealMoleFraction(std::vector<std::string> ids, std::vector<int> charges,
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

// A solid's source ln K law as a reaction record, for the shared law comparisons below.
inline ReactionSpec solid_law(const SolidSpec& s) { return {s.name, s.stoichiometry, s.equilibrium_constant, std::nullopt, s.correlation, s.reference}; }

// Exact declaration equality, shared by fitting admission/grouping and continuation.
inline bool same_reaction(const ReactionSpec& x, const ReactionSpec& y, bool compare_temperature_bounds = true, bool structural = false) {
    const auto optional_equal = [](const auto& a, const auto& b, const auto& equal) {
        return a.has_value() == b.has_value() && (!a || equal(*a, *b));
    };
    const auto target = [structural](double a, double b) { return structural || a == b; };
    return x.name == y.name && x.stoichiometry == y.stoichiometry
        && optional_equal(x.equilibrium_constant, y.equilibrium_constant, target)
        && optional_equal(x.affinity_j_per_mol, y.affinity_j_per_mol, target)
        && optional_equal(x.reference, y.reference, [](const auto& a, const auto& b) {
            return a.reference_pressure_pa == b.reference_pressure_pa && a.source_basis == b.source_basis;
        }) && optional_equal(x.correlation, y.correlation, [&](const auto& a, const auto& b) {
            return a.standard_state_id == b.standard_state_id && (structural
                || (std::tie(a.a, a.b, a.c, a.d, a.reference_temperature, a.e, a.f)
                    == std::tie(b.a, b.b, b.c, b.d, b.reference_temperature, b.e, b.f)
                    && (!compare_temperature_bounds || std::tie(a.temperature_min, a.temperature_max) == std::tie(b.temperature_min, b.temperature_max))));
        });
}

// Explicit finite TPD search options; unset supports select the ALG-STABILITY defaults.
enum class DiscoveryMode { Full, Fast, FastTested };
struct PreparedSuppliedPoints;
class PreparedSuppliedTrials final {
public:
    std::shared_ptr<const char> model;
    double temperature;
    std::vector<PhaseSpec> trials;
    std::shared_ptr<const PreparedSuppliedPoints> points;
    std::map<std::string, double> work;
};
struct StabilityOptions final {
    DiscoveryMode mode = DiscoveryMode::Full;
    bool evaluate_supplied_trials = false;
    std::shared_ptr<const PreparedSuppliedTrials> prepared_trials;
    std::optional<std::vector<PhaseSpec>> trials;
    std::optional<std::size_t> declaration_budget;
    std::optional<std::size_t> trial_budget;
    std::vector<double> standard_uncertainty_over_rt;
    bool retry_condensed_refusals = false;  // #21: the explicit test's guarded second pass, on request only
    bool plane_value_errors = false;
};
struct DiscoveryRecord;
struct StabilityReport;

// Critical-equation cell seeds (EqID critical_equations): molar density [mol/m^3] and an optional
// unnormalized amount-perturbation direction u; empty u selects the lowest mode of M at the seed.
struct CriticalDeclaration {
    double density_guess_mol_per_m3 = 0.0;
    std::vector<double> mode_guess;
};

struct SpinodalDeclaration final : CriticalDeclaration {
    std::vector<double> density_nodes;  // caller-declared connected scan, mol/m^3
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

enum class PhaseCalculation { Unknown, Declared, Discovery };
struct EquilibriumProblem final {
    PhaseCalculation phase_calculation = PhaseCalculation::Declared;
    std::vector<PhaseSpec> phases;
    FreeValue temperature;
    FreeValue pressure;
    FeedSpec feed;
    std::optional<std::vector<double>> prescribed_liquid_composition;
    std::vector<ReactionSpec> reactions;
    std::vector<double> standard_gibbs_j_per_mol;
    std::optional<NeutralReferenceDeclaration> source_reference;
    bool topology_declared = true;
    std::optional<StabilityOptions> discovery;
    std::optional<double> enthalpy_j;
    std::optional<double> entropy_j_per_k;
    std::optional<double> volume_m3;
    std::optional<double> internal_energy_j;
    std::optional<SpinodalDeclaration> spinodal;
    std::optional<CriticalDeclaration> critical;  // one free neutral phase, T and P free
    std::optional<EnvelopeCondition> envelope;    // feed plus a zero-pinned incipient phase
    std::optional<AffineLoading> loading;
    bool azeotrope = false;  // EqID azeotrope_condition: a binary bubble/dew cell whose feed composition is the seed
    // Caller row criterion for the solve; the named policy value is the only default.
    policy::SolverOptions solver_options;
    std::vector<SolidSpec> solids;

    bool has_reaction_reference() const {
        return std::any_of(reactions.begin(), reactions.end(),
            [](const ReactionSpec& reaction) { return reaction.reference.has_value(); });
    }
};

enum class RowKind {
    Chemical, Pressure, Moiety, Electro, Normalization, PrescribedComposition, Enthalpy, Entropy, InternalEnergy, Volume, ReactionMassAction,
    CriticalMode, CriticalNormalization, CriticalCubic, SpinodalMode, SpinodalNormalization, EnvelopeCondition, AzeotropeCondition, SolidAffinity, ReactionMolarity, SaturationPressure
};

// Caloric, critical and envelope rows: analytic in record_geometry, completed by row_series.
inline bool condition_row(RowKind kind) noexcept { return kind >= RowKind::Enthalpy && kind != RowKind::SolidAffinity; }

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
    case RowKind::InternalEnergy: return "extensive_internal_energy_specification";
    case RowKind::Volume: return "extensive_volume_specification";
    case RowKind::Pressure: return "pressure_equality";
    case RowKind::Moiety: return "material_balance";
    case RowKind::Electro: return "phase_electroneutrality";
    case RowKind::Normalization: return "phase_normalization";
    case RowKind::PrescribedComposition: return "prescribed_liquid_composition";
    case RowKind::Enthalpy: return "extensive_enthalpy_specification";
    case RowKind::Entropy: return "extensive_entropy_specification";
    case RowKind::CriticalMode: case RowKind::CriticalNormalization: case RowKind::CriticalCubic:
        return "critical_equations";
    case RowKind::SpinodalMode: case RowKind::SpinodalNormalization: return "diffusive_spinodal";
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
    std::optional<std::size_t> held_column{};
};

// max |v_i|, infinite if any entry is nonfinite.
template <typename Base>
Base maximum_abs(const std::vector<Base>& values) noexcept {
    Base result = 0.0;
    for (const Base value : values) {
        const Base magnitude = wide_math::abs(value);
        if (!wide_math::isfinite(magnitude)) return std::numeric_limits<Base>::infinity();
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

// EqID solid_standard_gibbs: ln K^rho_s(T, P0) = ln_k + scale + C_s, with C_s the reference contraction at the source
// pressure P0. slope and curvature (d ln K_src/dT, d2/dT2) exist only for a correlation: a scalar K has no T law.
struct SolidLogK final {
    double ln_k = 0.0, scale = 0.0;
    std::optional<double> slope, curvature;
    NeutralReferenceReactionDiagnostic reference;
};

// The reference columns a caller asks for at (T, P): temperature derivatives through temperature_order, pressure jets,
// first-order jets in the listed active coordinates (ascending), the order of the consumed mixed reference columns,
// the parameter coordinates whose second derivatives are needed, and numerical error estimates. Each owner reads only
// the fields it uses.
struct ReferenceJets final {
    std::size_t temperature_order = 0;
    std::size_t pressure_order = 0;
    std::vector<std::size_t> parameters;
    std::size_t reference_order = 0;
    std::vector<std::size_t> pairs;
    bool errors = false;
    void request_every_parameter(const HelmholtzModel& model) {  // for a caller that reads them all
        parameters.resize(model.active_parameters().size());
        std::iota(parameters.begin(), parameters.end(), std::size_t(0));
    }
};

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
    std::vector<SolidLogK> solid;  // per solid: the curvature and reference columns the actions move; empty unless requested
    // parameter_pair_derivative[i][j][k] is d2(g_k/RT)/dtheta_a dtheta_b, a and b from pair_coordinates[i], [j].
    std::vector<std::size_t> pair_coordinates;
    std::vector<std::vector<std::vector<double>>> parameter_pair_derivative;
    std::vector<NeutralReferenceReactionDiagnostic> reference_diagnostics;
};

struct ReferenceCalorics;
struct ParameterStandardRefused final {};
struct ParameterStandard final {
    std::vector<std::size_t> positions;
    std::variant<StandardTerm, ParameterStandardRefused> outcome;
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
    bool parameter_actions = false;
    std::shared_ptr<const char> state_model;  // Owns complete fixed-parameter derivative columns.
    // The referenced standard term for the parameter positions shared by a first-action batch.
    mutable std::optional<ParameterStandard> parameter_standard;
    mutable std::map<std::string, ResidualEnergyAction> energy_actions;
    // An anchored PH/PS cell's g jets at this centre; row_series adds pressure or parameter columns on demand.
    mutable std::shared_ptr<const ReferenceCalorics> calorics;
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
    std::vector<std::vector<double>> standard_minus_plane; // species, Taylor coefficient; requested projections only
    std::array<std::vector<std::vector<double>>, 3> reaction_shift; // L, L_T, L_TT Taylor coefficients
    std::vector<std::vector<double>> solid_amount, solid_affinity;  // [solid][k]: n_s and A_s/RT
    PropertyActionStatus energy = PropertyActionStatus::Available;
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

struct SpinodalBracket final { double lower, upper, lower_eigenvalue, upper_eigenvalue; };
struct SpinodalSearch final {
    std::vector<SpinodalBracket> brackets;
    std::vector<DensityDiagnosticPoint> finite_limits;
    std::vector<double> unresolved_densities;
    double lower = NAN, upper = NAN, seed = NAN, geometric_cap = NAN;
    DerivativeOrder evaluation_order = DerivativeOrder::Hessian;
    std::size_t evaluations = 0;
};

struct SpinodalKey;

struct CompiledSystem final {
    // Internal split context; ordinary declarations retain the inventory-limited optimizer box.
    bool reaction_relaxed_split = false;
    PhaseCalculation phase_calculation = PhaseCalculation::Unknown;
    std::vector<core::NumericalRank> rank_records;
    std::optional<core::NumericalRank> last_rank;
    std::vector<double> optimizer_identity() const;
    const HelmholtzModel* model = nullptr;
    std::shared_ptr<const char> eos_identity;
    const IdealMoleFraction* ideal = nullptr;
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
    std::vector<Layout> shape;  // per phase: its solved-chart columns
    std::vector<std::optional<std::size_t>> charge_row;  // per phase: ordinal of its charge row and multiplier, if any
    std::optional<StandardTerm> fixed_standard;  // the prepared standard term when T and P are both fixed
    struct Declaration;
    std::shared_ptr<const Declaration> original;

    // EqIDs solid_standard_gibbs and solid_saturation_index: g_s/RT and A_s/RT = lambda.B nu_s - g_s/RT.
    double solid_gibbs(const PhasePoint& prepared, std::size_t solid) const;
    double solid_affinity(const PhasePoint& prepared, std::size_t solid) const;

    // State-only preparation differentiates V, support amounts and T when free (phase Hessians);
    // fixed controls and active parameters are padding. Request parameter_actions for fitted row series.
    EPC_API PhasePoint prepare(const std::vector<double>& point, bool state_only = false,
        std::optional<core::Convergence>* first_numerical_failure = nullptr, std::optional<std::size_t> chemical_phase = std::nullopt,
        const std::vector<bool>& site_phases = {}, bool parameter_actions = false) const;
    template <typename Base>
    PhasePointValues<Base> prepare_values(const std::vector<Base>&, bool,
        std::optional<core::Convergence>* = nullptr, std::optional<std::size_t> = std::nullopt, const std::vector<bool>& site_phases = {}, bool parameter_actions = false, bool gradient_only = false, bool rank_only = false) const;
    template <typename Base> PointEvaluationValues<Base> evaluate_at(const PhasePointValues<Base>&, std::optional<std::size_t> = std::nullopt) const;
    // The chemical row of `species` in `phase` at mole fraction `fraction`; at fraction 1, its finite part.
    template <typename Base> Base chemical_row(const PhasePointValues<Base>&, std::size_t phase, std::size_t species, const Base& fraction) const;
    template <typename Base> std::vector<Base> jacobian_at(const PhasePointValues<Base>&, bool constraints_only = false) const;
    template <typename Base> RowUncertaintyValues<Base> row_uncertainty_at(const PhasePointValues<Base>&,
        const std::vector<std::vector<Base>>&) const;
    EPC_API PointEvaluation evaluate(const std::vector<double>& point) const;
    // Residuals that classify a solve: for a free-amount (Gibbs) cell the chemical and
    // pressure rows are scaled by N_p x_i / N_feed and N_p / N_feed, the first-order
    // conditions of the constrained Gibbs minimum in the solved chart; a pinned
    // (feasibility) cell classifies its rows unscaled. Trace moiety rows are relative to their
    // size in both routes (trace_moiety_relative); charged Vapor also retains the log-charge criterion.
    EPC_API std::vector<double> classified_residuals(const PhasePoint& prepared, const std::vector<double>& rows) const;
    // Moiety rows (Jacobian rows for width n) divided by their size, the larger of the feed inventory and the
    // row's term magnitude per mol of feed, when it is at most EQUILIBRIUM_TRACE_MOIETY_FRACTION (#633).
    // Sizes are measured at `sizes` when given, so a polish trial is judged with its step's divisors (#1094).
    template <typename Base> std::vector<Base> trace_moiety_relative(const PhasePointValues<Base>&, std::vector<Base>, std::size_t = 1,
                                                                     const PhasePointValues<Base>* sizes = nullptr) const;
    std::vector<double> trace_moiety_relative(const PhasePoint& prepared, std::vector<double> values,
                                              std::size_t width = 1) const;
    EPC_API std::vector<double> numerical_rows(const PhasePoint&, std::vector<double>, std::size_t width = 1,
                                               const PhasePoint* sizes = nullptr) const;
    EPC_API double classified_residual_norm(const std::vector<double>& point) const;
    EPC_API RowUncertainty row_uncertainty(const PhasePoint&) const;
    EPC_API PointEvaluation evaluate_at(const PhasePoint& prepared, std::optional<std::size_t> chemical_phase = std::nullopt) const;
    EPC_API std::vector<double> generic_point() const;
    // Quantity and basis of each solved-chart coordinate, in coordinate order (unknown_count entries).
    EPC_API std::vector<std::string> coordinate_labels() const;
    // Each jacobian(point) entry against central differences of evaluate() with h_j = relative_step max(1, |x_j|).
    EPC_API std::vector<core::DifferenceCheck> jacobian_difference(const std::vector<double>& point, double relative_step,
        double tau, double kappa) const;
    // First solved-chart coordinate of each declared phase: its log compositions, then a free amount and log packing.
    EPC_API std::vector<std::size_t> phase_offsets() const;
    EPC_API std::vector<double> jacobian(const std::vector<double>& point) const;
    EPC_API std::vector<double> jacobian_at(const PhasePoint& prepared, bool constraints_only = false) const;
    EPC_API std::vector<double> lever_jacobian(const std::vector<double>& point) const;
    EPC_API std::vector<double> objective_gradient(const std::vector<double>& point) const;
    std::vector<double> objective_gradient_at(const PhasePoint& prepared) const;
    EPC_API std::vector<double> hessian(
        const std::vector<double>& point,
        const std::vector<double>& constraint_multipliers = {}
    ) const;
    std::vector<double> hessian_at(
        const PhasePoint& prepared, const std::vector<double>& constraint_multipliers, bool physical = false
    ) const;
    // Multiplier-weighted constraint curvature only; the objective curvature is separate.
    EPC_API std::vector<double> entropy_hessian_at(const PhasePoint&, std::vector<double>* magnitude = nullptr) const;
    EPC_API std::vector<double> constraint_hessian_at(
        const PhasePoint& prepared, const std::vector<double>& constraint_multipliers, bool physical = false
    ) const;
    // Row and physical Taylor coefficients through `order` (at most 3; 2 in a critical cell) along the state
    // coefficients `state` (state[0] is the prepared point) and the control direction.
    EPC_API RowSeries row_series(
        const PhasePoint& prepared, const std::vector<std::vector<double>>& state,
        const std::vector<double>& control, std::size_t order, bool plane = false,
        bool parameter_pairs = false,  // two or more directions carry parameter increments
        const std::vector<std::size_t>& parameter_union = {}
    ) const;

    std::size_t loading_column() const noexcept { return primal_size - free_intensive_count - 1; }
    template<class Scalar> std::vector<Scalar> feed_at(const std::vector<Scalar>& point) const {
        std::vector<Scalar> feed(n_feed.begin(), n_feed.end());
        if (problem.loading) for (std::size_t i = 0; i < feed.size(); ++i)
            feed[i] = problem.loading->initial[i] + point[loading_column()] * problem.loading->direction[i];
        return feed;
    }
    double feed_total_at(const PhasePoint& at) const {
        if (!problem.loading) return feed_total;
        const auto feed = feed_at(at.point);
        return std::accumulate(feed.begin(), feed.end(), 0.0);
    }

    std::size_t component_count() const noexcept { return model ? model->parameters().component_count() : n_feed.size(); }

    template <class Base>
    void pin_liquid_composition(std::vector<Base>& point) const {
        if (!problem.prescribed_liquid_composition) return;
        for (std::size_t i = 0; i < shape[0].x.size(); ++i)
            point[shape[0].x[i]] = Base(std::log((*problem.prescribed_liquid_composition)[problem.phases[0].support[i]]));
    }

    std::size_t unknown_count() const noexcept {
        return ideal ? primal_size : primal_size + rank + charge_row_count;
    }


    bool has_caloric_condition() const noexcept {
        return problem.enthalpy_j.has_value() || problem.entropy_j_per_k.has_value() || problem.internal_energy_j.has_value();
    }

    // PH/PS completes record-less species from reactions or carries a source reference (EqID reference_species_calorics).
    bool anchored_calorics() const;

    // Apparent Cp is defined at driven T and P, or at free T with fixed P and an enthalpy target alone.
    bool admits_apparent_cp() const noexcept {
        const auto& p = problem;
        return !p.pressure.free && !p.pressure.absent
            && ((!p.temperature.free && !p.temperature.absent)
                || (p.temperature.free && p.enthalpy_j && !p.entropy_j_per_k && !p.volume_m3 && !p.internal_energy_j));
    }

    bool feasibility_cell() const noexcept {
        return ideal || pinned_count != 0 || has_caloric_condition() || problem.volume_m3 || (problem.critical.has_value() || problem.spinodal.has_value());
    }
    // One minimized phase: no optimizer ln eta bounds, prepare refuses eta outside the domain (#314).
    bool log_packing_unbounded() const noexcept { return problem.phases.size() == 1 && !feasibility_cell(); }

    std::size_t mode_row() const noexcept { return rows.size() - mode_count - 1 - problem.critical.has_value(); }
    const CriticalDeclaration& mode_seed() const { return problem.critical ? *problem.critical : *problem.spinodal; }
    std::size_t mode_offset() const noexcept { return primal_size - free_intensive_count - mode_count; }

    // Control direction [T, P, n_feed (C), active parameters (A)]; the extended layout appends the
    // (a, b, c, d) increments of every reaction law, then H/S, an admitted volume target and an admitted U target.
    std::size_t control_width() const noexcept { return 2 + n_feed.size() + (model ? model->active_parameters().size() : 0); }
    std::size_t extended_control_width() const noexcept { return control_width() + 4 * problem.reactions.size() + problem.solids.size() + 2 + problem.volume_m3.has_value() + problem.internal_energy_j.has_value(); }

    // A free amount is a coordinate seeded at 1 mol; a pinned amount is its fraction of the feed.
    double phase_amount(const PhaseSpec& phase) const noexcept { return phase.is_free() ? 1.0 : phase.amount.fraction * feed_total; }

    // The declared value, else the default seed of an unseeded free intensive: 300 K and 1e5 Pa.
    double seed_temperature() const noexcept { return problem.temperature.is_positive() ? problem.temperature.value : 300.0; }
    double seed_pressure() const noexcept { return problem.pressure.is_positive() ? problem.pressure.value : 1.0e5; }

    const std::optional<SpinodalSearch>& spinodal_diagnostics() const noexcept { return spinodal_diagnostics_; }

private:
    friend SpinodalKey;
    mutable core::Retained<SpinodalKey, SpinodalSearch> spinodal_selection_;
    mutable std::optional<SpinodalSearch> spinodal_diagnostics_;
    mutable bool spinodal_corrector_ = false;
};

EPC_API CompiledSystem compile_problem(
    const HelmholtzModel& model, const EquilibriumProblem& problem
);
EPC_API CompiledSystem compile_problem(const IdealMoleFraction&, const EquilibriumProblem&);
EPC_API CompiledSystem compile_problem(const HelmholtzModel&, const EquilibriumProblem&, const std::vector<std::string>&);
EPC_API CompiledSystem refresh_problem(const HelmholtzModel&, const EquilibriumProblem&, const std::vector<std::string>&, const CompiledSystem&);

// Width, a nonzero charge-neutral nu, exactly one positive K or law, a neutral reference with P0, and V_s > 0.
void admit_solid_declarations(const HelmholtzModel&, const std::optional<NeutralReferenceDeclaration>&, const std::vector<SolidSpec>&);
// Reads jets.temperature_order, parameters and pairs.
std::vector<SolidLogK> solid_log_k(const HelmholtzModel&, double temperature, double pressure,
    const std::optional<NeutralReferenceDeclaration>&, const std::vector<SolidSpec>&, const ReferenceJets& jets);
// The classified declaration of a solved point: a free solid is present exactly when its amount is positive.

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
EPC_API std::vector<InventoryBound> inventory_bounds(const CompiledSystem& compiled,
    const std::vector<std::size_t>& species, bool fraction, double floor);

// Candidate-topology outcome of the existing producer loop: the first vanished free phase
// in declaration order, else the first coincident pair in lexicographic ordinal order.
enum class Attempt : std::uint8_t { Cold, Warm, Newton, PhysicalAmounts };
enum class Retry : std::uint8_t { None, PhysicalAmounts, Source };
enum class Candidate : std::uint8_t { Absent, Refused, Unstable, Admitted, Published };
enum class TopologyEvent { None, PhaseDisappeared, PhasesCoalesced, BranchEnd };
enum class BranchEndCondition { LoadingFold, VaporMechanicalLimit, VaporKindRefused };
struct BubbleBranchEnd final {
    BranchEndCondition condition;
    std::pair<double, double> feed_bracket_mol;
    std::vector<double> limiting_feed_mol;
    double limiting_pressure_pa;
};

// Spectrum of M, the cubic form and the augmented-Jacobian diagnostics of a critical-equation
// candidate; it certifies no stability.
struct HelmholtzModeRecord {
    std::vector<double> eigenvalues, mode;
    double zero_band = NAN, condition_estimate = NAN, backward_error = NAN;
    std::size_t rank = 0;
};

struct CriticalRecord final : HelmholtzModeRecord { double cubic = NAN; };
struct SpinodalRecord final : HelmholtzModeRecord {
    SpinodalSearch search;
    double normalized_null_residual = NAN;
};

enum class Quantity { Linear, Positive, LogPositive };
enum class ComparisonMode { Strict, EstimatedError };
struct NumericalEstimate final {
    double absolute = NAN, chart = NAN;
    bool available = false, resolved = false;
    std::string cause = "not_evaluated";
};
struct ComparisonEligibility final {
    bool declaration = true, branch = true, phase_matching = true, inventory = true;
    bool residual_a = true, residual_b = true, admission_a = true, admission_b = true;
    bool valid() const { return declaration && branch && phase_matching && inventory && residual_a && residual_b && admission_a && admission_b; }
};
struct EquilibriumComparison final {
    std::string judgment = "failed", cause;
    double difference = 0, chart_difference = 0, radius_sum = 0;
};
template <typename Base>
Base log_difference(Base a, Base b) {
    const Base relative = (a - b) / b;
    return wide_math::abs(wide_math::isfinite(relative) && wide_math::abs(relative) < policy::LOG_RATIO_NEAR_UNITY.value ? wide_math::log1p(relative) : wide_math::log(a) - wide_math::log(b));
}
template <typename Base>
NumericalEstimate project_numerical_estimate(double value, Base mapped, Base radius, bool positive, bool constant) {
    core::Scope scope("project_numerical_estimate");
    const double rounding = constant ? 0.0 : std::numeric_limits<double>::epsilon() * (positive ? 1.0 : wide_math::abs(value));
    const Base displacement = positive ? log_difference(Base(value), mapped) : wide_math::abs(Base(value) - mapped);
    const double chart = static_cast<double>(radius + displacement + rounding);
    const double absolute = positive ? value * std::expm1(chart) : chart;
    const bool available = std::isfinite(absolute) && std::isfinite(chart);
    return {absolute, chart, available, available, available ? "" : "nonfinite_projection"};
}
inline EquilibriumComparison compare_equilibrium_quantity(Quantity quantity, double a, double b, double allowance,
    const NumericalEstimate& ea, const NumericalEstimate& eb, ComparisonEligibility eligibility, ComparisonMode mode) {
    core::Scope scope("compare_equilibrium_quantity");
    EquilibriumComparison out; out.difference = std::abs(a - b);
    if (!std::isfinite(a) || !std::isfinite(b) || !std::isfinite(allowance) || allowance < 0) { out.cause = "nonfinite_operand"; return out; }
    if (quantity == Quantity::Positive && (!(a > 0) || !(b > 0))) { out.cause = "positive_domain"; return out; }
    out.chart_difference = quantity == Quantity::Positive ? log_difference(a, b) : out.difference;
    if (!eligibility.valid()) { out.cause = "prerequisite"; return out; }
    if (out.difference <= allowance) { out.judgment = "ordinary"; out.cause = "ordinary_allowance"; return out; }
    if (mode == ComparisonMode::Strict) { out.cause = "strict_comparison"; return out; }
    const auto usable = [](const NumericalEstimate& e) {
        return e.available && e.resolved && std::isfinite(e.absolute) && std::isfinite(e.chart) && e.absolute >= 0 && e.chart >= 0;
    };
    if (!usable(ea) || !usable(eb)) { out.cause = "unavailable_or_unresolved"; return out; }
    if (ea.absolute <= allowance && eb.absolute <= allowance) { out.cause = "ordinary_achievable"; return out; }
    out.radius_sum = quantity == Quantity::Linear ? ea.absolute + eb.absolute : ea.chart + eb.chart;
    if (!std::isfinite(out.radius_sum)) { out.cause = "nonfinite_sum"; return out; }
    if (out.chart_difference <= out.radius_sum) {
        out.judgment = "estimated_error"; out.cause = "consistent_within_estimated_numerical_error";
    } else out.cause = "outside_estimates";
    return out;
}

struct EquilibriumNumericalEstimates final {
    std::string basis = "canonical_defining_rows", cause = "not_evaluated";
    std::string source_reference_method;
    std::size_t rank = 0, dimension = 0, precision_bits = std::numeric_limits<double>::digits;
    double threshold = 0, condition = INFINITY, inverse_norm = INFINITY, backward_error = INFINITY, physical_radius = INFINITY;
    std::vector<double> residuals, row_terms, row_inner, row_floor, physical_scale;
    equilibrium::NumericalEstimate temperature, pressure, loading;
    std::vector<equilibrium::NumericalEstimate> densities, packing, composition, amounts;
};

template <typename Base>
struct PhysicalProjection final {
    using Vector = Eigen::Matrix<Base, Eigen::Dynamic, 1>;
    Eigen::Matrix<Base, Eigen::Dynamic, Eigen::Dynamic> gradient;
    Vector values, scale;
};
template <typename Base> PhysicalProjection<Base> physical_projection(const CompiledSystem&, const PhasePointValues<Base>&);
// A pressure becomes a coordinate or output scale only when finite and positive (#1008).
template <typename Base> Base pressure_scale(Base pressure) {
    if (!(wide_math::isfinite(pressure) && pressure > 0)) throw std::domain_error("nonpositive_or_nonfinite_pressure");
    return pressure;
}
template <typename Base, typename Vector>
Base physical_radius(const CompiledSystem& c, const PhasePointValues<Base>& at, const Vector& error) {
    Base radius = error[1] / pressure_scale(at.pressure);
    for (Eigen::Index i = 2; i < error.size(); ++i)
        radius = std::max(radius, error[i] / ((i - 2) % (3 + c.component_count()) == 2 ? Base(c.feed_total) : Base(1)));
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
struct CriticalSearch;
struct LinearSolveDiagnostics final {
    std::string system = "one_component_two_density_reduction", reason;
    double residual_norm = std::numeric_limits<double>::infinity(), evaluation_floor = std::numeric_limits<double>::infinity();
    double jacobian_inverse_norm = std::numeric_limits<double>::infinity(), condition_estimate = std::numeric_limits<double>::infinity();
    double forward_error_estimate = std::numeric_limits<double>::infinity();
    std::size_t precision_bits = std::numeric_limits<double>::digits, refinement_iterations = 0;
    std::string point_kind = "last_linearization";
};

struct IpoptIteration final {
    int mode = 0, iteration = 0, line_search_trials = 0;
    double objective = 0, inf_pr = 0, inf_du = 0, mu = 0, step_norm = 0, regularization = 0, alpha_du = 0, alpha_pr = 0;
};
// One Ipopt intermediate callback; `run` counts the Ipopt runs seen by the collecting thread.
struct IpoptIterationRecord final {
    std::size_t run = 0;
    IpoptIteration iteration;
};
// The thread's open DEBUG collector, else null; closed it costs one pointer test per callback (#976).
EPC_API std::vector<IpoptIterationRecord>*& ipopt_iteration_log() noexcept;
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
    Attempt attempt = Attempt::Cold;
    Retry retry = Retry::None;
    Candidate candidate = Candidate::Absent;
    std::shared_ptr<const CriticalSearch> critical_search;
    std::shared_ptr<const EquilibriumResult> rejected_reduction;
    std::optional<core::Convergence> terminal_numerical_failure, last_full_evaluation, limit_record;
    std::vector<double> phase_pressures;
    std::optional<core::Convergence> first_numerical_failure;
    std::vector<core::NumericalRank> rank_records;
    std::optional<IpoptDiagnostics> ipopt_diagnostics, recovered_ipopt_diagnostics;
    std::optional<core::NumericalRank> last_rank;
    std::map<std::string, double> solver_settings;
    std::string diagnostic_options_file, linear_solver;  // linear_solver: empty when no Ipopt solve ran
    std::shared_ptr<const OptimizerState> optimizer_state;
    std::shared_ptr<const char> model_identity, ideal_identity;
    std::shared_ptr<const EquilibriumProblem> declaration;
    std::shared_ptr<DiscoveryRecord> discovery;
    std::shared_ptr<const StabilityReport> recovery_stability;  // #987: the trial search run before recovery
    std::shared_ptr<const StabilityReport> vanished_phase_stability;  // a vanished phase's trials, without it
    std::vector<double> vanished_composition;
    double vanished_amount = 0.0;
    bool success = false;
    std::string message;
    std::string solver_status;
    std::optional<double> temperature;
    std::optional<double> pressure;
    std::vector<double> amounts;
    std::optional<double> electric_potential_difference_v;  // second liquid minus first
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
    // ALG-TERMINATION: the species a default cold trial attempt left at their composition bounds, rows unconverged.
    std::vector<std::size_t> bound_species;
    struct HeldSpecies { std::size_t phase, species; double fraction, chemical_row, complementarity; std::optional<double> response{}; };
    std::vector<HeldSpecies> held_species;
    std::size_t unknown_count = 0;
    std::size_t jacobian_rank = 0;
    std::size_t generic_nullity = 0;
    // Absent without a finite final row residual; otherwise residual <= caller row_tolerance.
    std::optional<bool> requested_tolerance_met;
    TopologyEvent topology_event = TopologyEvent::None;
    std::optional<std::size_t> topology_phase;
    std::optional<std::pair<std::size_t, std::size_t>> topology_pair;
    std::vector<NeutralReferenceReactionDiagnostic> reference_diagnostics;
    std::optional<BubbleBranchEnd> branch_end;
    std::optional<CriticalRecord> critical;
    std::optional<SpinodalRecord> spinodal;
    std::vector<double> solid_amounts, solid_affinity_decades;  // per declared solid: mol, and SI_s in decades
    // A numerical failure after the driver began, with its work record; solve_compiled rethrows it, refuse clears it.
    std::exception_ptr failure;

    bool has_state() const noexcept {
        return success && temperature && pressure
            && !mole_fractions.empty() && (!amounts.empty()
                || (declaration && declaration->prescribed_liquid_composition
                    && molar_densities.size() == 2
                    && mole_fractions.size() == 2 * declaration->prescribed_liquid_composition->size()));
    }

    double residual_norm() const noexcept { return maximum_abs(residuals); }
};

CompiledSystem classify_result(CompiledSystem compiled, const EquilibriumResult& result);

// Ordered critical starts and every returned candidate; finite search certifies no stability or uniqueness.
struct CriticalSearch final {
    std::vector<std::array<double, 2>> seeds;  // attempted T [K], rho [mol/m^3]; original declaration first
    std::vector<EquilibriumResult> results;
    std::vector<std::array<std::size_t, 3>> added_work;  // cumulative iterations, evaluations, Hessian callbacks
    std::size_t scheduled = 0;  // generated seeds, excluding the original declaration
    bool budget_exhausted = false;
    std::optional<std::size_t> selected;
};

// `counted`: the estimate's roots and energy evaluations count in the enclosing calculation work.
EPC_API EquilibriumNumericalEstimates equilibrium_numerical_estimates(const CompiledSystem&, const EquilibriumResult&, const PhasePoint* = nullptr, bool counted = false, bool plane_rows = false,
    const std::vector<double>& prescribed_errors = {}, bool replay = false);
struct DriverLimits;
EPC_API EquilibriumNumericalEstimates equilibrium_numerical_estimates(const CompiledSystem&, const EquilibriumResult&, const std::vector<WiderScalar>&, const DriverLimits* = nullptr);
EPC_API EquilibriumNumericalEstimates equilibrium_numerical_estimates(const HelmholtzModel&, const EquilibriumResult&);
EPC_API EquilibriumNumericalEstimates equilibrium_numerical_estimates(const IdealMoleFraction&, const EquilibriumResult&);

// EqID reaction_activity_quotient; per declared phase/reaction, in component order.
struct ReactionActivities final {
    std::size_t phase_index = 0, reaction_index = 0;
    std::string basis, cause;
    std::optional<double> reference_pressure_pa, log_quotient;
    // Pressure basis (EOS models): P0, ln Q_P = ln prod (f_i/P0)^nu_i, and its partial actions along T, P and e_j - x (amount j).
    std::optional<double> standard_pressure_pa, log_pressure_quotient;
    std::vector<std::optional<double>> log_pressure_quotient_actions;
    std::string pressure_cause;
    SpeciesActivities species;
};
EPC_API std::vector<ReactionActivities> reaction_activities(const HelmholtzModel&, const EquilibriumResult&);
EPC_API std::vector<ReactionActivities> reaction_activities(const IdealMoleFraction&, const EquilibriumResult&);
// Per declared solid: 1 when present, absent otherwise.
EPC_API std::vector<std::optional<double>> solid_activities(const EquilibriumResult&);
// Finite non-negative amounts and compositions, positive densities and normalized phases, or throws.
void require_solution_shape(const EquilibriumResult& result, std::size_t phases, std::size_t components);

struct StabilityPhaseComparison final {
    std::size_t phase = 0;
    std::vector<std::size_t> added_species, missing_species;
    double composition_max_distance = 0.0, volume_log_distance = 0.0;
};
struct StabilityTrial final {
    PhaseSpec phase;
    std::string seed, evidence, message;
    EquilibriumResult optimizer;
    std::vector<double> raw_coordinates, raw_residuals, composition, split_coordinates;
    std::vector<double> reduced_split_coordinates;  // the split without `optimizer.bound_species` in the added phase
    std::vector<std::vector<double>> basis, row_transformation;
    // Quantities carry units in their names; errors are local estimates, never enclosures.
    std::map<std::string, double> quantities;
    std::vector<std::map<std::string, double>> split_checks;
    std::vector<StabilityPhaseComparison> present_phase_comparisons;
    std::string diagnostic_cause;
    std::optional<std::size_t> donor;
    bool resolved_descent = false;
    bool pivot_only = false;
};
enum class StabilityDescent { ResolvedDescent, NoWitness, Undecided };
struct StabilityReport final {
    EPC_API std::map<std::string, std::size_t> evidence_counts(const std::vector<std::size_t>& multiplicities = {}) const;
    EPC_API std::map<std::string, std::size_t> unresolved_reasons() const;
    EPC_API StabilityDescent descent() const;
    EPC_API const char* descent_name() const;
    EPC_API std::vector<std::size_t> undecided() const;
    EPC_API std::vector<StabilityTrial> witnesses() const;
    EPC_API std::vector<StabilityTrial> supplied_minima() const;
    std::optional<core::Convergence> first_numerical_failure;
    std::optional<core::NumericalRank> last_rank;
    bool completed = false;
    bool no_descent_in_trial_set = false;  // #987: completed, trials, none descending or unresolved
    std::string message;
    std::vector<StabilityTrial> trials;
    std::vector<double> candidate_coordinates, standard_over_rt, plane_over_rt;
    std::vector<std::map<std::string, double>> present_phases;
    std::map<std::string, double> work;
};
struct StabilityGridCoverage final {
    PhaseSpec phase;
    std::vector<std::size_t> omitted_species, nodes;
    std::size_t planned = 0, visited = 0, admitted = 0;
    std::optional<double> rounding_bound, h, mapping_allowance;
    std::optional<std::size_t> worst_gap, nearest_admitted;
    std::string cause;
    bool uninformative = false;
};
struct StabilityGridResult final {
    std::optional<std::uint64_t> raw_attempts;
    std::uint64_t ceiling;
    std::vector<StabilityGridCoverage> coverage;
    StabilityReport raw, refined;
    double elapsed_seconds = 0.0;
};
EPC_API StabilityGridResult stability_grid(const CompiledSystem&, const EquilibriumResult&, const std::vector<PhaseSpec>&,
    std::uint64_t denominator, std::uint64_t ceiling = 2048, bool refine = false, std::optional<std::size_t> trial_budget = {});

struct StabilityMeasure final {
    StabilityReport report;
    StabilityOptions options;
    std::optional<double> raw_minimum, value, uncertainty;
    std::optional<std::size_t> minimizing_trial;
    std::vector<std::size_t> minima;
    std::string cause, event;
    bool available = false, coincident_supports = false;
};
EPC_API std::pair<StabilityReport, bool> supplied_trial_measure(const CompiledSystem&, const EquilibriumResult&,
    const std::shared_ptr<const PreparedSuppliedTrials>&, const std::vector<double>& standard_uncertainty_over_rt);
EPC_API StabilityMeasure stability_measure(const CompiledSystem&, const EquilibriumResult&, const StabilityOptions& = {});

struct DiscoveryRecord final {
    DiscoveryMode mode = DiscoveryMode::Full;
    std::string termination;
    std::vector<std::string> unsearched;
    EquilibriumProblem declaration;
    EquilibriumResult candidate;
    std::vector<StabilityReport> rounds;
    std::vector<EquilibriumResult> declarations;
    std::vector<std::string> declaration_roles;
    std::vector<std::vector<double>> declaration_attempt_coordinates;
    std::map<std::string, double> work;
};
EPC_API void admit_trial_feed(CompiledSystem&);
EPC_API void normalize_trial_feed(CompiledSystem&);
EPC_API void admit_inventory_bound(const InventoryBound&);
EPC_API std::shared_ptr<const PreparedSuppliedTrials> prepare_supplied_trials(const HelmholtzModel&, double, const std::vector<PhaseSpec>&);
EPC_API CompiledSystem compile_stability_trial(const HelmholtzModel&, double, double, const PhaseSpec&,
    const std::vector<double>&, const std::vector<double>&, std::vector<std::vector<double>>&);
EPC_API double present_phase_residual_distance(double, double);
EPC_API bool descent_resolved(double, double, double);
EPC_API std::vector<double> reconstruct_trial_composition(std::vector<double>, const std::vector<std::size_t>&, const std::vector<int>&, PhaseKind = PhaseKind::Liquid, bool* pivot_only = nullptr, bool keep_charge = false);
// n dln(phi_i)/dn_j at fixed T, P on the support (row-major), from an IsothermalPhaseHessian evaluation at unit amount and v.
EPC_API std::vector<double> tangent_plane_jacobian(const HelmholtzDerivatives&, double, const std::vector<std::size_t>&);
EPC_API std::pair<EquilibriumProblem, std::vector<double>> remove_discovery_phase(const CompiledSystem&,
    const std::vector<double>&, std::map<std::string, double>* = nullptr);
EPC_API StabilityReport stability_test(const CompiledSystem&, const EquilibriumResult&, const StabilityOptions& = {});
EPC_API std::pair<EquilibriumResult, std::vector<StabilityReport>> cooling_onset(const CompiledSystem&);
EPC_API EquilibriumResult discover_phases(const HelmholtzModel&, const EquilibriumProblem&);

// Set only for an attempt whose declared failure is reported directly to a named caller.
struct EquilibriumAttemptReporting final {
    const std::vector<std::string>* composition_bound_species = nullptr;
};
EPC_API EquilibriumResult solve_equilibrium(
    const HelmholtzModel& model,
    const EquilibriumProblem& problem,
    const std::optional<EquilibriumResult>& start = std::nullopt, const std::vector<std::string>& component_ids = {},
    std::optional<std::size_t> work_limit = std::nullopt, std::optional<std::size_t> branch_end_work_limit = std::nullopt
);
EquilibriumResult solve_equilibrium_retaining_compiled(const HelmholtzModel&, const EquilibriumProblem&,
    std::optional<CompiledSystem>&, const std::optional<EquilibriumResult>& = std::nullopt,
    const std::vector<std::string>& component_ids = {}, const CompiledSystem* prepared = nullptr);
EPC_API EquilibriumResult solve_equilibrium(const HelmholtzModel&, const EquilibriumProblem&,
    const std::optional<EquilibriumResult>&, const std::vector<std::string>&, const CompiledSystem& prepared);
EPC_API EquilibriumResult solve_equilibrium(const IdealMoleFraction&, const EquilibriumProblem&,
    const std::optional<EquilibriumResult>& start = std::nullopt);
// The same solve of a compiled problem (its model must outlive it), inside the caller's open span if given.
// `recover`: the public no-start solve and cooling onset only (ALG-EQUILIBRIUM neutral recovery, #565).
// try_solve_compiled returns a numerical failure in its result; solve_compiled throws it.
struct DriverLimits;
EPC_API EquilibriumResult try_solve_compiled(const CompiledSystem& compiled, const std::optional<EquilibriumResult>& start = std::nullopt,
                                 core::Scope* outer = nullptr, const DriverLimits* limits = nullptr,
                                 std::shared_ptr<const OptimizerState> warm = {}, bool recover = false, bool select_endpoint = false, EquilibriumAttemptReporting reporting = {});
EPC_API EquilibriumResult solve_compiled(const CompiledSystem& compiled, const std::optional<EquilibriumResult>& start = std::nullopt,
                                 core::Scope* outer = nullptr, const DriverLimits* limits = nullptr,
                                 std::shared_ptr<const OptimizerState> warm = {}, bool recover = false, bool select_endpoint = false, EquilibriumAttemptReporting reporting = {});

EPC_API EquilibriumResult solve_saturation(const CompiledSystem&);
bool one_component_two_density_eligible(const CompiledSystem&);
struct SensitivityInfo final {
    std::size_t rank = 0;
    double condition_estimate = NAN, threshold = 0.0;
    bool nonsingular = false;
};

EPC_API SensitivityInfo inspect_sensitivity(
    const CompiledSystem& compiled,
    const std::vector<double>& point,
    bool estimate_condition = false
);
struct SolvedStateActionDiagnostics;
// All operands are generated from one exact trial declaration in declared phase order.
struct SaturationOperand final {
    equilibrium::Quantity quantity = equilibrium::Quantity::Linear;
    double value = NAN;
    equilibrium::NumericalEstimate error;
};
struct SaturationCandidate final {
    EquilibriumResult state;
    std::vector<SaturationOperand> operands;
    bool admitted = false;
};
struct SaturationComparison final {
    std::size_t left = 0, right = 0;
    std::vector<equilibrium::EquilibriumComparison> operands;
    std::string identity;
};
struct SaturationNode final {
    double loading = 0, mismatch = NAN, slope = NAN;
    EquilibriumResult state;
    EquilibriumNumericalEstimates estimates;
    std::string action_cause, prior_start;
    std::size_t driver_attempts = 0;
    double action_seconds = 0, estimate_seconds = 0;
    EquilibriumNumericalEstimates local_estimates;
    std::shared_ptr<SolvedStateActionDiagnostics> action_diagnostics;
};
enum class RefusedOperation : std::uint8_t {
    declaration, central_solve, source_recovery, loading_search, search_node_solve, search_node_action, search_node_estimate,
    search_candidate_solve, search_candidate_estimate, search_candidate_admission, search_root_selection,
    solved_state_action, thermochemistry, density_closure, observation_form
};
struct SaturationSearch final {
    std::optional<double> loading_cap;
    std::optional<core::Convergence> refusal_convergence;
    std::optional<RefusedOperation> refused_operation;
    std::vector<double> grid;
    std::vector<SaturationNode> nodes, rejected_starts;
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
void compare_loading_candidates(SaturationSearch&, equilibrium::ComparisonEligibility);
EPC_API void admit_loading_search(SaturationSearch&, std::optional<double> loading_cap = {});
EPC_API SaturationSearch search_loading(const HelmholtzModel&, EquilibriumProblem, const AffineLoading&,
    const std::vector<double>&, const DriverLimits* = nullptr,
    std::chrono::steady_clock::time_point = std::chrono::steady_clock::time_point::max(), std::optional<double> loading_cap = {},
    const SaturationSearch* earlier = nullptr);

// ---- Solved-state derivative actions: native records, no public promise.

enum class SolvedStateObservableKind {
    Temperature,
    Pressure,
    SpeciesStandardMinusPlane,
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
    Loading,
    PhaseBulkRelativePermittivity,
    AssociationSiteFraction
};

struct SolvedStateObservable final {
    SolvedStateObservableKind kind = SolvedStateObservableKind::Temperature;
    std::size_t phase = 0;
    std::size_t component = 0;
    PropertyObservable property = PropertyObservable::IdealEnthalpy;
    std::optional<std::size_t> site{};
    bool plane_uncertainty = false;
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
    double volume_m3 = 0.0;
    double internal_energy_j = 0.0;
    std::vector<double> solid_log_k;  // per declared solid: increment of its source ln K at the fixed T, or empty
};

EPC_API void validate_phase_site(const Parameters&, const PhaseSpec&, std::optional<std::size_t>);
EPC_API HelmholtzDerivatives phase_site_values(const HelmholtzModel&, const EquilibriumResult&, std::size_t phase);

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

// A source or solid reference, or a reference column a motion reads, is unavailable; the one producer of the label.
struct EPC_EXCEPTION ReferenceUnavailable final : std::invalid_argument {
    static std::string label(const std::string& cause) { return "reference_unavailable: " + cause; }
    explicit ReferenceUnavailable(const std::string& cause) : std::invalid_argument(label(cause)) {}
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
        std::optional<std::size_t> parameter, second_parameter;
        Column column; // raw derivatives of the record-anchored density-standard g, J/mol
    };
    std::vector<Derivative> gibbs;
    std::vector<std::vector<double>> reaction_map; // anchored inverse of nu_D by reaction
    std::vector<std::vector<double>> reaction_defect; // I - nu_D nu_D^+, admits combined source motions
    bool thermal_boundary = false;
    SolvedStateActionStatus status = SolvedStateActionStatus::Available;
    std::string message;
};
// T columns through max(jets.temperature_order, 1); reads pressure, parameters, pairs and errors.
EPC_API ReferenceCalorics reference_calorics(const CompiledSystem &compiled, double temperature, double pressure,
                                             const ReferenceJets &jets);
// Per unit x of one phase: ideal h = sum x(g - T g_T + RT), s = sum x(R - g_T - R ln(x/v rho0)) and cp = sum x(R - T g_TT)
// from the Taylor polynomial of the anchored g at its centre in (dT, dP, dtheta), exact through the computed order.
// Unavailable columns are skipped; their consumers refuse the motions that read them.
template <typename Scalar>
IdealMixtureValue<Scalar> anchored_ideal(const ReferenceCalorics& c, const std::vector<std::size_t>& support,
    const std::vector<Scalar>& x, const Scalar& volume, const Scalar& T, const Scalar& dT,
    const Scalar& dP = Scalar(0.0), const std::vector<Scalar>& dtheta = {}) {
    const Scalar R = Scalar(GAS_CONSTANT_J_PER_MOL_K);
    IdealMixtureValue<Scalar> out;
    for (const auto i : support) {
        std::array<Scalar, 3> g{Scalar(0.0), Scalar(0.0), Scalar(0.0)};  // g, g_T, g_TT
        for (const auto& e : c.gibbs) {
            if (e.column.status != SolvedStateActionStatus::Available || e.second_parameter || (e.parameter && *e.parameter >= dtheta.size())) continue;
            Scalar w = Scalar(e.column.values[i]);
            for (std::size_t k = 1; k <= e.pressure_order; ++k) w *= dP / Scalar(double(k));
            if (e.parameter) w *= dtheta[*e.parameter];
            for (std::size_t m = 0; m < 3 && m <= e.temperature_order; ++m) {
                Scalar term = w;
                for (std::size_t k = 1; k <= e.temperature_order - m; ++k) term *= dT / Scalar(double(k));
                g[m] += term;
            }
        }
        out.enthalpy += x[i] * (g[0] - T * g[1] + R * T);
        out.entropy += x[i] * (R - g[1] - R * scalar_log(x[i] / volume));
        out.isobaric_heat_capacity += x[i] * (R - T * g[2]);
    }
    return out;
}

struct SolvedStateActionResult final {
    SolvedStateObservable observable;
    std::string unit;
    std::size_t order = 0;
    std::optional<double> value;
    std::optional<double> action;
    std::optional<double> numerical_error;
    SolvedStateActionStatus status = SolvedStateActionStatus::CentralSolveUnavailable;
    std::optional<PropertyActionDiagnostics> property_diagnostics;
    std::optional<std::string> diagnostic_message;
    std::optional<core::Convergence> refusal;
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

EPC_API std::string solved_state_observable_id(const SolvedStateObservable& observable);
// Each observable's first-order action along `direction` against central differences of re-solved neighbors at
// +-step (verification only). Supported: fixed T or P, Amounts-basis feed, declared enthalpy or entropy.
EPC_API std::vector<core::DifferenceCheck> solved_state_difference(const CompiledSystem& compiled,
    const EquilibriumResult& central, const std::vector<SolvedStateObservable>& observables,
    const SolvedStateActionDirection& direction, double step, double tau, double kappa);

// Exact first-, second- or third-order mixed directional actions of physical solved-state
// observables at one accepted result; no nonlinear solve, no compilation, no mutation.
EPC_API SolvedStateActionBatch solved_state_actions(
    const CompiledSystem& compiled,
    const EquilibriumResult& central,
    const SolvedStateActionRequest& request
);
// The ASL reader resets this workspace before replacing its compiled system or result.
struct SolvedStateActionWorkspace;
EPC_API SolvedStateActionBatch solved_state_actions_prepared(const CompiledSystem& compiled, const EquilibriumResult& central,
    const SolvedStateActionRequest& request, std::shared_ptr<SolvedStateActionWorkspace>& workspace);
// Independent first-order directions at one accepted center; each response retains its own status.
struct SolvedStateCaloricCenter;
EPC_API std::vector<SolvedStateActionBatch> solved_state_first_actions(
    const CompiledSystem& compiled, const EquilibriumResult& central,
    const std::vector<SolvedStateObservable>& observables,
    const std::vector<SolvedStateActionDirection>& directions, std::shared_ptr<const SolvedStateCaloricCenter>* caloric = nullptr
);

}  // namespace epcsaft::equilibrium
