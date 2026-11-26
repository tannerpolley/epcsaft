#pragma once
#include "core/export.hpp"

#include "equilibrium/ipopt.hpp"
#include "equilibrium/problem.hpp"

#include <IpTNLP.hpp>
#include <Eigen/Dense>

#include <cstddef>
#include <chrono>
#include <functional>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace epcsaft::core { class Scope; }


namespace epcsaft::equilibrium {

struct HelmholtzSpectrum final {
    Eigen::MatrixXd matrix;
    Eigen::VectorXd values, mode;
    double band = NAN;
};
EPC_API HelmholtzSpectrum helmholtz_spectrum(const HelmholtzDerivatives&, const std::vector<double>&,
    const std::vector<std::size_t>&, const std::vector<double>& normalization = {});
EPC_API SpinodalSearch spinodal_search(const CompiledSystem&, const DriverLimits* = nullptr, double anchor = NAN, const SpinodalSearch* carried = nullptr, SpinodalSearch* attempted = nullptr);
EPC_API CompiledSystem spinodal_corrector(const CompiledSystem&, const DriverLimits*, double anchor, const SpinodalSearch* carried, SpinodalSearch& attempted);
EPC_API EquilibriumResult spinodal_search_failure(const CompiledSystem&, EquilibriumResult, const std::string&, const SpinodalSearch* attempted = nullptr);
// Returned work is accumulated unless the retry caller requests invocation work.
enum class Work : std::uint8_t { Accumulated, Invocation };
struct StageWork final {  // the work of the sub-solves that seed a solve, failed ones included
    std::size_t iterations = 0, evaluations = 0, hessian_calls = 0, stage = 0;
    std::optional<core::Convergence> first_numerical_failure;
    std::vector<core::NumericalRank> rank_records;
    std::optional<IpoptDiagnostics> ipopt_diagnostics, recovered_ipopt_diagnostics;
    template<class R> void observe(const R& r) {
        core::retain_first(first_numerical_failure, r.first_numerical_failure);
        for (const auto& rank : r.rank_records) core::retain_rank(rank_records, rank, {}, stage);
        if constexpr (std::is_same_v<R, EquilibriumResult>)
            if (r.last_rank) core::retain_rank(rank_records, *r.last_rank, {}, stage);
        ++stage;
        if (!recovered_ipopt_diagnostics && ipopt_diagnostics && ipopt_diagnostics->first_restoration_iteration && r.ipopt_diagnostics)
            recovered_ipopt_diagnostics = ipopt_diagnostics;
        if (!recovered_ipopt_diagnostics) recovered_ipopt_diagnostics = r.recovered_ipopt_diagnostics;
        if (r.ipopt_diagnostics) ipopt_diagnostics = r.ipopt_diagnostics;
    }
    void publish_observations(EquilibriumResult& r) && {
        r.first_numerical_failure = std::move(first_numerical_failure);
        r.rank_records = std::move(rank_records);
        r.ipopt_diagnostics = std::move(ipopt_diagnostics);
        r.recovered_ipopt_diagnostics = std::move(recovered_ipopt_diagnostics);
    }
    static void publish_observations(const EquilibriumResult& source, EquilibriumResult& destination) {
        StageWork observations; observations.observe(source); std::move(observations).publish_observations(destination);
        if (!destination.declaration) destination.declaration = source.declaration;
        if (source.last_rank) destination.last_rank = source.last_rank;
    }
    void report_observations(EquilibriumResult& r) const {
        StageWork combined = *this; combined.observe(r); std::move(combined).publish_observations(r);
    }
    void add(const EquilibriumResult& r) { iterations += r.iterations; evaluations += r.evaluations; hessian_calls += r.hessian_calls; observe(r); }
    void report(EquilibriumResult& r) && {
        r.iterations += iterations;
        r.evaluations += evaluations;
        r.hessian_calls += hessian_calls;
        observe(r);
        std::move(*this).publish_observations(r);
    }
};
// Hands the exception a failed result carries to `handle`, for handlers that read its type.
template<class Handle> void with_failure(const EquilibriumResult& failed, Handle&& handle) {
    try { std::rethrow_exception(failed.failure); } catch (const std::exception& error) { handle(error); }
}
// Shared general continuation; an original start skips initialization and reduced attempts.
EquilibriumResult continue_compiled(const CompiledSystem&, const std::optional<EquilibriumResult>&,
    core::Scope&, const DriverLimits*, std::shared_ptr<const OptimizerState>, bool recover,
    const EquilibriumResult& structural, StageWork stages = {}, std::optional<std::vector<double>> original = std::nullopt,
    std::shared_ptr<const EquilibriumResult> rejected = {}, bool select_endpoint = false, EquilibriumAttemptReporting reporting = {});
// `speciate` false: a supplied composition replaces the start's, so a reacting phase is not speciated.
// A `packed` phase's packing guess is supplied by the caller and replaces its density root.
EPC_API std::vector<double> boundary_start(const CompiledSystem& compiled, bool speciate = true, const DriverLimits* limits = nullptr,
                                   StageWork* stages = nullptr, const std::vector<bool>& packed = {}, bool correct_dense = false);
struct OneComponentTwoDensityPoint final {
    PhasePoint phase;
    std::vector<double> rows, jacobian;
    double floor;
};
EPC_API OneComponentTwoDensityPoint one_component_two_density_point(const CompiledSystem&, const std::vector<double>&, std::optional<core::Convergence>* = nullptr, const PhasePoint* = nullptr, std::optional<double> log_coefficient = std::nullopt);
EPC_API EquilibriumResult one_component_two_density_attempt(const CompiledSystem&, const std::vector<double>&,
    const EquilibriumResult&, core::Scope&, const DriverLimits*);
EPC_API EquilibriumResult publish_candidate(const CompiledSystem&, const PhasePoint&, EquilibriumResult,
    core::Scope&, const DriverLimits* = nullptr);
bool same_result_origin(const CompiledSystem&, const EquilibriumResult&);
// An empty result carrying the model identity and declaration that same_result_origin admits.
EPC_API EquilibriumResult result_origin(const CompiledSystem&);
bool newton_eligible(const CompiledSystem&), recovery_eligible(const CompiledSystem&);
bool discovery_newton_eligible(const CompiledSystem&);
EPC_API EquilibriumResult newton_attempt(const CompiledSystem&, const std::vector<double>& initial,
    const EquilibriumResult& structural, const DriverLimits* limits = nullptr, std::vector<double>* attempt = nullptr);
bool supports_split_seed(const CompiledSystem& compiled);
std::string declaration_difference(const EquilibriumProblem& a, const EquilibriumProblem& b, bool endpoint = true, bool structural = false);
// The same declared thermodynamic state: declaration_difference without seeds or labels, plus feed values
// and fixed T and P. Empty when equal; shared by previous-result starts and stability tests.
std::string declared_state_difference(const EquilibriumProblem& a, const EquilibriumProblem& b, bool structural = false);

// The producer's candidate-topology predicate on full coordinates: first vanished free phase
// in declaration order, else the first coincident pair; shared by the driver and continuation.
// It owns the chart reading: normalized phase compositions, phase amounts and ln eta.
struct TopologyCheck final {
    std::vector<std::vector<double>> compositions;
    std::vector<double> amounts, log_packing;
    TopologyEvent event = TopologyEvent::None;
    std::optional<std::size_t> phase;
    std::optional<std::pair<std::size_t, std::size_t>> pair;
    std::vector<double> distances;
    double distance = 0.0;
};

EPC_API TopologyCheck check_topology(const CompiledSystem& compiled, const std::vector<double>& coordinates);

// One density-branch admission owner for public solves and anchored continuation.
struct BranchCheck final {
    bool accepted = false;
    bool terminal = false;
    std::string message;
};

EPC_API BranchCheck admit_null_mode(const CompiledSystem&, EquilibriumResult&);
BranchCheck check_density_branches(
    const CompiledSystem& compiled, const EquilibriumResult& result,
    const std::vector<double>& anchors, bool initial, core::Scope& scope, const DriverLimits* limits = nullptr
);
// Cold root on the declared branch; throws when the stored negative-slope sample places a sole root on the other side (#423).
EPC_API DensityRoot declared_root(const HelmholtzModel& model, double T, double P, const std::vector<double>& x, PhaseKind kind, bool* same_side = nullptr);

struct PhaseAmountStepError final : std::logic_error {
    using std::logic_error::logic_error;
};

struct EPC_EXCEPTION BudgetExhausted final : std::runtime_error {
    using std::runtime_error::runtime_error;
};

struct DriverLimits final {
    std::size_t iterations;
    std::chrono::steady_clock::time_point deadline;
    mutable std::size_t consumed = 0;  // summed over the solve and its seed sub-solves
    // Recovery work (#565) also caps objective/polish evaluations and Hessian callbacks; a spent cap is final.
    bool recovery = false;
    std::size_t evaluation_cap = 0, hessian_cap = 0;
    mutable std::size_t evaluations = 0, hessian_calls = 0;
    mutable bool exhausted = false;
    std::size_t attempt_cap = 0;
    mutable std::size_t attempts = 0;
    // #21: a guarded pass refuses an optimizer or polish point for which this returns a cause.
    std::function<std::optional<std::string>(const PhasePoint&)> admit;
    static DriverLimits recovery_allowance(std::chrono::steady_clock::time_point deadline, std::size_t iterations) {
        return {iterations, deadline, 0, true, static_cast<std::size_t>(policy::EQUILIBRIUM_RECOVERY_EVALUATIONS.value),
            static_cast<std::size_t>(policy::EQUILIBRIUM_RECOVERY_HESSIAN_CALLS.value)};
    }
    using Limit = core::WorkLimit;
    mutable Limit refused = Limit::none;
    std::size_t attempt_iterations = policy::EQUILIBRIUM_MAX_ITERATIONS.value;

    bool deny(Limit gate) const {
        if (refused == Limit::none) refused = gate;
        return true;
    }
    bool iterations_spent() const { return consumed >= iterations && deny(Limit::iterations); }
    void retain_denial(core::Convergence& record) const { if (record.limit == Limit::none) record.limit = refused; }
    bool caps_spent() const { return (evaluations >= evaluation_cap && deny(Limit::evaluations)) || (hessian_calls >= hessian_cap && deny(Limit::hessian_calls)); }
    bool expired() const {
        exhausted = exhausted || (recovery && iterations_spent()) || (attempt_cap && attempts > attempt_cap && deny(Limit::attempts));
        return exhausted || (std::chrono::steady_clock::now() >= deadline && deny(Limit::deadline));
    }
    void check(const char* cause = "elapsed work budget exhausted") const { if (expired()) throw BudgetExhausted(cause); }
    // Before one evaluation or Hessian callback: false, refusing it, once any recovery cap is spent.
    bool spend(std::size_t& used, std::size_t cap) const {
        exhausted = exhausted || (recovery && ((used >= cap && deny(&used == &evaluations ? Limit::evaluations
            : &used == &hessian_calls ? Limit::hessian_calls : Limit::iterations)) || iterations_spent()));
        return !exhausted && ++used;
    }
};
bool supports_work_limit(const EquilibriumProblem&);
EquilibriumResult solve_from_feed_increment(const CompiledSystem&, const EquilibriumResult&, core::Scope&, std::optional<std::size_t>);
// ALG-EQUILIBRIUM source recovery: the finite pure-temperature schedule or mixture source/continuation;
// all child calculations spend `limits`. A refusal has no state and retains the recovery cause.
EPC_API EquilibriumResult recover_from_source(const CompiledSystem& target, const DriverLimits& limits);
EPC_API EquilibriumResult recover_from_source(const CompiledSystem&, EquilibriumResult,
    const DriverLimits* = nullptr, std::optional<core::Convergence>* = nullptr);
EPC_API std::string start_difference(const std::optional<EquilibriumResult>&, const CompiledSystem&);
EPC_API bool retain_limit(std::optional<core::Convergence>&, const core::Convergence&);
inline const core::Convergence& solve_refusal(const EquilibriumResult& result) {
    if (result.limit_record) return *result.limit_record;
    return result.terminal_numerical_failure ? *result.terminal_numerical_failure : result.convergence;
}
// #987: the finite trial search of a coalesced fixed-T/P flash before recovery, charged to `limits`.
EPC_API StabilityReport recovery_stability(const CompiledSystem& compiled, const DriverLimits& limits);
// The public explicit test with `recovery` funding its second pass over condensed-domain refusals (#21).
EPC_API StabilityReport recovered_stability_test(const CompiledSystem&, const EquilibriumResult&, const StabilityOptions&, DriverLimits recovery);
bool reservoir_applicable(const CompiledSystem&, std::size_t phase);
// ALG-STABILITY: one trial of the vanished phase's kind and support from the declaration without it, charged to `limits`.
StabilityReport vanished_phase_stability(const CompiledSystem&, std::size_t phase, const DriverLimits&,
    const std::vector<double>& first, EquilibriumResult& result);
EPC_API EquilibriumResult refine_selected_endpoint(const CompiledSystem&, EquilibriumResult, core::Scope&, const DriverLimits* = nullptr);
enum class GibbsDifferenceStatus { Decrease, NoDecrease, Overlap, Unavailable };
struct GibbsDifference final {
  GibbsDifferenceStatus status = GibbsDifferenceStatus::Unavailable;
  std::map<std::string, double> quantities;
  std::vector<double> coordinates; // repaired split, candidate-coordinate layout
  std::string cause;
};
EPC_API bool capped_gibbs_descent(double integral, double integral_error, double endpoint, double endpoint_error);
enum class CappedGibbsDescent { StandardsDiffer, Unresolved, Resolved };
EPC_API CappedGibbsDescent capped_gibbs_descent(double integral, double integral_error, double endpoint, double endpoint_error,
    const std::vector<double>& split_standard, const std::vector<double>& parent_standard);
// `donor` indexes the parent; the split's donor is its first phase, the trial its last.
EPC_API GibbsDifference conserved_split_difference(
    const CompiledSystem &, const EquilibriumResult &,
    const CompiledSystem &, const EquilibriumResult &, std::size_t donor,
    const StabilityOptions &);
// The driver internals the stability search (phase_stability.cpp) shares.
void refuse(EquilibriumResult& result, Candidate as = Candidate::Refused);
inline bool held_refusal(EquilibriumResult& result, core::Stop stop, const char* cause, std::string message) {
    result.convergence.stop = stop;
    result.convergence.cause = cause;
    result.message = std::move(message);
    return false;
}
double largest_difference(const std::vector<double>& left, const std::vector<double>& right);
EPC_API std::vector<double> add_multipliers(const CompiledSystem& compiled, const std::vector<double>& point, EquilibriumResult* observations = nullptr,
    const PhasePoint* cached_prepared = nullptr, const PointEvaluation* cached_rows = nullptr);
EPC_API EquilibriumResult structural_result(const CompiledSystem& compiled, core::Scope& scope);
// The lower composition bound of each species in the optimizer's log chart (ALG-EQUILIBRIUM).
std::vector<double> composition_box(const CompiledSystem& compiled);
EPC_API EquilibriumResult run_driver(
    const CompiledSystem& compiled,
    std::vector<double> initial, EquilibriumResult result, core::Scope& scope,
    std::vector<double>* attempt = nullptr, const DriverLimits* limits = nullptr, Attempt kind = Attempt::Cold,
    std::shared_ptr<const OptimizerState> warm = {}, bool predicted = false, bool refine_trial = false, Retry runs = Retry::PhysicalAmounts,
    EquilibriumAttemptReporting reporting = {}
);
EPC_API Retry retry_after(Attempt, bool reaction, bool non_error, bool fatal, Ipopt::SolverReturn, TopologyEvent, core::Stop, Candidate);
// One physical-mol attempt; null restart rebuilds the original declaration, otherwise coordinates are immutable.
EPC_API EquilibriumResult recover_physical_amounts(
    const CompiledSystem& compiled, EquilibriumResult result, bool eligible, core::Scope& scope,
    const std::vector<double>* restart = nullptr, const DriverLimits* limits = nullptr,
    std::vector<double>* attempt = nullptr, bool admit_density = false,
    std::size_t* declarations = nullptr, std::size_t declaration_budget = 0, bool* invoked = nullptr, Work counts = Work::Accumulated,
    EquilibriumAttemptReporting reporting = {}
);
// Internal corrector route for continuation: the full compiled coordinate vector is the
// seed of the same Ipopt lifecycle; a numerical failure ends in its result. Not bound; solve_equilibrium
// remains the public solve.
EPC_API EquilibriumResult solve_equilibrium_from_coordinates(
    const CompiledSystem&, const std::vector<double>&, std::vector<double>* = nullptr,
    const DriverLimits* = nullptr, std::shared_ptr<const OptimizerState> = {});
EPC_API EquilibriumResult solve_equilibrium_from_coordinates(
    const HelmholtzModel&, const EquilibriumProblem&, const std::vector<double>&, std::vector<double>* = nullptr,
    const DriverLimits* = nullptr, std::shared_ptr<const OptimizerState> = {});

// Declared seeds carry no thermodynamic declarations: all share the one problem.
struct EquilibriumStart final {
    std::vector<double> coordinates, amounts;
    std::vector<std::vector<double>> compositions;
    std::vector<std::optional<double>> packing_fractions;
    std::optional<double> temperature, pressure;
    std::optional<EquilibriumResult> previous;
};
struct EquilibriumWorkLimits final {
    std::optional<std::size_t> attempts, iterations, attempt_iterations;
    std::optional<double> seconds, attempt_seconds;
};
struct EquilibriumStartsResult final {
    std::optional<std::size_t> selected;
    std::vector<std::optional<EquilibriumResult>> results;
    std::vector<std::string> causes;
    std::vector<double> attempt_seconds;
    std::size_t iterations = 0;
    double seconds = 0.0, overrun_seconds = 0.0;
    bool budget_exhausted = false;
};
EPC_API EquilibriumStartsResult solve_equilibrium_starts(const HelmholtzModel& model,
    const EquilibriumProblem& problem, const std::vector<EquilibriumStart>& starts,
    const EquilibriumWorkLimits& limits);

// ---------------------------------------------------------------- continuation

enum class ContinuationStop {
    reached_target, initial_failure, invalid_path, singular_jacobian, branch_lost,
    phase_disappeared, phases_coalesced, left_domain, infeasible, nonfinite, stalled,
    step_too_small, budget_exhausted, critical_endpoint, turning_point, held_species
};

EPC_API const char* continuation_stop_name(ContinuationStop stop) noexcept;

struct ContinuationOptions final {
    double initial_step = 0.0;
    double minimum_step = 0.0;
    std::size_t max_steps = 0;
    std::size_t max_retries_per_step = 0;
    // Envelope mode: largest |delta s| of the local coordinate and the near-critical |ln K| straddle threshold.
    double envelope_max_log_step = 0.0;
    double envelope_straddle = 0.0;
};

struct ContinuationRequest final {
    EquilibriumProblem start_problem;
    EquilibriumProblem target_problem;
    std::optional<EquilibriumResult> start_result;
    ContinuationOptions options;
    std::vector<double> parameter_targets;  // Physical endpoints in model.active_parameters() order.
};

struct ContinuationResult final {
    StageWork observations;
    bool success = false;
    std::string message;
    ContinuationStop stop = ContinuationStop::initial_failure;
    std::vector<double> path_values;
    std::vector<std::vector<double>> parameter_values;
    std::vector<EquilibriumResult> states;
    double reached_path = 0.0;
    std::optional<double> attempted_path;
    std::optional<std::size_t> event_phase;
    std::optional<std::pair<std::size_t, std::size_t>> event_pair;
    std::size_t attempts = 0;
    std::size_t accepted_steps = 0;
    std::size_t retries = 0;
    std::optional<EquilibriumResult> last_result;
    std::optional<ContinuationStop> cause_stop;
    std::optional<core::Convergence> cause_convergence;
    std::optional<std::string> cause_message;
    std::vector<EnvelopeCondition> conditions;        // envelope: the selected condition at each accepted state
    std::vector<EquilibriumResult> critical_points;   // envelope: critical-equation candidates in trace order
    std::vector<std::vector<double>> critical_seeds;  // their seeds: ln(T/1 K), ln(P/1 Pa), ln(rho/(1 mol/m^3))
};

// Declared-branch continuation along a straight path in fixed intensives and resolved feed
// amounts, with a tangent predictor and the one Ipopt driver as corrector. A start problem with an
// envelope condition traces that envelope instead (ALG-EQUILIBRIUM); its target problem is unused.
EPC_API ContinuationResult continue_equilibrium(
    const HelmholtzModel& model, const ContinuationRequest& request
);

template<class Base, class Direction>
std::vector<Base> shifted_coordinates(std::vector<Base> point, const Direction& direction, Base scale) {
    for (std::size_t i = 0; i < point.size(); ++i) point[i] += scale * direction[i];
    return point;
}

enum class NewtonTrial { Refused, Accepted, Stop, Interrupted };

template<class Base, class Direction, class Trial>
NewtonTrial damped_newton_step(const std::vector<Base>& point, const Direction& direction, Base contraction_limit,
    const Trial& trial, bool fourth_trial = true) {
    const auto coordinates = [&](Base alpha) { return shifted_coordinates(point, direction, alpha); };
    const auto contracts = [&](Base norm) { return norm < contraction_limit; };
    constexpr double damping[] = {1.0, 0.5, 0.25, 0.125};
    for (std::uint32_t i = 0; i < (fourth_trial ? 4u : 3u); ++i)
        if (auto outcome = trial(damping[i], i, coordinates, contracts); outcome != NewtonTrial::Refused) return outcome;
    return NewtonTrial::Refused;
}

// J dy = -F_theta for one column F_theta at a prepared point; empty when J is singular or the solve is not admitted.
EPC_API std::vector<double> solve_implicit(
    const CompiledSystem& compiled,
    const PhasePoint& prepared,
    const std::vector<double>& column, std::optional<core::Convergence>* failure = nullptr, std::vector<double> jacobian = {}
);
EPC_API bool check_held_responses(const CompiledSystem&, const std::function<Eigen::MatrixXd()>&,
    const Eigen::VectorXd&, EquilibriumResult&);

EPC_API std::vector<double> sensitivity(
    const CompiledSystem& compiled,
    const std::vector<double>& point,
    const std::vector<double>& lever, std::optional<core::Convergence>* failure = nullptr
);

class GibbsMinimizer final : public Ipopt::TNLP {
public:
    EPC_API GibbsMinimizer(const CompiledSystem& compiled, std::vector<double> start,
                   core::Scope* scope = nullptr, const DriverLimits* limits = nullptr, bool physical_amounts = false);
    EPC_API ~GibbsMinimizer() override;
    EPC_API bool get_nlp_info(Ipopt::Index& n, Ipopt::Index& m,
                      Ipopt::Index& nnz_jac_g, Ipopt::Index& nnz_h_lag,
                      IndexStyleEnum& index_style) override;
    bool get_bounds_info(Ipopt::Index n, Ipopt::Number* x_l, Ipopt::Number* x_u,
                         Ipopt::Index m, Ipopt::Number* g_l, Ipopt::Number* g_u) override;
    bool get_starting_point(Ipopt::Index n, bool init_x, Ipopt::Number* x,
                            bool init_z, Ipopt::Number* z_L, Ipopt::Number* z_U,
                            Ipopt::Index m, bool init_lambda, Ipopt::Number* lambda) override;
    EPC_API bool eval_f(Ipopt::Index n, const Ipopt::Number* x, bool new_x,
                Ipopt::Number& value) override;
    EPC_API bool eval_grad_f(Ipopt::Index n, const Ipopt::Number* x, bool new_x,
                     Ipopt::Number* grad_f) override;
    EPC_API bool eval_g(Ipopt::Index n, const Ipopt::Number* x, bool new_x,
                Ipopt::Index m, Ipopt::Number* g) override;
    EPC_API bool eval_jac_g(Ipopt::Index n, const Ipopt::Number* x, bool new_x,
                    Ipopt::Index m, Ipopt::Index nele_jac, Ipopt::Index* iRow,
                    Ipopt::Index* jCol, Ipopt::Number* values) override;
    EPC_API bool eval_h(Ipopt::Index n, const Ipopt::Number* x, bool new_x,
                Ipopt::Number obj_factor, Ipopt::Index m, const Ipopt::Number* lambda,
                bool new_lambda, Ipopt::Index nele_hess, Ipopt::Index* iRow,
                Ipopt::Index* jCol, Ipopt::Number* values) override;
    bool intermediate_callback(Ipopt::AlgorithmMode mode, Ipopt::Index iter,
                                Ipopt::Number obj_value, Ipopt::Number inf_pr,
                                Ipopt::Number inf_du, Ipopt::Number mu,
                                Ipopt::Number d_norm, Ipopt::Number regularization_size,
                                Ipopt::Number alpha_du, Ipopt::Number alpha_pr,
                                Ipopt::Index ls_trials, const Ipopt::IpoptData* ip_data,
                                Ipopt::IpoptCalculatedQuantities* ip_cq) override;
    void finalize_solution(Ipopt::SolverReturn status, Ipopt::Index n,
                           const Ipopt::Number* x, const Ipopt::Number* z_L,
                           const Ipopt::Number* z_U, Ipopt::Index m,
                           const Ipopt::Number* g, const Ipopt::Number* lambda,
                           Ipopt::Number objective,
                           const Ipopt::IpoptData* ip_data,
                           Ipopt::IpoptCalculatedQuantities* ip_cq) override;

    bool initialize_warm(const OptimizerState&, double cap, bool predicted);
    bool physical_bounds(std::vector<double>& lower, std::vector<double>& upper);
    std::shared_ptr<const OptimizerState> optimizer_state() const { return endpoint_; }
    const std::vector<double>& solution() const noexcept { return solution_; }
    const std::vector<double>& constraint_multipliers() const noexcept {
        return constraint_multipliers_;
    }
    // Ipopt's lower-bound multipliers at the endpoint, per physical column of the extensive objective.
    const std::vector<double>& lower_multipliers() const noexcept { return lower_multipliers_; }
    const std::optional<PhasePoint>& last_prepared() const noexcept { return prepared_; }
    const std::optional<PointEvaluation>& last_rows() const noexcept { return rows_; }
    std::vector<double> last_classified() const { return rows_ && prepared_ ? compiled_.classified_residuals(*prepared_, rows_->values) : std::vector<double>{}; }
    Ipopt::SolverReturn status() const noexcept { return status_; }
    std::size_t iterations() const noexcept { return iterations_; }
    std::size_t evaluations() const noexcept { return evaluations_; }
    std::size_t hessian_calls() const noexcept { return hessian_calls_; }
    std::size_t prepared_points() const noexcept { return prepared_points_; }  // distinct callback points

    double regularization_size() const noexcept { return regularization_size_; }
    const std::optional<std::string>& step_rule_error() const noexcept { return step_rule_error_; }
    bool amount_limited() const { return !amount_columns_.empty(); }
    void begin_ipopt() { observations_.ipopt_diagnostics.emplace(); }
    void finish_ipopt(int initialized, int returned, bool predicted, bool physical) noexcept {
        auto& d = *observations_.ipopt_diagnostics;
        d.initialization_return = initialized; d.application_return = returned;
        d.warm_initialized = warm_; d.predicted = predicted; d.physical_amounts = physical;
    }
    void publish_observations(EquilibriumResult& result) const { observations_.report_observations(result); }

private:
    std::vector<std::size_t> constraint_rows_;
    // One latest-point phase evaluation per solve; every derived quantity is cached
    // with it and invalidated when the callback point changes.
    void retain_callback_failure(const std::exception& error, const char* operation) {
        if (dynamic_cast<const BudgetExhausted*>(&error)) return;
        core::retain_failure(observations_.first_numerical_failure, error, operation);
    }
    bool step_admitted(const Ipopt::Number* x);
    const PhasePoint& prepared(Ipopt::Index n, const Ipopt::Number* x, bool hessian = false);
    const PointEvaluation& rows_at(Ipopt::Index n, const Ipopt::Number* x);

    const CompiledSystem& compiled_;
    std::vector<double> start_, lower_multipliers_, warm_lower_, warm_upper_, warm_lambda_;
    std::shared_ptr<const OptimizerState> endpoint_;
    bool warm_ = false, restored_ = false;

    std::optional<PhasePoint> prepared_;
    std::optional<PointEvaluation> rows_;
    std::optional<std::vector<double>> objective_hessian_;
    std::optional<std::vector<double>> constraint_hessian_;
    std::vector<double> constraint_hessian_multipliers_;
    StageWork observations_;
    std::vector<double> solution_;
    std::vector<double> constraint_multipliers_;
    Ipopt::SolverReturn status_ = Ipopt::UNASSIGNED;
    std::size_t iterations_ = 0;
    std::size_t evaluations_ = 0;
    std::size_t hessian_calls_ = 0;
    std::size_t prepared_points_ = 0;
    core::Scope* scope_ = nullptr;
    double regularization_size_ = 0.0;
    const DriverLimits* limits_ = nullptr;
    std::size_t consumed_before_ = limits_ ? limits_->consumed : 0;  // the earlier stages' iterations
    std::vector<double> base_, lower_, z_lower_, z_upper_;
    std::vector<std::size_t> amount_columns_;  // the free phase amounts
    std::vector<double> chart_;  // per-column unit of the TNLP coordinates: physical = chart_ * TNLP
    std::optional<std::string> step_rule_error_;
};

}  // namespace epcsaft::equilibrium
