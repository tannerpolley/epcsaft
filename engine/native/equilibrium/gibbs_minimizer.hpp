#pragma once

#include "equilibrium/problem.hpp"

#include <IpTNLP.hpp>

#include <cstddef>
#include <chrono>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace epcsaft::core { class Scope; }

namespace epcsaft::equilibrium {

// Work actually performed on distinct callback points within one solve.
struct EvaluationWork final {
    std::size_t prepared_points = 0;
    std::size_t row_evaluations = 0;
    std::size_t jacobians = 0;
    std::size_t gradients = 0;
    std::size_t objective_hessians = 0;
    std::size_t constraint_hessians = 0;
};

struct NlpShape final {
    Ipopt::Index variables = 0;
    Ipopt::Index constraints = 0;
    Ipopt::Index jacobian_entries = 0;
    Ipopt::Index hessian_entries = 0;
};

// `speciate` false: a previous result replaces the start's compositions, so a reacting phase is not speciated.
struct DriverLimits;
std::vector<double> boundary_start(const CompiledSystem& compiled, bool speciate = true, const DriverLimits* limits = nullptr);
bool supports_split_seed(const CompiledSystem& compiled);
std::string declaration_difference(const EquilibriumProblem& a, const EquilibriumProblem& b, bool same_seeds = true);
inline bool admitted_packing_seed(double eta) {
    return std::isfinite(eta) && eta >= policy::PACKING_FRACTION_FLOOR.value
        && eta <= policy::PACKING_FRACTION_CEILING.value;
}

// The producer's candidate-topology predicate on full coordinates: first vanished free phase
// in declaration order, else the first coincident pair; shared by the driver and continuation.
struct TopologyCheck final {
    TopologyEvent event = TopologyEvent::None;
    std::optional<std::size_t> phase;
    std::optional<std::pair<std::size_t, std::size_t>> pair;
    std::vector<double> distances;
    double distance = 0.0;
};

TopologyCheck check_topology(const CompiledSystem& compiled, const std::vector<double>& coordinates);

// One density-branch admission owner for public solves and anchored continuation.
struct BranchCheck final {
    bool accepted = false;
    bool terminal = false;
    std::string message;
};

BranchCheck check_density_branches(
    const HelmholtzModel& model, const CompiledSystem& compiled, const EquilibriumResult& result,
    const std::vector<double>& anchors, bool initial, core::Scope& scope, const DriverLimits* limits = nullptr
);

struct DriverLimits final {
    std::size_t iterations;
    std::chrono::steady_clock::time_point deadline;
    mutable std::size_t consumed = 0;
    bool expired() const { return std::chrono::steady_clock::now() >= deadline; }
    void check() const { if (expired()) throw std::runtime_error("elapsed work budget exhausted"); }
};
// Internal corrector route for continuation: the full compiled coordinate vector is the
// seed of the same Ipopt lifecycle. Not bound; solve_equilibrium remains the public solve.
EquilibriumResult solve_equilibrium_from_coordinates(
    const HelmholtzModel& model,
    const EquilibriumProblem& problem,
    const std::vector<double>& coordinates, std::vector<double>* attempt = nullptr, const DriverLimits* limits = nullptr
);

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
EquilibriumStartsResult solve_equilibrium_starts(const HelmholtzModel& model,
    const EquilibriumProblem& problem, const std::vector<EquilibriumStart>& starts,
    const EquilibriumWorkLimits& limits);

// ---------------------------------------------------------------- continuation

enum class ContinuationStop {
    reached_target, initial_failure, invalid_path, singular_jacobian, branch_lost,
    phase_disappeared, phases_coalesced, left_domain, infeasible, nonfinite, stalled,
    step_too_small, budget_exhausted
};

const char* continuation_stop_name(ContinuationStop stop) noexcept;

struct ContinuationOptions final {
    double initial_step = 0.0;
    double minimum_step = 0.0;
    std::size_t max_steps = 0;
    std::size_t max_retries_per_step = 0;
};

struct ContinuationRequest final {
    EquilibriumProblem start_problem;
    EquilibriumProblem target_problem;
    std::optional<EquilibriumResult> start_result;
    ContinuationOptions options;
    std::vector<double> parameter_targets;  // Physical endpoints in model.active_parameters() order.
};

struct ContinuationResult final {
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
};

// Declared-branch continuation along a straight path in fixed intensives and resolved feed
// amounts, with a tangent predictor and the one Ipopt driver as corrector.
ContinuationResult continue_equilibrium(
    const HelmholtzModel& model, const ContinuationRequest& request
);

std::vector<double> solve_implicit(
    const CompiledSystem& compiled,
    const std::vector<double>& point,
    const std::vector<double>& parameter_rows
);

std::vector<double> sensitivity(
    const CompiledSystem& compiled,
    const std::vector<double>& point,
    const std::vector<double>& lever
);

struct SensitivityInfo final {
    std::size_t equations = 0;
    std::size_t parameters = 0;
    std::size_t rank = 0;
    double condition_estimate = 0.0;
    bool nonsingular = false;
};

SensitivityInfo inspect_sensitivity(
    const CompiledSystem& compiled,
    const std::vector<double>& point,
    const std::vector<double>& parameter_rows
);

class GibbsMinimizer final : public Ipopt::TNLP {
public:
    GibbsMinimizer(const CompiledSystem& compiled, std::vector<double> start,
                   core::Scope* scope = nullptr, const DriverLimits* limits = nullptr);
    bool get_nlp_info(Ipopt::Index& n, Ipopt::Index& m,
                      Ipopt::Index& nnz_jac_g, Ipopt::Index& nnz_h_lag,
                      IndexStyleEnum& index_style) override;
    bool get_bounds_info(Ipopt::Index n, Ipopt::Number* x_l, Ipopt::Number* x_u,
                         Ipopt::Index m, Ipopt::Number* g_l, Ipopt::Number* g_u) override;
    bool get_starting_point(Ipopt::Index n, bool init_x, Ipopt::Number* x,
                            bool init_z, Ipopt::Number* z_L, Ipopt::Number* z_U,
                            Ipopt::Index m, bool init_lambda, Ipopt::Number* lambda) override;
    bool eval_f(Ipopt::Index n, const Ipopt::Number* x, bool new_x,
                Ipopt::Number& value) override;
    bool eval_grad_f(Ipopt::Index n, const Ipopt::Number* x, bool new_x,
                     Ipopt::Number* grad_f) override;
    bool eval_g(Ipopt::Index n, const Ipopt::Number* x, bool new_x,
                Ipopt::Index m, Ipopt::Number* g) override;
    bool eval_jac_g(Ipopt::Index n, const Ipopt::Number* x, bool new_x,
                    Ipopt::Index m, Ipopt::Index nele_jac, Ipopt::Index* iRow,
                    Ipopt::Index* jCol, Ipopt::Number* values) override;
    bool eval_h(Ipopt::Index n, const Ipopt::Number* x, bool new_x,
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

    const std::vector<double>& solution() const noexcept { return solution_; }
    const std::vector<double>& constraint_multipliers() const noexcept {
        return constraint_multipliers_;
    }
    const std::optional<PointEvaluation>& last_rows() const noexcept { return rows_; }
    std::vector<double> last_classified() const { return rows_ && prepared_ ? compiled_.classified_residuals(*prepared_, rows_->values) : std::vector<double>{}; }
    Ipopt::SolverReturn status() const noexcept { return status_; }
    std::size_t iterations() const noexcept { return iterations_; }
    std::size_t evaluations() const noexcept { return evaluations_; }
    std::size_t hessian_calls() const noexcept { return hessian_calls_; }

    double regularization_size() const noexcept { return regularization_size_; }
    const EvaluationWork& work() const noexcept { return work_; }
    const std::optional<std::string>& preparation_error() const noexcept {
        return preparation_error_;
    }

    NlpShape shape() const {
        const Ipopt::Index variables = static_cast<Ipopt::Index>(
            compiled_.pinned_count == 0
                ? compiled_.primal_size
                : compiled_.primal_size + compiled_.rank + compiled_.charge_row_count);
        const Ipopt::Index constraints = static_cast<Ipopt::Index>(constraint_rows().size());
        return {variables, constraints, variables * constraints,
                variables * (variables + 1) / 2};
    }

    std::size_t primal_size() const noexcept { return compiled_.primal_size; }

    bool is_square() const noexcept { return compiled_.is_square(); }

private:
    std::vector<std::size_t> constraint_rows() const;
    // One latest-point phase evaluation per solve; every derived quantity is cached
    // with it and invalidated when the callback point changes.
    const PhasePoint& prepared(Ipopt::Index n, const Ipopt::Number* x);
    const PointEvaluation& rows_at(Ipopt::Index n, const Ipopt::Number* x);

    const CompiledSystem& compiled_;
    std::vector<double> start_;
    std::optional<PhasePoint> prepared_;
    std::optional<PointEvaluation> rows_;
    std::optional<std::vector<double>> jacobian_;
    std::optional<std::vector<double>> gradient_;
    std::optional<std::vector<double>> objective_hessian_;
    std::optional<std::vector<double>> constraint_hessian_;
    std::vector<double> constraint_hessian_multipliers_;
    std::optional<std::string> preparation_error_;
    EvaluationWork work_;
    std::vector<double> solution_;
    std::vector<double> constraint_multipliers_;
    Ipopt::SolverReturn status_ = Ipopt::UNASSIGNED;
    std::size_t iterations_ = 0;
    std::size_t evaluations_ = 0;
    std::size_t hessian_calls_ = 0;
    core::Scope* scope_ = nullptr;
    double regularization_size_ = 0.0;
    const DriverLimits* limits_ = nullptr;
};

}  // namespace epcsaft::equilibrium
