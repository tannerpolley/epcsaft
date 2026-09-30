#include "regression/fit.hpp"

#include "core/numerical_policy.hpp"
#include "core/trace.hpp"

#include <ceres/ceres.h>
#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>

namespace epcsaft::regression {
namespace {

void require(bool condition, const char* what) {
    if (!condition) throw std::invalid_argument(std::string("fit ") + what);
}

bool finite(const std::vector<double>& values) {
    return std::all_of(values.begin(), values.end(), [](double value) { return std::isfinite(value); });
}

// Chain of one row: optimizer_jacobian[c] = sqrt(w) (dr/dy) sum_{r: tie r = c} dy/dtheta_r scale_c.
void chain(const FitRequest& request, const std::vector<std::size_t>& tie, std::size_t row,
           double sqrt_weight, double prediction, const double* gradient, double* out) {
    const auto& coordinates = request.coordinates.coordinates;
    const double slope = sqrt_weight * residual_slope(request.observations[row], prediction);
    std::fill(out, out + coordinates.size(), 0.0);
    for (std::size_t r = 0; r < tie.size(); ++r)
        out[tie[r]] += slope * gradient[r] * coordinates[tie[r]].scale;
}

struct Evaluation final {
    std::vector<double> predictions, residuals, gradient;  // gradient: n x R row-major
    std::vector<ObservationStatus> statuses;
    std::vector<double> weighted, optimizer_gradient;
    double cost = 0.0;
};

// Each input row names the first row with the same state; this is the sole membership map.
std::vector<std::size_t> state_groups(const std::vector<Observation>& observations) {
    std::vector<std::size_t> groups(observations.size());
    // ponytail: quadratic declaration comparison; index exact state records if fit datasets grow large.
    for (std::size_t i = 0; i < observations.size(); ++i) {
        groups[i] = i;
        for (std::size_t j = 0; j < i; ++j)
            if (same_state(observations[i], observations[j])) {
                groups[i] = groups[j];
                break;
            }
    }
    return groups;
}

struct Trial final {
    std::unique_ptr<const HelmholtzModel> model;  // declared first: compiled centers are destroyed first
    std::vector<Observation> observations;
    std::vector<Prediction> centers;
    std::vector<double> optimizer;
    Evaluation evaluation;
    bool jacobian = false;
};

// One successful point owns its immutable model and compiled central states through later promotion.
class Evaluator final {
public:
    Evaluator(const FitRequest& request, FitResult& result, const std::vector<std::size_t>& groups,
              const std::vector<double>& sqrt_weights, Trial prepared = {})
        : request_(request), result_(result), groups_(groups), sqrt_weights_(sqrt_weights),
          active_(request.coordinates.active_parameters()), tie_(request.coordinates.tie()),
          starts_(groups.size()), prepared_(std::move(prepared)) {}

    bool evaluate(const std::vector<double>& optimizer, bool jacobian, Evaluation& out) {
        core::Scope scope("fit_evaluation");
        scope.probe("optimizer", optimizer);
        scope.probe("jacobian", jacobian ? 1.0 : 0.0);
        const bool same_point = cached_ && cached_->optimizer == optimizer;
        const bool hit = same_point && (!jacobian || cached_->jacobian);
        scope.probe("reuse", hit ? (cached_->jacobian ? 2.0 : 1.0) : same_point ? 3.0 : 0.0);
        if (hit) {
            out = cached_->evaluation;
            return true;
        }
        ++trials_;
        (jacobian ? result_.jacobian_evaluations : result_.residual_evaluations) += 1;
        const std::size_t n = request_.observations.size(), width = tie_.size();
        Trial fresh = std::move(prepared_);
        Trial& trial = same_point ? *cached_ : fresh;
        try {
            auto observations = request_.observations;
            request_.coordinates.apply(observations, optimizer);
            auto physical = request_.coordinates.physical(optimizer);
            physical.resize(active_.size());
            // Reuse validation only at exact EOS values and rewritten reaction declarations.
            if (!same_point && (!trial.model || trial.model->active_parameter_values() != physical
                || !std::equal(observations.begin(), observations.end(), trial.observations.begin(), trial.observations.end(),
                    [](const auto& x, const auto& y) { return pressure_closed(x.kind) || same_state(x, y); }))) {
                trial.model = std::make_unique<HelmholtzModel>(request_.coordinates.apply(request_.parameters, optimizer),
                                                              active_, request_.thermochemistry);
                trial.centers.assign(groups_.size(), Prediction{});
            }
            trial.observations = std::move(observations);
            if (!evaluate_rows(optimizer, jacobian, out, trial, n, width)) return false;
            out.weighted.resize(n);
            for (std::size_t i = 0; i < n; ++i) out.weighted[i] = sqrt_weights_[i] * out.residuals[i];
            // One half of the row-ordered squared sum; reject its first nonfinite operand or sum.
            double sum = 0.0;
            for (std::size_t i = 0; i < n; ++i) {
                const double weighted = out.weighted[i];
                sum += weighted * weighted;
                if (!std::isfinite(weighted) || !std::isfinite(sum))
                    return fail(i, optimizer, ObservationStatus::Nonfinite,
                        std::isfinite(weighted) ? "accumulated weighted cost not finite" : "weighted residual not finite", jacobian);
            }
            out.cost = 0.5 * sum;
            const std::size_t columns = request_.coordinates.coordinates.size();
            if (jacobian) out.optimizer_gradient.resize(n * columns);
            for (std::size_t i = 0; jacobian && i < n; ++i) {
                chain(request_, tie_, i, sqrt_weights_[i], out.predictions[i],
                      out.gradient.data() + i * width, out.optimizer_gradient.data() + i * columns);
                if (!std::all_of(out.optimizer_gradient.begin() + static_cast<std::ptrdiff_t>(i * columns),
                                 out.optimizer_gradient.begin() + static_cast<std::ptrdiff_t>((i + 1) * columns),
                                 [](double value) { return std::isfinite(value); }))
                    return fail(i, optimizer, ObservationStatus::Nonfinite, "optimizer Jacobian not finite", true);
            }
        } catch (const std::exception& error) {
            return fail(n, optimizer, ObservationStatus::InadmissibleTrial, error.what(), jacobian);
        }
        // Neither a failed row nor an overflowing weighted/coordinate chain publishes starts or a cache.
        for (std::size_t g = 0; g < groups_.size(); ++g)
            if (trial.centers[g].solved) starts_[g] = trial.centers[g].solved;
        trial.optimizer = optimizer;
        trial.evaluation = out;
        trial.jacobian = jacobian;
        if (!same_point) {
            cached_.reset();
            cached_.emplace(std::move(fresh));
        }
        return true;
    }

    bool gradient_at(const std::vector<double>& point, std::vector<double>& out) {
        Evaluation trial;
        if (!evaluate(point, true, trial)) return false;
        const std::size_t width = request_.coordinates.coordinates.size();
        out.assign(width, 0.0);
        for (std::size_t i = 0; i < trial.weighted.size(); ++i) {
            for (std::size_t c = 0; c < width; ++c)
                out[c] += trial.optimizer_gradient[i * width + c] * trial.weighted[i];
            if (!finite(out))
                return fail(i, point, ObservationStatus::Nonfinite, "weighted objective gradient not finite", true);
        }
        return true;
    }

    bool fail(std::size_t observation, const std::vector<double>& optimizer, ObservationStatus status,
              std::string message, bool jacobian) {
        result_.failures.push_back({trials_, observation, optimizer, status, std::move(message), jacobian});
        ++result_.trial_failures;
        return false;
    }

private:
    bool evaluate_rows(const std::vector<double>& optimizer, bool jacobian, Evaluation& out,
                       Trial& trial, std::size_t n, std::size_t width) {
        const auto& model = *trial.model;
        out = Evaluation{};
        out.predictions.assign(n, 0.0);
        out.residuals.assign(n, 0.0);
        if (jacobian) out.gradient.assign(n * width, 0.0);
        std::vector<std::vector<PredictionGradient>> gradients(n);
        std::vector<std::size_t> position(n, 0);
        for (std::size_t i = 0; i < n; ++i) {
            const auto& observation = trial.observations[i];
            const std::size_t g = groups_[i];
            auto& prediction = trial.centers[g];
            if (!prediction.value) {
                core::Scope scope("fit_group_prediction");
                scope.probe("anchor", static_cast<double>(g));
                scope.probe("targets", static_cast<double>(std::count(groups_.begin(), groups_.end(), g)));
                prediction = predict(model, observation, starts_[g], std::move(prediction.compiled));
                if (prediction.status != ObservationStatus::Available && starts_[g])
                    prediction = predict(model, observation, std::nullopt, std::move(prediction.compiled));
            } else if (!pressure_closed(observation.kind)) project_prediction(model, observation, prediction);
            out.statuses.push_back(prediction.status);
            if (prediction.status != ObservationStatus::Available)
                return fail(i, optimizer, prediction.status, prediction.message, jacobian);
            out.predictions[i] = *prediction.value;
            out.residuals[i] = *prediction.residual;
            if (!jacobian) continue;
            if (gradients[g].empty()) {
                core::Scope scope("fit_group_differentiation");
                scope.probe("anchor", static_cast<double>(g));
                std::vector<const Observation*> rows;
                for (std::size_t row = 0; row < n; ++row)
                    if (groups_[row] == g) rows.push_back(&trial.observations[row]);
                gradients[g] = pressure_closed(observation.kind)
                    ? std::vector<PredictionGradient>{differentiate(model, observation, prediction)}
                    : differentiate_group(model, rows, prediction, request_.coordinates);
            }
            const auto& gradient = gradients[g][position[g]++];
            if (gradient.status != ObservationStatus::Available)
                return fail(i, optimizer, gradient.status, gradient.message, true);
            std::copy(gradient.value.begin(), gradient.value.end(),
                      out.gradient.begin() + static_cast<std::ptrdiff_t>(i * width));
        }
        return true;
    }

    const FitRequest& request_;
    FitResult& result_;
    const std::vector<std::size_t>& groups_;
    const std::vector<double>& sqrt_weights_;
    std::vector<ActiveParameter> active_;
    std::vector<std::size_t> tie_;
    std::vector<std::optional<equilibrium::EquilibriumResult>> starts_;
    Trial prepared_;
    std::optional<Trial> cached_;
    std::size_t trials_ = 0;
};

std::vector<int> select_held(const std::vector<double>& point, const std::vector<double>& gradient,
                            const std::vector<double>& lower, const std::vector<double>& upper) {
    std::vector<int> held;
    for (std::size_t c = 0; c < point.size(); ++c)
        if ((point[c] - lower[c] < policy::FIT_ACTIVE_BOUND_DISTANCE.value && gradient[c] > 0.0)
            || (upper[c] - point[c] < policy::FIT_ACTIVE_BOUND_DISTANCE.value && gradient[c] < 0.0))
            held.push_back(static_cast<int>(c));
    return held;
}

// The optimizer's chart of the box, in place of Ceres' bounds: held coordinates are fixed; a step leaving
// the box is shortened along its own direction to the first bound it meets (clamping one component keeps
// the rest of a step computed as if that bound could be crossed), clamping coordinates already on it. As
// the callback it ends a round on accepted new bound contact or changed held membership, including release.
class BoxChart final : public ceres::Manifold, public ceres::IterationCallback {
public:
    BoxChart(const std::vector<double>& u, const std::vector<double>& lower, const std::vector<double>& upper,
             const std::vector<int>& held,
             Evaluator& evaluator)
        : u_(u), lower_(lower), upper_(upper), held_(held), evaluator_(evaluator) {
        for (int c = 0; c < int(lower.size()); ++c)
            if (std::find(held.begin(), held.end(), c) == held.end()) free_.push_back(c);
        for (const int c : free_) at_bound_.push_back(side(u[c], c));
    }
    int side(double value, int c) const { return value == lower_[c] ? -1 : value == upper_[c] ? 1 : 0; }  // Plus writes bounds exactly
    ceres::CallbackReturnType operator()(const ceres::IterationSummary& summary) override {
        core::Scope scope("fit_iteration");
        scope.probe("accepted", summary.step_is_successful ? 1.0 : 0.0, summary.iteration);
        scope.probe("cost", summary.cost, summary.iteration);
        scope.probe("ratio", summary.relative_decrease, summary.iteration);
        scope.probe("trust_radius", summary.trust_region_radius, summary.iteration);
        if (summary.iteration > 0 && !summary.step_is_successful) return ceres::SOLVER_CONTINUE;
        std::vector<double> gradient;
        // Ceres evaluates the accepted Jacobian and updates u_ before user callbacks: this is a cache hit.
        if (!evaluator_.gradient_at(u_, gradient)) return ceres::SOLVER_ABORT;
        const auto next = select_held(u_, gradient, lower_, upper_);
        if (scope.on()) {
            std::vector<double> mask(u_.size(), 0.0), projected(u_.size());
            for (const int c : next) mask[c] = 1.0;
            for (std::size_t c = 0; c < u_.size(); ++c)
                projected[c] = u_[c] - std::clamp(u_[c] - gradient[c], lower_[c], upper_[c]);
            scope.probe("optimizer", u_, summary.iteration);
            scope.probe("held", mask, summary.iteration);
            scope.probe("weighted_gradient", gradient, summary.iteration);
            scope.probe("projected_gradient", projected, summary.iteration);
        }
        bool contact = false;
        for (std::size_t k = 0; summary.iteration > 0 && summary.step_is_successful && k < free_.size(); ++k) {
            contact = contact || (side(u_[free_[k]], free_[k]) != 0 && side(u_[free_[k]], free_[k]) != at_bound_[k]);
            at_bound_[k] = side(u_[free_[k]], free_[k]);
        }
        return summary.iteration > 0 && (contact || next != held_)
            ? ceres::SOLVER_TERMINATE_SUCCESSFULLY : ceres::SOLVER_CONTINUE;
    }
    int AmbientSize() const override { return static_cast<int>(lower_.size()); }
    int TangentSize() const override { return static_cast<int>(free_.size()); }
    bool Plus(const double* x, const double* delta, double* out) const override {
        double alpha = 1.0;
        std::size_t limit = free_.size();
        for (std::size_t k = 0; k < free_.size(); ++k) {
            const int c = free_[k];
            const double bound = delta[k] < 0.0 ? lower_[c] : upper_[c];
            if (std::abs(bound - x[c]) >= policy::FIT_ACTIVE_BOUND_DISTANCE.value
                && (x[c] + delta[k] - bound) * delta[k] > 0.0 && (bound - x[c]) / delta[k] < alpha) {
                alpha = (bound - x[c]) / delta[k];
                limit = k;
            }
        }
        std::copy(x, x + lower_.size(), out);
        for (std::size_t k = 0; k < free_.size(); ++k)
            out[free_[k]] = std::clamp(x[free_[k]] + alpha * delta[k], lower_[free_[k]], upper_[free_[k]]);
        if (limit < free_.size()) out[free_[limit]] = delta[limit] < 0.0 ? lower_[free_[limit]] : upper_[free_[limit]];
        return true;
    }
    bool PlusJacobian(const double*, double* jacobian) const override {
        std::fill(jacobian, jacobian + lower_.size() * free_.size(), 0.0);
        for (std::size_t k = 0; k < free_.size(); ++k) jacobian[free_[k] * free_.size() + k] = 1.0;
        return true;
    }
    bool Minus(const double* y, const double* x, double* out) const override {
        for (std::size_t k = 0; k < free_.size(); ++k) out[k] = y[free_[k]] - x[free_[k]];
        return true;
    }
    bool MinusJacobian(const double*, double* jacobian) const override {
        std::fill(jacobian, jacobian + lower_.size() * free_.size(), 0.0);
        for (std::size_t k = 0; k < free_.size(); ++k) jacobian[k * lower_.size() + free_[k]] = 1.0;
        return true;
    }

private:
    const std::vector<double>& u_;
    std::vector<double> lower_, upper_;
    std::vector<int> free_;
    std::vector<int> at_bound_;
    std::vector<int> held_;
    Evaluator& evaluator_;
};

class Cost final : public ceres::CostFunction {
public:
    Cost(const FitRequest& request, Evaluator& evaluator)
        : request_(request), evaluator_(evaluator) {
        set_num_residuals(static_cast<int>(request.observations.size()));
        mutable_parameter_block_sizes()->push_back(static_cast<int>(request.coordinates.coordinates.size()));
    }

    bool Evaluate(double const* const* parameters, double* residuals, double** jacobians) const override {
        const std::size_t width = request_.coordinates.coordinates.size();
        const std::vector<double> optimizer(parameters[0], parameters[0] + width);
        if (!request_.coordinates.within_bounds(optimizer)) return false;
        const bool jacobian = jacobians != nullptr && jacobians[0] != nullptr;
        Evaluation trial;
        if (!evaluator_.evaluate(optimizer, jacobian, trial)) return false;
        std::copy(trial.weighted.begin(), trial.weighted.end(), residuals);
        if (jacobian) std::copy(trial.optimizer_gradient.begin(), trial.optimizer_gradient.end(), jacobians[0]);
        return true;
    }

private:
    const FitRequest& request_;
    Evaluator& evaluator_;
};

const char* termination_name(ceres::TerminationType termination) {
    switch (termination) {
    case ceres::CONVERGENCE: return "convergence";
    case ceres::NO_CONVERGENCE: return "no_convergence";
    case ceres::FAILURE: return "failure";
    case ceres::USER_SUCCESS: return "user_success";
    case ceres::USER_FAILURE: return "user_failure";
    }
    return "failure";
}

void covariance(FitResult& result, const FitRequest& request, std::size_t n) {
    auto& cov = result.covariance;
    const auto& set = request.coordinates;
    const auto tie = set.tie();
    const auto width = static_cast<Eigen::Index>(set.coordinates.size());
    cov.assumption = "unit-variance weighted residuals, independent errors and a regular local "
                     "linearization at an interior point; multiply by variance_factor for an "
                     "unknown common variance; no claim under bounds, rank deficiency or model discrepancy";
    Eigen::MatrixXd jacobian(static_cast<Eigen::Index>(n), width);
    for (Eigen::Index i = 0; i < jacobian.rows(); ++i)
        for (Eigen::Index c = 0; c < width; ++c)
            jacobian(i, c) = result.optimizer_jacobian[static_cast<std::size_t>(i * width + c)];
    const Eigen::JacobiSVD<Eigen::MatrixXd> svd(jacobian, Eigen::ComputeThinU | Eigen::ComputeThinV);
    const auto& singular = svd.singularValues();
    cov.singular_values.assign(singular.data(), singular.data() + singular.size());
    const double floor = policy::FIT_RANK_RECIPROCAL_CONDITION.value * (singular.size() > 0 ? singular(0) : 0.0);
    cov.rank = static_cast<std::size_t>(std::count_if(cov.singular_values.begin(), cov.singular_values.end(),
        [&](double value) { return value > floor && value > 0.0; }));
    // A bound is active in optimizer space, within the fit's own active distance of the values it bounds.
    const auto [lower, upper] = set.optimizer_bounds();
    for (std::size_t c = 0; c < set.coordinates.size(); ++c)
        if (std::min(result.optimizer[c] - lower[c], upper[c] - result.optimizer[c]) < policy::FIT_ACTIVE_BOUND_DISTANCE.value)
            cov.active_bounds.push_back(c);
    // The variance factor counts only rows that enter the objective (w > 0).
    double weighted_sum = 0.0;
    std::size_t weighted_rows = 0;
    for (std::size_t i = 0; i < n; ++i) {
        weighted_sum += result.weighted_residuals[i] * result.weighted_residuals[i];
        weighted_rows += request.weights.empty() || request.weights[i] > 0.0 ? 1 : 0;
    }
    if (weighted_rows > cov.rank) cov.variance_factor = weighted_sum / static_cast<double>(weighted_rows - cov.rank);
    // An unfinished iterate keeps its rank, singular values and bound diagnostics but is not a
    // parameter-covariance estimate.
    if (result.status != "converged") { cov.status = "not_converged"; cov.variance_factor.reset(); return; }
    if (cov.rank < set.coordinates.size()) { cov.status = "rank_deficient"; return; }
    if (!cov.active_bounds.empty()) { cov.status = "active_bound"; return; }
    const Eigen::MatrixXd optimizer = svd.matrixV() * singular.cwiseAbs2().cwiseInverse().asDiagonal()
        * svd.matrixV().transpose();
    Eigen::MatrixXd map = Eigen::MatrixXd::Zero(static_cast<Eigen::Index>(tie.size()), width);
    for (std::size_t r = 0; r < tie.size(); ++r)
        map(static_cast<Eigen::Index>(r), static_cast<Eigen::Index>(tie[r])) = set.coordinates[tie[r]].scale;
    const Eigen::MatrixXd physical = map * optimizer * map.transpose();
    cov.optimizer.assign(optimizer.data(), optimizer.data() + optimizer.size());
    cov.physical.assign(physical.data(), physical.data() + physical.size());
    cov.status = finite(cov.optimizer) && finite(cov.physical) ? "available" : "unavailable";
    if (cov.status == "unavailable") { cov.optimizer.clear(); cov.physical.clear(); }
}

}  // namespace

FitResult fit(const FitRequest& request) {
    FitResult result;
    const auto& set = request.coordinates;
    const std::size_t n = request.observations.size();
    std::vector<double> sqrt_weights(n, 1.0);
    const auto groups = state_groups(request.observations);
    Trial prepared;
    try {
        validate(set, request.parameters, request.observations, request.thermochemistry);
        prepared.model = std::make_unique<HelmholtzModel>(
            request.parameters, set.active_parameters(), request.thermochemistry);
        const auto& model = *prepared.model;
        prepared.observations = request.observations;
        require(n > 0, "needs observations");
        core::Scope validation("fit_validation");
        prepared.centers.resize(groups.size());
        for (std::size_t i = 0; i < n; ++i) {
            auto& center = prepared.centers[groups[i]].compiled;
            const auto& observation = request.observations[i];
            if (!center && !pressure_closed(observation.kind) && observation.kind != ObservationKind::LogEquilibriumConstant)
                center = equilibrium::compile_problem(model, observation.problem);
            validate(observation, model, center ? &*center : nullptr);
        }
        require(request.weights.empty() || request.weights.size() == n, "weights width");
        for (std::size_t i = 0; i < request.weights.size(); ++i) {
            require(std::isfinite(request.weights[i]) && request.weights[i] >= 0.0, "weight");
            sqrt_weights[i] = std::sqrt(request.weights[i]);
        }
        const auto& controls = request.controls;
        require(std::isfinite(controls.function_tolerance) && controls.function_tolerance >= 0.0
                && std::isfinite(controls.gradient_tolerance) && controls.gradient_tolerance >= 0.0
                && std::isfinite(controls.parameter_tolerance) && controls.parameter_tolerance >= 0.0
                && controls.maximum_iterations > 0 && std::isfinite(controls.maximum_elapsed_time_seconds)
                && controls.maximum_elapsed_time_seconds > 0.0, "controls");
    } catch (const std::invalid_argument& error) {
        result.message = error.what();
        return result;
    }
    const std::size_t width = set.coordinates.size();
    Evaluator evaluator(request, result, groups, sqrt_weights, std::move(prepared));
    std::vector<double> optimizer = set.start(request.parameters, request.observations, request.thermochemistry);
    Evaluation initial;
    if (!evaluator.evaluate(optimizer, true, initial)) {
        result.status = "initial_evaluation_failed";
        result.message = result.failures.back().message;
        result.statuses = std::move(initial.statuses);  // the rows evaluated up to the failure
        return result;
    }
    result.initial_cost = initial.cost;

    Cost cost(request, evaluator);
    ceres::Problem::Options problem_options;
    problem_options.cost_function_ownership = ceres::DO_NOT_TAKE_OWNERSHIP;
    const auto [lower, upper] = set.optimizer_bounds();
    ceres::Solver::Options options;
    options.minimizer_type = ceres::TRUST_REGION;
    options.trust_region_strategy_type = ceres::LEVENBERG_MARQUARDT;
    options.linear_solver_type = ceres::DENSE_QR;
    options.dense_linear_algebra_library_type = ceres::EIGEN;
    options.num_threads = 1;
    options.logging_type = ceres::SILENT;
    options.minimizer_progress_to_stdout = false;
    const auto& controls = request.controls;
    options.function_tolerance = controls.function_tolerance > 0.0 ? controls.function_tolerance : policy::FIT_FUNCTION_TOLERANCE.value;
    options.gradient_tolerance = controls.gradient_tolerance > 0.0 ? controls.gradient_tolerance : policy::FIT_GRADIENT_TOLERANCE.value;
    options.parameter_tolerance = controls.parameter_tolerance > 0.0 ? controls.parameter_tolerance : policy::FIT_PARAMETER_TOLERANCE.value;
    options.update_state_every_iteration = true;
    // Bounds as an active set (a projected step alone stalls against a bound it wants to cross): a
    // coordinate on a bound (within FIT_ACTIVE_BOUND_DISTANCE) whose objective gradient points out of the
    // box is held where it is for a round. A round ends on Ceres convergence, accepted new bound contact,
    // or changed accepted held membership (including release); the set is rebuilt from the cached
    // Jacobian, g = J^T r through the same evaluation, checks and failure records as a Ceres trial.
    std::vector<double> gradient;
    ceres::Solver::Summary summary;
    std::vector<int> held;
    int remaining = controls.maximum_iterations;
    double seconds = controls.maximum_elapsed_time_seconds;
    bool settled = false, all_held = false, exhausted = false;
    int round = 0;
    for (bool first = true;; first = false) {
        if (!evaluator.gradient_at(optimizer, gradient)) {
            summary.termination_type = ceres::FAILURE;
            summary.message = result.failures.back().message;
            break;
        }
        const auto next = select_held(optimizer, gradient, lower, upper);
        // Unchanged membership retains Ceres' stopping criteria; it does not certify the free gradient.
        all_held = next.size() == width;
        settled = (!first && summary.termination_type == ceres::CONVERGENCE && next == held) || all_held;
        exhausted = !settled && (remaining <= 0 || seconds <= 0.0);
        if (settled || exhausted) break;
        held = next;
        core::Scope round_scope("fit_round");
        round_scope.probe("round", round++);
        if (round_scope.on()) {
            std::vector<double> mask(width, 0.0);
            for (const int c : held) mask[c] = 1.0;
            round_scope.probe("held", mask);
        }
        ceres::Problem problem(problem_options);
        problem.AddResidualBlock(&cost, nullptr, optimizer.data());
        auto* box = new BoxChart(optimizer, lower, upper, held, evaluator);  // owned by the problem
        problem.SetManifold(optimizer.data(), box);
        options.callbacks = {box};
        options.max_num_iterations = remaining;
        options.max_solver_time_in_seconds = seconds;
        ceres::Solve(options, &problem, &summary);
        round_scope.probe("termination", static_cast<double>(summary.termination_type));
        for (std::size_t k = first ? 0 : 1; k < summary.iterations.size(); ++k) {
            result.iteration_seconds.push_back(summary.iterations[k].iteration_time_in_seconds);
            result.iteration_costs.push_back(summary.iterations[k].cost);
        }
        remaining -= static_cast<int>(summary.iterations.size()) - 1;
        seconds -= summary.total_time_in_seconds;
        if (summary.termination_type != ceres::CONVERGENCE && summary.termination_type != ceres::USER_SUCCESS) break;
    }
    result.termination = settled ? "convergence" : exhausted ? "no_convergence" : termination_name(summary.termination_type);
    result.iterations = result.iteration_costs.size();
    result.message = all_held ? "every coordinate is held at an active bound"
        : exhausted ? (remaining <= 0 ? "iteration limit reached" : "time limit reached") : summary.message;
    result.status = settled ? "converged" : summary.termination_type == ceres::FAILURE
        || summary.termination_type == ceres::USER_FAILURE ? "failed" : "limit_reached";

    result.optimizer = optimizer;
    Evaluation final;
    if (!set.within_bounds(optimizer) || !evaluator.evaluate(optimizer, true, final)) {
        result.status = "failed";
        result.message = "the final point could not be evaluated: " + result.message;
        result.statuses = std::move(final.statuses);
        return result;
    }
    result.physical = set.physical(optimizer);
    result.parameters = set.apply(request.parameters, optimizer);
    result.predictions = final.predictions;
    result.residuals = final.residuals;
    result.statuses = final.statuses;
    result.physical_jacobian = final.gradient;
    result.weighted_residuals = final.weighted;
    result.optimizer_jacobian = final.optimizer_gradient;
    result.final_cost = final.cost;
    result.usable = result.status != "failed" && finite(result.physical);
    if (result.usable) covariance(result, request, n);
    return result;
}

}  // namespace epcsaft::regression
