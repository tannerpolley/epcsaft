#include "regression/fit.hpp"

#include "core/numerical_policy.hpp"

#include <ceres/ceres.h>
#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
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

struct Evaluation final {
    std::vector<double> predictions, residuals, gradient;  // gradient: n x R row-major
    std::vector<ObservationStatus> statuses;
};

// One trial: apply the optimizer point, rebuild the immutable model and evaluate every
// observation through the shared owners; the first typed failure rejects the trial.
class Evaluator final {
public:
    Evaluator(const FitRequest& request, FitResult& result)
        : request_(request), result_(result), active_(request.coordinates.active_parameters()),
          starts_(request.observations.size()) {}

    bool evaluate(const std::vector<double>& optimizer, bool jacobian, Evaluation& out) {
        // The last Jacobian evaluation serves every later request at the same point unchanged.
        if (jacobian && cached_ && cached_->first == optimizer) {
            out = cached_->second;
            return true;
        }
        ++trials_;
        (jacobian ? result_.jacobian_evaluations : result_.residual_evaluations) += 1;
        const std::size_t n = request_.observations.size(), width = active_.size();
        try {
            if (!evaluate_rows(optimizer, jacobian, out, n, width)) return false;
        } catch (const std::exception& error) {  // the model or an owner refused the trial by name
            return fail(n, optimizer, ObservationStatus::InadmissibleTrial, error.what(), jacobian);
        }
        if (jacobian) cached_.emplace(optimizer, out);
        return true;
    }

private:
    bool evaluate_rows(const std::vector<double>& optimizer, bool jacobian, Evaluation& out,
                       std::size_t n, std::size_t width) {
        const HelmholtzModel model(request_.coordinates.apply(request_.parameters, optimizer),
                                   active_, request_.thermochemistry);
        out = Evaluation{};
        out.predictions.assign(n, 0.0);
        out.residuals.assign(n, 0.0);
        if (jacobian) out.gradient.assign(n * width, 0.0);
        // Warm starts come only from a trial in which every observation succeeded (whether or not
        // the optimizer accepts its step), and a warm start that fails is retried cold.
        auto fresh = starts_;
        for (std::size_t i = 0; i < n; ++i) {
            const auto& observation = request_.observations[i];
            auto prediction = predict(model, observation, starts_[i]);
            if (prediction.status != ObservationStatus::Available && starts_[i].has_value())
                prediction = predict(model, observation, std::nullopt);
            out.statuses.push_back(prediction.status);
            if (prediction.solved.has_value()) fresh[i] = prediction.solved;
            if (prediction.status != ObservationStatus::Available)
                return fail(i, optimizer, prediction.status, prediction.message, jacobian);
            out.predictions[i] = *prediction.value;
            out.residuals[i] = *prediction.residual;
            if (!jacobian) continue;
            const auto gradient = differentiate(model, observation, prediction);
            if (gradient.status != ObservationStatus::Available)
                return fail(i, optimizer, gradient.status, gradient.message, true);
            std::copy(gradient.value.begin(), gradient.value.end(),
                      out.gradient.begin() + static_cast<std::ptrdiff_t>(i * width));
        }
        starts_ = std::move(fresh);
        return true;
    }

public:
    bool fail(std::size_t observation, const std::vector<double>& optimizer, ObservationStatus status,
              std::string message, bool jacobian) {
        result_.failures.push_back({trials_, observation, optimizer, status, std::move(message), jacobian});
        ++result_.trial_failures;
        return false;
    }

private:

    const FitRequest& request_;
    FitResult& result_;
    std::vector<ActiveParameter> active_;
    std::vector<std::optional<equilibrium::EquilibriumResult>> starts_;
    std::optional<std::pair<std::vector<double>, Evaluation>> cached_;
    std::size_t trials_ = 0;
};

// The objective the solver accumulates: one half of the running sum of squared weighted residuals
// in row order. Returns the first row at which the weighted residual, its square or the running sum
// stops being representable, else none; `cost` holds the objective only when every row is valid.
// The cause names the row's own residual when that is what overflowed, else the running sum.
const char* cost_cause(const double* weighted, std::size_t row, const char* prefix) {
    static const char* const names[4] = {"weighted residual not finite", "accumulated weighted cost not finite",
                                         "final weighted residual not finite", "final accumulated weighted cost not finite"};
    return names[(prefix != nullptr ? 2 : 0) + (std::isfinite(weighted[row]) ? 1 : 0)];
}

std::optional<std::size_t> objective(const double* weighted, std::size_t n, double& cost) {
    double sum = 0.0;
    cost = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        sum += weighted[i] * weighted[i];
        if (!std::isfinite(weighted[i]) || !std::isfinite(sum)) return i;
    }
    cost = 0.5 * sum;
    return std::nullopt;
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

// The optimizer's chart of the box, in place of Ceres' bounds: held coordinates are fixed; a step leaving
// the box is shortened along its own direction to the first bound it meets (clamping one component keeps
// the rest of a step computed as if that bound could be crossed), clamping coordinates already on it. As
// the iteration callback it ends a round when an accepted step puts a free coordinate on a new bound.
class BoxChart final : public ceres::Manifold, public ceres::IterationCallback {
public:
    BoxChart(const std::vector<double>& u, const std::vector<double>& lower, const std::vector<double>& upper,
             const std::vector<int>& held) : u_(u), lower_(lower), upper_(upper) {
        for (int c = 0; c < int(lower.size()); ++c)
            if (std::find(held.begin(), held.end(), c) == held.end()) free_.push_back(c);
        for (const int c : free_) at_bound_.push_back(side(u[c], c));
    }
    int side(double value, int c) const { return value == lower_[c] ? -1 : value == upper_[c] ? 1 : 0; }  // Plus writes bounds exactly
    ceres::CallbackReturnType operator()(const ceres::IterationSummary& summary) override {
        bool contact = false;
        for (std::size_t k = 0; summary.iteration > 0 && summary.step_is_successful && k < free_.size(); ++k) {
            contact = contact || (side(u_[free_[k]], free_[k]) != 0 && side(u_[free_[k]], free_[k]) != at_bound_[k]);
            at_bound_[k] = side(u_[free_[k]], free_[k]);
        }
        return contact ? ceres::SOLVER_TERMINATE_SUCCESSFULLY : ceres::SOLVER_CONTINUE;
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
};

class Cost final : public ceres::CostFunction {
public:
    Cost(const FitRequest& request, Evaluator& evaluator, std::vector<double> sqrt_weights)
        : request_(request), evaluator_(evaluator), sqrt_weights_(std::move(sqrt_weights)),
          tie_(request.coordinates.tie()) {
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
        const std::size_t n = trial.residuals.size();
        for (std::size_t i = 0; i < n; ++i) residuals[i] = sqrt_weights_[i] * trial.residuals[i];
        double cost = 0.0;
        if (const auto row = objective(residuals, n, cost))
            return evaluator_.fail(*row, optimizer, ObservationStatus::Nonfinite, cost_cause(residuals, *row, nullptr), jacobian);
        for (std::size_t i = 0; jacobian && i < n; ++i) {
            chain(request_, tie_, i, sqrt_weights_[i], trial.predictions[i],
                  trial.gradient.data() + i * tie_.size(), jacobians[0] + i * width);
            if (!std::all_of(jacobians[0] + i * width, jacobians[0] + (i + 1) * width,
                             [](double value) { return std::isfinite(value); }))
                return evaluator_.fail(i, optimizer, ObservationStatus::Nonfinite, "optimizer Jacobian not finite", true);
        }
        return true;
    }

private:
    const FitRequest& request_;
    Evaluator& evaluator_;
    std::vector<double> sqrt_weights_;
    std::vector<std::size_t> tie_;
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
    try {
        validate(set, request.parameters);
        const HelmholtzModel model(request.parameters, set.active_parameters(), request.thermochemistry);
        require(n > 0, "needs observations");
        for (const auto& observation : request.observations) validate(observation, model);
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
    Evaluator evaluator(request, result);
    std::vector<double> optimizer = set.start(request.parameters);
    Evaluation initial;
    if (!evaluator.evaluate(optimizer, true, initial)) {
        result.status = "initial_evaluation_failed";
        result.message = result.failures.back().message;
        result.statuses = std::move(initial.statuses);  // the rows evaluated up to the failure
        return result;
    }
    std::vector<double> initial_weighted(n);
    for (std::size_t i = 0; i < n; ++i) initial_weighted[i] = sqrt_weights[i] * initial.residuals[i];
    if (const auto row = objective(initial_weighted.data(), n, result.initial_cost)) {
        evaluator.fail(*row, optimizer, ObservationStatus::Nonfinite, cost_cause(initial_weighted.data(), *row, nullptr), true);
        result.status = "initial_evaluation_failed";
        result.message = result.failures.back().message;
        result.statuses = std::move(initial.statuses);
        return result;
    }

    Cost cost(request, evaluator, sqrt_weights);
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
    // box is held where it is for a round. A round ends when Ceres converges or an accepted step puts a
    // free coordinate on a bound (any other stop ends the fit); the set is rebuilt there from the cached
    // Jacobian, g = J^T r through the same evaluation, checks and failure records as a Ceres trial.
    std::vector<double> residuals(n), jacobian(n * width), gradient;
    const auto gradient_at = [&](const std::vector<double>& point, std::vector<double>& out) {
        const double* parameters[] = {point.data()};
        double* blocks[] = {jacobian.data()};
        if (!cost.Evaluate(parameters, residuals.data(), blocks)) return false;
        out.assign(width, 0.0);
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t c = 0; c < width; ++c) out[c] += jacobian[i * width + c] * residuals[i];
        return true;
    };
    ceres::Solver::Summary summary;
    std::vector<int> held;
    int remaining = controls.maximum_iterations;
    double seconds = controls.maximum_elapsed_time_seconds;
    bool settled = false, all_held = false, exhausted = false;
    for (bool first = true;; first = false) {
        if (!gradient_at(optimizer, gradient)) break;
        std::vector<int> next;
        for (std::size_t c = 0; c < width; ++c)
            if ((optimizer[c] - lower[c] < policy::FIT_ACTIVE_BOUND_DISTANCE.value && gradient[c] > 0.0)
                || (upper[c] - optimizer[c] < policy::FIT_ACTIVE_BOUND_DISTANCE.value && gradient[c] < 0.0))
                next.push_back(static_cast<int>(c));
        // Converged with an unchanged set, or every coordinate held: a KKT point of the bounded problem.
        all_held = next.size() == width;
        settled = (!first && summary.termination_type == ceres::CONVERGENCE && next == held) || all_held;
        exhausted = !settled && (remaining <= 0 || seconds <= 0.0);
        if (settled || exhausted) break;
        held = next;
        ceres::Problem problem(problem_options);
        problem.AddResidualBlock(&cost, nullptr, optimizer.data());
        auto* box = new BoxChart(optimizer, lower, upper, held);  // owned by the problem, outlives the solve
        problem.SetManifold(optimizer.data(), box);
        options.callbacks = {box};
        options.max_num_iterations = remaining;
        options.max_solver_time_in_seconds = seconds;
        ceres::Solve(options, &problem, &summary);
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
    const auto tie = set.tie();
    result.physical = set.physical(optimizer);
    result.parameters = set.apply(request.parameters, optimizer);
    result.predictions = final.predictions;
    result.residuals = final.residuals;
    result.statuses = final.statuses;
    result.physical_jacobian = final.gradient;
    result.weighted_residuals.resize(n);
    result.optimizer_jacobian.assign(n * width, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        result.weighted_residuals[i] = sqrt_weights[i] * final.residuals[i];
        chain(request, tie, i, sqrt_weights[i], final.predictions[i],
              final.gradient.data() + i * tie.size(), result.optimizer_jacobian.data() + i * width);
    }
    if (const auto row = objective(result.weighted_residuals.data(), n, result.final_cost)) {
        evaluator.fail(*row, optimizer, ObservationStatus::Nonfinite, cost_cause(result.weighted_residuals.data(), *row, "final"), true);
        result.status = "failed";
        result.message = result.failures.back().message + ": " + result.message;
    }
    for (std::size_t i = 0; i < n; ++i) {
        if (std::all_of(result.optimizer_jacobian.begin() + static_cast<std::ptrdiff_t>(i * width),
                        result.optimizer_jacobian.begin() + static_cast<std::ptrdiff_t>((i + 1) * width),
                        [](double value) { return std::isfinite(value); })) continue;
        evaluator.fail(i, optimizer, ObservationStatus::Nonfinite, "final optimizer Jacobian not finite", true);
        result.status = "failed";
        result.message = result.failures.back().message + ": " + result.message;
    }
    result.usable = result.status != "failed" && finite(result.physical);
    if (result.usable) covariance(result, request, n);
    return result;
}

}  // namespace epcsaft::regression
