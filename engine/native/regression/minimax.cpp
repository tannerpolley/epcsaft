#include "regression/fit.hpp"
#include "equilibrium/ipopt.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <set>
#include <stdexcept>

namespace epcsaft::regression {
namespace {
using namespace Ipopt;
void require(bool ok, const char* cause) { if (!ok) throw std::invalid_argument(cause); }
bool positive(double x) { return std::isfinite(x) && x > 0.0; }
class Minimax final : public TNLP {
public:
    Minimax(const FitRequest& request, FitResult& result, const ObservationFunction& evaluate,
            const std::vector<double>& start, const std::vector<CheckResult>& frozen)
        : request_(request), result_(result), evaluate_(evaluate), spec_(*request.minimax), u_(start.size()) {
        result_.minimax.emplace();
        auto& out = *result_.minimax;
        out.frozen_checks = frozen;
        require(u_ > 0 && start.size() == request.coordinates.coordinates.size(), "minimax start width");
        require(!spec_.targets.empty() && positive(spec_.q) && positive(spec_.tolerance)
            && positive(spec_.activity_distance) && positive(spec_.margin), "minimax controls");
        require(frozen.size() == request.checks.size() && spec_.group_scales.size() == frozen.size(), "minimax groups/scales");
        std::tie(out.lower, out.upper) = request.coordinates.optimizer_bounds();
        out.point = start;
        ObservationEvaluation initial;
        require(evaluate_(start, true, initial), "minimax initial evaluation unavailable");
        double t = 0.0;
        std::set<std::size_t> targets;
        for (const auto r : spec_.targets) {
            require(targets.insert(r).second && r < request.observations.size(), "minimax target index");
            const auto e = check_error(CheckMetric::AbsoluteLog10Ratio, request.observations[r], initial.predictions.at(r)).first;
            require(std::isfinite(e), "minimax nonpositive target");
            t = std::max(t, std::abs(e));
            paired("target", r, request.checks.size(), CheckMetric::AbsoluteLog10Ratio, 0.0, spec_.q, u_);
        }
        out.point.push_back(t + spec_.margin * spec_.q);
        variable_scales_ = std::vector<double>(u_, 1.0);
        variable_scales_.push_back(spec_.q);
        for (std::size_t g = 0; g < frozen.size(); ++g) {
            const auto& group = request.checks[g];
            if (group.report_only) continue;
            require(frozen[g].status != "unavailable" && positive(spec_.group_scales[g]), "minimax guard unavailable/scale");
            std::vector<std::size_t> members;
            for (const auto r : group.rows) {
                if (frozen[g].maximum_threshold)
                    paired("maximum", r, g, group.metric, *frozen[g].maximum_threshold, spec_.group_scales[g], {});
                if (frozen[g].mean_threshold) {
                    const auto s = out.point.size();
                    const auto e = check_error(group.metric, request.observations.at(r), initial.predictions.at(r)).first;
                    require(std::isfinite(e), "minimax initial guard error");
                    out.point.push_back(std::abs(e) + spec_.margin * spec_.group_scales[g]);
                    variable_scales_.push_back(spec_.group_scales[g]);
                    out.slack_groups.push_back(g);
                    out.slack_rows.push_back(r);
                    members.push_back(s);
                    paired("link", r, g, group.metric, 0.0, spec_.group_scales[g], s);
                }
            }
            if (!members.empty())
                out.rows.push_back({"mean", 0, g, group.metric, 0.0, *frozen[g].mean_threshold, spec_.group_scales[g], {}, members});
        }
        out.lower.resize(out.point.size(), 0.0);
        out.upper.resize(out.point.size(), INFINITY);
    }
    void paired(const char* kind, std::size_t r, std::size_t g, CheckMetric metric,
                double threshold, double scale, std::optional<std::size_t> auxiliary) {
        for (double sign : {1.0, -1.0})
            result_.minimax->rows.push_back({kind, r, g, metric, sign, threshold, scale, auxiliary, {}});
    }
    bool assemble(const Number* x, bool columns) {
        auto& out = *result_.minimax;
        ObservationEvaluation evaluation;
        if (!evaluate_(std::vector<double>(x, x + u_), columns, evaluation)) return false;
        out.values.resize(out.rows.size());
        if (columns) out.jacobian.assign(out.rows.size() * out.point.size(), 0.0);
        for (std::size_t r = 0; r < out.rows.size(); ++r) {
            const auto& row = out.rows[r];
            double value = -row.threshold;
            if (row.members.empty()) {
                const auto [error, slope] = check_error(row.metric, request_.observations.at(row.row), evaluation.predictions.at(row.row));
                if (!std::isfinite(error) || !std::isfinite(slope)) return false;
                value += row.sign * error;
                for (std::size_t c = 0; columns && c < u_; ++c)
                    out.jacobian[r * out.point.size() + c] = row.sign * slope
                        * evaluation.optimizer_gradient.at(row.row * u_ + c) / row.scale;
                if (row.auxiliary) {
                    value -= x[*row.auxiliary];
                    if (columns) out.jacobian[r * out.point.size() + *row.auxiliary] = -1.0 / row.scale;
                }
            } else for (const auto s : row.members) {
                value += x[s] / row.members.size();
                if (columns) out.jacobian[r * out.point.size() + s] = 1.0 / (row.members.size() * row.scale);
            }
            out.values[r] = value / row.scale;
            if (!std::isfinite(out.values[r])) return false;
        }
        evaluation_ = std::move(evaluation);
        return true;
    }
    bool get_nlp_info(Index& n, Index& m, Index& nj, Index& nh, IndexStyleEnum& style) override {
        n = result_.minimax->point.size();
        m = result_.minimax->rows.size();
        nj = n * m;
        nh = 0;
        style = C_STYLE;
        return true;
    }
    bool get_bounds_info(Index n, Number* xl, Number* xu, Index m, Number* gl, Number* gu) override {
        const auto& out = *result_.minimax;
        for (Index c = 0; c < n; ++c) {
            xl[c] = std::isfinite(out.lower[c]) ? out.lower[c] : -1e19;
            xu[c] = std::isfinite(out.upper[c]) ? out.upper[c] : 1e19;
        }
        std::fill(gl, gl + m, -1e19);
        std::fill(gu, gu + m, 0.0);
        return true;
    }
    bool get_starting_point(Index, bool ix, Number* x, bool iz, Number*, Number*, Index, bool il, Number*) override {
        if (!ix || iz || il) return false;
        const auto& point = result_.minimax->point;
        std::copy(point.begin(), point.end(), x);
        return true;
    }
    bool eval_f(Index, const Number* x, bool, Number& f) override { f = x[u_] / spec_.q; return std::isfinite(f); }
    bool eval_grad_f(Index n, const Number*, bool, Number* gradient) override {
        std::fill(gradient, gradient + n, 0.0);
        gradient[u_] = 1.0 / spec_.q;
        return true;
    }
    bool eval_g(Index, const Number* x, bool, Index, Number* values) override {
        if (!assemble(x, false)) return false;
        const auto& v = result_.minimax->values;
        std::copy(v.begin(), v.end(), values);
        return true;
    }
    bool eval_jac_g(Index n, const Number* x, bool, Index m, Index, Index* ir, Index* jc, Number* data) override {
        if (!data) {
            for (Index r = 0; r < m; ++r) for (Index c = 0; c < n; ++c) {
                ir[r * n + c] = r;
                jc[r * n + c] = c;
            }
            return true;
        }
        if (!assemble(x, true)) return false;
        const auto& j = result_.minimax->jacobian;
        std::copy(j.begin(), j.end(), data);
        return true;
    }
    bool eval_h(Index, const Number*, bool, Number, Index, const Number*, bool, Index entries, Index*, Index*, Number*) override {
        return entries == 0;
    }
    bool intermediate_callback(AlgorithmMode, Index iteration, Number, Number, Number, Number, Number,
                               Number, Number, Number, Index, const IpoptData*, IpoptCalculatedQuantities*) override {
        result_.iterations = iteration;
        return std::chrono::duration<double>(std::chrono::steady_clock::now() - begin_).count()
            < request_.controls.maximum_elapsed_time_seconds;
    }
    void finalize_solution(SolverReturn status, Index n, const Number* x, const Number* zl, const Number* zu,
                           Index m, const Number*, const Number* lambda, Number, const IpoptData*, IpoptCalculatedQuantities*) override {
        auto& out = *result_.minimax;
        out.solver_status = equilibrium::solver_status_name(status);
        out.point.assign(x, x + n);
        out.z_lower.assign(zl, zl + n);
        out.z_upper.assign(zu, zu + n);
        out.multipliers.assign(lambda, lambda + m);
        finalized_ = true;
    }
    void publish() {
        const auto& out = *result_.minimax;
        publish_candidate(result_, std::move(evaluation_), std::vector<double>(out.point.begin(), out.point.begin() + u_));
        result_.initial_cost = result_.final_cost = NAN;
        result_.weighted_residuals.clear(); result_.covariance.status = "minimax_unavailable";
        result_.usable = true;
    }
    void point_values() {
        auto& out = *result_.minimax;
        require(assemble(out.point.data(), true), "minimax point unavailable");
        publish();
        result_.status = "point_evaluated";
        out.availability = "values_only";
    }
    void finish() {
        auto& out = *result_.minimax;
        if (!finalized_ || !assemble(out.point.data(), true)) { out.availability = "point_or_columns_unavailable"; return; }
        for (const auto* operands : {&out.multipliers, &out.z_lower, &out.z_upper, &out.jacobian})
            require(std::all_of(operands->begin(), operands->end(), [](double x) { return std::isfinite(x); }), "minimax nonfinite KKT operands");
        const std::size_t n = out.point.size(), physical = evaluation_.chart.physical.size();
        std::vector<double> stationarity(n, 0.0), full(physical, 0.0);
        stationarity[u_] = 1.0 / spec_.q;
        out.t = out.point[u_];
        out.availability = "available";
        for (const auto r : spec_.targets)
            out.target_maximum = std::max(out.target_maximum, std::abs(check_error(CheckMetric::AbsoluteLog10Ratio,
                request_.observations.at(r), evaluation_.predictions.at(r)).first));
        out.gap = out.t - out.target_maximum;
        for (std::size_t r = 0; r < out.rows.size(); ++r) {
            const auto& row = out.rows[r];
            const double lambda = out.multipliers[r], multiplier = spec_.q * lambda / row.scale;
            out.physical_values.push_back(row.scale * out.values[r]);
            out.physical_multipliers.push_back(multiplier);
            out.primal = std::max(out.primal, out.values[r]); out.dual_sign = std::max(out.dual_sign, -lambda);
            out.complementarity = std::max(out.complementarity, std::abs(lambda * out.values[r]));
            for (std::size_t c = 0; c < n; ++c) stationarity[c] += lambda * out.jacobian[r * n + c];
            double root_error = 0.0;
            if (row.members.empty()) {
                const double slope = check_error(row.metric, request_.observations.at(row.row), evaluation_.predictions.at(row.row)).second;
                for (std::size_t p = 0; p < physical; ++p)
                    full[p] += multiplier * row.sign * slope * evaluation_.gradient.at(row.row * physical + p);
                for (std::size_t h = 0; h < evaluation_.chart.holds.size(); ++h) {
                    const auto& hold = evaluation_.chart.holds[h];
                    root_error += std::abs(slope * evaluation_.gradient.at(row.row * physical
                        + request_.coordinates.active_parameters().size() + h) / row.scale)
                        * hold.uncertainty;
                }
            }
            out.root_errors.push_back(root_error);
        }
        for (std::size_t c = 0; c < n; ++c) {
            const double zl = out.z_lower[c], zu = out.z_upper[c];
            out.physical_z_lower.push_back(spec_.q * zl);
            out.physical_z_upper.push_back(spec_.q * zu);
            const double dl = out.point[c] - out.lower[c], du = out.upper[c] - out.point[c];
            out.lower_distances.push_back(dl / variable_scales_[c]); out.upper_distances.push_back(du / variable_scales_[c]);
            if (std::isfinite(dl)) {
                out.primal = std::max(out.primal, -dl / variable_scales_[c]);
                out.complementarity = std::max(out.complementarity, std::abs(zl * dl));
                if (dl / variable_scales_[c] <= spec_.activity_distance) out.active_lower.push_back(c);
            }
            if (std::isfinite(du)) {
                out.primal = std::max(out.primal, -du / variable_scales_[c]);
                out.complementarity = std::max(out.complementarity, std::abs(zu * du));
                if (du / variable_scales_[c] <= spec_.activity_distance) out.active_upper.push_back(c);
            }
            out.dual_sign = std::max({out.dual_sign, -zl, -zu});
            out.stationarity = std::max(out.stationarity, std::abs(stationarity[c] - zl + zu));
        }
        const auto independent = request_.coordinates.active_parameters().size();
        for (std::size_t h = 0; h < evaluation_.chart.holds.size(); ++h) {
            const auto& hold = evaluation_.chart.holds[h];
            const auto& declaration = request_.holds.at(h);
            const double diameter = evaluation_.chart.physical.at(independent + h);
            const bool regular = std::isfinite(hold.derivative) && hold.derivative != 0.0
                && diameter > declaration.lower && diameter < declaration.upper;
            const double nu = regular ? -full.at(independent + h) / hold.derivative : NAN;
            out.hold_multipliers.push_back(nu); out.hold_causes.push_back(regular ? "available" : "singular_or_bracket_boundary");
            for (std::size_t p = 0; p < physical; ++p) full[p] += nu * hold.parameter_derivative.at(p);
        }
        for (std::size_t c = 0; c < u_; ++c) {
            double value = spec_.q * (-out.z_lower[c] + out.z_upper[c]);
            for (std::size_t p = 0; p < physical; ++p) value += full[p] * evaluation_.chart.map.at(p * u_ + c);
            out.full_stationarity.push_back(value);
        }
        for (std::size_t h = 0; h < evaluation_.chart.holds.size(); ++h) out.full_stationarity.push_back(full[independent + h]);
        out.feasible = out.primal <= spec_.tolerance;
        out.stationary = std::max({out.stationarity, out.dual_sign, out.complementarity}) <= spec_.tolerance;
        publish();
        result_.status = out.solver_status == "SUCCESS" ? "converged" : "limit_or_failure";
        result_.termination = out.solver_status;
        result_.message = "local constrained estimate; global feasibility and optimality unresolved";
    }
    bool finalized_ = false;
    std::chrono::steady_clock::time_point begin_ = std::chrono::steady_clock::now();
private:
    const FitRequest& request_;
    FitResult& result_;
    const ObservationFunction& evaluate_;
    const MinimaxRequest& spec_;
    std::size_t u_;
    std::vector<double> variable_scales_;
    ObservationEvaluation evaluation_;
};
}  // namespace

FitResult solve_minimax(const FitRequest& request, const std::vector<double>& start,
    const std::vector<CheckResult>& frozen, const ObservationFunction& evaluate) {
    FitResult result; result.solver_options = request.controls.solver_options.resolved({});
    try {
        require(request.minimax.has_value(), "missing minimax declaration");
        SmartPtr<Minimax> nlp = new Minimax(request, result, evaluate, start, frozen);
        const auto [initialized, returned] = equilibrium::run_ipopt(nlp, request.minimax->tolerance,
            request.controls.maximum_iterations, request.controls.maximum_elapsed_time_seconds, false, 1.0, nullptr, 0.0, {},
            &result.minimax->linear_solver);
        result.minimax->application_status = equilibrium::application_return_name(returned);
        if (initialized != Solve_Succeeded || (returned <= Invalid_Problem_Definition && returned != Invalid_Number_Detected)) {
            result.message = "Ipopt " + equilibrium::application_return_name(initialized) + "; " + result.minimax->application_status;
            result.minimax->availability = "fatal_initialization_or_application";
        } else nlp->finish();
    } catch (const std::exception& error) { result.message = error.what(); }
    return result;
}
FitResult fit_minimax(const FitRequest& request, const std::optional<std::vector<double>>& evaluation_point) {
    FitResult failures; failures.solver_options = request.controls.solver_options.resolved({});
    try {
        for (const auto& coordinate : request.coordinates.coordinates)
            require(!coordinate.reaction || coordinate.reaction->coefficient != ReactionCoefficient::EnthalpyShift,
                "reaction enthalpy shift is not supported in minimax");
        validate_fit_request(request, !evaluation_point);
        const auto frozen = evaluate_checks(request, request.baseline.value_or(request.parameters));
        std::vector<std::size_t> selected = request.minimax->targets;
        for (const auto& group : request.checks) if (!group.report_only) selected.insert(selected.end(), group.rows.begin(), group.rows.end());
        ObservationEvaluator evaluator(request, failures, selected);
        const auto point = request.coordinates.start(request.parameters, request.observations, request.thermochemistry);
        const ObservationFunction callback = [&](const auto& u, bool columns, auto& out) { return evaluator.evaluate(u, columns, out); };
        FitResult result; result.solver_options = request.controls.solver_options.resolved({});
        if (evaluation_point) {
            Minimax nlp(request, result, callback, *evaluation_point, frozen);
            nlp.point_values();
        } else result = solve_minimax(request, point, frozen, callback);
        result.failures = std::move(failures.failures); result.trial_failures = failures.trial_failures;
        result.residual_evaluations = failures.residual_evaluations; result.jacobian_evaluations = failures.jacobian_evaluations;
        if (result.usable || evaluation_point) finish_checks(result, request);
        return result;
    } catch (const std::exception& error) { failures.message = error.what(); return failures; }
}
}  // namespace epcsaft::regression
