#pragma once

#include "regression/coordinates.hpp"
#include "regression/observations.hpp"

#include <cstddef>
#include <optional>
#include <string>
#include <vector>

namespace epcsaft::regression {

// Ceres trust-region controls; a zero tolerance selects its named policy default. The
// iteration and time budgets are explicit caller controls, not tolerances.
struct FitControls final {
    double function_tolerance = 0.0;
    double gradient_tolerance = 0.0;
    double parameter_tolerance = 0.0;
    int maximum_iterations = 100;
    double maximum_elapsed_time_seconds = 600.0;
};

enum class CheckMetric { AbsoluteDifference, AbsoluteLog10Ratio };
struct CheckGroup final {
    std::string name;
    std::vector<std::size_t> rows;
    CheckMetric metric = CheckMetric::AbsoluteDifference;
    std::optional<double> maximum, mean;
    bool baseline = false, report_only = false;
};
struct CheckResult final {
    std::string name, unit, status;
    std::vector<std::size_t> rows, failed_rows;
    std::vector<double> errors, slopes;
    double maximum = 0.0, mean = 0.0;
    std::size_t argmax = 0;
    std::optional<double> maximum_threshold, mean_threshold, baseline_maximum, baseline_mean;
    bool report_only = false;
};
struct FitRequest final {
    Parameters parameters;
    std::optional<ThermochemistryRecord> thermochemistry;
    CoordinateSet coordinates;
    std::vector<Observation> observations;
    // Relative weights w_i >= 0 (empty: all one). The optimizer residual is sqrt(w_i) r_i with
    // r_i in the declared form; an absolute uncertainty enters only through a Difference scale.
    std::vector<double> weights;
    FitControls controls;
    std::vector<HydrationHold> holds;
    std::vector<CheckGroup> checks;
    std::optional<Parameters> baseline;
    bool singular_vectors = false;
};

// One rejected trial: the observation that failed, its typed cause and the optimizer point. A
// trial the model or an owner refuses by exception (InadmissibleTrial) carries the observation
// count as its index and the verbatim message, so an internal error stays visible as its cause.
struct TrialFailure final {
    std::size_t trial = 0;
    std::size_t observation = 0;
    std::vector<double> optimizer;
    ObservationStatus status = ObservationStatus::Nonfinite;
    std::string message;
    bool jacobian = false;
};

// (J^T J)^-1 of the weighted optimizer Jacobian, mapped to physical rows by the evaluated D.
// Available only at an interior point of full rank; otherwise the cause is the status.
struct FitCovariance final {
    std::string status = "not_evaluated";  // available | not_converged | rank_deficient | active_bound | unavailable
    std::string assumption;
    std::size_t rank = 0;
    std::vector<double> singular_values, right_vectors, left_vectors;
    double projected_gradient_norm = 0.0;
    std::vector<std::size_t> objective_rows;  // left-vector row labels, positive-weight observations
    std::vector<std::size_t> active_bounds;  // optimizer coordinates at a bound
    std::vector<double> optimizer;  // U x U, symmetric
    std::vector<double> physical;   // R x R, symmetric
    std::optional<double> variance_factor;  // sum(w r^2) / (n_w - rank) over the n_w rows with w > 0
};

struct FitResult final {
    std::string status = "invalid_setup";  // converged | limit_reached | initial_evaluation_failed | failed
    std::string termination = "not_run";
    std::string message;
    bool usable = false;
    std::vector<double> optimizer;
    std::vector<double> physical;  // one value per physical row of the coordinate set
    Parameters parameters;
    std::vector<double> predictions;
    std::vector<double> residuals;           // declared form, unweighted
    std::vector<double> weighted_residuals;  // sqrt(w) r, the optimizer residual
    std::vector<double> physical_jacobian;   // n x R: d(value)/d(theta)
    std::vector<double> optimizer_jacobian;  // n x U: sqrt(w) (dr/dy) (dy/dtheta) D
    std::vector<ObservationStatus> statuses;  // NotEvaluated for skipped rows or rows after a failed prediction
    std::vector<TrialFailure> failures;
    std::size_t iterations = 0;
    std::size_t residual_evaluations = 0;
    std::size_t jacobian_evaluations = 0;
    std::size_t trial_failures = 0;
    std::vector<double> iteration_seconds, iteration_costs;  // per Ceres iteration, 0 = the start
    // 0.5 sum(w r^2), accumulated in row order; on a failed evaluation both are exactly zero, never a
    // partial or infinite sum, and only status and usable say whether a zero is a cost.
    double initial_cost = 0.0;
    double final_cost = 0.0;
    FitCovariance covariance;
    EvaluatedCoordinates chart;
    std::vector<CheckResult> checks;
    std::string acceptance = "not_evaluated";
};

struct ObservationEvaluation final {
    std::vector<double> predictions, residuals, gradient, weighted, optimizer_gradient;
    std::vector<ObservationStatus> statuses;
    double cost = 0.0;
    EvaluatedCoordinates chart;
};
class ObservationEvaluator final {
public:
    struct Impl;
    ObservationEvaluator(const FitRequest& request, FitResult& failures, const std::vector<std::size_t>& rows);
    bool evaluate(const std::vector<double>& optimizer, bool columns, ObservationEvaluation& out);
private:
    std::shared_ptr<Impl> impl_;
};
std::vector<CheckResult> evaluate_checks(const FitRequest& request, const Parameters& stored, FitResult* final = nullptr);
std::vector<Observation> observations_at(const FitRequest& request, const std::vector<double>& optimizer,
    const std::vector<double>& errors);  // signed native differences or decades, for manufactured recovery
FitResult fit(const FitRequest& request);
std::vector<FitResult> fit_sequence(FitRequest request, const std::vector<std::vector<double>>& starts,
    int iterations, double seconds);  // shared budgets, ordered results; unstarted points are absent

}  // namespace epcsaft::regression
