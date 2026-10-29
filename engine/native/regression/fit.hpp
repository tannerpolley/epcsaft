#pragma once

#include "regression/coordinates.hpp"
#include "regression/observations.hpp"

#include <cstddef>
#include <functional>
#include <optional>
#include <string>
#include <vector>

namespace epcsaft::regression {

// Ceres trust-region controls; omitted tolerances select their named policies. The
// iteration and time budgets are caller inputs with no default (ALG-FIT): a fit refuses its
// request until the caller sets both.
struct FitControls final {
    policy::SolverOptions solver_options;
    int maximum_iterations = 0;
    double maximum_elapsed_time_seconds = 0.0;
    int workers = 0;  // processes evaluating state groups, this one included; 0: one per hardware thread
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
// Guarded minimax: physical epigraph/slacks; positive row scales and exact first actions.
struct MinimaxRequest final {
    std::vector<std::size_t> targets;
    std::vector<double> group_scales;
    double q = 0.0, tolerance = 0.0, activity_distance = 0.0, margin = 0.0;
};
struct MinimaxRow final {
    std::string kind;
    std::size_t row = 0, group = 0;
    CheckMetric metric = CheckMetric::AbsoluteDifference;
    double sign = 0.0, threshold = 0.0, scale = 0.0;
    std::optional<std::size_t> auxiliary;
    std::vector<std::size_t> members;
};
struct MinimaxResult final {
    std::string availability = "not_evaluated", application_status, solver_status, linear_solver;
    bool feasible = false, stationary = false;
    double t = 0.0, target_maximum = 0.0, gap = 0.0;
    double primal = 0.0, stationarity = 0.0, dual_sign = 0.0, complementarity = 0.0;
    std::vector<double> point, lower, upper, values, physical_values, jacobian;
    std::vector<double> multipliers, physical_multipliers, z_lower, z_upper, physical_z_lower, physical_z_upper;
    std::vector<double> lower_distances, upper_distances, hold_multipliers, full_stationarity, root_errors;
    std::vector<std::size_t> slack_groups, slack_rows, active_lower, active_upper;
    std::vector<std::string> hold_causes;
    std::vector<MinimaxRow> rows;
    std::vector<CheckResult> frozen_checks;
    std::optional<equilibrium::IpoptDiagnostics> ipopt_diagnostics;
};
enum class MeasuredInputKind { Temperature, Pressure, MoleFraction, Molality };
struct MeasuredInput final {
    MeasuredInputKind kind = MeasuredInputKind::Temperature;
    std::size_t species = 0;  // MoleFraction only
    double sigma = 0.0, lower = 0.0, upper = 0.0;
    std::string unit;
};
struct MeasurementGroup final {
    std::vector<std::size_t> rows;
    std::vector<MeasuredInput> inputs;
    std::vector<std::size_t> dependent_species;
    double feed_total_mol = 0.0;  // required for an equilibrium feed chart
};
struct MeasurementUncertainty final {
    std::vector<MeasurementGroup> groups;
    std::vector<double> output_sigma;
    // At most one, row-major m x m: positive-weight outputs in row order, then active inputs in column order.
    std::optional<std::vector<double>> covariance, correlation;
};
// Under a declared matrix, penalties are squared whitened input components, input_cost is conditional on the
// outputs and jacobian is whitened; corrections stay standardized physical corrections.
struct InputReconciliation final {
    std::vector<double> values, corrections, penalties, jacobian;
    std::vector<std::size_t> row_groups, column_groups, column_inputs;
    double output_cost = 0.0, input_cost = 0.0;
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
    std::optional<MinimaxRequest> minimax;
    std::optional<MeasurementUncertainty> uncertainty;
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
    std::vector<ResidualEnergyAction> parameter_actions{};
    std::optional<core::Convergence> first_numerical_failure;
};

// (J^T J)^-1 of the weighted (whitened, under a declared matrix) optimizer Jacobian, mapped to physical rows by the evaluated D.
// Available only at an interior point of full rank; otherwise the cause is the status.
struct FitCovariance final {
    std::string status = "not_evaluated";  // available | not_converged | rank_deficient | active_bound | unavailable | minimax_unavailable
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
    policy::SolverOptions solver_options;
    // Ceres: converged | limit_reached | initial_evaluation_failed | failed | invalid_setup;
    // minimax: converged | limit_or_failure | point_evaluated | invalid_setup.
    std::string status = "invalid_setup";
    std::string termination = "not_run";
    std::string message;
    bool usable = false;  // a finite candidate point with values was published; status says how it terminated
    std::vector<double> optimizer;
    std::vector<double> physical;  // one value per physical row of the coordinate set
    Parameters parameters;
    std::vector<std::optional<Prediction>> saturation_predictions;
    std::vector<double> predictions;
    std::vector<double> residuals;           // declared form, unweighted
    std::vector<double> weighted_residuals;  // sqrt(w) r, the optimizer residual; whitened under a declared matrix
    std::vector<double> physical_jacobian;   // n x R: d(value)/d(theta)
    std::vector<double> optimizer_jacobian;  // n x U: Ceres weighted residual columns; minimax raw value columns
    std::vector<ObservationStatus> statuses;  // NotEvaluated for skipped rows or rows after a failed prediction
    std::vector<TrialFailure> failures;
    std::size_t iterations = 0;
    std::size_t residual_evaluations = 0;
    std::size_t jacobian_evaluations = 0;
    std::size_t trial_failures = 0;
    std::vector<double> iteration_seconds, iteration_costs;  // per Ceres iteration, 0 = the start
    // Ceres: 0.5 sum(w r^2); minimax: unavailable (NaN). Failed Ceres evaluations are zero, never a
    // partial or infinite sum, and only status and usable say whether a zero is a cost.
    double initial_cost = 0.0;
    double final_cost = 0.0;
    FitCovariance covariance;
    std::vector<core::NumericalRank> rank_records;
    std::optional<core::Convergence> first_numerical_failure;
    EvaluatedCoordinates chart;
    std::vector<CheckResult> checks;
    std::string acceptance = "not_evaluated";
    std::optional<MinimaxResult> minimax;
    std::optional<InputReconciliation> reconciliation;
};

struct ObservationEvaluation final {
    policy::SolverOptions solver_options;
    std::vector<std::optional<Prediction>> saturation_predictions;
    std::vector<double> predictions, residuals, gradient, weighted, optimizer_gradient;
    // Public point columns: optimizer_gradient is dy/dv; this field is dr/dv at unit weight (with uncertainty, the
    // complete standardized and, under a declared matrix, whitened objective columns).
    std::vector<double> residual_optimizer_gradient;
    std::vector<ObservationStatus> statuses;
    double cost = 0.0;
    std::optional<core::Convergence> first_numerical_failure;
    std::size_t first_numerical_failure_observation = 0;
    EvaluatedCoordinates chart;
    std::optional<InputReconciliation> reconciliation;
    std::vector<double> input_gradient;
    std::vector<std::vector<double>> state_densities, reference_values;
};
class ObservationEvaluator final {
public:
    struct Impl;
    ObservationEvaluator(const FitRequest& request, FitResult& failures, const std::vector<std::size_t>& rows);
    bool evaluate(const std::vector<double>& optimizer, bool columns, ObservationEvaluation& out);
private:
    std::shared_ptr<Impl> impl_;
};
std::vector<CheckResult> evaluate_checks(const FitRequest& request, const Parameters& stored, FitResult* final = nullptr,
    const std::vector<Observation>* reconciled = nullptr);
std::vector<Observation> observations_at(const FitRequest& request, const std::vector<double>& optimizer,
    const std::vector<double>& errors);  // signed native differences or decades, for manufactured recovery
// The Python host's worker command (empty: a fit evaluates every row itself) and cancellation poll; a true
// poll stops a least-squares fit by throwing FitInterrupted (#380).
struct FitHost final {
    std::vector<std::string> worker;
    bool (*interrupted)() = nullptr;
};
void set_fit_host(FitHost host, const void* extension);
struct __attribute__((visibility("default"))) FitInterrupted final {};
[[noreturn]] void serve_fit(int fd);  // a worker process: its rows of every evaluation, over fd
FitResult fit(const FitRequest& request);
std::vector<FitResult> fit_sequence(FitRequest request, const std::vector<std::vector<double>>& starts,
    int iterations, double seconds);  // shared budgets, ordered results; unstarted points are absent

// Shared unweighted check algebra and final publication, used by both native estimators.
std::pair<double, double> check_error(CheckMetric metric, const Observation& row, double prediction);
void publish_candidate(FitResult& result, ObservationEvaluation evaluation, const std::vector<double>& point);
void finish_checks(FitResult& result, const FitRequest& request, const std::vector<Observation>* reconciled = nullptr);
using ObservationFunction = std::function<bool(const std::vector<double>&, bool, ObservationEvaluation&)>;
FitResult solve_minimax(const FitRequest& request, const std::vector<double>& start,
    const std::vector<CheckResult>& frozen, const ObservationFunction& evaluate);
FitResult fit_minimax(const FitRequest& request, const std::optional<std::vector<double>>& point = {});
void validate_fit_request(const FitRequest& request, bool optimizing = true);

}  // namespace epcsaft::regression
