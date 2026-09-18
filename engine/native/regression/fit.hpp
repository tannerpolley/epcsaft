#pragma once

#include "regression/coordinates.hpp"
#include "regression/observations.hpp"

#include <cstddef>
#include <optional>
#include <string>
#include <vector>

namespace epcsaft::regression {

// Ceres trust-region controls; the named policy entries are the only defaults.
struct FitControls final {
    double function_tolerance = 0.0;
    double gradient_tolerance = 0.0;
    double parameter_tolerance = 0.0;
    int maximum_iterations = 100;
    double maximum_elapsed_time_seconds = 600.0;
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
};

// One rejected trial: the observation that failed, its typed cause and the optimizer point.
struct TrialFailure final {
    std::size_t trial = 0;
    std::size_t observation = 0;
    std::vector<double> optimizer;
    ObservationStatus status = ObservationStatus::Nonfinite;
    std::string message;
    bool jacobian = false;
};

// (J^T J)^-1 of the weighted optimizer Jacobian, mapped to physical rows by D = S diag(scale).
// Available only at an interior point of full rank; otherwise the cause is the status.
struct FitCovariance final {
    std::string status = "not_evaluated";  // available | rank_deficient | active_bound | unavailable
    std::string assumption;
    std::size_t rank = 0;
    std::vector<double> singular_values;
    std::vector<std::size_t> active_bounds;  // optimizer coordinates at a bound
    std::vector<double> optimizer;  // U x U row-major
    std::vector<double> physical;   // R x R row-major
    std::optional<double> variance_factor;  // sum(w r^2) / (n - rank) when n > rank
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
    std::vector<double> optimizer_jacobian;  // n x U: sqrt(w) (dr/dy) (dy/dtheta) S diag(scale)
    std::vector<ObservationStatus> statuses;
    std::vector<TrialFailure> failures;
    std::size_t iterations = 0;
    std::size_t residual_evaluations = 0;
    std::size_t jacobian_evaluations = 0;
    std::size_t trial_failures = 0;
    double initial_cost = 0.0;  // 0.5 sum(w r^2)
    double final_cost = 0.0;
    FitCovariance covariance;
};

FitResult fit(const FitRequest& request);

}  // namespace epcsaft::regression
