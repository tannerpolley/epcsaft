#pragma once

#include "owner_evaluator.hpp"

#include <cstddef>
#include <string>
#include <vector>

namespace epcsaft_engine::regression {

struct EvaluatorHandle final {
  interfaces::OwnerEvaluatorHandle evaluator;
};

struct SolverControls final {
  double function_tolerance;
  double gradient_tolerance;
  double parameter_tolerance;
  int maximum_iterations;
  double maximum_elapsed_time_seconds;
  std::size_t maximum_uncached_owner_calls;
};

struct LeastSquaresRequest final {
  std::vector<EvaluatorHandle> evaluators;
  std::vector<double> observations;
  std::vector<double> residual_scales;
  std::string residual_scale_metadata;
  std::vector<double> initial_parameters;
  std::vector<double> lower_bounds;
  std::vector<double> upper_bounds;
  std::vector<double> affine_origins;
  std::vector<double> affine_scales;
  SolverControls controls;
};

struct EvaluatorRows final {
  std::string evaluator_fingerprint;
  std::size_t evaluator_index;
  std::size_t row_begin;
  std::size_t row_end;
};

struct FailureEvent final {
  std::size_t evaluator_index;
  std::size_t row_begin;
  std::size_t row_end;
  std::string evaluation_phase;
  std::vector<double> physical_parameters;
  std::string code;
  std::string reason;
  std::size_t count;
};

struct CeresOptionsRecord final {
  std::string minimizer_type;
  std::string trust_region_strategy_type;
  std::string linear_solver_type;
  std::string dense_linear_algebra_library_type;
  bool jacobi_scaling;
  bool use_inner_iterations;
  bool use_nonmonotonic_steps;
  int maximum_line_search_step_size_iterations;
  int maximum_consecutive_invalid_steps;
  int number_of_threads;
  double initial_trust_region_radius;
  double maximum_trust_region_radius;
  double minimum_trust_region_radius;
  double minimum_relative_decrease;
  double minimum_lm_diagonal;
  double maximum_lm_diagonal;
  bool use_mixed_precision_solves;
  int maximum_refinement_iterations;
  double function_tolerance;
  double gradient_tolerance;
  double parameter_tolerance;
  int maximum_iterations;
  double maximum_elapsed_time_seconds;
};

struct SolverCounters final {
  std::size_t owner_call_limit{0};
  std::size_t uncached_residual_owner_attempts{0};
  std::size_t uncached_jacobian_owner_attempts{0};
  std::size_t owner_cache_hits{0};
  std::size_t typed_trial_failures{0};
  int ceres_residual_evaluations{0};
  int ceres_jacobian_evaluations{0};
  int accepted_steps{0};
  int rejected_steps{0};
  int iterations{0};
  bool iteration_limit_reached{false};
  bool elapsed_time_limit_reached{false};
  bool owner_call_limit_reached{false};
};

struct CovarianceResult final {
  std::string status{"not_evaluated"};
  std::string algorithm{"dense_svd"};
  int null_space_rank{-1};
  bool apply_loss_function{true};
  int number_of_threads{1};
  double minimum_reciprocal_condition_number{1.0e-14};
  std::vector<double> optimizer_coordinates;
  std::vector<double> physical_coordinates;
  bool active_bound{false};
};

struct LeastSquaresResult final {
  std::string status{"invalid_setup"};
  std::string termination{"not_run"};
  std::string message;
  bool usable{false};
  std::string parameter_contract_fingerprint;
  std::string residual_scale_metadata;
  std::vector<EvaluatorRows> evaluator_rows;
  std::vector<double> parameters;
  std::vector<double> predictions;
  std::vector<double> residuals;
  std::vector<double> optimizer_jacobian;
  std::vector<FailureEvent> failures;
  CeresOptionsRecord requested_options;
  CeresOptionsRecord resolved_options;
  SolverCounters counters;
  CovarianceResult covariance;
};

[[nodiscard]] LeastSquaresResult solve(const LeastSquaresRequest &request);

} // namespace epcsaft_engine::regression
