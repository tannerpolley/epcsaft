#include "ceres_least_squares.hpp"

#include <ceres/ceres.h>
#include <ceres/covariance.h>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <string>
#include <utility>
#include <variant>
#include <vector>

namespace epcsaft_engine::regression {
namespace {

constexpr std::size_t kMaximumFailureCodeLength = 128;
constexpr std::size_t kMaximumFailureReasonLength = 512;

bool finite(const std::vector<double> &values) {
  return std::all_of(values.begin(), values.end(),
                     [](double value) { return std::isfinite(value); });
}

bool finite_or_infinite(const std::vector<double> &values) {
  return std::all_of(values.begin(), values.end(),
                     [](double value) { return !std::isnan(value); });
}

std::string bounded(std::string value, std::size_t maximum) {
  if (value.size() > maximum) {
    value.resize(maximum);
  }
  return value;
}

interfaces::EvaluationFailure failure(std::string code, std::string reason) {
  return {bounded(std::move(code), kMaximumFailureCodeLength),
          bounded(std::move(reason), kMaximumFailureReasonLength)};
}

CeresOptionsRecord requested_options(const SolverControls &controls) {
  return {
      "trust_region",
      "levenberg_marquardt",
      "dense_qr",
      "eigen",
      false,
      false,
      false,
      0,
      5,
      1,
      1.0e4,
      1.0e16,
      1.0e-32,
      1.0e-3,
      1.0e-6,
      1.0e32,
      false,
      0,
      controls.function_tolerance,
      controls.gradient_tolerance,
      controls.parameter_tolerance,
      controls.maximum_iterations,
      controls.maximum_elapsed_time_seconds,
  };
}

struct EvaluatorCache final {
  bool valid{false};
  bool exact_jacobian{false};
  std::vector<double> parameters;
  interfaces::EvaluationResult result{failure("cache_empty", "cache is empty")};
};

struct CombinedEvaluation final {
  std::vector<double> predictions;
  std::vector<double> jacobian;
};

class EvaluationState final {
public:
  EvaluationState(const LeastSquaresRequest &request,
                  std::vector<EvaluatorRows> rows,
                  std::size_t parameter_count)
      : request_(request), rows_(std::move(rows)),
        parameter_count_(parameter_count), caches_(request.evaluators.size()) {
    counters_.owner_call_limit = request.controls.maximum_uncached_owner_calls;
  }

  std::variant<CombinedEvaluation, interfaces::EvaluationFailure>
  evaluate(const std::vector<double> &parameters, bool exact_jacobian,
           const std::string &phase, bool bypass_cache) {
    CombinedEvaluation combined;
    combined.predictions.reserve(request_.observations.size());
    if (exact_jacobian) {
      combined.jacobian.reserve(request_.observations.size() * parameter_count_);
    }
    for (std::size_t index = 0; index < request_.evaluators.size(); ++index) {
      const auto result = evaluate_one(index, parameters, exact_jacobian,
                                       bypass_cache);
      if (const auto *owner_failure =
              std::get_if<interfaces::EvaluationFailure>(&result)) {
        failures_.push_back({
            index,
            rows_[index].row_begin,
            rows_[index].row_end,
            phase,
            parameters,
            owner_failure->code,
            owner_failure->reason,
            failures_.size() + 1,
        });
        if (phase == "ceres_trial") {
          ++counters_.typed_trial_failures;
          if (exact_jacobian) {
            required_jacobian_failed_ = true;
          }
        }
        return *owner_failure;
      }
      const auto &evaluation = std::get<interfaces::Evaluation>(result);
      combined.predictions.insert(combined.predictions.end(),
                                  evaluation.predictions.begin(),
                                  evaluation.predictions.end());
      if (exact_jacobian) {
        combined.jacobian.insert(combined.jacobian.end(),
                                 evaluation.jacobian.begin(),
                                 evaluation.jacobian.end());
      }
    }
    return combined;
  }

  SolverCounters &counters() noexcept { return counters_; }
  const std::vector<FailureEvent> &failures() const noexcept { return failures_; }
  bool required_jacobian_failed() const noexcept {
    return required_jacobian_failed_;
  }

private:
  interfaces::EvaluationResult evaluate_one(std::size_t index,
                                             const std::vector<double> &parameters,
                                             bool exact_jacobian,
                                             bool bypass_cache) {
    auto &cache = caches_[index];
    if (!bypass_cache && cache.valid &&
        cache.exact_jacobian == exact_jacobian &&
        cache.parameters == parameters) {
      ++counters_.owner_cache_hits;
      return cache.result;
    }

    const std::size_t attempts = counters_.uncached_residual_owner_attempts +
                                 counters_.uncached_jacobian_owner_attempts;
    interfaces::EvaluationResult result;
    if (attempts >= request_.controls.maximum_uncached_owner_calls) {
      counters_.owner_call_limit_reached = true;
      result = failure("owner_call_limit", "uncached owner-call limit reached");
    } else {
      if (exact_jacobian) {
        ++counters_.uncached_jacobian_owner_attempts;
      } else {
        ++counters_.uncached_residual_owner_attempts;
      }
      result = request_.evaluators[index]->evaluate(parameters, exact_jacobian);
      result = validate(index, exact_jacobian, std::move(result));
    }
    cache = {true, exact_jacobian, parameters, result};
    return result;
  }

  interfaces::EvaluationResult
  validate(std::size_t index, bool exact_jacobian,
           interfaces::EvaluationResult result) const {
    if (auto *owner_failure =
            std::get_if<interfaces::EvaluationFailure>(&result)) {
      if (owner_failure->code.empty()) {
        owner_failure->code = "owner_failure";
      }
      if (owner_failure->reason.empty()) {
        owner_failure->reason = "owner returned an empty failure reason";
      }
      owner_failure->code =
          bounded(std::move(owner_failure->code), kMaximumFailureCodeLength);
      owner_failure->reason =
          bounded(std::move(owner_failure->reason), kMaximumFailureReasonLength);
      return result;
    }
    const auto &evaluation = std::get<interfaces::Evaluation>(result);
    const std::size_t outputs =
        request_.evaluators[index]->output_count();
    if (evaluation.predictions.size() != outputs ||
        !finite(evaluation.predictions)) {
      return failure("owner_contract_failure",
                     "owner returned invalid predictions");
    }
    const std::size_t expected = exact_jacobian ? outputs * parameter_count_ : 0;
    if (evaluation.jacobian.size() != expected ||
        !finite(evaluation.jacobian)) {
      return failure("owner_contract_failure",
                     "owner returned an invalid Jacobian");
    }
    return result;
  }

  const LeastSquaresRequest &request_;
  std::vector<EvaluatorRows> rows_;
  std::size_t parameter_count_;
  std::vector<EvaluatorCache> caches_;
  SolverCounters counters_;
  std::vector<FailureEvent> failures_;
  bool required_jacobian_failed_{false};
};

class OwnerCostFunction final : public ceres::CostFunction {
public:
  OwnerCostFunction(EvaluationState &state, const LeastSquaresRequest &request,
                    std::vector<double> residual_scales)
      : state_(state), request_(request),
        residual_scales_(std::move(residual_scales)) {
    set_num_residuals(static_cast<int>(request.observations.size()));
    mutable_parameter_block_sizes()->push_back(
        static_cast<int>(request.initial_parameters.size()));
  }

  bool Evaluate(double const *const *parameters, double *residuals,
                double **jacobians) const override {
    const std::size_t parameter_count = request_.initial_parameters.size();
    std::vector<double> physical(parameter_count);
    for (std::size_t index = 0; index < parameter_count; ++index) {
      physical[index] = request_.affine_origins[index] +
                        request_.affine_scales[index] * parameters[0][index];
    }
    for (std::size_t index = 0; index < parameter_count; ++index) {
      if (!std::isfinite(physical[index]) ||
          physical[index] < request_.lower_bounds[index] ||
          physical[index] > request_.upper_bounds[index]) {
        return false;
      }
    }
    const bool exact_jacobian = jacobians != nullptr && jacobians[0] != nullptr;
    const auto result =
        state_.evaluate(physical, exact_jacobian, phase_, false);
    const auto *evaluation = std::get_if<CombinedEvaluation>(&result);
    if (evaluation == nullptr) {
      return false;
    }
    for (std::size_t row = 0; row < request_.observations.size(); ++row) {
      residuals[row] = (evaluation->predictions[row] -
                        request_.observations[row]) /
                       residual_scales_[row];
      if (!std::isfinite(residuals[row])) {
        return false;
      }
    }
    if (exact_jacobian) {
      for (std::size_t row = 0; row < request_.observations.size(); ++row) {
        for (std::size_t column = 0; column < parameter_count; ++column) {
          jacobians[0][row * parameter_count + column] =
              evaluation->jacobian[row * parameter_count + column] *
              request_.affine_scales[column] / residual_scales_[row];
          if (!std::isfinite(
                  jacobians[0][row * parameter_count + column])) {
            return false;
          }
        }
      }
    }
    return true;
  }

  void set_phase(std::string phase) { phase_ = std::move(phase); }

private:
  EvaluationState &state_;
  const LeastSquaresRequest &request_;
  std::vector<double> residual_scales_;
  std::string phase_{"ceres_trial"};
};

std::string termination_name(ceres::TerminationType termination) {
  switch (termination) {
  case ceres::CONVERGENCE:
    return "convergence";
  case ceres::NO_CONVERGENCE:
    return "no_convergence";
  case ceres::FAILURE:
    return "failure";
  case ceres::USER_SUCCESS:
    return "user_success";
  case ceres::USER_FAILURE:
    return "user_failure";
  }
  return "failure";
}

std::string minimizer_name(ceres::MinimizerType value) {
  return value == ceres::TRUST_REGION ? "trust_region" : "line_search";
}

std::string trust_region_name(ceres::TrustRegionStrategyType value) {
  return value == ceres::LEVENBERG_MARQUARDT ? "levenberg_marquardt"
                                             : "dogleg";
}

std::string linear_solver_name(ceres::LinearSolverType value) {
  return value == ceres::DENSE_QR ? "dense_qr" : "different";
}

std::string dense_library_name(ceres::DenseLinearAlgebraLibraryType value) {
  return value == ceres::EIGEN ? "eigen" : "different";
}

bool same_fixed_options(const ceres::Solver::Summary &summary) {
  return summary.minimizer_type == ceres::TRUST_REGION &&
         summary.trust_region_strategy_type == ceres::LEVENBERG_MARQUARDT &&
         summary.linear_solver_type_given == ceres::DENSE_QR &&
         summary.linear_solver_type_used == ceres::DENSE_QR &&
         summary.dense_linear_algebra_library_type == ceres::EIGEN &&
         summary.num_threads_given == 1 && summary.num_threads_used == 1 &&
         !summary.inner_iterations_given && !summary.inner_iterations_used &&
         !summary.mixed_precision_solves_used;
}

void apply_options(const CeresOptionsRecord &record,
                   ceres::Solver::Options &options) {
  options.minimizer_type = ceres::TRUST_REGION;
  options.trust_region_strategy_type = ceres::LEVENBERG_MARQUARDT;
  options.linear_solver_type = ceres::DENSE_QR;
  options.dense_linear_algebra_library_type = ceres::EIGEN;
  options.jacobi_scaling = record.jacobi_scaling;
  options.use_inner_iterations = record.use_inner_iterations;
  options.use_nonmonotonic_steps = record.use_nonmonotonic_steps;
  options.max_num_line_search_step_size_iterations =
      record.maximum_line_search_step_size_iterations;
  options.max_num_consecutive_invalid_steps =
      record.maximum_consecutive_invalid_steps;
  options.num_threads = record.number_of_threads;
  options.initial_trust_region_radius = record.initial_trust_region_radius;
  options.max_trust_region_radius = record.maximum_trust_region_radius;
  options.min_trust_region_radius = record.minimum_trust_region_radius;
  options.min_relative_decrease = record.minimum_relative_decrease;
  options.min_lm_diagonal = record.minimum_lm_diagonal;
  options.max_lm_diagonal = record.maximum_lm_diagonal;
  options.use_mixed_precision_solves = record.use_mixed_precision_solves;
  options.max_num_refinement_iterations = record.maximum_refinement_iterations;
  options.function_tolerance = record.function_tolerance;
  options.gradient_tolerance = record.gradient_tolerance;
  options.parameter_tolerance = record.parameter_tolerance;
  options.max_num_iterations = record.maximum_iterations;
  options.max_solver_time_in_seconds = record.maximum_elapsed_time_seconds;
}

std::vector<double> expand_residual_scales(const LeastSquaresRequest &request) {
  if (request.residual_scales.empty()) {
    return std::vector<double>(request.observations.size(), 1.0);
  }
  if (request.residual_scales.size() == 1) {
    return std::vector<double>(request.observations.size(),
                               request.residual_scales.front());
  }
  return request.residual_scales;
}

} // namespace

LeastSquaresResult solve(const LeastSquaresRequest &request) {
  LeastSquaresResult result;
  result.requested_options = requested_options(request.controls);
  result.resolved_options = result.requested_options;
  result.counters.owner_call_limit =
      request.controls.maximum_uncached_owner_calls;
  result.residual_scale_metadata = request.residual_scale_metadata;

  const std::size_t parameter_count = request.initial_parameters.size();
  if (request.evaluators.empty() || parameter_count == 0 ||
      request.observations.empty()) {
    result.message = "at least one evaluator, parameter, and observation is required";
    return result;
  }
  if (parameter_count > static_cast<std::size_t>(std::numeric_limits<int>::max()) ||
      request.observations.size() >
          static_cast<std::size_t>(std::numeric_limits<int>::max())) {
    result.message = "Ceres dimensions exceed the supported integer range";
    return result;
  }
  const auto vectors_have_parameter_count = [&](const auto &values) {
    return values.size() == parameter_count && finite(values);
  };
  if (!vectors_have_parameter_count(request.initial_parameters) ||
      request.lower_bounds.size() != parameter_count ||
      request.upper_bounds.size() != parameter_count ||
      !finite_or_infinite(request.lower_bounds) ||
      !finite_or_infinite(request.upper_bounds) ||
      !vectors_have_parameter_count(request.affine_origins) ||
      !vectors_have_parameter_count(request.affine_scales) ||
      !finite(request.observations)) {
    result.message = "parameter, bound, affine, and observation inputs must be finite and dimensionally consistent";
    return result;
  }
  for (std::size_t index = 0; index < parameter_count; ++index) {
    if (request.lower_bounds[index] > request.upper_bounds[index] ||
        request.initial_parameters[index] < request.lower_bounds[index] ||
        request.initial_parameters[index] > request.upper_bounds[index] ||
        !(request.affine_scales[index] > 0.0)) {
      result.message = "bounds, initial parameters, or affine scales are invalid";
      return result;
    }
  }
  const auto &controls = request.controls;
  if (!std::isfinite(controls.function_tolerance) ||
      !std::isfinite(controls.gradient_tolerance) ||
      !std::isfinite(controls.parameter_tolerance) ||
      controls.function_tolerance < 0.0 || controls.gradient_tolerance < 0.0 ||
      controls.parameter_tolerance < 0.0 || controls.maximum_iterations <= 0 ||
      !std::isfinite(controls.maximum_elapsed_time_seconds) ||
      !(controls.maximum_elapsed_time_seconds > 0.0) ||
      controls.maximum_uncached_owner_calls == 0) {
    result.message = "solver tolerances and resource limits are invalid";
    return result;
  }
  if (!request.residual_scales.empty() &&
      request.residual_scales.size() != 1 &&
      request.residual_scales.size() != request.observations.size()) {
    result.message = "residual scales must be absent, scalar, or diagonal";
    return result;
  }
  if (!finite(request.residual_scales) ||
      std::any_of(request.residual_scales.begin(), request.residual_scales.end(),
                  [](double value) { return !(value > 0.0); })) {
    result.message = "residual scales must be finite and positive";
    return result;
  }

  std::size_t output_count = 0;
  std::string parameter_contract_fingerprint;
  result.evaluator_rows.reserve(request.evaluators.size());
  for (std::size_t index = 0; index < request.evaluators.size(); ++index) {
    const auto &evaluator = request.evaluators[index];
    if (!evaluator || evaluator->parameter_count() != parameter_count ||
        evaluator->output_count() == 0 ||
        evaluator->parameter_contract_fingerprint().empty() ||
        evaluator->evaluator_fingerprint().empty()) {
      result.message = "evaluator dimensions and fingerprints are required";
      return result;
    }
    if (index == 0) {
      parameter_contract_fingerprint =
          evaluator->parameter_contract_fingerprint();
    } else if (evaluator->parameter_contract_fingerprint() !=
               parameter_contract_fingerprint) {
      result.message = "evaluator parameter contracts do not match";
      return result;
    }
    const std::size_t next = output_count + evaluator->output_count();
    if (next < output_count) {
      result.message = "evaluator output count overflow";
      return result;
    }
    result.evaluator_rows.push_back(
        {evaluator->evaluator_fingerprint(), index, output_count, next});
    output_count = next;
  }
  if (output_count != request.observations.size()) {
    result.message = "concatenated evaluator outputs do not match observations";
    return result;
  }
  result.parameter_contract_fingerprint = parameter_contract_fingerprint;

  const std::vector<double> residual_scales = expand_residual_scales(request);
  EvaluationState state(request, result.evaluator_rows, parameter_count);
  const auto initial = state.evaluate(request.initial_parameters, true,
                                      "initial_evaluation", false);
  if (std::holds_alternative<interfaces::EvaluationFailure>(initial)) {
    result.status = "initial_evaluation_failed";
    result.message = "initial owner evaluation failed";
    result.failures = state.failures();
    result.counters = state.counters();
    return result;
  }

  std::vector<double> optimizer_parameters(parameter_count);
  std::vector<double> optimizer_lower(parameter_count);
  std::vector<double> optimizer_upper(parameter_count);
  for (std::size_t index = 0; index < parameter_count; ++index) {
    optimizer_parameters[index] =
        (request.initial_parameters[index] - request.affine_origins[index]) /
        request.affine_scales[index];
    optimizer_lower[index] =
        (request.lower_bounds[index] - request.affine_origins[index]) /
        request.affine_scales[index];
    optimizer_upper[index] =
        (request.upper_bounds[index] - request.affine_origins[index]) /
        request.affine_scales[index];
  }
  if (!finite(optimizer_parameters) || !finite_or_infinite(optimizer_lower) ||
      !finite_or_infinite(optimizer_upper)) {
    result.message = "affine optimizer coordinates and bounds must not be NaN";
    result.failures = state.failures();
    result.counters = state.counters();
    return result;
  }
  for (std::size_t index = 0; index < parameter_count; ++index) {
    if (optimizer_parameters[index] < optimizer_lower[index] ||
        optimizer_parameters[index] > optimizer_upper[index]) {
      result.message = "initial optimizer coordinates are outside transformed bounds";
      result.failures = state.failures();
      result.counters = state.counters();
      return result;
    }
  }

  const auto &initial_evaluation = std::get<CombinedEvaluation>(initial);
  for (std::size_t row = 0; row < output_count; ++row) {
    const double residual =
        (initial_evaluation.predictions[row] - request.observations[row]) /
        residual_scales[row];
    if (!std::isfinite(residual)) {
      result.message = "initial scaled residual is not finite";
      result.failures = state.failures();
      result.counters = state.counters();
      return result;
    }
    for (std::size_t column = 0; column < parameter_count; ++column) {
      const double derivative =
          initial_evaluation.jacobian[row * parameter_count + column] *
          request.affine_scales[column] / residual_scales[row];
      if (!std::isfinite(derivative)) {
        result.message = "initial optimizer Jacobian is not finite";
        result.failures = state.failures();
        result.counters = state.counters();
        return result;
      }
    }
  }

  OwnerCostFunction cost_function(state, request, residual_scales);
  ceres::Problem::Options problem_options;
  problem_options.cost_function_ownership = ceres::DO_NOT_TAKE_OWNERSHIP;
  ceres::Problem problem(problem_options);
  problem.AddResidualBlock(&cost_function, nullptr, optimizer_parameters.data());
  for (std::size_t index = 0; index < parameter_count; ++index) {
    problem.SetParameterLowerBound(optimizer_parameters.data(),
                                   static_cast<int>(index),
                                   optimizer_lower[index]);
    problem.SetParameterUpperBound(optimizer_parameters.data(),
                                   static_cast<int>(index),
                                   optimizer_upper[index]);
  }

  ceres::Solver::Options options;
  apply_options(result.requested_options, options);
  ceres::Solver::Summary summary;
  ceres::Solve(options, &problem, &summary);

  result.termination = termination_name(summary.termination_type);
  result.resolved_options.minimizer_type =
      minimizer_name(summary.minimizer_type);
  result.resolved_options.trust_region_strategy_type =
      trust_region_name(summary.trust_region_strategy_type);
  result.resolved_options.linear_solver_type =
      linear_solver_name(summary.linear_solver_type_used);
  result.resolved_options.dense_linear_algebra_library_type =
      dense_library_name(summary.dense_linear_algebra_library_type);
  result.resolved_options.number_of_threads = summary.num_threads_used;
  result.resolved_options.use_inner_iterations = summary.inner_iterations_used;
  result.resolved_options.use_mixed_precision_solves =
      summary.mixed_precision_solves_used;
  auto &counters = state.counters();
  counters.ceres_residual_evaluations =
      std::max(0, summary.num_residual_evaluations);
  counters.ceres_jacobian_evaluations =
      std::max(0, summary.num_jacobian_evaluations);
  counters.accepted_steps = std::max(0, summary.num_successful_steps);
  counters.rejected_steps = std::max(0, summary.num_unsuccessful_steps);
  counters.iterations = static_cast<int>(summary.iterations.size());
  counters.iteration_limit_reached =
      summary.termination_type == ceres::NO_CONVERGENCE &&
      counters.iterations >= controls.maximum_iterations;
  counters.elapsed_time_limit_reached =
      summary.termination_type == ceres::NO_CONVERGENCE &&
      summary.total_time_in_seconds >= controls.maximum_elapsed_time_seconds;

  if (!same_fixed_options(summary)) {
    result.status = "failed";
    result.message = "Ceres did not retain the requested solver configuration";
  } else if (state.required_jacobian_failed()) {
    result.status = "failed";
    result.message = "a required owner Jacobian evaluation failed";
  } else if (summary.termination_type == ceres::CONVERGENCE) {
    result.status = "converged";
  } else if (summary.termination_type == ceres::NO_CONVERGENCE ||
             summary.termination_type == ceres::USER_SUCCESS) {
    result.status = "limit_reached";
  } else {
    result.status = "failed";
    result.message = "Ceres terminated without a usable point";
  }

  if (result.status == "converged" || result.status == "limit_reached") {
    std::vector<double> physical(parameter_count);
    for (std::size_t index = 0; index < parameter_count; ++index) {
      physical[index] = request.affine_origins[index] +
                        request.affine_scales[index] * optimizer_parameters[index];
    }
    const bool physical_valid =
        finite(physical) &&
        std::equal(physical.begin(), physical.end(), request.lower_bounds.begin(),
                   [](double value, double lower) { return value >= lower; }) &&
        std::equal(physical.begin(), physical.end(), request.upper_bounds.begin(),
                   [](double value, double upper) { return value <= upper; });
    if (!physical_valid) {
      result.status = "failed";
      result.message = "final physical parameters are nonfinite or outside bounds";
    }
    const auto final = physical_valid
                           ? state.evaluate(physical, true, "final_validation", true)
                           : std::variant<CombinedEvaluation,
                                          interfaces::EvaluationFailure>(
                                 failure("invalid_final_parameters",
                                         "final physical parameters are invalid"));
    if (const auto *evaluation = std::get_if<CombinedEvaluation>(&final)) {
      result.parameters = std::move(physical);
      result.predictions = evaluation->predictions;
      result.physical_jacobian = evaluation->jacobian;
      result.optimizer_jacobian.resize(output_count * parameter_count);
      result.residuals.resize(output_count);
      for (std::size_t row = 0; row < output_count; ++row) {
        result.residuals[row] =
            (result.predictions[row] - request.observations[row]) /
            residual_scales[row];
        for (std::size_t column = 0; column < parameter_count; ++column) {
          result.optimizer_jacobian[row * parameter_count + column] =
              evaluation->jacobian[row * parameter_count + column] *
              request.affine_scales[column] / residual_scales[row];
        }
      }
      if (finite(result.residuals) && finite(result.optimizer_jacobian)) {
        result.usable = true;
      } else {
        result.status = "failed";
        result.message = "final scaled residuals or optimizer Jacobian are not finite";
        result.parameters.clear();
        result.predictions.clear();
        result.residuals.clear();
        result.physical_jacobian.clear();
        result.optimizer_jacobian.clear();
      }
    } else {
      result.status = "failed";
      if (physical_valid) {
        result.message = "fresh final owner value/Jacobian evaluation failed";
      }
    }
  }

  if (result.usable) {
    cost_function.set_phase("covariance");
    ceres::Covariance::Options covariance_options;
    covariance_options.algorithm_type = ceres::DENSE_SVD;
    covariance_options.null_space_rank = -1;
    covariance_options.apply_loss_function = true;
    covariance_options.num_threads = 1;
    covariance_options.min_reciprocal_condition_number = 1.0e-14;
    ceres::Covariance covariance(covariance_options);
    for (std::size_t index = 0; index < parameter_count; ++index) {
      result.covariance.active_bound =
          result.covariance.active_bound ||
          result.parameters[index] == request.lower_bounds[index] ||
          result.parameters[index] == request.upper_bounds[index];
    }
    result.covariance.optimizer_coordinates.resize(parameter_count *
                                                    parameter_count);
    if (!covariance.Compute({optimizer_parameters.data()}, &problem) ||
        !covariance.GetCovarianceBlock(
            optimizer_parameters.data(), optimizer_parameters.data(),
            result.covariance.optimizer_coordinates.data()) ||
        !finite(result.covariance.optimizer_coordinates)) {
      result.covariance.status = "unavailable";
      result.covariance.optimizer_coordinates.clear();
    } else {
      result.covariance.physical_coordinates.resize(parameter_count *
                                                     parameter_count);
      for (std::size_t row = 0; row < parameter_count; ++row) {
        for (std::size_t column = 0; column < parameter_count; ++column) {
          result.covariance.physical_coordinates[row * parameter_count + column] =
              request.affine_scales[row] *
              result.covariance.optimizer_coordinates[row * parameter_count +
                                                       column] *
              request.affine_scales[column];
        }
      }
      if (!finite(result.covariance.physical_coordinates)) {
        result.covariance.status = "unavailable";
        result.covariance.optimizer_coordinates.clear();
        result.covariance.physical_coordinates.clear();
      } else {
        result.covariance.status =
            result.covariance.active_bound ? "active_bound" : "available";
      }
    }
  }

  result.failures = state.failures();
  result.counters = state.counters();
  if (result.message.empty()) {
    result.message = result.usable ? "usable fitted point" : "unusable fitted point";
  }
  return result;
}

} // namespace epcsaft_engine::regression
