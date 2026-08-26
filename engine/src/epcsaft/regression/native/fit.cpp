#include "fit.hpp"

#include <Eigen/Dense>
#include <Eigen/SVD>
#include <coin/IpIpoptApplication.hpp>
#include <coin/IpSolveStatistics.hpp>
#include <coin/IpTNLP.hpp>
#include <pybind11/stl.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace py = pybind11;

namespace epcsaft::regression {
namespace {

constexpr double kOptimizerBasinAffineDistanceMax = 1.0e-6;
constexpr double kOptimizerBasinCostDifferenceFactor = 1.0e-8;

enum class ResidualKind {
  LogRatio,
  ScaledDifference,
  LeftCensored,
  RightCensored,
};

struct Target final {
  double observed;
  ResidualKind residual;
  double scale;
  double multiplier;
};

struct Observation final {
  double multiplier;
  std::vector<Target> targets;
};

struct OwnerFailure final {
  std::string code;
  std::string reason;
  std::size_t count{0};
};

struct OwnerEvaluation final {
  std::vector<double> values;
  std::vector<double> derivatives;
  std::string receipt_json;
};

struct Coordinate final {
  double origin;
  double scale;
  double lower_affine;
  double upper_affine;
};

struct Start final {
  std::string identity;
  std::vector<double> affine;
};

struct Controls final {
  std::vector<Start> starts;
  double affine_bound_margin;
  int maximum_iterations;
  double optimizer_tolerance;
};

struct FullSpaceSetup final {
  std::vector<double> initial;
  std::vector<double> lower;
  std::vector<double> upper;
  std::vector<double> constraint_lower;
  std::vector<double> constraint_upper;
  std::vector<int> jacobian_rows;
  std::vector<int> jacobian_columns;
};

struct Diagnostics final {
  std::vector<double> singular_values;
  std::vector<double> affine_covariance;
  int rank{0};
  int free_subspace_rank{0};
  double condition_number{std::numeric_limits<double>::infinity()};
  std::size_t degrees_of_freedom{0};
  std::size_t information_row_count{0};
  std::vector<std::size_t> free_coordinate_indices;
  double residual_variance{std::numeric_limits<double>::quiet_NaN()};
  std::string covariance_status{"not_evaluated"};
};

struct BatchOutcome final {
  std::size_t index{0};
  std::size_t row_count{0};
  std::string status{"not_evaluated"};
  std::string failure_code;
};

struct StartResult final {
  std::size_t index{0};
  std::string start_identity;
  std::vector<double> initial_affine;
  std::vector<double> final_affine;
  std::vector<double> predictions;
  std::vector<double> residuals;
  std::vector<double> affine_jacobian;
  std::vector<OwnerFailure> owner_failures;
  std::vector<std::string> owner_receipts_json;
  std::vector<BatchOutcome> batch_outcomes;
  Diagnostics diagnostics;
  std::vector<std::size_t> active_bounds;
  std::string solver_status;
  std::string failure_reason;
  std::string trial_failure_code;
  int iterations{0};
  std::size_t non_evaluable_count{0};
  double cost{std::numeric_limits<double>::infinity()};
  bool usable{false};
  std::size_t optimizer_basin{std::numeric_limits<std::size_t>::max()};
};

ResidualKind residual_kind(const std::string &value) {
  if (value == "log_ratio") {
    return ResidualKind::LogRatio;
  }
  if (value == "scaled_difference") {
    return ResidualKind::ScaledDifference;
  }
  if (value == "left_censored") {
    return ResidualKind::LeftCensored;
  }
  if (value == "right_censored") {
    return ResidualKind::RightCensored;
  }
  throw std::invalid_argument("unknown residual kind: " + value);
}

template <typename Value>
Value required(const py::dict &mapping, const char *name) {
  if (!mapping.contains(name)) {
    throw std::invalid_argument(std::string("missing field: ") + name);
  }
  return py::cast<Value>(mapping[name]);
}

std::vector<Observation> parse_observations(const py::iterable &values) {
  std::vector<Observation> observations;
  for (const py::handle item : values) {
    const auto mapping = py::cast<py::dict>(item);
    Observation observation{
        required<double>(mapping, "multiplier"),
        {},
    };
    if (!(std::isfinite(observation.multiplier) &&
          observation.multiplier > 0.0)) {
      throw std::invalid_argument(
          "observation multiplier must be positive and finite");
    }
    for (const py::handle target_item :
         required<py::iterable>(mapping, "targets")) {
      const auto target_mapping = py::cast<py::dict>(target_item);
      Target target{
          required<double>(target_mapping, "observed"),
          residual_kind(required<std::string>(target_mapping, "residual")),
          required<double>(target_mapping, "scale"),
          target_mapping.contains("multiplier")
              ? py::cast<double>(target_mapping["multiplier"])
              : 1.0,
      };
      if (!(std::isfinite(target.observed) && std::isfinite(target.scale) &&
            target.scale > 0.0 && std::isfinite(target.multiplier) &&
            target.multiplier > 0.0)) {
        throw std::invalid_argument(
            "target values, scales, and multipliers must be finite");
      }
      if (target.residual == ResidualKind::LogRatio && target.observed <= 0.0) {
        throw std::invalid_argument(
            "log-ratio observations must be strictly positive");
      }
      observation.targets.push_back(target);
    }
    if (observation.targets.empty()) {
      throw std::invalid_argument("observation targets must be nonempty");
    }
    observations.push_back(std::move(observation));
  }
  if (observations.empty()) {
    throw std::invalid_argument("observations must be nonempty");
  }
  return observations;
}

std::vector<Coordinate> parse_coordinates(const py::iterable &values) {
  std::vector<Coordinate> coordinates;
  for (const py::handle item : values) {
    const auto mapping = py::cast<py::dict>(item);
    Coordinate coordinate{
        required<double>(mapping, "origin"),
        required<double>(mapping, "scale"),
        required<double>(mapping, "lower_affine"),
        required<double>(mapping, "upper_affine"),
    };
    if (!(std::isfinite(coordinate.origin) && std::isfinite(coordinate.scale) &&
          coordinate.scale > 0.0 && std::isfinite(coordinate.lower_affine) &&
          std::isfinite(coordinate.upper_affine) &&
          coordinate.lower_affine < coordinate.upper_affine)) {
      throw std::invalid_argument("coordinate origin, scale, and affine bounds "
                                  "must be finite and consistent");
    }
    coordinates.push_back(coordinate);
  }
  if (coordinates.empty()) {
    throw std::invalid_argument("coordinates must be nonempty");
  }
  return coordinates;
}

Controls parse_controls(const py::dict &mapping) {
  const auto raw_starts = required<py::iterable>(mapping, "starts");
  std::vector<Start> starts;
  for (const py::handle item : raw_starts) {
    const auto start = py::cast<py::dict>(item);
    starts.push_back({required<std::string>(start, "identity"),
                      required<std::vector<double>>(start, "affine")});
  }
  Controls controls{
      std::move(starts),
      required<double>(mapping, "affine_bound_margin"),
      required<int>(mapping, "maximum_iterations"),
      required<double>(mapping, "optimizer_tolerance"),
  };
  if (controls.starts.empty() || controls.maximum_iterations <= 0 ||
      !(std::isfinite(controls.affine_bound_margin) &&
        controls.affine_bound_margin >= 0.0 &&
        std::isfinite(controls.optimizer_tolerance) &&
        controls.optimizer_tolerance > 0.0)) {
    throw std::invalid_argument("solver controls must be positive and finite");
  }
  return controls;
}

std::vector<FullSpaceSetup> parse_full_space_setups(
    const py::iterable &values) {
  std::vector<FullSpaceSetup> setups;
  for (const py::handle item : values) {
    const auto mapping = py::cast<py::dict>(item);
    FullSpaceSetup setup{
        required<std::vector<double>>(mapping, "initial"),
        required<std::vector<double>>(mapping, "lower"),
        required<std::vector<double>>(mapping, "upper"),
        required<std::vector<double>>(mapping, "constraint_lower"),
        required<std::vector<double>>(mapping, "constraint_upper"),
        required<std::vector<int>>(mapping, "jacobian_rows"),
        required<std::vector<int>>(mapping, "jacobian_columns"),
    };
    if (setup.initial.empty() || setup.lower.size() != setup.initial.size() ||
        setup.upper.size() != setup.initial.size() ||
        setup.constraint_lower.empty() ||
        setup.constraint_upper.size() != setup.constraint_lower.size() ||
        setup.jacobian_rows.size() != setup.jacobian_columns.size()) {
      throw std::invalid_argument("full-space setup dimensions are invalid");
    }
    for (std::size_t index = 0; index < setup.initial.size(); ++index) {
      if (!std::isfinite(setup.initial[index]) ||
          !std::isfinite(setup.lower[index]) ||
          !std::isfinite(setup.upper[index]) ||
          setup.lower[index] >= setup.upper[index] ||
          setup.initial[index] < setup.lower[index] ||
          setup.initial[index] > setup.upper[index]) {
        throw std::invalid_argument("full-space local bounds are invalid");
      }
    }
    for (std::size_t index = 0; index < setup.constraint_lower.size(); ++index) {
      if (!std::isfinite(setup.constraint_lower[index]) ||
          !std::isfinite(setup.constraint_upper[index]) ||
          setup.constraint_lower[index] > setup.constraint_upper[index]) {
        throw std::invalid_argument("full-space constraint bounds are invalid");
      }
    }
    setups.push_back(std::move(setup));
  }
  return setups;
}

class OwnerEvaluationError final : public std::runtime_error {
public:
  OwnerEvaluationError(std::string code, std::string reason)
      : std::runtime_error(std::move(reason)), code_(std::move(code)) {}

  const std::string &code() const noexcept { return code_; }

private:
  std::string code_;
};

class OwnerCost final {
public:
  OwnerCost(PyObject *owner_evaluator, std::size_t start_index,
            std::size_t observation_index, std::size_t active_parameter_count,
            std::vector<std::size_t> identity_coordinates,
            std::vector<double> coordinate_origins,
            std::vector<double> coordinate_scales, Observation observation,
            std::vector<std::string> failure_codes)
      : owner_evaluator_(owner_evaluator), start_index_(start_index),
        observation_index_(observation_index),
        active_parameter_count_(active_parameter_count),
        identity_coordinates_(std::move(identity_coordinates)),
        coordinate_origins_(std::move(coordinate_origins)),
        coordinate_scales_(std::move(coordinate_scales)),
        observation_(std::move(observation)),
        failure_codes_(std::move(failure_codes)) {}

  std::size_t row_count() const noexcept { return observation_.targets.size(); }

  bool evaluate(double const *const *parameters, double *residuals,
                double **jacobians, double *predictions, bool receipt_requested,
                std::string *receipt_json) const {
    try {
      std::vector<double> active(active_parameter_count_);
      for (std::size_t index = 0; index < active_parameter_count_; ++index) {
        const std::size_t coordinate = identity_coordinates_[index];
        active[index] =
            coordinate_origins_[coordinate] +
            parameters[coordinate][0] * coordinate_scales_[coordinate];
      }
      const OwnerEvaluation evaluation =
          evaluate_owner(active, receipt_requested);
      if (receipt_json != nullptr) {
        *receipt_json = evaluation.receipt_json;
      }
      for (std::size_t row = 0; row < observation_.targets.size(); ++row) {
        const auto &target = observation_.targets[row];
        const double predicted = evaluation.values[row];
        if (!std::isfinite(predicted) ||
            (target.residual == ResidualKind::LogRatio && predicted <= 0.0)) {
          return reject("observation_out_of_support",
                        "owner returned an inadmissible target value");
        }
        if (predictions != nullptr) {
          predictions[row] = predicted;
        }
        const double multiplier = observation_.multiplier * target.multiplier;
        double factor = multiplier / target.scale;
        switch (target.residual) {
        case ResidualKind::LogRatio:
          factor = multiplier / predicted;
          residuals[row] = multiplier * std::log(predicted / target.observed);
          break;
        case ResidualKind::ScaledDifference:
          residuals[row] =
              multiplier * (predicted - target.observed) / target.scale;
          break;
        case ResidualKind::LeftCensored:
          if (predicted <= target.observed) {
            factor = 0.0;
            residuals[row] = 0.0;
          } else {
            residuals[row] =
                multiplier * (predicted - target.observed) / target.scale;
          }
          break;
        case ResidualKind::RightCensored:
          if (predicted >= target.observed) {
            factor = 0.0;
            residuals[row] = 0.0;
          } else {
            residuals[row] =
                multiplier * (predicted - target.observed) / target.scale;
          }
          break;
        }
        if (jacobians == nullptr) {
          continue;
        }
        for (std::size_t coordinate = 0; coordinate < coordinate_scales_.size();
             ++coordinate) {
          if (jacobians[coordinate] == nullptr) {
            continue;
          }
          double derivative = 0.0;
          for (std::size_t identity = 0; identity < active_parameter_count_;
               ++identity) {
            if (identity_coordinates_[identity] == coordinate) {
              derivative +=
                  evaluation
                      .derivatives[row * active_parameter_count_ + identity];
            }
          }
          jacobians[coordinate][row] =
              factor * derivative * coordinate_scales_[coordinate];
        }
      }
      if (!std::all_of(residuals, residuals + observation_.targets.size(),
                       [](double value) { return std::isfinite(value); })) {
        return reject("native_residual_failure",
                      "residual evaluation was nonfinite");
      }
      return true;
    } catch (const OwnerEvaluationError &error) {
      return reject(error.code(), error.what());
    } catch (const std::exception &error) {
      return reject("owner_contract_failure", error.what());
    }
  }

  std::size_t rejection_count() const noexcept {
    return std::accumulate(owner_failures_.begin(), owner_failures_.end(),
                           std::size_t{0},
                           [](std::size_t count, const OwnerFailure &failure) {
                             return count + failure.count;
                           });
  }
  const std::vector<OwnerFailure> &owner_failures() const noexcept {
    return owner_failures_;
  }

  bool information_bearing(std::size_t row, double predicted) const noexcept {
    const auto &target = observation_.targets[row];
    if (target.residual == ResidualKind::LeftCensored) {
      return predicted > target.observed;
    }
    if (target.residual == ResidualKind::RightCensored) {
      return predicted < target.observed;
    }
    return true;
  }

  bool censor_kink(std::size_t row, double predicted) const noexcept {
    const auto &target = observation_.targets[row];
    if (target.residual != ResidualKind::LeftCensored &&
        target.residual != ResidualKind::RightCensored) {
      return false;
    }
    const double tolerance = 1.0e-12 * std::max({1.0, std::abs(predicted),
                                                 std::abs(target.observed)});
    return std::abs(predicted - target.observed) <= tolerance;
  }

private:
  OwnerEvaluation evaluate_owner(const std::vector<double> &active,
                                 bool receipt_requested) const {
    py::gil_scoped_acquire acquire;
    const auto evaluator =
        py::reinterpret_borrow<py::function>(py::handle(owner_evaluator_));
    const py::dict result = py::cast<py::dict>(
        evaluator(start_index_, observation_index_, active, receipt_requested));
    const std::string status = required<std::string>(result, "status");
    if (status == "rejected") {
      const auto failure = required<py::dict>(result, "failure");
      const std::string code = required<std::string>(failure, "code");
      const std::string reason = required<std::string>(failure, "reason");
      if (std::find(failure_codes_.begin(), failure_codes_.end(), code) ==
          failure_codes_.end()) {
        throw OwnerEvaluationError("owner_contract_failure",
                                   "owner returned unknown failure code '" +
                                       code + "': " + reason);
      }
      throw OwnerEvaluationError(code, reason);
    }
    if (status != "accepted") {
      throw OwnerEvaluationError("owner_contract_failure",
                                 "unknown owner evaluation status");
    }
    OwnerEvaluation evaluation{
        required<std::vector<double>>(result, "predictions"),
        required<std::vector<double>>(result, "jacobian"),
        result.contains("receipt_json")
            ? py::cast<std::string>(result["receipt_json"])
            : std::string{},
    };
    if (evaluation.values.size() != observation_.targets.size() ||
        evaluation.derivatives.size() !=
            observation_.targets.size() * active_parameter_count_) {
      throw OwnerEvaluationError(
          "owner_contract_failure",
          "owner block has the wrong value or Jacobian shape");
    }
    if (!std::all_of(evaluation.derivatives.begin(),
                     evaluation.derivatives.end(),
                     [](double value) { return std::isfinite(value); })) {
      throw OwnerEvaluationError(
          "owner_contract_failure",
          "owner block contains a nonfinite Jacobian derivative");
    }
    if (!receipt_requested && !evaluation.receipt_json.empty()) {
      throw OwnerEvaluationError(
          "owner_contract_failure",
          "owner returned a receipt for an ordinary trial evaluation");
    }
    return evaluation;
  }

  bool reject(const std::string &code, const std::string &reason) const {
    const auto match =
        std::find_if(owner_failures_.begin(), owner_failures_.end(),
                     [&](const OwnerFailure &failure) {
                       return failure.code == code && failure.reason == reason;
                     });
    if (match == owner_failures_.end()) {
      owner_failures_.push_back({code, reason, 1});
    } else {
      ++match->count;
    }
    return false;
  }

  PyObject *owner_evaluator_;
  std::size_t start_index_;
  std::size_t observation_index_;
  std::size_t active_parameter_count_;
  std::vector<std::size_t> identity_coordinates_;
  std::vector<double> coordinate_origins_;
  std::vector<double> coordinate_scales_;
  Observation observation_;
  std::vector<std::string> failure_codes_;
  mutable std::vector<OwnerFailure> owner_failures_;
};

Diagnostics diagnose(const std::vector<double> &jacobian,
                     const std::vector<double> &residuals, std::size_t rows,
                     std::size_t columns,
                     const std::vector<std::size_t> &active_bounds,
                     const std::vector<bool> &information_rows,
                     bool censor_kink) {
  Diagnostics result;
  for (std::size_t coordinate = 0; coordinate < columns; ++coordinate) {
    if (std::find(active_bounds.begin(), active_bounds.end(), coordinate) ==
        active_bounds.end()) {
      result.free_coordinate_indices.push_back(coordinate);
    }
  }
  result.information_row_count = static_cast<std::size_t>(
      std::count(information_rows.begin(), information_rows.end(), true));
  if (result.information_row_count == 0) {
    result.covariance_status = "unavailable_no_information_rows";
    return result;
  }
  Eigen::MatrixXd information_matrix(
      static_cast<Eigen::Index>(result.information_row_count),
      static_cast<Eigen::Index>(columns));
  Eigen::Index information_row = 0;
  for (std::size_t row = 0; row < rows; ++row) {
    if (!information_rows[row]) {
      continue;
    }
    for (std::size_t column = 0; column < columns; ++column) {
      information_matrix(information_row, static_cast<Eigen::Index>(column)) =
          jacobian[row * columns + column];
    }
    ++information_row;
  }
  const Eigen::JacobiSVD<Eigen::MatrixXd> svd(
      information_matrix, Eigen::ComputeThinU | Eigen::ComputeThinV);
  result.singular_values.assign(svd.singularValues().data(),
                                svd.singularValues().data() +
                                    svd.singularValues().size());
  if (result.singular_values.empty()) {
    result.condition_number = 1.0;
    result.covariance_status = "unavailable_rank_deficient";
    return result;
  }
  const double threshold =
      100.0 * std::numeric_limits<double>::epsilon() *
      static_cast<double>(std::max(result.information_row_count, columns)) *
      result.singular_values.front();
  result.rank = static_cast<int>(std::count_if(
      result.singular_values.begin(), result.singular_values.end(),
      [threshold](double value) { return value > threshold; }));
  if (result.rank == static_cast<int>(columns)) {
    result.condition_number =
        result.singular_values.front() /
        result.singular_values[static_cast<std::size_t>(result.rank - 1)];
  }
  Eigen::MatrixXd free_matrix(
      static_cast<Eigen::Index>(result.information_row_count),
      static_cast<Eigen::Index>(result.free_coordinate_indices.size()));
  for (std::size_t column = 0; column < result.free_coordinate_indices.size();
       ++column) {
    free_matrix.col(static_cast<Eigen::Index>(column)) = information_matrix.col(
        static_cast<Eigen::Index>(result.free_coordinate_indices[column]));
  }
  if (result.free_coordinate_indices.empty()) {
    result.degrees_of_freedom = result.information_row_count;
    result.covariance_status = "unavailable_active_bounds";
    return result;
  }
  Eigen::JacobiSVD<Eigen::MatrixXd> free_svd(
      free_matrix, Eigen::ComputeThinU | Eigen::ComputeThinV);
  if (!result.free_coordinate_indices.empty()) {
    const double free_threshold =
        100.0 * std::numeric_limits<double>::epsilon() *
        static_cast<double>(std::max(result.information_row_count,
                                     result.free_coordinate_indices.size())) *
        free_svd.singularValues()[0];
    result.free_subspace_rank = static_cast<int>(std::count_if(
        free_svd.singularValues().data(),
        free_svd.singularValues().data() + free_svd.singularValues().size(),
        [free_threshold](double value) { return value > free_threshold; }));
  }
  result.degrees_of_freedom =
      result.information_row_count >
              static_cast<std::size_t>(result.free_subspace_rank)
          ? result.information_row_count -
                static_cast<std::size_t>(result.free_subspace_rank)
          : 0;
  if (!active_bounds.empty()) {
    result.covariance_status = "unavailable_active_bounds";
    return result;
  }
  if (censor_kink) {
    result.covariance_status = "unavailable_censor_kink";
    return result;
  }
  if (result.free_subspace_rank !=
      static_cast<int>(result.free_coordinate_indices.size())) {
    result.covariance_status = "unavailable_rank_deficient";
    return result;
  }
  if (result.degrees_of_freedom == 0) {
    result.covariance_status = "unavailable_insufficient_degrees_of_freedom";
    return result;
  }
  double information_sum_squares = 0.0;
  for (std::size_t row = 0; row < rows; ++row) {
    if (information_rows[row]) {
      information_sum_squares += residuals[row] * residuals[row];
    }
  }
  result.residual_variance =
      information_sum_squares / static_cast<double>(result.degrees_of_freedom);
  const Eigen::VectorXd inverse_variances =
      free_svd.singularValues().array().square().inverse() *
      result.residual_variance;
  const Eigen::MatrixXd covariance = free_svd.matrixV() *
                                     inverse_variances.asDiagonal() *
                                     free_svd.matrixV().transpose();
  result.affine_covariance.reserve(result.free_coordinate_indices.size() *
                                   result.free_coordinate_indices.size());
  for (std::size_t row = 0; row < result.free_coordinate_indices.size();
       ++row) {
    for (std::size_t column = 0; column < result.free_coordinate_indices.size();
         ++column) {
      result.affine_covariance.push_back(covariance(
          static_cast<Eigen::Index>(row), static_cast<Eigen::Index>(column)));
    }
  }
  result.covariance_status = "available";
  return result;
}

std::string ipopt_status_name(Ipopt::ApplicationReturnStatus status) {
  switch (status) {
  case Ipopt::Solve_Succeeded:
    return "solve_succeeded";
  case Ipopt::Solved_To_Acceptable_Level:
    return "solved_to_acceptable_level";
  case Ipopt::Infeasible_Problem_Detected:
    return "infeasible_problem_detected";
  case Ipopt::Search_Direction_Becomes_Too_Small:
    return "search_direction_too_small";
  case Ipopt::Diverging_Iterates:
    return "diverging_iterates";
  case Ipopt::Feasible_Point_Found:
    return "feasible_point_found";
  case Ipopt::Maximum_Iterations_Exceeded:
    return "maximum_iterations_exceeded";
  case Ipopt::Restoration_Failed:
    return "restoration_failed";
  case Ipopt::Error_In_Step_Computation:
    return "error_in_step_computation";
  case Ipopt::Not_Enough_Degrees_Of_Freedom:
    return "not_enough_degrees_of_freedom";
  case Ipopt::Invalid_Problem_Definition:
    return "invalid_problem_definition";
  case Ipopt::Invalid_Option:
    return "invalid_option";
  case Ipopt::Invalid_Number_Detected:
    return "invalid_number_detected";
  case Ipopt::Unrecoverable_Exception:
    return "unrecoverable_exception";
  case Ipopt::NonIpopt_Exception_Thrown:
    return "non_ipopt_exception_thrown";
  case Ipopt::Insufficient_Memory:
    return "insufficient_memory";
  case Ipopt::Internal_Error:
    return "internal_error";
  default:
    return "ipopt_status_" + std::to_string(static_cast<int>(status));
  }
}

bool ipopt_solution_usable(Ipopt::ApplicationReturnStatus status) {
  return status == Ipopt::Solve_Succeeded ||
         status == Ipopt::Solved_To_Acceptable_Level ||
         status == Ipopt::Search_Direction_Becomes_Too_Small ||
         status == Ipopt::Feasible_Point_Found ||
         status == Ipopt::Maximum_Iterations_Exceeded;
}

class ReducedSpaceNlp final : public Ipopt::TNLP {
public:
  ReducedSpaceNlp(const std::vector<std::unique_ptr<OwnerCost>> &costs,
                  const std::vector<Coordinate> &coordinates,
                  std::vector<double> initial)
      : costs_(costs), coordinates_(coordinates), initial_(std::move(initial)),
        solution_(initial_) {}

  bool get_nlp_info(Ipopt::Index &n, Ipopt::Index &m, Ipopt::Index &nnz_jac_g,
                    Ipopt::Index &nnz_h_lag,
                    IndexStyleEnum &index_style) override {
    n = static_cast<Ipopt::Index>(coordinates_.size());
    m = 0;
    nnz_jac_g = 0;
    nnz_h_lag = 0;
    index_style = TNLP::C_STYLE;
    return true;
  }

  bool get_bounds_info(Ipopt::Index n, Ipopt::Number *x_l, Ipopt::Number *x_u,
                       Ipopt::Index m, Ipopt::Number *,
                       Ipopt::Number *) override {
    if (n != static_cast<Ipopt::Index>(coordinates_.size()) || m != 0) {
      return false;
    }
    for (std::size_t index = 0; index < coordinates_.size(); ++index) {
      x_l[index] = coordinates_[index].lower_affine;
      x_u[index] = coordinates_[index].upper_affine;
    }
    return true;
  }

  bool get_starting_point(Ipopt::Index n, bool init_x, Ipopt::Number *x,
                          bool init_z, Ipopt::Number *, Ipopt::Number *,
                          Ipopt::Index m, bool init_lambda,
                          Ipopt::Number *) override {
    if (n != static_cast<Ipopt::Index>(initial_.size()) || m != 0 || !init_x ||
        init_z || init_lambda) {
      return false;
    }
    std::copy(initial_.begin(), initial_.end(), x);
    return true;
  }

  bool eval_f(Ipopt::Index n, const Ipopt::Number *x, bool,
              Ipopt::Number &objective) override {
    if (!evaluate(n, x)) {
      return false;
    }
    objective = cached_objective_;
    return true;
  }

  bool eval_grad_f(Ipopt::Index n, const Ipopt::Number *x, bool,
                   Ipopt::Number *gradient) override {
    if (!evaluate(n, x)) {
      return false;
    }
    std::copy(cached_gradient_.begin(), cached_gradient_.end(), gradient);
    return true;
  }

  bool eval_g(Ipopt::Index, const Ipopt::Number *, bool, Ipopt::Index m,
              Ipopt::Number *) override {
    return m == 0;
  }

  bool eval_jac_g(Ipopt::Index, const Ipopt::Number *, bool, Ipopt::Index m,
                  Ipopt::Index nonzero_count, Ipopt::Index *, Ipopt::Index *,
                  Ipopt::Number *) override {
    return m == 0 && nonzero_count == 0;
  }

  bool eval_h(Ipopt::Index, const Ipopt::Number *, bool, Ipopt::Number,
              Ipopt::Index m, const Ipopt::Number *, bool,
              Ipopt::Index nonzero_count, Ipopt::Index *, Ipopt::Index *,
              Ipopt::Number *) override {
    return m == 0 && nonzero_count == 0;
  }

  void finalize_solution(Ipopt::SolverReturn, Ipopt::Index n,
                         const Ipopt::Number *x, const Ipopt::Number *,
                         const Ipopt::Number *, Ipopt::Index m,
                         const Ipopt::Number *, const Ipopt::Number *,
                         Ipopt::Number, const Ipopt::IpoptData *,
                         Ipopt::IpoptCalculatedQuantities *) override {
    if (n == static_cast<Ipopt::Index>(solution_.size()) && m == 0 &&
        x != nullptr) {
      std::copy(x, x + n, solution_.begin());
    }
  }

  const std::vector<double> &solution() const noexcept { return solution_; }

private:
  bool evaluate(Ipopt::Index n, const Ipopt::Number *x) const {
    if (n != static_cast<Ipopt::Index>(coordinates_.size())) {
      return false;
    }
    if (cached_x_.size() == coordinates_.size() &&
        std::equal(x, x + n, cached_x_.begin())) {
      return cached_evaluable_;
    }
    cached_x_.assign(x, x + n);
    cached_objective_ = 0.0;
    cached_gradient_.assign(coordinates_.size(), 0.0);
    std::vector<const double *> blocks;
    blocks.reserve(coordinates_.size());
    for (Ipopt::Index index = 0; index < n; ++index) {
      blocks.push_back(x + index);
    }
    for (const auto &cost : costs_) {
      const std::size_t rows = cost->row_count();
      std::vector<double> residuals(rows);
      std::vector<std::vector<double>> jacobian_blocks(
          coordinates_.size(), std::vector<double>(rows));
      std::vector<double *> jacobian_pointers;
      jacobian_pointers.reserve(coordinates_.size());
      for (auto &block : jacobian_blocks) {
        jacobian_pointers.push_back(block.data());
      }
      if (!cost->evaluate(blocks.data(), residuals.data(),
                          jacobian_pointers.data(), nullptr, false, nullptr)) {
        cached_evaluable_ = false;
        return false;
      }
      for (std::size_t row = 0; row < rows; ++row) {
        cached_objective_ += 0.5 * residuals[row] * residuals[row];
        for (std::size_t coordinate = 0; coordinate < coordinates_.size();
             ++coordinate) {
          cached_gradient_[coordinate] +=
              residuals[row] * jacobian_blocks[coordinate][row];
        }
      }
    }
    cached_evaluable_ =
        std::isfinite(cached_objective_) &&
        std::all_of(cached_gradient_.begin(), cached_gradient_.end(),
                    [](double value) { return std::isfinite(value); });
    return cached_evaluable_;
  }

  const std::vector<std::unique_ptr<OwnerCost>> &costs_;
  const std::vector<Coordinate> &coordinates_;
  std::vector<double> initial_;
  std::vector<double> solution_;
  mutable std::vector<double> cached_x_;
  mutable std::vector<double> cached_gradient_;
  mutable double cached_objective_{0.0};
  mutable bool cached_evaluable_{false};
};

class SparseFullSpaceNlp final : public Ipopt::TNLP {
public:
  SparseFullSpaceNlp(PyObject *evaluator, std::size_t start_index,
                     std::vector<double> initial, std::vector<double> lower,
                     std::vector<double> upper, std::size_t residual_count,
                     std::vector<double> constraint_lower,
                     std::vector<double> constraint_upper,
                     std::vector<int> jacobian_rows,
                     std::vector<int> jacobian_columns)
      : evaluator_(evaluator), start_index_(start_index),
        initial_(std::move(initial)),
        lower_(std::move(lower)), upper_(std::move(upper)),
        residual_count_(residual_count),
        constraint_lower_(std::move(constraint_lower)),
        constraint_upper_(std::move(constraint_upper)),
        jacobian_rows_(std::move(jacobian_rows)),
        jacobian_columns_(std::move(jacobian_columns)), solution_(initial_) {
    Py_INCREF(evaluator_);
  }

  ~SparseFullSpaceNlp() override {
    py::gil_scoped_acquire acquire;
    Py_DECREF(evaluator_);
  }

  bool get_nlp_info(Ipopt::Index &n, Ipopt::Index &m,
                    Ipopt::Index &nnz_jac_g, Ipopt::Index &nnz_h_lag,
                    IndexStyleEnum &index_style) override {
    n = static_cast<Ipopt::Index>(initial_.size());
    m = static_cast<Ipopt::Index>(constraint_lower_.size());
    nnz_jac_g = static_cast<Ipopt::Index>(jacobian_rows_.size());
    nnz_h_lag = 0;
    index_style = TNLP::C_STYLE;
    return true;
  }

  bool get_bounds_info(Ipopt::Index n, Ipopt::Number *x_l,
                       Ipopt::Number *x_u, Ipopt::Index m,
                       Ipopt::Number *g_l, Ipopt::Number *g_u) override {
    if (!dimensions(n, m)) {
      return false;
    }
    std::copy(lower_.begin(), lower_.end(), x_l);
    std::copy(upper_.begin(), upper_.end(), x_u);
    std::copy(constraint_lower_.begin(), constraint_lower_.end(), g_l);
    std::copy(constraint_upper_.begin(), constraint_upper_.end(), g_u);
    return true;
  }

  bool get_starting_point(Ipopt::Index n, bool init_x, Ipopt::Number *x,
                          bool init_z, Ipopt::Number *, Ipopt::Number *,
                          Ipopt::Index m, bool init_lambda,
                          Ipopt::Number *) override {
    if (!dimensions(n, m) || !init_x || init_z || init_lambda) {
      return false;
    }
    std::copy(initial_.begin(), initial_.end(), x);
    return true;
  }

  bool eval_f(Ipopt::Index n, const Ipopt::Number *x, bool,
              Ipopt::Number &objective) override {
    if (!evaluate(n, x)) {
      return false;
    }
    objective = objective_;
    return true;
  }

  bool eval_grad_f(Ipopt::Index n, const Ipopt::Number *x, bool,
                   Ipopt::Number *gradient) override {
    if (!evaluate(n, x)) {
      return false;
    }
    std::copy(gradient_.begin(), gradient_.end(), gradient);
    return true;
  }

  bool eval_g(Ipopt::Index n, const Ipopt::Number *x, bool, Ipopt::Index m,
              Ipopt::Number *constraints) override {
    if (!dimensions(n, m) || !evaluate(n, x)) {
      return false;
    }
    std::copy(constraints_.begin(), constraints_.end(), constraints);
    return true;
  }

  bool eval_jac_g(Ipopt::Index n, const Ipopt::Number *x, bool,
                  Ipopt::Index m, Ipopt::Index nonzero_count,
                  Ipopt::Index *rows, Ipopt::Index *columns,
                  Ipopt::Number *values) override {
    if (!dimensions(n, m) ||
        nonzero_count != static_cast<Ipopt::Index>(jacobian_rows_.size())) {
      return false;
    }
    if (values == nullptr) {
      for (std::size_t index = 0; index < jacobian_rows_.size(); ++index) {
        rows[index] = jacobian_rows_[index];
        columns[index] = jacobian_columns_[index];
      }
      return true;
    }
    if (!evaluate(n, x)) {
      return false;
    }
    std::copy(constraint_jacobian_.begin(), constraint_jacobian_.end(), values);
    return true;
  }

  bool eval_h(Ipopt::Index, const Ipopt::Number *, bool, Ipopt::Number,
              Ipopt::Index, const Ipopt::Number *, bool,
              Ipopt::Index nonzero_count, Ipopt::Index *, Ipopt::Index *,
              Ipopt::Number *) override {
    return nonzero_count == 0;
  }

  void finalize_solution(Ipopt::SolverReturn, Ipopt::Index n,
                         const Ipopt::Number *x, const Ipopt::Number *,
                         const Ipopt::Number *, Ipopt::Index m,
                         const Ipopt::Number *, const Ipopt::Number *,
                         Ipopt::Number, const Ipopt::IpoptData *,
                         Ipopt::IpoptCalculatedQuantities *) override {
    if (dimensions(n, m) && x != nullptr && evaluate(n, x)) {
      std::copy(x, x + n, solution_.begin());
    }
  }

  const std::vector<double> &solution() const noexcept { return solution_; }
private:
  bool dimensions(Ipopt::Index n, Ipopt::Index m) const noexcept {
    return n == static_cast<Ipopt::Index>(initial_.size()) &&
           m == static_cast<Ipopt::Index>(constraint_lower_.size());
  }

  bool evaluate(Ipopt::Index n, const Ipopt::Number *x) const {
    if (n != static_cast<Ipopt::Index>(initial_.size()) || x == nullptr) {
      return false;
    }
    if (cached_x_.size() == initial_.size() &&
        std::equal(x, x + n, cached_x_.begin())) {
      return !callback_failed_;
    }
    cached_x_.assign(x, x + n);
    try {
      py::gil_scoped_acquire acquire;
      const py::dict result = py::cast<py::dict>(
          py::reinterpret_borrow<py::function>(py::handle(evaluator_))(
              start_index_, cached_x_));
      const auto residuals = required<std::vector<double>>(result, "residuals");
      const auto residual_jacobian =
          required<std::vector<double>>(result, "residual_jacobian");
      constraints_ = required<std::vector<double>>(result, "constraints");
      constraint_jacobian_ =
          required<std::vector<double>>(result, "constraint_jacobian");
      if (residuals.size() != residual_count_ ||
          residual_jacobian.size() != residual_count_ * initial_.size() ||
          constraints_.size() != constraint_lower_.size() ||
          constraint_jacobian_.size() != jacobian_rows_.size()) {
        throw std::invalid_argument("full-space callback shape changed");
      }
      objective_ = 0.0;
      gradient_.assign(initial_.size(), 0.0);
      for (std::size_t row = 0; row < residual_count_; ++row) {
        objective_ += 0.5 * residuals[row] * residuals[row];
        for (std::size_t column = 0; column < initial_.size(); ++column) {
          gradient_[column] += residuals[row] *
                               residual_jacobian[row * initial_.size() + column];
        }
      }
      const auto finite = [](const std::vector<double> &values) {
        return std::all_of(values.begin(), values.end(),
                           [](double value) { return std::isfinite(value); });
      };
      if (!std::isfinite(objective_) || !finite(gradient_) ||
          !finite(constraints_) || !finite(constraint_jacobian_)) {
        throw std::domain_error("full-space callback returned nonfinite values");
      }
      callback_failed_ = false;
      return true;
    } catch (const std::exception &) {
      callback_failed_ = true;
      return false;
    }
  }

  PyObject *evaluator_;
  std::size_t start_index_;
  std::vector<double> initial_;
  std::vector<double> lower_;
  std::vector<double> upper_;
  std::size_t residual_count_;
  std::vector<double> constraint_lower_;
  std::vector<double> constraint_upper_;
  std::vector<int> jacobian_rows_;
  std::vector<int> jacobian_columns_;
  std::vector<double> solution_;
  mutable std::vector<double> cached_x_;
  mutable std::vector<double> gradient_;
  mutable std::vector<double> constraints_;
  mutable std::vector<double> constraint_jacobian_;
  mutable double objective_{0.0};
  mutable bool callback_failed_{false};
};

template <typename Problem>
Ipopt::ApplicationReturnStatus solve_ipopt(Problem *problem,
                                           const Controls &controls,
                                           int &iterations,
                                           std::vector<double> &solution,
                                           bool full_space = false) {
  Ipopt::SmartPtr<Ipopt::TNLP> nlp = problem;
  Ipopt::SmartPtr<Ipopt::IpoptApplication> application =
      IpoptApplicationFactory();
  application->Options()->SetStringValue("option_file_name", "");
  application->Options()->SetIntegerValue("print_level", 0);
  application->Options()->SetStringValue("sb", "yes");
  application->Options()->SetIntegerValue("max_iter",
                                          controls.maximum_iterations);
  const double ipopt_tolerance = 0.001 * controls.optimizer_tolerance;
  application->Options()->SetNumericValue("tol", ipopt_tolerance);
  application->Options()->SetIntegerValue("acceptable_iter", 0);
  application->Options()->SetStringValue("hessian_approximation",
                                         "limited-memory");
  application->Options()->SetStringValue("mu_strategy", "adaptive");
  application->Options()->SetNumericValue("mu_min", 0.01 * ipopt_tolerance);
  application->Options()->SetStringValue(
      "nlp_scaling_method", full_space ? "gradient-based" : "none");
  application->Options()->SetNumericValue("bound_relax_factor", 0.0);
  application->Options()->SetStringValue("honor_original_bounds", "yes");
  application->Options()->SetNumericValue("bound_push", 1.0e-12);
  application->Options()->SetStringValue("check_derivatives_for_naninf", "yes");
  const auto initialize_status = application->Initialize();
  if (initialize_status != Ipopt::Solve_Succeeded) {
    solution = problem->solution();
    return initialize_status;
  }
  const auto status = application->OptimizeTNLP(nlp);
  const Ipopt::SmartPtr<Ipopt::SolveStatistics> statistics =
      application->Statistics();
  if (IsValid(statistics)) {
    iterations = statistics->IterationCount();
  }
  solution = problem->solution();
  return status;
}

StartResult solve_start(std::size_t index, PyObject *owner_evaluator,
                        PyObject *full_space_evaluator,
                        const FullSpaceSetup *full_space_setup,
                        std::size_t active_parameter_count,
                        const std::vector<std::size_t> &identity_coordinates,
                        const std::vector<Coordinate> &coordinates,
                        const std::vector<Observation> &observations,
                        const Controls &controls, const Start &start,
                        const std::vector<std::string> &failure_codes) {
  StartResult result;
  result.index = index;
  result.start_identity = start.identity;
  result.initial_affine = start.affine;
  result.final_affine = start.affine;
  std::vector<double> origins;
  std::vector<double> scales;
  origins.reserve(coordinates.size());
  scales.reserve(coordinates.size());
  for (const auto &coordinate : coordinates) {
    origins.push_back(coordinate.origin);
    scales.push_back(coordinate.scale);
  }
  std::vector<std::unique_ptr<OwnerCost>> costs;
  costs.reserve(observations.size());
  for (std::size_t observation_index = 0;
       observation_index < observations.size(); ++observation_index) {
    costs.push_back(std::make_unique<OwnerCost>(
        owner_evaluator, index, observation_index, active_parameter_count,
        identity_coordinates, origins, scales, observations[observation_index],
        failure_codes));
  }
  Ipopt::ApplicationReturnStatus solver_status;
  if (full_space_setup == nullptr) {
    auto *problem = new ReducedSpaceNlp(costs, coordinates, result.final_affine);
    solver_status =
        solve_ipopt(problem, controls, result.iterations, result.final_affine);
  } else {
    std::vector<double> initial = result.final_affine;
    std::vector<double> lower;
    std::vector<double> upper;
    lower.reserve(coordinates.size() + full_space_setup->lower.size());
    upper.reserve(coordinates.size() + full_space_setup->upper.size());
    for (const auto &coordinate : coordinates) {
      lower.push_back(coordinate.lower_affine);
      upper.push_back(coordinate.upper_affine);
    }
    initial.insert(initial.end(), full_space_setup->initial.begin(),
                   full_space_setup->initial.end());
    lower.insert(lower.end(), full_space_setup->lower.begin(),
                 full_space_setup->lower.end());
    upper.insert(upper.end(), full_space_setup->upper.begin(),
                 full_space_setup->upper.end());
    std::size_t residual_count = 0;
    for (const auto &observation : observations) {
      residual_count += observation.targets.size();
    }
    std::vector<double> solution;
    auto *problem = new SparseFullSpaceNlp(
        full_space_evaluator, index, std::move(initial), std::move(lower),
        std::move(upper), residual_count,
        full_space_setup->constraint_lower, full_space_setup->constraint_upper,
        full_space_setup->jacobian_rows,
        full_space_setup->jacobian_columns);
    solver_status = solve_ipopt(problem, controls, result.iterations, solution, true);
    result.final_affine.assign(solution.begin(),
                               solution.begin() + coordinates.size());
  }
  result.solver_status = ipopt_status_name(solver_status);
  const bool solver_usable = ipopt_solution_usable(solver_status);
  for (std::size_t coordinate = 0; coordinate < coordinates.size();
       ++coordinate) {
    double &affine = result.final_affine[coordinate];
    if (std::abs(affine - coordinates[coordinate].lower_affine) <=
        controls.affine_bound_margin) {
      affine = coordinates[coordinate].lower_affine;
    } else if (std::abs(affine - coordinates[coordinate].upper_affine) <=
               controls.affine_bound_margin) {
      affine = coordinates[coordinate].upper_affine;
    }
  }
  std::vector<double *> blocks;
  blocks.reserve(result.final_affine.size());
  for (double &value : result.final_affine) {
    blocks.push_back(&value);
  }
  std::size_t residual_count = 0;
  for (const auto &observation : observations) {
    residual_count += observation.targets.size();
  }
  result.predictions.assign(residual_count,
                            std::numeric_limits<double>::quiet_NaN());
  result.residuals.assign(residual_count,
                          std::numeric_limits<double>::quiet_NaN());
  result.affine_jacobian.assign(residual_count * coordinates.size(), 0.0);
  std::vector<bool> information_rows(residual_count, true);
  result.batch_outcomes.reserve(costs.size());
  for (std::size_t batch_index = 0; batch_index < costs.size(); ++batch_index) {
    result.batch_outcomes.push_back(
        {batch_index, costs[batch_index]->row_count(), "not_evaluated", ""});
  }
  bool censor_kink = false;
  std::vector<const double *> final_blocks;
  final_blocks.reserve(blocks.size());
  for (const double *block : blocks) {
    final_blocks.push_back(block);
  }
  std::size_t offset = 0;
  bool final_evaluable = true;
  const bool request_receipts = solver_usable;
  for (std::size_t batch_index = 0; batch_index < costs.size(); ++batch_index) {
    const auto &cost = costs[batch_index];
    const std::size_t rows = cost->row_count();
    std::vector<std::vector<double>> jacobian_blocks(coordinates.size(),
                                                     std::vector<double>(rows));
    std::vector<double *> jacobian_pointers;
    jacobian_pointers.reserve(coordinates.size());
    for (auto &block : jacobian_blocks) {
      jacobian_pointers.push_back(block.data());
    }
    std::string receipt_json;
    if (!cost->evaluate(final_blocks.data(), result.residuals.data() + offset,
                        jacobian_pointers.data(),
                        result.predictions.data() + offset, request_receipts,
                        &receipt_json)) {
      result.batch_outcomes[batch_index].status = "rejected";
      const auto &failures = cost->owner_failures();
      result.batch_outcomes[batch_index].failure_code =
          failures.empty() ? "owner_contract_failure" : failures.back().code;
      final_evaluable = false;
      break;
    }
    result.batch_outcomes[batch_index].status = "evaluated";
    if (!receipt_json.empty()) {
      result.owner_receipts_json.push_back(std::move(receipt_json));
    }
    for (std::size_t row = 0; row < rows; ++row) {
      information_rows[offset + row] =
          cost->information_bearing(row, result.predictions[offset + row]);
      censor_kink = censor_kink ||
                    cost->censor_kink(row, result.predictions[offset + row]);
      for (std::size_t coordinate = 0; coordinate < coordinates.size();
           ++coordinate) {
        result
            .affine_jacobian[(offset + row) * coordinates.size() + coordinate] =
            jacobian_blocks[coordinate][row];
      }
    }
    offset += rows;
  }
  result.usable = solver_usable && final_evaluable;
  if (!result.usable) {
    result.owner_receipts_json.clear();
  }
  result.non_evaluable_count = 0;
  for (const auto &cost : costs) {
    result.non_evaluable_count += cost->rejection_count();
    for (const auto &failure : cost->owner_failures()) {
      const auto match = std::find_if(
          result.owner_failures.begin(), result.owner_failures.end(),
          [&](const OwnerFailure &existing) {
            return existing.code == failure.code &&
                   existing.reason == failure.reason;
          });
      if (match == result.owner_failures.end()) {
        result.owner_failures.push_back(failure);
      } else {
        match->count += failure.count;
      }
    }
  }
  if (!result.owner_failures.empty()) {
    result.failure_reason = result.owner_failures.front().reason;
    result.trial_failure_code = result.owner_failures.front().code;
  }
  if (result.usable) {
    result.cost = 0.5 * std::inner_product(result.residuals.begin(),
                                           result.residuals.end(),
                                           result.residuals.begin(), 0.0);
    result.failure_reason.clear();
    for (std::size_t coordinate = 0; coordinate < coordinates.size();
         ++coordinate) {
      const double affine = result.final_affine[coordinate];
      if (std::abs(affine - coordinates[coordinate].lower_affine) <=
              controls.affine_bound_margin ||
          std::abs(affine - coordinates[coordinate].upper_affine) <=
              controls.affine_bound_margin) {
        result.active_bounds.push_back(coordinate);
      }
    }
    result.diagnostics =
        diagnose(result.affine_jacobian, result.residuals, residual_count,
                 coordinates.size(), result.active_bounds, information_rows,
                 censor_kink);
    result.trial_failure_code.clear();
  } else {
    if (result.failure_reason.empty()) {
      result.failure_reason = result.solver_status;
    }
    if (result.trial_failure_code.empty()) {
      result.trial_failure_code = "optimizer_failure";
    }
    result.cost = std::numeric_limits<double>::infinity();
  }
  return result;
}

void assign_optimizer_basins(std::vector<StartResult> &results,
                             const std::vector<Coordinate> &coordinates) {
  std::vector<std::size_t> representatives;
  for (auto &result : results) {
    if (!result.usable) {
      continue;
    }
    for (std::size_t basin = 0; basin < representatives.size(); ++basin) {
      const auto &other = results[representatives[basin]];
      bool same = true;
      for (std::size_t coordinate = 0; coordinate < coordinates.size();
           ++coordinate) {
        const double span = coordinates[coordinate].upper_affine -
                            coordinates[coordinate].lower_affine;
        const double distance = std::abs(result.final_affine[coordinate] -
                                         other.final_affine[coordinate]) /
                                span;
        same = same && distance <= kOptimizerBasinAffineDistanceMax;
      }
      const double cost_scale =
          std::max({1.0, std::abs(result.cost), std::abs(other.cost)});
      if (same && std::abs(result.cost - other.cost) <=
                      kOptimizerBasinCostDifferenceFactor * cost_scale) {
        result.optimizer_basin = basin;
        break;
      }
    }
    if (result.optimizer_basin == std::numeric_limits<std::size_t>::max()) {
      result.optimizer_basin = representatives.size();
      representatives.push_back(result.index);
    }
  }
}

py::dict diagnostics_mapping(const Diagnostics &value) {
  py::dict result;
  result["singular_values"] = value.singular_values;
  result["rank"] = value.rank;
  result["free_subspace_rank"] = value.free_subspace_rank;
  result["condition_number"] = value.condition_number;
  result["affine_covariance"] = value.affine_covariance;
  result["degrees_of_freedom"] = value.degrees_of_freedom;
  result["information_row_count"] = value.information_row_count;
  result["free_coordinate_indices"] = value.free_coordinate_indices;
  result["residual_variance"] = value.residual_variance;
  result["covariance_status"] = value.covariance_status;
  return result;
}

py::dict owner_failure_mapping(const OwnerFailure &value) {
  py::dict result;
  result["code"] = value.code;
  result["reason"] = value.reason;
  result["count"] = value.count;
  return result;
}

py::dict batch_outcome_mapping(const BatchOutcome &value) {
  py::dict result;
  result["index"] = value.index;
  result["row_count"] = value.row_count;
  result["status"] = value.status;
  result["failure_code"] = value.failure_code;
  return result;
}

py::dict result_mapping(const StartResult &value,
                        const std::vector<Coordinate> &coordinates) {
  py::dict result;
  std::vector<double> initial_physical;
  std::vector<double> final_physical;
  std::vector<double> physical_jacobian = value.affine_jacobian;
  for (std::size_t coordinate = 0; coordinate < coordinates.size();
       ++coordinate) {
    initial_physical.push_back(coordinates[coordinate].origin +
                               value.initial_affine[coordinate] *
                                   coordinates[coordinate].scale);
    final_physical.push_back(coordinates[coordinate].origin +
                             value.final_affine[coordinate] *
                                 coordinates[coordinate].scale);
    for (std::size_t row = 0; row < value.residuals.size(); ++row) {
      physical_jacobian[row * coordinates.size() + coordinate] /=
          coordinates[coordinate].scale;
    }
  }
  result["index"] = value.index;
  result["start_identity"] = value.start_identity;
  result["initial_affine"] = value.initial_affine;
  result["initial_physical"] = initial_physical;
  result["final_physical"] = final_physical;
  result["final_affine"] = value.final_affine;
  result["predictions"] = value.predictions;
  result["residuals"] = value.residuals;
  result["physical_jacobian"] = physical_jacobian;
  result["affine_jacobian"] = value.affine_jacobian;
  result["diagnostics"] = diagnostics_mapping(value.diagnostics);
  result["active_bounds"] = value.active_bounds;
  result["solver_status"] = value.solver_status;
  result["failure_reason"] = value.failure_reason;
  result["trial_failure_code"] = value.trial_failure_code;
  result["iterations"] = value.iterations;
  result["non_evaluable_count"] = value.non_evaluable_count;
  py::list owner_failures;
  for (const auto &failure : value.owner_failures) {
    owner_failures.append(owner_failure_mapping(failure));
  }
  result["owner_failures"] = owner_failures;
  result["owner_receipts_json"] = value.owner_receipts_json;
  py::list batch_outcomes;
  for (const auto &outcome : value.batch_outcomes) {
    batch_outcomes.append(batch_outcome_mapping(outcome));
  }
  result["batch_outcomes"] = batch_outcomes;
  result["cost"] = value.cost;
  result["usable"] = value.usable;
  if (value.usable) {
    result["optimizer_basin"] = value.optimizer_basin;
  } else {
    result["optimizer_basin"] = py::none();
  }
  return result;
}

py::dict fit_resolved(const py::function &owner_evaluator,
                      std::size_t active_parameter_count,
                      const std::vector<std::size_t> &identity_coordinates,
                      const py::iterable &coordinate_values,
                      const py::iterable &observation_values,
                      const py::dict &control_values,
                      const std::vector<std::string> &failure_codes,
                      const py::object &full_space_evaluator,
                      const py::iterable &full_space_values) {
  const auto coordinates = parse_coordinates(coordinate_values);
  const auto observations = parse_observations(observation_values);
  const auto controls = parse_controls(control_values);
  const auto full_space_setups = parse_full_space_setups(full_space_values);
  if (active_parameter_count == 0 ||
      active_parameter_count != identity_coordinates.size() ||
      std::any_of(identity_coordinates.begin(), identity_coordinates.end(),
                  [&](std::size_t coordinate) {
                    return coordinate >= coordinates.size();
                  })) {
    throw std::invalid_argument("identity coordinates must bind every active "
                                "identity to one coordinate");
  }
  if (controls.starts.size() == 0) {
    throw std::invalid_argument("at least one explicit start is required");
  }
  if (failure_codes.empty()) {
    throw std::invalid_argument(
        "the owner failure vocabulary must be nonempty");
  }
  const bool full_space = !full_space_evaluator.is_none();
  if (full_space != !full_space_setups.empty() ||
      (full_space &&
       (full_space_setups.size() != controls.starts.size() ||
        !PyCallable_Check(full_space_evaluator.ptr())))) {
    throw std::invalid_argument(
        "full-space evaluator and one setup per start must be supplied together");
  }
  for (const auto &setup : full_space_setups) {
    const std::size_t variable_count = coordinates.size() + setup.initial.size();
    for (std::size_t index = 0; index < setup.jacobian_rows.size(); ++index) {
      if (setup.jacobian_rows[index] < 0 ||
          setup.jacobian_rows[index] >=
              static_cast<int>(setup.constraint_lower.size()) ||
          setup.jacobian_columns[index] < 0 ||
          setup.jacobian_columns[index] >= static_cast<int>(variable_count)) {
        throw std::invalid_argument("full-space sparse index is invalid");
      }
    }
  }
  std::vector<std::string> start_identities;
  start_identities.reserve(controls.starts.size());
  for (const auto &start : controls.starts) {
    if (start.identity.empty() || start.affine.size() != coordinates.size() ||
        !std::all_of(start.affine.begin(), start.affine.end(),
                     [](double value) { return std::isfinite(value); })) {
      throw std::invalid_argument("each explicit start needs an identity and "
                                  "one finite affine value per coordinate");
    }
    for (std::size_t coordinate = 0; coordinate < coordinates.size();
         ++coordinate) {
      if (start.affine[coordinate] < coordinates[coordinate].lower_affine ||
          start.affine[coordinate] > coordinates[coordinate].upper_affine) {
        throw std::invalid_argument(
            "explicit affine start lies outside its bounds");
      }
    }
    start_identities.push_back(start.identity);
  }
  std::sort(start_identities.begin(), start_identities.end());
  if (std::adjacent_find(start_identities.begin(), start_identities.end()) !=
      start_identities.end()) {
    throw std::invalid_argument("explicit start identities must be unique");
  }
  const auto &start_values = controls.starts;
  std::vector<StartResult> results(start_values.size());
  {
    py::gil_scoped_release release;
    for (std::size_t index = 0; index < start_values.size(); ++index) {
      try {
        results[index] = solve_start(
            index, owner_evaluator.ptr(),
            full_space ? full_space_evaluator.ptr() : nullptr,
            full_space ? &full_space_setups[index] : nullptr,
            active_parameter_count, identity_coordinates, coordinates,
            observations, controls, start_values[index], failure_codes);
      } catch (const std::exception &error) {
        results[index].index = index;
        results[index].start_identity = start_values[index].identity;
        results[index].initial_affine = start_values[index].affine;
        results[index].final_affine = start_values[index].affine;
        results[index].failure_reason = error.what();
        results[index].trial_failure_code = "optimizer_failure";
        for (std::size_t batch = 0; batch < observations.size(); ++batch) {
          results[index].batch_outcomes.push_back(
              {batch, observations[batch].targets.size(), "not_evaluated",
               "optimizer_failure"});
        }
      }
    }
  }
  assign_optimizer_basins(results, coordinates);
  py::list starts_result;
  for (const auto &result : results) {
    starts_result.append(result_mapping(result, coordinates));
  }
  py::dict output;
  output["starts"] = starts_result;
  return output;
}

} // namespace
} // namespace epcsaft::regression

void bind_parameter_fit(py::module_ &module) {
  module.def("_fit_resolved", &epcsaft::regression::fit_resolved,
             py::arg("owner_evaluator"), py::arg("active_parameter_count"),
             py::arg("identity_coordinates"), py::arg("coordinates"),
             py::arg("observations"), py::arg("controls"),
             py::arg("failure_codes"),
             py::arg("full_space_evaluator") = py::none(),
             py::arg("full_space_setups") = py::tuple());
}
