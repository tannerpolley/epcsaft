#include "fit.hpp"

#include <Eigen/Dense>
#include <Eigen/SVD>
#include <ceres/ceres.h>
#include <pybind11/stl.h>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <limits>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>
#include <vector>

namespace py = pybind11;

namespace epcsaft::regression {
namespace {

enum class ResidualKind { LogRatio, ScaledDifference };

struct Target final {
  double observed;
  ResidualKind residual;
  double scale;
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
};

struct Coordinate final {
  double scale;
  double lower;
  double upper;
  double input;
};

struct Controls final {
  std::size_t start_count;
  bool include_input;
  std::size_t workers;
  int maximum_iterations;
  double function_tolerance;
  double gradient_tolerance;
  double parameter_tolerance;
};

struct Diagnostics final {
  std::vector<double> singular_values;
  int rank{0};
  double condition_number{std::numeric_limits<double>::infinity()};
};

struct StartResult final {
  std::size_t index{0};
  std::vector<double> initial_affine;
  std::vector<double> final_affine;
  std::vector<double> predictions;
  std::vector<double> residuals;
  std::vector<double> affine_jacobian;
  std::vector<OwnerFailure> owner_failures;
  Diagnostics diagnostics;
  std::vector<std::size_t> active_bounds;
  std::string solver_status;
  std::string failure_reason;
  int iterations{0};
  std::size_t non_evaluable_count{0};
  double cost{std::numeric_limits<double>::infinity()};
  bool usable{false};
  std::size_t basin{std::numeric_limits<std::size_t>::max()};
};

ResidualKind residual_kind(const std::string &value) {
  if (value == "log_ratio") {
    return ResidualKind::LogRatio;
  }
  if (value == "scaled_difference") {
    return ResidualKind::ScaledDifference;
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
      throw std::invalid_argument("observation multiplier must be positive and finite");
    }
    for (const py::handle target_item :
         required<py::iterable>(mapping, "targets")) {
      const auto target_mapping = py::cast<py::dict>(target_item);
      Target target{
          required<double>(target_mapping, "observed"),
          residual_kind(required<std::string>(target_mapping, "residual")),
          required<double>(target_mapping, "scale"),
      };
      if (!(std::isfinite(target.observed) && std::isfinite(target.scale) &&
            target.scale > 0.0)) {
        throw std::invalid_argument("target values and scales must be finite");
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
        required<double>(mapping, "scale"),
        required<double>(mapping, "lower"),
        required<double>(mapping, "upper"),
        required<double>(mapping, "input"),
    };
    if (!(std::isfinite(coordinate.scale) && coordinate.scale > 0.0 &&
          std::isfinite(coordinate.lower) && std::isfinite(coordinate.upper) &&
          coordinate.lower < coordinate.upper &&
          std::isfinite(coordinate.input) &&
          coordinate.input >= coordinate.lower &&
          coordinate.input <= coordinate.upper)) {
      throw std::invalid_argument(
          "coordinate scale, bounds, and input must be finite and consistent");
    }
    coordinates.push_back(coordinate);
  }
  if (coordinates.empty()) {
    throw std::invalid_argument("coordinates must be nonempty");
  }
  return coordinates;
}

Controls parse_controls(const py::dict &mapping) {
  Controls controls{
      required<std::size_t>(mapping, "start_count"),
      required<bool>(mapping, "include_input"),
      required<std::size_t>(mapping, "workers"),
      required<int>(mapping, "maximum_iterations"),
      required<double>(mapping, "function_tolerance"),
      required<double>(mapping, "gradient_tolerance"),
      required<double>(mapping, "parameter_tolerance"),
  };
  if (controls.start_count == 0 || controls.workers == 0 ||
      controls.maximum_iterations <= 0 ||
      !(std::isfinite(controls.function_tolerance) &&
        controls.function_tolerance > 0.0 &&
        std::isfinite(controls.gradient_tolerance) &&
        controls.gradient_tolerance > 0.0 &&
        std::isfinite(controls.parameter_tolerance) &&
        controls.parameter_tolerance > 0.0)) {
    throw std::invalid_argument("solver controls must be positive and finite");
  }
  return controls;
}

class OwnerEvaluationError final : public std::runtime_error {
public:
  OwnerEvaluationError(std::string code, std::string reason)
      : std::runtime_error(std::move(reason)), code_(std::move(code)) {}

  const std::string &code() const noexcept { return code_; }

private:
  std::string code_;
};

class OwnerCost final : public ceres::CostFunction {
public:
  OwnerCost(PyObject *owner_evaluator, std::size_t observation_index,
          std::size_t active_parameter_count,
          std::vector<std::size_t> identity_coordinates,
          std::vector<double> coordinate_scales, Observation observation)
      : owner_evaluator_(owner_evaluator), observation_index_(observation_index),
        active_parameter_count_(active_parameter_count),
        identity_coordinates_(std::move(identity_coordinates)),
        coordinate_scales_(std::move(coordinate_scales)),
        observation_(std::move(observation)) {
    set_num_residuals(static_cast<int>(observation_.targets.size()));
    mutable_parameter_block_sizes()->assign(coordinate_scales_.size(), 1);
  }

  bool Evaluate(double const *const *parameters, double *residuals,
                double **jacobians) const override {
    return evaluate(parameters, residuals, jacobians, nullptr);
  }

  bool evaluate(double const *const *parameters, double *residuals,
                double **jacobians, double *predictions) const {
    try {
      std::vector<double> active(active_parameter_count_);
      for (std::size_t index = 0; index < active_parameter_count_; ++index) {
        const std::size_t coordinate = identity_coordinates_[index];
        active[index] =
            parameters[coordinate][0] * coordinate_scales_[coordinate];
      }
      const OwnerEvaluation evaluation = evaluate_owner(active);
      for (std::size_t row = 0; row < observation_.targets.size(); ++row) {
        const auto &target = observation_.targets[row];
        const double predicted = evaluation.values[row];
        if (!std::isfinite(predicted) ||
            (target.residual == ResidualKind::LogRatio && predicted <= 0.0)) {
          return reject("owner-inadmissible-value",
                        "owner returned an inadmissible target value");
        }
        if (predictions != nullptr) {
          predictions[row] = predicted;
        }
        const double factor = target.residual == ResidualKind::LogRatio
                                  ? observation_.multiplier / predicted
                                  : observation_.multiplier / target.scale;
        residuals[row] = target.residual == ResidualKind::LogRatio
                             ? observation_.multiplier *
                                   std::log(predicted / target.observed)
                             : observation_.multiplier *
                                   (predicted - target.observed) / target.scale;
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
                  evaluation.derivatives[row * active_parameter_count_ + identity];
            }
          }
          jacobians[coordinate][row] = factor * derivative *
                                       coordinate_scales_[coordinate];
        }
      }
      if (!std::all_of(residuals, residuals + observation_.targets.size(),
                       [](double value) { return std::isfinite(value); })) {
        return reject("nonfinite-residual", "residual evaluation was nonfinite");
      }
      return true;
    } catch (const OwnerEvaluationError &error) {
      return reject(error.code(), error.what());
    } catch (const std::exception &error) {
      return reject("owner-evaluation-error", error.what());
    }
  }

  std::size_t rejection_count() const noexcept {
    return std::accumulate(
        owner_failures_.begin(), owner_failures_.end(), std::size_t{0},
        [](std::size_t count, const OwnerFailure &failure) {
          return count + failure.count;
        });
  }
  const std::vector<OwnerFailure> &owner_failures() const noexcept {
    return owner_failures_;
  }

private:
  OwnerEvaluation evaluate_owner(const std::vector<double> &active) const {
    py::gil_scoped_acquire acquire;
    const auto evaluator = py::reinterpret_borrow<py::function>(
        py::handle(owner_evaluator_));
    const py::dict result = py::cast<py::dict>(
        evaluator(observation_index_, active));
    const std::string status = required<std::string>(result, "status");
    if (status == "rejected") {
      const auto failure = required<py::dict>(result, "failure");
      const std::string code = required<std::string>(failure, "code");
      const std::string reason = required<std::string>(failure, "reason");
      throw OwnerEvaluationError(code, reason);
    }
    if (status != "accepted") {
      throw OwnerEvaluationError("owner-status", "unknown owner evaluation status");
    }
    OwnerEvaluation evaluation{
        required<std::vector<double>>(result, "predictions"),
        required<std::vector<double>>(result, "jacobian"),
    };
    if (evaluation.values.size() != observation_.targets.size() ||
        evaluation.derivatives.size() !=
            observation_.targets.size() * active_parameter_count_) {
      throw OwnerEvaluationError(
          "owner-shape-mismatch", "owner block has the wrong value or Jacobian shape");
    }
    return evaluation;
  }

  bool reject(const std::string &code, const std::string &reason) const {
    const auto match = std::find_if(
        owner_failures_.begin(), owner_failures_.end(),
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
  std::size_t observation_index_;
  std::size_t active_parameter_count_;
  std::vector<std::size_t> identity_coordinates_;
  std::vector<double> coordinate_scales_;
  Observation observation_;
  mutable std::vector<OwnerFailure> owner_failures_;
};

double radical_inverse(std::size_t index, std::size_t base) {
  double result = 0.0;
  double factor = 1.0 / static_cast<double>(base);
  while (index > 0) {
    result += factor * static_cast<double>(index % base);
    index /= base;
    factor /= static_cast<double>(base);
  }
  return result;
}

std::size_t prime(std::size_t ordinal) {
  std::size_t candidate = 2;
  for (std::size_t found = 0;; ++candidate) {
    bool is_prime = true;
    for (std::size_t divisor = 2; divisor * divisor <= candidate; ++divisor) {
      if (candidate % divisor == 0) {
        is_prime = false;
        break;
      }
    }
    if (is_prime && found++ == ordinal) {
      return candidate;
    }
  }
}

std::vector<std::vector<double>>
starts(const std::vector<Coordinate> &coordinates, const Controls &controls) {
  std::vector<std::vector<double>> result;
  result.reserve(controls.start_count);
  if (controls.include_input) {
    std::vector<double> input;
    input.reserve(coordinates.size());
    for (const auto &coordinate : coordinates) {
      input.push_back(coordinate.input / coordinate.scale);
    }
    result.push_back(std::move(input));
  }
  for (std::size_t sequence = 1; result.size() < controls.start_count;
       ++sequence) {
    std::vector<double> point;
    point.reserve(coordinates.size());
    for (std::size_t coordinate = 0; coordinate < coordinates.size();
         ++coordinate) {
      const double fraction = radical_inverse(sequence, prime(coordinate));
      const double physical = coordinates[coordinate].lower +
                              fraction * (coordinates[coordinate].upper -
                                          coordinates[coordinate].lower);
      point.push_back(physical / coordinates[coordinate].scale);
    }
    result.push_back(std::move(point));
  }
  return result;
}

Diagnostics diagnose(const std::vector<double> &jacobian, std::size_t rows,
                     std::size_t columns) {
  using RowMajor =
      Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
  const Eigen::Map<const RowMajor> matrix(jacobian.data(),
                                          static_cast<Eigen::Index>(rows),
                                          static_cast<Eigen::Index>(columns));
  const Eigen::JacobiSVD<Eigen::MatrixXd> svd(matrix, Eigen::ComputeThinU |
                                                          Eigen::ComputeThinV);
  Diagnostics result;
  result.singular_values.assign(svd.singularValues().data(),
                                svd.singularValues().data() +
                                    svd.singularValues().size());
  if (result.singular_values.empty()) {
    result.condition_number = 1.0;
    return result;
  }
  const double threshold = 100.0 * std::numeric_limits<double>::epsilon() *
                           static_cast<double>(std::max(rows, columns)) *
                           result.singular_values.front();
  result.rank = static_cast<int>(std::count_if(
      result.singular_values.begin(), result.singular_values.end(),
      [threshold](double value) { return value > threshold; }));
  if (result.rank == static_cast<int>(columns)) {
    result.condition_number =
        result.singular_values.front() /
        result.singular_values[static_cast<std::size_t>(result.rank - 1)];
  }
  return result;
}

ceres::Solver::Options solver_options(const Controls &controls) {
  ceres::Solver::Options options;
  options.linear_solver_type = ceres::DENSE_QR;
  options.max_num_iterations = controls.maximum_iterations;
  options.function_tolerance = controls.function_tolerance;
  options.gradient_tolerance = controls.gradient_tolerance;
  options.parameter_tolerance = controls.parameter_tolerance;
  options.logging_type = ceres::SILENT;
  options.num_threads = 1;
  return options;
}

StartResult solve_start(std::size_t index, PyObject *owner_evaluator,
                        std::size_t active_parameter_count,
                        const std::vector<std::size_t> &identity_coordinates,
                        const std::vector<Coordinate> &coordinates,
                        const std::vector<Observation> &observations,
                        const Controls &controls,
                        const std::vector<double> &initial) {
  StartResult result;
  result.index = index;
  result.initial_affine = initial;
  result.final_affine = initial;
  std::vector<double> scales;
  scales.reserve(coordinates.size());
  for (const auto &coordinate : coordinates) {
    scales.push_back(coordinate.scale);
  }
  std::vector<std::unique_ptr<OwnerCost>> costs;
  costs.reserve(observations.size());
  ceres::Problem::Options problem_options;
  problem_options.cost_function_ownership = ceres::DO_NOT_TAKE_OWNERSHIP;
  ceres::Problem problem(problem_options);
  std::vector<double *> blocks;
  blocks.reserve(result.final_affine.size());
  for (double &value : result.final_affine) {
    blocks.push_back(&value);
  }
  for (std::size_t observation_index = 0;
       observation_index < observations.size(); ++observation_index) {
    costs.push_back(std::make_unique<OwnerCost>(
        owner_evaluator, observation_index, active_parameter_count,
        identity_coordinates, scales, observations[observation_index]));
    problem.AddResidualBlock(costs.back().get(), nullptr, blocks);
  }
  for (std::size_t coordinate = 0; coordinate < coordinates.size();
       ++coordinate) {
    problem.SetParameterLowerBound(blocks[coordinate], 0,
                                   coordinates[coordinate].lower /
                                       coordinates[coordinate].scale);
    problem.SetParameterUpperBound(blocks[coordinate], 0,
                                   coordinates[coordinate].upper /
                                       coordinates[coordinate].scale);
  }
  ceres::Solver::Summary summary;
  ceres::Solve(solver_options(controls), &problem, &summary);
  result.solver_status = summary.BriefReport();
  result.iterations = static_cast<int>(summary.iterations.size());
  result.cost = summary.final_cost;
  std::size_t residual_count = 0;
  for (const auto &observation : observations) {
    residual_count += observation.targets.size();
  }
  result.predictions.assign(residual_count,
                            std::numeric_limits<double>::quiet_NaN());
  result.residuals.assign(residual_count,
                          std::numeric_limits<double>::quiet_NaN());
  result.affine_jacobian.assign(residual_count * coordinates.size(), 0.0);
  std::vector<const double *> final_blocks;
  final_blocks.reserve(blocks.size());
  for (const double *block : blocks) {
    final_blocks.push_back(block);
  }
  std::size_t offset = 0;
  bool final_evaluable = true;
  for (const auto &cost : costs) {
    const std::size_t rows = static_cast<std::size_t>(cost->num_residuals());
    std::vector<std::vector<double>> jacobian_blocks(coordinates.size(),
                                                     std::vector<double>(rows));
    std::vector<double *> jacobian_pointers;
    jacobian_pointers.reserve(coordinates.size());
    for (auto &block : jacobian_blocks) {
      jacobian_pointers.push_back(block.data());
    }
    if (!cost->evaluate(final_blocks.data(), result.residuals.data() + offset,
                        jacobian_pointers.data(),
                        result.predictions.data() + offset)) {
      final_evaluable = false;
      break;
    }
    for (std::size_t row = 0; row < rows; ++row) {
      for (std::size_t coordinate = 0; coordinate < coordinates.size();
           ++coordinate) {
        result
            .affine_jacobian[(offset + row) * coordinates.size() + coordinate] =
            jacobian_blocks[coordinate][row];
      }
    }
    offset += rows;
  }
  result.usable = summary.IsSolutionUsable() && final_evaluable;
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
  }
  if (result.usable) {
    result.failure_reason.clear();
    result.diagnostics =
        diagnose(result.affine_jacobian, residual_count, coordinates.size());
    for (std::size_t coordinate = 0; coordinate < coordinates.size();
         ++coordinate) {
      const double physical =
          result.final_affine[coordinate] * coordinates[coordinate].scale;
      const double tolerance =
          1.0e-10 * std::max(1.0, coordinates[coordinate].upper -
                                      coordinates[coordinate].lower);
      if (std::abs(physical - coordinates[coordinate].lower) <= tolerance ||
          std::abs(physical - coordinates[coordinate].upper) <= tolerance) {
        result.active_bounds.push_back(coordinate);
      }
    }
  } else {
    if (result.failure_reason.empty()) {
      result.failure_reason = summary.message;
    }
    result.cost = std::numeric_limits<double>::infinity();
  }
  return result;
}

void assign_basins(std::vector<StartResult> &results,
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
        const double span =
            coordinates[coordinate].upper - coordinates[coordinate].lower;
        const double distance = std::abs(result.final_affine[coordinate] -
                                         other.final_affine[coordinate]) *
                                coordinates[coordinate].scale / span;
        same = same && distance <= 1.0e-6;
      }
      const double cost_scale =
          std::max({1.0, std::abs(result.cost), std::abs(other.cost)});
      if (same && std::abs(result.cost - other.cost) <= 1.0e-8 * cost_scale) {
        result.basin = basin;
        break;
      }
    }
    if (result.basin == std::numeric_limits<std::size_t>::max()) {
      result.basin = representatives.size();
      representatives.push_back(result.index);
    }
  }
}

py::dict diagnostics_mapping(const Diagnostics &value) {
  py::dict result;
  result["singular_values"] = value.singular_values;
  result["rank"] = value.rank;
  result["condition_number"] = value.condition_number;
  return result;
}

py::dict owner_failure_mapping(const OwnerFailure &value) {
  py::dict result;
  result["code"] = value.code;
  result["reason"] = value.reason;
  result["count"] = value.count;
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
    initial_physical.push_back(value.initial_affine[coordinate] *
                               coordinates[coordinate].scale);
    final_physical.push_back(value.final_affine[coordinate] *
                             coordinates[coordinate].scale);
    for (std::size_t row = 0; row < value.residuals.size(); ++row) {
      physical_jacobian[row * coordinates.size() + coordinate] /=
          coordinates[coordinate].scale;
    }
  }
  result["index"] = value.index;
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
  result["iterations"] = value.iterations;
  result["non_evaluable_count"] = value.non_evaluable_count;
  py::list owner_failures;
  for (const auto &failure : value.owner_failures) {
    owner_failures.append(owner_failure_mapping(failure));
  }
  result["owner_failures"] = owner_failures;
  result["cost"] = value.cost;
  result["usable"] = value.usable;
  if (value.usable) {
    result["basin"] = value.basin;
  } else {
    result["basin"] = py::none();
  }
  return result;
}

py::dict fit_resolved(const py::function &owner_evaluator,
                      std::size_t active_parameter_count,
                      const std::vector<std::size_t> &identity_coordinates,
                      const py::iterable &coordinate_values,
                      const py::iterable &observation_values,
                      const py::dict &control_values) {
  const auto coordinates = parse_coordinates(coordinate_values);
  const auto observations = parse_observations(observation_values);
  const auto controls = parse_controls(control_values);
  if (active_parameter_count == 0 ||
      active_parameter_count != identity_coordinates.size() ||
      std::any_of(identity_coordinates.begin(), identity_coordinates.end(),
                  [&](std::size_t coordinate) {
                    return coordinate >= coordinates.size();
                  })) {
    throw std::invalid_argument("identity coordinates must bind every active "
                                "identity to one coordinate");
  }
  const auto start_values = starts(coordinates, controls);
  std::vector<StartResult> results(start_values.size());
  std::atomic<std::size_t> next{0};
  const std::size_t worker_count =
      std::min(controls.workers, start_values.size());
  {
    py::gil_scoped_release release;
    auto solve_next = [&]() {
      for (;;) {
        const std::size_t index = next.fetch_add(1);
        if (index >= start_values.size()) {
          return;
        }
        try {
          results[index] = solve_start(
              index, owner_evaluator.ptr(), active_parameter_count,
              identity_coordinates, coordinates, observations, controls,
              start_values[index]);
        } catch (const std::exception &error) {
          results[index].index = index;
          results[index].initial_affine = start_values[index];
          results[index].final_affine = start_values[index];
          results[index].failure_reason = error.what();
        }
      }
    };
    std::vector<std::thread> workers;
    workers.reserve(worker_count);
    for (std::size_t worker = 0; worker < worker_count; ++worker) {
      workers.emplace_back(solve_next);
    }
    for (auto &worker : workers) {
      worker.join();
    }
  }
  assign_basins(results, coordinates);
  py::list starts_result;
  for (const auto &result : results) {
    starts_result.append(result_mapping(result, coordinates));
  }
  py::dict output;
  output["starts"] = starts_result;
  output["coordinate_count"] = coordinates.size();
  output["residual_count"] =
      std::accumulate(observations.begin(), observations.end(), std::size_t{0},
                      [](std::size_t count, const Observation &observation) {
                        return count + observation.targets.size();
                      });
  output["owner_callbacks_serialized"] = true;
  output["start_order"] = "halton";
  return output;
}

} // namespace
} // namespace epcsaft::regression

void bind_parameter_fit(py::module_ &module) {
  module.def("_fit_resolved", &epcsaft::regression::fit_resolved,
             py::arg("owner_evaluator"), py::arg("active_parameter_count"),
             py::arg("identity_coordinates"), py::arg("coordinates"),
             py::arg("observations"), py::arg("controls"));
}
