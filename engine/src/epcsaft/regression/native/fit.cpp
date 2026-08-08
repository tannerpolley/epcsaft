#include "fit.hpp"

#include "eos.hpp"
#include "phase_block.hpp"
#include "pure_saturation.hpp"

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
#include <mutex>
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
enum class ObservationKind { PureSaturation, FixedStatePressure };

struct Target final {
  std::size_t value_index;
  double observed;
  ResidualKind residual;
  double scale;
  double value_factor;
};

struct Observation final {
  ObservationKind kind;
  double temperature_k;
  double liquid_density_upper_mol_m3;
  double molar_density_mol_m3;
  std::vector<double> composition;
  double multiplier;
  std::vector<Target> targets;
};

struct EosEvaluation final {
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
  std::vector<double> residuals;
  std::vector<double> affine_jacobian;
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

std::mutex &saturation_solver_mutex() {
  static std::mutex mutex;
  return mutex;
}

ResidualKind residual_kind(const std::string &value) {
  if (value == "log_ratio") {
    return ResidualKind::LogRatio;
  }
  if (value == "scaled_difference") {
    return ResidualKind::ScaledDifference;
  }
  throw std::invalid_argument("unknown residual kind: " + value);
}

ObservationKind observation_kind(const std::string &value) {
  if (value == "pure_saturation") {
    return ObservationKind::PureSaturation;
  }
  if (value == "fixed_state_pressure") {
    return ObservationKind::FixedStatePressure;
  }
  throw std::invalid_argument("unknown observation family: " + value);
}

std::size_t target_index(ObservationKind kind, const std::string &value) {
  if (kind == ObservationKind::FixedStatePressure) {
    if (value == "pressure_pa") {
      return 0;
    }
    throw std::invalid_argument("unknown fixed-state pressure target: " +
                                value);
  }
  if (value == "saturation_pressure_pa") {
    return 0;
  }
  if (value == "vapor_molar_density_mol_per_m3") {
    return 1;
  }
  if (value == "liquid_molar_density_mol_per_m3") {
    return 2;
  }
  throw std::invalid_argument("unknown pure-saturation target: " + value);
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
    const ObservationKind kind =
        observation_kind(required<std::string>(mapping, "family"));
    Observation observation{
        kind,
        required<double>(mapping, "temperature_k"),
        std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::quiet_NaN(),
        {},
        required<double>(mapping, "multiplier"),
        {},
    };
    if (!(std::isfinite(observation.temperature_k) &&
          observation.temperature_k > 0.0 &&
          std::isfinite(observation.multiplier) &&
          observation.multiplier > 0.0)) {
      throw std::invalid_argument(
          "observation temperature and multiplier must be positive and finite");
    }
    if (kind == ObservationKind::PureSaturation) {
      observation.liquid_density_upper_mol_m3 =
          required<double>(mapping, "liquid_density_upper_mol_m3");
      if (!(std::isfinite(observation.liquid_density_upper_mol_m3) &&
            observation.liquid_density_upper_mol_m3 > 0.0)) {
        throw std::invalid_argument(
            "pure-saturation density bound must be positive and finite");
      }
    } else {
      observation.molar_density_mol_m3 =
          required<double>(mapping, "molar_density_mol_m3");
      observation.composition =
          required<std::vector<double>>(mapping, "composition");
      const double total = std::accumulate(observation.composition.begin(),
                                           observation.composition.end(), 0.0);
      if (!(std::isfinite(observation.molar_density_mol_m3) &&
            observation.molar_density_mol_m3 > 0.0 &&
            !observation.composition.empty() &&
            std::all_of(observation.composition.begin(),
                        observation.composition.end(),
                        [](double value) {
                          return std::isfinite(value) && value > 0.0;
                        }) &&
            std::abs(total - 1.0) <= 1.0e-12)) {
        throw std::invalid_argument(
            "fixed-state density and normalized composition are invalid");
      }
    }
    for (const py::handle target_item :
         required<py::iterable>(mapping, "targets")) {
      const auto target_mapping = py::cast<py::dict>(target_item);
      Target target{
          target_index(kind, required<std::string>(target_mapping, "label")),
          required<double>(target_mapping, "observed"),
          residual_kind(required<std::string>(target_mapping, "residual")),
          required<double>(target_mapping, "scale"),
          required<double>(target_mapping, "value_factor"),
      };
      if (!(std::isfinite(target.observed) && std::isfinite(target.scale) &&
            target.scale > 0.0 && std::isfinite(target.value_factor) &&
            target.value_factor > 0.0)) {
        throw std::invalid_argument("target values and scales must be finite");
      }
      if ((kind == ObservationKind::FixedStatePressure ||
           target.residual == ResidualKind::LogRatio) &&
          target.observed <= 0.0) {
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

template <typename Block>
EosEvaluation fixed_state_pressure_evaluation(const Block &block,
                                              double temperature_k,
                                              std::size_t width,
                                              std::size_t volume_coordinate,
                                              std::size_t parameter_coordinate,
                                              std::size_t parameter_count) {
  if (block.hessian.size() != width * width) {
    throw std::runtime_error("EOS fixed-state Hessian has the wrong shape");
  }
  EosEvaluation result{{block.pressure_pa}, {}};
  result.derivatives.reserve(parameter_count);
  for (std::size_t parameter = 0; parameter < parameter_count; ++parameter) {
    result.derivatives.push_back(
        -epcsaft::detail::GAS_CONSTANT * temperature_k *
        block.hessian[volume_coordinate * width + parameter_coordinate +
                      parameter]);
  }
  return result;
}

class EosCost final : public ceres::CostFunction {
public:
  EosCost(const NativeModel &model,
          std::vector<NativeActiveCoordinate> active_coordinates,
          std::vector<std::size_t> identity_coordinates,
          std::vector<double> coordinate_scales, Observation observation)
      : model_(model), active_coordinates_(std::move(active_coordinates)),
        identity_coordinates_(std::move(identity_coordinates)),
        coordinate_scales_(std::move(coordinate_scales)),
        observation_(std::move(observation)) {
    set_num_residuals(static_cast<int>(observation_.targets.size()));
    mutable_parameter_block_sizes()->assign(coordinate_scales_.size(), 1);
  }

  bool Evaluate(double const *const *parameters, double *residuals,
                double **jacobians) const override {
    try {
      auto active = active_coordinates_;
      for (std::size_t index = 0; index < active.size(); ++index) {
        const std::size_t coordinate = identity_coordinates_[index];
        active[index].value =
            parameters[coordinate][0] * coordinate_scales_[coordinate];
      }
      const EosEvaluation evaluation = evaluate_eos(active);
      for (std::size_t row = 0; row < observation_.targets.size(); ++row) {
        const auto &target = observation_.targets[row];
        const double predicted =
            evaluation.values[target.value_index] * target.value_factor;
        if (!std::isfinite(predicted) ||
            ((observation_.kind == ObservationKind::FixedStatePressure ||
              target.residual == ResidualKind::LogRatio) &&
             predicted <= 0.0)) {
          return reject("EOS returned an inadmissible target value");
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
          for (std::size_t identity = 0; identity < active_coordinates_.size();
               ++identity) {
            if (identity_coordinates_[identity] == coordinate) {
              derivative +=
                  evaluation
                      .derivatives[target.value_index * active_coordinates_.size() +
                                   identity];
            }
          }
          jacobians[coordinate][row] = factor * target.value_factor *
                                       derivative *
                                       coordinate_scales_[coordinate];
        }
      }
      if (!std::all_of(residuals, residuals + observation_.targets.size(),
                       [](double value) { return std::isfinite(value); })) {
        return reject("residual evaluation was nonfinite");
      }
      return true;
    } catch (const std::exception &error) {
      return reject(error.what());
    }
  }

  std::size_t rejection_count() const noexcept { return rejection_count_; }
  const std::string &first_failure() const noexcept { return first_failure_; }

private:
  EosEvaluation
  evaluate_eos(const std::vector<NativeActiveCoordinate> &active) const {
    if (observation_.kind == ObservationKind::PureSaturation) {
      SaturationSensitivityResult result;
      {
        std::lock_guard<std::mutex> lock(saturation_solver_mutex());
        result = solve_local_saturation_active(
            model_, observation_.temperature_k,
            observation_.liquid_density_upper_mol_m3, active);
      }
      if (!result.solve.accepted) {
        throw std::runtime_error(result.solve.failure_reason);
      }
      return {
          {
              result.solve.evaluation.saturation_pressure_pa,
              1.0 / result.solve.evaluation.vapor.volume_m3,
              1.0 / result.solve.evaluation.liquid.volume_m3,
          },
          std::move(result.derivatives),
      };
    }

    const std::size_t components = model_.component_count();
    const std::size_t parameters = active.size();
    const double volume = 1.0 / observation_.molar_density_mol_m3;
    const std::size_t width = components + 1 + parameters;
    const auto amounts = components == 1 ? std::vector<double>{1.0}
                                         : observation_.composition;
    const auto block = model_.evaluate_phase_derivatives(
        observation_.temperature_k, amounts, volume, active, 2);
    return fixed_state_pressure_evaluation(block, observation_.temperature_k,
                                           width, components, components + 1,
                                           parameters);
  }

  bool reject(const std::string &reason) const {
    ++rejection_count_;
    if (first_failure_.empty()) {
      first_failure_ =
          reason.empty() ? "EOS evaluation was non-evaluable" : reason;
    }
    return false;
  }

  NativeModel model_;
  std::vector<NativeActiveCoordinate> active_coordinates_;
  std::vector<std::size_t> identity_coordinates_;
  std::vector<double> coordinate_scales_;
  Observation observation_;
  mutable std::size_t rejection_count_{0};
  mutable std::string first_failure_;
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

StartResult solve_start(std::size_t index, const NativeModel &model,
                        const std::vector<NativeActiveCoordinate> &active_coordinates,
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
  std::vector<std::unique_ptr<EosCost>> costs;
  costs.reserve(observations.size());
  ceres::Problem::Options problem_options;
  problem_options.cost_function_ownership = ceres::DO_NOT_TAKE_OWNERSHIP;
  ceres::Problem problem(problem_options);
  std::vector<double *> blocks;
  blocks.reserve(result.final_affine.size());
  for (double &value : result.final_affine) {
    blocks.push_back(&value);
  }
  for (const auto &observation : observations) {
    costs.push_back(std::make_unique<EosCost>(
        model, active_coordinates, identity_coordinates, scales, observation));
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
  for (const auto &cost : costs) {
    result.non_evaluable_count += cost->rejection_count();
    if (result.failure_reason.empty() && !cost->first_failure().empty()) {
      result.failure_reason = cost->first_failure();
    }
  }

  std::size_t residual_count = 0;
  for (const auto &observation : observations) {
    residual_count += observation.targets.size();
  }
  result.residuals.resize(residual_count);
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
    if (!cost->Evaluate(final_blocks.data(), result.residuals.data() + offset,
                        jacobian_pointers.data())) {
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
  result["residuals"] = value.residuals;
  result["physical_jacobian"] = physical_jacobian;
  result["affine_jacobian"] = value.affine_jacobian;
  result["diagnostics"] = diagnostics_mapping(value.diagnostics);
  result["active_bounds"] = value.active_bounds;
  result["solver_status"] = value.solver_status;
  result["failure_reason"] = value.failure_reason;
  result["iterations"] = value.iterations;
  result["non_evaluable_count"] = value.non_evaluable_count;
  result["cost"] = value.cost;
  result["usable"] = value.usable;
  if (value.usable) {
    result["basin"] = value.basin;
  } else {
    result["basin"] = py::none();
  }
  return result;
}

py::dict fit_resolved(const NativeModel &model,
                      const std::vector<NativeActiveCoordinate> &active_coordinates,
                      const std::vector<std::size_t> &identity_coordinates,
                      const py::iterable &coordinate_values,
                      const py::iterable &observation_values,
                      const py::dict &control_values) {
  const auto coordinates = parse_coordinates(coordinate_values);
  const auto observations = parse_observations(observation_values);
  const auto controls = parse_controls(control_values);
  if (active_coordinates.empty() ||
      active_coordinates.size() != identity_coordinates.size() ||
      std::any_of(identity_coordinates.begin(), identity_coordinates.end(),
                  [&](std::size_t coordinate) {
                    return coordinate >= coordinates.size();
                  })) {
    throw std::invalid_argument("identity coordinates must bind every active "
                                "identity to one coordinate");
  }
  if (std::any_of(observations.begin(), observations.end(),
                  [&](const Observation &observation) {
                    return observation.kind ==
                               ObservationKind::FixedStatePressure &&
                           observation.composition.size() !=
                               model.component_count();
                  })) {
    throw std::invalid_argument(
        "fixed-state composition must match the EOS component count");
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
              index, model, active_coordinates, identity_coordinates, coordinates,
              observations, controls, start_values[index]);
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
  output["single_threaded_solves"] = true;
  output["start_order"] = "halton";
  return output;
}

} // namespace
} // namespace epcsaft::regression

void bind_parameter_fit(py::module_ &module) {
  module.def("_fit_resolved", &epcsaft::regression::fit_resolved,
             py::arg("model"), py::arg("active_coordinates"),
             py::arg("identity_coordinates"), py::arg("coordinates"),
             py::arg("observations"), py::arg("controls"));
}
