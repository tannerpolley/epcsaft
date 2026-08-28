#include "fit_core.hpp"
#include "pure_saturation.hpp"

#include <Eigen/Dense>
#include <Eigen/SVD>
#include <coin/IpIpoptApplication.hpp>
#include <coin/IpSolveStatistics.hpp>
#include <coin/IpTNLP.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

namespace epcsaft::regression {
namespace {

constexpr double kOptimizerBasinAffineDistanceMax = 1.0e-6;
constexpr double kOptimizerBasinCostDifferenceFactor = 1.0e-8;

struct OwnerEvaluation final {
  std::vector<double> values;
  std::vector<double> derivatives;
  std::optional<epcsaft_engine::equilibrium::GrepeResult> owner_result;
};

class OwnerEvaluationError final : public std::runtime_error {
public:
  OwnerEvaluationError(std::string code, std::string reason)
      : std::runtime_error(std::move(reason)), code_(std::move(code)) {}

  const std::string &code() const noexcept { return code_; }

private:
  std::string code_;
};

NativeModel with_active_values(const NativeModel &base,
                               const std::vector<std::string> &identities,
                               const std::vector<double> &active) {
  auto coordinates = base.active_coordinates();
  for (auto &coordinate : coordinates) {
    const auto found =
        std::find(identities.begin(), identities.end(), coordinate.identity);
    if (found == identities.end()) {
      throw OwnerEvaluationError("parameter_topology_mismatch",
                                 "native active coordinate identity changed");
    }
    const std::size_t index =
        static_cast<std::size_t>(found - identities.begin());
    coordinate.value = active[index];
  }
  return base.with_active_coordinates(coordinates);
}

OwnerEvaluation evaluate_direct(const DirectEosExperiment &experiment,
                                const std::vector<double> &active) {
  std::vector<std::string> identities;
  for (const auto &coordinate : experiment.model.active_coordinates()) {
    identities.push_back(coordinate.identity);
  }
  const NativeModel model =
      with_active_values(experiment.model, identities, active);
  std::vector<double> values;
  std::vector<double> derivatives;
  const std::size_t parameter_count = active.size();
  if (experiment.kind == DirectEosKind::FixedState) {
    if (!(experiment.molar_density_mol_m3 > 0.0) ||
        experiment.composition.size() != model.component_count()) {
      throw OwnerEvaluationError("parameter_topology_mismatch",
                                 "fixed-state native experiment is invalid");
    }
    const double volume = 1.0 / experiment.molar_density_mol_m3;
    const auto block = model.evaluate_phase_derivatives(
        experiment.temperature_k, experiment.composition, volume,
        model.active_coordinates(), 2);
    const std::size_t component_count = experiment.composition.size();
    const std::size_t width = component_count + 1 + parameter_count;
    values.push_back(block.helmholtz_over_rt_reference_amount);
    values.push_back(block.pressure_pa);
    values.insert(values.end(), block.gradient.begin(),
                  block.gradient.begin() +
                      static_cast<std::ptrdiff_t>(component_count));
    derivatives.resize(values.size() * parameter_count);
    for (std::size_t parameter = 0; parameter < parameter_count; ++parameter) {
      const std::size_t column = component_count + 1 + parameter;
      derivatives[parameter] = block.gradient[column];
      derivatives[parameter_count + parameter] =
          -8.31446261815324 * experiment.temperature_k *
          block.hessian[component_count * width + column];
      for (std::size_t component = 0; component < component_count; ++component) {
        derivatives[(component + 2) * parameter_count + parameter] =
            block.hessian[component * width + column];
      }
    }
  } else {
    const auto saturation = solve_local_saturation_active(
        model, experiment.temperature_k, model.active_coordinates());
    if (!saturation.solve.accepted) {
      throw OwnerEvaluationError(
          "physical_rejection",
          saturation.solve.rejection_status + ": " +
              (saturation.solve.failure_reason.empty()
                   ? saturation.solve.solver_status
                   : saturation.solve.failure_reason));
    }
    const auto &solve = saturation.solve;
    values = {solve.evaluation.saturation_pressure_pa,
              solve.evaluation.vapor.amount_mol / solve.evaluation.vapor.volume_m3,
              solve.evaluation.liquid.amount_mol / solve.evaluation.liquid.volume_m3};
    derivatives = saturation.derivatives;
    const auto &molar_masses = model.molar_masses_kg_per_mol();
    if (molar_masses.size() == 1 && std::isfinite(molar_masses.front())) {
      values.push_back(values[1] * molar_masses.front());
      values.push_back(values[2] * molar_masses.front());
      for (std::size_t parameter = 0; parameter < parameter_count; ++parameter) {
        derivatives.push_back(
            derivatives[parameter_count + parameter] * molar_masses.front());
      }
      for (std::size_t parameter = 0; parameter < parameter_count; ++parameter) {
        derivatives.push_back(
            derivatives[2 * parameter_count + parameter] * molar_masses.front());
      }
    }
  }
  OwnerEvaluation selected;
  selected.values.reserve(experiment.selected_rows.size());
  selected.derivatives.reserve(experiment.selected_rows.size() * parameter_count);
  for (const std::size_t row : experiment.selected_rows) {
    if (row >= values.size() || (row + 1) * parameter_count > derivatives.size()) {
      throw OwnerEvaluationError("owner_contract_failure",
                                 "direct EOS selected row is unavailable");
    }
    selected.values.push_back(values[row]);
    selected.derivatives.insert(
        selected.derivatives.end(),
        derivatives.begin() + static_cast<std::ptrdiff_t>(row * parameter_count),
        derivatives.begin() +
            static_cast<std::ptrdiff_t>((row + 1) * parameter_count));
  }
  return selected;
}

OwnerEvaluation select_grepe(const GrepeExperiment &experiment,
                             std::vector<double> active,
                             epcsaft_engine::equilibrium::GrepeResult native,
                             bool retain_result) {
  const auto &equilibrium = native.equilibrium;
  if (!equilibrium.accepted) {
    std::string code = native.failure_code;
    if (code.empty()) {
      code = equilibrium.provider_domain_status == "failed"
                 ? "provider_domain_rejection"
                 : equilibrium.physical_status == "failed"
                       ? "physical_rejection"
                       : "numerical_convergence_failure";
    }
    throw OwnerEvaluationError(
        code, equilibrium.failure_reason.empty() ? equilibrium.solver_status
                                                 : equilibrium.failure_reason);
  }
  const std::size_t input_width = equilibrium.input_order.size();
  if (input_width == 0 ||
      equilibrium.output_input_derivatives.size() !=
          equilibrium.output_values.size() * input_width) {
    throw OwnerEvaluationError("owner_contract_failure",
                               "GREPE returned an invalid ordered Jacobian");
  }
  std::vector<std::size_t> active_columns;
  active_columns.reserve(experiment.model.active_coordinates().size());
  for (const auto &coordinate : experiment.model.active_coordinates()) {
    const auto found = std::find(equilibrium.input_order.begin(),
                                 equilibrium.input_order.end(),
                                 coordinate.identity);
    if (found == equilibrium.input_order.end()) {
      throw OwnerEvaluationError("first_order_failure",
                                 "GREPE omitted an active derivative column");
    }
    active_columns.push_back(static_cast<std::size_t>(
        found - equilibrium.input_order.begin()));
  }
  OwnerEvaluation selected;
  selected.values.reserve(experiment.selected_rows.size());
  selected.derivatives.reserve(experiment.selected_rows.size() * active.size());
  for (const std::size_t row : experiment.selected_rows) {
    if (row >= equilibrium.output_values.size()) {
      throw OwnerEvaluationError("owner_contract_failure",
                                 "GREPE selected output row is unavailable");
    }
    selected.values.push_back(equilibrium.output_values[row]);
    for (const std::size_t column : active_columns) {
      selected.derivatives.push_back(
          equilibrium.output_input_derivatives[row * input_width + column]);
    }
  }
  if (selected.derivatives.size() != selected.values.size() * active.size()) {
    throw OwnerEvaluationError("owner_contract_failure",
                               "GREPE active Jacobian order changed");
  }
  if (retain_result) {
    selected.owner_result = std::move(native);
  }
  return selected;
}

OwnerEvaluation evaluate_grepe(const GrepeExperiment &experiment,
                               const std::vector<double> &active,
                               bool retain_result) {
  std::vector<std::string> identities;
  for (const auto &coordinate : experiment.model.active_coordinates()) {
    identities.push_back(coordinate.identity);
  }
  auto native = epcsaft_engine::equilibrium::solve_grepe(
      with_active_values(experiment.model, identities, active),
      with_active_values(experiment.incipient_model, identities, active),
      experiment.request);
  return select_grepe(experiment, active, std::move(native), retain_result);
}

OwnerEvaluation finalize_grepe(
    const GrepeExperiment &experiment, const std::vector<double> &active,
    const epcsaft_engine::equilibrium::GrepeLocalStructure &structure,
    const epcsaft_engine::equilibrium::GrepeLocalPoint &point,
    const std::vector<double> &local, const std::vector<double> &lower,
    const std::vector<double> &upper) {
  auto native = epcsaft_engine::equilibrium::finalize_grepe_local(
      with_active_values(experiment.model,
                         structure.active_parameter_identities, active),
      with_active_values(experiment.incipient_model,
                         structure.active_parameter_identities, active),
      experiment.request, structure, point, local, lower, upper);
  return select_grepe(experiment, active, std::move(native), true);
}

class NativeCost final {
public:
  NativeCost(const Experiment &experiment, std::size_t active_parameter_count,
            std::vector<std::size_t> identity_coordinates,
            std::vector<double> coordinate_origins,
            std::vector<double> coordinate_scales, Observation observation,
            std::vector<std::string> failure_codes)
      : experiment_(experiment), active_parameter_count_(active_parameter_count),
        identity_coordinates_(std::move(identity_coordinates)),
        coordinate_origins_(std::move(coordinate_origins)),
        coordinate_scales_(std::move(coordinate_scales)),
        observation_(std::move(observation)),
        failure_codes_(std::move(failure_codes)) {}

  std::size_t row_count() const noexcept { return observation_.targets.size(); }

  bool evaluate(double const *const *parameters, double *residuals,
                double **jacobians, double *predictions, bool retain_result,
                std::optional<epcsaft_engine::equilibrium::GrepeResult>
                    *owner_result) const {
    try {
      std::vector<double> active(active_parameter_count_);
      for (std::size_t index = 0; index < active_parameter_count_; ++index) {
        const std::size_t coordinate = identity_coordinates_[index];
        active[index] =
            coordinate_origins_[coordinate] +
            parameters[coordinate][0] * coordinate_scales_[coordinate];
      }
      const OwnerEvaluation evaluation = std::visit(
          [&](const auto &experiment) -> OwnerEvaluation {
            using Value = std::decay_t<decltype(experiment)>;
            if constexpr (std::is_same_v<Value, DirectEosExperiment>) {
              return evaluate_direct(experiment, active);
            } else {
              return evaluate_grepe(experiment, active, retain_result);
            }
          },
          experiment_);
      if (owner_result != nullptr) {
        *owner_result = evaluation.owner_result;
      }
      return write_evaluation(evaluation, residuals, jacobians, predictions);
    } catch (const OwnerEvaluationError &error) {
      return reject(error.code(), error.what());
    } catch (const std::exception &error) {
      return reject("owner_contract_failure", error.what());
    }
  }

  bool finalize(double const *const *parameters,
                const epcsaft_engine::equilibrium::GrepeLocalStructure &structure,
                const epcsaft_engine::equilibrium::GrepeLocalPoint &point,
                const std::vector<double> &local,
                const std::vector<double> &lower,
                const std::vector<double> &upper, double *residuals,
                double **jacobians, double *predictions,
                std::optional<epcsaft_engine::equilibrium::GrepeResult>
                    *owner_result) const {
    try {
      const auto *experiment = std::get_if<GrepeExperiment>(&experiment_);
      if (experiment == nullptr) {
        return reject("owner_contract_failure",
                      "full-space finalization requires GREPE");
      }
      std::vector<double> active(active_parameter_count_);
      for (std::size_t index = 0; index < active_parameter_count_; ++index) {
        const std::size_t coordinate = identity_coordinates_[index];
        active[index] = coordinate_origins_[coordinate] +
                        parameters[coordinate][0] * coordinate_scales_[coordinate];
      }
      const OwnerEvaluation evaluation =
          finalize_grepe(*experiment, active, structure, point, local, lower,
                         upper);
      if (owner_result != nullptr) {
        *owner_result = evaluation.owner_result;
      }
      return write_evaluation(evaluation, residuals, jacobians, predictions);
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
  bool write_evaluation(const OwnerEvaluation &evaluation, double *residuals,
                        double **jacobians, double *predictions) const {
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
        factor = multiplier / (predicted * target.scale);
        residuals[row] =
            multiplier * std::log(predicted / target.observed) / target.scale;
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
            derivative += evaluation.derivatives[
                row * active_parameter_count_ + identity];
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
  }

  bool reject(const std::string &raw_code, const std::string &reason) const {
    const std::string code =
        std::find(failure_codes_.begin(), failure_codes_.end(), raw_code) ==
                failure_codes_.end()
            ? "owner_contract_failure"
            : raw_code;
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

  const Experiment &experiment_;
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

struct RetainedLocalState final {
  std::vector<double> values;
  std::vector<double> lower;
  std::vector<double> upper;
};

class ReducedSpaceNlp final : public Ipopt::TNLP {
public:
  ReducedSpaceNlp(const std::vector<std::unique_ptr<NativeCost>> &costs,
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
  std::vector<RetainedLocalState> retained_local_states() const { return {}; }
  std::vector<epcsaft_engine::equilibrium::GrepeLocalPoint>
  retained_points() const { return {}; }

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

  const std::vector<std::unique_ptr<NativeCost>> &costs_;
  const std::vector<Coordinate> &coordinates_;
  std::vector<double> initial_;
  std::vector<double> solution_;
  mutable std::vector<double> cached_x_;
  mutable std::vector<double> cached_gradient_;
  mutable double cached_objective_{0.0};
  mutable bool cached_evaluable_{false};
};

struct FullSpaceBlock final {
  const GrepeExperiment *experiment;
  const Observation *observation;
  std::size_t structure_index;
  std::size_t local_offset;
  std::size_t local_width;
  std::size_t constraint_offset;
  std::vector<double> lower;
  std::vector<double> upper;
  std::vector<double> constraints;
  std::vector<double> constraint_local_jacobian;
  std::vector<double> constraint_parameter_jacobian;
  std::vector<double> outputs;
  std::vector<double> output_local_jacobian;
  std::vector<double> output_parameter_jacobian;
};

struct FullSpaceStructure final {
  std::string fingerprint;
  std::size_t experiment_index;
  epcsaft_engine::equilibrium::GrepeLocalStructure value;
};

class FullSpaceCampaign final {
public:
  explicit FullSpaceCampaign(const std::vector<Experiment> &experiments) {
    experiment_structures_.reserve(experiments.size());
    for (std::size_t index = 0; index < experiments.size(); ++index) {
      const auto &experiment = std::get<GrepeExperiment>(experiments[index]);
      const auto retained = std::find_if(
          structures_.begin(), structures_.end(), [&](const auto &value) {
            return value.fingerprint == experiment.structure_fingerprint;
          });
      if (retained != structures_.end()) {
        experiment_structures_.push_back(
            static_cast<std::size_t>(retained - structures_.begin()));
        continue;
      }
      experiment_structures_.push_back(structures_.size());
      structures_.push_back({
          experiment.structure_fingerprint,
          index,
          epcsaft_engine::equilibrium::compile_grepe_local_structure(
              experiment.model, experiment.incipient_model, experiment.request),
      });
    }
  }

  std::size_t structure_index(std::size_t experiment) const {
    return experiment_structures_.at(experiment);
  }

  const epcsaft_engine::equilibrium::GrepeLocalStructure &
  structure(std::size_t index) const {
    return structures_.at(index).value;
  }

  std::vector<epcsaft_engine::equilibrium::GrepeLocalPoint>
  evaluate_points(const std::vector<Experiment> &experiments,
                  const std::vector<double> &active) const {
    std::vector<epcsaft_engine::equilibrium::GrepeLocalPoint> result;
    result.reserve(structures_.size());
    for (const auto &retained : structures_) {
      const auto &experiment =
          std::get<GrepeExperiment>(experiments[retained.experiment_index]);
      result.push_back(epcsaft_engine::equilibrium::evaluate_grepe_local_point(
          with_active_values(experiment.model,
                             retained.value.active_parameter_identities, active),
          with_active_values(experiment.incipient_model,
                             retained.value.active_parameter_identities, active),
          experiment.request, retained.value));
    }
    return result;
  }

private:
  std::vector<FullSpaceStructure> structures_;
  std::vector<std::size_t> experiment_structures_;
};

class SparseFullSpaceNlp final : public Ipopt::TNLP {
public:
  SparseFullSpaceNlp(const std::vector<Experiment> &experiments,
                     const std::vector<Observation> &observations,
                     const FullSpaceCampaign &campaign,
                     const std::vector<std::size_t> &identity_coordinates,
                     const std::vector<Coordinate> &coordinates,
                     std::vector<double> affine)
      : experiments_(experiments), campaign_(campaign),
        identity_coordinates_(identity_coordinates), coordinates_(coordinates),
        initial_(std::move(affine)) {
    lower_.reserve(initial_.size());
    upper_.reserve(initial_.size());
    for (const auto &coordinate : coordinates_) {
      lower_.push_back(coordinate.lower_affine);
      upper_.push_back(coordinate.upper_affine);
    }
    const std::vector<double> active = active_values(initial_.data());
    points_ = campaign_.evaluate_points(experiments_, active);
    std::size_t constraint_offset = 0;
    for (std::size_t index = 0; index < experiments.size(); ++index) {
      const auto &experiment = std::get<GrepeExperiment>(experiments[index]);
      const std::size_t structure_index = campaign_.structure_index(index);
      const auto &structure = campaign_.structure(structure_index);
      const auto native = epcsaft_engine::equilibrium::initialize_grepe_local(
          with_active_values(experiment.model,
                             structure.active_parameter_identities, active),
          with_active_values(experiment.incipient_model,
                             structure.active_parameter_identities, active),
          experiment.request, structure, points_[structure_index]);
      const auto &equilibrium = native.equilibrium;
      if (equilibrium.local_variables.empty() ||
          equilibrium.local_lower_bounds.size() !=
              equilibrium.local_variables.size() ||
          equilibrium.local_upper_bounds.size() !=
              equilibrium.local_variables.size()) {
        throw OwnerEvaluationError(
            native.failure_code.empty() ? "physical_rejection"
                                        : native.failure_code,
            equilibrium.failure_reason.empty() ? equilibrium.solver_status
                                               : equilibrium.failure_reason);
      }
      const std::size_t local_offset = initial_.size();
      const std::size_t local_width = equilibrium.local_variables.size();
      const std::size_t output_count = experiment.request.outputs.size();
      FullSpaceBlock block{
          &experiment,
          &observations[index],
          structure_index,
          local_offset,
          local_width,
          constraint_offset,
          equilibrium.local_lower_bounds,
          equilibrium.local_upper_bounds,
          std::vector<double>(local_width),
          std::vector<double>(local_width * local_width),
          std::vector<double>(
              local_width * structure.active_parameter_identities.size()),
          std::vector<double>(output_count),
          std::vector<double>(output_count * local_width),
          std::vector<double>(output_count
                              * structure.active_parameter_identities.size()),
      };
      blocks_.push_back(std::move(block));
      initial_.insert(initial_.end(), equilibrium.local_variables.begin(),
                      equilibrium.local_variables.end());
      lower_.insert(lower_.end(), equilibrium.local_lower_bounds.begin(),
                    equilibrium.local_lower_bounds.end());
      upper_.insert(upper_.end(), equilibrium.local_upper_bounds.begin(),
                    equilibrium.local_upper_bounds.end());
      for (std::size_t row = 0; row < local_width; ++row) {
        for (std::size_t coordinate = 0; coordinate < coordinates_.size();
             ++coordinate) {
          jacobian_rows_.push_back(
              static_cast<int>(constraint_offset + row));
          jacobian_columns_.push_back(static_cast<int>(coordinate));
        }
        for (std::size_t local = 0; local < local_width; ++local) {
          jacobian_rows_.push_back(
              static_cast<int>(constraint_offset + row));
          jacobian_columns_.push_back(
              static_cast<int>(local_offset + local));
        }
      }
      constraint_offset += local_width;
      residual_count_ += observations[index].targets.size();
    }
    constraint_count_ = constraint_offset;
    solution_ = initial_;
    source_theta_ = active;
    constraints_.reserve(constraint_count_);
    constraint_jacobian_.reserve(jacobian_rows_.size());
  }

  bool get_nlp_info(Ipopt::Index &n, Ipopt::Index &m,
                    Ipopt::Index &nnz_jac_g, Ipopt::Index &nnz_h_lag,
                    IndexStyleEnum &index_style) override {
    n = static_cast<Ipopt::Index>(initial_.size());
    m = static_cast<Ipopt::Index>(constraint_count_);
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
    std::fill(g_l, g_l + m, 0.0);
    std::fill(g_u, g_u + m, 0.0);
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
  std::vector<RetainedLocalState> retained_local_states() const {
    std::vector<RetainedLocalState> states;
    states.reserve(blocks_.size());
    for (const auto &block : blocks_) {
      states.push_back({
          std::vector<double>(
              solution_.begin() + static_cast<std::ptrdiff_t>(block.local_offset),
              solution_.begin() + static_cast<std::ptrdiff_t>(
                                      block.local_offset + block.local_width)),
          block.lower,
          block.upper,
      });
    }
    return states;
  }
  std::vector<epcsaft_engine::equilibrium::GrepeLocalPoint>
  retained_points() const { return points_; }
private:
  bool dimensions(Ipopt::Index n, Ipopt::Index m) const noexcept {
    return n == static_cast<Ipopt::Index>(initial_.size()) &&
           m == static_cast<Ipopt::Index>(constraint_count_);
  }

  std::vector<double> active_values(const double *x) const {
    std::vector<double> result(identity_coordinates_.size());
    for (std::size_t identity = 0; identity < result.size(); ++identity) {
      const std::size_t coordinate = identity_coordinates_[identity];
      result[identity] = coordinates_[coordinate].origin +
                         x[coordinate] * coordinates_[coordinate].scale;
    }
    return result;
  }

  bool evaluate(Ipopt::Index n, const Ipopt::Number *x) const {
    if (n != static_cast<Ipopt::Index>(initial_.size()) || x == nullptr) {
      return false;
    }
    if (cached_x_.size() == initial_.size() &&
        std::equal(x, x + n, cached_x_.begin())) {
      return evaluable_;
    }
    cached_x_.assign(x, x + n);
    try {
      objective_ = 0.0;
      gradient_.assign(initial_.size(), 0.0);
      constraints_.clear();
      constraint_jacobian_.clear();
      const std::vector<double> active = active_values(x);
      if (source_theta_ != active) {
        points_ = campaign_.evaluate_points(experiments_, active);
        source_theta_ = active;
      }
      for (auto &block : blocks_) {
        const auto &structure = campaign_.structure(block.structure_index);
        const std::vector<double> local(
            x + static_cast<std::ptrdiff_t>(block.local_offset),
            x + static_cast<std::ptrdiff_t>(block.local_offset +
                                            block.local_width));
        epcsaft_engine::equilibrium::evaluate_grepe_local_blocks(
            with_active_values(block.experiment->model,
                               structure.active_parameter_identities, active),
            with_active_values(block.experiment->incipient_model,
                               structure.active_parameter_identities, active),
            block.experiment->request, structure, points_[block.structure_index],
            local, block.lower, block.upper,
            {
                block.constraints.data(),
                block.constraints.size(),
                block.constraint_local_jacobian.data(),
                block.constraint_local_jacobian.size(),
                block.constraint_parameter_jacobian.data(),
                block.constraint_parameter_jacobian.size(),
                block.outputs.data(),
                block.outputs.size(),
                block.output_local_jacobian.data(),
                block.output_local_jacobian.size(),
                block.output_parameter_jacobian.data(),
                block.output_parameter_jacobian.size(),
            });
        for (std::size_t target_index = 0;
             target_index < block.observation->targets.size(); ++target_index) {
          const std::size_t output_row =
              block.experiment->selected_rows[target_index];
          const auto &target = block.observation->targets[target_index];
          const double predicted = block.outputs[output_row];
          const double multiplier =
              block.observation->multiplier * target.multiplier;
          double factor = multiplier / target.scale;
          double residual;
          switch (target.residual) {
          case ResidualKind::LogRatio:
            factor = multiplier / (predicted * target.scale);
            residual = multiplier * std::log(predicted / target.observed) /
                       target.scale;
            break;
          case ResidualKind::LeftCensored:
            if (predicted <= target.observed) {
              factor = residual = 0.0;
            } else {
              residual = multiplier * (predicted - target.observed) /
                         target.scale;
            }
            break;
          case ResidualKind::RightCensored:
            if (predicted >= target.observed) {
              factor = residual = 0.0;
            } else {
              residual = multiplier * (predicted - target.observed) /
                         target.scale;
            }
            break;
          default:
            residual = multiplier * (predicted - target.observed) /
                       target.scale;
            break;
          }
          objective_ += 0.5 * residual * residual;
          for (std::size_t identity = 0;
               identity < identity_coordinates_.size(); ++identity) {
            const std::size_t parameter_width =
                block.output_parameter_jacobian.size() / block.outputs.size();
            const std::size_t coordinate = identity_coordinates_[identity];
            gradient_[coordinate] +=
                residual * factor *
                block.output_parameter_jacobian[
                    output_row * parameter_width + identity] *
                coordinates_[coordinate].scale;
          }
          for (std::size_t local_column = 0;
               local_column < block.local_width; ++local_column) {
            gradient_[block.local_offset + local_column] +=
                residual * factor *
                block.output_local_jacobian[
                    output_row * block.local_width + local_column];
          }
        }
        constraints_.insert(constraints_.end(),
                            block.constraints.begin(), block.constraints.end());
        const std::size_t input_width =
            block.constraint_parameter_jacobian.size() / block.local_width;
        for (std::size_t row = 0; row < block.local_width; ++row) {
          for (std::size_t coordinate = 0; coordinate < coordinates_.size();
               ++coordinate) {
            double derivative = 0.0;
            for (std::size_t identity = 0;
                 identity < identity_coordinates_.size(); ++identity) {
              if (identity_coordinates_[identity] != coordinate) {
                continue;
              }
              derivative +=
                  block.constraint_parameter_jacobian[row * input_width + identity];
            }
            constraint_jacobian_.push_back(
                derivative * coordinates_[coordinate].scale);
          }
          constraint_jacobian_.insert(
              constraint_jacobian_.end(),
              block.constraint_local_jacobian.begin() +
                  static_cast<std::ptrdiff_t>(row * block.local_width),
              block.constraint_local_jacobian.begin() +
                  static_cast<std::ptrdiff_t>((row + 1) * block.local_width));
        }
      }
      const auto finite = [](const std::vector<double> &values) {
        return std::all_of(values.begin(), values.end(),
                           [](double value) { return std::isfinite(value); });
      };
      if (constraints_.size() != constraint_count_ ||
          constraint_jacobian_.size() != jacobian_rows_.size() ||
          !std::isfinite(objective_) || !finite(gradient_) ||
          !finite(constraints_) || !finite(constraint_jacobian_)) {
        throw std::domain_error("native full-space values are invalid");
      }
      evaluable_ = true;
      return true;
    } catch (const std::exception &error) {
      evaluable_ = false;
      return false;
    }
  }

  const std::vector<Experiment> &experiments_;
  const FullSpaceCampaign &campaign_;
  const std::vector<std::size_t> &identity_coordinates_;
  const std::vector<Coordinate> &coordinates_;
  mutable std::vector<FullSpaceBlock> blocks_;
  std::vector<double> initial_;
  std::vector<double> lower_;
  std::vector<double> upper_;
  std::size_t residual_count_{0};
  std::size_t constraint_count_{0};
  std::vector<int> jacobian_rows_;
  std::vector<int> jacobian_columns_;
  std::vector<double> solution_;
  mutable std::vector<double> cached_x_;
  mutable std::vector<double> source_theta_;
  mutable std::vector<epcsaft_engine::equilibrium::GrepeLocalPoint> points_;
  mutable std::vector<double> gradient_;
  mutable std::vector<double> constraints_;
  mutable std::vector<double> constraint_jacobian_;
  mutable double objective_{0.0};
  mutable bool evaluable_{false};
};

template <typename Problem>
Ipopt::ApplicationReturnStatus solve_ipopt(Problem *problem,
                                           const Controls &controls,
                                           int &iterations,
                                           std::vector<double> &solution,
                                           bool full_space = false,
                                           std::vector<RetainedLocalState>
                                               *retained_states = nullptr,
                                           std::vector<epcsaft_engine::equilibrium::GrepeLocalPoint>
                                               *retained_points = nullptr) {
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
  if (retained_states != nullptr) {
    *retained_states = problem->retained_local_states();
  }
  if (retained_points != nullptr) {
    *retained_points = problem->retained_points();
  }
  return status;
}

StartResult solve_start(std::size_t index, const FitRequest &request,
                        bool full_space,
                        const FullSpaceCampaign *full_space_campaign,
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
  std::vector<std::unique_ptr<NativeCost>> costs;
  costs.reserve(observations.size());
  for (std::size_t observation_index = 0;
       observation_index < observations.size(); ++observation_index) {
    costs.push_back(std::make_unique<NativeCost>(
        request.experiments[observation_index], request.active_parameter_count,
        identity_coordinates, origins, scales, observations[observation_index],
        failure_codes));
  }
  Ipopt::ApplicationReturnStatus solver_status;
  std::vector<RetainedLocalState> retained_states;
  std::vector<epcsaft_engine::equilibrium::GrepeLocalPoint> retained_points;
  if (!full_space) {
    auto *problem = new ReducedSpaceNlp(costs, coordinates, result.final_affine);
    solver_status =
        solve_ipopt(problem, controls, result.iterations, result.final_affine);
  } else {
    if (full_space_campaign == nullptr) {
      throw std::logic_error("full-space GREPE campaign is unavailable");
    }
    std::vector<double> solution;
    auto *problem = new SparseFullSpaceNlp(
        request.experiments, observations, *full_space_campaign,
        identity_coordinates, coordinates, result.final_affine);
    solver_status = solve_ipopt(problem, controls, result.iterations, solution,
                                true, &retained_states, &retained_points);
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
  const bool retain_owner_results = solver_usable;
  result.owner_results.resize(costs.size());
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
    const bool evaluated = full_space
        ? solver_usable && batch_index < retained_states.size()
            && cost->finalize(
                final_blocks.data(),
                full_space_campaign->structure(
                    full_space_campaign->structure_index(batch_index)),
                retained_points[full_space_campaign->structure_index(batch_index)],
                retained_states[batch_index].values,
                retained_states[batch_index].lower,
                retained_states[batch_index].upper,
                result.residuals.data() + offset, jacobian_pointers.data(),
                result.predictions.data() + offset,
                &result.owner_results[batch_index])
        : cost->evaluate(
            final_blocks.data(), result.residuals.data() + offset,
            jacobian_pointers.data(), result.predictions.data() + offset,
            retain_owner_results, &result.owner_results[batch_index]);
    if (!evaluated) {
      result.batch_outcomes[batch_index].status = "rejected";
      const auto &failures = cost->owner_failures();
      result.batch_outcomes[batch_index].failure_code =
          failures.empty() ? "owner_contract_failure" : failures.back().code;
      final_evaluable = false;
      break;
    }
    result.batch_outcomes[batch_index].status = "evaluated";
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
    result.owner_results.clear();
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

bool full_space_eligible(const FitRequest &request) {
  return std::all_of(
      request.experiments.begin(), request.experiments.end(),
      [](const Experiment &value) {
        const auto *experiment = std::get_if<GrepeExperiment>(&value);
        if (experiment == nullptr || experiment->request.temperature.solved ||
            experiment->request.boundary.has_value()) {
          return false;
        }
        const bool finite_only = experiment->request.phases.size() == 1 &&
            !experiment->request.pressure.solved &&
            !experiment->request.phases[0].ideal_gas &&
            experiment->request.phases[0].finite &&
            experiment->request.phases[0].support_name == "all_components" &&
            !experiment->request.phases[0].start.has_value() &&
            experiment->request.continuation_phases.size() == 1 &&
            experiment->request.continuation_phases[0].identity ==
                experiment->request.phases[0].identity;
        const bool finite_incipient = experiment->request.phases.size() == 2 &&
            experiment->request.pressure.solved &&
            std::count_if(
                experiment->request.phases.begin(),
                experiment->request.phases.end(),
                [](const auto &phase) { return phase.finite; }) == 1 &&
            std::all_of(
                experiment->request.phases.begin(),
                experiment->request.phases.end(),
                [](const auto &phase) { return !phase.start.has_value(); });
        if (!finite_only && !finite_incipient) {
          return false;
        }
        return std::all_of(
            experiment->selected_rows.begin(), experiment->selected_rows.end(),
            [&](std::size_t row) {
              return row < experiment->request.outputs.size() &&
                     (experiment->request.outputs[row].selector ==
                          "phase.amount" ||
                      experiment->request.outputs[row].selector ==
                          "phase.mole_fraction" ||
                      (finite_incipient && (
                          experiment->request.outputs[row].selector ==
                              "system.pressure" ||
                          experiment->request.outputs[row].selector ==
                              "phase.partial_pressure")));
            });
      });
}

void validate_request(const FitRequest &request) {
  if (request.active_parameter_count == 0 ||
      request.active_parameter_count != request.identity_coordinates.size() ||
      request.observations.size() != request.experiments.size() ||
      std::any_of(request.identity_coordinates.begin(),
                  request.identity_coordinates.end(),
                  [&](std::size_t coordinate) {
                    return coordinate >= request.coordinates.size();
                  })) {
    throw std::invalid_argument("identity coordinates must bind every active "
                                "identity to one coordinate");
  }
  if (request.coordinates.empty() || request.observations.empty() ||
      request.controls.starts.empty()) {
    throw std::invalid_argument("at least one explicit start is required");
  }
  if (request.failure_codes.empty()) {
    throw std::invalid_argument(
        "the owner failure vocabulary must be nonempty");
  }
  if (request.controls.maximum_iterations <= 0 ||
      !std::isfinite(request.controls.affine_bound_margin) ||
      request.controls.affine_bound_margin < 0.0 ||
      !std::isfinite(request.controls.optimizer_tolerance) ||
      request.controls.optimizer_tolerance <= 0.0) {
    throw std::invalid_argument("solver controls must be positive and finite");
  }
  for (const auto &coordinate : request.coordinates) {
    if (!std::isfinite(coordinate.origin) || !std::isfinite(coordinate.scale) ||
        coordinate.scale <= 0.0 || !std::isfinite(coordinate.lower_affine) ||
        !std::isfinite(coordinate.upper_affine) ||
        coordinate.lower_affine >= coordinate.upper_affine) {
      throw std::invalid_argument("coordinate values and bounds are invalid");
    }
  }
  std::vector<std::string> start_identities;
  start_identities.reserve(request.controls.starts.size());
  for (const auto &start : request.controls.starts) {
    if (start.identity.empty() ||
        start.affine.size() != request.coordinates.size() ||
        !std::all_of(start.affine.begin(), start.affine.end(),
                     [](double value) { return std::isfinite(value); })) {
      throw std::invalid_argument("each explicit start needs an identity and "
                                  "one finite affine value per coordinate");
    }
    for (std::size_t coordinate = 0; coordinate < request.coordinates.size();
         ++coordinate) {
      if (start.affine[coordinate] <
              request.coordinates[coordinate].lower_affine ||
          start.affine[coordinate] >
              request.coordinates[coordinate].upper_affine) {
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
  for (std::size_t index = 0; index < request.observations.size(); ++index) {
    const auto &observation = request.observations[index];
    if (observation.targets.empty()) {
      throw std::invalid_argument("observation targets must be nonempty");
    }
    std::visit(
        [&](const auto &experiment) {
          const auto &active = experiment.model.active_coordinates();
          if (experiment.selected_rows.size() != observation.targets.size()) {
            throw std::invalid_argument(
                "each experiment must select one row per observation target");
          }
          if (active.size() != request.active_parameter_count) {
            throw std::invalid_argument(
                "experiment active coordinates must match the fitted active width");
          }
          using Value = std::decay_t<decltype(experiment)>;
          if constexpr (std::is_same_v<Value, GrepeExperiment>) {
            if (experiment.structure_fingerprint.size() != 71 ||
                experiment.structure_fingerprint.rfind("sha256:", 0) != 0 ||
                experiment.structure_fingerprint.find_first_not_of(
                    "0123456789abcdef", 7) != std::string::npos) {
              throw std::invalid_argument(
                  "GREPE structure fingerprint must be canonical SHA-256");
            }
            auto base = active.begin();
            for (const auto &incipient :
                 experiment.incipient_model.active_coordinates()) {
              base = std::find_if(
                  base, active.end(), [&](const auto &coordinate) {
                    return coordinate.identity == incipient.identity;
                  });
              if (base == active.end()) {
                throw std::invalid_argument(
                    "GREPE incipient active identities must preserve base order");
              }
              ++base;
            }
          }
        },
        request.experiments[index]);
    if (!std::isfinite(observation.multiplier) ||
        observation.multiplier <= 0.0) {
      throw std::invalid_argument(
          "observation multipliers must be positive and finite");
    }
    for (const auto &target : observation.targets) {
      if (!std::isfinite(target.observed) || !std::isfinite(target.scale) ||
          target.scale <= 0.0 || !std::isfinite(target.multiplier) ||
          target.multiplier <= 0.0 ||
          (target.residual == ResidualKind::LogRatio &&
           target.observed <= 0.0)) {
        throw std::invalid_argument(
            "target values, scales, and multipliers are invalid");
      }
    }
  }
}

std::vector<double> start_active_values(const FitRequest &request,
                                        const Start &start) {
  std::vector<double> active(request.active_parameter_count);
  for (std::size_t identity = 0; identity < active.size(); ++identity) {
    const std::size_t coordinate = request.identity_coordinates[identity];
    active[identity] = request.coordinates[coordinate].origin +
                       start.affine[coordinate] *
                           request.coordinates[coordinate].scale;
  }
  return active;
}

} // namespace

FitResult fit_native(const FitRequest &request) {
  validate_request(request);
  FitResult output;
  const bool full_space = full_space_eligible(request);
  output.evaluation_strategy =
      full_space ? "simultaneous_sparse_full_space" : "reduced_space";
  std::optional<FullSpaceCampaign> full_space_campaign;
  if (full_space) {
    full_space_campaign.emplace(request.experiments);
  }

  for (const auto &start : request.controls.starts) {
    const std::vector<double> active = start_active_values(request, start);
    for (std::size_t index = 0; index < request.experiments.size(); ++index) {
      const auto *direct =
          std::get_if<DirectEosExperiment>(&request.experiments[index]);
      if (direct == nullptr) {
        continue;
      }
      try {
        const auto evaluation = evaluate_direct(*direct, active);
        if (evaluation.values.size() != request.observations[index].targets.size() ||
            evaluation.derivatives.size() !=
                evaluation.values.size() * request.active_parameter_count ||
            !std::all_of(evaluation.values.begin(), evaluation.values.end(),
                         [](double value) { return std::isfinite(value); }) ||
            !std::all_of(evaluation.derivatives.begin(),
                         evaluation.derivatives.end(),
                         [](double value) { return std::isfinite(value); })) {
          throw OwnerEvaluationError(
              "owner_contract_failure",
              "native direct EOS preflight returned an incomplete block");
        }
      } catch (const OwnerEvaluationError &error) {
        output.preflight_failure = PreflightFailure{
            error.code(), start.identity, "observation-" + std::to_string(index),
            error.what()};
        return output;
      }
    }
  }

  output.starts.resize(request.controls.starts.size());
  for (std::size_t index = 0; index < request.controls.starts.size(); ++index) {
    try {
      output.starts[index] = solve_start(
          index, request, full_space,
          full_space ? &*full_space_campaign : nullptr,
          request.identity_coordinates,
          request.coordinates, request.observations, request.controls,
          request.controls.starts[index], request.failure_codes);
    } catch (const std::exception &error) {
      auto &result = output.starts[index];
      result.index = index;
      result.start_identity = request.controls.starts[index].identity;
      result.initial_affine = request.controls.starts[index].affine;
      result.final_affine = request.controls.starts[index].affine;
      result.failure_reason = error.what();
      result.trial_failure_code = "optimizer_failure";
      for (std::size_t batch = 0; batch < request.observations.size(); ++batch) {
        result.batch_outcomes.push_back(
            {batch, request.observations[batch].targets.size(), "not_evaluated",
             "optimizer_failure"});
      }
    }
  }
  assign_optimizer_basins(output.starts, request.coordinates);
  return output;
}

} // namespace epcsaft::regression
