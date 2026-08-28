#pragma once

#include "eos.hpp"
#include "grepe.hpp"

#include <cstddef>
#include <limits>
#include <optional>
#include <string>
#include <variant>
#include <vector>

namespace epcsaft::regression {

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

enum class DirectEosKind { FixedState, PureSaturation };

struct DirectEosExperiment final {
  NativeModel model;
  DirectEosKind kind;
  double temperature_k;
  double molar_density_mol_m3;
  std::vector<double> composition;
  std::vector<std::size_t> selected_rows;
};

struct GrepeExperiment final {
  std::string structure_fingerprint;
  NativeModel model;
  NativeModel incipient_model;
  epcsaft_engine::equilibrium::GrepeRequest request;
  std::vector<std::size_t> selected_rows;
};

using Experiment = std::variant<DirectEosExperiment, GrepeExperiment>;

struct FitRequest final {
  std::size_t active_parameter_count;
  std::vector<std::size_t> identity_coordinates;
  std::vector<Coordinate> coordinates;
  std::vector<Observation> observations;
  Controls controls;
  std::vector<std::string> failure_codes;
  std::vector<Experiment> experiments;
  std::string evaluation_strategy{"automatic"};
};

struct OwnerFailure final {
  std::string code;
  std::string reason;
  std::size_t count{0};
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

struct ParameterDomainConstraintDiagnostic final {
  std::string identity;
  double temperature_k;
  double initial_value;
  double final_value;
  double upper;
  std::vector<double> affine_jacobian;
};

struct KktDiagnostics final {
  bool retained{false};
  std::vector<double> coordinate_lower_bound_multipliers;
  std::vector<double> coordinate_upper_bound_multipliers;
  std::vector<double> parameter_constraint_multipliers;
  double objective{std::numeric_limits<double>::quiet_NaN()};
  double primal_feasibility_inf_norm{std::numeric_limits<double>::quiet_NaN()};
  double stationarity_inf_norm{std::numeric_limits<double>::quiet_NaN()};
  double complementarity_inf_norm{std::numeric_limits<double>::quiet_NaN()};
  double dual_sign_violation_inf_norm{std::numeric_limits<double>::quiet_NaN()};
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
  std::vector<ParameterDomainConstraintDiagnostic>
      parameter_domain_constraints;
  std::vector<BatchOutcome> batch_outcomes;
  std::vector<std::optional<epcsaft_engine::equilibrium::GrepeResult>>
      owner_results;
  Diagnostics diagnostics;
  KktDiagnostics kkt_diagnostics;
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

struct PreflightFailure final {
  std::string code;
  std::string start_identity;
  std::string batch_identity;
  std::string diagnostic;
};

struct FitResult final {
  std::vector<StartResult> starts;
  std::optional<PreflightFailure> preflight_failure;
  std::string evaluation_strategy{"reduced_space"};
};

[[nodiscard]] FitResult fit_native(const FitRequest &request);

} // namespace epcsaft::regression
