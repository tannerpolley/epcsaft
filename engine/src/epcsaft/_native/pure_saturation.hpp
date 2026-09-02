#pragma once

#include <array>
#include <limits>
#include <optional>
#include <string>
#include <vector>

#include "eos.hpp"
#include "phase_block.hpp"

namespace epcsaft {

struct SaturationPhaseEvaluation {
  double amount_mol = 0.0;
  double volume_m3 = 0.0;
  NativePurePhaseBlock eos{};
};

struct NlpEvaluation {
  std::array<double, 3> constraints{};
  std::array<double, 9> jacobian{};
  std::array<double, 6> lagrangian_hessian_lower{};
  SaturationPhaseEvaluation vapor{};
  SaturationPhaseEvaluation liquid{};
  double saturation_pressure_pa = 0.0;
};

[[nodiscard]] NlpEvaluation
evaluate_saturation_nlp(const NativeModel &model, double temperature_k,
                        const std::array<double, 3> &variables,
                        const std::array<double, 3> &multipliers);

struct SaturationBounds {
  std::array<double, 3> solver_lower{};
  std::array<double, 3> solver_upper{};
  std::array<double, 2> vapor_density_mol_m3{};
  std::array<double, 2> liquid_density_mol_m3{};
};

struct SaturationSolveResult {
  bool accepted = false;
  bool candidate_available = false;
  bool optimizer_converged = false;
  std::string optimizer_status;
  SaturationBounds bounds{};
  double feasibility_measure = std::numeric_limits<double>::infinity();
  double feasibility_tolerance = std::numeric_limits<double>::infinity();
  double pressure_relative_residual = std::numeric_limits<double>::infinity();
  double vapor_pressure_relative_residual =
      std::numeric_limits<double>::infinity();
  double liquid_pressure_relative_residual =
      std::numeric_limits<double>::infinity();
  double chemical_potential_absolute_residual =
      std::numeric_limits<double>::infinity();
  double phase_density_distance = std::numeric_limits<double>::quiet_NaN();
  double vapor_residual_density_derivative_m3_per_mol =
      std::numeric_limits<double>::quiet_NaN();
  double liquid_residual_density_derivative_m3_per_mol =
      std::numeric_limits<double>::quiet_NaN();
  int physical_status_count = 0;
  std::string physical_status = "not_adjudicated";
  std::string rejection_status = "solve-failure";
  std::string failure_reason;
  NlpEvaluation evaluation{};
};

struct SaturationSensitivityResult {
  SaturationSolveResult solve;
  std::vector<std::string> parameter_identities;
  std::vector<double> derivatives;
  std::string sensitivity_status = "unavailable";
  std::string sensitivity_failure_reason = "not_requested";
  std::optional<std::size_t> sensitivity_rank;
  std::optional<std::size_t> sensitivity_dimension;
  std::optional<double> sensitivity_condition_number;
  std::optional<double> sensitivity_backward_error;
};

[[nodiscard]] SaturationSolveResult
solve_local_saturation(const NativeModel &model, double temperature_k,
                       double liquid_density_upper_mol_m3,
                       std::optional<std::array<double, 3>> explicit_start =
                           std::nullopt);

[[nodiscard]] SaturationSensitivityResult solve_local_saturation_active(
    const NativeModel &model, double temperature_k,
    const std::vector<NativeActiveCoordinate> &active_parameters);

} // namespace epcsaft
