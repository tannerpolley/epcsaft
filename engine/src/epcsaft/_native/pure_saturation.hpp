#pragma once

#include <array>
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
  double objective = 0.0;
  std::array<double, 3> objective_gradient{};
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

struct SaturationSolveResult {
  bool accepted = false;
  bool solver_converged = false;
  bool numerical_converged = false;
  bool physical_accepted = false;
  std::string failure_reason;
  std::string solver_status;
  int iterations = 0;
  int attempts = 0;
  int confirmation_solves = 0;
  std::array<double, 3> solver_lower_bounds{};
  std::array<double, 3> solver_upper_bounds{};
  struct AttemptRecord {
    std::string role;
    std::array<double, 3> initial_guess{};
    bool solver_converged = false;
    std::string solver_status;
    int iterations = 0;
    double constraint_violation = 0.0;
    std::string callback_error;
  };
  std::vector<AttemptRecord> attempt_log;
  double solver_constraint_violation = 0.0;
  double confirmation_max_relative_difference = 0.0;
  double pressure_relative_residual = 0.0;
  double chemical_potential_absolute_residual = 0.0;
  double phase_density_distance = 0.0;
  NlpEvaluation evaluation{};
};

struct SaturationSensitivityResult {
  SaturationSolveResult solve;
  std::vector<std::string> parameter_identities;
  std::vector<double> derivatives;
};

[[nodiscard]] SaturationSolveResult
solve_local_saturation(const NativeModel &model, double temperature_k,
                       double liquid_density_upper_mol_m3);

[[nodiscard]] SaturationSensitivityResult solve_local_saturation_active(
    const NativeModel &model, double temperature_k,
    double liquid_density_upper_mol_m3,
    const std::vector<NativeParameterRecord> &active_parameters);

} // namespace epcsaft
