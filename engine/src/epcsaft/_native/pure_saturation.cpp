#include "pure_saturation.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

namespace epcsaft {

namespace {

constexpr double kGasConstantJPerMolK = 8.31446261815324;
constexpr double kVaporDensityLower = 1.0e-8;
constexpr double kVaporDensityUpper = 5.0e3;
constexpr double kLiquidDensityLower = 5.001e3;
constexpr double kPressureLowerPa = 1.0e-12;
constexpr double kPressureUpperPa = 1.0e8;
constexpr double kSolverConstraintTolerance = 1.0e-10;
constexpr double kConstraintTolerance = 1.0e-8;
constexpr double kPhaseDistanceTolerance = 1.0e-3;

void require_finite(double value, const char *name) {
  if (!std::isfinite(value)) {
    throw std::invalid_argument(std::string(name) + " must be finite");
  }
}

struct PhaseDerivatives {
  SaturationPhaseEvaluation phase;
  double pressure_q;
  double pressure_qq;
  double chemical_potential_q;
  double chemical_potential_qq;
};

PhaseDerivatives evaluate_phase_derivatives(const NativeModel &model,
                                            double temperature_k,
                                            double log_density) {
  const double volume_m3 = std::exp(-log_density);
  require_finite(volume_m3, "phase volume");
  if (volume_m3 <= 0.0) {
    throw std::invalid_argument("phase volume must be positive");
  }
  SaturationPhaseEvaluation phase{
      1.0,
      volume_m3,
      model.evaluate_pure_phase(temperature_k, 1.0, volume_m3),
  };
  const double rt = kGasConstantJPerMolK * temperature_k;
  const double h_nv = phase.eos.hessian[1];
  const double h_vv = phase.eos.hessian[3];
  const double t_nvv = phase.eos.third[3];
  const double t_vvv = phase.eos.third[7];
  return {
      std::move(phase),
      rt * h_vv * volume_m3,
      -rt * (t_vvv * volume_m3 * volume_m3 + h_vv * volume_m3),
      -h_nv * volume_m3,
      h_nv * volume_m3 + t_nvv * volume_m3 * volume_m3,
  };
}

struct AttemptResult {
  bool solver_converged = false;
  std::string solver_status;
  std::string callback_error;
  int iterations = 0;
  double solver_constraint_violation = std::numeric_limits<double>::infinity();
  std::array<double, 3> initial_guess{};
  std::array<double, 3> variables{};
};

std::array<double, 3> solver_lower_bounds() {
  return {
      std::log(kVaporDensityLower),
      std::log(kLiquidDensityLower),
      std::log(kPressureLowerPa),
  };
}

std::array<double, 3> solver_upper_bounds(double liquid_density_upper_mol_m3) {
  return {
      std::log(kVaporDensityUpper),
      std::log(liquid_density_upper_mol_m3),
      std::log(kPressureUpperPa),
  };
}

NlpEvaluation
evaluate_saturation_state(const NativeModel &model, double temperature_k,
                          const std::array<double, 3> &variables) {
  auto phase = [&](double log_density) {
    const double volume = std::exp(-log_density);
    const auto block = model.evaluate_fixed_topology_association(
        temperature_k, 1.0, volume, {});
    SaturationPhaseEvaluation value{
        1.0,
        volume,
        {
            block.helmholtz_over_rt_reference_amount,
            {block.gradient[0], block.gradient[1]},
            {block.hessian[0], block.hessian[1], block.hessian[2],
             block.hessian[3]},
            {},
            block.pressure_pa,
            block.chemical_potential_over_rt,
        },
    };
    const double pressure_q =
        kGasConstantJPerMolK * temperature_k * block.hessian[3] * volume;
    const double chemical_potential_q = -block.hessian[1] * volume;
    return std::tuple{std::move(value), pressure_q, chemical_potential_q};
  };
  auto [vapor, vapor_pressure_q, vapor_mu_q] = phase(variables[0]);
  auto [liquid, liquid_pressure_q, liquid_mu_q] = phase(variables[1]);
  const double pressure = std::exp(variables[2]);
  const double vapor_ratio = vapor.eos.pressure_pa / pressure;
  const double liquid_ratio = liquid.eos.pressure_pa / pressure;
  return {
      0.0,
      {},
      {
          vapor_ratio - 1.0,
          liquid_ratio - 1.0,
          vapor.eos.chemical_potential_over_rt -
              liquid.eos.chemical_potential_over_rt,
      },
      {
          vapor_pressure_q / pressure,
          0.0,
          -vapor_ratio,
          0.0,
          liquid_pressure_q / pressure,
          -liquid_ratio,
          vapor_mu_q,
          -liquid_mu_q,
          0.0,
      },
      {},
      std::move(vapor),
      std::move(liquid),
      pressure,
  };
}

AttemptResult run_saturation_newton(const NativeModel &model,
                                    double temperature_k,
                                    const std::array<double, 3> &initial,
                                    double liquid_density_upper_mol_m3) {
  AttemptResult result;
  result.initial_guess = initial;
  result.variables = initial;
  const auto lower = solver_lower_bounds();
  const auto upper = solver_upper_bounds(liquid_density_upper_mol_m3);
  try {
    for (int iteration = 0; iteration < 80; ++iteration) {
      const auto evaluation = evaluate_saturation_state(
          model, temperature_k, result.variables);
      const double residual = std::max({std::abs(evaluation.constraints[0]),
                                        std::abs(evaluation.constraints[1]),
                                        std::abs(evaluation.constraints[2])});
      result.iterations = iteration;
      result.solver_constraint_violation = residual;
      if (residual <= kSolverConstraintTolerance) {
        result.solver_converged = true;
        result.solver_status = "solve_succeeded";
        return result;
      }
      std::vector<double> right_hand_side{-evaluation.constraints[0],
                                          -evaluation.constraints[1],
                                          -evaluation.constraints[2]};
      const std::vector<double> jacobian(evaluation.jacobian.begin(),
                                         evaluation.jacobian.end());
      const auto step =
          detail::solve_linear(jacobian, std::move(right_hand_side));
      bool accepted = false;
      for (double scale = 1.0; scale >= 1.0e-6; scale *= 0.5) {
        std::array<double, 3> trial{};
        bool bounded = true;
        for (std::size_t coordinate = 0; coordinate < trial.size();
             ++coordinate) {
          trial[coordinate] =
              result.variables[coordinate] + scale * step[coordinate];
          bounded = bounded && trial[coordinate] > lower[coordinate] &&
                    trial[coordinate] < upper[coordinate];
        }
        if (!bounded) {
          continue;
        }
        const auto candidate =
            evaluate_saturation_state(model, temperature_k, trial);
        const double candidate_residual =
            std::max({std::abs(candidate.constraints[0]),
                      std::abs(candidate.constraints[1]),
                      std::abs(candidate.constraints[2])});
        if (candidate_residual < residual) {
          result.variables = trial;
          accepted = true;
          break;
        }
      }
      if (!accepted) {
        result.solver_status = "line_search_failed";
        return result;
      }
    }
    result.solver_status = "maximum_iterations_exceeded";
  } catch (const std::exception &error) {
    result.callback_error = error.what();
    result.solver_status = "evaluation_failed";
  }
  return result;
}

double relative_difference(double first, double second) {
  return std::abs(first - second) /
         std::max({std::abs(first), std::abs(second), 1.0e-300});
}

bool inside_bounds(const std::array<double, 3> &variables,
                   double liquid_density_upper_mol_m3) {
  constexpr double kMargin = 1.0e-7;
  const std::array<double, 3> lower = solver_lower_bounds();
  const std::array<double, 3> upper =
      solver_upper_bounds(liquid_density_upper_mol_m3);
  for (std::size_t index = 0; index < variables.size(); ++index) {
    if (variables[index] <= lower[index] + kMargin ||
        variables[index] >= upper[index] - kMargin) {
      return false;
    }
  }
  return true;
}

bool physical_acceptance(const AttemptResult &attempt,
                         const NlpEvaluation &evaluation,
                         double liquid_density_upper_mol_m3,
                         SaturationSolveResult &result) {
  const double vapor_density = 1.0 / evaluation.vapor.volume_m3;
  const double liquid_density = 1.0 / evaluation.liquid.volume_m3;
  result.pressure_relative_residual = std::max(
      std::abs(evaluation.constraints[0]), std::abs(evaluation.constraints[1]));
  result.chemical_potential_absolute_residual =
      std::abs(evaluation.constraints[2]);
  result.phase_density_distance =
      (liquid_density - vapor_density) / liquid_density;
  const bool finite =
      std::isfinite(evaluation.saturation_pressure_pa) &&
      std::isfinite(vapor_density) && std::isfinite(liquid_density) &&
      std::isfinite(result.pressure_relative_residual) &&
      std::isfinite(result.chemical_potential_absolute_residual) &&
      std::isfinite(result.phase_density_distance);
  const bool locally_stable =
      evaluation.jacobian[0] > 0.0 && evaluation.jacobian[4] > 0.0;
  return attempt.solver_converged &&
         std::isfinite(attempt.solver_constraint_violation) &&
         attempt.solver_constraint_violation <= kSolverConstraintTolerance &&
         attempt.callback_error.empty() && finite &&
         evaluation.saturation_pressure_pa > 0.0 && vapor_density > 0.0 &&
         liquid_density > vapor_density &&
         result.pressure_relative_residual <= kConstraintTolerance &&
         result.chemical_potential_absolute_residual <= kConstraintTolerance &&
         result.phase_density_distance > kPhaseDistanceTolerance &&
         locally_stable &&
         inside_bounds(attempt.variables, liquid_density_upper_mol_m3);
}

SaturationSolveResult solve_saturation(const NativeModel &model,
                                       double temperature_k,
                                       double liquid_density_upper_mol_m3) {
  const std::array<double, 7> vapor_seeds{1.0e-5, 1.0e-3, 0.1,   10.0,
                                          100.0,  1000.0, 3000.0};
  const std::array<double, 3> liquid_seeds{18'000.0, 12'000.0, 24'000.0};
  SaturationSolveResult last;
  last.failure_reason = "no Newton seed produced a certified local boundary";
  last.solver_lower_bounds = solver_lower_bounds();
  last.solver_upper_bounds = solver_upper_bounds(liquid_density_upper_mol_m3);
  for (double vapor_density : vapor_seeds) {
    for (double liquid_density : liquid_seeds) {
      const std::array<double, 3> initial{
          std::log(vapor_density),
          std::log(liquid_density),
          std::log(
              std::max(vapor_density * kGasConstantJPerMolK * temperature_k,
                       kPressureLowerPa)),
      };
      AttemptResult attempt = run_saturation_newton(
          model, temperature_k, initial, liquid_density_upper_mol_m3);
      ++last.attempts;
      last.attempt_log.push_back({
          "primary",
          attempt.initial_guess,
          attempt.solver_converged,
          attempt.solver_status,
          attempt.iterations,
          attempt.solver_constraint_violation,
          attempt.callback_error,
      });
      if (!attempt.solver_converged) {
        last.failure_reason = attempt.callback_error.empty()
                                  ? attempt.solver_status
                                  : attempt.callback_error;
        continue;
      }
      NlpEvaluation evaluation = evaluate_saturation_state(
          model, temperature_k, attempt.variables);
      SaturationSolveResult candidate = last;
      candidate.solver_converged = true;
      candidate.solver_status = attempt.solver_status;
      candidate.iterations = attempt.iterations;
      candidate.solver_constraint_violation =
          attempt.solver_constraint_violation;
      candidate.evaluation = evaluation;
      candidate.physical_accepted = physical_acceptance(
          attempt, evaluation, liquid_density_upper_mol_m3, candidate);
      if (!candidate.physical_accepted) {
        last = candidate;
        last.failure_reason =
            "Newton solution failed local physical acceptance";
        continue;
      }
      const std::array<std::array<double, 3>, 2> perturbations{{
          {0.05, -0.05, 0.02},
          {-0.05, 0.05, -0.02},
      }};
      bool confirmed = true;
      double maximum_difference = 0.0;
      for (const auto &perturbation : perturbations) {
        std::array<double, 3> seed{};
        for (std::size_t index = 0; index < seed.size(); ++index) {
          seed[index] = attempt.variables[index] + perturbation[index];
        }
        AttemptResult confirmation = run_saturation_newton(
            model, temperature_k, seed, liquid_density_upper_mol_m3);
        ++candidate.attempts;
        candidate.attempt_log.push_back({
            "confirmation",
            confirmation.initial_guess,
            confirmation.solver_converged,
            confirmation.solver_status,
            confirmation.iterations,
            confirmation.solver_constraint_violation,
            confirmation.callback_error,
        });
        if (!confirmation.solver_converged) {
          confirmed = false;
          break;
        }
        const auto check = evaluate_saturation_state(
            model, temperature_k, confirmation.variables);
        maximum_difference =
            std::max({maximum_difference,
                      relative_difference(evaluation.saturation_pressure_pa,
                                          check.saturation_pressure_pa),
                      relative_difference(evaluation.vapor.volume_m3,
                                          check.vapor.volume_m3),
                      relative_difference(evaluation.liquid.volume_m3,
                                          check.liquid.volume_m3)});
        ++candidate.confirmation_solves;
      }
      candidate.confirmation_max_relative_difference = maximum_difference;
      candidate.numerical_converged = confirmed && maximum_difference <= 1.0e-7;
      if (!candidate.numerical_converged) {
        last = candidate;
        last.failure_reason =
            "perturbed Newton solves did not confirm the boundary";
        continue;
      }
      candidate.accepted = true;
      candidate.failure_reason.clear();
      return candidate;
    }
  }
  return last;
}

} // namespace

NlpEvaluation
evaluate_saturation_nlp(const NativeModel &model, double temperature_k,
                        const std::array<double, 3> &variables,
                        const std::array<double, 3> &multipliers) {
  require_finite(temperature_k, "temperature");
  if (temperature_k <= 0.0) {
    throw std::invalid_argument("temperature must be positive");
  }
  for (double value : variables) {
    require_finite(value, "NLP variable");
  }
  for (double value : multipliers) {
    require_finite(value, "constraint multiplier");
  }

  PhaseDerivatives vapor =
      evaluate_phase_derivatives(model, temperature_k, variables[0]);
  PhaseDerivatives liquid =
      evaluate_phase_derivatives(model, temperature_k, variables[1]);
  const double saturation_pressure_pa = std::exp(variables[2]);
  require_finite(saturation_pressure_pa, "saturation pressure");
  if (saturation_pressure_pa <= 0.0) {
    throw std::invalid_argument("saturation pressure must be positive");
  }

  const double vapor_pressure_ratio =
      vapor.phase.eos.pressure_pa / saturation_pressure_pa;
  const double liquid_pressure_ratio =
      liquid.phase.eos.pressure_pa / saturation_pressure_pa;
  const double vapor_pressure_q_ratio =
      vapor.pressure_q / saturation_pressure_pa;
  const double liquid_pressure_q_ratio =
      liquid.pressure_q / saturation_pressure_pa;

  NlpEvaluation result{
      0.0,
      {0.0, 0.0, 0.0},
      {
          vapor_pressure_ratio - 1.0,
          liquid_pressure_ratio - 1.0,
          vapor.phase.eos.chemical_potential_over_rt -
              liquid.phase.eos.chemical_potential_over_rt,
      },
      {
          vapor_pressure_q_ratio,
          0.0,
          -vapor_pressure_ratio,
          0.0,
          liquid_pressure_q_ratio,
          -liquid_pressure_ratio,
          vapor.chemical_potential_q,
          -liquid.chemical_potential_q,
          0.0,
      },
      {
          multipliers[0] * vapor.pressure_qq / saturation_pressure_pa +
              multipliers[2] * vapor.chemical_potential_qq,
          0.0,
          multipliers[1] * liquid.pressure_qq / saturation_pressure_pa -
              multipliers[2] * liquid.chemical_potential_qq,
          -multipliers[0] * vapor_pressure_q_ratio,
          -multipliers[1] * liquid_pressure_q_ratio,
          multipliers[0] * vapor_pressure_ratio +
              multipliers[1] * liquid_pressure_ratio,
      },
      std::move(vapor.phase),
      std::move(liquid.phase),
      saturation_pressure_pa,
  };
  return result;
}

SaturationSolveResult
solve_local_saturation(const NativeModel &model, double temperature_k,
                       double liquid_density_upper_mol_m3) {
  require_finite(temperature_k, "temperature");
  if (temperature_k <= 0.0) {
    throw std::invalid_argument("temperature must be positive");
  }
  require_finite(liquid_density_upper_mol_m3, "liquid density upper bound");
  if (liquid_density_upper_mol_m3 <= kLiquidDensityLower) {
    throw std::invalid_argument(
        "liquid density upper bound must exceed 5001 mol m^-3");
  }
  return solve_saturation(model, temperature_k,
                          liquid_density_upper_mol_m3);
}

SaturationSensitivityResult solve_local_saturation_active(
    const NativeModel &model, double temperature_k,
    double liquid_density_upper_mol_m3,
    const std::vector<NativeParameterRecord> &active_parameters) {
  if (active_parameters.empty()) {
    throw std::invalid_argument(
        "active saturation parameter set must be nonempty");
  }
  SaturationSensitivityResult result;
  const NativeModel active_model = model.with_active_parameters(active_parameters);
  result.solve = solve_local_saturation(
      active_model, temperature_k, liquid_density_upper_mol_m3);
  if (!result.solve.accepted) {
    throw std::runtime_error(
        "saturation sensitivity requires an accepted solve");
  }
  const std::size_t count = active_parameters.size();
  result.parameter_identities.reserve(count);
  for (const auto &parameter : active_parameters) {
    result.parameter_identities.push_back(parameter.identity);
  }

  auto phase_partials = [&](double volume) {
    std::vector<double> pressure(count);
    std::vector<double> chemical_potential(count);
    const auto block =
        active_model.evaluate_fixed_topology_association_active_parameters(
            temperature_k, 1.0, volume, active_parameters);
    const std::size_t width = count + 2;
    for (std::size_t parameter = 0; parameter < count; ++parameter) {
      const std::size_t coordinate = parameter + 2;
      pressure[parameter] = -kGasConstantJPerMolK * temperature_k *
                            block.hessian[width + coordinate];
      chemical_potential[parameter] = block.hessian[coordinate];
    }
    return std::pair{std::move(pressure), std::move(chemical_potential)};
  };

  const auto vapor = phase_partials(result.solve.evaluation.vapor.volume_m3);
  const auto liquid = phase_partials(result.solve.evaluation.liquid.volume_m3);
  const double saturation_pressure =
      result.solve.evaluation.saturation_pressure_pa;
  const std::vector<double> jacobian(result.solve.evaluation.jacobian.begin(),
                                     result.solve.evaluation.jacobian.end());
  result.derivatives.resize(3 * count);
  const double vapor_density = 1.0 / result.solve.evaluation.vapor.volume_m3;
  const double liquid_density = 1.0 / result.solve.evaluation.liquid.volume_m3;
  for (std::size_t parameter = 0; parameter < count; ++parameter) {
    std::vector<double> right_hand_side{
        -vapor.first[parameter] / saturation_pressure,
        -liquid.first[parameter] / saturation_pressure,
        -vapor.second[parameter] + liquid.second[parameter],
    };
    const auto state =
        detail::solve_linear(jacobian, std::move(right_hand_side));
    result.derivatives[parameter] = saturation_pressure * state[2];
    result.derivatives[count + parameter] = vapor_density * state[0];
    result.derivatives[2 * count + parameter] = liquid_density * state[1];
  }
  return result;
}

} // namespace epcsaft
