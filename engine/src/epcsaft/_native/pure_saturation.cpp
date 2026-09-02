#include "pure_saturation.hpp"
#include "physical_constants.hpp"
#include "sensitivity_linear_solve.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>
#include <IpIpoptApplication.hpp>
#include <IpTNLP.hpp>
namespace epcsaft {
namespace {
using epcsaft::detail::GAS_CONSTANT;
constexpr double kVaporDensityLower = 1.0e-8;
constexpr double kVaporDensityUpper = 5.0e3;
constexpr double kLiquidDensityLower = 5.001e3;
constexpr double kPressureLowerPa = 1.0e-12;
constexpr double kPressureUpperPa = 1.0e8;
// Dimensionless closure bound: 2e-8 is the retained roundoff/conditioning scale.
constexpr double kFeasibilityTolerance = 2.0e-8;
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
  SaturationPhaseEvaluation phase{1.0, volume_m3,
                                 model.evaluate_pure_phase(temperature_k, 1.0, volume_m3)};
  const double rt = GAS_CONSTANT * temperature_k;
  const double h_nv = phase.eos.hessian[1];
  const double h_vv = phase.eos.hessian[3];
  const double t_nvv = phase.eos.third[3];
  const double t_vvv = phase.eos.third[7];
  return {std::move(phase), rt * h_vv * volume_m3,
          -rt * (t_vvv * volume_m3 * volume_m3 + h_vv * volume_m3),
          -h_nv * volume_m3, h_nv * volume_m3 + t_nvv * volume_m3 * volume_m3};
}
struct AttemptResult {
  bool solution_available = false;
  bool optimizer_converged = false;
  std::string optimizer_status;
  std::string callback_error;
  std::array<double, 3> variables{};
};
SaturationBounds saturation_bounds(double liquid_density_upper_mol_m3) {
  return {
      {std::log(kVaporDensityLower), std::log(kLiquidDensityLower),
       std::log(kPressureLowerPa)},
      {std::log(kVaporDensityUpper), std::log(liquid_density_upper_mol_m3),
       std::log(kPressureUpperPa)},
      {kVaporDensityLower, kVaporDensityUpper},
      {kLiquidDensityLower, liquid_density_upper_mol_m3},
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
        GAS_CONSTANT * temperature_k * block.hessian[3] * volume;
    const double chemical_potential_q = -block.hessian[1] * volume;
    return std::tuple{std::move(value), pressure_q, chemical_potential_q};
  };
  auto [vapor, vapor_pressure_q, vapor_mu_q] = phase(variables[0]);
  auto [liquid, liquid_pressure_q, liquid_mu_q] = phase(variables[1]);
  const double pressure = std::exp(variables[2]);
  const double vapor_ratio = vapor.eos.pressure_pa / pressure;
  const double liquid_ratio = liquid.eos.pressure_pa / pressure;
  return {{vapor_ratio - 1.0, liquid_ratio - 1.0,
           vapor.eos.chemical_potential_over_rt - liquid.eos.chemical_potential_over_rt},
          {vapor_pressure_q / pressure, 0.0, -vapor_ratio, 0.0,
           liquid_pressure_q / pressure, -liquid_ratio, vapor_mu_q, -liquid_mu_q, 0.0},
          {}, std::move(vapor), std::move(liquid), pressure};
}
std::string ipopt_status_name(Ipopt::ApplicationReturnStatus status) {
  return status == Ipopt::Solve_Succeeded
             ? "solve_succeeded"
             : "ipopt_status_" + std::to_string(static_cast<int>(status));
}
class SaturationTnlp final : public Ipopt::TNLP {
public:
  SaturationTnlp(const NativeModel &model, double temperature_k,
                 std::array<double, 3> initial,
                 SaturationBounds bounds,
                 bool exact_hessian)
      : model_(model), temperature_k_(temperature_k), initial_(initial),
        bounds_(std::move(bounds)), exact_hessian_(exact_hessian) {}
  bool get_nlp_info(Ipopt::Index &n, Ipopt::Index &m, Ipopt::Index &nnz_jac_g,
                    Ipopt::Index &nnz_h_lag, IndexStyleEnum &style) override {
    n = m = 3; nnz_jac_g = 9; nnz_h_lag = 6; style = TNLP::C_STYLE; return true;
  }
  bool get_bounds_info(Ipopt::Index, Ipopt::Number *x_l, Ipopt::Number *x_u,
                       Ipopt::Index m, Ipopt::Number *g_l, Ipopt::Number *g_u) override {
    std::copy(bounds_.solver_lower.begin(), bounds_.solver_lower.end(), x_l);
    std::copy(bounds_.solver_upper.begin(), bounds_.solver_upper.end(), x_u);
    std::fill(g_l, g_l + m, 0.0); std::fill(g_u, g_u + m, 0.0); return true;
  }
  bool get_starting_point(Ipopt::Index, bool init_x, Ipopt::Number *x, bool init_z,
                          Ipopt::Number *, Ipopt::Number *, Ipopt::Index, bool init_lambda,
                          Ipopt::Number *) override {
    if (!init_x || init_z || init_lambda) return false;
    std::copy(initial_.begin(), initial_.end(), x);
    return true;
  }
  bool eval_f(Ipopt::Index, const Ipopt::Number *, bool, Ipopt::Number &objective) override { objective = 0.0; return true; }
  bool eval_grad_f(Ipopt::Index n, const Ipopt::Number *, bool, Ipopt::Number *gradient) override { std::fill(gradient, gradient + n, 0.0); return true; }
  bool eval_g(Ipopt::Index, const Ipopt::Number *x, bool, Ipopt::Index,
              Ipopt::Number *constraints) override {
    NlpEvaluation evaluation; if (!evaluate(x, evaluation)) return false;
    std::copy(evaluation.constraints.begin(), evaluation.constraints.end(), constraints); return true;
  }
  bool eval_jac_g(Ipopt::Index, const Ipopt::Number *x, bool, Ipopt::Index,
                  Ipopt::Index count, Ipopt::Index *rows, Ipopt::Index *columns,
                  Ipopt::Number *values) override {
    if (count != 9) return false;
    if (values == nullptr) { for (Ipopt::Index i = 0; i < 9; ++i) { rows[i] = i / 3; columns[i] = i % 3; } return true; }
    NlpEvaluation evaluation; if (!evaluate(x, evaluation)) return false;
    std::copy(evaluation.jacobian.begin(), evaluation.jacobian.end(), values); return true;
  }
  bool eval_h(Ipopt::Index, const Ipopt::Number *x, bool, Ipopt::Number,
              Ipopt::Index, const Ipopt::Number *lambda, bool, Ipopt::Index count,
              Ipopt::Index *rows, Ipopt::Index *columns, Ipopt::Number *values) override {
    if (!exact_hessian_ || count != 6) return false;
    if (values == nullptr) { std::size_t i = 0; for (Ipopt::Index r = 0; r < 3; ++r) for (Ipopt::Index c = 0; c <= r; ++c) { rows[i] = r; columns[i++] = c; } return true; }
    try { const auto e = evaluate_saturation_nlp(model_, temperature_k_, {x[0], x[1], x[2]}, {lambda[0], lambda[1], lambda[2]}); std::copy(e.lagrangian_hessian_lower.begin(), e.lagrangian_hessian_lower.end(), values); return true; }
    catch (const std::exception &error) { callback_error_ = error.what(); return false; }
  }
  void finalize_solution(Ipopt::SolverReturn, Ipopt::Index n, const Ipopt::Number *x,
                         const Ipopt::Number *, const Ipopt::Number *, Ipopt::Index,
                         const Ipopt::Number *, const Ipopt::Number *, Ipopt::Number,
                         const Ipopt::IpoptData *, Ipopt::IpoptCalculatedQuantities *) override {
    if (x) {
      std::copy(x, x + n, solution_.begin());
      solution_available_ = true;
    }
  }
  [[nodiscard]] const std::array<double, 3> &solution() const { return solution_; }
  [[nodiscard]] bool solution_available() const { return solution_available_; }
  [[nodiscard]] const std::string &callback_error() const { return callback_error_; }
private:
  bool evaluate(const Ipopt::Number *x, NlpEvaluation &result) {
    try { result = evaluate_saturation_state(model_, temperature_k_, {x[0], x[1], x[2]}); return true; }
    catch (const std::exception &error) { callback_error_ = error.what(); return false; }
  }
  const NativeModel &model_;
  double temperature_k_;
  std::array<double, 3> initial_;
  SaturationBounds bounds_;
  bool exact_hessian_;
  std::array<double, 3> solution_{};
  bool solution_available_ = false;
  std::string callback_error_;
};
AttemptResult run_saturation_ipopt(const NativeModel &model, double temperature_k,
                                   const std::array<double, 3> &initial,
                                   const SaturationBounds &bounds) {
  const bool exact_hessian = model.site_components().empty();
  auto *raw_problem = new SaturationTnlp(
      model, temperature_k, initial, bounds, exact_hessian);
  Ipopt::SmartPtr<Ipopt::TNLP> problem = raw_problem;
  Ipopt::SmartPtr<Ipopt::IpoptApplication> application = IpoptApplicationFactory();
  application->Options()->SetIntegerValue("print_level", 0);
  application->Options()->SetStringValue("sb", "yes");
  application->Options()->SetIntegerValue("max_iter", 50);
  application->Options()->SetNumericValue("tol", kFeasibilityTolerance);
  application->Options()->SetNumericValue("acceptable_tol", kFeasibilityTolerance);
  application->Options()->SetNumericValue("acceptable_constr_viol_tol",
                                          kFeasibilityTolerance);
  application->Options()->SetNumericValue("constr_viol_tol", kFeasibilityTolerance);
  application->Options()->SetStringValue("jacobian_approximation", "exact");
  application->Options()->SetStringValue("hessian_approximation", exact_hessian ? "exact" : "limited-memory");
  application->Options()->SetStringValue("nlp_scaling_method", "none");
  application->Options()->SetNumericValue("bound_relax_factor", 0.0);
  application->Options()->SetStringValue("honor_original_bounds", "yes");
  AttemptResult result;
  const auto initialized = application->Initialize();
  if (initialized != Ipopt::Solve_Succeeded) {
    result.optimizer_status = "initialization_" + ipopt_status_name(initialized);
    return result;
  }
  const auto status = application->OptimizeTNLP(problem);
  result.optimizer_status = ipopt_status_name(status);
  result.optimizer_converged = status == Ipopt::Solve_Succeeded;
  result.solution_available = raw_problem->solution_available();
  if (result.solution_available) result.variables = raw_problem->solution();
  result.callback_error = raw_problem->callback_error();
  return result;
}
bool inside_bounds(const std::array<double, 3> &variables,
                   const SaturationBounds &bounds) {
  for (std::size_t index = 0; index < variables.size(); ++index) {
    if (variables[index] < bounds.solver_lower[index] ||
        variables[index] > bounds.solver_upper[index]) {
      return false;
    }
  }
  return true;
}
double feasibility_measure(const NlpEvaluation &evaluation) {
  return std::max({std::abs(evaluation.constraints[0]), std::abs(evaluation.constraints[1]),
                   std::abs(evaluation.constraints[2])});
}
bool assess_candidate(const AttemptResult &attempt,
                      const NlpEvaluation &evaluation,
                      const SaturationBounds &bounds,
                      SaturationSolveResult &result) {
  const double vapor_density = 1.0 / evaluation.vapor.volume_m3;
  const double liquid_density = 1.0 / evaluation.liquid.volume_m3;
  const double saturation_pressure = evaluation.saturation_pressure_pa;
  const double vapor_pressure_difference =
      std::abs(evaluation.vapor.eos.pressure_pa - saturation_pressure);
  const double liquid_pressure_difference =
      std::abs(evaluation.liquid.eos.pressure_pa - saturation_pressure);
  result.vapor_pressure_relative_residual =
      vapor_pressure_difference / saturation_pressure;
  result.liquid_pressure_relative_residual =
      liquid_pressure_difference / saturation_pressure;
  result.pressure_relative_residual =
      std::max(result.vapor_pressure_relative_residual,
               result.liquid_pressure_relative_residual);
  result.chemical_potential_absolute_residual =
      std::abs(evaluation.constraints[2]);
  result.phase_density_distance =
      (liquid_density - vapor_density) / liquid_density;
  result.vapor_residual_density_derivative_m3_per_mol =
      evaluation.jacobian[0] / vapor_density;
  result.liquid_residual_density_derivative_m3_per_mol =
      evaluation.jacobian[4] / liquid_density;
  result.feasibility_measure = feasibility_measure(evaluation);
  result.feasibility_tolerance = kFeasibilityTolerance;
  const bool finite_values =
      std::isfinite(evaluation.saturation_pressure_pa) &&
      std::isfinite(vapor_density) && std::isfinite(liquid_density) &&
      std::isfinite(result.feasibility_measure) &&
      std::isfinite(result.pressure_relative_residual) &&
      std::isfinite(result.chemical_potential_absolute_residual) &&
      std::isfinite(result.phase_density_distance);
  const bool finite_stability =
      std::isfinite(result.vapor_residual_density_derivative_m3_per_mol) &&
      std::isfinite(result.liquid_residual_density_derivative_m3_per_mol);
  const bool locally_stable = finite_stability &&
      evaluation.jacobian[0] > 0.0 && evaluation.jacobian[4] > 0.0;
  const bool in_bounds = inside_bounds(attempt.variables, bounds);
  result.candidate_available = finite_values &&
      in_bounds;
  const bool physical[] = {
      finite_values,
      evaluation.saturation_pressure_pa > 0.0,
      vapor_density > 0.0,
      liquid_density > vapor_density,
      result.phase_density_distance > kPhaseDistanceTolerance,
      locally_stable,
      in_bounds,
  };
  result.physical_status_count = 0;
  for (const bool status : physical) result.physical_status_count += status;
  result.accepted = result.candidate_available &&
      result.feasibility_measure <= kFeasibilityTolerance &&
      result.physical_status_count == 7;
  result.physical_status = result.accepted ? "accepted" : "rejected";
  if (result.accepted) {
    result.rejection_status = "none";
  } else if (finite_values && liquid_density <= vapor_density) {
    result.rejection_status = "root-swap";
  } else if (finite_values && result.phase_density_distance <= kPhaseDistanceTolerance) {
    result.rejection_status = "coalescence";
  } else if (finite_values && (!locally_stable || !in_bounds)) {
    result.rejection_status = "branch-loss";
  } else {
    result.rejection_status = "solve-failure";
  }
  if (!result.accepted) {
    result.failure_reason = result.rejection_status == "solve-failure"
        ? (result.candidate_available ? "coexistence feasibility or local physical checks failed"
                                      : "candidate was not finite and in-domain")
        : result.rejection_status;
  } else {
    result.failure_reason.clear();
  }
  return result.candidate_available;
}
SaturationSolveResult solve_saturation(const NativeModel &model,
                                       double temperature_k,
                                       double liquid_density_upper_mol_m3,
                                       std::optional<std::array<double, 3>> explicit_start) {
  const std::array<double, 7> vapor_seeds{1.0e-5, 1.0e-3, 0.1,   10.0,
                                          100.0,  1000.0, 3000.0};
  const std::array<double, 3> liquid_seeds{18'000.0, 12'000.0, 24'000.0};
  const SaturationBounds bounds = saturation_bounds(liquid_density_upper_mol_m3);
  SaturationSolveResult result;
  result.bounds = bounds;
  result.failure_reason = "no finite in-domain saturation candidate";
  std::vector<std::array<double, 3>> starts;
  if (explicit_start.has_value()) {
    const auto &start = *explicit_start;
    for (double value : start) require_finite(value, "saturation start");
    if (start[0] <= 0.0 || start[1] <= 0.0 || start[2] <= 0.0)
      throw std::invalid_argument("saturation start values must be positive");
    starts.push_back({std::log(start[0]), std::log(start[1]), std::log(start[2])});
  }
  for (double vapor_density : vapor_seeds) {
    for (double liquid_density : liquid_seeds) {
      starts.push_back({std::log(vapor_density), std::log(liquid_density),
                        std::log(std::max(vapor_density * GAS_CONSTANT * temperature_k,
                                          kPressureLowerPa))});
    }
  }
  bool have_candidate = false;
  double best_feasibility = std::numeric_limits<double>::infinity();
  int best_physical_count = -1;
  std::size_t best_start = 0;
  bool explicit_certified = false;
  SaturationSolveResult last_failed_attempt = result;
  for (std::size_t start_index = 0; start_index < starts.size(); ++start_index) {
    if (explicit_start.has_value() && start_index > 0 && explicit_certified)
      break;
    const AttemptResult attempt = run_saturation_ipopt(
        model, temperature_k, starts[start_index], bounds);
    SaturationSolveResult candidate;
    candidate.optimizer_converged = attempt.optimizer_converged;
    candidate.optimizer_status = attempt.optimizer_status;
    candidate.bounds = bounds;
    candidate.feasibility_tolerance = kFeasibilityTolerance;
    candidate.failure_reason = attempt.callback_error.empty()
        ? "no finite in-domain saturation candidate"
        : attempt.callback_error;
    candidate.evaluation = {};
    bool available = false;
    try {
      if (attempt.solution_available && inside_bounds(attempt.variables, bounds)) {
        const auto evaluation = evaluate_saturation_state(
            model, temperature_k, attempt.variables);
        candidate.evaluation = evaluation;
        available = assess_candidate(attempt, evaluation,
                                     bounds, candidate);
      }
    } catch (const std::exception &error) {
      candidate.failure_reason = error.what();
    }
    if (!available) {
      last_failed_attempt = std::move(candidate);
      continue;
    }
    if (explicit_start.has_value() && start_index == 0)
      explicit_certified = candidate.accepted;
    const bool better = !have_candidate ||
        (candidate.feasibility_measure < best_feasibility) ||
        (candidate.feasibility_measure == best_feasibility &&
         (candidate.physical_status_count > best_physical_count ||
          (candidate.physical_status_count == best_physical_count &&
           start_index < best_start)));
    if (better) {
      have_candidate = true;
      best_feasibility = candidate.feasibility_measure;
      best_physical_count = candidate.physical_status_count;
      best_start = start_index;
      result = std::move(candidate);
    }
  }
  if (!have_candidate) result = std::move(last_failed_attempt);
  return result;
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
                       double liquid_density_upper_mol_m3,
                       std::optional<std::array<double, 3>> explicit_start) {
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
                          liquid_density_upper_mol_m3, explicit_start);
}
SaturationSensitivityResult solve_local_saturation_active(
    const NativeModel &model, double temperature_k,
    const std::vector<NativeActiveCoordinate> &active_parameters) {
  if (active_parameters.empty()) {
    throw std::invalid_argument(
        "active saturation parameter set must be nonempty");
  }
  SaturationSensitivityResult result;
  const NativeModel active_model = model.with_active_coordinates(active_parameters);
  const double liquid_density_upper_mol_m3 =
      active_model.molar_density_upper_bound(temperature_k, {1.0});
  result.solve = solve_local_saturation(active_model, temperature_k,
                                        liquid_density_upper_mol_m3);
  if (!result.solve.accepted) {
    result.sensitivity_failure_reason = "value_not_certified";
    return result;
  }
  const std::size_t sensitivity_dimension = 3;
  const std::size_t count = active_parameters.size();
  result.parameter_identities.reserve(count);
  for (const auto &parameter : active_parameters) {
    result.parameter_identities.push_back(parameter.identity);
  }
  auto phase_partials = [&](double volume) {
    std::vector<double> pressure(count);
    std::vector<double> chemical_potential(count);
    const auto block = active_model.evaluate_phase_derivatives(
        temperature_k, {1.0}, volume, active_parameters, 2);
    const std::size_t width = count + 2;
    for (std::size_t parameter = 0; parameter < count; ++parameter) {
      const std::size_t coordinate = parameter + 2;
      pressure[parameter] = -GAS_CONSTANT * temperature_k *
                            block.hessian[width + coordinate];
      chemical_potential[parameter] = block.hessian[coordinate];
    }
    return std::pair{std::move(pressure), std::move(chemical_potential)};
  };
  std::pair<std::vector<double>, std::vector<double>> vapor;
  std::pair<std::vector<double>, std::vector<double>> liquid;
  try {
    vapor = phase_partials(result.solve.evaluation.vapor.volume_m3);
    liquid = phase_partials(result.solve.evaluation.liquid.volume_m3);
  } catch (const std::exception &error) {
    result.sensitivity_status = "unavailable";
    result.sensitivity_failure_reason = error.what();
    return result;
  }
  const double saturation_pressure =
      result.solve.evaluation.saturation_pressure_pa;
  result.derivatives.resize(3 * count);
  const double vapor_density = 1.0 / result.solve.evaluation.vapor.volume_m3;
  const double liquid_density = 1.0 / result.solve.evaluation.liquid.volume_m3;
  std::vector<double> residual_state_jacobian(
      result.solve.evaluation.jacobian.begin(),
      result.solve.evaluation.jacobian.end());
  std::vector<double> residual_parameter_jacobian(3 * count);
  for (std::size_t parameter = 0; parameter < count; ++parameter) {
    residual_parameter_jacobian[parameter] =
        vapor.first[parameter] / saturation_pressure;
    residual_parameter_jacobian[count + parameter] =
        liquid.first[parameter] / saturation_pressure;
    residual_parameter_jacobian[2 * count + parameter] =
        vapor.second[parameter] - liquid.second[parameter];
  }
  const auto sensitivity =
      epcsaft_engine::interfaces::solve_rank_revealing_linear_system(
          residual_state_jacobian, residual_parameter_jacobian,
          sensitivity_dimension, count);
  result.sensitivity_rank = sensitivity.rank;
  result.sensitivity_dimension = sensitivity_dimension;
  if (std::isfinite(sensitivity.condition_number))
    result.sensitivity_condition_number = sensitivity.condition_number;
  if (std::isfinite(sensitivity.backward_error))
    result.sensitivity_backward_error = sensitivity.backward_error;
  const bool sensitivity_available = sensitivity.solution.size() ==
      sensitivity_dimension * count;
  result.sensitivity_status = sensitivity_available ? "available" : "unavailable";
  if (!sensitivity_available) {
    result.sensitivity_failure_reason = sensitivity.rank !=
            sensitivity_dimension
        ? "singular_residual_jacobian"
        : sensitivity.condition_number >
                epcsaft_engine::interfaces::
                    kSolvedStateSensitivityConditionNumberMax
            ? "ill_conditioned_residual_jacobian"
            : "linear_solve_backward_error";
    result.derivatives.clear();
    return result;
  }
  result.sensitivity_failure_reason.clear();
  for (std::size_t parameter = 0; parameter < count; ++parameter) {
    result.derivatives[parameter] = saturation_pressure *
        sensitivity.solution[2 * count + parameter];
    result.derivatives[count + parameter] = vapor_density *
        sensitivity.solution[parameter];
    result.derivatives[2 * count + parameter] = liquid_density *
        sensitivity.solution[count + parameter];
  }
  return result;
}
} // namespace epcsaft
