#include "pure_saturation.hpp"

#include <array>
#include <utility>
#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

namespace {

py::dict phase_to_dict(const epcsaft::SaturationPhaseEvaluation &phase) {
  py::dict result;
  result["amount_mol"] = phase.amount_mol;
  result["volume_m3"] = phase.volume_m3;
  result["molar_density_mol_m3"] = phase.amount_mol / phase.volume_m3;
  result["pressure_pa"] = phase.eos.pressure_pa;
  result["chemical_potential_over_rt"] = phase.eos.chemical_potential_over_rt;
  return result;
}

py::list attempt_log_to_list(
    const std::vector<epcsaft::SaturationSolveResult::AttemptRecord>
        &attempts) {
  py::list result;
  for (const auto &attempt : attempts) {
    py::dict item;
    item["role"] = attempt.role;
    item["initial_guess"] = attempt.initial_guess;
    item["solver_converged"] = attempt.solver_converged;
    item["solver_status"] = attempt.solver_status;
    item["iterations"] = attempt.iterations;
    item["constraint_violation"] = attempt.constraint_violation;
    item["callback_error"] = attempt.callback_error;
    result.append(std::move(item));
  }
  return result;
}

py::dict evaluate_saturation_nlp(const epcsaft::NativeModel &model,
                                 double temperature_k,
                                 const std::array<double, 3> &variables,
                                 const std::array<double, 3> &multipliers) {
  const auto evaluation = epcsaft::evaluate_saturation_nlp(
      model, temperature_k, variables, multipliers);
  py::dict result;
  result["objective"] = evaluation.objective;
  result["objective_gradient"] = evaluation.objective_gradient;
  result["constraints"] = evaluation.constraints;
  result["jacobian"] = evaluation.jacobian;
  result["lagrangian_hessian_lower"] = evaluation.lagrangian_hessian_lower;
  return result;
}

py::dict saturation_to_dict(const epcsaft::SaturationSolveResult &solve,
                            double temperature_k,
                            bool exact_derivatives) {
  py::dict diagnostics;
  diagnostics["solver_converged"] = solve.solver_converged;
  diagnostics["solver_status"] = solve.solver_status;
  diagnostics["iterations"] = solve.iterations;
  diagnostics["attempts"] = solve.attempts;
  diagnostics["attempt_log"] = attempt_log_to_list(solve.attempt_log);
  diagnostics["solver_lower_bounds"] = solve.solver_lower_bounds;
  diagnostics["solver_upper_bounds"] = solve.solver_upper_bounds;
  diagnostics["solver_constraint_violation"] =
      solve.solver_constraint_violation;
  diagnostics["numerical_converged"] = solve.numerical_converged;
  diagnostics["confirmation_solves"] = solve.confirmation_solves;
  diagnostics["confirmation_max_relative_difference"] =
      solve.confirmation_max_relative_difference;
  diagnostics["physical_accepted"] = solve.physical_accepted;
  diagnostics["pressure_relative_residual"] = solve.pressure_relative_residual;
  diagnostics["vapor_pressure_relative_residual"] =
      solve.vapor_pressure_relative_residual;
  diagnostics["liquid_pressure_relative_residual"] =
      solve.liquid_pressure_relative_residual;
  diagnostics["chemical_potential_absolute_residual"] =
      solve.chemical_potential_absolute_residual;
  diagnostics["phase_density_distance"] = solve.phase_density_distance;
  diagnostics["vapor_residual_density_derivative_m3_per_mol"] =
      solve.vapor_residual_density_derivative_m3_per_mol;
  diagnostics["liquid_residual_density_derivative_m3_per_mol"] =
      solve.liquid_residual_density_derivative_m3_per_mol;
  diagnostics["rejection_status"] = solve.rejection_status;
  diagnostics["exact_derivatives"] = exact_derivatives;
  diagnostics["globality_certificate"] = false;
  diagnostics["failure_reason"] = solve.failure_reason;

  py::dict result;
  result["accepted"] = solve.accepted;
  result["temperature_k"] = temperature_k;
  result["diagnostics"] = diagnostics;
  if (solve.accepted) {
    result["saturation_pressure_pa"] = solve.evaluation.saturation_pressure_pa;
    result["vapor"] = phase_to_dict(solve.evaluation.vapor);
    result["liquid"] = phase_to_dict(solve.evaluation.liquid);
  }
  return result;
}

py::dict solve_saturation(const epcsaft::NativeModel &model,
                          double temperature_k,
                          double liquid_density_upper_mol_m3) {
  epcsaft::SaturationSolveResult solve;
  {
    py::gil_scoped_release release;
    solve = epcsaft::solve_local_saturation(model, temperature_k,
                                            liquid_density_upper_mol_m3);
  }
  return saturation_to_dict(solve, temperature_k, solve.accepted);
}

py::dict solve_saturation_active(
    const epcsaft::NativeModel &model, double temperature_k,
    const std::vector<epcsaft::NativeActiveCoordinate> &active_parameters) {
  epcsaft::SaturationSensitivityResult sensitivity;
  {
    py::gil_scoped_release release;
    sensitivity = epcsaft::solve_local_saturation_active(
        model, temperature_k, active_parameters);
  }
  const bool exact_derivatives =
      sensitivity.solve.accepted &&
      sensitivity.parameter_identities.size() == active_parameters.size() &&
      sensitivity.derivatives.size() == 3 * active_parameters.size();
  py::dict result =
      saturation_to_dict(sensitivity.solve, temperature_k, exact_derivatives);
  result["parameter_identities"] = sensitivity.parameter_identities;
  result["derivatives"] = sensitivity.derivatives;
  return result;
}

} // namespace

void bind_pure_saturation(py::module_ &module) {
  module.def("_evaluate_saturation_nlp", &evaluate_saturation_nlp,
             py::arg("model"), py::arg("temperature_k"), py::arg("variables"),
             py::arg("multipliers"));
  module.def("_solve_saturation", &solve_saturation, py::arg("model"),
             py::arg("temperature_k"), py::arg("liquid_density_upper_mol_m3"));
  module.def("_solve_saturation_active", &solve_saturation_active,
             py::arg("model"), py::arg("temperature_k"),
             py::arg("active_parameters"));
}
