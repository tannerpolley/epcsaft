#include "pure_saturation.hpp"

#include <array>
#include <limits>
#include <optional>
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

py::dict evaluate_saturation_nlp(const epcsaft::NativeModel &model,
                                 double temperature_k,
                                 const std::array<double, 3> &variables,
                                 const std::array<double, 3> &multipliers) {
  const auto evaluation = epcsaft::evaluate_saturation_nlp(
      model, temperature_k, variables, multipliers);
  py::dict result;
  result["constraints"] = evaluation.constraints;
  result["jacobian"] = evaluation.jacobian;
  result["lagrangian_hessian_lower"] = evaluation.lagrangian_hessian_lower;
  return result;
}

py::dict saturation_to_dict(const epcsaft::SaturationSolveResult &solve,
                            double temperature_k, bool exact_derivatives,
                            const std::string &sensitivity_status = "not_requested",
                            const std::string &sensitivity_failure_reason = "not_requested",
                            std::size_t sensitivity_rank = 0,
                            std::size_t sensitivity_dimension = 0,
                            double sensitivity_condition_number =
                                std::numeric_limits<double>::infinity(),
                            double sensitivity_backward_error =
                                std::numeric_limits<double>::infinity()) {
  py::dict diagnostics;
  diagnostics["candidate_available"] = solve.candidate_available;
  diagnostics["optimizer_converged"] = solve.optimizer_converged;
  diagnostics["optimizer_status"] = solve.optimizer_status;
  diagnostics["feasibility_measure"] = solve.feasibility_measure;
  diagnostics["feasibility_tolerance"] = solve.feasibility_tolerance;
  diagnostics["physical_accepted"] = solve.physical_accepted;
  diagnostics["physical_status"] = solve.physical_status;
  diagnostics["physical_status_count"] = solve.physical_status_count;
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
  diagnostics["sensitivity_status"] = sensitivity_status;
  diagnostics["sensitivity_failure_reason"] = sensitivity_failure_reason;
  diagnostics["sensitivity_rank"] = sensitivity_rank;
  diagnostics["sensitivity_dimension"] = sensitivity_dimension;
  diagnostics["sensitivity_condition_number"] = sensitivity_condition_number;
  diagnostics["sensitivity_backward_error"] = sensitivity_backward_error;
  diagnostics["globality_certificate"] = false;
  diagnostics["failure_reason"] = solve.failure_reason;

  py::dict result;
  result["accepted"] = solve.accepted;
  result["temperature_k"] = temperature_k;
  result["diagnostics"] = diagnostics;
  if (solve.candidate_available) {
    result["saturation_pressure_pa"] = solve.evaluation.saturation_pressure_pa;
    result["vapor"] = phase_to_dict(solve.evaluation.vapor);
    result["liquid"] = phase_to_dict(solve.evaluation.liquid);
  }
  return result;
}

py::dict solve_saturation(const epcsaft::NativeModel &model,
                          double temperature_k,
                          double liquid_density_upper_mol_m3,
                          std::optional<std::array<double, 3>> explicit_start) {
  epcsaft::SaturationSolveResult solve;
  {
    py::gil_scoped_release release;
    solve = epcsaft::solve_local_saturation(model, temperature_k,
                                            liquid_density_upper_mol_m3,
                                            explicit_start);
  }
  return saturation_to_dict(solve, temperature_k, false);
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
      sensitivity.sensitivity_status == "available" &&
      sensitivity.parameter_identities.size() == active_parameters.size() &&
      sensitivity.derivatives.size() == 3 * active_parameters.size();
  py::dict result =
      saturation_to_dict(sensitivity.solve, temperature_k, exact_derivatives,
                         sensitivity.sensitivity_status,
                         sensitivity.sensitivity_failure_reason,
                         sensitivity.sensitivity_rank,
                         sensitivity.sensitivity_dimension,
                         sensitivity.sensitivity_condition_number,
                         sensitivity.sensitivity_backward_error);
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
             py::arg("temperature_k"), py::arg("liquid_density_upper_mol_m3"),
             py::arg("explicit_start") = std::nullopt);
  module.def("_solve_saturation_active", &solve_saturation_active,
             py::arg("model"), py::arg("temperature_k"),
             py::arg("active_parameters"));
}
