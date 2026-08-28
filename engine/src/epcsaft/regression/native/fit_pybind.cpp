#include "fit_core.hpp"
#include "fit_pybind.hpp"
#include "grepe_pybind.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace py = pybind11;
using namespace pybind11::literals;

namespace epcsaft::regression {
namespace {

template <typename Value>
Value required(const py::dict &mapping, const char *name) {
  if (!mapping.contains(name)) {
    throw std::invalid_argument(std::string("missing field: ") + name);
  }
  return py::cast<Value>(mapping[name]);
}

ResidualKind residual_kind(const std::string &value) {
  if (value == "log_ratio") {
    return ResidualKind::LogRatio;
  }
  if (value == "scaled_difference") {
    return ResidualKind::ScaledDifference;
  }
  if (value == "left_censored") {
    return ResidualKind::LeftCensored;
  }
  if (value == "right_censored") {
    return ResidualKind::RightCensored;
  }
  throw std::invalid_argument("unknown residual kind: " + value);
}

std::vector<Coordinate> parse_coordinates(const py::iterable &values) {
  std::vector<Coordinate> result;
  for (const py::handle item : values) {
    const auto value = py::cast<py::dict>(item);
    result.push_back({required<double>(value, "origin"),
                      required<double>(value, "scale"),
                      required<double>(value, "lower_affine"),
                      required<double>(value, "upper_affine")});
  }
  return result;
}

std::vector<Observation> parse_observations(const py::iterable &values) {
  std::vector<Observation> result;
  for (const py::handle item : values) {
    const auto value = py::cast<py::dict>(item);
    Observation observation{required<double>(value, "multiplier"), {}};
    for (const py::handle target_item :
         required<py::iterable>(value, "targets")) {
      const auto target = py::cast<py::dict>(target_item);
      observation.targets.push_back(
          {required<double>(target, "observed"),
           residual_kind(required<std::string>(target, "residual")),
           required<double>(target, "scale"),
           target.contains("multiplier")
               ? py::cast<double>(target["multiplier"])
               : 1.0});
    }
    result.push_back(std::move(observation));
  }
  return result;
}

Controls parse_controls(const py::dict &value) {
  std::vector<Start> starts;
  for (const py::handle item : required<py::iterable>(value, "starts")) {
    const auto start = py::cast<py::dict>(item);
    starts.push_back({required<std::string>(start, "identity"),
                      required<std::vector<double>>(start, "affine")});
  }
  return {std::move(starts), required<double>(value, "affine_bound_margin"),
          required<int>(value, "maximum_iterations"),
          required<double>(value, "optimizer_tolerance")};
}

epcsaft_engine::equilibrium::GrepeIntensive
parse_intensive(const py::dict &value) {
  epcsaft_engine::equilibrium::GrepeIntensive result;
  result.solved = required<bool>(value, "solved");
  result.value = required<double>(value, "value");
  if (result.solved) {
    result.bounds = required<std::array<double, 2>>(value, "bounds");
    result.starts = required<std::vector<double>>(value, "starts");
  }
  return result;
}

std::vector<Experiment> parse_experiments(const py::iterable &values) {
  std::vector<Experiment> result;
  for (const py::handle item : values) {
    const auto value = py::cast<py::dict>(item);
    const std::string kind = required<std::string>(value, "kind");
    if (kind == "direct_eos") {
      const std::string descriptor =
          required<std::string>(value, "descriptor_identity");
      const auto state = required<py::dict>(value, "state");
      result.emplace_back(DirectEosExperiment{
          required<NativeModel>(value, "model"),
          descriptor == "eos.fixed-state" ? DirectEosKind::FixedState
                                            : DirectEosKind::PureSaturation,
          required<double>(state, "temperature_k"),
          state.contains("molar_density_mol_per_m3")
              ? py::cast<double>(state["molar_density_mol_per_m3"])
              : 0.0,
          value.contains("composition") && !value["composition"].is_none()
              ? py::cast<std::vector<double>>(value["composition"])
              : std::vector<double>{1.0},
          required<std::vector<std::size_t>>(value, "selected_rows"),
      });
      if (descriptor != "eos.fixed-state" &&
          descriptor != "eos.pure-saturation") {
        throw std::invalid_argument("unknown direct EOS descriptor identity");
      }
      continue;
    }
    if (kind == "grepe") {
      const std::string fingerprint =
          required<std::string>(value, "source_reference_fingerprint");
      const auto temperature =
          parse_intensive(required<py::dict>(value, "temperature"));
      const auto pressure =
          parse_intensive(required<py::dict>(value, "pressure"));
      result.emplace_back(GrepeExperiment{
          required<std::string>(value, "structure_fingerprint"),
          required<NativeModel>(value, "model"),
          required<NativeModel>(value, "incipient_model"),
          epcsaft_engine::equilibrium::parse_grepe_request(
              value["problem"], true, fingerprint, &temperature, &pressure),
          required<std::vector<std::size_t>>(value, "selected_rows"),
      });
      continue;
    }
    throw std::invalid_argument("unknown native regression experiment kind");
  }
  return result;
}

py::dict diagnostics_mapping(const Diagnostics &value) {
  py::dict result;
  result["singular_values"] = value.singular_values;
  result["rank"] = value.rank;
  result["free_subspace_rank"] = value.free_subspace_rank;
  result["condition_number"] = value.condition_number;
  result["affine_covariance"] = value.affine_covariance;
  result["degrees_of_freedom"] = value.degrees_of_freedom;
  result["information_row_count"] = value.information_row_count;
  result["free_coordinate_indices"] = value.free_coordinate_indices;
  result["residual_variance"] = value.residual_variance;
  result["covariance_status"] = value.covariance_status;
  return result;
}

py::dict start_mapping(const StartResult &value,
                       const std::vector<Coordinate> &coordinates) {
  py::dict result;
  std::vector<double> initial_physical;
  std::vector<double> final_physical;
  std::vector<double> physical_jacobian = value.affine_jacobian;
  for (std::size_t coordinate = 0; coordinate < coordinates.size();
       ++coordinate) {
    initial_physical.push_back(coordinates[coordinate].origin +
                               value.initial_affine[coordinate] *
                                   coordinates[coordinate].scale);
    final_physical.push_back(coordinates[coordinate].origin +
                             value.final_affine[coordinate] *
                                 coordinates[coordinate].scale);
    for (std::size_t row = 0; row < value.residuals.size(); ++row) {
      physical_jacobian[row * coordinates.size() + coordinate] /=
          coordinates[coordinate].scale;
    }
  }
  result["index"] = value.index;
  result["start_identity"] = value.start_identity;
  result["initial_affine"] = value.initial_affine;
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
  result["trial_failure_code"] = value.trial_failure_code;
  result["iterations"] = value.iterations;
  result["non_evaluable_count"] = value.non_evaluable_count;
  py::list parameter_constraints;
  for (const auto &constraint : value.parameter_domain_constraints) {
    parameter_constraints.append(py::dict(
        "identity"_a = constraint.identity,
        "temperature_k"_a = constraint.temperature_k,
        "initial_value"_a = constraint.initial_value,
        "final_value"_a = constraint.final_value,
        "upper"_a = constraint.upper,
        "affine_jacobian"_a = constraint.affine_jacobian));
  }
  result["parameter_domain_constraints"] = parameter_constraints;
  py::list failures;
  for (const auto &failure : value.owner_failures) {
    failures.append(py::dict("code"_a = failure.code, "reason"_a = failure.reason,
                             "count"_a = failure.count));
  }
  result["owner_failures"] = failures;
  result["owner_receipts_json"] = py::list();
  py::list outcomes;
  for (const auto &outcome : value.batch_outcomes) {
    outcomes.append(py::dict("index"_a = outcome.index,
                             "row_count"_a = outcome.row_count,
                             "status"_a = outcome.status,
                             "failure_code"_a = outcome.failure_code));
  }
  result["batch_outcomes"] = outcomes;
  py::list owner_results;
  for (const auto &owner : value.owner_results) {
    owner_results.append(owner.has_value() ? py::cast(*owner) : py::none());
  }
  result["owner_results"] = owner_results;
  result["cost"] = value.cost;
  result["usable"] = value.usable;
  result["optimizer_basin"] =
      value.usable ? py::cast(value.optimizer_basin) : py::none();
  return result;
}

py::dict fit_resolved(
    std::size_t active_parameter_count,
    const std::vector<std::size_t> &identity_coordinates,
    const py::iterable &coordinate_values,
    const py::iterable &observation_values, const py::dict &control_values,
    const std::vector<std::string> &failure_codes,
    const py::iterable &experiment_values,
    const std::string &evaluation_strategy) {
  FitRequest request{active_parameter_count,
                     identity_coordinates,
                     parse_coordinates(coordinate_values),
                     parse_observations(observation_values),
                     parse_controls(control_values),
                     failure_codes,
                     parse_experiments(experiment_values),
                     evaluation_strategy};
  FitResult native;
  {
    py::gil_scoped_release release;
    native = fit_native(request);
  }
  py::dict result;
  result["evaluation_strategy"] = native.evaluation_strategy;
  if (native.preflight_failure.has_value()) {
    const auto &failure = *native.preflight_failure;
    result["preflight_failure"] =
        py::dict("code"_a = failure.code,
                 "start_identity"_a = failure.start_identity,
                 "batch_identity"_a = failure.batch_identity,
                 "diagnostic"_a = failure.diagnostic);
  } else {
    result["preflight_failure"] = py::none();
  }
  py::list starts;
  for (const auto &start : native.starts) {
    starts.append(start_mapping(start, request.coordinates));
  }
  result["starts"] = starts;
  return result;
}

} // namespace
} // namespace epcsaft::regression

void bind_parameter_fit(py::module_ &module) {
  module.def("_fit_resolved", &epcsaft::regression::fit_resolved,
             py::arg("active_parameter_count"),
             py::arg("identity_coordinates"), py::arg("coordinates"),
             py::arg("observations"), py::arg("controls"),
             py::arg("failure_codes"), py::arg("experiments"),
             py::arg("evaluation_strategy") = "automatic");
}
