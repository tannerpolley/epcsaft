#include "ceres_least_squares.hpp"
#include "eos_regression_evaluator.hpp"
#include "grepe_pybind.hpp"
#include "../../../cpp/equilibrium/owner_evaluator.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <cstddef>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace py = pybind11;
using namespace pybind11::literals;

namespace {

template <typename Value>
Value required(const py::dict& mapping, const char* name) {
    if (!mapping.contains(name)) {
        throw std::invalid_argument(std::string("missing field: ") + name);
    }
    return py::cast<Value>(mapping[name]);
}

std::string optional_text(const py::dict& mapping, const char* name) {
    return mapping.contains(name) && !mapping[name].is_none()
        ? py::cast<std::string>(mapping[name]) : std::string{};
}

struct CompiledOwner final {
    explicit CompiledOwner(
        std::shared_ptr<const epcsaft_engine::interfaces::OwnerEvaluator> owner,
        py::dict sidecar
    ) : owner(std::move(owner)), sidecar(std::move(sidecar)) {}

    std::shared_ptr<const epcsaft_engine::interfaces::OwnerEvaluator> owner;
    py::dict sidecar;
};

std::vector<epcsaft::NativeActiveCoordinate> active_coordinates(
    const py::iterable& values
) {
    std::vector<epcsaft::NativeActiveCoordinate> result;
    for (const py::handle item : values) {
        result.push_back(py::cast<epcsaft::NativeActiveCoordinate>(item));
    }
    return result;
}

epcsaft::EosRegressionSpec parse_eos_spec(
    const py::dict& value,
    const std::vector<epcsaft::NativeActiveCoordinate>& active
) {
    const std::string family = required<std::string>(value, "family");
    epcsaft::EosRegressionSpec result;
    if (family == "fixed-state") {
        result.family = epcsaft::EosRegressionFamily::FixedState;
    } else if (family == "pure-saturation") {
        result.family = epcsaft::EosRegressionFamily::PureSaturation;
    } else {
        throw std::invalid_argument("unsupported EOS evaluator family");
    }
    result.temperature_k = required<double>(value, "temperature_k");
    result.row_identities = required<std::vector<std::string>>(value, "row_identities");
    result.active_parameter_identities = required<std::vector<std::string>>(
        value, "active_parameter_identities");
    result.active_parameter_units = required<std::vector<std::string>>(
        value, "active_parameter_units");
    result.active_parameter_bases = required<std::vector<std::string>>(
        value, "active_parameter_bases");
    result.source_spec_fingerprint = required<std::string>(
        value, "source_spec_fingerprint");
    result.evaluator_fingerprint = required<std::string>(
        value, "evaluator_fingerprint");
    if (result.family == epcsaft::EosRegressionFamily::FixedState) {
        result.molar_density_mol_m3 = required<double>(value, "molar_density_mol_m3");
        result.mole_fractions = required<std::vector<double>>(value, "mole_fractions");
    }
    if (result.active_parameter_identities.size() != active.size()) {
        throw std::invalid_argument("active parameter metadata does not match coordinates");
    }
    return result;
}

py::dict eos_sidecar(const epcsaft::EosRegressionEvaluatorSidecar& value) {
    py::dict result;
    result["owner_kind"] = value.owner_kind;
    result["parameter_fingerprint"] = value.parameter_fingerprint;
    result["parameter_contract_fingerprint"] = value.parameter_contract_fingerprint;
    result["source_spec_fingerprint"] = value.source_spec_fingerprint;
    result["evaluator_fingerprint"] = value.evaluator_fingerprint;
    result["active_parameter_identities"] = value.active_parameter_identities;
    result["active_parameter_values"] = value.active_parameter_values;
    result["active_parameter_units"] = value.active_parameter_units;
    result["active_parameter_bases"] = value.active_parameter_bases;
    py::list rows;
    for (const auto& row : value.output_rows) {
        const char* kind = "unknown";
        switch (row.kind) {
        case epcsaft::EosRegressionRowKind::FixedHelmholtz: kind = "fixed_helmholtz"; break;
        case epcsaft::EosRegressionRowKind::FixedPressure: kind = "fixed_pressure"; break;
        case epcsaft::EosRegressionRowKind::FixedChemicalPotential: kind = "fixed_chemical_potential"; break;
        case epcsaft::EosRegressionRowKind::SaturationPressure: kind = "saturation_pressure"; break;
        case epcsaft::EosRegressionRowKind::VaporMolarDensity: kind = "vapor_molar_density"; break;
        case epcsaft::EosRegressionRowKind::LiquidMolarDensity: kind = "liquid_molar_density"; break;
        case epcsaft::EosRegressionRowKind::VaporMassDensity: kind = "vapor_mass_density"; break;
        case epcsaft::EosRegressionRowKind::LiquidMassDensity: kind = "liquid_mass_density"; break;
        }
        rows.append(py::dict("identity"_a = row.identity, "unit"_a = row.unit,
                             "basis"_a = row.basis, "kind"_a = kind));
    }
    result["output_rows"] = rows;
    result["derivative_provenance"] = value.derivative_provenance;
    result["acceptance_policy"] = value.acceptance_policy;
    result["topology_fingerprint"] = value.topology_fingerprint;
    result["component_order_fingerprint"] = value.component_order_fingerprint;
    result["domain_fingerprint"] = value.domain_fingerprint;
    result["artifact_fingerprint"] = value.artifact_fingerprint;
    result["candidate_qualification"] = value.candidate_qualification;
    result["candidate_source_fingerprint"] = value.candidate_source_fingerprint;
    result["temperature_domain_k"] = py::make_tuple(value.temperature_min_k, value.temperature_max_k);
    result["pressure_domain_pa"] = py::make_tuple(value.pressure_min_pa, value.pressure_max_pa);
    result["ion_mole_fraction_max"] = value.ion_mole_fraction_max;
    return result;
}

epcsaft_engine::equilibrium::EquilibriumEvaluatorSetup parse_equilibrium_setup(
    const py::dict& value
) {
    epcsaft_engine::equilibrium::EquilibriumEvaluatorSetup result;
    result.owner_kind = optional_text(value, "owner_kind");
    result.parameter_baseline_digest = required<std::string>(value, "parameter_baseline_digest");
    result.source_spec_digest = required<std::string>(value, "source_spec_digest");
    result.anchor_identity = required<std::string>(value, "anchor_identity");
    result.anchor_reference_fingerprint = required<std::string>(value, "anchor_reference_fingerprint");
    result.derivative_provenance = optional_text(value, "derivative_provenance");
    result.topology_policy = optional_text(value, "topology_policy");
    result.branch_policy = optional_text(value, "branch_policy");
    result.conditioning_policy = optional_text(value, "conditioning_policy");
    for (const py::handle item : required<py::iterable>(value, "active_parameters")) {
        const py::dict parameter = py::cast<py::dict>(item);
        result.active_parameters.push_back({
            required<std::string>(parameter, "identity"),
            required<std::string>(parameter, "unit"),
            required<std::string>(parameter, "basis")});
    }
    for (const py::handle item : required<py::iterable>(value, "output_rows")) {
        const py::dict row = py::cast<py::dict>(item);
        result.output_rows.push_back({
            required<std::string>(row, "identity"),
            required<std::string>(row, "unit"),
            required<std::string>(row, "basis")});
    }
    return result;
}

py::dict equilibrium_sidecar(
    const epcsaft_engine::equilibrium::EquilibriumEvaluatorSidecar& value
) {
    py::dict result;
    result["owner_kind"] = value.owner_kind;
    result["parameter_baseline_digest"] = value.parameter_baseline_digest;
    result["source_spec_digest"] = value.source_spec_digest;
    result["anchor_identity"] = value.anchor_identity;
    result["anchor_reference_fingerprint"] = value.anchor_reference_fingerprint;
    py::list parameters;
    for (const auto& parameter : value.active_parameters) {
        parameters.append(py::dict("identity"_a = parameter.identity,
                                   "unit"_a = parameter.unit,
                                   "basis"_a = parameter.basis));
    }
    result["active_parameters"] = parameters;
    py::list rows;
    for (const auto& row : value.output_rows) {
        rows.append(py::dict("identity"_a = row.identity, "unit"_a = row.unit,
                             "basis"_a = row.basis));
    }
    result["output_rows"] = rows;
    result["derivative_provenance"] = value.derivative_provenance;
    py::list phases;
    for (const auto& phase : value.phases) {
        phases.append(py::dict("identity"_a = phase.identity, "role"_a = phase.role,
                               "support"_a = phase.support,
                               "support_indices"_a = phase.support_indices,
                               "finite"_a = phase.finite, "ideal_gas"_a = phase.ideal_gas));
    }
    result["phases"] = phases;
    result["temperature_role"] = value.temperature_role;
    result["pressure_role"] = value.pressure_role;
    result["reaction_identities"] = value.reaction_identities;
    result["reaction_phase_identities"] = value.reaction_phase_identities;
    result["topology_policy"] = value.topology_policy;
    result["branch_policy"] = value.branch_policy;
    result["conditioning_policy"] = value.conditioning_policy;
    result["parameter_contract_fingerprint"] = value.parameter_contract_fingerprint;
    result["evaluator_fingerprint"] = value.evaluator_fingerprint;
    return result;
}

py::dict options_mapping(const epcsaft_engine::regression::CeresOptionsRecord& value) {
    py::dict result;
    result["minimizer_type"] = value.minimizer_type;
    result["trust_region_strategy_type"] = value.trust_region_strategy_type;
    result["linear_solver_type"] = value.linear_solver_type;
    result["dense_linear_algebra_library_type"] = value.dense_linear_algebra_library_type;
    result["jacobi_scaling"] = value.jacobi_scaling;
    result["use_inner_iterations"] = value.use_inner_iterations;
    result["use_nonmonotonic_steps"] = value.use_nonmonotonic_steps;
    result["maximum_line_search_step_size_iterations"] = value.maximum_line_search_step_size_iterations;
    result["maximum_consecutive_invalid_steps"] = value.maximum_consecutive_invalid_steps;
    result["number_of_threads"] = value.number_of_threads;
    result["initial_trust_region_radius"] = value.initial_trust_region_radius;
    result["maximum_trust_region_radius"] = value.maximum_trust_region_radius;
    result["minimum_trust_region_radius"] = value.minimum_trust_region_radius;
    result["minimum_relative_decrease"] = value.minimum_relative_decrease;
    result["minimum_lm_diagonal"] = value.minimum_lm_diagonal;
    result["maximum_lm_diagonal"] = value.maximum_lm_diagonal;
    result["use_mixed_precision_solves"] = value.use_mixed_precision_solves;
    result["maximum_refinement_iterations"] = value.maximum_refinement_iterations;
    result["function_tolerance"] = value.function_tolerance;
    result["gradient_tolerance"] = value.gradient_tolerance;
    result["parameter_tolerance"] = value.parameter_tolerance;
    result["maximum_iterations"] = value.maximum_iterations;
    result["maximum_elapsed_time_seconds"] = value.maximum_elapsed_time_seconds;
    return result;
}

py::dict solve_mapping(const epcsaft_engine::regression::LeastSquaresResult& value) {
    py::dict result;
    result["status"] = value.status;
    result["termination"] = value.termination;
    result["message"] = value.message;
    result["usable"] = value.usable;
    result["parameter_contract_fingerprint"] = value.parameter_contract_fingerprint;
    result["residual_scale_metadata"] = value.residual_scale_metadata;
    py::list rows;
    for (const auto& row : value.evaluator_rows) {
        rows.append(py::dict("evaluator_fingerprint"_a = row.evaluator_fingerprint,
                             "evaluator_index"_a = row.evaluator_index,
                             "row_begin"_a = row.row_begin, "row_end"_a = row.row_end));
    }
    result["evaluator_rows"] = rows;
    result["parameters"] = value.parameters;
    result["predictions"] = value.predictions;
    result["residuals"] = value.residuals;
    result["physical_jacobian"] = value.physical_jacobian;
    result["optimizer_jacobian"] = value.optimizer_jacobian;
    py::list failures;
    for (const auto& failure : value.failures) {
        failures.append(py::dict(
            "evaluator_index"_a = failure.evaluator_index,
            "row_begin"_a = failure.row_begin, "row_end"_a = failure.row_end,
            "evaluation_phase"_a = failure.evaluation_phase,
            "physical_parameters"_a = failure.physical_parameters,
            "code"_a = failure.code, "reason"_a = failure.reason,
            "count"_a = failure.count));
    }
    result["failures"] = failures;
    result["requested_options"] = options_mapping(value.requested_options);
    result["resolved_options"] = options_mapping(value.resolved_options);
    const auto& counters = value.counters;
    result["counters"] = py::dict(
        "owner_call_limit"_a = counters.owner_call_limit,
        "uncached_residual_owner_attempts"_a = counters.uncached_residual_owner_attempts,
        "uncached_jacobian_owner_attempts"_a = counters.uncached_jacobian_owner_attempts,
        "owner_cache_hits"_a = counters.owner_cache_hits,
        "typed_trial_failures"_a = counters.typed_trial_failures,
        "ceres_residual_evaluations"_a = counters.ceres_residual_evaluations,
        "ceres_jacobian_evaluations"_a = counters.ceres_jacobian_evaluations,
        "accepted_steps"_a = counters.accepted_steps,
        "rejected_steps"_a = counters.rejected_steps,
        "iterations"_a = counters.iterations,
        "iteration_limit_reached"_a = counters.iteration_limit_reached,
        "elapsed_time_limit_reached"_a = counters.elapsed_time_limit_reached,
        "owner_call_limit_reached"_a = counters.owner_call_limit_reached);
    const auto& covariance = value.covariance;
    result["covariance"] = py::dict(
        "status"_a = covariance.status, "algorithm"_a = covariance.algorithm,
        "null_space_rank"_a = covariance.null_space_rank,
        "apply_loss_function"_a = covariance.apply_loss_function,
        "number_of_threads"_a = covariance.number_of_threads,
        "minimum_reciprocal_condition_number"_a = covariance.minimum_reciprocal_condition_number,
        "optimizer_coordinates"_a = covariance.optimizer_coordinates,
        "physical_coordinates"_a = covariance.physical_coordinates,
        "active_bound"_a = covariance.active_bound);
    return result;
}

std::vector<std::shared_ptr<const epcsaft_engine::interfaces::OwnerEvaluator>> owners(
    const py::iterable& values
) {
    std::vector<std::shared_ptr<const epcsaft_engine::interfaces::OwnerEvaluator>> result;
    for (const py::handle item : values) {
        result.push_back(py::cast<std::shared_ptr<CompiledOwner>>(item)->owner);
    }
    return result;
}

std::shared_ptr<CompiledOwner> compile_eos(
    const epcsaft::NativeModel& model, const py::dict& specification,
    const py::iterable& active_values
) {
    auto active = active_coordinates(active_values);
    auto compiled = epcsaft::compile_eos_regression_evaluator(
        model, parse_eos_spec(specification, active), active);
    auto owner = std::dynamic_pointer_cast<const epcsaft::EosOwnerEvaluator>(compiled.evaluator);
    if (!owner) {
        throw std::runtime_error("EOS owner evaluator type changed");
    }
    return std::make_shared<CompiledOwner>(
        compiled.evaluator, eos_sidecar(owner->sidecar()));
}

std::shared_ptr<CompiledOwner> compile_equilibrium(
    const epcsaft::NativeModel& model, const epcsaft::NativeModel& incipient_model,
    const py::handle& problem, const std::string& source_reference_fingerprint,
    const py::dict& setup_mapping
) {
    auto request = epcsaft_engine::equilibrium::parse_grepe_request(
        problem, true, source_reference_fingerprint);
    // Regression owns a fixed local branch: compilation is never allowed to
    // generate starts or discover a topology.
    request.declared_starts_only = true;
    auto evaluator = epcsaft_engine::equilibrium::compile_equilibrium_owner(
        model, incipient_model, std::move(request), parse_equilibrium_setup(setup_mapping));
    return std::make_shared<CompiledOwner>(
        evaluator, equilibrium_sidecar(evaluator->sidecar()));
}

py::dict solve_least_squares(
    const py::iterable& evaluator_values, const std::vector<double>& observations,
    const std::vector<double>& residual_scales, const std::string& residual_scale_metadata,
    const std::vector<double>& initial_parameters, const std::vector<double>& lower_bounds,
    const std::vector<double>& upper_bounds, const std::vector<double>& affine_origins,
    const std::vector<double>& affine_scales, const py::dict& control_values
) {
    epcsaft_engine::regression::LeastSquaresRequest request;
    request.evaluators = owners(evaluator_values);
    request.observations = observations;
    request.residual_scales = residual_scales;
    request.residual_scale_metadata = residual_scale_metadata;
    request.initial_parameters = initial_parameters;
    request.lower_bounds = lower_bounds;
    request.upper_bounds = upper_bounds;
    request.affine_origins = affine_origins;
    request.affine_scales = affine_scales;
    request.controls = {
        required<double>(control_values, "function_tolerance"),
        required<double>(control_values, "gradient_tolerance"),
        required<double>(control_values, "parameter_tolerance"),
        required<int>(control_values, "maximum_iterations"),
        required<double>(control_values, "maximum_elapsed_time_seconds"),
        required<std::size_t>(control_values, "maximum_uncached_owner_calls")};
    epcsaft_engine::regression::LeastSquaresResult solved;
    {
        py::gil_scoped_release release;
        solved = epcsaft_engine::regression::solve(request);
    }
    return solve_mapping(solved);
}

}  // namespace

void bind_regression(py::module_& module) {
    py::class_<CompiledOwner, std::shared_ptr<CompiledOwner>>(module, "CompiledOwner")
        .def_property_readonly("parameter_count", [](const CompiledOwner& value) {
            return value.owner->parameter_count();
        })
        .def_property_readonly("output_count", [](const CompiledOwner& value) {
            return value.owner->output_count();
        })
        .def_property_readonly("parameter_contract_fingerprint", [](const CompiledOwner& value) {
            return value.owner->parameter_contract_fingerprint();
        })
        .def_property_readonly("evaluator_fingerprint", [](const CompiledOwner& value) {
            return value.owner->evaluator_fingerprint();
        })
        .def_property_readonly("sidecar", [](const CompiledOwner& value) {
            return value.sidecar;
        });
    module.def("_compile_eos_regression_evaluator", &compile_eos,
               py::arg("model"), py::arg("specification"), py::arg("active_parameters"));
    module.def("_compile_equilibrium_regression_evaluator", &compile_equilibrium,
               py::arg("model"), py::arg("incipient_model"), py::arg("problem"),
               py::arg("source_reference_fingerprint"), py::arg("setup"));
    module.def("_solve_least_squares", &solve_least_squares,
               py::arg("evaluators"), py::arg("observations"),
               py::arg("residual_scales"), py::arg("residual_scale_metadata"),
               py::arg("initial_parameters"), py::arg("lower_bounds"),
               py::arg("upper_bounds"), py::arg("affine_origins"),
               py::arg("affine_scales"), py::arg("controls"));
}
