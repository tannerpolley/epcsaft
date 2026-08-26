#include "chemical_equilibrium.hpp"
#include "provider.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numeric>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "native_transport.hpp"

namespace py = pybind11;

namespace epcsaft_engine::equilibrium {
namespace {

struct ChemicalProviderMetadata {
    std::vector<std::string> component_ids;
    std::vector<int> charges;
};

const epcsaft_native_transport& checked_chemical_transport(
    const py::capsule& capsule
) {
    const char* name = capsule.name();
    if (name == nullptr
        || std::string_view(name) != EPCSAFT_NATIVE_TRANSPORT_CAPSULE_NAME) {
        throw py::value_error("expected private Engine native transport");
    }
    const auto* transport = capsule.get_pointer<epcsaft_native_transport>();
    if (transport == nullptr || transport->model_context == nullptr
        || transport->component_count < 2
        || transport->component_ids == nullptr
        || transport->component_charges == nullptr
        || (transport->evaluate_mixture_phase == nullptr
            && transport->evaluate_electrolyte_phase == nullptr)
        || transport->evaluate_molar_volume_bounds == nullptr
        || transport->evaluate_packing_fraction == nullptr) {
        throw py::value_error(
            "Provider capsule is missing the chemical-equilibrium phase callbacks"
        );
    }
    if (!std::isfinite(transport->source_temperature_min_k)
        || !std::isfinite(transport->source_temperature_max_k)
        || transport->source_temperature_min_k <= 0.0
        || transport->source_temperature_max_k
            < transport->source_temperature_min_k
        || (!std::isnan(transport->total_ion_mole_fraction_max)
            && (!std::isfinite(transport->total_ion_mole_fraction_max)
                || transport->total_ion_mole_fraction_max < 0.0
                || transport->total_ion_mole_fraction_max > 1.0))) {
        throw py::value_error("Provider source-domain metadata is invalid");
    }
    return *transport;
}

const epcsaft_native_transport& checked_reactive_finite_transport(
    const py::capsule& capsule
) {
    const char* name = capsule.name();
    if (name == nullptr
        || std::string_view(name) != EPCSAFT_NATIVE_TRANSPORT_CAPSULE_NAME) {
        throw py::value_error("expected private Engine native transport");
    }
    const auto* transport = capsule.get_pointer<epcsaft_native_transport>();
    if (transport == nullptr || transport->model_context == nullptr
        || transport->component_count < 2
        || transport->component_ids == nullptr
        || transport->component_charges == nullptr
        || transport->evaluate_mixture_phase == nullptr
        || transport->evaluate_molar_volume_bounds == nullptr
        || transport->evaluate_packing_fraction == nullptr
        || transport->evaluate_inverse_packing_geometry == nullptr) {
        throw py::value_error(
            "Provider capsule is missing the reactive finite-phase callbacks"
        );
    }
    const bool has_temperature_domain =
        std::isfinite(transport->source_temperature_min_k)
        && std::isfinite(transport->source_temperature_max_k);
    if ((std::isfinite(transport->source_temperature_min_k)
            != std::isfinite(transport->source_temperature_max_k))
        || (has_temperature_domain
            && (transport->source_temperature_min_k <= 0.0
                || transport->source_temperature_max_k
                    < transport->source_temperature_min_k))) {
        throw py::value_error("Provider source-domain metadata is invalid");
    }
    return *transport;
}

ChemicalProviderMetadata chemical_provider_metadata(const epcsaft_native_transport& transport) {
    ChemicalProviderMetadata result;
    result.component_ids.reserve(transport.component_count);
    result.charges.reserve(transport.component_count);
    for (std::size_t component = 0; component < transport.component_count; ++component) {
        if (transport.component_ids[component] == nullptr
            || transport.component_ids[component][0] == '\0') {
            throw py::value_error("Provider component identity is incomplete");
        }
        result.component_ids.emplace_back(transport.component_ids[component]);
        result.charges.push_back(static_cast<int>(transport.component_charges[component]));
    }
    return result;
}

DenseMatrix dense_matrix(const py::handle& value, const char* field) {
    const std::vector<std::vector<double>> rows = py::cast<std::vector<std::vector<double>>>(
        value
    );
    if (rows.empty()) {
        return {};
    }
    const std::size_t columns = rows.front().size();
    if (columns == 0) {
        throw py::value_error(std::string(field) + " rows must not be empty");
    }
    DenseMatrix result{rows.size(), columns, {}};
    result.values.reserve(rows.size() * columns);
    for (const std::vector<double>& row : rows) {
        if (row.size() != columns) {
            throw py::value_error(std::string(field) + " must be rectangular");
        }
        result.values.insert(result.values.end(), row.begin(), row.end());
    }
    return result;
}

DenseMatrix phase_direction_matrix(
    const std::vector<std::vector<double>>& rows,
    std::size_t expected_columns,
    const char* field
) {
    if (rows.empty()) {
        return {0, expected_columns, {}};
    }
    const std::size_t columns = rows.front().size();
    DenseMatrix result{rows.size(), columns, {}};
    result.values.reserve(rows.size() * columns);
    for (const std::vector<double>& row : rows) {
        if (row.size() != columns) {
            throw py::value_error(std::string(field) + " must be rectangular");
        }
        result.values.insert(result.values.end(), row.begin(), row.end());
    }
    if (columns != expected_columns) {
        throw py::value_error(
            std::string(field) + " columns must match the phase-species dimension"
        );
    }
    return result;
}

std::vector<EquilibriumConstantRecord> equilibrium_constant_records(
    const py::handle& value
) {
    std::vector<EquilibriumConstantRecord> result;
    for (const py::handle item : py::cast<py::tuple>(value)) {
        const py::dict record = py::cast<py::dict>(item);
        result.push_back({
            py::cast<std::string>(record["source_id"]),
            py::cast<std::string>(record["reference_id"]),
            py::cast<std::string>(record["reaction_orientation"]),
            py::cast<std::string>(record["conversion_id"]),
            py::cast<bool>(record["dimensionless"]),
            py::cast<double>(record["temperature_k"]),
            py::cast<double>(record["pressure_pa"]),
        });
    }
    return result;
}

ReactionSystemInput reaction_system_input(const py::dict& spec) {
    ReactionSystemInput input;
    input.species_ids = py::cast<std::vector<std::string>>(spec["species_ids"]);
    input.charges = py::cast<std::vector<int>>(spec["charges"]);
    if (spec.contains("provider_fingerprint")) {
        input.provider_fingerprint = py::cast<std::string>(spec["provider_fingerprint"]);
    }
    input.molar_masses_kg_per_mol = py::cast<std::vector<double>>(
        spec["molar_masses_kg_per_mol"]
    );
    input.balance_matrix = dense_matrix(spec["balance_matrix"], "balance matrix");
    input.conserved_totals = py::cast<std::vector<double>>(spec["conserved_totals"]);
    input.reaction_matrix = dense_matrix(spec["reaction_matrix"], "reaction matrix");
    if (input.reaction_matrix.rows == 0) {
        input.reaction_matrix.columns = input.species_ids.size();
    }
    input.feed_amounts = py::cast<std::vector<double>>(spec["feed_amounts"]);
    input.ln_k = py::cast<std::vector<double>>(spec["ln_k"]);
    input.equilibrium_constant_records = equilibrium_constant_records(
        spec["equilibrium_constant_records"]
    );
    input.temperature_k = py::cast<double>(spec["temperature_k"]);
    input.pressure_pa = py::cast<double>(spec["pressure_pa"]);
    return input;
}

ChemicalContinuationStart continuation_start(const py::dict& spec) {
    const py::dict payload = py::cast<py::dict>(spec["continuation_start"]);
    ChemicalContinuationStart start;
    start.amounts = py::cast<std::vector<double>>(payload["amounts"]);
    start.volume_m3 = py::cast<double>(payload["volume_m3"]);
    if (start.amounts.empty()
        || std::any_of(start.amounts.begin(), start.amounts.end(), [](double amount) {
            return !std::isfinite(amount) || amount <= 0.0;
        })
        || !std::isfinite(start.volume_m3)
        || start.volume_m3 <= 0.0) {
        throw py::value_error(
            "immutable continuation start must be strictly positive and finite"
        );
    }
    return start;
}

py::dict chemical_result(const ChemicalSolveResult& evaluation) {
    py::dict result;
    const auto optional_float = [](double value) -> py::object {
        return std::isfinite(value) ? py::cast(value) : py::none();
    };
    const auto optional_index = [](long value) -> py::object {
        return value >= 0 ? py::cast(value) : py::none();
    };
    result["accepted"] = evaluation.accepted;
    result["solver_status"] = evaluation.solver_status;
    result["callback_error"] = evaluation.callback_error;
    result["chemical_certification_level"] =
        evaluation.chemical_certification_level;
    result["boundary_status"] = evaluation.boundary_status;
    result["structural_zero_species_indices"] =
        evaluation.structural_zero_species_indices;
    result["numerical_status"] = evaluation.numerical_status;
    result["physical_status"] = evaluation.physical_status;
    result["provider_domain_status"] = evaluation.provider_domain_status;
    result["local_minimum_status"] = evaluation.local_minimum_status;
    result["negative_curvature_recovery_status"] =
        evaluation.negative_curvature_recovery_status;
    result["negative_curvature_recovery_attempts"] =
        evaluation.negative_curvature_recovery_attempts;
    result["negative_curvature_recovery_selected_sign"] =
        evaluation.negative_curvature_recovery_selected_sign;
    result["trace_status"] = evaluation.trace_status;
    result["globality_certificate"] = "not_guaranteed";
    result["amounts"] = evaluation.amounts;
    result["volume_m3"] = evaluation.volume_m3;
    result["local_variables"] = evaluation.local_variables;
    result["local_lower_bounds"] = evaluation.local_lower_bounds;
    result["local_upper_bounds"] = evaluation.local_upper_bounds;
    result["balance_inf_norm"] = evaluation.balance_inf_norm;
    result["charge_inf_norm"] = evaluation.charge_inf_norm;
    result["pressure_relative_residual"] = evaluation.pressure_relative_residual;
    result["reaction_affinity_inf_norm"] = evaluation.reaction_affinity_inf_norm;
    result["reaction_affinity_residuals"] =
        evaluation.reaction_affinity_residuals;
    result["packing_fraction"] = evaluation.packing_fraction;
    result["packing_fraction_min"] = optional_float(
        evaluation.packing_fraction_min
    );
    result["packing_fraction_max"] = optional_float(
        evaluation.packing_fraction_max
    );
    result["total_ion_fraction"] = optional_float(
        evaluation.total_ion_fraction
    );
    result["total_ion_fraction_max"] = optional_float(
        evaluation.total_ion_fraction_max
    );
    result["minimum_amount_mol"] = optional_float(
        evaluation.minimum_amount_mol
    );
    result["trace_floor_mol"] = optional_float(evaluation.trace_floor_mol);
    result["kkt_stationarity_inf_norm"] = evaluation.kkt_stationarity_inf_norm;
    result["physical_stationarity_residuals"] =
        evaluation.physical_stationarity_residuals;
    result["physical_equality_multipliers"] =
        evaluation.physical_equality_multipliers;
    result["physical_constraint_jacobian"] =
        evaluation.physical_constraint_jacobian;
    result["physical_constraint_shape"] = py::make_tuple(
        evaluation.physical_constraint_rows,
        evaluation.physical_constraint_columns
    );
    result["physical_lagrangian_gradient"] =
        evaluation.physical_lagrangian_gradient;
    result["physical_objective_gradient"] =
        evaluation.physical_objective_gradient;
    result["physical_to_chart_jacobian"] =
        evaluation.physical_to_chart_jacobian;
    result["physical_to_chart_shape"] = py::make_tuple(
        evaluation.physical_to_chart_rows,
        evaluation.physical_to_chart_columns
    );
    result["chart_physical_pullback_residual_inf_norm"] = optional_float(
        evaluation.chart_physical_pullback_residual_inf_norm
    );
    result["complementarity_inf_norm"] = evaluation.complementarity_inf_norm;
    result["kkt_dimension"] = evaluation.kkt_dimension;
    result["kkt_rank"] = evaluation.kkt_rank;
    result["condition_number_inf"] =
        optional_float(evaluation.condition_number_inf);
    result["failure_kind"] = evaluation.failure_kind;
    result["failure_reason"] = evaluation.failure_reason;
    result["active_lower_bounds"] = evaluation.active_lower_bounds;
    result["active_upper_bounds"] = evaluation.active_upper_bounds;
    result["active_constraint_bounds"] =
        evaluation.active_constraint_bounds;
    result["reduced_hessian_status"] = evaluation.reduced_hessian_status;
    result["reduced_hessian"] = evaluation.reduced_hessian;
    result["reduced_hessian_nullspace_basis"] =
        evaluation.reduced_hessian_nullspace_basis;
    result["reduced_hessian_nullspace_shape"] = py::make_tuple(
        evaluation.reduced_hessian_nullspace_rows,
        evaluation.reduced_hessian_nullspace_columns
    );
    result["reduced_hessian_eigenvalues"] =
        evaluation.reduced_hessian_eigenvalues;
    result["reduced_hessian_spectrum_status"] =
        evaluation.reduced_hessian_spectrum_status;
    result["reduced_hessian_raw_inertia"] = py::make_tuple(
        evaluation.reduced_hessian_raw_positive_eigenvalues,
        evaluation.reduced_hessian_raw_zero_eigenvalues,
        evaluation.reduced_hessian_raw_negative_eigenvalues
    );
    result["reduced_hessian_inertia"] = py::make_tuple(
        evaluation.reduced_hessian_positive_eigenvalues,
        evaluation.reduced_hessian_zero_eigenvalues,
        evaluation.reduced_hessian_negative_eigenvalues
    );
    result["reduced_hessian_scale"] =
        optional_float(evaluation.reduced_hessian_scale);
    result["reduced_hessian_eigenvalue_tolerance"] =
        optional_float(evaluation.reduced_hessian_eigenvalue_tolerance);
    result["objective_gradient"] = evaluation.objective_gradient;
    result["constraint_values"] = evaluation.constraint_values;
    result["constraint_jacobian"] = evaluation.constraint_jacobian;
    result["lagrangian_gradient"] = evaluation.lagrangian_gradient;
    result["equality_multipliers"] = evaluation.equality_multipliers;
    result["chart_stationarity_inf_norm"] =
        optional_float(evaluation.chart_stationarity_inf_norm);
    result["lagrangian_hessian"] = evaluation.lagrangian_hessian;
    result["covariant_lagrangian_hessian"] =
        evaluation.covariant_lagrangian_hessian;
    result["derivative_audits_performed"] =
        evaluation.derivative_audits_performed;
    result["derivative_check_step"] = optional_float(
        evaluation.derivative_check_step
    );
    result["objective_gradient_check_relative_error"] = optional_float(
        evaluation.objective_gradient_check_relative_error
    );
    result["physical_objective_gradient_check_relative_error"] = optional_float(
        evaluation.physical_objective_gradient_check_relative_error
    );
    result["constraint_jacobian_check_relative_error"] = optional_float(
        evaluation.constraint_jacobian_check_relative_error
    );
    result["lagrangian_hessian_check_relative_error"] = optional_float(
        evaluation.lagrangian_hessian_check_relative_error
    );
    result["reduced_hessian_check_relative_error"] = optional_float(
        evaluation.reduced_hessian_check_relative_error
    );
    result["derivative_check_worst_entry"] =
        evaluation.derivative_check_worst_entry;
    result["derivative_check_worst_relative_error"] = optional_float(
        evaluation.derivative_check_worst_relative_error
    );
    result["derivative_check_worst_analytic_value"] = optional_float(
        evaluation.derivative_check_worst_analytic_value
    );
    result["derivative_check_worst_finite_difference_value"] = optional_float(
        evaluation.derivative_check_worst_finite_difference_value
    );
    result["derivative_check_worst_step"] = optional_float(
        evaluation.derivative_check_worst_step
    );
    const std::size_t derivative_dimension =
        evaluation.objective_gradient.size();
    result["derivative_coordinate_order"] =
        evaluation.derivative_coordinate_order.size() == derivative_dimension
        ? py::cast(evaluation.derivative_coordinate_order)
        : py::cast(std::vector<std::string>{});
    result["derivative_objective_basis"] =
        evaluation.derivative_objective_basis;
    result["derivative_constraint_basis"] =
        evaluation.derivative_constraint_basis;
    result["derivative_constraint_order"] =
        evaluation.derivative_constraint_order.size()
            == evaluation.constraint_values.size()
        ? py::cast(evaluation.derivative_constraint_order)
        : py::cast(std::vector<std::string>{});
    result["objective_gradient_shape"] = py::make_tuple(
        derivative_dimension
    );
    result["constraint_jacobian_shape"] = py::make_tuple(
        evaluation.constraint_values.size(), derivative_dimension
    );
    result["lagrangian_hessian_shape"] = py::make_tuple(
        derivative_dimension, derivative_dimension
    );
    result["kkt_root_values"] = evaluation.kkt_root_values;
    result["kkt_root_jacobian"] = evaluation.kkt_root_jacobian;
    result["kkt_root_shape"] = py::make_tuple(
        evaluation.kkt_root_rows,
        evaluation.kkt_root_columns
    );
    result["kkt_root_status"] = evaluation.kkt_root_status;
    result["kkt_root_jacobian_check_relative_error"] = optional_float(
        evaluation.kkt_root_jacobian_check_relative_error
    );
    const auto criterion = [](
        const char* name,
        double value,
        double limit,
        const char* comparison,
        bool passed
    ) {
        py::dict record;
        record["name"] = name;
        record["status"] = std::isfinite(value) && passed ? "passed" : "failed";
        record["value"] = value;
        record["limit"] = limit;
        record["comparison"] = comparison;
        return record;
    };
    py::list numerical_criteria;
    numerical_criteria.append(criterion(
        "solver_termination", evaluation.solver_status == "solve_succeeded" ? 1.0 : 0.0,
        1.0, ">=", evaluation.solver_status == "solve_succeeded"
    ));
    numerical_criteria.append(criterion(
        "balance_inf_norm", evaluation.balance_inf_norm,
        kChemicalBalanceTolerance, "<=",
        evaluation.balance_inf_norm <= kChemicalBalanceTolerance
    ));
    numerical_criteria.append(criterion(
        "complementarity_inf_norm", evaluation.complementarity_inf_norm,
        kChemicalKktTolerance, "<=",
        evaluation.complementarity_inf_norm <= kChemicalKktTolerance
    ));
    numerical_criteria.append(criterion(
        "kkt_stationarity_inf_norm", evaluation.kkt_stationarity_inf_norm,
        kChemicalKktTolerance, "<=",
        evaluation.kkt_stationarity_inf_norm <= kChemicalKktTolerance
    ));
    numerical_criteria.append(criterion(
        "chart_stationarity_inf_norm", evaluation.chart_stationarity_inf_norm,
        kChemicalKktTolerance, "<=",
        evaluation.chart_stationarity_inf_norm <= kChemicalKktTolerance
    ));
    numerical_criteria.append(criterion(
        "chart_physical_pullback_residual_inf_norm",
        evaluation.chart_physical_pullback_residual_inf_norm,
        kChemicalChartPullbackTolerance,
        "<=",
        evaluation.chart_physical_pullback_residual_inf_norm
            <= kChemicalChartPullbackTolerance
    ));
    numerical_criteria.append(criterion(
        "derivative_evidence_finite", evaluation.callback_error.empty() ? 1.0 : 0.0,
        1.0, ">=", evaluation.callback_error.empty()
    ));
    if (evaluation.derivative_audits_performed) {
        numerical_criteria.append(criterion(
            "objective_gradient_spanning_relative_error",
            evaluation.objective_gradient_check_relative_error,
            kChemicalObjectiveDerivativeTolerance,
            "<=",
            evaluation.objective_gradient_check_relative_error
                <= kChemicalObjectiveDerivativeTolerance
        ));
        numerical_criteria.append(criterion(
            "physical_objective_gradient_spanning_relative_error",
            evaluation.physical_objective_gradient_check_relative_error,
            kChemicalObjectiveDerivativeTolerance,
            "<=",
            evaluation.physical_objective_gradient_check_relative_error
                <= kChemicalObjectiveDerivativeTolerance
        ));
        numerical_criteria.append(criterion(
            "constraint_jacobian_spanning_relative_error",
            evaluation.constraint_jacobian_check_relative_error,
            kChemicalConstraintDerivativeTolerance,
            "<=",
            evaluation.constraint_jacobian_check_relative_error
                <= kChemicalConstraintDerivativeTolerance
        ));
        numerical_criteria.append(criterion(
            "lagrangian_hessian_spanning_relative_error",
            evaluation.lagrangian_hessian_check_relative_error,
            kChemicalHessianDerivativeTolerance,
            "<=",
            evaluation.lagrangian_hessian_check_relative_error
                <= kChemicalHessianDerivativeTolerance
        ));
        numerical_criteria.append(criterion(
            "reduced_hessian_spanning_relative_error",
            evaluation.reduced_hessian_check_relative_error,
            kChemicalHessianDerivativeTolerance,
            "<=",
            evaluation.reduced_hessian_check_relative_error
                <= kChemicalHessianDerivativeTolerance
        ));
        numerical_criteria.append(criterion(
            "kkt_root_jacobian_spanning_relative_error",
            evaluation.kkt_root_jacobian_check_relative_error,
            kChemicalHessianDerivativeTolerance,
            "<=",
            evaluation.kkt_root_jacobian_check_relative_error
                <= kChemicalHessianDerivativeTolerance
        ));
    }
    numerical_criteria.append(criterion(
        "kkt_rank",
        static_cast<double>(evaluation.kkt_rank),
        static_cast<double>(evaluation.kkt_dimension),
        ">=",
        evaluation.kkt_rank == evaluation.kkt_dimension
    ));
    numerical_criteria.append(criterion(
        "kkt_condition_number_inf",
        evaluation.condition_number_inf,
        kChemicalConditionNumberMax,
        "<=",
        evaluation.condition_number_inf <= kChemicalConditionNumberMax
    ));
    numerical_criteria.append(criterion(
        "strict_local_minimum",
        evaluation.local_minimum_status == "passed" ? 1.0 : 0.0,
        1.0,
        ">=",
        evaluation.local_minimum_status == "passed"
    ));
    result["numerical_criteria"] = std::move(numerical_criteria);
    py::list physical_criteria;
    physical_criteria.append(criterion(
        "balance_inf_norm", evaluation.balance_inf_norm,
        kChemicalBalanceTolerance, "<=",
        evaluation.balance_inf_norm <= kChemicalBalanceTolerance
    ));
    physical_criteria.append(criterion(
        "charge_inf_norm", evaluation.charge_inf_norm,
        kChemicalBalanceTolerance, "<=",
        evaluation.charge_inf_norm <= kChemicalBalanceTolerance
    ));
    physical_criteria.append(criterion(
        "pressure_relative_residual",
        evaluation.pressure_relative_residual,
        kChemicalPressureTolerance,
        "<=",
        evaluation.pressure_relative_residual <= kChemicalPressureTolerance
    ));
    physical_criteria.append(criterion(
        "reaction_affinity_inf_norm",
        evaluation.reaction_affinity_inf_norm,
        kChemicalAffinityTolerance,
        "<=",
        evaluation.reaction_affinity_inf_norm <= kChemicalAffinityTolerance
    ));
    if (std::isfinite(evaluation.minimum_amount_mol)
        && std::isfinite(evaluation.trace_floor_mol)) {
        physical_criteria.append(criterion(
            "minimum_amount_mol",
            evaluation.minimum_amount_mol,
            evaluation.trace_floor_mol,
            ">",
            evaluation.minimum_amount_mol > evaluation.trace_floor_mol
        ));
    }
    if (std::isfinite(evaluation.packing_fraction_min)
        && std::isfinite(evaluation.packing_fraction_max)) {
        physical_criteria.append(criterion(
            "packing_fraction_lower_bound",
            evaluation.packing_fraction,
            evaluation.packing_fraction_min,
            ">=",
            evaluation.packing_fraction >= evaluation.packing_fraction_min
        ));
        physical_criteria.append(criterion(
            "packing_fraction_upper_bound",
            evaluation.packing_fraction,
            evaluation.packing_fraction_max,
            "<=",
            evaluation.packing_fraction <= evaluation.packing_fraction_max
        ));
    }
    if (std::isfinite(evaluation.total_ion_fraction_max)) {
        physical_criteria.append(criterion(
            "total_ion_fraction",
            evaluation.total_ion_fraction,
            evaluation.total_ion_fraction_max + 1.0e-12,
            "<=",
            evaluation.total_ion_fraction
                <= evaluation.total_ion_fraction_max + 1.0e-12
        ));
    }
    result["physical_criteria"] = std::move(physical_criteria);
    py::dict search;
    search["status"] = evaluation.search.status;
    search["primary_budget"] = evaluation.search.primary_budget;
    search["primary_attempt_count"] =
        evaluation.search.primary_attempt_count;
    search["generated_start_count"] =
        evaluation.search.generated_start_count;
    search["budget_truncated_start_count"] =
        evaluation.search.budget_truncated_start_count;
    search["duplicate_start_count"] =
        evaluation.search.duplicate_start_count;
    search["infeasible_start_count"] =
        evaluation.search.infeasible_start_count;
    search["evaluated_start_count"] =
        evaluation.search.evaluated_start_count;
    search["domain_rejected_start_count"] =
        evaluation.search.domain_rejected_start_count;
    search["construction_rejected_start_count"] =
        evaluation.search.construction_rejected_start_count;
    search["selected_basin_ordinal"] =
        optional_index(evaluation.search.selected_basin_ordinal);
    search["selected_objective"] =
        optional_float(evaluation.search.selected_objective);
    search["selection_label"] = evaluation.search.selection_label;
    py::list attempts;
    for (const ChemicalSearchAttempt& attempt : evaluation.search.attempts) {
        py::dict record;
        record["ordinal"] = attempt.ordinal;
        record["primary_ordinal"] = attempt.primary_ordinal;
        record["kind"] = attempt.kind;
        record["parent_ordinal"] = optional_index(attempt.parent_ordinal);
        record["start_identity"] = attempt.start_identity;
        record["start_construction_status"] =
            attempt.start_construction_status;
        record["retraction_status"] = attempt.retraction_status;
        record["provider_domain_status"] =
            attempt.provider_domain_status;
        record["solver_status"] = attempt.solver_status;
        record["callback_error"] = attempt.callback_error;
        record["terminal_status"] = attempt.terminal_status;
        record["amounts"] = attempt.amounts;
        record["volume_m3"] = optional_float(attempt.volume_m3);
        record["objective"] = optional_float(attempt.objective);
        record["balance_inf_norm"] =
            optional_float(attempt.balance_inf_norm);
        record["charge_inf_norm"] =
            optional_float(attempt.charge_inf_norm);
        record["pressure_relative_residual"] =
            optional_float(attempt.pressure_relative_residual);
        record["reaction_affinity_inf_norm"] =
            optional_float(attempt.reaction_affinity_inf_norm);
        record["kkt_stationarity_inf_norm"] =
            optional_float(attempt.kkt_stationarity_inf_norm);
        record["complementarity_inf_norm"] =
            optional_float(attempt.complementarity_inf_norm);
        record["kkt_dimension"] = attempt.kkt_dimension;
        record["kkt_rank"] = attempt.kkt_rank;
        record["condition_number_inf"] =
            optional_float(attempt.condition_number_inf);
        record["local_minimum_status"] =
            attempt.local_minimum_status;
        record["trace_status"] = attempt.trace_status;
        record["basin_ordinal"] = optional_index(attempt.basin_ordinal);
        record["recovery_seed_count"] = attempt.recovery_seed_count;
        record["recovery_solve_count"] = attempt.recovery_solve_count;
        attempts.append(std::move(record));
    }
    search["attempts"] = std::move(attempts);
    py::list basins;
    for (const ChemicalSearchBasin& basin : evaluation.search.basins) {
        py::dict record;
        record["ordinal"] = basin.ordinal;
        record["representative_attempt_ordinal"] =
            basin.representative_attempt_ordinal;
        record["amounts"] = basin.amounts;
        record["volume_m3"] = optional_float(basin.volume_m3);
        record["objective"] = optional_float(basin.objective);
        basins.append(std::move(record));
    }
    search["basins"] = std::move(basins);
    py::list prefixes;
    for (const ChemicalSearchBudgetPrefix& prefix
         : evaluation.search.budget_prefixes) {
        py::dict record;
        record["primary_budget"] = prefix.primary_budget;
        record["attempted_primary_ordinals"] =
            prefix.attempted_primary_ordinals;
        record["basin_ordinals"] = prefix.basin_ordinals;
        record["selected_basin_ordinal"] =
            optional_index(prefix.selected_basin_ordinal);
        record["selection_changed"] = prefix.selection_changed;
        prefixes.append(std::move(record));
    }
    search["budget_prefixes"] = std::move(prefixes);
    result["search"] = std::move(search);
    py::dict sensitivities;
    sensitivities["status"] = evaluation.sensitivities.status;
    sensitivities["failure_reason"] = evaluation.sensitivities.failure_reason;
    const std::size_t species_count = evaluation.amounts.size();
    if (evaluation.sensitivities.amount_derivatives.size()
            != evaluation.sensitivities.parameter_order.size() * species_count
        || evaluation.sensitivities.volume_derivatives.size()
            != evaluation.sensitivities.parameter_order.size()) {
        throw std::runtime_error(
            "internal chemical-sensitivity result dimensions are inconsistent"
        );
    }
    py::tuple parameter_order(evaluation.sensitivities.parameter_order.size());
    for (std::size_t parameter = 0;
         parameter < evaluation.sensitivities.parameter_order.size();
         ++parameter) {
        parameter_order[parameter] =
            evaluation.sensitivities.parameter_order[parameter];
    }
    sensitivities["parameter_order"] = std::move(parameter_order);
    py::tuple amount_derivatives(evaluation.sensitivities.parameter_order.size());
    for (std::size_t parameter = 0;
         parameter < evaluation.sensitivities.parameter_order.size();
         ++parameter) {
        py::tuple row(species_count);
        for (std::size_t species = 0; species < species_count; ++species) {
            row[species] = evaluation.sensitivities.amount_derivatives[
                parameter * species_count + species
            ];
        }
        amount_derivatives[parameter] = std::move(row);
    }
    sensitivities["amount_derivatives"] = std::move(amount_derivatives);
    sensitivities["volume_derivatives"] = py::tuple(py::cast(
        evaluation.sensitivities.volume_derivatives
    ));
    sensitivities["kkt_parameter_jacobian"] = py::tuple(py::cast(
        evaluation.sensitivities.kkt_parameter_jacobian
    ));
    sensitivities["kkt_dimension"] = evaluation.sensitivities.kkt_dimension;
    sensitivities["kkt_rank"] = evaluation.sensitivities.kkt_rank;
    sensitivities["condition_number_inf"] =
        evaluation.sensitivities.condition_number_inf;
    sensitivities["active_lower_bounds"] = py::tuple(py::cast(
        evaluation.sensitivities.active_lower_bounds
    ));
    sensitivities["active_upper_bounds"] = py::tuple(py::cast(
        evaluation.sensitivities.active_upper_bounds
    ));
    sensitivities["active_constraint_bounds"] = py::tuple(py::cast(
        evaluation.sensitivities.active_constraint_bounds
    ));
    sensitivities["active_trace_species"] = py::tuple(py::cast(
        evaluation.sensitivities.active_trace_species
    ));
    sensitivities["chart_topology"] =
        evaluation.sensitivities.chart_topology;
    sensitivities["parameter_fingerprint"] =
        evaluation.sensitivities.parameter_fingerprint;
    sensitivities["provider_parameter_status"] =
        evaluation.sensitivities.provider_parameter_status;
    sensitivities["provider_parameter_failure_reason"] =
        evaluation.sensitivities.provider_parameter_failure_reason;
    sensitivities["reference_parameter_status"] =
        evaluation.sensitivities.reference_parameter_status;
    sensitivities["reference_parameter_failure_reason"] =
        evaluation.sensitivities.reference_parameter_failure_reason;
    result["sensitivities"] = std::move(sensitivities);
    return result;
}

py::dict solve_manufactured(const py::dict& spec, double trace_floor) {
    const ReactionSystemInput input = reaction_system_input(spec);
    const CompiledReactionSystem compiled = compile_reaction_system(input);
    return chemical_result(
        solve_manufactured_ideal_reaction(
            compiled,
            input.temperature_k,
            input.pressure_pa,
            {},
            trace_floor
        )
    );
}

void validate_provider_identity(
    const ChemicalProviderMetadata& metadata,
    const ReactionSystemInput& input
) {
    if (input.species_ids != metadata.component_ids) {
        throw py::value_error("Provider capsule component order does not match the reaction system");
    }
    if (input.charges != metadata.charges) {
        throw py::value_error("Provider capsule charges do not match the reaction system");
    }
}

CompiledReactionSystem compile_provider_input(
    const epcsaft_native_transport& transport,
    const ChemicalProviderMetadata& metadata,
    const ReactionSystemInput& input,
    const std::vector<double>& packing_bounds
) {
    validate_provider_identity(metadata, input);
    if (input.temperature_k < transport.source_temperature_min_k
        || input.temperature_k > transport.source_temperature_max_k) {
        throw py::value_error("temperature is outside the Provider source domain");
    }
    const CompiledReactionSystem compiled = compile_reaction_system(input);
    double ionic_feed = 0.0;
    double total_feed = 0.0;
    for (std::size_t species = 0; species < input.feed_amounts.size(); ++species) {
        total_feed += input.feed_amounts[species];
        if (input.charges[species] != 0) {
            ionic_feed += input.feed_amounts[species];
        }
    }
    if (std::isfinite(transport.total_ion_mole_fraction_max)
        && (total_feed <= 0.0
            || ionic_feed / total_feed > transport.total_ion_mole_fraction_max)) {
        throw py::value_error("feed composition exceeds the Provider source domain");
    }
    if (packing_bounds.size() != 2) {
        throw py::value_error("packing_fraction_bounds must contain two values");
    }
    return compiled;
}

void set_provider_parameter_status(
    ChemicalSolveResult& evaluation,
    std::size_t active_parameter_count
) {
    if (active_parameter_count == 0) {
        evaluation.sensitivities.provider_parameter_status = "not_applicable";
        evaluation.sensitivities.provider_parameter_failure_reason.clear();
    } else if (evaluation.sensitivities.status == "available") {
        evaluation.sensitivities.provider_parameter_status = "available";
        evaluation.sensitivities.provider_parameter_failure_reason.clear();
    } else {
        evaluation.sensitivities.provider_parameter_status = "unavailable";
        evaluation.sensitivities.provider_parameter_failure_reason =
            evaluation.sensitivities.failure_reason;
    }
}

py::dict solve_provider_input(
    const epcsaft_native_transport& transport,
    const ChemicalProviderMetadata& metadata,
    const ReactionSystemInput& input,
    const ProviderContext& provider,
    const std::vector<double>& packing_bounds,
    double trace_floor,
    const std::vector<double>& ln_k_pressure_derivatives_per_pa = {},
    const std::vector<double>& ln_k_parameter_derivatives = {},
    const ProviderActiveParameterSet* active_parameters = nullptr,
    const ChemicalContinuationStart* continuation_start = nullptr
) {
    const CompiledReactionSystem compiled = compile_provider_input(
        transport, metadata, input, packing_bounds
    );
    std::vector<double> compiled_ln_k_pressure_derivatives_per_pa;
    if (!ln_k_pressure_derivatives_per_pa.empty()) {
        if (ln_k_pressure_derivatives_per_pa.size()
                != input.reaction_matrix.rows
            || compiled.supplied_reaction_transform.rows
                != compiled.reaction_matrix.rows
            || compiled.supplied_reaction_transform.columns
                != input.reaction_matrix.rows) {
            throw py::value_error(
                "source pressure derivatives do not match the compiled reaction basis"
            );
        }
        compiled_ln_k_pressure_derivatives_per_pa.assign(
            compiled.reaction_matrix.rows, 0.0
        );
        for (std::size_t reaction = 0;
             reaction < compiled.reaction_matrix.rows;
             ++reaction) {
            for (std::size_t supplied = 0;
                 supplied < input.reaction_matrix.rows;
                 ++supplied) {
                compiled_ln_k_pressure_derivatives_per_pa[reaction] +=
                    compiled.supplied_reaction_transform(reaction, supplied)
                    * ln_k_pressure_derivatives_per_pa[supplied];
            }
        }
    }
    std::vector<double> compiled_ln_k_parameter_derivatives;
    const std::size_t active_parameter_count = active_parameters == nullptr
        ? 0
        : active_parameters->identities.size();
    if (!ln_k_parameter_derivatives.empty()) {
        if (active_parameter_count == 0
            || ln_k_parameter_derivatives.size()
                != input.reaction_matrix.rows * active_parameter_count
            || compiled.supplied_reaction_transform.rows
                != compiled.reaction_matrix.rows
            || compiled.supplied_reaction_transform.columns
                != input.reaction_matrix.rows) {
            throw py::value_error(
                "source parameter derivatives do not match the compiled reaction basis"
            );
        }
        compiled_ln_k_parameter_derivatives.assign(
            compiled.reaction_matrix.rows * active_parameter_count, 0.0
        );
        for (std::size_t reaction = 0;
             reaction < compiled.reaction_matrix.rows;
             ++reaction) {
            for (std::size_t supplied = 0;
                 supplied < input.reaction_matrix.rows;
                 ++supplied) {
                for (std::size_t parameter = 0;
                     parameter < active_parameter_count;
                     ++parameter) {
                    compiled_ln_k_parameter_derivatives[
                        reaction * active_parameter_count + parameter
                    ] += compiled.supplied_reaction_transform(reaction, supplied)
                        * ln_k_parameter_derivatives[
                            supplied * active_parameter_count + parameter
                        ];
                }
            }
        }
    } else if (active_parameter_count != 0) {
        compiled_ln_k_parameter_derivatives.assign(
            compiled.reaction_matrix.rows * active_parameter_count, 0.0
        );
    }
    ChemicalSolveResult evaluation = solve_provider_reaction(
        compiled,
        provider,
        input.temperature_k,
        input.pressure_pa,
        packing_bounds[0],
        packing_bounds[1],
        transport.total_ion_mole_fraction_max,
        trace_floor,
        compiled_ln_k_pressure_derivatives_per_pa,
        compiled_ln_k_parameter_derivatives,
        active_parameters,
        continuation_start
    );
    set_provider_parameter_status(evaluation, active_parameter_count);
    py::dict result = chemical_result(evaluation);
    result["parameter_fingerprint"] = input.provider_fingerprint;
    result["packing_fraction_bounds"] = packing_bounds;
    result["provider_transport_capsule_name"] =
        EPCSAFT_NATIVE_TRANSPORT_CAPSULE_NAME;
    result["provider_evaluation_count"] = provider.evaluation_count();
    return result;
}

SourceReferenceTransferEvidence source_reference_transfer_evidence(
    const py::dict& source_standard_state,
    const char* key
) {
    if (!source_standard_state.contains(key)) {
        throw py::value_error(
            std::string("installed Provider source-reference transfer is missing: ")
            + key
        );
    }
    const py::dict payload = py::cast<py::dict>(source_standard_state[key]);
    SourceReferenceTransferEvidence result;
    result.status = py::cast<std::string>(payload["status"]);
    result.domain_status = py::cast<std::string>(payload["domain_status"]);
    result.convergence_status = py::cast<std::string>(
        payload["convergence_status"]
    );
    result.reference_state_id = py::cast<std::string>(
        payload["reference_state_id"]
    );
    result.activity_convention_id = py::cast<std::string>(
        payload["activity_convention_id"]
    );
    result.component_ids = py::cast<std::vector<std::string>>(
        payload["component_ids"]
    );
    result.neutral_basis_row_count = py::cast<std::size_t>(
        payload["neutral_basis_row_count"]
    );
    result.neutral_basis = py::cast<std::vector<double>>(
        payload["neutral_basis"]
    );
    result.log_fugacity_contractions = py::cast<std::vector<double>>(
        payload["log_fugacity_contractions"]
    );
    result.activity_scale_log_contractions = py::cast<std::vector<double>>(
        payload["activity_scale_log_contractions"]
    );
    result.transfer_log_contractions = py::cast<std::vector<double>>(
        payload["transfer_log_contractions"]
    );
    result.source_composition = py::cast<std::vector<double>>(
        payload["source_composition"]
    );
    result.derivative_availability = py::cast<std::vector<std::string>>(
        payload["derivative_availability"]
    );
    result.pressure_derivatives_per_pa = py::cast<std::vector<double>>(
        payload["pressure_derivatives_per_pa"]
    );
    result.active_parameter_identities = py::cast<std::vector<std::string>>(
        payload["active_parameter_identities"]
    );
    result.parameter_derivatives = py::cast<std::vector<double>>(
        payload["parameter_derivatives"]
    );
    result.temperature_k = py::cast<double>(payload["temperature_k"]);
    result.pressure_pa = py::cast<double>(payload["pressure_pa"]);
    result.source_reference_pressure_pa = py::cast<double>(
        payload["source_reference_pressure_pa"]
    );
    result.standard_molality_mol_per_kg = py::cast<double>(
        payload["standard_molality_mol_per_kg"]
    );
    result.reference_convergence_error = py::cast<double>(
        payload["reference_convergence_error"]
    );
    result.source_temperature_min_k = py::cast<double>(
        payload["source_temperature_min_k"]
    );
    result.source_temperature_max_k = py::cast<double>(
        payload["source_temperature_max_k"]
    );
    result.source_pressure_min_pa = py::cast<double>(
        payload["source_pressure_min_pa"]
    );
    result.source_pressure_max_pa = py::cast<double>(
        payload["source_pressure_max_pa"]
    );
    result.parameter_fingerprint = py::cast<std::string>(
        payload["parameter_fingerprint"]
    );
    result.topology_fingerprint = py::cast<std::string>(
        payload["topology_fingerprint"]
    );
    result.component_order_fingerprint = py::cast<std::string>(
        payload["component_order_fingerprint"]
    );
    result.reference_state_fingerprint = py::cast<std::string>(
        payload["reference_state_fingerprint"]
    );
    result.domain_fingerprint = py::cast<std::string>(
        payload["domain_fingerprint"]
    );
    result.artifact_fingerprint = py::cast<std::string>(
        payload["artifact_fingerprint"]
    );
    result.helmholtz_basis_id = py::cast<std::string>(
        payload["helmholtz_basis_id"]
    );
    const auto fingerprint_valid = [](const std::string& fingerprint) {
        return fingerprint.size() == 71
            && fingerprint.rfind("sha256:", 0) == 0
            && std::all_of(
                fingerprint.begin() + 7,
                fingerprint.end(),
                [](char value) {
                    return (value >= '0' && value <= '9')
                        || (value >= 'a' && value <= 'f');
                }
            );
    };
    if (!fingerprint_valid(result.parameter_fingerprint)
        || !fingerprint_valid(result.topology_fingerprint)
        || !fingerprint_valid(result.component_order_fingerprint)
        || !fingerprint_valid(result.reference_state_fingerprint)
        || !fingerprint_valid(result.domain_fingerprint)
        || !fingerprint_valid(result.artifact_fingerprint)) {
        throw py::value_error(
            "installed Provider source-reference fingerprints are incomplete"
        );
    }
    return result;
}

py::dict solve_provider_source(
    const py::capsule& capsule,
    const py::dict& spec,
    const py::dict& source_standard_state,
    const std::vector<double>& packing_bounds,
    double trace_floor,
    bool sensitivities_requested,
    const ProviderActiveParameterSet* active_parameters,
    const ChemicalContinuationStart* continuation_start
) {
    const epcsaft_native_transport& transport = checked_chemical_transport(capsule);
    const ChemicalProviderMetadata metadata = chemical_provider_metadata(transport);
    ReactionSystemInput input = reaction_system_input(spec);
    validate_provider_identity(metadata, input);
    const std::string source_standard_state_id = py::cast<std::string>(
        source_standard_state["id"]
    );
    const std::string source_activity_scale_id = py::cast<std::string>(
        source_standard_state["activity_scale_id"]
    );
    const std::vector<double> log_activity_scale_factors =
        py::cast<std::vector<double>>(
            source_standard_state["log_activity_scale_factors"]
        );
    const double source_reference_pressure_pa = py::cast<double>(
        source_standard_state["reference_pressure_pa"]
    );
    const std::string source_reference_activity_convention_id =
        py::cast<std::string>(
            source_standard_state["source_reference_activity_convention_id"]
        );
    const double source_reference_standard_molality_mol_per_kg =
        py::cast<double>(
            source_standard_state[
                "source_reference_standard_molality_mol_per_kg"
            ]
        );
    if (source_standard_state_id.empty() || source_activity_scale_id.empty()
        || source_reference_activity_convention_id.empty()) {
        throw py::value_error("source standard-state identity is incomplete");
    }
    if (!std::isfinite(source_reference_pressure_pa)
        || source_reference_pressure_pa <= 0.0) {
        throw py::value_error(
            "source standard-state reference pressure is invalid"
        );
    }
    if (!std::isfinite(source_reference_standard_molality_mol_per_kg)
        || source_reference_standard_molality_mol_per_kg <= 0.0) {
        throw py::value_error("source standard molality is invalid");
    }
    for (const EquilibriumConstantRecord& record : input.equilibrium_constant_records) {
        if (record.conversion_id
                != "source-standard-state-to-provider-neutral-reference"
            || record.reaction_orientation != "products_positive"
            || !record.dimensionless
            || record.reference_id.empty()
            || record.temperature_k != input.temperature_k
            || record.pressure_pa != source_reference_pressure_pa) {
            throw py::value_error(
                "source equilibrium-constant provenance is incompatible"
            );
        }
        if (source_standard_state_id != record.reference_id) {
            throw py::value_error(
                "source standard-state identity is inconsistent"
            );
        }
    }
    const ProviderContext provider(transport);
    input.provider_fingerprint = provider.fingerprint();
    const std::string provider_topology_fingerprint =
        provider.parameterized_topology_fingerprint();
    const SourceReferenceTransferEvidence reference =
        source_reference_transfer_evidence(
            source_standard_state, "provider_transfer"
        );
    const bool pressure_derivative_available = std::find(
        reference.derivative_availability.begin(),
        reference.derivative_availability.end(),
        "pressure"
    ) != reference.derivative_availability.end();
    const bool parameter_derivative_available = std::find(
        reference.derivative_availability.begin(),
        reference.derivative_availability.end(),
        "parameters"
    ) != reference.derivative_availability.end();
    if (sensitivities_requested != pressure_derivative_available
        || (active_parameters != nullptr) != parameter_derivative_available
        || (active_parameters != nullptr
            && reference.active_parameter_identities
                != active_parameters->identities)) {
        throw py::value_error(
            "installed Provider source-reference derivative receipt is incompatible"
        );
    }
    if (reference.reference_state_id != source_standard_state_id
        || reference.activity_convention_id
            != source_reference_activity_convention_id
        || reference.source_reference_pressure_pa
            != source_reference_pressure_pa
        || reference.standard_molality_mol_per_kg
            != source_reference_standard_molality_mol_per_kg) {
        throw py::value_error(
            "installed Provider source-reference declaration is inconsistent"
        );
    }
    if (reference.topology_fingerprint != provider_topology_fingerprint
        || reference.source_temperature_min_k != transport.source_temperature_min_k
        || reference.source_temperature_max_k != transport.source_temperature_max_k
        || !std::isfinite(reference.source_pressure_min_pa)
        || !std::isfinite(reference.source_pressure_max_pa)
        || reference.source_pressure_min_pa <= 0.0
        || reference.source_pressure_max_pa < reference.source_pressure_min_pa
        || input.pressure_pa < reference.source_pressure_min_pa
        || input.pressure_pa > reference.source_pressure_max_pa) {
        throw py::value_error(
            "installed Provider source-reference topology or domain is incompatible"
        );
    }
    const SourceStandardStateResult transformed = transform_source_standard_state(
        input.reaction_matrix,
        input.ln_k,
        log_activity_scale_factors,
        input.species_ids,
        input.charges,
        input.provider_fingerprint,
        input.temperature_k,
        input.pressure_pa,
        reference
    );
    input.ln_k = transformed.ln_k_provider_basis;
    for (EquilibriumConstantRecord& record : input.equilibrium_constant_records) {
        record.reference_id = "provider-helmholtz-coordinate-basis";
        record.conversion_id = "already-provider-basis";
        record.pressure_pa = input.pressure_pa;
    }
    py::dict result = solve_provider_input(
        transport,
        metadata,
        input,
        provider,
        packing_bounds,
        trace_floor,
        transformed.ln_k_pressure_derivatives_per_pa,
        transformed.ln_k_parameter_derivatives,
        active_parameters,
        continuation_start
    );
    py::dict sensitivities = py::cast<py::dict>(result["sensitivities"]);
    sensitivities["reference_parameter_status"] =
        active_parameters == nullptr ? "not_applicable" : "available";
    sensitivities["reference_parameter_failure_reason"] = "";
    result["standard_offsets"] = transformed.standard_offsets;
    result["ln_k_provider_basis"] = transformed.ln_k_provider_basis;
    result["reference_representation_residual_inf_norm"] =
        transformed.representation_residual_inf_norm;
    result["source_standard_state_id"] = source_standard_state_id;
    result["source_activity_scale_id"] = source_activity_scale_id;
    result["provider_reference_id"] = reference.helmholtz_basis_id;
    result["reference_derivative_availability"] =
        (pressure_derivative_available ? 2 : 0)
        | (parameter_derivative_available ? 8 : 0);
    result["reference_convergence_error"] =
        reference.reference_convergence_error;
    result["source_reference_state_fingerprint"] =
        reference.reference_state_fingerprint;
    return result;
}

py::dict solve_chemical_equilibrium(
    const py::object& capsule,
    const py::dict& spec,
    const py::object& source_standard_state,
    const py::object& packing_fraction_bounds,
    double trace_floor,
    bool sensitivities_requested,
    bool active_parameters_requested
) {
    if (capsule.is_none()) {
        if (!source_standard_state.is_none() || !packing_fraction_bounds.is_none()) {
            throw py::value_error(
                "ideal-gas chemical equilibrium cannot consume Provider metadata"
            );
        }
        return solve_manufactured(spec, trace_floor);
    }
    if (packing_fraction_bounds.is_none()) {
        throw py::value_error(
            "packing_fraction_bounds are required from the calling formulation"
        );
    }
    const std::vector<double> packing_bounds =
        py::cast<std::vector<double>>(packing_fraction_bounds);
    const ChemicalContinuationStart anchored = spec.contains("continuation_start")
        ? continuation_start(spec)
        : ChemicalContinuationStart{};
    const ChemicalContinuationStart* anchored_ptr =
        spec.contains("continuation_start") ? &anchored : nullptr;
    const py::capsule provider_capsule = py::cast<py::capsule>(capsule);
    if (source_standard_state.is_none()) {
        const epcsaft_native_transport& transport =
            checked_chemical_transport(provider_capsule);
        const ChemicalProviderMetadata metadata = chemical_provider_metadata(transport);
        ReactionSystemInput input = reaction_system_input(spec);
        const ProviderContext provider(transport);
        input.provider_fingerprint = provider.fingerprint();
        const ProviderActiveParameterSet resolved = provider.active_parameters();
        if (active_parameters_requested != !resolved.identities.empty()) {
            throw py::value_error(
                "Provider capsule active-parameter identity does not match the request"
            );
        }
        return solve_provider_input(
            transport,
            metadata,
            input,
            provider,
            packing_bounds,
            trace_floor,
            {},
            {},
            resolved.identities.empty() ? nullptr : &resolved,
            anchored_ptr
        );
    }
    const epcsaft_native_transport& transport =
        checked_chemical_transport(provider_capsule);
    const ProviderContext provider(transport);
    const ProviderActiveParameterSet resolved = provider.active_parameters();
    if (active_parameters_requested != !resolved.identities.empty()) {
        throw py::value_error(
            "Provider capsule active-parameter identity does not match the request"
        );
    }
    return solve_provider_source(
        provider_capsule,
        spec,
        py::cast<py::dict>(source_standard_state),
        packing_bounds,
        trace_floor,
        sensitivities_requested,
        resolved.identities.empty() ? nullptr : &resolved,
        anchored_ptr
    );
}

py::dict evaluate_provider_local(
    const py::capsule& capsule,
    const py::dict& spec,
    const std::vector<double>& packing_bounds,
    double trace_floor,
    const std::vector<double>& variables,
    const std::vector<double>& lower,
    const std::vector<double>& upper,
    bool active_parameters_requested
) {
    const epcsaft_native_transport& transport = checked_chemical_transport(capsule);
    const ChemicalProviderMetadata metadata = chemical_provider_metadata(transport);
    ReactionSystemInput input = reaction_system_input(spec);
    const ProviderContext provider(transport);
    input.provider_fingerprint = provider.fingerprint();
    const ProviderActiveParameterSet parameters = provider.active_parameters();
    if (active_parameters_requested != !parameters.identities.empty()) {
        throw py::value_error(
            "Provider capsule active-parameter identity does not match the request"
        );
    }
    const CompiledReactionSystem compiled = compile_provider_input(
        transport, metadata, input, packing_bounds
    );
    const std::size_t parameter_count = parameters.identities.size();
    ChemicalSolveResult evaluation = evaluate_provider_reaction_local(
        compiled, provider, input.temperature_k, input.pressure_pa,
        packing_bounds[0], packing_bounds[1],
        transport.total_ion_mole_fraction_max, trace_floor,
        variables, lower, upper,
        std::vector<double>(compiled.reaction_matrix.rows * parameter_count, 0.0),
        parameter_count == 0 ? nullptr : &parameters
    );
    set_provider_parameter_status(evaluation, parameter_count);
    return chemical_result(evaluation);
}

py::dict reactive_finite_phase_result(const ReactiveFinitePhaseResult& evaluation) {
    py::dict result;
    const auto optional_float = [](double value) -> py::object {
        return std::isfinite(value) ? py::cast(value) : py::none();
    };
    result["accepted"] = evaluation.accepted;
    result["solver_status"] = evaluation.solver_status;
    result["callback_error"] = evaluation.callback_error;
    result["callback_failure_kind"] = evaluation.callback_failure_kind;
    result["failure_reason"] = evaluation.failure_reason;
    result["numerical_status"] = evaluation.numerical_status;
    result["physical_status"] = evaluation.physical_status;
    result["provider_domain_status"] = evaluation.provider_domain_status;
    result["resolved_policy"] = evaluation.resolved_policy;
    result["phase_amounts_mol"] = evaluation.phase_amounts_mol;
    result["phase_volumes_m3"] = evaluation.phase_volumes_m3;
    result["phase_pressures_pa"] = evaluation.phase_pressures_pa;
    result["phase_packing_fractions"] = evaluation.phase_packing_fractions;
    result["phase_chemical_potentials_over_rt"] =
        evaluation.phase_chemical_potentials_over_rt;
    result["phase_log_fugacities"] = evaluation.phase_log_fugacities;
    result["objective_over_rt"] = optional_float(evaluation.objective_over_rt);
    result["balance_inf_norm"] = optional_float(evaluation.balance_inf_norm);
    result["pressure_relative_inf_norm"] =
        optional_float(evaluation.pressure_relative_inf_norm);
    result["reaction_affinity_inf_norm"] =
        optional_float(evaluation.reaction_affinity_inf_norm);
    result["transfer_chemical_potential_inf_norm"] =
        optional_float(evaluation.transfer_chemical_potential_inf_norm);
    result["kkt_stationarity_inf_norm"] =
        optional_float(evaluation.kkt_stationarity_inf_norm);
    result["phase_distinctness"] = optional_float(evaluation.phase_distinctness);
    result["reduced_hessian_min_eigenvalue"] =
        optional_float(evaluation.reduced_hessian_min_eigenvalue);
    result["declared_direction_rank"] = evaluation.declared_direction_rank;
    result["admissible_nullity"] = evaluation.admissible_nullity;
    result["optimizer_iterations"] = evaluation.optimizer_iterations;
    result["provider_evaluation_count"] = evaluation.provider_evaluation_count;
    result["derivative_audits_performed"] = evaluation.derivative_audits_performed;
    result["derivative_audit_status"] = evaluation.derivative_audit_status;
    result["objective_gradient_audit_max_relative_error"] = optional_float(
        evaluation.objective_gradient_audit_max_relative_error
    );
    result["constraint_jacobian_audit_max_relative_error"] = optional_float(
        evaluation.constraint_jacobian_audit_max_relative_error
    );
    result["lagrangian_hessian_audit_max_relative_error"] = optional_float(
        evaluation.lagrangian_hessian_audit_max_relative_error
    );
    result["parameter_order"] = evaluation.parameter_order;
    result["sensitivity_status"] = evaluation.sensitivity_status;
    result["sensitivity_failure_reason"] = evaluation.sensitivity_failure_reason;
    result["amount_derivatives"] = evaluation.amount_derivatives;
    result["volume_derivatives"] = evaluation.volume_derivatives;
    result["pressure_derivatives"] = evaluation.pressure_derivatives;
    result["chemical_potential_derivatives"] =
        evaluation.chemical_potential_derivatives;
    result["log_fugacity_derivatives"] = evaluation.log_fugacity_derivatives;
    return result;
}

py::dict solve_reactive_finite_phases(
    const py::capsule& capsule,
    const py::dict& spec,
    const py::object& source_standard_state,
    const std::vector<std::array<double, 2>>& packing_fraction_intervals,
    const std::vector<py::dict>& start_payloads,
    double trace_floor,
    bool active_parameters_requested,
    bool derivative_audits_requested,
    const py::object& phase_direction_system,
    std::vector<std::vector<std::size_t>> phase_support_indices
) {
    const epcsaft_native_transport& transport = checked_reactive_finite_transport(capsule);
    const ChemicalProviderMetadata metadata = chemical_provider_metadata(transport);
    ReactionSystemInput input = reaction_system_input(spec);
    validate_provider_identity(metadata, input);
    if (input.temperature_k < transport.source_temperature_min_k
        || input.temperature_k > transport.source_temperature_max_k) {
        throw py::value_error("temperature is outside the Provider source domain");
    }
    const ProviderContext provider(transport);
    input.provider_fingerprint = provider.fingerprint();
    const ProviderActiveParameterSet parameters = provider.active_parameters();
    if (active_parameters_requested != !parameters.identities.empty()) {
        throw py::value_error(
            "Provider capsule active-parameter identity does not match the request"
        );
    }
    CompiledReactionSystem compiled;
    const bool configured_directions = !phase_direction_system.is_none();
    if (source_standard_state.is_none()) {
        compiled = compile_reaction_system(input, !configured_directions);
    } else {
        const py::dict standard = py::cast<py::dict>(source_standard_state);
        const std::string standard_id = py::cast<std::string>(standard["id"]);
        const std::vector<double> log_activity_scale_factors =
            py::cast<std::vector<double>>(standard["log_activity_scale_factors"]);
        const double reference_pressure_pa = py::cast<double>(
            standard["reference_pressure_pa"]
        );
        if (standard_id.empty() || !std::isfinite(reference_pressure_pa)
            || reference_pressure_pa <= 0.0) {
            throw py::value_error("source standard-state declaration is invalid");
        }
        for (const EquilibriumConstantRecord& record : input.equilibrium_constant_records) {
            if (record.reference_id != standard_id
                || record.conversion_id
                    != "source-standard-state-to-provider-neutral-reference"
                || record.pressure_pa != reference_pressure_pa) {
                throw py::value_error(
                    "source equilibrium-constant declaration is incompatible"
                );
            }
        }
        const SourceReferenceTransferEvidence reference =
            source_reference_transfer_evidence(standard, "provider_transfer");
        if (!parameters.identities.empty()
            && reference.active_parameter_identities != parameters.identities) {
            throw py::value_error(
                "source-reference active parameters do not match the Provider"
            );
        }
        const SourceStandardStateResult transformed =
            transform_source_standard_state(
                input.reaction_matrix,
                input.ln_k,
                log_activity_scale_factors,
                input.species_ids,
                input.charges,
                input.provider_fingerprint,
                input.temperature_k,
                input.pressure_pa,
                reference
            );
        input.ln_k = transformed.ln_k_provider_basis;
        for (EquilibriumConstantRecord& record : input.equilibrium_constant_records) {
            record.reference_id = "provider-helmholtz-coordinate-basis";
            record.conversion_id = "already-provider-basis";
            record.pressure_pa = input.pressure_pa;
        }
        compiled = compile_reaction_system(input, !configured_directions);
        const std::size_t parameter_count = parameters.identities.size();
        if (parameter_count != 0) {
            if (transformed.ln_k_parameter_derivatives.size()
                != input.reaction_matrix.rows * parameter_count) {
                throw py::value_error(
                    "source-reference parameter derivatives are unavailable"
                );
            }
            compiled.g_ref_parameter_derivatives.assign(
                compiled.species_count * parameter_count, 0.0
            );
            for (std::size_t parameter = 0;
                 parameter < parameter_count;
                 ++parameter) {
                ReactionSystemInput derivative_input = input;
                for (std::size_t reaction = 0;
                     reaction < input.reaction_matrix.rows;
                     ++reaction) {
                    derivative_input.ln_k[reaction] =
                        transformed.ln_k_parameter_derivatives[
                            reaction * parameter_count + parameter
                        ];
                }
                const CompiledReactionSystem derivative_system =
                    compile_reaction_system(
                        derivative_input, !configured_directions
                    );
                if (derivative_system.species_count != compiled.species_count
                    || derivative_system.retained_species_indices
                        != compiled.retained_species_indices) {
                    throw py::value_error(
                        "source-reference derivative species face changed"
                    );
                }
                for (std::size_t species = 0;
                     species < compiled.species_count;
                     ++species) {
                    compiled.g_ref_parameter_derivatives[
                        species * parameter_count + parameter
                    ] = derivative_system.g_ref[species];
                }
            }
        }
    }
    std::vector<ReactiveFinitePhaseStart> starts;
    starts.reserve(start_payloads.size());
    for (const py::dict& payload : start_payloads) {
        starts.push_back({
            py::cast<std::vector<double>>(payload["amounts_mol"]),
            py::cast<double>(payload["molar_volume_m3_per_mol"]),
        });
    }
    if (phase_support_indices.empty()) {
        const std::vector<std::size_t> full_support = [&input]() {
            std::vector<std::size_t> indices(input.species_ids.size());
            std::iota(indices.begin(), indices.end(), 0);
            return indices;
        }();
        phase_support_indices.assign(starts.size(), full_support);
    }
    std::optional<PhaseDirectionCompilation> phase_directions;
    if (!phase_direction_system.is_none()) {
        const py::dict declaration = py::cast<py::dict>(phase_direction_system);
        std::vector<std::vector<double>> reference_rows;
        reference_rows.reserve(starts.size());
        for (const ReactiveFinitePhaseStart& start : starts) {
            reference_rows.push_back(start.amounts_mol);
        }
        for (std::size_t species = 0;
             species < input.species_ids.size();
             ++species) {
            double declared_total = 0.0;
            double scale = std::abs(input.feed_amounts[species]);
            for (const ReactiveFinitePhaseStart& start : starts) {
                if (start.amounts_mol.size() != input.species_ids.size()) {
                    throw py::value_error(
                        "phase reference amounts must match the species order"
                    );
                }
                declared_total += start.amounts_mol[species];
                scale += std::abs(start.amounts_mol[species]);
            }
            if (std::abs(declared_total - input.feed_amounts[species])
                > 1.0e-10 * std::max(1.0, scale)) {
                throw py::value_error(
                    "phase reference amounts must sum to the declared feed"
                );
            }
        }
        phase_directions = compile_phase_direction_system(
            py::cast<std::vector<std::string>>(declaration["phase_ids"]),
            py::cast<std::vector<std::string>>(declaration["species_ids"]),
            py::cast<std::vector<double>>(declaration["charges"]),
            phase_direction_matrix(
                py::cast<std::vector<std::vector<double>>>(
                    declaration["balance_matrix"]
                ),
                input.species_ids.size(),
                "balance_matrix"
            ),
            phase_direction_matrix(
                reference_rows,
                input.species_ids.size(),
                "phase_reference_amounts"
            ),
            phase_direction_matrix(
                py::cast<std::vector<std::vector<double>>>(
                    declaration["direction_rows"]
                ),
                starts.size() * input.species_ids.size(),
                "direction_rows"
            ),
            py::cast<std::vector<double>>(declaration["driving_constants"]),
            kPhaseDirectionRankTolerance,
            kPhaseDirectionConservationTolerance
        );
    }
    return reactive_finite_phase_result(
        solve_provider_reactive_finite_phases(
            compiled,
            provider,
            input.temperature_k,
            input.pressure_pa,
            packing_fraction_intervals,
            starts,
            trace_floor,
            parameters.identities.empty() ? nullptr : &parameters,
            500,
            derivative_audits_requested,
            phase_directions.has_value() ? &*phase_directions : nullptr,
            phase_support_indices
        )
    );
}

std::vector<std::vector<double>> matrix_rows(const DenseMatrix& matrix) {
    std::vector<std::vector<double>> result(
        matrix.rows,
        std::vector<double>(matrix.columns)
    );
    for (std::size_t row = 0; row < matrix.rows; ++row) {
        for (std::size_t column = 0; column < matrix.columns; ++column) {
            result[row][column] = matrix.values[row * matrix.columns + column];
        }
    }
    return result;
}

py::dict compile_phase_direction_system_binding(
    const std::vector<std::string>& phase_ids,
    const std::vector<std::string>& species_ids,
    const std::vector<double>& charges,
    const std::vector<std::vector<double>>& balance_rows,
    const std::vector<std::vector<double>>& phase_reference_rows,
    const std::vector<std::vector<double>>& direction_rows,
    const std::vector<double>& driving_constants,
    double rank_tolerance,
    double conservation_tolerance
) {
    const std::size_t dimension = phase_ids.size() * species_ids.size();
    const PhaseDirectionCompilation compilation =
        compile_phase_direction_system(
            phase_ids,
            species_ids,
            charges,
            phase_direction_matrix(
                balance_rows,
                species_ids.size(),
                "balance_matrix"
            ),
            phase_direction_matrix(
                phase_reference_rows,
                species_ids.size(),
                "phase_reference_amounts"
            ),
            phase_direction_matrix(
                direction_rows,
                dimension,
                "direction_rows"
            ),
            driving_constants,
            rank_tolerance,
            conservation_tolerance
        );
    py::dict result;
    result["phase_ids"] = compilation.phase_ids;
    result["species_ids"] = compilation.species_ids;
    result["charges"] = compilation.charges;
    result["coordinate_count"] = dimension;
    py::list coordinate_order;
    for (const std::string& phase_id : compilation.phase_ids) {
        for (const std::string& species_id : compilation.species_ids) {
            coordinate_order.append(py::make_tuple(phase_id, species_id));
        }
    }
    result["coordinate_order"] = std::move(coordinate_order);
    result["direction_matrix"] = matrix_rows(compilation.direction_matrix);
    result["directions"] = matrix_rows(compilation.direction_matrix);
    result["invariant_rows"] = matrix_rows(compilation.invariant_matrix);
    result["invariant_matrix"] = matrix_rows(compilation.invariant_matrix);
    result["phase_reference_amounts"] = compilation.phase_reference_amounts;
    result["driving_constants"] = compilation.driving_constants;
    result["invariant_totals"] = compilation.invariant_totals;
    result["reference_potentials"] = compilation.reference_potentials;
    result["g"] = compilation.reference_potentials;
    result["global_balance_residuals"] = compilation.global_balance_residuals;
    result["phase_charge_residuals"] = compilation.phase_charge_residuals;
    result["reference_phase_charge_residuals"] =
        compilation.reference_phase_charge_residuals;
    result["direction_rank"] = compilation.direction_rank;
    result["rank"] = compilation.direction_rank;
    result["invariant_rank"] = compilation.invariant_rank;
    result["direction_count"] = compilation.direction_matrix.rows;
    result["invariant_count"] = compilation.invariant_matrix.rows;
    result["direction_scale"] = compilation.direction_scale;
    result["rank_tolerance"] = compilation.rank_tolerance;
    result["conservation_tolerance"] = compilation.conservation_tolerance;
    result["invariant_orthonormal_residual_inf_norm"] =
        compilation.invariant_orthonormal_residual_inf;
    result["invariant_direction_residual_inf_norm"] =
        compilation.invariant_direction_residual_inf;
    result["potential_equation_residual_inf_norm"] =
        compilation.potential_equation_residual_inf;
    result["minimum_norm_potential_norm"] =
        compilation.minimum_norm_potential_norm;
    result["global_balance_residual_inf_norm"] =
        compilation.global_balance_residual_inf;
    result["phase_charge_residual_inf_norm"] =
        compilation.phase_charge_residual_inf;
    result["reference_neutrality_residual_inf_norm"] =
        compilation.reference_neutrality_residual_inf;
    return result;
}

}  // namespace

void bind_chemical_equilibrium(py::module_& module) {
    module.def(
        "_chemical_equilibrium",
        &solve_chemical_equilibrium,
        py::arg("capsule"),
        py::arg("spec"),
        py::arg("source_standard_state"),
        py::arg("packing_fraction_bounds"),
        py::arg("trace_floor"),
        py::arg("sensitivities_requested") = false,
        py::arg("active_parameters_requested") = false
    );
    module.def(
        "_evaluate_provider_reaction_local",
        &evaluate_provider_local,
        py::arg("capsule"),
        py::arg("spec"),
        py::arg("packing_fraction_bounds"),
        py::arg("trace_floor"),
        py::arg("variables"),
        py::arg("lower"),
        py::arg("upper"),
        py::arg("active_parameters_requested") = false
    );
    module.def(
        "_reactive_finite_phases",
        &solve_reactive_finite_phases,
        py::arg("capsule"),
        py::arg("spec"),
        py::arg("source_standard_state"),
        py::arg("packing_fraction_intervals"),
        py::arg("starts"),
        py::arg("trace_floor"),
        py::arg("active_parameters_requested") = false,
        py::arg("derivative_audits_requested") = false,
        py::arg("phase_direction_system") = py::none(),
        py::arg("phase_support_indices") =
            std::vector<std::vector<std::size_t>>{}
    );
    module.def(
        "_compile_phase_direction_system",
        &compile_phase_direction_system_binding,
        py::arg("phase_ids"),
        py::arg("species_ids"),
        py::arg("charges"),
        py::arg("balance_matrix"),
        py::arg("phase_reference_amounts"),
        py::arg("direction_rows"),
        py::arg("driving_constants"),
        py::arg("rank_tolerance") = kPhaseDirectionRankTolerance,
        py::arg("conservation_tolerance") = kPhaseDirectionConservationTolerance
    );
    // Keep the private seam discoverable under the concise plural spelling as
    // well; both names call the same compiler and have no separate behavior.
    module.def(
        "_compile_phase_directions",
        &compile_phase_direction_system_binding,
        py::arg("phase_ids"),
        py::arg("species_ids"),
        py::arg("charges"),
        py::arg("balance_matrix"),
        py::arg("phase_reference_amounts"),
        py::arg("directions"),
        py::arg("driving_constants"),
        py::arg("rank_tolerance") = kPhaseDirectionRankTolerance,
        py::arg("conservation_tolerance") = kPhaseDirectionConservationTolerance
    );
}

}  // namespace epcsaft_engine::equilibrium
