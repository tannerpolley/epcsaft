#include "chemical_equilibrium.hpp"
#include "provider.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numeric>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

namespace epcsaft_engine::equilibrium {
namespace {

struct ChemicalProviderMetadata {
    std::vector<std::string> component_ids;
    std::vector<int> charges;
};

ChemicalProviderMetadata chemical_provider_metadata(const ProviderContext& provider) {
    if (provider.component_count() < 2
        || (!provider.supports_mixture_phase()
            && !provider.supports_electrolyte_phase())) {
        throw py::value_error(
            "Provider model does not support chemical-equilibrium phase evaluation"
        );
    }
    ChemicalProviderMetadata result;
    result.component_ids = provider.component_ids();
    result.charges = provider.component_charges();
    for (const std::string& component_id : result.component_ids) {
        if (component_id.empty()) {
            throw py::value_error("Provider component identity is incomplete");
        }
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

void validate_provider_identity(
    const ChemicalProviderMetadata& metadata,
    const ReactionSystemInput& input
) {
    if (input.species_ids != metadata.component_ids) {
        throw py::value_error("Provider component order does not match the reaction system");
    }
    if (input.charges != metadata.charges) {
        throw py::value_error("Provider charges do not match the reaction system");
    }
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
    result["minimum_amount_mol"] = optional_float(evaluation.minimum_amount_mol);
    result["minimum_bound_margin"] = optional_float(evaluation.minimum_bound_margin);
    result["reduced_hessian_min_eigenvalue"] =
        optional_float(evaluation.reduced_hessian_min_eigenvalue);
    result["reduced_hessian_scale"] =
        optional_float(evaluation.reduced_hessian_scale);
    result["reduced_hessian_threshold"] =
        optional_float(evaluation.reduced_hessian_threshold);
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
    result["local_variables"] = evaluation.local_variables;
    result["local_lower_bounds"] = evaluation.local_lower_bounds;
    result["local_upper_bounds"] = evaluation.local_upper_bounds;
    result["input_order"] = evaluation.input_order;
    result["kkt_input_jacobian"] = evaluation.kkt_input_jacobian;
    result["input_amount_derivatives"] = evaluation.input_amount_derivatives;
    result["input_volume_derivatives"] = evaluation.input_volume_derivatives;
    result["input_pressure_derivatives"] = evaluation.input_pressure_derivatives;
    result["input_chemical_potential_derivatives"] =
        evaluation.input_chemical_potential_derivatives;
    result["input_log_fugacity_derivatives"] =
        evaluation.input_log_fugacity_derivatives;
    result["kkt_root_values"] = evaluation.kkt_root_values;
    result["kkt_root_jacobian"] = evaluation.kkt_root_jacobian;
    result["kkt_dimension"] = evaluation.kkt_dimension;
    result["kkt_rank"] = evaluation.kkt_rank;
    result["condition_number_inf"] = optional_float(evaluation.condition_number_inf);
    result["output_values"] = evaluation.output_values;
    result["output_input_derivatives"] = evaluation.output_input_derivatives;
    result["output_local_derivatives"] = evaluation.output_local_derivatives;
    result["output_parameter_partials"] = evaluation.output_parameter_partials;
    return result;
}

CompiledReactionSystem expand_reactive_species_face(
    CompiledReactionSystem compiled,
    const ReactionSystemInput& input
) {
    if (compiled.species_count == input.species_ids.size()) {
        return compiled;
    }
    const std::size_t full_count = input.species_ids.size();
    const std::size_t reduced_count = compiled.species_count;
    const auto expand_matrix = [full_count, reduced_count, &compiled](
                                   const DenseMatrix& reduced
                               ) {
        DenseMatrix full{
            reduced.rows,
            full_count,
            std::vector<double>(reduced.rows * full_count, 0.0),
        };
        for (std::size_t row = 0; row < reduced.rows; ++row) {
            for (std::size_t retained = 0; retained < reduced_count; ++retained) {
                full(row, compiled.retained_species_indices[retained]) =
                    reduced(row, retained);
            }
        }
        return full;
    };
    const auto expand_species = [full_count, reduced_count, &compiled](
                                    const std::vector<double>& reduced,
                                    std::size_t width
                                ) {
        std::vector<double> full(full_count * width, 0.0);
        for (std::size_t retained = 0; retained < reduced_count; ++retained) {
            const std::size_t species = compiled.retained_species_indices[retained];
            std::copy_n(
                reduced.begin() + static_cast<std::ptrdiff_t>(retained * width),
                width,
                full.begin() + static_cast<std::ptrdiff_t>(species * width)
            );
        }
        return full;
    };
    compiled.balance_matrix = expand_matrix(compiled.balance_matrix);
    compiled.reaction_matrix = expand_matrix(compiled.reaction_matrix);
    compiled.g_ref = expand_species(compiled.g_ref, 1);
    if (!compiled.g_ref_pressure_derivatives_per_pa.empty()) {
        compiled.g_ref_pressure_derivatives_per_pa = expand_species(
            compiled.g_ref_pressure_derivatives_per_pa, 1
        );
    }
    if (!compiled.g_ref_parameter_derivatives.empty()) {
        compiled.g_ref_parameter_derivatives = expand_species(
            compiled.g_ref_parameter_derivatives,
            compiled.g_ref_parameter_derivatives.size() / reduced_count
        );
    }
    compiled.species_count = full_count;
    compiled.charges = input.charges;
    compiled.molar_masses_kg_per_mol = input.molar_masses_kg_per_mol;
    compiled.feed_amounts = input.feed_amounts;
    return compiled;
}

py::dict solve_reactive_finite_phases(
    const epcsaft::NativeModel& model,
    const py::dict& spec,
    const py::object& source_standard_state,
    const std::vector<std::array<double, 2>>& packing_fraction_intervals,
    const std::vector<py::dict>& start_payloads,
    double trace_floor,
    bool active_parameters_requested,
    bool derivative_audits_requested,
    const py::object& phase_direction_system,
    std::vector<std::vector<std::size_t>> phase_support_indices,
    const std::vector<py::dict>& output_payloads,
    const std::vector<double>& local_variables,
    const std::vector<double>& local_lower_bounds,
    const std::vector<double>& local_upper_bounds,
    bool ideal_gas
) {
    const ProviderContext provider(model);
    if (!ideal_gas && !provider.supports_mixture_phase()) {
        throw py::value_error(
            "Provider model does not support reactive finite-phase evaluation"
        );
    }
    const ChemicalProviderMetadata metadata = chemical_provider_metadata(provider);
    ReactionSystemInput input = reaction_system_input(spec);
    if (!ideal_gas) {
        validate_provider_identity(metadata, input);
    }
    if (!ideal_gas
        && (input.temperature_k < provider.metadata().temperature_min_k
            || input.temperature_k > provider.metadata().temperature_max_k)) {
        throw py::value_error("temperature is outside the Provider source domain");
    }
    input.provider_fingerprint = provider.fingerprint();
    const ProviderActiveParameterSet parameters = provider.active_parameters();
    if (active_parameters_requested != !parameters.identities.empty()) {
        throw py::value_error(
            "Provider active-parameter identity does not match the request"
        );
    }
    CompiledReactionSystem compiled;
    double source_representation_residual = 0.0;
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
        source_representation_residual = transformed.representation_residual_inf_norm;
        input.ln_k = transformed.ln_k_provider_basis;
        for (EquilibriumConstantRecord& record : input.equilibrium_constant_records) {
            record.reference_id = "provider-helmholtz-coordinate-basis";
            record.conversion_id = "already-provider-basis";
            record.pressure_pa = input.pressure_pa;
        }
        compiled = compile_reaction_system(input, !configured_directions);
        const std::size_t parameter_count = parameters.identities.size();
        if (!transformed.ln_k_pressure_derivatives_per_pa.empty()) {
            ReactionSystemInput derivative_input = input;
            derivative_input.ln_k = transformed.ln_k_pressure_derivatives_per_pa;
            const CompiledReactionSystem derivative_system = compile_reaction_system(
                derivative_input, !configured_directions
            );
            if (derivative_system.species_count != compiled.species_count
                || derivative_system.retained_species_indices
                    != compiled.retained_species_indices) {
                throw py::value_error(
                    "source-reference pressure derivative species face changed"
                );
            }
            compiled.g_ref_pressure_derivatives_per_pa = derivative_system.g_ref;
        }
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
    compiled = expand_reactive_species_face(std::move(compiled), input);
    compiled.ideal_gas = ideal_gas;
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
    for (std::vector<std::size_t>& support : phase_support_indices) {
        support.erase(
            std::remove_if(
                support.begin(),
                support.end(),
                [&compiled](std::size_t species) {
                    return std::binary_search(
                        compiled.removed_species_indices.begin(),
                        compiled.removed_species_indices.end(),
                        species
                    );
                }
            ),
            support.end()
        );
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
            const double relative_tolerance = std::min(
                1.0e-12,
                256.0 * std::numeric_limits<double>::epsilon()
                    * static_cast<double>(std::max<std::size_t>(1, starts.size()))
            );
            const double tolerance = relative_tolerance
                * std::max(scale, std::numeric_limits<double>::min());
            if (std::abs(declared_total - input.feed_amounts[species])
                > tolerance) {
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
        for (std::size_t direction = 0;
             direction < phase_directions->direction_matrix.rows;
             ++direction) {
            double change = 0.0;
            double scale = 0.0;
            for (std::size_t phase = 0; phase < starts.size(); ++phase) {
                for (std::size_t species = 0;
                     species < input.species_ids.size();
                     ++species) {
                    const double contribution = phase_directions->direction_matrix(
                        direction,
                        phase * input.species_ids.size() + species
                    ) * input.molar_masses_kg_per_mol[species];
                    change += contribution;
                    scale += std::abs(contribution);
                }
            }
            if (std::abs(change)
                > kPhaseDirectionConservationTolerance * std::max(1.0, scale)) {
                throw py::value_error(
                    "phase direction violates global molar-mass conservation"
                );
            }
        }
    }
    std::vector<ReactiveFinitePhaseOutput> outputs;
    outputs.reserve(output_payloads.size());
    for (const py::dict& payload : output_payloads) {
        outputs.push_back({
            py::cast<std::string>(payload["selector"]),
            py::cast<std::size_t>(payload["phase_index"]),
            py::cast<std::vector<double>>(payload["coefficients"]),
            py::cast<std::vector<double>>(
                payload["solvent_mass_coefficients_kg_per_mol"]
            ),
        });
    }
    const ReactiveFinitePhaseResult evaluation = local_variables.empty()
        ? solve_provider_reactive_finite_phases(
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
            phase_support_indices,
            outputs
        )
        : evaluate_provider_reactive_finite_phases_local(
            compiled,
            provider,
            input.temperature_k,
            input.pressure_pa,
            packing_fraction_intervals,
            local_variables,
            local_lower_bounds,
            local_upper_bounds,
            parameters.identities.empty() ? nullptr : &parameters,
            phase_directions.has_value() ? &*phase_directions : nullptr,
            phase_support_indices,
            outputs
        );
    py::dict result = reactive_finite_phase_result(evaluation);
    result["source_reference_representation_residual_inf_norm"] =
        source_representation_residual;
    return result;
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
        "_reactive_finite_phases",
        &solve_reactive_finite_phases,
        py::arg("model"),
        py::arg("spec"),
        py::arg("source_standard_state"),
        py::arg("packing_fraction_intervals"),
        py::arg("starts"),
        py::arg("trace_floor"),
        py::arg("active_parameters_requested") = false,
        py::arg("derivative_audits_requested") = false,
        py::arg("phase_direction_system") = py::none(),
        py::arg("phase_support_indices") =
            std::vector<std::vector<std::size_t>>{},
        py::arg("outputs") = std::vector<py::dict>{},
        py::arg("local_variables") = std::vector<double>{},
        py::arg("local_lower_bounds") = std::vector<double>{},
        py::arg("local_upper_bounds") = std::vector<double>{},
        py::arg("ideal_gas") = false
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
