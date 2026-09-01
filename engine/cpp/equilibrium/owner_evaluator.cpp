#include "owner_evaluator.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace epcsaft_engine::equilibrium {
namespace {

constexpr const char* kHomogeneousKind =
    "D-034.fixed_TP_homogeneous_reactive";
constexpr const char* kBubbleVleKind =
    "D-035.fixed_topology_reactive_bubble_vle";

interfaces::EvaluationFailure failure(
    std::string code,
    std::string reason
) {
    return {std::move(code), std::move(reason)};
}

bool finite(const std::vector<double>& values) {
    return std::all_of(values.begin(), values.end(), [](double value) {
        return std::isfinite(value);
    });
}

bool unique_nonempty(const std::vector<std::string>& values) {
    if (!std::all_of(values.begin(), values.end(), [](const std::string& value) {
        return !value.empty();
    })) {
        return false;
    }
    std::vector<std::string> sorted = values;
    std::sort(sorted.begin(), sorted.end());
    return std::adjacent_find(sorted.begin(), sorted.end()) == sorted.end();
}

std::string failure_code(const GrepeResult& result) {
    if (!result.failure_code.empty()) {
        return result.failure_code;
    }
    const auto& equilibrium = result.equilibrium;
    if (equilibrium.eos_domain_status == "failed") {
        return "eos_domain_rejection";
    }
    if (equilibrium.physical_status == "failed") {
        return "physical_rejection";
    }
    if (equilibrium.sensitivity_status == "unavailable") {
        return "singular_sensitivity";
    }
    return "numerical_convergence_failure";
}

interfaces::EvaluationFailure invalid_result(
    const std::exception& error
) {
    const std::string reason = error.what();
    const bool domain = reason.find("EOS") != std::string::npos
        || reason.find("packing") != std::string::npos
        || reason.find("domain") != std::string::npos;
    return failure(domain ? "eos_domain_rejection" : "owner_failure", reason);
}

epcsaft::NativeModel model_at(
    const epcsaft::NativeModel& model,
    const std::vector<double>& values
) {
    std::vector<epcsaft::NativeActiveCoordinate> coordinates = model.active_coordinates();
    if (coordinates.size() != values.size()) {
        throw std::invalid_argument("active parameter dimensions changed");
    }
    if (coordinates.empty()) {
        return model;
    }
    for (std::size_t index = 0; index < coordinates.size(); ++index) {
        coordinates[index].value = values[index];
    }
    return model.with_active_coordinates(coordinates);
}

std::vector<double> values_for_model(
    const epcsaft::NativeModel& model,
    const std::vector<std::string>& parameter_ids,
    const std::vector<double>& values
) {
    std::vector<double> result;
    for (const auto& coordinate : model.active_coordinates()) {
        const auto found = std::find(
            parameter_ids.begin(), parameter_ids.end(), coordinate.identity
        );
        if (found == parameter_ids.end()) {
            throw std::invalid_argument("active parameter order changed");
        }
        result.push_back(values[static_cast<std::size_t>(found - parameter_ids.begin())]);
    }
    return result;
}

std::vector<std::string> identities(
    const std::vector<epcsaft::NativeActiveCoordinate>& coordinates
) {
    return epcsaft::active_parameter_identities(coordinates);
}

void hash_field(std::uint64_t& state, const std::string& field) {
    for (const unsigned char value : field) {
        state ^= value;
        state *= 1099511628211ULL;
    }
    state ^= 0xffU;
    state *= 1099511628211ULL;
}

std::string canonical_hash(const std::vector<std::string>& fields) {
    std::uint64_t state = 14695981039346656037ULL;
    for (const auto& field : fields) {
        hash_field(state, field);
    }
    std::ostringstream rendered;
    rendered << "fnv1a64:" << std::hex << std::setw(16) << std::setfill('0') << state;
    return rendered.str();
}

std::string canonical_double(double value) {
    std::ostringstream rendered;
    rendered << std::setprecision(17) << value;
    return rendered.str();
}

void append_doubles(std::vector<std::string>& fields, const std::vector<double>& values) {
    for (const double value : values) {
        fields.push_back(canonical_double(value));
    }
}

void append_indices(std::vector<std::string>& fields, const std::vector<std::size_t>& values) {
    for (const std::size_t value : values) {
        fields.push_back(std::to_string(value));
    }
}

void append_matrix(std::vector<std::string>& fields, const DenseMatrix& matrix) {
    fields.push_back(std::to_string(matrix.rows));
    fields.push_back(std::to_string(matrix.columns));
    append_doubles(fields, matrix.values);
}

std::string canonical_parameter_contract(
    const epcsaft::NativeModel& model,
    const std::vector<EvaluatorParameter>& parameters
) {
    std::vector<std::string> fields{"parameter-contract-v1", model.parameter_fingerprint()};
    for (const auto& parameter : parameters) {
        fields.insert(fields.end(), {parameter.identity, parameter.unit, parameter.basis});
    }
    return canonical_hash(fields);
}

std::string canonical_evaluator_fingerprint(
    const epcsaft::NativeModel& model,
    const GrepeRequest& request,
    const EquilibriumEvaluatorSidecar& sidecar
) {
    std::vector<std::string> fields{
        "equilibrium-owner-v1",
        model.parameter_fingerprint(),
        request.identity,
        sidecar.owner_kind,
        sidecar.source_spec_digest,
        sidecar.anchor_identity,
        sidecar.anchor_receipt_identity,
        sidecar.temperature_role,
        sidecar.pressure_role,
        sidecar.derivative_provenance,
        sidecar.topology_policy,
        sidecar.branch_policy,
        sidecar.conditioning_policy,
    };
    fields.insert(fields.end(), sidecar.active_parameter_identities.begin(),
        sidecar.active_parameter_identities.end());
    fields.insert(fields.end(), sidecar.active_parameter_units.begin(),
        sidecar.active_parameter_units.end());
    fields.insert(fields.end(), sidecar.active_parameter_bases.begin(),
        sidecar.active_parameter_bases.end());
    for (const auto& phase : sidecar.phases) {
        fields.insert(fields.end(), {
            phase.identity, phase.role, phase.support,
            phase.finite ? "finite" : "incipient",
            phase.ideal_gas ? "ideal-gas" : "eos"
        });
        append_indices(fields, phase.support_indices);
    }
    fields.insert(fields.end(), model.component_ids().begin(), model.component_ids().end());
    for (const int charge : model.charge_numbers()) {
        fields.push_back(std::to_string(charge));
    }
    fields.insert(fields.end(), sidecar.reaction_identities.begin(),
        sidecar.reaction_identities.end());
    fields.insert(fields.end(), sidecar.reaction_phase_identities.begin(),
        sidecar.reaction_phase_identities.end());
    fields.insert(fields.end(), request.reaction.species_ids.begin(),
        request.reaction.species_ids.end());
    for (const int charge : request.reaction.charges) {
        fields.push_back(std::to_string(charge));
    }
    fields.push_back(request.reaction.eos_fingerprint);
    append_doubles(fields, request.reaction.molar_masses_kg_per_mol);
    append_matrix(fields, request.reaction.balance_matrix);
    append_matrix(fields, request.reaction.reaction_matrix);
    append_doubles(fields, request.reaction.feed_amounts);
    append_doubles(fields, request.reaction.ln_k);
    append_doubles(fields, request.reaction.ln_k_parameter_derivatives);
    for (std::size_t index = 0; index < sidecar.output_row_identities.size(); ++index) {
        fields.insert(fields.end(), {
            sidecar.output_row_identities[index], sidecar.output_row_units[index],
            sidecar.output_row_bases[index]
        });
        const auto& output = request.outputs[index];
        fields.insert(fields.end(), {output.selector, std::to_string(output.phase_index)});
        append_doubles(fields, output.coefficients);
        append_doubles(fields, output.solvent_mass_coefficients_kg_per_mol);
    }
    fields.push_back(canonical_double(request.temperature.value));
    fields.push_back(canonical_double(request.pressure.value));
    append_doubles(fields, {
        request.temperature.bounds[0], request.temperature.bounds[1],
        request.pressure.bounds[0], request.pressure.bounds[1]
    });
    append_doubles(fields, request.temperature.starts);
    append_doubles(fields, request.pressure.starts);
    fields.push_back(canonical_double(request.continuation_pressure_pa));
    fields.push_back(canonical_double(request.trace_floor));
    fields.push_back(request.declared_starts_only ? "declared-starts-only" : "generated-starts");
    for (const auto& anchor : request.continuation_phases) {
        fields.push_back(anchor.identity);
        append_doubles(fields, anchor.mole_fractions);
        fields.push_back(canonical_double(anchor.molar_volume_m3_per_mol));
    }
    for (const auto& phase : request.phases) {
        fields.insert(fields.end(), {
            phase.model_reference_id,
            phase.ideal_gas ? "ideal-gas" : "eos",
            phase.finite ? "finite" : "incipient",
        });
        append_doubles(fields, {
            phase.packing_fraction_interval[0], phase.packing_fraction_interval[1]
        });
        if (phase.start.has_value()) {
            append_doubles(fields, phase.start->amounts_mol);
            fields.push_back(canonical_double(phase.start->molar_volume_m3_per_mol));
        }
    }
    return canonical_hash(fields);
}

EquilibriumEvaluatorSidecar make_sidecar(
    const epcsaft::NativeModel& model,
    const GrepeRequest& request,
    EquilibriumEvaluatorSetup setup,
    bool bubble_vle
) {
    const auto active = identities(model.active_coordinates());
    if (setup.active_parameters.size() != active.size()) {
        throw std::invalid_argument("active parameter metadata does not match model");
    }
    for (std::size_t index = 0; index < active.size(); ++index) {
        if (setup.active_parameters[index].identity != active[index]
            || setup.active_parameters[index].unit.empty()
            || setup.active_parameters[index].basis.empty()) {
            throw std::invalid_argument("active parameter metadata order changed");
        }
    }
    if (setup.output_rows.size() != request.outputs.size()
        || setup.output_rows.empty()) {
        throw std::invalid_argument("output metadata does not match request");
    }
    for (const auto& row : setup.output_rows) {
        if (row.identity.empty() || row.unit.empty() || row.basis.empty()) {
            throw std::invalid_argument("output metadata is incomplete");
        }
    }
    if (setup.parameter_baseline_digest.empty()
        || setup.source_spec_digest.empty()
        || setup.anchor_identity.empty()
        || setup.anchor_receipt_identity.empty()) {
        throw std::invalid_argument("equilibrium evaluator identity metadata is incomplete");
    }
    std::vector<std::string> output_identities;
    output_identities.reserve(setup.output_rows.size());
    for (const auto& row : setup.output_rows) {
        output_identities.push_back(row.identity);
    }
    if (!unique_nonempty(active) || !unique_nonempty(output_identities)) {
        throw std::invalid_argument("evaluator identities are invalid");
    }
    EquilibriumEvaluatorSidecar result;
    result.owner_kind = setup.owner_kind.empty()
        ? (bubble_vle ? kBubbleVleKind : kHomogeneousKind)
        : setup.owner_kind;
    result.parameter_baseline_digest = std::move(setup.parameter_baseline_digest);
    result.source_spec_digest = std::move(setup.source_spec_digest);
    result.anchor_identity = std::move(setup.anchor_identity);
    result.anchor_receipt_identity = std::move(setup.anchor_receipt_identity);
    result.active_parameter_identities = active;
    for (const auto& parameter : setup.active_parameters) {
        result.active_parameter_units.push_back(parameter.unit);
        result.active_parameter_bases.push_back(parameter.basis);
    }
    for (const auto& row : setup.output_rows) {
        result.output_row_identities.push_back(row.identity);
        result.output_row_units.push_back(row.unit);
        result.output_row_bases.push_back(row.basis);
    }
    result.derivative_provenance = setup.derivative_provenance.empty()
        ? "equilibrium-solved-state-exact-total-derivative"
        : std::move(setup.derivative_provenance);
    result.phases.reserve(request.phases.size());
    for (const auto& phase : request.phases) {
        result.phases.push_back({
            phase.identity,
            phase.role == PhaseRole::Liquid ? "liquid" : "vapor",
            phase.support_name,
            phase.support_indices,
            phase.finite,
            phase.ideal_gas,
        });
    }
    result.temperature_role = "fixed";
    result.pressure_role = bubble_vle ? "solved" : "fixed";
    for (const auto& record : request.reaction.equilibrium_constant_records) {
        result.reaction_identities.push_back(record.reaction_id);
    }
    result.reaction_phase_identities = request.reaction_phase_ids;
    result.topology_policy = setup.topology_policy.empty()
        ? "declared-topology-only"
        : std::move(setup.topology_policy);
    result.branch_policy = setup.branch_policy.empty()
        ? "immutable-declared-local-branch"
        : std::move(setup.branch_policy);
    result.conditioning_policy = setup.conditioning_policy.empty()
        ? "typed-refusal-on-singular-or-ill-conditioned-sensitivity"
        : std::move(setup.conditioning_policy);
    result.parameter_contract_fingerprint = canonical_parameter_contract(
        model, setup.active_parameters
    );
    result.evaluator_fingerprint = canonical_evaluator_fingerprint(
        model, request, result
    );
    return result;
}

interfaces::EvaluationResult map_result(
    const GrepeResult& result,
    const EquilibriumEvaluatorSidecar& sidecar,
    bool exact_jacobian
) {
    const auto& equilibrium = result.equilibrium;
    if (!equilibrium.accepted) {
        return failure(failure_code(result), equilibrium.failure_reason.empty()
            ? "declared equilibrium owner rejected the trial"
            : equilibrium.failure_reason);
    }
    const std::size_t output_count = sidecar.output_row_identities.size();
    const std::size_t parameter_count = sidecar.active_parameter_identities.size();
    if (equilibrium.parameter_order != sidecar.active_parameter_identities
        || equilibrium.output_values.size() != output_count
        || !finite(equilibrium.output_values)) {
        return failure("owner_contract_failure", "equilibrium output contract changed");
    }
    if (exact_jacobian) {
        const std::size_t width = equilibrium.input_order.size();
        if (equilibrium.sensitivity_status != "available"
            || equilibrium.output_input_derivatives.size() != output_count * width
            || !finite(equilibrium.output_input_derivatives)) {
            return failure(
                equilibrium.sensitivity_status == "unavailable"
                    ? "singular_sensitivity" : "owner_contract_failure",
                equilibrium.sensitivity_failure_reason.empty()
                    ? "exact equilibrium total derivatives are unavailable"
                    : equilibrium.sensitivity_failure_reason
            );
        }
    }
    interfaces::Evaluation evaluation;
    evaluation.predictions = equilibrium.output_values;
    if (!exact_jacobian) {
        return evaluation;
    }
    evaluation.jacobian.assign(output_count * parameter_count, 0.0);
    const std::size_t width = equilibrium.input_order.size();
    if (equilibrium.output_input_derivatives.size() != output_count * width) {
        return failure("owner_contract_failure", "equilibrium input derivative layout changed");
    }
    for (std::size_t output = 0; output < output_count; ++output) {
        for (std::size_t parameter = 0; parameter < parameter_count; ++parameter) {
            const auto found = std::find(
                equilibrium.input_order.begin(), equilibrium.input_order.end(),
                sidecar.active_parameter_identities[parameter]
            );
            if (found == equilibrium.input_order.end()) {
                return failure("owner_contract_failure", "equilibrium parameter order changed");
            }
            evaluation.jacobian[output * parameter_count + parameter] =
                equilibrium.output_input_derivatives[
                    output * width + static_cast<std::size_t>(found - equilibrium.input_order.begin())
                ];
        }
    }
    if (!finite(evaluation.jacobian)) {
        return failure("owner_contract_failure", "equilibrium Jacobian is nonfinite");
    }
    return evaluation;
}

} // namespace

EquilibriumOwnerEvaluator::EquilibriumOwnerEvaluator(
    epcsaft::NativeModel model,
    epcsaft::NativeModel incipient_model,
    GrepeRequest request,
    EquilibriumEvaluatorSidecar sidecar,
    bool bubble_vle
)
    : model_(std::move(model)), incipient_model_(std::move(incipient_model)),
      request_(std::move(request)), sidecar_(std::move(sidecar)), bubble_vle_(bubble_vle) {}

std::size_t EquilibriumOwnerEvaluator::parameter_count() const noexcept {
    return sidecar_.active_parameter_identities.size();
}

std::size_t EquilibriumOwnerEvaluator::output_count() const noexcept {
    return sidecar_.output_row_identities.size();
}

const std::string& EquilibriumOwnerEvaluator::parameter_contract_fingerprint()
    const noexcept {
    return sidecar_.parameter_contract_fingerprint;
}

const std::string& EquilibriumOwnerEvaluator::evaluator_fingerprint()
    const noexcept {
    return sidecar_.evaluator_fingerprint;
}

interfaces::EvaluationResult EquilibriumOwnerEvaluator::evaluate(
    const std::vector<double>& physical_parameters,
    bool exact_jacobian
) const noexcept {
    try {
        if (physical_parameters.size() != parameter_count()) {
            return failure("parameter_order_mismatch", "active parameter dimensions changed");
        }
        if (!finite(physical_parameters)) {
            return failure("invalid_parameter", "active parameters must be finite");
        }
        const epcsaft::NativeModel model = model_at(model_, physical_parameters);
        const std::vector<double> incipient_values = values_for_model(
            incipient_model_, sidecar_.active_parameter_identities, physical_parameters
        );
        const epcsaft::NativeModel incipient_model = model_at(
            incipient_model_, incipient_values
        );
        GrepeRequest request = request_;
        request.require_output_derivatives = exact_jacobian;
        const GrepeResult result = bubble_vle_
            ? solve_declared_bubble_vle(model, incipient_model, request)
            : solve_declared_homogeneous_reactive(model, request);
        return map_result(result, sidecar_, exact_jacobian);
    } catch (const std::exception& error) {
        return invalid_result(error);
    } catch (...) {
        return failure("owner_failure", "equilibrium owner evaluation failed");
    }
}

EquilibriumEvaluatorCompilation compile_equilibrium_owner(
    const epcsaft::NativeModel& model,
    const epcsaft::NativeModel& incipient_model,
    GrepeRequest request,
    EquilibriumEvaluatorSetup setup
) {
    const bool bubble_vle = admits_declared_bubble_vle(request);
    const bool homogeneous = admits_declared_homogeneous_reactive(request);
    if (!bubble_vle && !homogeneous) {
        throw std::invalid_argument("only admitted declared reactive topologies are supported");
    }
    const std::string expected_kind = bubble_vle ? kBubbleVleKind : kHomogeneousKind;
    if (!setup.owner_kind.empty() && setup.owner_kind != expected_kind) {
        throw std::invalid_argument("equilibrium evaluator owner kind does not match topology");
    }
    request.require_output_derivatives = true;
    const EquilibriumEvaluatorSidecar sidecar = make_sidecar(
        model, request, std::move(setup), bubble_vle
    );
    if (sidecar.owner_kind != expected_kind) {
        throw std::invalid_argument("equilibrium evaluator owner kind is invalid");
    }
    const auto active = model.active_coordinates();
    const auto incipient_active = incipient_model.active_coordinates();
    const auto primary_ids = identities(active);
    const auto incipient_ids = identities(incipient_active);
    if (bubble_vle
        && model.parameter_fingerprint() != incipient_model.parameter_fingerprint()) {
        throw std::invalid_argument("bubble evaluator model fingerprints do not match");
    }
    if (!unique_nonempty(primary_ids) || !unique_nonempty(incipient_ids)) {
        throw std::invalid_argument("equilibrium active parameter identities are invalid");
    }
    for (const std::string& identity : incipient_ids) {
        if (std::find(primary_ids.begin(), primary_ids.end(), identity) == primary_ids.end()) {
            throw std::invalid_argument("incipient active parameter order changed");
        }
    }
    auto evaluator = std::shared_ptr<const EquilibriumOwnerEvaluator>(
        new EquilibriumOwnerEvaluator(
            model,
            incipient_model,
            std::move(request),
            sidecar,
            bubble_vle
        )
    );
    return {std::move(evaluator)};
}

} // namespace epcsaft_engine::equilibrium
