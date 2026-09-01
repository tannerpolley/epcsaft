#include "owner_evaluator.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
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

bool all_neutral(
    const std::vector<std::size_t>& indices,
    const std::vector<int>& charges
) {
    return !indices.empty() && std::all_of(
        indices.begin(), indices.end(), [&](std::size_t index) {
            return index < charges.size() && charges[index] == 0;
        }
    );
}

bool unique_indices(const std::vector<std::size_t>& indices) {
    std::vector<std::size_t> sorted = indices;
    std::sort(sorted.begin(), sorted.end());
    return std::adjacent_find(sorted.begin(), sorted.end()) == sorted.end();
}

bool is_homogeneous_request(const GrepeRequest& request) {
    if (request.temperature.solved || request.pressure.solved
        || request.boundary.has_value() || request.directions.has_value()
        || request.phases.size() != 1 || !request.phases[0].finite
        || (request.feed.has_value() && request.feed->affine)) {
        return false;
    }
    if (!request.reaction_phase_ids.empty()
        && request.reaction_phase_ids != std::vector<std::string>{request.phases[0].identity}) {
        return false;
    }
    return std::isfinite(request.temperature.value)
        && request.temperature.value > 0.0
        && std::isfinite(request.pressure.value)
        && request.pressure.value > 0.0;
}

bool is_bubble_request(const GrepeRequest& request) {
    if (request.temperature.solved || !request.pressure.solved
        || request.boundary.has_value() || request.directions.has_value()
        || request.feed.has_value() || request.phases.size() != 2) {
        return false;
    }
    const auto finite_phase = std::find_if(
        request.phases.begin(), request.phases.end(), [](const GrepePhase& phase) {
            return phase.finite;
        }
    );
    const auto incipient_phase = std::find_if(
        request.phases.begin(), request.phases.end(), [](const GrepePhase& phase) {
            return !phase.finite;
        }
    );
    std::vector<std::size_t> expected_support(request.reaction.species_ids.size());
    std::iota(expected_support.begin(), expected_support.end(), 0);
    if (finite_phase == request.phases.end() || incipient_phase == request.phases.end()
        || finite_phase->role != PhaseRole::Liquid
        || incipient_phase->role != PhaseRole::Vapor
        || finite_phase->support_indices != expected_support
        || !unique_indices(incipient_phase->support_indices)
        || !all_neutral(incipient_phase->support_indices, request.reaction.charges)) {
        return false;
    }
    return (request.reaction_phase_ids.empty()
        || request.reaction_phase_ids == std::vector<std::string>{finite_phase->identity})
        && std::isfinite(request.temperature.value)
        && request.temperature.value > 0.0
        && std::isfinite(request.pressure.value)
        && request.pressure.value > 0.0;
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

std::vector<epcsaft::NativeActiveCoordinate> updated_coordinates(
    const std::vector<epcsaft::NativeActiveCoordinate>& baseline,
    const std::vector<epcsaft::NativeActiveCoordinate>& active,
    const std::vector<double>& values
) {
    std::vector<epcsaft::NativeActiveCoordinate> result = baseline;
    for (auto& coordinate : result) {
        const auto found = std::find_if(
            active.begin(), active.end(), [&](const auto& candidate) {
                return candidate.identity == coordinate.identity;
            }
        );
        if (found == active.end()) {
            throw std::invalid_argument("active parameter order changed");
        }
        coordinate.value = values[static_cast<std::size_t>(found - active.begin())];
    }
    return result;
}

epcsaft::NativeModel model_at(
    const epcsaft::NativeModel& model,
    const std::vector<epcsaft::NativeActiveCoordinate>& baseline,
    const std::vector<epcsaft::NativeActiveCoordinate>& active,
    const std::vector<double>& values
) {
    if (baseline.empty()) {
        return model;
    }
    return model.with_active_coordinates(updated_coordinates(baseline, active, values));
}

std::vector<std::string> identities(
    const std::vector<epcsaft::NativeActiveCoordinate>& coordinates
) {
    return epcsaft::active_parameter_identities(coordinates);
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
        || setup.parameter_contract_fingerprint.empty()
        || setup.evaluator_fingerprint.empty()) {
        throw std::invalid_argument("equilibrium evaluator fingerprints are incomplete");
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
    result.parameter_contract_fingerprint =
        std::move(setup.parameter_contract_fingerprint);
    result.evaluator_fingerprint = std::move(setup.evaluator_fingerprint);
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
    std::vector<epcsaft::NativeActiveCoordinate> active_parameters,
    std::vector<epcsaft::NativeActiveCoordinate> incipient_parameters,
    EquilibriumEvaluatorSidecar sidecar,
    bool bubble_vle
)
    : model_(std::move(model)), incipient_model_(std::move(incipient_model)),
      request_(std::move(request)), active_parameters_(std::move(active_parameters)),
      incipient_parameters_(std::move(incipient_parameters)),
      sidecar_(std::move(sidecar)), bubble_vle_(bubble_vle) {}

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
        const epcsaft::NativeModel model = model_at(
            model_, active_parameters_, active_parameters_, physical_parameters
        );
        const epcsaft::NativeModel incipient_model = model_at(
            incipient_model_, incipient_parameters_, active_parameters_, physical_parameters
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
    const bool bubble_vle = is_bubble_request(request);
    const bool homogeneous = is_homogeneous_request(request);
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
            std::vector<epcsaft::NativeActiveCoordinate>(active.begin(), active.end()),
            std::vector<epcsaft::NativeActiveCoordinate>(
                incipient_active.begin(), incipient_active.end()
            ),
            sidecar,
            bubble_vle
        )
    );
    return {std::move(evaluator), sidecar};
}

} // namespace epcsaft_engine::equilibrium
