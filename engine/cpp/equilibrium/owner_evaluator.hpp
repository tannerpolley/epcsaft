#pragma once

#include "../interfaces/owner_evaluator.hpp"
#include "grepe.hpp"

#include <cstddef>
#include <memory>
#include <string>
#include <vector>

namespace epcsaft_engine::equilibrium {

struct EvaluatorParameter final {
    std::string identity;
    std::string unit;
    std::string basis;
};

struct EvaluatorOutput final {
    std::string identity;
    std::string unit;
    std::string basis;
};

struct EvaluatorPhase final {
    std::string identity;
    std::string role;
    std::string support;
    std::vector<std::size_t> support_indices;
    bool finite;
    bool ideal_gas;
};

struct EquilibriumEvaluatorSetup final {
    std::string owner_kind;
    std::string parameter_baseline_digest;
    std::string source_spec_digest;
    std::string anchor_identity;
    std::string anchor_receipt_identity;
    std::string derivative_provenance;
    std::string topology_policy;
    std::string branch_policy;
    std::string conditioning_policy;
    std::vector<EvaluatorParameter> active_parameters;
    std::vector<EvaluatorOutput> output_rows;
};

struct EquilibriumEvaluatorSidecar final {
    std::string owner_kind;
    std::string parameter_baseline_digest;
    std::string source_spec_digest;
    std::string anchor_identity;
    std::string anchor_receipt_identity;
    std::vector<EvaluatorParameter> active_parameters;
    std::vector<EvaluatorOutput> output_rows;
    std::string derivative_provenance;
    std::vector<EvaluatorPhase> phases;
    std::string temperature_role;
    std::string pressure_role;
    std::vector<std::string> reaction_identities;
    std::vector<std::string> reaction_phase_identities;
    std::string topology_policy;
    std::string branch_policy;
    std::string conditioning_policy;
    std::string parameter_contract_fingerprint;
    std::string evaluator_fingerprint;
};

class EquilibriumOwnerEvaluator final
    : public interfaces::OwnerEvaluator {
public:
    [[nodiscard]] std::size_t parameter_count() const noexcept override;
    [[nodiscard]] std::size_t output_count() const noexcept override;
    [[nodiscard]] const std::string& parameter_contract_fingerprint()
        const noexcept override;
    [[nodiscard]] const std::string& evaluator_fingerprint()
        const noexcept override;
    [[nodiscard]] interfaces::EvaluationResult evaluate(
        const std::vector<double>& physical_parameters,
        bool exact_jacobian
    ) const noexcept override;

    [[nodiscard]] const EquilibriumEvaluatorSidecar& sidecar() const noexcept {
        return sidecar_;
    }

private:
    friend std::shared_ptr<const EquilibriumOwnerEvaluator> compile_equilibrium_owner(
        const epcsaft::NativeModel&,
        const epcsaft::NativeModel&,
        GrepeRequest,
        EquilibriumEvaluatorSetup
    );

    EquilibriumOwnerEvaluator(
        epcsaft::NativeModel model,
        epcsaft::NativeModel incipient_model,
        GrepeRequest request,
        EquilibriumEvaluatorSidecar sidecar,
        bool bubble_vle
    );

    epcsaft::NativeModel model_;
    epcsaft::NativeModel incipient_model_;
    GrepeRequest request_;
    EquilibriumEvaluatorSidecar sidecar_;
    bool bubble_vle_;
};

[[nodiscard]] std::shared_ptr<const EquilibriumOwnerEvaluator> compile_equilibrium_owner(
    const epcsaft::NativeModel& model,
    const epcsaft::NativeModel& incipient_model,
    GrepeRequest request,
    EquilibriumEvaluatorSetup setup
);

} // namespace epcsaft_engine::equilibrium
