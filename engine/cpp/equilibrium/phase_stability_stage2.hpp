#pragma once

#include <optional>
#include <string>
#include <vector>

#include "phase_stability_stage1.hpp"

namespace epcsaft_engine::equilibrium {

enum class PhaseStabilityStage2Outcome {
    NegativeWitness,
    NoNegativeWitnessDetected,
    Indeterminate,
};

struct PhaseStabilityTpdEvaluation {
    double value = 0.0;
    std::vector<double> gradient;
    std::vector<double> hessian;
};

struct PhaseStabilityStage2Result {
    PhaseStabilityStage2Outcome outcome = PhaseStabilityStage2Outcome::Indeterminate;
    std::string reason = "not_run";
    std::optional<PhaseStabilityPressureEnvelopeResult> reference_envelope;
    std::optional<PhaseStabilityStateEvaluation> reference;
    std::optional<PhaseStabilityStageICandidate> negative_witness;
    std::optional<double> minimum_tpd;
    PhaseStabilityStageTiming timing;
};

[[nodiscard]] PhaseStabilityTpdEvaluation evaluate_phase_stability_tpd(
    const PhaseStabilityStateEvaluation& reference,
    const std::vector<double>& feed,
    const PhaseStabilityStateEvaluation& trial,
    const std::vector<double>& independent
);

[[nodiscard]] PhaseStabilityStage2Result run_phase_stability_stage2(
    const PhaseStabilityStage1Result& step1,
    const PhaseStabilityStateEvaluator& evaluator,
    int search_work_budget
);

}  // namespace epcsaft_engine::equilibrium
