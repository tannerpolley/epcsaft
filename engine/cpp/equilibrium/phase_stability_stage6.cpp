#include "phase_stability_stage6.hpp"
#include "phase_stability_tolerances.hpp"

#include <algorithm>
#include <cmath>

namespace epcsaft_engine::equilibrium {

PhaseStabilityStage6Result run_phase_stability_stage6(
    const PhaseStabilityStage1Result& step1,
    const PhaseStabilityStage4Result& step4,
    const PhaseStabilityPersistentState& state,
    const PhaseStabilityPackingFractionEvaluator& packing_fraction
) {
    PhaseStabilityStage6Result result;
    result.timing.invocation_count = 1;
    if (step4.status != "complete" || !step1.coordinates
        || !packing_fraction) {
        result.reason = "invalid_step6_input";
        return result;
    }
    for (const PhaseStabilityMPoint& point : state.M) {
        if (point.reduced_gibbs_gradient.size() < state.feed.size()) {
            result.reason = "missing_step6_gradient";
            return result;
        }
        PhaseStabilityCandidateDecision decision;
        decision.insertion_id = point.insertion_id;
        double lower = point.reduced_gibbs;
        for (std::size_t index = 0; index < state.feed.size(); ++index) {
            lower += state.multipliers[index] * (
                state.feed[index]
                - point.independent_modified_fractions[index]
            );
        }
        decision.gap_passed = audit_phase_stability_tolerance(
            kPhaseStabilityPaperStep6Gap, state.upper_bound - lower
        ).passed;
        decision.derivative_passed = true;
        for (std::size_t index = 0; index < state.feed.size(); ++index) {
            if (point.independent_modified_fractions[index]
                <= step1.coordinates->independent_lower_bounds[index]
                    + kPhaseStabilityBoundActivity.atol) {
                continue;
            }
            decision.derivative_passed = decision.derivative_passed
                && audit_phase_stability_tolerance(
                    kPhaseStabilityPaperStep6Derivative,
                    point.reduced_gibbs_gradient[index]
                        - state.multipliers[index],
                    std::abs(state.multipliers[index])
                ).passed;
        }
        if (!decision.gap_passed || !decision.derivative_passed) {
            decision.reason = decision.gap_passed
                ? "derivative_failed"
                : "bound_gap_failed";
            result.decisions.push_back(std::move(decision));
            continue;
        }
        PhaseStabilityMPoint candidate = point;
        candidate.packing_fraction = packing_fraction(
            point.independent_modified_fractions, point.volume
        );
        decision.pairwise_distinct = std::all_of(
            result.candidates.begin(),
            result.candidates.end(),
            [&](const PhaseStabilityMPoint& retained) {
                if (std::abs(
                        candidate.packing_fraction
                        - retained.packing_fraction
                    ) >= kPhaseStabilityPaperStep6PackingDistinct.atol) {
                    return true;
                }
                for (std::size_t index = 0; index < state.feed.size(); ++index) {
                    if (std::abs(
                            candidate.independent_modified_fractions[index]
                            - retained.independent_modified_fractions[index]
                        ) >= kPhaseStabilityPaperStep6CompositionDistinct.atol) {
                        return true;
                    }
                }
                return false;
            }
        );
        decision.retained = decision.pairwise_distinct;
        decision.reason = decision.retained ? "retained" : "not_distinct";
        if (decision.retained) {
            result.candidates.push_back(std::move(candidate));
        }
        result.decisions.push_back(std::move(decision));
    }
    result.status = "complete";
    result.reason = result.candidates.size() >= 2
        ? "candidate_set"
        : "insufficient_candidates";
    result.timing.terminal_status = result.status;
    result.timing.terminal_reason = result.reason;
    result.timing.next_step = result.candidates.size() >= 2 ? 8 : 7;
    return result;
}

}  // namespace epcsaft_engine::equilibrium
