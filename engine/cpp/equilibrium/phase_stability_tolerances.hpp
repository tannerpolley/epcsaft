#pragma once

#include <cmath>
#include <stdexcept>

namespace epcsaft_engine::equilibrium {

enum class PhaseStabilityToleranceRelation {
    AbsAtMost,
    AtLeast,
    GreaterThan,
    LessThanNegative,
    SolverTarget,
};

struct PhaseStabilityTolerance {
    const char* name;
    const char* category;
    const char* failure_meaning;
    PhaseStabilityToleranceRelation relation;
    double atol;
    double rtol;
};

struct PhaseStabilityToleranceAudit {
    const PhaseStabilityTolerance* tolerance = nullptr;
    double residual = 0.0;
    double scale = 0.0;
    double threshold = 0.0;
    bool passed = false;
};

inline constexpr PhaseStabilityTolerance kPhaseStabilityPolytopeFeasibility{
    "polytope_feasibility", "representation",
    "modified composition is outside the complete Step-1 polytope",
    PhaseStabilityToleranceRelation::AbsAtMost, 1.0e-12, 0.0,
};
inline constexpr PhaseStabilityTolerance kPhaseStabilityCompositionSum{
    "composition_sum", "representation", "composition normalization unresolved",
    PhaseStabilityToleranceRelation::AbsAtMost, 1.0e-9, 0.0,
};
inline constexpr PhaseStabilityTolerance kPhaseStabilityChargeBalance{
    "charge_balance", "representation", "scaled electroneutrality unresolved",
    PhaseStabilityToleranceRelation::AbsAtMost, 1.0e-9, 0.0,
};
inline constexpr PhaseStabilityTolerance kPhaseStabilityReconstructedIon{
    "reconstructed_ion", "representation", "eliminated-ion reconstruction is negative",
    PhaseStabilityToleranceRelation::AbsAtMost, 1.0e-9, 0.0,
};
inline constexpr PhaseStabilityTolerance kPhaseStabilityBoundActivity{
    "bound_activity", "representation", "bound activity is numerically unresolved",
    PhaseStabilityToleranceRelation::AbsAtMost, 1.0e-8, 0.0,
};
inline constexpr PhaseStabilityTolerance kPhaseStabilityRootPressure{
    "root_pressure", "root", "relative pressure closure failed",
    PhaseStabilityToleranceRelation::AbsAtMost, 1.0e-8, 0.0,
};
inline constexpr PhaseStabilityTolerance kPhaseStabilityRootLogVolumeWidth{
    "root_log_volume_width", "root", "log-volume root bracket is unresolved",
    PhaseStabilityToleranceRelation::AbsAtMost, 1.0e-9, 0.0,
};
inline constexpr PhaseStabilityTolerance kPhaseStabilityRootStationary{
    "root_stationary", "root", "pressure stationary-point refinement failed",
    PhaseStabilityToleranceRelation::AbsAtMost, 1.0e-9, 0.0,
};
inline constexpr PhaseStabilityTolerance kPhaseStabilityRootBoundary{
    "root_boundary", "root", "root is indistinguishable from Provider boundary",
    PhaseStabilityToleranceRelation::AbsAtMost, 1.0e-8, 0.0,
};
inline constexpr PhaseStabilityTolerance kPhaseStabilityRootDuplicate{
    "root_duplicate", "root", "root identity is numerically unresolved",
    PhaseStabilityToleranceRelation::AbsAtMost, 1.0e-8, 0.0,
};
inline constexpr PhaseStabilityTolerance kPhaseStabilityJointVolumeConsistency{
    "joint_volume_consistency", "representation",
    "Provider state does not match the requested joint-search molar volume",
    PhaseStabilityToleranceRelation::AbsAtMost, 1.0e-10, 0.0,
};
inline constexpr PhaseStabilityTolerance kPhaseStabilityMechanicalMargin{
    "mechanical_margin", "topology", "mechanical stability is marginal",
    PhaseStabilityToleranceRelation::GreaterThan, 1.0e-6, 0.0,
};
inline constexpr PhaseStabilityTolerance kPhaseStabilityStableObjectiveTie{
    "stable_objective_tie", "topology", "stable-root objective ordering is unresolved",
    PhaseStabilityToleranceRelation::AbsAtMost, 1.0e-8, 1.0e-9,
};
inline constexpr PhaseStabilityTolerance kPhaseStabilityTpdNegativeMargin{
    "tpd_negative_margin", "stage_i", "TPD is not materially negative",
    PhaseStabilityToleranceRelation::LessThanNegative, 1.0e-8, 0.0,
};
inline constexpr PhaseStabilityTolerance kPhaseStabilityTpdReferenceZero{
    "tpd_reference_zero", "stage_i", "TPD is nonzero at its tangent point",
    PhaseStabilityToleranceRelation::AbsAtMost, 1.0e-12, 0.0,
};
inline constexpr PhaseStabilityTolerance kPhaseStabilityLpPrimal{
    "lp_primal", "stage_ii_lp", "upper LP primal feasibility failed",
    PhaseStabilityToleranceRelation::AbsAtMost, 1.0e-9, 1.0e-8,
};
inline constexpr PhaseStabilityTolerance kPhaseStabilityLpDual{
    "lp_dual", "stage_ii_lp", "upper LP dual feasibility failed",
    PhaseStabilityToleranceRelation::AbsAtMost, 1.0e-9, 1.0e-8,
};
inline constexpr PhaseStabilityTolerance kPhaseStabilityLpComplementarity{
    "lp_complementarity", "stage_ii_lp", "upper LP complementarity failed",
    PhaseStabilityToleranceRelation::AbsAtMost, 1.0e-8, 0.0,
};
inline constexpr PhaseStabilityTolerance kPhaseStabilityLpActiveCut{
    "lp_active_cut", "stage_ii_lp_diagnostic", "cut is not diagnostically active",
    PhaseStabilityToleranceRelation::AbsAtMost, 1.0e-7, 0.0,
};
inline constexpr PhaseStabilityTolerance kPhaseStabilityStage2KktPrimal{
    "stage2_kkt_primal", "stage_ii_kkt",
    "Step-5 original-coordinate primal feasibility failed",
    PhaseStabilityToleranceRelation::AbsAtMost, 1.0e-8, 0.0,
};
inline constexpr PhaseStabilityTolerance kPhaseStabilityStage2KktDualSign{
    "stage2_kkt_dual_sign", "stage_ii_kkt",
    "Step-5 bound or inequality multiplier sign failed",
    PhaseStabilityToleranceRelation::AbsAtMost, 1.0e-9, 0.0,
};
inline constexpr PhaseStabilityTolerance kPhaseStabilityStage2KktPullback{
    "stage2_kkt_pullback", "stage_ii_kkt",
    "Step-5 physical-to-modified dual reconstruction failed",
    PhaseStabilityToleranceRelation::AbsAtMost, 1.0e-9, 1.0e-9,
};
inline constexpr PhaseStabilityTolerance kPhaseStabilityStage2KktStationarity{
    "stage2_kkt_stationarity", "stage_ii_kkt",
    "Step-5 original-coordinate stationarity failed",
    PhaseStabilityToleranceRelation::AbsAtMost, 1.0e-6, 0.0,
};
inline constexpr PhaseStabilityTolerance kPhaseStabilityStage2KktComplementarity{
    "stage2_kkt_complementarity", "stage_ii_kkt",
    "Step-5 bound or inequality complementarity failed",
    PhaseStabilityToleranceRelation::AbsAtMost, 1.0e-8, 0.0,
};
inline constexpr PhaseStabilityTolerance kPhaseStabilityStage2KktRankPivot{
    "stage2_kkt_rank_pivot", "stage_ii_kkt",
    "Step-5 active-constraint Jacobian is rank deficient",
    PhaseStabilityToleranceRelation::GreaterThan, 1.0e-12, 0.0,
};
inline constexpr PhaseStabilityTolerance kPhaseStabilityFarkasRowSign{
    "farkas_row_sign", "stage_iii_lp_certificate",
    "phase-feasibility Farkas row multiplier has the wrong bound sign",
    PhaseStabilityToleranceRelation::AbsAtMost, 1.0e-10, 0.0,
};
inline constexpr PhaseStabilityTolerance kPhaseStabilityFarkasDual{
    "farkas_dual", "stage_iii_lp_certificate",
    "phase-feasibility Farkas column dual is incompatible with variable bounds",
    PhaseStabilityToleranceRelation::AbsAtMost, 1.0e-10, 0.0,
};
inline constexpr PhaseStabilityTolerance kPhaseStabilityFarkasContradiction{
    "farkas_contradiction", "stage_iii_lp_certificate",
    "phase-feasibility Farkas contradiction margin is not strict",
    PhaseStabilityToleranceRelation::GreaterThan, 1.0e-9, 1.0e-9,
};
inline constexpr PhaseStabilityTolerance kPhaseStabilityStage6Gap{
    "step6_gap", "stage_ii_step6", "same-major upper/lower gap failed",
    PhaseStabilityToleranceRelation::AbsAtMost, 1.0e-8, 0.0,
};
inline constexpr PhaseStabilityTolerance kPhaseStabilityPaperStep6Gap{
    "paper_step6_gap", "stage_ii_paper",
    "PHASE_STABILITY Eq. (66) upper/lower agreement failed",
    PhaseStabilityToleranceRelation::AbsAtMost, 5.0e-2, 0.0,
};
inline constexpr PhaseStabilityTolerance kPhaseStabilityPaperStep6Derivative{
    "paper_step6_derivative", "stage_ii_paper",
    "PHASE_STABILITY Eq. (66) derivative agreement failed",
    PhaseStabilityToleranceRelation::AbsAtMost, 1.0e-8, 0.5,
};
inline constexpr PhaseStabilityTolerance kPhaseStabilityPaperStep6PackingDistinct{
    "paper_step6_packing_distinct", "phase_identity",
    "PHASE_STABILITY Eq. (66) packing fractions are not distinct",
    PhaseStabilityToleranceRelation::AtLeast, 1.0e-3, 0.0,
};
inline constexpr PhaseStabilityTolerance kPhaseStabilityPaperStep6CompositionDistinct{
    "paper_step6_composition_distinct", "phase_identity",
    "PHASE_STABILITY Eq. (66) modified compositions are not distinct",
    PhaseStabilityToleranceRelation::AtLeast, 1.0e-3, 0.0,
};
inline constexpr PhaseStabilityTolerance kPhaseStabilityMRepresentationEquivalent{
    "m_representation_equivalent", "phase_identity",
    "Step-5 terminal differs in modified composition or relative molar volume",
    PhaseStabilityToleranceRelation::AbsAtMost, 1.0e-8, 0.0,
};
inline constexpr PhaseStabilityTolerance kPhaseStabilityStage3ModifiedBalance{
    "stage3_modified_balance", "stage_iii_physical", "modified material balance failed",
    PhaseStabilityToleranceRelation::AbsAtMost, 1.0e-8, 0.0,
};
inline constexpr PhaseStabilityTolerance kPhaseStabilityStage3ExplicitBalance{
    "stage3_explicit_balance", "stage_iii_physical", "explicit material balance failed",
    PhaseStabilityToleranceRelation::AbsAtMost, 1.0e-8, 0.0,
};
inline constexpr PhaseStabilityTolerance kPhaseStabilityStage3Charge{
    "stage3_charge", "stage_iii_physical", "scaled phase electroneutrality failed",
    PhaseStabilityToleranceRelation::AbsAtMost, 1.0e-9, 0.0,
};
inline constexpr PhaseStabilityTolerance kPhaseStabilityStage3Pressure{
    "stage3_pressure", "stage_iii_physical", "relative phase pressure closure failed",
    PhaseStabilityToleranceRelation::AbsAtMost, 1.0e-8, 0.0,
};
inline constexpr PhaseStabilityTolerance kPhaseStabilityStage3Stationarity{
    "stage3_stationarity", "stage_iii_kkt", "Stage-III stationarity failed",
    PhaseStabilityToleranceRelation::AbsAtMost, 1.0e-7, 0.0,
};
inline constexpr PhaseStabilityTolerance kPhaseStabilityStage3DualSign{
    "stage3_dual_sign", "stage_iii_kkt", "Stage-III multiplier sign failed",
    PhaseStabilityToleranceRelation::AbsAtMost, 1.0e-9, 0.0,
};
inline constexpr PhaseStabilityTolerance kPhaseStabilityStage3Complementarity{
    "stage3_complementarity", "stage_iii_kkt", "Stage-III complementarity failed",
    PhaseStabilityToleranceRelation::AbsAtMost, 1.0e-8, 0.0,
};
inline constexpr PhaseStabilityTolerance kPhaseStabilityPaperFreeEnergyGap{
    "paper_free_energy_gap", "stage_iii_paper",
    "PHASE_STABILITY Eq. (68) free-energy convergence failed",
    PhaseStabilityToleranceRelation::AbsAtMost, 1.0e-4, 0.0,
};
inline constexpr PhaseStabilityTolerance kPhaseStabilityPaperPotentialRatio{
    "paper_potential_ratio", "stage_iii_paper",
    "PHASE_STABILITY Eq. (69) modified-potential convergence failed",
    PhaseStabilityToleranceRelation::AbsAtMost, 1.0e-2, 0.0,
};
inline constexpr PhaseStabilityTolerance kPhaseStabilityPhaseRetirementMargin{
    "phase_retirement_margin", "stage_iii_kkt", "phase retirement evidence is marginal",
    PhaseStabilityToleranceRelation::GreaterThan, 1.0e-8, 0.0,
};
inline constexpr PhaseStabilityTolerance kPhaseStabilityPhaseMerge{
    "phase_merge", "phase_identity", "phases are not numerical duplicates",
    PhaseStabilityToleranceRelation::AbsAtMost, 1.0e-5, 0.0,
};
inline constexpr PhaseStabilityTolerance kPhaseStabilityPaperPhaseCoalescence{
    "paper_phase_coalescence", "phase_identity",
    "three-or-more-phase solve has no paper-scale coalescent pair",
    PhaseStabilityToleranceRelation::AbsAtMost, 1.0e-2, 0.0,
};
inline constexpr PhaseStabilityTolerance kPhaseStabilityIpoptTarget{
    "ipopt_target", "solver", "Ipopt target tolerance was not met",
    PhaseStabilityToleranceRelation::SolverTarget, 1.0e-10, 0.0,
};
inline constexpr PhaseStabilityTolerance kPhaseStabilityStage3IpoptTarget{
    "stage3_ipopt_target", "stage_iii_solver",
    "Stage-III Ipopt target tolerance was not met",
    PhaseStabilityToleranceRelation::SolverTarget, 1.0e-8, 0.0,
};
inline constexpr PhaseStabilityTolerance kPhaseStabilityIpoptConstraint{
    "ipopt_constraint", "solver", "Ipopt constraint target was not met",
    PhaseStabilityToleranceRelation::SolverTarget, 1.0e-10, 0.0,
};
[[nodiscard]] inline PhaseStabilityToleranceAudit audit_phase_stability_tolerance(
    const PhaseStabilityTolerance& tolerance,
    double residual,
    double scale = 0.0
) {
    if (!std::isfinite(residual) || !std::isfinite(scale) || scale < 0.0) {
        throw std::invalid_argument(
            "PHASE_STABILITY tolerance residual must be finite and scale must be "
            "finite and nonnegative"
        );
    }
    const double threshold = tolerance.atol + tolerance.rtol * scale;
    bool passed = false;
    switch (tolerance.relation) {
        case PhaseStabilityToleranceRelation::AbsAtMost:
        case PhaseStabilityToleranceRelation::SolverTarget:
            passed = std::abs(residual) <= threshold;
            break;
        case PhaseStabilityToleranceRelation::AtLeast:
            passed = residual >= threshold;
            break;
        case PhaseStabilityToleranceRelation::GreaterThan:
            passed = residual > threshold;
            break;
        case PhaseStabilityToleranceRelation::LessThanNegative:
            passed = residual < -threshold;
            break;
    }
    return {&tolerance, residual, scale, threshold, passed};
}

}  // namespace epcsaft_engine::equilibrium
