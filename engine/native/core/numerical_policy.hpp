#pragma once

#include "core/physical_constants.hpp"

#include <cmath>
#include <limits>

// Every number in the engine that is not a physical constant, with the reason it
// has the value it has. Four categories:
//
//   1  physical constant        -- measured or SI-defined; lives in
//                                  physical_constants.hpp, never here
//   2  thermodynamic constraint -- forced by the equations; changing it makes
//                                  the model wrong, not merely different
//   3  algorithm tolerance      -- a convergence or discretisation choice with a
//                                  stated source and a stated consequence
//   4  unjustified legacy       -- a number nobody can source; this category is empty
//
// Category 4 is empty and must stay empty. An entry that cannot be justified
// into category 3 with a source does not get carried; a provisional category-3
// choice may remain only when its rationale, predicate, failure meaning and
// unresolved accuracy are explicit. The audit found the drifted values below
// by comparing against the paper, and they are pinned here so a solver stage
// cannot quietly move them again.
//
// Stages add their own entries as they are written. Nothing outside this header
// may define a tolerance.

namespace epcsaft::policy {

enum class Category {
    ThermodynamicConstraint = 2,
    AlgorithmTolerance = 3,
};

struct Tolerance final {
    const char* name;
    double value;
    Category category;
    const char* source;
    const Tolerance* upstream = nullptr;
    double factor = 1.0;
};

constexpr Tolerance derived(const Tolerance& upstream, double factor,
                            const char* name, const char* why) {
    if (!(factor >= 1.0)) throw "derived tolerance tighter than its input";
    return {name, upstream.value * factor, Category::AlgorithmTolerance, why,
            &upstream, factor};
}

// ---------------------------------------------------------------- constraints

// Dimensionless packing fraction. Reject eta >= the equal-sphere
// face-centred-cubic geometric limit as a conservative PC-SAFT admission
// bound; this is a domain constraint, not a mathematical singularity claim.
inline constexpr Tolerance CLOSE_PACKING_FRACTION{"CLOSE_PACKING_FRACTION",
    PI / (3.0 * 1.4142135623730951),
    Category::ThermodynamicConstraint,
    "equal-sphere face-centred-cubic admission bound; eta >= value is a domain failure",
};

// A site fraction is a probability. Outside (0, 1] the mass-action system has
// left the domain and the solve has failed, whatever its residual says.
inline constexpr Tolerance SITE_FRACTION_FLOOR{"SITE_FRACTION_FLOOR",
    0.0,
    Category::ThermodynamicConstraint,
    "X_a is the unbonded fraction of site a, so 0 < X_a <= 1",
};
inline constexpr Tolerance SITE_FRACTION_CEILING{"SITE_FRACTION_CEILING",
    1.0,
    Category::ThermodynamicConstraint,
    "X_a is the unbonded fraction of site a, so 0 < X_a <= 1",
};

// ---------------------------------------------------------------- association

// The mass-action residual is dimensionless and O(1), so this is about 1000 ulp
// -- tight enough that the site fractions do not limit the derivative accuracy,
// loose enough to be reachable in double precision.
inline constexpr Tolerance ASSOCIATION_RESIDUAL{"ASSOCIATION_RESIDUAL",
    2.0e-13,
    Category::AlgorithmTolerance,
    "carried unchanged; Newton reaches it in 3-8 steps from a warm start",
};

inline constexpr Tolerance ASSOCIATION_ACTION_RESIDUAL{"ASSOCIATION_ACTION_RESIDUAL",
    1.0e-10, Category::AlgorithmTolerance,
    "approved P4 normalized mass-action coefficient residual; exceeding it makes the action unavailable"};

// A budget, not a convergence criterion. Newton with the analytic Jacobian
// converges in single digits; reaching this many steps means the step selection
// is failing and the caller should hear about it rather than get a wrong answer
// from a truncated solve.
inline constexpr Tolerance ASSOCIATION_MAX_ITERATIONS{"ASSOCIATION_MAX_ITERATIONS",
    50.0,
    Category::AlgorithmTolerance,
    "budget; measured 3-8 Newton steps, the substitution solver it replaces "
    "needed 126-31056",
};

// Below this the Jacobian's smallest pivot no longer separates from its largest
// and the Newton step is not trustworthy; the solver falls back to substitution.
inline constexpr Tolerance ASSOCIATION_JACOBIAN_CONDITION{"ASSOCIATION_JACOBIAN_CONDITION",
    1.0e-13,
    Category::AlgorithmTolerance,
    "relative pivot floor in the dense solve",
};

inline constexpr Tolerance ASSOCIATION_JACOBIAN_RANK{"ASSOCIATION_JACOBIAN_RANK",
    std::numeric_limits<double>::epsilon(), Category::AlgorithmTolerance,
    "FullPivLU relative rank base; multiply by the site count"};

// ------------------------------------------------------------------- density

// The pressure row is solved in reduced units where the residual is O(1).
inline constexpr Tolerance DENSITY_ROOT_RESIDUAL = derived(ASSOCIATION_RESIDUAL, 5.0, "DENSITY_ROOT_RESIDUAL",
    "bracketed Newton on the reduced pressure residual; the tolerance is on "
    "the reduced residual; the delivered relative pressure error is |r|/z, "
    "about 1e-8 on a dense liquid root, measured worst 1.7e-8"
);

inline constexpr Tolerance DENSITY_MAX_ITERATIONS{"DENSITY_MAX_ITERATIONS",
    60.0,
    Category::AlgorithmTolerance,
    "budget; 5-8 steps from an anchor, 20-40 for a cold bracket",
};

// The coarse sweep that brackets the roots before Newton refines them. The tree
// this replaces used 360 intervals here and a second 64-interval scan subdivided
// to depth 8 elsewhere, for a root that a bracket plus Newton finds in single
// digits. Two scanners became one; this is its only discretisation number.
inline constexpr Tolerance DENSITY_BRACKET_INTERVALS{"DENSITY_BRACKET_INTERVALS",
    40.0,
    Category::AlgorithmTolerance,
    "coarse bracket over log volume; refined by Newton, not by subdivision",
};

inline constexpr Tolerance DENSITY_LOWER_FRACTION{"DENSITY_LOWER_FRACTION",
    1.0e-8,
    Category::AlgorithmTolerance,
    "Stage 2 brief; logarithmic density scan starts at a dilute ideal-gas fraction",
};

inline constexpr Tolerance DENSITY_UPPER_FRACTION{"DENSITY_UPPER_FRACTION",
    0.999,
    Category::AlgorithmTolerance,
    "Stage 2 brief; logarithmic density scan stops below close packing",
};

inline constexpr Tolerance DENSITY_BISECTION_FRACTION{"DENSITY_BISECTION_FRACTION",
    0.5,
    Category::AlgorithmTolerance,
    "Stage 2 brief; safeguarded Newton uses the bracket midpoint",
};
inline constexpr Tolerance DENSITY_ANCHOR_LOG_WINDOW{"DENSITY_ANCHOR_LOG_WINDOW", 0.009950330853168083,
    Category::AlgorithmTolerance, "approved P2 local density envelope log1p(0.01)"};
inline constexpr Tolerance DENSITY_ANCHOR_SCAN_INTERVALS{"DENSITY_ANCHOR_SCAN_INTERVALS", 5.0,
    Category::AlgorithmTolerance, "approved P2 centered local anchor uniqueness scan"};

inline constexpr Tolerance COMPOSITION_SUM_TOLERANCE{"COMPOSITION_SUM_TOLERANCE",
    1.0e-12,
    Category::AlgorithmTolerance,
    "carried composition validation threshold from electrolyte_reference.hpp",
};

inline constexpr Tolerance REFERENCE_MOLALITY_ONE{"REFERENCE_MOLALITY_ONE",
    1.0e-6,
    Category::AlgorithmTolerance,
    "electrolyte_reference.hpp molality ladder",
};
inline constexpr Tolerance REFERENCE_MOLALITY_TWO{"REFERENCE_MOLALITY_TWO",
    1.0e-8,
    Category::AlgorithmTolerance,
    "electrolyte_reference.hpp molality ladder",
};
inline constexpr Tolerance REFERENCE_MOLALITY_THREE{"REFERENCE_MOLALITY_THREE",
    1.0e-10,
    Category::AlgorithmTolerance,
    "electrolyte_reference.hpp molality ladder",
};
inline constexpr Tolerance REFERENCE_MOLALITY_FOUR{"REFERENCE_MOLALITY_FOUR",
    1.0e-12,
    Category::AlgorithmTolerance,
    "electrolyte_reference.hpp molality ladder",
};
inline constexpr Tolerance REFERENCE_CONVERGENCE_WITNESS{"REFERENCE_CONVERGENCE_WITNESS",
    1.0e-14,
    Category::AlgorithmTolerance,
    "electrolyte_reference.hpp convergence witness",
};
inline constexpr Tolerance REFERENCE_CONVERGENCE_TOLERANCE{"REFERENCE_CONVERGENCE_TOLERANCE",
    5.0e-5,
    Category::AlgorithmTolerance,
    "electrolyte_reference.hpp reference convergence threshold",
};

// ---------------------------------------------------------------- stability
//
// Pereira, Jackson, Galindo and Adjiman, Comput. Chem. Eng. 36 (2012) 99-118.
// The two entries marked "drifted" were changed in the tree this replaces, with
// measured consequences; they are pinned at the published values here so the
// stability stage starts from the paper.

// Gap between the upper and lower bounds within which a candidate phase is
// accepted. The old tree used 5e-2, five times the paper's value, so candidates
// qualified before the bounds had closed: the third stage was entered early,
// took 14-42 percent of all evaluations on the measured feeds, and had 60-90
// percent of its solutions rejected again by the free-energy test.
inline constexpr Tolerance STABILITY_CANDIDATE_GAP{"STABILITY_CANDIDATE_GAP",
    1.0e-2,
    Category::AlgorithmTolerance,
    "Pereira 2012 eps_b; drifted to 5e-2 in the superseded tree",
};

// Lower bound on a composition coordinate in the inner problem. The old tree
// scaled this by the feed (1e-10 f_k), which put solutions closer to the bound
// than the paper's trace treatment assumes.
inline constexpr Tolerance STABILITY_COMPOSITION_FLOOR{"STABILITY_COMPOSITION_FLOOR",
    1.0e-8,
    Category::AlgorithmTolerance,
    "dimensionless x_i coordinate; require x_i >= floor for log(x_i); Pereira 2012; provisional domain cutoff, not an accuracy bound",
};
inline constexpr Tolerance STABILITY_PHASE_AMOUNT_FLOOR{"STABILITY_PHASE_AMOUNT_FLOOR",
    1.0e-8,
    Category::AlgorithmTolerance,
    "dimensionless n_phase/feed_total; provisional free-phase disappearance topology cutoff; not a stability result or accuracy bound",
};

inline constexpr Tolerance STABILITY_MULTIPLIER_GAP{"STABILITY_MULTIPLIER_GAP",
    0.5,
    Category::AlgorithmTolerance,
    "Pereira 2012 eps_lambda",
};

inline constexpr Tolerance STABILITY_COMPOSITION_GAP{"STABILITY_COMPOSITION_GAP",
    1.0e-3,
    Category::AlgorithmTolerance,
    "Pereira 2012 eps_x",
};

inline constexpr Tolerance STABILITY_PACKING_GAP{"STABILITY_PACKING_GAP",
    1.0e-3,
    Category::AlgorithmTolerance,
    "Pereira 2012 eps_eta",
};

// ------------------------------------------------------------- regression

// Ceres trust-region stopping criteria for the bounded least-squares fit (ALG-FIT). The
// caller's FitControls override them; these are the only defaults and are not accuracy claims.
inline constexpr Tolerance FIT_FUNCTION_TOLERANCE{"FIT_FUNCTION_TOLERANCE",
    1.0e-12,
    Category::AlgorithmTolerance,
    "relative cost decrease below which the trust-region iteration stops",
};
inline constexpr Tolerance FIT_GRADIENT_TOLERANCE{"FIT_GRADIENT_TOLERANCE",
    1.0e-12,
    Category::AlgorithmTolerance,
    "infinity norm of the projected optimizer-coordinate gradient at which the fit stops",
};
inline constexpr Tolerance FIT_PARAMETER_TOLERANCE{"FIT_PARAMETER_TOLERANCE",
    1.0e-10,
    Category::AlgorithmTolerance,
    "relative optimizer-coordinate step below which the fit stops",
};
// Singular values of the weighted optimizer Jacobian below this fraction of the largest are
// treated as a null direction: covariance is withheld as rank deficient rather than inverted.
// The published object is (J^T J)^-1, whose condition is cond(J)^2 = 1e28 at the threshold, so
// a retained direction can carry about one percent SVD error in its smallest singular value.
inline constexpr Tolerance FIT_RANK_RECIPROCAL_CONDITION{"FIT_RANK_RECIPROCAL_CONDITION",
    1.0e-14,
    Category::AlgorithmTolerance,
    "reciprocal condition number separating an identifiable direction from a null direction",
};

// ------------------------------------------------------------- equilibrium

inline constexpr Tolerance EQUILIBRIUM_ROW_RESIDUAL{"EQUILIBRIUM_ROW_RESIDUAL",
    1.0e-10,
    Category::AlgorithmTolerance,
    "dimensionless equilibrium-row infinity norm; two decades above the "
    "association residual entering the assembly floor",
};

// The optimizer works below the caller's row criterion: the recomputed row maximum tracks
// the Ipopt tolerance and overshoots it by a hair (measured 1.0002-1.02x), so classifying
// at the same value is systematically marginal. Inner Ipopt tolerances are
// min(row_tolerance, EQUILIBRIUM_ROW_RESIDUAL) times this factor; the caller criterion alone
// classifies rows. Owner decision 2026-09-16 (E4); provisional, not an accuracy bound.
inline constexpr Tolerance EQUILIBRIUM_INNER_TOLERANCE_FACTOR{"EQUILIBRIUM_INNER_TOLERANCE_FACTOR",
    0.1,
    Category::AlgorithmTolerance,
    "provisional decade below the caller row criterion for the Ipopt tol, acceptable_tol, "
    "dual_inf_tol and constr_viol_tol settings; the caller criterion classifies rows",
};

inline constexpr Tolerance EQUILIBRIUM_MAX_ITERATIONS{"EQUILIBRIUM_MAX_ITERATIONS",
    200.0,
    Category::AlgorithmTolerance,
    "budget; analytic equilibrium Jacobian with Ipopt",
};

inline constexpr Tolerance PHASE_COINCIDENCE_DISTANCE{"PHASE_COINCIDENCE_DISTANCE",
    1.0e-6,
    Category::AlgorithmTolerance,
    "phase amount columns become numerically dependent below this composition "
    "or log-packing distance",
};

inline constexpr Tolerance PACKING_FRACTION_CEILING{"PACKING_FRACTION_CEILING",
    0.74,
    Category::AlgorithmTolerance,
    "fluid-coordinate ceiling below close packing",
};

// Floor of the log-packing chart, shared by the driver's variable bounds and the
// continuation predictor's domain test so both admit the same coordinates.
inline constexpr Tolerance PACKING_FRACTION_FLOOR{"PACKING_FRACTION_FLOOR",
    1.0e-12,
    Category::AlgorithmTolerance,
    "fluid-coordinate floor of the log-packing chart",
};

// A budget with a name. The old tree's global search stopped at 3250 samples and
// reported that as a one-phase result, so a decision that ran out of work was
// indistinguishable from a decision that converged. Reaching this is a named
// stop the caller can see, not an answer.
inline constexpr Tolerance STABILITY_SEARCH_BUDGET{"STABILITY_SEARCH_BUDGET",
    300.0,
    Category::AlgorithmTolerance,
    "budget; the floor analysis puts a binary decision at 40-150 evaluations",
};

// Below this reduced screening length the closed form of the finite-size correction
// loses precision to cancellation -- log(1+r) - r + r^2/2 is a difference of terms that
// agree to leading order -- so the series is used instead. Both forms are exact; this is
// where double precision stops preferring one.
//
// It is also a recording-time predicate: the branch taken is baked into a tape, so a tape
// recorded on one side of it is wrong on the other. That is why the tape cache is keyed on
// the predicate signature and not on topology alone.
inline constexpr Tolerance SCREENING_SERIES_SWITCH{"SCREENING_SERIES_SWITCH",
    1.0e-3,
    Category::AlgorithmTolerance,
    "cancellation floor of the closed-form chi; carried from the superseded tree",
};

// --------------------------------------------------------------- convergence

inline constexpr Tolerance FLOOR_MARGIN{
    "FLOOR_MARGIN",
    10.0,
    Category::AlgorithmTolerance,
    "provisional guessed decade of slack across first-order bound norm inequalities; "
    "predicate: residual <= FLOOR_MARGIN * floor; failure meaning: classified "
    "converged_to_floor; numerical adequacy/calibration remains open, not an accuracy bound",
};

inline constexpr Tolerance STALL_WINDOW{
    "STALL_WINDOW",
    5.0,
    Category::AlgorithmTolerance,
    "provisional guessed window below Ipopt acceptable_iter default 15; predicate: "
    "residual history has STALL_WINDOW entries without contraction; failure meaning: "
    "classified stalled; numerical adequacy/calibration remains open, not an accuracy bound",
};

}  // namespace epcsaft::policy
