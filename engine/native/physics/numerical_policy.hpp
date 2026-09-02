#pragma once

#include "physics/physical_constants.hpp"

#include <cmath>

// Every number in the engine that is not a physical constant, with the reason it
// has the value it has. Four categories:
//
//   1  physical constant        -- measured or SI-defined; lives in
//                                  physical_constants.hpp, never here
//   2  thermodynamic constraint -- forced by the equations; changing it makes
//                                  the model wrong, not merely different
//   3  algorithm tolerance      -- a convergence or discretisation choice with a
//                                  stated source and a stated consequence
//   4  unjustified legacy       -- a number nobody can source
//
// Category 4 is empty and must stay empty. An entry that cannot be justified
// into category 3 with a source does not get carried; the audit found the
// drifted values below by comparing against the paper, and they are pinned here
// so a solver stage cannot quietly move them again.
//
// Stages add their own entries as they are written. Nothing outside this header
// may define a tolerance.

namespace epcsaft::policy {

enum class Category {
    ThermodynamicConstraint = 2,
    AlgorithmTolerance = 3,
};

struct Tolerance final {
    double value;
    Category category;
    const char* source;
};

// ---------------------------------------------------------------- constraints

// Above this the hard-sphere term's logarithm has no real value: the packing
// fraction of identical spheres in a face-centred cubic lattice. Not a
// tolerance -- a state at or above it is not a fluid.
inline const Tolerance CLOSE_PACKING_FRACTION{
    PI / (3.0 * 1.4142135623730951),
    Category::ThermodynamicConstraint,
    "face-centred cubic close packing, pi/(3 sqrt 2)",
};

// A site fraction is a probability. Outside (0, 1] the mass-action system has
// left the domain and the solve has failed, whatever its residual says.
inline const Tolerance SITE_FRACTION_FLOOR{
    0.0,
    Category::ThermodynamicConstraint,
    "X_a is the unbonded fraction of site a, so 0 < X_a <= 1",
};
inline const Tolerance SITE_FRACTION_CEILING{
    1.0,
    Category::ThermodynamicConstraint,
    "X_a is the unbonded fraction of site a, so 0 < X_a <= 1",
};

// ---------------------------------------------------------------- association

// The mass-action residual is dimensionless and O(1), so this is about 1000 ulp
// -- tight enough that the site fractions do not limit the derivative accuracy,
// loose enough to be reachable in double precision.
inline const Tolerance ASSOCIATION_RESIDUAL{
    2.0e-13,
    Category::AlgorithmTolerance,
    "carried unchanged; Newton reaches it in 3-8 steps from a warm start",
};

// A budget, not a convergence criterion. Newton with the analytic Jacobian
// converges in single digits; reaching this many steps means the step selection
// is failing and the caller should hear about it rather than get a wrong answer
// from a truncated solve.
inline const Tolerance ASSOCIATION_MAX_ITERATIONS{
    50.0,
    Category::AlgorithmTolerance,
    "budget; measured 3-8 Newton steps, the substitution solver it replaces "
    "needed 126-31056",
};

// Below this the Jacobian's smallest pivot no longer separates from its largest
// and the Newton step is not trustworthy; the solver falls back to substitution.
inline const Tolerance ASSOCIATION_JACOBIAN_CONDITION{
    1.0e-13,
    Category::AlgorithmTolerance,
    "relative pivot floor in the dense solve",
};

// ------------------------------------------------------------------- density

// The pressure row is solved in reduced units where the residual is O(1).
inline const Tolerance DENSITY_ROOT_RESIDUAL{
    1.0e-12,
    Category::AlgorithmTolerance,
    "bracketed Newton on the pressure row",
};

inline const Tolerance DENSITY_MAX_ITERATIONS{
    60.0,
    Category::AlgorithmTolerance,
    "budget; 5-8 steps from an anchor, 20-40 for a cold bracket",
};

// The coarse sweep that brackets the roots before Newton refines them. The tree
// this replaces used 360 intervals here and a second 64-interval scan subdivided
// to depth 8 elsewhere, for a root that a bracket plus Newton finds in single
// digits. Two scanners became one; this is its only discretisation number.
inline const Tolerance DENSITY_BRACKET_INTERVALS{
    40.0,
    Category::AlgorithmTolerance,
    "coarse bracket over log volume; refined by Newton, not by subdivision",
};

inline const Tolerance DENSITY_LOWER_FRACTION{
    1.0e-8,
    Category::AlgorithmTolerance,
    "Stage 2 brief; logarithmic density scan starts at a dilute ideal-gas fraction",
};

inline const Tolerance DENSITY_UPPER_FRACTION{
    0.999,
    Category::AlgorithmTolerance,
    "Stage 2 brief; logarithmic density scan stops below close packing",
};

inline const Tolerance DENSITY_BISECTION_FRACTION{
    0.5,
    Category::AlgorithmTolerance,
    "Stage 2 brief; safeguarded Newton uses the bracket midpoint",
};

inline const Tolerance ELECTROLYTE_SERIES_THRESHOLD{
    1.0e-3,
    Category::AlgorithmTolerance,
    "electrolyte.hpp finite-size chi series switch",
};

inline const Tolerance COMPOSITION_SUM_TOLERANCE{
    1.0e-12,
    Category::AlgorithmTolerance,
    "carried composition validation threshold from electrolyte_reference.hpp",
};

inline const Tolerance REFERENCE_MOLALITY_ONE{
    1.0e-6,
    Category::AlgorithmTolerance,
    "electrolyte_reference.hpp molality ladder",
};
inline const Tolerance REFERENCE_MOLALITY_TWO{
    1.0e-8,
    Category::AlgorithmTolerance,
    "electrolyte_reference.hpp molality ladder",
};
inline const Tolerance REFERENCE_MOLALITY_THREE{
    1.0e-10,
    Category::AlgorithmTolerance,
    "electrolyte_reference.hpp molality ladder",
};
inline const Tolerance REFERENCE_MOLALITY_FOUR{
    1.0e-12,
    Category::AlgorithmTolerance,
    "electrolyte_reference.hpp molality ladder",
};
inline const Tolerance REFERENCE_CONVERGENCE_WITNESS{
    1.0e-14,
    Category::AlgorithmTolerance,
    "electrolyte_reference.hpp convergence witness",
};
inline const Tolerance REFERENCE_CONVERGENCE_TOLERANCE{
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
inline const Tolerance STABILITY_CANDIDATE_GAP{
    1.0e-2,
    Category::AlgorithmTolerance,
    "Pereira 2012 eps_b; drifted to 5e-2 in the superseded tree",
};

// Lower bound on a composition coordinate in the inner problem. The old tree
// scaled this by the feed (1e-10 f_k), which put solutions closer to the bound
// than the paper's trace treatment assumes.
inline const Tolerance STABILITY_COMPOSITION_FLOOR{
    1.0e-8,
    Category::AlgorithmTolerance,
    "Pereira 2012; drifted to 1e-10 f_k in the superseded tree",
};

inline const Tolerance STABILITY_MULTIPLIER_GAP{
    0.5,
    Category::AlgorithmTolerance,
    "Pereira 2012 eps_lambda",
};

inline const Tolerance STABILITY_COMPOSITION_GAP{
    1.0e-3,
    Category::AlgorithmTolerance,
    "Pereira 2012 eps_x",
};

inline const Tolerance STABILITY_PACKING_GAP{
    1.0e-3,
    Category::AlgorithmTolerance,
    "Pereira 2012 eps_eta",
};

// A budget with a name. The old tree's global search stopped at 3250 samples and
// reported that as a one-phase result, so a decision that ran out of work was
// indistinguishable from a decision that converged. Reaching this is a named
// stop the caller can see, not an answer.
inline const Tolerance STABILITY_SEARCH_BUDGET{
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
inline const Tolerance SCREENING_SERIES_SWITCH{
    1.0e-3,
    Category::AlgorithmTolerance,
    "cancellation floor of the closed-form chi; carried from the superseded tree",
};

}  // namespace epcsaft::policy
