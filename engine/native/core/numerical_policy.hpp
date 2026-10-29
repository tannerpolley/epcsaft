#pragma once

#include "core/physical_constants.hpp"

#include <limits>
#include <array>
#include <cmath>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>

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
// unresolved accuracy are explicit.
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
    "Newton reaches it in 3-8 steps from a warm start",
};

inline constexpr Tolerance ASSOCIATION_ACTION_RESIDUAL{"ASSOCIATION_ACTION_RESIDUAL",
    1.0e-10, Category::AlgorithmTolerance,
    "approved P4 normalized mass-action coefficient residual; exceeding it makes the action unavailable"};

// Strong association starts from X = 1 far above the root, where a full Newton step leaves (0, 1].
inline constexpr Tolerance ASSOCIATION_BOUNDARY_FRACTION{"ASSOCIATION_BOUNDARY_FRACTION", 0.9,
    Category::AlgorithmTolerance, "fraction-to-boundary: a damped Newton step keeps X_a >= 0.1 X_a"};
// A budget, not a convergence criterion. Newton with the analytic Jacobian
// converges in single digits; reaching this many steps means the step selection
// is failing and the caller should hear about it rather than get a wrong answer
// from a truncated solve.
inline constexpr Tolerance ASSOCIATION_MAX_ITERATIONS{"ASSOCIATION_MAX_ITERATIONS",
    50.0,
    Category::AlgorithmTolerance,
    "budget; measured 3-8 Newton steps",
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
    "the reduced residual; publish the bracket endpoint with smaller |P_calc-P_target|; "
    "relative pressure error is |r|/Z_target, with Z_target=P_target/(rho R T)"
);

inline constexpr Tolerance DENSITY_MAX_ITERATIONS{"DENSITY_MAX_ITERATIONS",
    60.0,
    Category::AlgorithmTolerance,
    "budget for each root refinement, spinodal bisection and slope-minimum search",
};

// The coarse sweep that brackets the roots before Newton refines them; a bracket
// plus Newton finds the root in single digits. This is the scan's only
// discretisation number.
inline constexpr Tolerance DENSITY_BRACKET_INTERVALS{"DENSITY_BRACKET_INTERVALS",
    40.0,
    Category::AlgorithmTolerance,
    "coarse log-density scan; slope zeros separate stable brackets before safeguarded Newton",
};

inline constexpr Tolerance DENSITY_LOWER_FRACTION{"DENSITY_LOWER_FRACTION",
    1.0e-8,
    Category::AlgorithmTolerance,
    "fraction of the ideal-gas density P/RT where the log scan starts; a root lies above it unless Z > 1e8",
};

inline constexpr Tolerance DENSITY_UPPER_FRACTION{"DENSITY_UPPER_FRACTION",
    0.999,
    Category::AlgorithmTolerance,
    "fraction of the model's close-packing molar density where the log scan stops",
};

inline constexpr Tolerance DENSITY_BISECTION_FRACTION{"DENSITY_BISECTION_FRACTION",
    0.5,
    Category::AlgorithmTolerance,
    "bracket midpoint of the safeguarded Newton and the fold bisection",
};
inline constexpr Tolerance DENSITY_ANCHOR_LOG_WINDOW{"DENSITY_ANCHOR_LOG_WINDOW", 0.009950330853168083,
    Category::AlgorithmTolerance, "approved P2 local density envelope log1p(0.01)"};
inline constexpr Tolerance DENSITY_ANCHOR_SCAN_INTERVALS{"DENSITY_ANCHOR_SCAN_INTERVALS", 5.0,
    Category::AlgorithmTolerance, "approved P2 centered local anchor uniqueness scan"};
inline constexpr Tolerance DENSITY_ANCHOR_REGULARITY{"DENSITY_ANCHOR_REGULARITY",
    DENSITY_ANCHOR_LOG_WINDOW.value / 100.0, Category::AlgorithmTolerance,
    "admitted relative density error of an anchored root, one hundredth of DENSITY_ANCHOR_LOG_WINDOW"};
// The structured search's checks of the single mechanical loop (ALG-DENSITY).
inline constexpr Tolerance DENSITY_SLOPE_SEPARATION{"DENSITY_SLOPE_SEPARATION", 1.0e-6, Category::AlgorithmTolerance,
    "relative density separation above which two samples' P_rho order is inspected; not a rounding bound"};
// ALG-DENSITY: q=RT(1+u), c=-RT(v+w), w=2u/rho imply slope norm <= 2RT+|q|
// and curvature norm <= |c|+4(RT+|q|)/rho. Higham (1993), pp. 784-785, gamma_n:
// 9 local conversion factors + 5 norm/bound factors + 2 comparison factors, using epsilon
// conservatively in place of unit roundoff. Conditional on Taylor coefficients and normal arithmetic;
// this does not bound the complete EOS or implicit association error.
inline constexpr Tolerance DENSITY_STRUCTURE_ROUNDOFF{"DENSITY_STRUCTURE_ROUNDOFF",
    16.0 * std::numeric_limits<double>::epsilon() / (1.0 - 16.0 * std::numeric_limits<double>::epsilon()),
    Category::AlgorithmTolerance, "ALG-DENSITY local derivative assembly and comparison gamma_16"};


inline constexpr Tolerance COMPOSITION_SUM_TOLERANCE{"COMPOSITION_SUM_TOLERANCE",
    1.0e-12,
    Category::AlgorithmTolerance,
    "absolute defect of a unit sum or charge balance of O(1) terms, about 4500 ulp at 1",
};

// #214 NaCl scalar at 298.15 K, 1 bar meets 1e-4 in Phi; below floor unresolved, no error bound.
inline constexpr Tolerance OSMOTIC_MOLALITY_FLOOR{"OSMOTIC_MOLALITY_FLOOR",
    1.0e-8,
    Category::AlgorithmTolerance,
    "owner decision #214: common-pressure scalar cancellation reporting floor, mol/kg",
};

inline constexpr Tolerance EQUILIBRIUM_COMPOSITION_BOX{"EQUILIBRIUM_COMPOSITION_BOX",
    1.0e-8,
    Category::AlgorithmTolerance,
    "Pereira 2012 inner-problem composition bound; a regularising barrier, not a solution domain or accuracy bound",
};
// Largest fraction of a species' attainable inventory u_i the box may hold in the phases: the
// box is min(BOX, this u_i/N_feed), unchanged wherever the uniform box displaces less. One
// decade below the 10 % at which the retained trace-Li probe still polished strictly (#172).
inline constexpr Tolerance EQUILIBRIUM_BOX_INVENTORY_FRACTION{"EQUILIBRIUM_BOX_INVENTORY_FRACTION",
    1.0e-2,
    Category::AlgorithmTolerance,
    "largest fraction of a species' attainable inventory held by the composition box; not an accuracy bound",
};
// Lower bound of the inventory-scaled box of all-free cells: it binds only below an attainable
// inventory of 1e-18 N_feed, under the balance rows' rounding resolution. No root, action,
// continuation step or seed is cut at it (#266, #322): x_B = K/(1+K) = 1e-24 solves, has exact
// actions and continues. Feasibility cells are bounded by the smallest normal double instead.
inline constexpr Tolerance EQUILIBRIUM_BOX_FLOOR{"EQUILIBRIUM_BOX_FLOOR",
    1.0e-20,
    Category::AlgorithmTolerance,
    "dimensionless lower bound of the all-free optimizer log-composition box; not a solution domain or accuracy bound",
};
inline constexpr Tolerance PHASE_DISAPPEARANCE_AMOUNT{"PHASE_DISAPPEARANCE_AMOUNT",
    1.0e-8,
    Category::AlgorithmTolerance,
    "dimensionless n_phase/feed_total; provisional free-phase disappearance topology cutoff and interior guard of the "
    "split seed; not a stability result or accuracy bound",
};
inline constexpr Tolerance EQUILIBRIUM_AMOUNT_STEP_RETENTION{"EQUILIBRIUM_AMOUNT_STEP_RETENTION",
    0.1,
    Category::AlgorithmTolerance,
    "dimensionless kappa_N: an optimizer trial point keeps at least kappa_N of each free phase amount's distance to "
    "its relaxed lower bound at the search base (ALG-EQUILIBRIUM); provisional numerical choice, not an amount floor, "
    "stability result or accuracy bound",
};

// ------------------------------------------------------------- regression

// Also the conservation test of ALG-STABILITY's amount correction, finite split and phase removal,
// whose component and moiety errors are divided by max(1 mol, N_feed).
inline constexpr Tolerance ADDITION_MATERIAL_BALANCE{"ADDITION_MATERIAL_BALANCE", 1.0e-12,
    Category::AlgorithmTolerance, "#210: component feed balance divided by max(1 mol, endpoint feed amounts); "
    "admission of two declared endpoints and of stability inventory conservation, independent of equilibrium stopping error"};

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
// A coordinate within this many of its own scale units of a bound is on it and is held while its gradient
// points out of the box: an interior start a rounding step from a bound does not stall the free step.
inline constexpr Tolerance FIT_ACTIVE_BOUND_DISTANCE{"FIT_ACTIVE_BOUND_DISTANCE",
    1.0e-4,
    Category::AlgorithmTolerance,
    "optimizer-coordinate distance to a bound below which the coordinate is on the bound",
};
// Inside |y - y_obs| < LOG_RATIO_NEAR_UNITY * y_obs the log-ratio residual is evaluated as
// log1p((y - y_obs) / y_obs), which keeps its leading digits near unity; outside, as the
// difference of logarithms, which stays finite across the double range.
inline constexpr Tolerance LOG_RATIO_NEAR_UNITY{"LOG_RATIO_NEAR_UNITY",
    0.5,
    Category::AlgorithmTolerance,
    "relative half-band around unity inside which ln(y / y_obs) is evaluated by log1p",
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

// ALG-FIT/#245: maximum physical-chart change for an optional density initializer.
inline constexpr Tolerance FIT_DENSITY_PARAMETER_FRACTION{"FIT_DENSITY_PARAMETER_FRACTION", 0.10,
    Category::AlgorithmTolerance, "abs(new-old)/max(1,abs(new),abs(old)); includes held diameters"};

// ------------------------------------------------------------- equilibrium

__attribute__((visibility("default"))) inline constexpr Tolerance EQUILIBRIUM_ROW_RESIDUAL{"EQUILIBRIUM_ROW_RESIDUAL",
    1.0e-10,
    Category::AlgorithmTolerance,
    "dimensionless equilibrium-row infinity norm; two decades above the "
    "association residual entering the assembly floor",
};
// #811: empirical retention policy for trial polish; ln(1 + COMPOSITION_BOX / ROW_RESIDUAL) = ln 101.
inline constexpr Tolerance EQUILIBRIUM_TRIAL_LOG_STEP{"EQUILIBRIUM_TRIAL_LOG_STEP",
    4.61512051684126, Category::AlgorithmTolerance,
    "largest log mole fraction or log packing change per optimized-trial polish proposal; "
    "an empirical retention policy, not an accuracy bound"};
inline constexpr Tolerance EQUILIBRIUM_ACTION_BACKWARD_ERROR{"EQUILIBRIUM_ACTION_BACKWARD_ERROR", 1.0e-10,
    Category::AlgorithmTolerance, "normwise backward error of each implicit Taylor solve; catches a stale "
    "Jacobian or jet mismatch"};
inline constexpr Tolerance EQUILIBRIUM_ACTION_FORWARD_ERROR{"EQUILIBRIUM_ACTION_FORWARD_ERROR", 1.0e-6,
    Category::AlgorithmTolerance, "first-order relative forward-error bound 2 cond_inf(J) eta of each Taylor "
    "coefficient; owner decision 2026-09-23 (#87)"};

// The optimizer works below the caller's row criterion: the recomputed row maximum tracks
// the Ipopt tolerance and overshoots it by a hair (measured 1.0002-1.02x), so classifying
// at the same value is systematically marginal. Inner Ipopt tolerances are
// the requested row tolerance times this factor; the caller criterion alone
// classifies rows. Owner decision 2026-09-16 (E4); provisional, not an accuracy bound. Pinned
// square systems use the caller criterion itself; their Newton polish absorbs the overshoot
// (owner decision 2026-09-23).
inline constexpr Tolerance EQUILIBRIUM_INNER_TOLERANCE_FACTOR{"EQUILIBRIUM_INNER_TOLERANCE_FACTOR",
    0.1,
    Category::AlgorithmTolerance,
    "provisional decade below the caller row criterion for the Ipopt tol, acceptable_tol, "
    "constr_viol_tol and the objective-normalized dual stop; the caller criterion classifies rows",
};

inline constexpr Tolerance EQUILIBRIUM_POLISH_STEPS{"EQUILIBRIUM_POLISH_STEPS",
    3.0,
    Category::AlgorithmTolerance,
    "raw-contracting Newton steps at positive compositions, more while raw rows exceed the caller "
    "criterion; return the lowest classified residual among the optimizer start and accepted steps",
};

inline constexpr Tolerance TWO_DENSITY_FORWARD_ERROR{"TWO_DENSITY_FORWARD_ERROR", 1e-10,
    Category::AlgorithmTolerance,
    "#854 owner-approved log-density error target; exceeding it engages wider evaluation"};

inline constexpr Tolerance MIXTURE_FORWARD_ERROR{"MIXTURE_FORWARD_ERROR", 1e-7,
    Category::AlgorithmTolerance,
    "#864 physical relative pressure/log-density, absolute mole-fraction and feed-scaled amount error target; exceeding it engages wider evaluation"};

inline constexpr Tolerance EQUILIBRIUM_POLISH_CORRECTION{"EQUILIBRIUM_POLISH_CORRECTION",
    1.0e-12, Category::AlgorithmTolerance,
    "#272 scaled undamped Newton correction and forward-error estimate; local termination, not a state-error certificate",
};
// ALG-EQUILIBRIUM polish (#633): a moiety no larger than this fraction of the feed, by inventory and by the
// terms its balance row sums, is judged relative to that size. The threshold of #272's trace-species definition.
inline constexpr Tolerance EQUILIBRIUM_TRACE_MOIETY_FRACTION{"EQUILIBRIUM_TRACE_MOIETY_FRACTION",
    1.0e-6, Category::AlgorithmTolerance,
    "largest moiety size, per mol of feed, whose balance row is judged relative to that size; not an accuracy bound",
};
// A reacting one-phase start is the ideal-solution speciation of its feed with frozen nonideal terms
// (ALG-EQUILIBRIUM, #154/#157); its element-potential Newton iteration stops when every log-amount step,
// trace species included, is below this. The frozen terms bound the start's accuracy, not this value.
inline constexpr Tolerance EQUILIBRIUM_START_SPECIATION_STEP{"EQUILIBRIUM_START_SPECIATION_STEP",
    1.0e-6,
    Category::AlgorithmTolerance,
    "largest log-amount Newton step at which the ideal start speciation stops; a start, not an accuracy bound",
};

inline constexpr Tolerance EQUILIBRIUM_START_PRESSURE_STEP{"EQUILIBRIUM_START_PRESSURE_STEP", 1.0e-3, Category::AlgorithmTolerance, "boundary start step; not an accuracy bound"},
    EQUILIBRIUM_START_PRESSURE_PASSES{"EQUILIBRIUM_START_PRESSURE_PASSES", 8.0, Category::AlgorithmTolerance, "boundary start passes; not an accuracy bound"};
inline constexpr Tolerance EQUILIBRIUM_MAX_ITERATIONS{"EQUILIBRIUM_MAX_ITERATIONS",
    200.0,
    Category::AlgorithmTolerance,
    "budget; analytic equilibrium Jacobian with Ipopt",
};
// ALG-STABILITY (#647): direct substitutions before second-order steps in a Michelsen tangent-plane trial.
inline constexpr Tolerance TANGENT_PLANE_SUBSTITUTIONS{"TANGENT_PLANE_SUBSTITUTIONS", 5.0, Category::AlgorithmTolerance,
    "count; Michelsen 1982 eq. 19 steps before eqs. 25-27, as in FeOs stability_analysis"};
// ALG-EQUILIBRIUM neutral recovery (#565): sources at T_seed r^j, j = 1..8; straight continuation steps 0.05 to 1e-6.
inline constexpr Tolerance COOLING_RATIO{"COOLING_RATIO", 0.95, Category::AlgorithmTolerance,
    "T_(j+1)/T_j of #565 recovery sources and #610 PH/PS starts"};
inline constexpr Tolerance EQUILIBRIUM_RECOVERY_SOURCES{"EQUILIBRIUM_RECOVERY_SOURCES", 8.0, Category::AlgorithmTolerance,
    "#565 lead-approved finite source schedule, ten retries per continuation step"};
inline constexpr Tolerance EQUILIBRIUM_RECOVERY_WORK{"EQUILIBRIUM_RECOVERY_WORK", 512.0, Category::AlgorithmTolerance,
    "#565 lead-approved added iterations and continuation attempts; 4x evaluations, 2x Hessian callbacks"};

inline constexpr Tolerance EQUILIBRIUM_CALORIC_SEED_ITERATIONS{"EQUILIBRIUM_CALORIC_SEED_ITERATIONS", 20.0,
    Category::AlgorithmTolerance, "budget of the discardable PH/PS fixed-T seed; converged seeds took 4-14 (#579)"};

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

// Floor of the log-packing chart, shared by the compiled preparation's domain test, the
// optimizer bounds of multi-phase solves and feasibility cells, and the continuation
// predictor's domain test so all admit the same coordinates.
inline constexpr Tolerance PACKING_FRACTION_FLOOR{"PACKING_FRACTION_FLOOR",
    1.0e-12,
    Category::AlgorithmTolerance,
    "fluid-coordinate floor of the log-packing chart",
};

// Below this reduced screening length the closed form of the finite-size correction
// loses precision to cancellation -- log(1+r) - r + r^2/2 is a difference of terms that
// agree to leading order -- so the series through r^24 is used instead. Against 50-digit
// references (2026-09-23), the worst relative errors of chi and its first and second
// Taylor derivatives are 3e-14, 7e-13 and 1.1e-11 on either side of 0.2; the former
// switch at 1e-3 with six series terms left the closed form's second derivative wrong
// by a factor of 3.3 just above it. Both branches are on the tape behind a conditional
// expression, so one tape serves either side.
inline constexpr Tolerance SCREENING_SERIES_SWITCH{"SCREENING_SERIES_SWITCH",
    0.2,
    Category::AlgorithmTolerance,
    "crossover of closed-form cancellation and the 24-term series truncation in chi and its derivatives",
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

inline constexpr Tolerance EQUILIBRIUM_WARM_PUSH{"EQUILIBRIUM_WARM_PUSH", 1.0e-6,
    Category::AlgorithmTolerance, "Ipopt 3.11.9 manual, Warm-Starting Capabilities via AMPL; #505"};
inline constexpr Tolerance EQUILIBRIUM_WARM_MU{"EQUILIBRIUM_WARM_MU", 1.0e-6,
    Category::AlgorithmTolerance, "Ipopt 3.11.9 manual, monotone initial barrier for similar active sets; #505"};
// Every entry, readable as epcsaft.tolerances(); the architecture check requires listing and a consumer.
inline constexpr const Tolerance* ENTRIES[] = {
    &EQUILIBRIUM_WARM_PUSH, &EQUILIBRIUM_WARM_MU,
    &CLOSE_PACKING_FRACTION, &SITE_FRACTION_FLOOR, &SITE_FRACTION_CEILING, &ASSOCIATION_RESIDUAL,
    &ASSOCIATION_ACTION_RESIDUAL, &ASSOCIATION_MAX_ITERATIONS, &ASSOCIATION_BOUNDARY_FRACTION,
    &ASSOCIATION_JACOBIAN_CONDITION, &ASSOCIATION_JACOBIAN_RANK, &DENSITY_ROOT_RESIDUAL, &DENSITY_MAX_ITERATIONS,
    &DENSITY_BRACKET_INTERVALS, &DENSITY_LOWER_FRACTION, &DENSITY_UPPER_FRACTION, &DENSITY_BISECTION_FRACTION,
    &DENSITY_ANCHOR_LOG_WINDOW, &DENSITY_ANCHOR_SCAN_INTERVALS, &DENSITY_ANCHOR_REGULARITY, &DENSITY_SLOPE_SEPARATION, &DENSITY_STRUCTURE_ROUNDOFF,
    &COMPOSITION_SUM_TOLERANCE,
    &OSMOTIC_MOLALITY_FLOOR,
    &EQUILIBRIUM_COMPOSITION_BOX, &EQUILIBRIUM_BOX_INVENTORY_FRACTION,
    &EQUILIBRIUM_BOX_FLOOR,
    &PHASE_DISAPPEARANCE_AMOUNT, &EQUILIBRIUM_AMOUNT_STEP_RETENTION, &FIT_FUNCTION_TOLERANCE, &FIT_GRADIENT_TOLERANCE, &FIT_PARAMETER_TOLERANCE, &FIT_ACTIVE_BOUND_DISTANCE,
    &LOG_RATIO_NEAR_UNITY, &FIT_RANK_RECIPROCAL_CONDITION, &FIT_DENSITY_PARAMETER_FRACTION, &ADDITION_MATERIAL_BALANCE, &EQUILIBRIUM_ROW_RESIDUAL,
    &EQUILIBRIUM_TRIAL_LOG_STEP, &EQUILIBRIUM_ACTION_BACKWARD_ERROR, &EQUILIBRIUM_ACTION_FORWARD_ERROR, &EQUILIBRIUM_INNER_TOLERANCE_FACTOR,
    &EQUILIBRIUM_POLISH_STEPS, &TWO_DENSITY_FORWARD_ERROR, &MIXTURE_FORWARD_ERROR, &EQUILIBRIUM_POLISH_CORRECTION, &EQUILIBRIUM_TRACE_MOIETY_FRACTION, &EQUILIBRIUM_START_SPECIATION_STEP, &EQUILIBRIUM_START_PRESSURE_STEP, &EQUILIBRIUM_START_PRESSURE_PASSES, &EQUILIBRIUM_MAX_ITERATIONS, &TANGENT_PLANE_SUBSTITUTIONS, &EQUILIBRIUM_CALORIC_SEED_ITERATIONS, &COOLING_RATIO, &EQUILIBRIUM_RECOVERY_SOURCES, &EQUILIBRIUM_RECOVERY_WORK, &PHASE_COINCIDENCE_DISTANCE, &PACKING_FRACTION_CEILING,
    &PACKING_FRACTION_FLOOR, &SCREENING_SERIES_SWITCH, &FLOOR_MARGIN, &STALL_WINDOW};

// Six stopping controls; defaults remain the policies above. Presence survives inheritance.
enum SolverOption : std::size_t { association, density, row, function, gradient, parameter };
inline constexpr std::array<const Tolerance*, 6> SOLVER_POLICIES = {&ASSOCIATION_RESIDUAL,
    &DENSITY_ROOT_RESIDUAL, &EQUILIBRIUM_ROW_RESIDUAL, &FIT_FUNCTION_TOLERANCE,
    &FIT_GRADIENT_TOLERANCE, &FIT_PARAMETER_TOLERANCE};
inline constexpr std::array<const char*, 6> SOLVER_UNITS = {"1 (mass action)", "1 (reduced pressure)",
    "1 (classified rows)", "1 (relative cost)", "scaled objective / optimizer coordinate", "1 (relative step)"};

struct SolverOptions final {
    std::array<std::optional<double>, 6> values{};
    std::array<std::shared_ptr<const Tolerance>, 6> owners{};

    double value(std::size_t i) const { return values.at(i).value_or(SOLVER_POLICIES.at(i)->value); }
    const Tolerance* criterion(std::size_t i) const { return owners.at(i) ? owners[i].get() : SOLVER_POLICIES.at(i); }
    static void validate(std::size_t i, double v) {
        if (!std::isfinite(v) || !(v > 0.0))
            throw std::invalid_argument(std::string(SOLVER_POLICIES.at(i)->name) + "=" + std::to_string(v)
                + " [" + SOLVER_UNITS.at(i) + "] requires 0 < value <= double max, finite");
    }
    void set(std::size_t i, double v) {
        validate(i, v);
        if (values.at(i)) throw std::invalid_argument(std::string("duplicate ") + SOLVER_POLICIES[i]->name);
        values[i] = v;
    }
    void set(std::string_view name, double v, unsigned scope = 63) {
        for (std::size_t i = 0; i < values.size(); ++i) if (name == SOLVER_POLICIES[i]->name) {
            if (!(scope & (1u << i))) throw std::invalid_argument(std::string(name) + " outside owner scope");
            set(i, v); return;
        }
        throw std::invalid_argument(std::string(name) + " is unknown or fixed");
    }
    SolverOptions resolved(const SolverOptions& inherited) const {
        SolverOptions result = inherited;
        for (std::size_t i = 0; i < values.size(); ++i) {
            if (values[i]) result.values[i] = values[i];
            const double v = result.value(i);
            validate(i, v);
            if (v == SOLVER_POLICIES[i]->value) result.owners[i].reset();
            else if (!result.owners[i] || result.owners[i]->value != v) {
                auto criterion = *SOLVER_POLICIES[i];
                criterion.value = v; criterion.source = "caller-supplied stopping tolerance";
                criterion.upstream = nullptr; criterion.factor = 1.0;
                result.owners[i] = std::make_shared<const Tolerance>(criterion);
            }
        }
        return result;
    }
    bool operator==(const SolverOptions& other) const {
        for (std::size_t i = 0; i < values.size(); ++i) if (value(i) != other.value(i)) return false;
        return true;
    }
};

}  // namespace epcsaft::policy
