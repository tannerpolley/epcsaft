# Parameter regression formulation

Required general fitting mathematics; pure saturation is one observation
construction. Installed availability is recorded in the
[GitHub issue roadmap](https://github.com/tannerpolley/ePC-SAFT/milestones).

Use the plan's adoption gate to adapt justified existing native fitting and
evaluator logic. Implement only observation forms consumed by current scope;
the broader algebra below defines meaning rather than requiring unused nodes.

An observation combines a state/equilibrium specification, an expression, a
reported value and an explicit discrepancy/uncertainty model. Expressions use
four node kinds: Leaf (property, coordinate or constant), weighted Sum,
Product with declared exponents, and Log. This algebra covers partial pressure,
finite-dose heat, extraction fractions and ratios of selectivities without
application-specific property code. Every operation must respect its domain:
positive log arguments, nonzero denominators and valid real powers.

## Built-in observations (R1)

The implemented inventory is finite; each kind names one shared owner and one
predicted quantity in SI molar units. Phases and species are indices into the
declared problem; a mass basis multiplies by the fixed molar masses.

| Kind | Fixed inputs | Owner and value | Unit |
|---|---|---|---|
| Density | $T$, $P$, normalized $\mathbf{x}$, declared branch | pressure closure `state_at_pressure`; the resolved branch and stable-root count are retained (one stable root is returned as Unspecified) | mol m$^{-3}$ or kg m$^{-3}$ |
| PhaseDensity | a declared specification, phase $p$ | the driver's phase molar density | mol m$^{-3}$ or kg m$^{-3}$ |
| EquilibriumPressure | $T$ fixed, one free and one pinned-zero amount (bubble, dew, pure saturation) | the driver's free pressure | Pa |
| EquilibriumTemperature | $P$ fixed, the same phase specification with $T$ free; a cold guess must lie above the answer (a guess below reports a coalesced-phase topology event) | the driver's free temperature | K |
| PhaseMoleFraction | any declared specification, phase $p$, nonempty species set $S$ within its support | $\sum_{i\in S} x_{p,i}$ | 1 |
| PhaseFraction | free amounts | $\beta_p = n_p / \sum_q n_q$ over the solved amounts, never the feed total | 1 |
| PhaseAmountRatio | declared numerator phase/species set $(p,S)$ and denominator $(q,S')$ | $A/B$ with $A=\sum_S n_{p,i}$, $B=\sum_{S'}n_{q,i}$ and $B\ne0$ | 1 |
| PartialPressure | any declared specification, phase $p$, nonempty species set $S$ within its support | $P\sum_{i\in S} x_{p,i}$ with the solved system pressure; `LogRatio` gives $\ln(p_S/p_{S,obs})$ | Pa |

An unavailable prediction is typed, never a zero residual: a vanished or
coalesced phase is `topology_event` before any solve status; a rejected or
state-less solve is `central_solve_unavailable` with its convergence record;
a missing density root is `density_unavailable` with the thrown cause; a mass
basis without a recorded ionic molar mass is `missing_molar_mass`.

Residual forms: `Difference` $r=(y-y_{obs})/s$ with a positive
nondimensionalisation scale $s$ in the observation's unit (not an
uncertainty), and `LogRatio` $r=\ln(y/y_{obs})$ for positive quantities, which
carries no scale, evaluated as $\ln y-\ln y_{obs}$ across the dynamic range
and as $\mathrm{log1p}$ inside the named `LOG_RATIO_NEAR_UNITY` band; a nonpositive log argument is
`nonpositive_log_argument` and a residual that overflows is `nonfinite`
while the independently valid value is retained. A mass basis is defined only
for the two density kinds and is refused by name elsewhere.
A constant factor multiplying an amount ratio is folded into its observed value;
a Difference residual must scale its nondimensionalization scale by that same
factor. Phase volumes are not inferred.
The two forms are not commensurable by themselves; relative weights
$w_i\ge 0$ enter the objective as $\sqrt{w_i}\,r_i$ and are declared
separately from any absolute uncertainty.

## Parameter coordinates (R1)

An EOS physical coordinate is exactly the Properties active-parameter coordinate:
one component slot for $m$, $\sigma$ and $\varepsilon/k$, one reciprocal pair
for $k_{ij}$ and $l_{ij}$, one site pair or diagonal self-pair for
$\varepsilon^{AB}/k$ and $\kappa^{AB}$, one reciprocal pair with a positive
reference temperature for the $k_{ij}$ slope $b_{ij}$ of EqID:
`cross_dispersion_energy`, in native units (1, Å, K). Every slot
absent from the active list is fixed. An optimizer coordinate is affine,
$\theta = o + s\,u$ with $s>0$, with inclusive bounds in physical units mapped
by the same transform; one optimizer coordinate may drive several tied slots of
the same family through the fitting layer's linear map $\theta = o + S u$; tied
values must agree at the start. Physical rows contain all EOS slots in coordinate declaration order, followed by
one row per reaction coordinate in declaration order. The model's active list
contains the EOS rows only. Optimizer starts, bounds, columns and covariance
keep coordinate declaration order; the tie map relates those two orders.

A reaction coordinate is a named scalar ln K (coefficient a) or a stored
`ReactionLogPolynomial` coefficient a, b, c, d, in units 1, K, 1, K⁻¹. It uses
the same affine map and inclusive bounds. Name tying requires identical
stoichiometry, constant representation and starts, correlation coefficients,
Tref, interval and standard-state identity, and source basis/pressure. Referenced
matches also share the solvent and common-molality standard; their dilution-path
ratios may differ. The derivative is EqID:
`reaction_constant_coordinate_derivative` evaluated at the solved temperature.
Source/reference declarations are preserved in every trial. Scalar constants
are rewritten as exp(ln K); correlations keep their native stored coefficients.
Affinity targets, supplied Gibbs values, dependent reaction sets, thermochemistry-
fixed unreferenced constants and referenced free-temperature requests cannot be
fitted. No coefficient Hessians or caloric coefficient actions are claimed.

## Prediction Jacobians and fitting (R2, R3)

$\partial y/\partial\theta$ comes from the shared owners only: the Properties
pressure-chart density action for Density and the solved-state actions of the
equilibrium driver for every other kind (a linear sum over $S$, the quotient
rule for $\beta_p$ and amount ratios ($d(A/B)=(dA\,B-A\,dB)/B^2$), the product rule for a mass basis and for
$\mathrm{d}(P\sum_S x) = \mathrm{d}P\sum_S x + P\sum_S\mathrm{d}x$). One
first-action batch per state carries every declared observable and active
direction at one factorization of the solved Jacobian, reusing the problem compiled for the solve; a referenced
problem builds its standard-state parameter jets once for all directions, since
they enter the rows only through $\sum_a (\partial s/\partial\theta_a)(\theta_a-\theta_{a,c})$. The optimizer
Jacobian is $\sqrt{w}\,(\mathrm{d}r/\mathrm{d}y)\,(\partial y/\partial\theta)\,S\,\mathrm{diag}(s)$
with $\mathrm{d}r/\mathrm{d}y = 1/s$ or $1/y$. The bounded fit is a
Levenberg--Marquardt trust region (Ceres, dense QR, one thread) over $u$ in a
chart of the box instead of Ceres' bounds (whose projected line search would
re-solve and differentiate all states at every trial and whose per-coordinate
clamp spoils a step computed as if a bound could be crossed): a step leaving the
box is shortened along its direction to the first bound (the trust-region ratio
keeps the full step's predicted decrease, so it is conservative). Bounds form an
active set in optimizer space, against the physical bounds mapped by the affine
transform (rounding is monotone, so a physical value inside its bounds is inside
these): a coordinate within the named `FIT_ACTIVE_BOUND_DISTANCE` of a bound
whose objective gradient $J^\top r$ points out of the box is on it and is held
where it is for a round, never moved onto the bound; a free coordinate on a
bound it moves toward is clamped instead of shortening the step. A round ends
when Ceres converges, an accepted step puts a free coordinate exactly on a
bound it was not on at the previous accepted point, or the held membership
selected by the full accepted weighted optimizer gradient $J^\top r$ changes,
including release when a held coordinate's gradient turns inward. Rejected
trials cannot change the set. The accepted-point Jacobian is reused without
another thermodynamic solve; a failed gradient stops with its typed cause.
Any other Ceres stop ends the fit. The set is rebuilt at that point, and the
fit converges when a converged round leaves the set unchanged or every
coordinate is held. Unchanged membership retains Ceres' function, parameter
and gradient stopping criteria; it does not alone certify a small free
projected gradient. Every coordinate held satisfies the bound-gradient signs
to the declared distance. A bound within the same
distance is active for the covariance, and physical values are clamped to their
bounds against the affine map's rounding. Iterations
and time are shared by the rounds; a
trial in which any observation or action is unavailable is rejected with its
typed cause and observation index retained, never zeroed. Exact raw equilibrium
declarations, including feed basis, phase supports and seeds, reaction/reference inputs and row tolerance, share one compiled central
solve within a fitting evaluation; projection, residual and weight metadata
remain independent. Density observations remain separate pressure closures.
Validation preparation is reused only at exactly matching resolved EOS values
and reaction declarations. Scalar K and all correlation coefficients belong
to that comparison; a changed reaction coordinate invalidates the compiled
problem even when EOS values match. The last fully successful optimizer point
retains its immutable model and central states: a subsequent equilibrium-observation Jacobian request at that
exact point adds actions without compiling or solving equilibrium again.
Density differentiation retains its pressure-closure work for each direction
batch. Weighted objective and optimizer-Jacobian
finiteness are checked before publishing the trial; failed promotion retains
only the successful residual trial. Warm starts are taken
only from a trial in which every observation succeeded and are retried cold on
failure. The conditional covariance is $(J^\top J)^{-1}$ of the weighted
optimizer Jacobian, mapped to physical rows by $D = S\,\mathrm{diag}(s)$; it is
withheld as `rank_deficient` when a singular value falls below the named
reciprocal-condition policy, as `active_bound` at a bound, and as
`not_converged` for an iterate the optimizer did not accept as converged (its
rank, singular values and bound diagnostics stay). Rank and the conditional
covariance describe local conditioning at the fitted point under the stated
assumptions; a converged multistart spread describes reproducibility across
those starts; neither establishes structural or global identifiability or
predictive accuracy. The variance
factor $\sum w r^2/(n_w-\mathrm{rank})$ over the $n_w$ rows with positive
weight is reported separately for an unknown common variance and is not
applied; a bound is detected in optimizer space against the exact values
handed to the optimizer. The objective is one half of the sum of squared weighted residuals; because the
terms are non-negative its overflow is independent of summation order, and the
reported row is where the row-order running sum stops being representable. That
row, or a row whose own weighted residual is not finite, is a typed `nonfinite`
trial failure at every evaluation, and no infinite or partial sum is reported as
a cost.
The result carries a `status` field (converged, limit_reached,
initial_evaluation_failed, failed) and the solver's own termination message,
which distinguishes an iteration limit from a time limit; a nonzero count of
typed trial failures is a diagnostic of rejected trials, not a failed fit.
Each Ceres iteration's wall time and cost (iteration 0 is the start) and each
failure's residual/Jacobian flag are reported, so a slow or failing fit is
visible from its result.

Fixed-coordinate leaves are valid constants with zero derivative. Observable
and parameter order, units, bases and chain rules must remain explicit through
physical and optimizer coordinates. The common state/equilibrium owner supplies
values and required actions; the fitting layer owns no duplicate EOS or equilibrium.
The [action table](../../science/algorithms.md#alg-actions) determines required
energy order; first residual Jacobians can require third energy derivatives
when the observed quantity is a heat capacity.

Fit the declared discrepancy under stated bounds, transforms, starts, loss and
numerical controls. Retain raw predictions, physical Jacobians, weighted
optimizer-coordinate Jacobians, failed row identities, conditioning and the
termination reason (the status and the solver's termination message). A failed or unsupported state must not become a successful zero
residual or silently disappear from the objective.

Covariance interpretation depends on the data model. With known absolute
independent Gaussian uncertainty and a regular locally linear fit, the inverse
weighted information matrix estimates parameter covariance. With an unknown
common variance estimated from residuals, a separately estimated variance
factor is required. Correlated errors, constraints, rank deficiency, nonlinear
losses and model discrepancy require their stated treatment; no universal
variance multiplier or unconditional covariance promise applies.

A manufactured identifiable case can verify parameter recovery. Agreement with
published property data does not imply unique recovery of a published physical
parameter tuple. Record identifiability, fitted-domain and holdout evidence;
parameter adoption remains an application decision. Use the retained
[literature reproduction classes](literature-reproduction-contract.md), including
the distinction between measured values, evaluated references and source model
curves. Historical Ceres callable names and stage line/time budgets are not
current interfaces or universal acceptance criteria.
