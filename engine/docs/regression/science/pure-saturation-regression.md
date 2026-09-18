# Regression formulation

Required general fitting mathematics; pure saturation is one observation
construction. Installed availability is recorded in the
[GitHub issue roadmap](https://github.com/tannerpolley/ePC-SAFT-greenfield/milestones).

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

An unavailable prediction is typed, never a zero residual: a vanished or
coalesced phase is `topology_event` before any solve status; a rejected or
state-less solve is `central_solve_unavailable` with its convergence record;
a missing density root is `density_unavailable` with the thrown cause; a mass
basis without a recorded ionic molar mass is `missing_molar_mass`.

Residual forms: `Difference` $r=(y-y_{obs})/s$ with a positive
nondimensionalisation scale $s$ in the observation's unit (not an
uncertainty), and `LogRatio` $r=\ln(y/y_{obs})$ for positive quantities, which
carries no scale; a nonpositive log argument is `nonpositive_log_argument`.
The two forms are not commensurable by themselves; relative weights
$w_i\ge 0$ enter the objective as $\sqrt{w_i}\,r_i$ and are declared
separately from any absolute uncertainty.

## Parameter coordinates (R1)

A physical coordinate is exactly the Properties active-parameter coordinate:
one component slot for $m$, $\sigma$ and $\varepsilon/k$, one reciprocal pair
for $k_{ij}$ and $l_{ij}$, one site pair or diagonal self-pair for
$\varepsilon^{AB}/k$ and $\kappa^{AB}$, in native units (1, Å, K). Every slot
absent from the active list is fixed. An optimizer coordinate is affine,
$\theta = o + s\,u$ with $s>0$, with inclusive bounds in physical units mapped
by the same transform; one optimizer coordinate may drive several tied slots of
the same family, a regression-owned linear map $\theta = o + S u$ whose tied
values must agree at the start. The model's active list is the coordinate
set's physical list in order, so every parameter direction of the shared
owners is indexed by these rows.

## Prediction Jacobians and fitting (R2, R3)

$\partial y/\partial\theta$ comes from the shared owners only: the Properties
pressure-chart density action for Density and the solved-state actions of the
equilibrium driver for every other kind (a linear sum over $S$, the quotient
rule for $\beta_p$, the product rule for a mass basis). The optimizer
Jacobian is $\sqrt{w}\,(\mathrm{d}r/\mathrm{d}y)\,(\partial y/\partial\theta)\,S\,\mathrm{diag}(s)$
with $\mathrm{d}r/\mathrm{d}y = 1/s$ or $1/y$. The bounded fit is a
Levenberg--Marquardt trust region (Ceres, dense QR, one thread) over $u$; a
trial in which any observation or action is unavailable is rejected with its
typed cause and observation index retained, never zeroed. Warm starts are taken
only from a trial in which every observation succeeded and are retried cold on
failure. The conditional covariance is $(J^\top J)^{-1}$ of the weighted
optimizer Jacobian, mapped to physical rows by $D = S\,\mathrm{diag}(s)$; it is
withheld as `rank_deficient` when a singular value falls below the named
reciprocal-condition policy, and as `active_bound` at a bound. The variance
factor $\sum w r^2/(n_w-\mathrm{rank})$ over the $n_w$ rows with positive
weight is reported separately for an unknown common variance and is not
applied; a bound is detected in optimizer space against the exact values
handed to the optimizer.

Fixed-coordinate leaves are valid constants with zero derivative. Observable
and parameter order, units, bases and chain rules must remain explicit through
physical and optimizer coordinates. The common state/equilibrium owner supplies
values and required actions; Regression owns no duplicate EOS or equilibrium.
The [action table](../../science/algorithms.md#alg-actions) determines required
energy order; first residual Jacobians can require third energy derivatives
when the observed quantity is a heat capacity.

Fit the declared discrepancy under stated bounds, transforms, starts, loss and
numerical controls. Retain raw predictions, physical Jacobians, weighted
optimizer-coordinate Jacobians, failed row identities, conditioning and
termination. A failed or unsupported state must not become a successful zero
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
