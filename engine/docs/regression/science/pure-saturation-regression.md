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
finite-addition enthalpy, extraction fractions and ratios of selectivities without
application-specific property code. Every operation must respect its domain:
positive log arguments, nonzero denominators and valid real powers.

## Built-in observations (R1)

The implemented inventory is finite; each kind names one shared owner and one
predicted quantity in SI molar units. Phases and species are indices into the
declared problem; a mass basis multiplies by the fixed molar masses.

| Kind | Fixed inputs | Owner and value | Unit |
|---|---|---|---|
| Density | $T$, $P$, normalized $\mathbf{x}$, declared branch | pressure closure `state_at_pressure`; the resolved branch and stable-root count are retained (one stable root is returned as Unspecified) | mol m$^{-3}$ or kg m$^{-3}$ |
| LogEquilibriumConstant | fixed positive `problem.temperature.value`, named `reaction`, declared source record | `evaluate_reaction_properties` for a/b/c/d, or ln(scalar K); no equilibrium solve | 1 |
| IsobaricHeatCapacity | $T$, $P$, normalized $\mathbf{x}$, declared branch, fixed thermochemistry record | pressure closure `state_at_pressure`; total molar $c_p$ | J mol$^{-1}$ K$^{-1}$ |
| EquilibriumIsobaricHeatCapacity | driven $T/P$, declared equilibrium/feed, fixed record thermochemistry and source references | Equilibrium's extensive H thermal A1 divided by the feed sum and thermal increment (EqID: `equilibrium_isobaric_heat_capacity`) | J mol-feed$^{-1}$ K$^{-1}$ |
| AdditionEnthalpy | declared initial and final `Problem`s, positive `feed_amount` $n_{add}$, and feed `temperature`, `pressure`, `composition` $\mathbf y$, and `branch` | Parameter regression; total-system enthalpy increment less feed enthalpy (EqID: `finite_addition_enthalpy`) | J mol$^{-1}$ of total added feed |
| PhaseDensity | a declared specification, phase $p$ | the driver's phase molar density | mol m$^{-3}$ or kg m$^{-3}$ |
| EquilibriumPressure | $T$ fixed, one free and one pinned-zero amount (bubble, dew, pure saturation) | the driver's free pressure | Pa |
| EquilibriumTemperature | $P$ fixed, the same phase specification with $T$ free; a cold guess must lie above the answer (a guess below reports a coalesced-phase topology event) | the driver's free temperature | K |
| PhaseMoleFraction | any declared specification, phase $p$, nonempty species set $S$ within its support | $\sum_{i\in S} x_{p,i}$ | 1 |
| PhaseFraction | free amounts | $\beta_p = n_p / \sum_q n_q$ over the solved amounts, never the feed total | 1 |
| PhaseAmountRatio | declared numerator phase/species set $(p,S)$ and denominator $(q,S')$ | $A/B$ with $A=\sum_S n_{p,i}$, $B=\sum_{S'}n_{q,i}$ and $B\ne0$ | 1 |
| PartialPressure | any declared specification, phase $p$, nonempty species set $S$ within its support | $P\sum_{i\in S} x_{p,i}$ with the solved system pressure; `LogRatio` gives $\ln(p_S/p_{S,obs})$ | Pa |

`AdditionEnthalpy` evaluates two separately declared endpoints; their
specifications and temperatures may differ. Both use the same parameters,
thermochemistry, ordered reactions and source reference. The compiled feed
amounts must satisfy $\mathbf n^{(2)}_{feed}-\mathbf n^{(1)}_{feed}
=n_{add}\mathbf y$ componentwise within $10^{-12}\max(1,|n|)$. Each feed species
needs an ideal-gas record, and the feed shares the endpoints' parameters and
record. The endpoint value sums $n_p h_p$ over every phase not pinned at zero,
including a phase with nonzero vapor amount. When both endpoint problems are at
the same fixed pressure and the feed enters at that pressure, the result is the
isobaric heat of addition; for free-pressure endpoints it is an enthalpy
increment, whose calorimetric interpretation depends on the caller's energy
balance. Referenced free-temperature endpoints consume the same native A1
calorics at each solved temperature, with fixed feed-state actions.

An unavailable prediction is typed, never a zero residual: a vanished or
coalesced phase is `topology_event` before any solve status; a rejected or
state-less solve is `central_solve_unavailable` with its convergence record;
a missing density root is `density_unavailable` with the thrown cause; a mass
basis without a recorded ionic molar mass is `missing_molar_mass`. Missing,
outside-interval or nonfinite thermochemistry for total $c_p$ is
`thermochemistry_unavailable`, with the Ideal status name as its cause. Every
declared component requires its own thermal interval, including zero fractions.

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

`AdditionEnthalpy` uses a signed value, which is negative for exothermic
addition. Exothermic data use `Difference`; `LogRatio` requires positive
predicted and observed values. Any other constant reporting basis is folded into the
observed value and, for `Difference`, its scale; a state-dependent
per-absorbed-amount basis is unavailable.

A direct source-basis ln K row uses `ObservationKind.LogEquilibriumConstant`,
`Observation.reaction` for the reaction name and the fixed temperature in its
problem. Its problem supplies the same reaction/source declarations as the
other observations; phase, pressure and feed data are not consumed for a direct
prediction. EqID: `source_log_equilibrium_constant_observation` gives its value,
exact coefficient columns and zero EOS columns. The same #187 naming, tying and
fixed-target rejections apply. Use Difference: ln K is already logarithmic.
For absolute standard uncertainty σ in ln K, choose scale=1 and weight=1/σ².

The `AdditionEnthalpy` prediction derivative uses the endpoint actions for
both extensive enthalpies and phase amounts, then subtracts the feed residual
enthalpy action at fixed feed state. It consumes `PhaseAmount` and
`TotalEnthalpy` central values from the solved-state owner rather than
reassembling phase calorics in parameter regression. The feed's ideal enthalpy is fixed
thermochemistry; its EOS-coordinate residual action comes from the Properties
pressure chart. Build the direction once per endpoint using that endpoint's
solved temperature and reaction-row mapping, and keep both endpoint solves
through Jacobian evaluation. Reaction-coefficient enthalpy actions include the
temperature-slope derivative $\phi'_k(T)$ (EqID:
`reaction_constant_coordinate_derivative`); EOS-coordinate actions consume
$C_{T\theta}$ (EqID:
`standard_state_transfer_temperature_parameter_derivative`). A refused
column returns `action_unavailable` for that named endpoint and parameter or
reaction, while other directions and central values remain available.

Admission rejects nonfinite or nonpositive $n_{add}$, an invalid feed state,
the first component that violates the material balance, or endpoints with
different ordered reactions or source references. Endpoint solve failures keep
their `initial:` or `final:` prefix and driver cause. A missing feed density
root returns `density_unavailable` with its convergence record; a missing,
out-of-interval or nonfinite feed record returns `thermochemistry_unavailable`
with the `feed:` prefix. Refused endpoint caloric values retain the endpoint
prefix and `thermochemistry_unavailable`; unavailable feed actions retain
`feed:` and `action_unavailable`. No failure becomes a zero residual.

The existing constructor can map a declared `Problem`; then set the direct kind
and named reaction before passing the row to `fit`:

```python
row = regression.observation(parameters, "pressure", problem,
    observed=ln_k_observed, form="difference", scale=1.0, mass_basis=False)
row.kind = regression.ObservationKind.LogEquilibriumConstant
row.reaction = "named_reaction"
```

## Parameter coordinates (R1)

An EOS physical coordinate is exactly the Properties active-parameter coordinate:
one component slot for $m$, $\sigma$ and $\varepsilon/k$, one reciprocal pair
for $k_{ij}$ and $l_{ij}$, one site pair or diagonal self-pair for
$\varepsilon^{AB}/k$ and $\kappa^{AB}$, one reciprocal pair with a positive
reference temperature for the $k_{ij}$ slope $b_{ij}$ of EqID:
`cross_dispersion_energy`, in native units (1, Å, K). Every slot
absent from the active list is fixed. EqID: `affine_parameter_coordinate_map`
owns the optimizer mapping. A `CoordinateTarget` declares a family, indices
and a finite nonzero coefficient. Global dimensionless Born scalars own index
`[0]`; their Python targets use `(family, None, coefficient)` with no component
name. Other targets own component or reciprocal-pair indices. One coordinate
may drive targets across families; its positive scale and inclusive bounds
apply to the base value q. A reciprocal pair is one physical row. Independent
EOS rows follow coordinate and target declaration order, followed by held Born
rows and reaction rows. The active list contains the EOS and held rows only;
`tie()` gives row-to-column identity, not derivative magnitudes.
The public point reader's `optimizer_gradient` is dy/dv;
`residual_optimizer_gradient` is dr/dv at unit row weight. These coincide for
Difference scale 1; otherwise Difference divides by its declared scale and
natural-log LogRatio divides by the native prediction. Ceres applies its row
weights to residual columns; minimax applies native check slopes to prediction columns.

Native `CoordinateSet.start()` divides the first target value by its coefficient
and inverts the base origin/scale. Every stored target and mirrored entry must
agree with its forward-mapped value within eight machine epsilons times the
larger magnitude; coefficient-one ties require exact stored equality. Zero
has zero allowance. Nonfinite values, inconsistent equalities and starts outside
bounds are rejected without averaging or projection. Bounds are enforced in
optimizer space. The existing boundary-rounding clamp applies only to q before
multiplication; exact columns describe admitted points, including feasible
one-sided boundary actions. Targets have no independent offsets or clipping.

A reaction coordinate is a named scalar ln K (coefficient a) or a stored
`ReactionLogPolynomial` coefficient a, b, c, d, in units 1, K, 1, K⁻¹. It uses
the same affine map and inclusive bounds. Name tying requires identical
stoichiometry, constant representation and starts, correlation coefficients,
Tref, interval and standard-state identity, and source basis/pressure. Referenced
matches also share the solvent, the common-molality standard and the effective
species reference of every reacting species, since no conversion between
conventions is declared (EqID: `neutral_reaction_reference_transfer`); their
dilution-path ratios may differ. The derivative is EqID:
`reaction_constant_coordinate_derivative` evaluated at the solved temperature.
Source/reference declarations are preserved in every trial. Scalar constants
are rewritten as exp(ln K); correlations keep their native stored coefficients.
Affinity targets, supplied Gibbs values, dependent reaction sets, thermochemistry-
fixed unreferenced constants cannot be fitted. Coefficient motions are evaluated
by the native law at T(t), including the solved-temperature chain and thermal
products. Apparent-Cp fitting consumes Equilibrium H A2(T,coordinate), including all four
phase-product terms. Physical increments are divided out once before affine
ties, residual slopes and weights. No outer coefficient Hessian is consumed.

The following declaration uses the same native chain for the fixed-pivot
direction in `reaction_constant_coordinate_derivative`:

```python
regression.reaction_enthalpy_shift(
    name, pivot_temperature=..., origin=..., scale=..., bounds=...
)
```

Pivot temperature is required in K;
origin, scale, bounds and the physical result are in J/mol. Ordinary least squares
and coordinate-aware `regression.differentiate` admit the declaration; sequence
starts and minimax refuse it. The point reader's residual columns have unit
weights. `result.chart.reaction_laws[name]` gives the fitted native correlation;
`regression.reaction_log10_reciprocal(law)` supplies $(A,B,C)$ for the declared
negative base-10 reciprocal form. Replay preserves caller-owned source metadata.

A plain coefficient and one enthalpy shift may target the same reaction. Each
trial assigns the plain coefficient first, then adds the pivot direction to the
fresh supplied law. The pivot invariant compares laws at the same plain value.
The plain physical value remains distinct from the completed law's coefficient.
Two shifts on one reaction and duplicate plain coefficients are refused.
Replay attaches `result.chart.reaction_laws` and fitted EOS parameters directly,
without applying the coordinates a second time; preserve native double values.

For a nonzero starting shift in a shift-only fit, pre-apply the known amount
$h_0$ to the supplied native law, start the fitted increment $h$ at zero, and
translate total-shift bounds $[l,u]$ to $[l-h_0,u-h_0]$. Serialize $h_0$, $h$ and
the total $h_0+h$ alongside the completed law. The Engine cannot infer an
unreported pre-applied amount. This recipe does not preserve an unshifted plain
coefficient meaning in a combined pair with nonzero $h_0$: that would also
require translating the plain coefficient's origin and bounds. The required
combined chart uses $h_0=0$.

## Prediction Jacobians and fitting (R2, R3)

$\partial y/\partial\theta$ comes from the shared owners only: the Properties
pressure-chart density action for Density, the pressure-chart total caloric
action for IsobaricHeatCapacity (one direction per call), and the solved-state
actions of the equilibrium driver for the equilibrium kinds (a linear sum over $S$, the quotient
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
box is shortened along its direction to the first bound. The trust-region ratio
keeps the full step's predicted decrease while the evaluated step may be shortened
or clamped, and is not shown to be conservative: in the 149-row reacting MEA fit
(#455) it rejected five shortened steps that lowered the cost (full-step ratio at
most $5.6\times10^{-4}$, shortened-step ratio above 0.9997, threshold $10^{-3}$)
under the earlier release rule, and no such step under the current one. Bounds form an
active set in optimizer space, against the physical bounds mapped by the affine
transform (rounding is monotone, so a physical value inside its bounds is inside
these): a coordinate within the named `FIT_ACTIVE_BOUND_DISTANCE` of a bound
whose objective gradient $J^\top r$ points out of the box is on it and is held
where it is for a round, never moved onto the bound; a free coordinate on a
bound it moves toward is clamped instead of shortening the step. A round ends
when Ceres converges, an accepted step puts a free coordinate exactly on a
bound it was not on at the previous accepted point, or the full accepted
weighted optimizer gradient $J^\top r$ selects a coordinate that is not held.
Rejected trials cannot change the set. The accepted-point Jacobian is reused
without another thermodynamic solve; a failed gradient stops with its typed cause.
Any other Ceres stop ends the fit. The set is rebuilt at that point: selected
coordinates join it, and the held coordinates no longer selected (gradient turned
inward) are released together, but only after a Ceres convergence report or when
no coordinate would otherwise stay free. On a held coordinate $\pm g_c$ is its
bound multiplier only where the free gradient vanishes, so release waits until
Ceres stops on the current set (the drop timing of Nocedal and Wright, 2006,
Algorithm 16.3). The convergence report (function, gradient, parameter or
minimum-radius test) is a stopping trigger, not a stationarity certificate;
releasing several coordinates at once is a choice, revisiting a face is not
excluded, and a budget stop leaves before another set update, so a needed release
can be missed (the fit reports its limit). The fit converges when a converged
round leaves the set unchanged or every coordinate is held. Unchanged membership
retains Ceres' function, parameter and gradient stopping criteria; it does not
alone certify a small free projected gradient. Every coordinate held satisfies the bound-gradient signs
to the declared distance. A bound within the same
distance is active for the covariance, and physical values are clamped to their
bounds against the affine map's rounding. Iterations
and time are shared by the rounds; a
trial in which any positive-weight objective observation or action is unavailable is rejected with its
typed cause and observation index retained, never zeroed. Exact raw equilibrium
declarations, including feed basis, phase supports and seeds, reaction/reference inputs and row tolerance, share one compiled central
solve within a fitting evaluation; projection, residual and weight metadata
remain independent. Density and heat-capacity observations remain separate pressure closures.
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
only from a trial in which every positive-weight objective observation succeeded and are retried cold on
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
energy order. The total $c_p$ first parameter Jacobian uses third energy
derivatives, the pressure density lift and implicit site coefficients through
the existing Properties action (EqID: `fixed_pressure_heat_capacity_parameter_derivative`).
The ideal-gas record is fixed and is not fitted; errors in it can be absorbed
into fitted residual parameters. No parameter Hessians or reacting-equilibrium
heat capacities are claimed.

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


## Declared measurement uncertainty

`MeasurementUncertainty` selects errors-in-variables estimation within the
existing native Ceres fit. It declares input and output standard deviations in
physical units and, optionally, one known joint covariance or correlation; no
variance is fitted. The objective and exact
columns are EqIDs `errors_in_variables_objective` and
`errors_in_variables_residual_derivative`, specializing Zavala (2008), thesis
chapter 4, §4.2, Eq. (4.2), printed p. 53. This defines a conditional synthetic
measurement model and establishes no physical validation or global optimum.

Each `MeasurementGroup` explicitly names output rows sharing one physical input
measurement. Distinct groups remain independent at equal nominal states. A row
belongs to at most one group; rows in a group must agree on state, branch and
reference declarations. One positive-sigma input creates one correction and
one prior residual, regardless of output count. A group containing only
zero-weight rows creates neither. Zero input sigma fixes that input exactly.
An absent or otherwise empty declaration follows the ordinary fit. Supplying
only output sigmas gives the equivalent fixed-input Difference weights
$w_i=(a_i/\sigma_{y,i})^2$, retaining absolute-sigma covariance meaning.
Positive external weights must equal one; arbitrary relative weights and
LogRatio do not express the selected physical-unit Gaussian likelihood.

Measured T is in K, P in Pa, independent fractions in mol/mol and salt molality
in mol/kg solvent. Bounds are explicit physical domains, independent of sigma.
The affine fraction chart, EqID `independent_feed_measurement_coordinates`,
resolves the original feed basis once and holds its declared total in mol.
Dependent species solve normalization and charge; the full box is checked
before fitting. Structural zeros, references, stoichiometry, species identities
and phase support stay fixed. No trial input is clipped or renormalized.
Salt molality uses the existing formula-composition owner and a fixed pure
solvent. Changing m does not change the reference ladder.

| Observations | Admitted uncertain inputs |
|---|---|
| Density, IsobaricHeatCapacity | T, P, independent fractions |
| PhaseDensity, EquilibriumPressure, EquilibriumTemperature, PhaseMoleFraction, PhaseFraction, PartialPressure, PhaseAmountRatio, EquilibriumIsobaricHeatCapacity | Fixed T/P and independent feed fractions of one ordinary equilibrium problem |
| LogEquilibriumConstant | Fixed source-law T |
| MeanIonicActivityCoefficient, OsmoticCoefficient | T, P, salt molality in one pure solvent |

Uncertain free solved T/P, compound addition enthalpy, solvation/transfer,
permittivity, caloric enthalpy rows, critical/envelope declarations, references,
hydration targets and mixed-solvent ionic inputs are unsupported in this mode.
Fixed rows can accompany admitted uncertain rows with declared output sigma.
Missing/nonpositive objective sigma, negative/nonfinite input sigma, unsupported
units, overlapping groups, conflicting states and nonphysical feed boxes are
refused with their cause.

`covariance` (physical units) or `correlation` declares one dense row-major
matrix; declaring both is refused. Its m rows are the positive-weight outputs
in observation order, then the positive-sigma inputs in
`column_groups`/`column_inputs` order, so correlations may cross measurement
groups. Zero-weight outputs and zero-sigma inputs occupy no row; all sigmas
remain required. With D the diagonal of those sigmas, R = D⁻¹CD⁻¹. A wrong size,
an empty or nonfinite matrix, a diagonal or symmetry error above 64 m ε, or a
matrix that is not positive definite is refused before fitting. The objective
is ½‖L⁻¹e‖² with R = LLᵀ (EqID `errors_in_variables_objective`).

The existing point-local parameter map D, residual slope and weighting act once
on parameter and input columns. Equilibrium actions share one defining-Jacobian
factorization. Pressure density actions include the molar-mass composition
derivative for kg/m³. Ionic T/P actions differentiate finite, terminal and pure
states; m differentiates only the finite state and molality conversion. Raw
reference actions are the terminal-pair limit, admitted separately for each
consuming sigma when $\sigma_q$ times the limit change $e_f$ of EqID:
`reference_limit_terminal_pair` is at most $5\times10^{-5}$. Failed actions retain their cause.

`InputReconciliation` returns computed inputs, standardized corrections,
penalties, separate output/input costs and the complete standardized Jacobian.
With a matrix, `corrections` remain the standardized physical corrections, while
`penalties` are the squared whitened input components, `output_cost` is half the
squared whitened output block (the marginal output term) and `input_cost` the
remaining conditional term, which depends on output errors under output/input
correlation; their sum is the objective. The Jacobian, `weighted_residuals`,
fit `optimizer_jacobian`, projected gradient, SVD and covariance are whitened;
the public point `optimizer_gradient` stays the raw dy/d(v,u) and
`residual_optimizer_gradient` holds the complete whitened columns. Group labels
identify measurements, not owners of whitened rows.
Its first rows are output rows; subsequent rows and columns follow
`column_groups`/`column_inputs`, after the model-parameter columns. `row_groups`
uses the observation count for an ungrouped output. The authoritative request
retains measured values, sigmas, units, bounds and membership. Ordinary model
parameter vectors and output-row arrays retain their existing order and width.
Objective predictions use reconciled inputs; stored-parameter checks use the
original measured inputs. Report-only failures do not change objective acceptance.

At a regular converged interior point the model covariance is the vv block of
$(J^TJ)^{-1}$ for the full augmented J, then mapped by D. This eliminates input
directions; using only output parameter columns would condition on fitted
inputs. Known absolute sigmas or covariance have no common variance-factor
estimate; the covariance assumption text names the joint Gaussian covariance
when a matrix is declared. Rank
deficiency, unfinished fitting or any active model/input bound withholds
covariance. This is a local Gaussian approximation. Ordered starts reset all
input corrections to zero and keep the shared work budget; minimax with an
uncertainty declaration is refused.

## Electrolyte observations and dependent Born coordinates

[ALG-FIT](../../science/algorithms.md#alg-fit) owns native execution and exact
derivatives. Gamma and Phi are dimensionless on the declared formula-unit
molality basis; solvation and transfer are signed J/mol in the Figiel x-process
convention. Transfer is reference A to composition B. Permittivity is dimensionless;
a declared reference composition makes its prediction a ratio, with the exact
quotient derivative of both pressure-closed states. Liquid branch, solvent
mass, charge-neutral positive integer stoichiometry and selected-ion membership
are admitted before fitting. Phi keeps the single-solvent molality floor.

Each request-local reference shares T/P, solvent fractions, ions and stoichiometry.
Finite states also retain molality or full composition; gamma and Phi share the
same finite activity evaluation. An unchanged optimizer point retains the model,
references and finite states through value-to-column promotion. An active-parameter
change discards them. Values and admitted first columns stay with Properties;
parameter regression performs projection, residual transforms and the chart.

A hydration hold declares a solvation observation/target, Born bracket, energy
and diameter tolerances, and iteration limit. Newton uses exact H_d and bisection
when its step leaves the bracket. Acceptance requires both |H| and the exact
Newton correction or refined bracket width within the inner tolerances.
A held slot cannot be independently free or coupled to another held reference.
The evaluated chart returns physical rows (independent EOS, held Born, then
reaction rows), D, achieved errors, H_d and raw parameter columns. The fitter
and `evaluate_coordinates` use this same chart. Raw columns hold other physical
parameters fixed; the complete optimizer chain includes target coefficients, scales and holds.
Covariance maps as D C_u D^T and remains conditional on fixed hydration targets.

`CheckGroup` declares row membership, absolute native differences or absolute
log10 ratios in decades, max/mean limits, baseline floors and report-only use.
Baseline limits come from a separate fixed parameter declaration. Final checks
use the stored fitted diameters, preserve failed members, and expose signed
errors/slopes, magnitudes, argmax, thresholds and baseline metrics. They change
no least-squares weights or objective. Zero-weight rows are admitted but evaluated
only finally; a logarithm-domain failure preserves its independently valid value.
Optional ordered optimizer starts run through the native `fit_sequence` owner,
with per-fit controls and shared iteration/elapsed limits. Each point reapplies
the same chart; returned candidates retain their checks and causes. Points not
started after exhaustion are absent. This is a bounded local sequence, with no
global-optimum claim.

The native `observations_at` operation constructs manufactured recovery data
from the same evaluated chart and observation evaluator; signed input offsets
are native-unit differences or decades. It performs no fit. This is numerical
verification input generation, not a calibrated parameter dataset.

Singular vectors label independent optimizer and objective-row orders. SVD,
rank, bounds and projected gradient are local diagnostics. Gamma and Phi are
related by Gibbs-Duhem, so independent-error covariance is conditional sensitivity,
not measurement uncertainty or structural identifiability. Profiles and native
constrained minimax estimation remain outside this slice's delivered capability.

## Native constrained estimation

An explicit `MinimaxRequest` declares target rows, objective scale $q$, one
positive scale per check group, positive scaled tolerance, bound-activity
distance and starting margin. The existing fit controls declare iteration
and elapsed limits; the existing native ordered sequence owns shared budgets.
Python transports declarations and reads results. [ALG-FIT](../../science/algorithms.md#alg-fit)
owns the method and numerical lifecycle; the EOS/reference and eliminated-hold
equations remain unchanged.

For targets, $e_i=\log_{10}(y_i/y_{obs,i})$ in decades, with positive arguments.
For guards, $e_j$ is either that unweighted decade error or $y_j-y_{obs,j}$ in
the observation's native unit. Frozen source-baseline rules and floors retain
the existing check meaning. The physical nonlinear program is

$$
\min t,\quad t\geq0,\quad \pm e_i-t\leq0,\quad
\pm e_j-B_g^{max}\leq0,
$$

with a separate $s_{gj}\geq0$ per mean-group member,

$$
\pm e_j-s_{gj}\leq0,\qquad
\frac{1}{n_g}\sum_j s_{gj}-B_g^{mean}\leq0.
$$

The optimizer objective is $t/q$, and each inequality is divided by its
positive declared scale $\sigma$. The target scale is $q$. Errors/slacks retain
decades, J/mol or dimensionless units in the returned point; optimizer
coordinates retain the existing affine convention. Exact first error columns
are $(dy/du)/(\ln(10)y)$ for decades and $dy/du$ for differences. Ceres weights
and residual scales do not enter these constraints. Report-only failures leave
the constrained feasible set unchanged; unavailable enforced members cannot
pass by omission. Final checks use $|e|$, not loose optimizer slacks.

For scaled $c\leq0$ and $l\leq v\leq h$, the diagnostic Lagrangian is
$L=t/q+\lambda^Tc+z_L^T(l-v)+z_U^T(v-h)$, with nonnegative multipliers.
`MinimaxResult` retains the positive-part feasibility norm, infinity norm of
$\nabla F+J_c^T\lambda-z_L+z_U$, negative-part multiplier norm and infinity norm
of complementary products. Auxiliary bound distances are divided by their
error scale; independent distances use $u$. Physical row multipliers are
$q\lambda/\sigma$, physical bound multipliers are $qz$; their units follow the
physical objective in decades per constraint/variable unit. Constraint labels
retain every target/max/link pair and mean aggregate, with original row/group
indices. Slack labels retain group and member indices.

Eliminated holds have no Ipopt equality multiplier. A regular interior hold
permits reconstructed $\nu=-L_d/H_d$ and held/full independent stationarity;
independent tied coordinates retain their affine scales in that comparison.
Bracket-boundary or singular holds retain their cause without a full-space
certificate. Missing final points, exact actions or multipliers make the
associated diagnostics unavailable. `FitResult` retains the common candidate,
failures, chart and checks. Its minimax optimizer Jacobian contains raw value
columns; weighted residuals and least-squares costs/covariance are unavailable.
`evaluate_minimax` returns point values/columns/checks only, with no solver or
KKT conclusion. Same-point root refinement tests numerical error, not physical
uncertainty. Manufactured minimax/equality examples and Engine synthetic rows
establish numerical verification only.
