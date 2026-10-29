# Scientific algorithms

Accepted requirements, reconciled against the 2026-09-06 transfer audit.
These methods consume the [equation definitions](../equations.md),
[EOS conventions](full-eos.md) and [equilibrium formulation](../equilibrium/formulation.md).
Algorithm IDs below are stable references. Current implementation and required
numerical falsifiers are mapped in the owning [GitHub issues](https://github.com/tannerpolley/ePC-SAFT/issues).
No section claims that all required actions currently exist. The implementation
plan owns current versus deferred scope. PH/PS uses the neutral caloric specifications below; neutral critical-equation
candidates and phase envelopes use ALG-EQUILIBRIUM; UV remains future work. The
explicit native fixed-T/P discovery and tangent-plane testing path in
ALG-STABILITY has retained numerical qualification on the stated neutral and
admitted ionic/reacting cases, supported by the independent delivered review
within the retained limits. Reuse accepted numerical owners and independent cases.
Declared-branch continuation and mechanical/local checks remain required.

See [ideal mole-fraction and molarity specifications](#ideal-mole-fraction-and-molarity-specifications)
for the ideal-liquid specialization.

## ALG-ADMISSION

Inputs are ordered components, explicit charges, source-backed coefficients,
units, fixed topology, source domains, state/specification and requested actions.

```text
resolve only selected components, pairs and sites
check required records, units, topology and mathematical domain
check physical charge/balance constraints for ordinary states
intersect known source applicability domains; retain unknown applicability
label requested extrapolation separately from physical validation
reject missing physical inputs; never substitute a coefficient or charge default
```

Unconstrained composition differentiation uses an explicitly declared extension
of the energy. It does not authorize charged physical states. Parameters fixed
by the request have zero derivatives; dependent combining-rule outputs carry
the chain from independent coefficients. Discrete topology is not fitted.

## ALG-CALORICS

Evaluate one explicitly selected source correlation on its declared domain.
Use its own integration constants and pressure standard; do not infer a common
anchor, source interval or continuity at a join. For piecewise inputs, require
adjacent intervals to meet at one explicitly owned endpoint, reject overlaps and
gaps, and return the owned branch value with cross-join derivatives unavailable
when the branch is discontinuous. Differentiate smooth correlations on the same
selected interval; reselect the source interval when a temperature leaves it. Missing
thermal inputs make only consumers of those inputs unavailable.

Ideal and Total calorics consume one component interval vector in installed EOS
component order and one explicit $P_0$. Ideal $h,u,s,c_p,c_v$ follow their
registered equations; Total values add the existing residual fields at the same
state. Missing thermal input leaves residual outputs intact and makes only Ideal
and Total consumers unavailable. An owned endpoint has a value, but a derivative
that would cross an undeclared join is unavailable. No total value falls back to
its residual part.

For the centered polynomial, evaluate $C_p=\sum_k a_k(T-T_{ref})^k$ and anchor
$h$ and $s$ at the supplied $T_{ref}$. The constant term uses the exact
$a_0\log(T/T_{ref})$ identity. No insertion derivative at a zero composition
amount is promised.

For each nonzero degree $k\geq1$, let $\Delta=T-T_{ref}$ and
$u=\Delta/T_{ref}$. Accumulate the entropy contribution as
$a_k|\Delta|^k$ times a dimensionless moment, multiplying the coefficient
before each power. This avoids the artificial overflow of $T_{ref}^k$ or
$|\Delta|^k$ before a small coefficient is applied. It does not cure an
ill-conditioned sum of distinct polynomial terms or a genuinely unrepresentable
result. The constant logarithm uses `log1p(u)` near the anchor and a difference
of positive logarithms away from it, avoiding overflow/underflow of $T/T_{ref}$.

For $0\leq u<1$, use
\[
\frac{I_k(u)}{u^k}=\frac{u}{1+u}\sum_{j\geq0}t_j,\qquad
 t_0=\frac1{k+1},\qquad
 \frac{t_{j+1}}{t_j}=\frac{u}{1+u}\frac{j+1}{k+j+2},\qquad
 I_k(u)=\int_0^u\frac{z^k}{1+z}\,dz.
\]
The ratio is at most $1/2$. With $b$ binary significand digits, retain the
initial term and $b+2$ updates; the relative geometric tail is below $2^{-b-2}$.
For $u\geq1$, use $J_0=\log(1+u)$ and
$J_p=1/p-J_{p-1}/u$ through $p=k$. Its error multiplier $1/u\leq1$
avoids the unstable near-anchor forward recurrence. Form $1/u=T_{ref}/\Delta$
and the initial logarithm from positive temperatures so an overflowing ratio
need not be formed in this branch.

For $u<0$, write $v=-u$ and $c_k=\exp(-1/k)$. Below $c_k$, use
\[
\frac{I_k(-v)}{v^k}=(-1)^{k+1}v
 \sum_{j=0}^{N-1}\frac{v^j}{k+j+1},\qquad
 N=(b+2+\lceil3\log_2(k+1)\rceil)(k+1).
\]
Evaluate this positive series by backward Horner summation. The omitted sum
is bounded by $v^N/[(k+N+1)(1-v)]$; since $v\leq\exp(-1/k)$,
$v^N\leq\exp(-N/k)$. Differentiating the positive-term geometric majorant
bounds differentiated tails as well; the degree-dependent margin suppresses
these through the third order tested here. Floating-point rounding remains
separate from this truncation estimate.

At and above $c_k$, use the finite logarithmic tail
\[
\frac{I_k(-v)}{v^k}=(-1)^{k+1}
 \frac{-\log(T/T_{ref})-\sum_{p=1}^k v^p/p}{v^k}.
\]
Here $v^k\geq e^{-1}$ and the numerator is at least $e^{-2}/2$, as seen
by retaining terms $p=k+1,\ldots,2k$ in the logarithmic remainder. Thus neither
the normalization nor a vanishing remainder amplifies roundoff as in a
fixed crossover independent of degree. Use the positive-temperature logarithm
even when subtraction rounds $v$ to one. Finite-precision summation error still
depends on degree and conditioning; the evidence qualifies its recorded cases.

Source coefficients, degree and interval ownership are immutable during a tape's
lifetime. Every temperature-dependent representation choice uses scalar
conditional expressions with finite clamped arguments in its inactive branches.
No recording-time value extraction or temperature-dependent loop count is used.
The same tape replays on its selected smooth source interval, including across
these arithmetic representation changes and the reference. Selecting another
physical source interval remains a separate admission operation.

## ALG-PARAMETERS

The solvent membership and normalization for `solvent-mole-weighting` are
fixed parameter entries, retained through selection, reordering, fitting-worker
transport and optimizer identity. `bulk_relative_permittivity` owns the sum or
quotient; it reuses the existing component laws and automatic differentiation.


`BornDiameter`, `DebyeHuckelDiameter` and `PackingDiameter` own positive
Angstrom values in charged component slots. `IonSuppressionCoefficient` owns
one dimensionless charged-slot value; a common C ties all declared charged
slots in the existing affine chart. `DielectricVolume` (m^3/mol) and
`DielectricTemperature` (K) own one neutral-component slot each and are admitted
only under the density-temperature permittivity (EqID:
`density_correlated_bulk_permittivity`); an L/mol coordinate takes target
coefficient 1e-3 (#441). Their values use `ParameterRecord<Scalar>`
with the other active families. Charges, masses, formulation choices and Born
correction/solvation factors remain fixed. No cross-family factor tie is added.

Resolve source values into documented native units before constructing active
coordinates. Native units are dimensionless for segment count, `k_ij`, `l_ij`,
and association volume; Angstrom for sigma; and Kelvin for dispersion and
association energy-over-k. An active scalar owns one native slot for a scalar
family, both reciprocal row-major off-diagonal slots for an unordered
unlike-component `k_ij` or `l_ij` pair, or the corresponding association
site-pair slots. A distinct association site pair owns both reciprocal slots;
a same-site association pair may own its one diagonal slot when active.
Component diagonal interaction slots remain fixed. The coordinate is linear in
every owned slot, with unit derivative; no implicit log or normalization is
applied.

For `p_native = s_p p_source`, apply
`dF/dp_source = s_p dF/dp_native` and
`d2F/(dp_source dq_source) = s_p s_q
d2F/(dp_native dq_native)`. Active pair derivatives sum the two raw
reciprocal-slot derivatives. The pure packing coefficient has no `k_ij` or
`l_ij` dependency under `packing_coefficient_parameter_derivative`; their
energy and cross-diameter dependencies remain on the common Helmholtz tape.

## ALG-ASSOCIATION

Optional `trace(diagnostics=True)` observes existing mass-action iterations as
`association_iteration`: residual infinity norm, accepted Newton damping (zero
for non-Newton updates), pivot ratio (NaN without Newton), and Newton selection.
Its reason retains the analytic start, Newton cause or `Newton`; terminal convergence remains
on the span alone. `association_derivative_failure` observes the existing LU
at derivative admission: site count, rank, largest/smallest absolute pivots and
pivot ratio, with the singular or ill-conditioned cause. It performs no second
factorization.

The shared model retains at most two most-recently-used thermodynamic recordings,
keyed by composition support, electrolyte predicates, selected thermal intervals
and residual/ideal output availability. Parameters and thermal source data remain
immutable within that model. Switching recordings never reuses a numerical state:
each caller refreshes the point, site fractions, Taylor coefficients and required
site factorization under the existing model lock. A miss records the shared EOS
before replacing the older entry. This bounded reuse avoids rebuilding on every
aqueous/organic phase switch; it does not change the association solve below.

Use EqID: `association_mass_action`, EqID: `association_helmholtz` and the
implicit mass-action identities in full EOS §5. Declare initialization and
iteration budget; use Newton's method on the admitted site domain, with
step control that keeps site fractions in (0,1]. Accept only a finite residual
meeting the declared criterion. Return nonconvergence, domain escape or a
singular defining Jacobian explicitly; the last iterate is not success.

For the first update, form $q_a=\sum_b K_{ab}$ at $X=1$, with
$K_{ab}=\widehat\rho x_{i(b)}\nu_b\Delta_{ab}$. Keep that start when
$\max q_a\le1/\texttt{ASSOCIATION_JACOBIAN_CONDITION}=10^{13}$.
Above this scale, let $r_a$ be the common computed $q_b$ of all partners with
$K_{ab}>0$, if they are exactly equal; otherwise use $r_a=q_a$ (also for an
empty row). With $c=1+(q_a-r_a)$ and $d=\operatorname{hypot}(c,2\sqrt{r_a})$,
the initial estimate is $2/(c+d)$ for $c\ge0$ and $(d-c)/(2r_a)$ otherwise.
This is the positive root of the two-population mass-action reduction;
it is an estimate for a heterogeneous graph. A finite estimate in $(0,1]$
replaces the first update and consumes one of the existing fifty iterations.
An inadmissible estimate is discarded and the original first Newton/substitution
step runs from $X=1$, without clamping or an extra budget. Later iterations,
residual criterion and implicit derivative admission are unchanged. Detailed
capture names the admitted update `analytic association start`. The estimate
is outside the thermodynamic recording and is never differentiated.

```text
solve site mass-action values outside the reusable thermodynamic recording
form partials of the defining equations at the converged value
solve the implicit linear equations for each requested derivative action
include those derivatives in the common energy differentiation
verify derivative backward error and report conditioning
```

An equivalent stationary elimination is allowed only on its regular domain
and must reproduce those identities. A weighted association Hessian may lose
rank at zero species amount although the unweighted mass-action Jacobian does
not. Do not add composition floors to hide this defect. AD dynamic parameters
supply updated values, not automatically derivatives with respect to themselves;
association values and independent physical coefficients require their actual
implicit/explicit chains.

Phase Hessian orders differentiate only V, the positive amounts and T when it is
free (`PhaseHessian`; `IsothermalPhaseHessian` omits T). Fixed controls, absent
species and active parameters stay in the physical array layout as uncomputed
padding, not physical zero derivatives; full `Hessian` requests keep every
coordinate. Site partials of the requested columns use the smaller of forward
coordinate sweeps and reverse site-output sweeps, never the full tape Jacobian;
one unweighted factorization and its multiplier solve supply the implicit
gradient and Hessian, and each Hessian entry contracts the row coordinate's
site tangent in O(sites). The association gradient floor is the maximum over
the requested columns only, so it bounds exactly the gradient entries computed.

Optimize each recorded Helmholtz tape once with
`optimize("no_cumulative_sum_op")` before evaluation. CppAD retains conditional
expressions and their selected branches while skipping unused branch work;
cumulative-sum operators are disabled. Atomic dependencies retain CppAD's
dependency rules; no user-defined atomic functions occur in this recording.
Optimization discards Taylor coefficients, so a new tape's first evaluation
recomputes order zero. See CppAD's [optimization documentation](https://cppad.readthedocs.io/stable-2025/optimize.html)
and [conditional expressions](https://cppad.readthedocs.io/stable-2025/CondExp.html#optimize).

Reuse order-zero Taylor coefficients only on the same tape at a byte-identical
full input, including site fractions and active parameters. A changed input
recomputes order zero; a replacement tape starts without cached coefficients.
Positive-order sweeps recompute the requested directions while preserving order
zero. This changes neither the implicit derivative equations nor their admission
criteria.

## ALG-DENSITY

Detailed capture emits `density_scan_input` (T/K, P/Pa, lower/upper density in
mol/m3, interval count), `density_scan_point` (density, reduced residual,
pressure slope, linked), and `density_scan_skip` (density, attempted node, zero
link) with the originating association cause. `density_scan_failure` retains a
terminating non-association cause. `anchor_failure` retains the offered density
and original exception; `anchor_cold_retry` also gives the reference rung.
These records observe the existing scan and retry decisions.

At each failed grid node, the shared cold/local scan keeps the attempted ordinal,
density in mol/m3 at seventeen significant digits, and the entire original cause.
It continues the declared grid only for a primal `ConvergenceFailure` whose
criterion is `ASSOCIATION_RESIDUAL`. Any other `std::exception` propagates
unchanged when no earlier interruption exists; after an interruption it is
appended and terminates the scan with the combined refusal. An interrupted
scan throws `density scan interrupted by association` before slope searches,
cell processing or branch selection. It returns neither density nor stable-root
count; capture-disabled calls retain every visited failed-node cause in the
exception. Unvisited nodes after a terminal exception are not claimed evaluated.
The same rule applies to the local anchor scan. Pressure actions and fitting
predictions retain the complete cause as density unavailable, and required
equilibrium initialization fails. Existing optional attempts can discard a
refused scan without certifying its roots. An uninterrupted finite scan keeps
its existing search meaning, not exhaustive root enumeration.

Inputs include T, composition, target pressure, admissible packing interval,
root-selection request, numerical policy and declared search budget.

State construction evaluates the full Helmholtz gradient and every row of the
T, V and active-parameter Hessian columns. Uncaptured construction omits unused
composition columns; active Core capture retains the complete Hessian payload.

Density closure evaluates normalized Taylor coefficients of the shared residual
Helmholtz equations in volume at fixed temperature and composition. Order two
supplies pressure and slope; order three adds pressure curvature. Association
uses the primal mass-action solve and one checked Jacobian factorization for
its implicit coefficients. The stationary association potential eliminates the
highest site coefficient. Temperature and composition stay double in this
chart; tape calls retain their differentiated scalar dependencies. Gradient and
Hessian diagnostic requests still use the tape. Per-thread buffers hold scratch
only, with no model or request cache.

Every unanchored request checks the liquid side first, then locates the
vapor root when requested or when a sibling is needed; other liquid requests
check vapor existence. An explicit vapor request without an actual sibling
stops liquid location once a strict pressure bracket has positive slope and
positive curvature at both endpoints, retaining the deciding endpoint's fold
separation check and the same locating proposals. Under the assumed shape,
this convex bracket contains one stable crossing. The vapor bracket must lie
strictly below it. Only the vapor root is refined and published when both
sides exist; the unrefined liquid bracket is private existence evidence.
If vapor is absent, liquid location resumes fully before publishing its sole
root. Selected-liquid, unspecified and actual-sibling requests retain full
liquid location. Anchored requests retain their independent local search.
The [#512 retained evaluation counts](../../../build/issue-512/pressure-evaluations.csv)
qualify reduced work; quiet-host timing remains open. Removing refinement
samples does not guarantee detection of every otherwise observable shape
contradiction or evaluation failure. Sampled curvature and slope order test consistency with
one interior minimum of the pressure slope. They do not prove this shape.
On a concave pressure interval, the pressure tangents bound pressure from above;
on a convex interval they bound it from below. Their envelope can exclude the
corresponding crossing only with the existing pressure separation. A negative
slope is a loop witness. A positive intersection of tangents to the slope cannot
exclude a negative dip. This unresolved exit, unresolved order comparisons,
degenerate or failed samples, and an exhausted structured-search budget hand the
request to the inherited scan. A resolved sampled contradiction instead abandons
the shape inference and inserts a finite local subdivision into that scan.
The original evaluation failure remains in the observation; the scan retains
its association-connectivity and failure rules.

For density-ordered samples $a,b$, let $r=RT$, $q=P_\rho$ and $c=P_{\rho\rho}$.
The single `DENSITY_STRUCTURE_ROUNDOFF` policy owns
$g=16\epsilon/(1-16\epsilon)$, using binary64 epsilon conservatively in place of
unit roundoff. At an existing order exit, compute

$$E_q=g(2r+|q|),\qquad E_c=g\left(|c|+\frac{4(r+|q|)}{\rho}\right).$$

A contradiction is resolved when $c_a>E_{c,a}$ and $c_b<-E_{c,b}$, or when both
curvatures have the same resolved sign and the existing density-spacing filter
applies: $q_b-q_a>E_{q,a}+E_{q,b}$ on the concave side, or
$q_a-q_b>E_{q,a}+E_{q,b}$ on the convex side. Equality, overlapping intervals,
unresolved curvature signs and arithmetic outside the finite normal range
remain recoverable at that same exit; structured inference does not continue.
The unchanged `DENSITY_SLOPE_SEPARATION` controls slope inspection, not resolution.

With $a_2=2\,\mathrm{energy.c[2]}$, $a_3=6\,\mathrm{energy.c[3]}$,
$u=a_2/\rho^2$, $v=a_3/\rho^4$ and $w=2a_2/\rho^3$, the existing conversion is
$q=r(1+u)$, $c=-r(v+w)$. Its cancellation-aware norms satisfy
$r(1+|u|)\le2r+|q|$ and
$r(|v|+|w|)\le|c|+2r|w|\le|c|+4(r+|q|)/\rho$.
The rounding-factor count reserves nine for conversion, five for replacing and
evaluating norms, and two for comparison. The arithmetic model and product
bound $\gamma_n=nu/(1-nu)$ follow
[Higham (1993), pp. 784–785, Eq. (1.2) and the bound preceding Eq. (2.3)](https://nhigham.com/wp-content/uploads/2023/10/high93s.pdf);
these norm inequalities and their application here are the #513 derivation.
The lead's 2026-10-15 build amendment permits one `Sample` boolean, set from the
existing conversion intermediates only when all are finite and normal. Only
order-exit classification reads it; false keeps the recoverable scan fallback.
No coefficients are reconstructed and no additional pressure evaluation occurs.
The radii concern local assembly conditional on the returned Taylor
coefficients; they do not enclose complete EOS or implicit association error.

After a resolved contradiction, retain its density-ordered sample pair and all
already evaluated samples. Intersect that interval with the admitted cold domain
or the anchor's original local window. An empty intersection still abandons the
shape inference and scans the unchanged window. Otherwise start with
`DENSITY_BRACKET_INTERVALS` equal log-density cells, including the exact interval
ends and existing samples, then bisect each cell geometrically. Old nodes are
reused exactly. No enumeration cache, subdivision or pressure evaluation is
added before a resolved contradiction.

Evaluate all local mesh nodes before processing cells. Reuse the same scan-node
ownership, fold bisection and safeguarded Newton refinement. Local fold bisection
finishes the sampled slope transition, retaining its final endpoints and stable
side, rather than stopping after a pressure bracket. Retain signed pressure
crossings, signed slope transitions and exact-node owners in density order;
negative-slope pressure crossings are evidence only. A negative-slope cell that
rises through the target remains unresolved. Exact positive-slope node roots and
strict upward crossings use the existing closure criterion. Do not use the
single-minimum quadratic exclusion or the exact-midpoint shortcut locally.

Require agreement on two successive nested partitions, first possible at the
40/80-cell levels before inserted existing samples. Each fine crossing or fold
must have the same coarse parent cell or exact boundary owner, with the same
ordered signs and multiplicity. Counts or nearby root densities alone do not
establish identity. If the partitions disagree, subdivide again. A resolved
zero-crossing interval may agree. Keep the ordinary cold/local mesh outside the
interval, insert the exact interval ends, and process each complementary cell
once. Boundary-node roots have one owner; approximate density merging is not
used. Keep a negative-slope witness from the combined finite search.

Enumeration has a chosen aggregate cap of
`DENSITY_BRACKET_INTERVALS * DENSITY_MAX_ITERATIONS` = 2,400 **new pressure
evaluations**, checked before every uncached evaluation in the shared evaluator,
including fold and Newton helpers and unsuccessful attempts. A finer mesh is not
started unless all missing nodes can fit in the remainder. Existing per-bracket
limits also apply. This resource cap is composed from existing policy values;
it is neither a convergence guarantee nor a bound on pre-contradiction work,
complementary scanning, association iterations or wall time.

Unresolved helper work, cap exhaustion, adjacency before agreement or pressure
closure, nonfinite evaluations, association interruptions, tangent-pressure
ambiguity, or nonunique anchor admission raises `DensityStructureError`. It
retains the original T/P/sample operands and rounding radii, intersected interval,
last attempted cell count, resolved stable count, enumeration evaluations and
remaining cause. A successful local partition never erases an unresolved
complement. There is no cold retry or conditional single-loop polishing after
contradicted enumeration. Ordinary rows retain the sole-root polishing path.

Cold Vapor and Liquid requests select the lowest and highest resolved stable
roots; multiple stable roots make Unspecified ambiguous. A sole root retains
count one and Unspecified. A requested actual sibling contains the opposite
resolved root and its own convergence record. An anchor requires exactly one
stable root in its whole original window, keeps its requested branch and has no
global root count. Required species/reference admission retains its existing
witness decisions, stated as admission by a finite search rather than proof of
branch absence.

This is finite-search numerical resolution. Two nested meshes can both miss an
arbitrarily narrow pair of folds or crossings; exact point derivatives do not
bound derivatives between samples. Root counts describe the resolved stable
roots within the admitted bounds and inherited complementary scan. They do not
certify all EOS roots, Gibbs stability, composition stability, or physical phase
validity for cold decane or any additional phase.

Required shared density closure, the first neutral-reference anchor, a
parameter-observation column and a required boundary-phase vapor anchor retain
the complete typed refusal before any cold retry. Pressure-action batches carry
one native boolean derived from the exception's C++ type beside the complete
cause. Optional initializers and the Newton attempt retain their existing
fallback and mandatory final admission.

Under the same conditional shape, an anchored root's window edge facing the
slope minimum certifies positive slope throughout the window when that edge
has the root's curvature sign and positive slope. The root, both final refinement
endpoints and their bracket stay inside the local window. Every unproved case
uses the inherited local scan. This is local uniqueness only; caller branch
admission and corrected-density checks still apply. The existence count and an
actually located sibling are separate: an existence-only result never invents
another root density. Sole-root metadata and the scan's negative-slope witness
retain their existing meanings. Finite sampling cannot certify arbitrary
multi-loop completeness; resolved contradictions use the local subdivision above.

A neutral reference trace at system pressure first tries the compile-time
pure-solvent liquid density as an anchor for its pure-solvent point: the
anchored solve must find a locally unique stable root in its window. An anchored
failure other than `DensityStructureError` falls back to the cold scan and its
branch admission; later trace points
anchor on the previous point as before. Like those points, an anchored success is
admitted by local uniqueness and skips the global scan and its failure checks.

Identical unanchored T/P/composition requests for vapor and liquid may share one
bounded search and its stable-root list: the sibling is the opposite extreme
stable root, or the same unbranched root when only one exists. Anchored requests
cannot supply a sibling. Boundary initialization shares its two feed requests;
changed phase compositions still resolve their own roots with the same checks.

```text
scan deterministically in log density within the admitted packing domain,
    keeping the reduced pressure residual and analytic pressure slope
at a positive sampled local minimum of the slope, search by golden section:
    a negative-slope sample exposes a loop and joins the connected cells;
    a positive quadratic-model minimum estimate stops the search;
    exhausted work or a collapsed bracket without resolution is a failure
own an exact stable scan-node root once, keeping its adjacent cells connected
per cell, discard a negative-slope interval unless its pressure rises through
    the target; bisect such a cell toward the crossing until a stable sample
    splits it into two fold cells (an exact stable root there is its one root);
    exhausted work or adjacent densities without one is a failure
bisect any slope sign
    transition to adjacent representable densities or an exact slope zero,
    stopping early at a stable sample whose pressure passes the target
if no stable sample brackets the target and it is within the existing
    reduced-residual tolerance of that turning-point pressure, report an
    unresolved spinodal
trim the cell to its positive-slope side and refine each strict pressure
    crossing once with the existing safeguarded Newton method
select by declared branch policy and return the stable-root count
```

Outside a resolved contradicted interval, the inherited scan assumes one mechanical loop per isotherm, with one interior
minimum of `P_rho`. A second high-density loop can leave the stable liquid
interval inside one cell whose end slopes are both negative, as for cold
propylene (#391). When that cell's pressure rises through the target, the mean
value theorem requires a positive slope somewhere inside it, not one stable
interval. Midpoint bisection toward the crossing, stopping at floating-point
adjacency, exposes a positive slope or fails the search rather than returning
an incomplete root list. `DENSITY_MAX_ITERATIONS` bounds these hidden-crossing
midpoints per suspect cell, not the whole density search; the fold cells it
splits off have their own bisection budgets. An exact stable root at a midpoint
is taken as the cell's one root, which assumes the cell holds one stable
interval. This does not qualify arbitrary multi-loop isotherms: a
negative-slope cell whose end pressures lie on one side of the target can still
hide a pair of crossings. In the single-loop case, once a negative-slope sample
exposes the loop, its two
turning points partition the connected domain into increasing, decreasing
and increasing pressure intervals. Each increasing interval contains at most
one pressure root. Exact stable scan nodes and strict crossings have separate
ownership; two approximate Newton results are never counted or merged to
decide root identity. Association failures break scan connectivity.

Fold bisection uses the existing `DENSITY_MAX_ITERATIONS` budget and stops at
floating-point adjacency, retaining the stable-side endpoint. Pressure increases
monotonically across the stable side, so a stable bisection sample whose
residual has passed the target by more than `DENSITY_ROOT_RESIDUAL` already
brackets that side's unique root; the fold itself is resolved only when no such
sample appears, which is exactly the case needed to show that no root exists.
A cell whose stable end is already past the target by more than
`DENSITY_ROOT_RESIDUAL` R T rho_upper in pressure (vapor side above it, liquid
side below it) is rootless on that side
without any fold bisection, since the stable pressure only moves further away
toward the fold. It does not
estimate a fold from pressure-residual coincidence. Without such a stable sample, a
target indistinguishable from a fold pressure at `DENSITY_ROOT_RESIDUAL` is
explicitly unresolved,
including a true tangent root; this is numerical unavailability, not physical
nonexistence. The root residual tolerance is unchanged. Refinement retains both bracket endpoint
pressure evaluations and chooses the endpoint with smaller `|P_calc - P_target|`
(the upper endpoint on a tie) as the safeguarded Newton origin and publication
point. Its own reduced residual, exact derivative, association floor and
amplification accompany the result. For an approximate root, both
`|r| <= tolerance` and `(rho_upper - rho_lower) |dr/d rho| <= tolerance`
apply at that endpoint, with a finite, nonzero derivative. Here
`r = (P_calc - P_target)/(rho R T)` and
`Z_target = P_target/(rho R T)`, so relative pressure error is `|r|/Z_target`;
the reduced-residual tolerance alone does not bound that relative error.
Once the residual passes, the next sample is a quarter radius past the Newton
root toward the far bracket end; midpoint bisection remains the fallback when
that sample leaves the strict bracket. Selecting the better endpoint as origin
can change the Newton trajectory and evaluation count. Thus the represented
crossing lies within the returned local-sensitivity radius
`tolerance |d rho/d r|`; an exact numerical zero is a separate stop. This
bracket witness is needed because a small residual alone does not make the
first-order radius a density-error bound near a fold. It does not bound
floating-point error in EOS evaluation itself.

The quadratic slope estimate is evaluated only on a sampled local-minimum
bracket. It is a bound for the interpolating quadratic, not for the EOS:
a narrower unsampled dip with greater curvature is not excluded. Near a
critical point the pressure is locally cubic in density, not in log density,
so that limit does not make the interpolation exact. Search completeness
therefore remains bounded by the declared scan and supported shape; no
global stability or complete root enumeration is claimed.

An unspecified request with multiple resolved stable roots is ambiguous.
A genuine sole root remains unlabelled. Its side is read from the stored
negative-slope sample of the same scan (single-loop case, #423): below the sample
it is vapor-side, above it liquid-side, and with no sample it stays unlabelled; no
evaluation, density cutoff or tolerance is added. Liquid references (the
neutral trace, the ion-free ionic solvent reference, its pure solvent and the finite
electrolyte-activity state) refuse a consumed root on the established vapor side;
pure-species overrides retain their own branch admission, and unused reference states do not cause refusal. Fit observations keep the unlabelled root and its count. An anchored request uses the same
partition and counting method within its existing local window and keeps the
requested branch only after local uniqueness is established. It returns the sole
stable root already refined by that local scan, with its own convergence record,
combined anchor/refinement iterations and all density evaluations; it does not
publish the seed merely because its reduced residual passes. That window
does not establish uniqueness outside itself or global branch persistence;
a failed anchor must not silently switch to an unrelated cold root.

Fixed-composition density diagnostics use caller-declared, strictly increasing
positive nodes inside the packing admission cap. They reuse the pressure and
exact pressure-slope evaluator and the bounded adaptive slope search above.
The declared strict pressure-crossing brackets are recorded before adaptation;
extra adaptive samples and resolved crossings never change that count.
Each unavailable node carries its cause and breaks connectivity. Exact node
zeros have one owner; roots are not merged by approximate density distance.

The diagnostic pressure roots include negative-slope crossings. A spinodal
solves $P_\rho=0$ on the connected declared isotherm. Sign brackets are bisected
with caller-declared density/log-width limits and evaluation/time budgets;
exact numerical zeros are reported separately. Pressure residuals have units
Pa; spinodal residuals are $P_\rho/(RT)$, dimensionless. Budget exhaustion
retains resolved rows and the current unresolved bracket. A pressure within
the existing reduced-residual tolerance of a spinodal is reported as an
unresolved tangent pressure root. Negative-slope diagnostic roots do not enter
ordinary stable-branch selection or pressure inverse actions. The inverse
$d\rho/dP=1/P_\rho$ is singular at a spinodal. Canonical pressure, packing and
association equations and derivative orders remain unchanged.

An optional finite-domain search declares a known finite anchor, larger search
cap and `Gradient`, `PressureSlope` or `Hessian` evaluation request. Native
doubling and bisection retain the last finite/first unavailable densities and
cause separately from the geometric packing cap. It reports only the first
finite-to-unavailable transition reached from the anchor; later finite islands
are not excluded. Reaching the declared search cap without failure means only
that this interval remained finite. A finite-evaluation boundary is not a
spinodal or a thermodynamic singularity by definition. These finite scans do
not establish all roots or physical nonexistence.

## ALG-EQUILIBRIUM

All admitted specifications use the one equilibrium driver required by the
[project map](../../../docs/scientific/CONTEXT.md). Compile balances, supports,
pins and free intensives from data. An omitted phase support is every component,
restricted in a vapor to the zero-charge components. Explicit charged Vapor support
is admitted only under `DensityTemperature`, retaining its kind, ionic coordinates
and intensive electroneutrality even at zero pinned amount (#771's scoped amendment
to #438). `DipolarPerturbation` charged Vapor is refused pending separate numerical
qualification despite its vacuum limit; density-independent forms retain their
condensed-domain refusal. This is a numerical declaration capability, with physical
vapor applicability unvalidated. Supports are admitted
before any reference calculation; an omitted support keeps its full declared width for
supplied guesses and starts. A compile-time declaration error propagates to the caller
from every entry, including declared starts and continuation endpoints. The balance basis $B$ spans the left null space
of the reaction stoichiometry in echelon form (unit pivots, coefficients below
$\sqrt{\epsilon}$ set to zero), so each row counts one conserved moiety (#172). The pivots
are the species taking part in the most reactions, ties by declaration index, each accepted
while the remaining species keep the full reaction rank on the exact stoichiometry: shared species ($\mathrm{H^+}$, a carrier) become pivots and each reaction's
peripheral species stays off them, so a trace moiety such as $\mathrm{Li^+}+\mathrm{LiA}$
is one row. The rule is structural, so the basis is the same along a feed path. An
orthonormal SVD basis, or echelon pivots holding both $\mathrm{Li^+}$ and $\mathrm{LiA}$,
mixes that moiety into rows of major species, whose rounding then limits its balance and
makes the action Jacobian ill-conditioned in proportion to $N_{\mathrm{feed}}/u_{\mathrm{Li}}$
(equilibrated condition 9.6e10 and 1.1e10 against 345 and 276 at 1e-8 and 3e-10 mol/L Li).
Before positive charts or rank checks, cache the existing amount bounds on the
operative supplied basis, including stability trials. For ordinary reacting
declarations, an apparent LP zero proposes a full-component nonnegative conserved
functional. Verify its zero feed contraction and every original reaction
contraction with native error-free dot products. Only a verified certificate
removes a species; an unresolved apparent zero cannot prove absence. Reduce state
balances by independent echelon row operations while retaining their full feed
columns and the discarded inventory restrictions. Full reaction/reference
completion precedes restriction. For a pinned-zero support intersecting the
certified absent set, evaluate admission on its effective support after removal:
require a nonempty support of individually uncharged species, exactly two
opposite liquid/vapor phases, fixed positive T, free P, one zero pin and one
free amount. Apply the existing support and supplied-guess restriction, retaining
original declared widths for physical starts. Otherwise keep
`zero_inventory_incipient_support` with its species, phase and domain reason.
The formulation's emerging-phase limit on the same exact face owns this domain;
trial chemistry does not establish a global physical-inventory certificate.
Numerical qualification checks the reduced Jacobian at the solved endpoint,
separate mechanical branches and the unchanged raw/classified residual criteria.

### One-component two-density reduction

The general driver and native prepared `CompiledSystem.solve_saturation()` entry
reuse this single reduction/refinement owner. The prepared entry uses the same
eligibility predicate and falls back to the general entry for an explicit
options file or an ineligible declaration. Previous-result, optimizer-transfer
and bounded calls retain their existing general entry points. Rejection retains
an owned failed result and continues once from the original start, adding its
work once without replacing the fallback's own endpoint or linear diagnostics.
The bulk snapshot retains optional absence until a record is evaluated and an
owning reference to the complete result. Phase pressures are evaluated natively
at the actually exported densities. Direct pure-state evaluations retain inner
first observations in evaluation order, including rejected trials and wider
arithmetic; typed exceptions also retain their terminal record. Allocation
failure propagates.

An ordinary compiled declaration with one neutral nonideal component, fixed
positive T, free P, two liquid/vapor supports {0}, one free amount and one
zero-pinned amount uses the exact two-density reduction inside the existing
driver. Other Newton exclusions remain effective. The existing start owner,
full-row assembly, topology/local checks, branch admission and immutable-start
Ipopt fallback retain their responsibilities. Complete preparation remains the
owner of solved-state actions; the limited isothermal preparation supplies only
energy and amount/volume derivatives and their homogeneity reconstruction.

At unit phase amount, let $E=A^r/(RT)$, $v=1/\rho$,
$Z=1-vE_v$, $K=P_\rho/(RT)$, $u=\ln\rho$ and $r=\rho_V/\rho_L$.
The reduced pressure row is $(P_L-P_V)/(RT\rho_L)$ and the chemical row is the
difference of the full chemical rows after eliminating their common standard
term and multiplier. EqID `one_component_two_density` defines the rows and exact Jacobian. Their log-density Jacobian retains the derivative of the
liquid pressure scale. The full rows still determine the caller criterion.

The prospective double endpoint retains the compulsory fresh full correction.
Its corrected endpoint starts the polish residual history. Further exact
undamped double corrections continue within the same work limits until the
next log-density correction is at most `EQUILIBRIUM_POLISH_CORRECTION`, or
the reduced residual reaches its evaluation floor or stops decreasing.
The extended-precision engagement check follows this double polish.
The same equilibrated full-pivot factorization supplies correction and inverse
row sums by solving unit right-hand sides and undoing both scales. Define
$I=\|J^{-1}\|_\infty$, $\kappa=\|J\|_\infty I$ and
$E_u=I(\|F\|_\infty+\eta_q)$, in the original two rows and log densities.
The residual floor is
$$
\eta_q=\epsilon_q\max\left\{
1+|v_LE_{v,L}|+r(1+|v_VE_{v,V}|),\
|u_L-u_V|+|E_L|+|E_V|+|v_LE_{v,L}|+|v_VE_{v,V}|
\right\}
+\max\{v_Lg_L+rv_Vg_V,\ v_Lg_L+v_Vg_V\}.
$$
Here $g_p$ estimates the error in the unit-amount residual-energy volume
gradient, with units mol/m³. Its pressure-row and chemical-row factors differ.
The estimate is local and first order, with no rigorous global-error claim.

The double association path retains its admitted site solve, factorization and
existing gradient floor. The extended instantiation refines values and implicit
coefficients in affine volume and amount directions against the shared wide
mass-action residual. It reuses the admitted double full-pivot factorization as
a preconditioner and retains the achieved linear defect. Values and coefficients
stay wide. Each coefficient must meet working epsilon times its coefficient
scale within the cumulative association budget, including the starting solve.
Nonfinite or out-of-domain sites, rank loss, unchanged coefficients above that
floor, or exhausted work refuse the calculation with their cause.

For $q=(V,n_1,\ldots,n_c)$, use dimensionless affine paths
$q(t)=q_0+t s_i e_i$ and sums of two such directions, with $s_V=V_0$ and
$s_{n_i}=1$ mol. Let $w_a=n_{i(a)}\nu_a$, $r_k=[t^k]F(q(t),\widetilde X(t))$
be freshly evaluated achieved residuals after the last coefficient update,
and $J_k=[t^k]F_X(q(t),\widetilde X(t))$. The coupled remaining correction is
$$
J_0e_0=-r_0,\qquad J_0e_1=-r_1-J_1e_0,\qquad
J_0e_2=-r_2-J_1e_1-J_2e_0.
$$
The stationary association potential has derivative
$\mathcal A_{X_a}=-w_a r_a/\widetilde X_a$. For
$b_a=w_a r_a/\widetilde X_a$, the local gradient estimate is
$$
E_1=\sum_a\left(|b_{a,1}e_{a,0}|+|b_{a,0}e_{a,1}|\right),\qquad
b_{a,1}=\frac{w_{a,1}r_{a,0}+w_{a,0}r_{a,1}
-w_{a,0}r_{a,0}\widetilde X_{a,1}/\widetilde X_{a,0}}{\widetilde X_{a,0}}.
$$
Here $w_{a,1}=\nu_a s_{n_i}\mathbf1_{i(a)=i}$ on an amount path and zero on
a volume path; $w_{a,2}=0$. The gradient error is $E_1/|s_i|$. Second-coefficient
verification uses $E_2=\sum_a\sum_{k=0}^2|b_{a,k}e_{a,2-k}|$,
$2E_2/s_i^2$ for a diagonal Hessian entry and
$(E_2^{i+j}+E_2^i+E_2^j)/|s_i s_j|$ for a polarized mixed entry.
These are local linearized-defect estimates, excluding higher correction terms
and parameter/source uncertainty; they are not rigorous enclosures.

Properties retains each achieved direction-specific gradient error beside its
matching derivative. Chemical defining rows consume $e_{\mathcal A_{n_i}}$;
pressure rows consume $|V|e_{\mathcal A_V}$; balance and normalization rows have
zero association inner error. No separate first-order energy-value error is
added to these rows. The cheap double precision decision instead consumes the
existing cached double gradient floor, with the same pressure factor $|V|$,
without coefficient replay, another site solve, EOS evaluation or density work.
The shared row floor remains the componentwise
$\max\{\gamma_n|J|h,\eta_{\rm inner}+\epsilon S\}$, with its existing row scaling
and physical projection. Arithmetic error remains separate from achieved inner
error; a tiny inner estimate is not the total floating-point error.

Refinement engages once when the double estimate or its actual correction
exceeds `TWO_DENSITY_FORWARD_ERROR`, the 1e-10 log-density target. It repeats the shared Helmholtz,
temperature/geometry, polar and association arithmetic in long double, using
the existing Newton damping and bounds and remaining work limits. No fluid or
temperature selector is used. After a fresh full correction, the reduced rows
must meet the unchanged criterion and either both correction and estimate meet
EQUILIBRIUM_POLISH_CORRECTION, or the residual is at most the estimated floor
and the correction is at most $I\eta_q$. The estimated log-density intervals
must remain separated, $u_L-u_V>2E_u$. Existing branch admission still decides
whether the two phases are numerically resolved.

The exported double coordinates are prepared and checked again, with no
subsequent double Newton correction. Re-evaluation at the actually exported
double densities gives
$E_u=I(\|F_{\rm exported}\|_\infty+\eta_q)
+\|u_{\rm exported}-u_{\rm internal}\|_\infty$;
$\operatorname{expm1}(E_u)$ is the relative-density error estimate.
The result's linear_solve_diagnostics reports the row norm, floor, inverse norm,
condition estimate, coordinate forward-error estimate, working significand bits
and extended correction count. A rejected attempt is identified explicitly on
failure and is retained in trace evidence; its estimate is absent from a
successful fallback state's diagnostics. If long double has no more bits than
double, the admitted double result remains, with its actual estimate, 53 bits,
zero extended corrections and the reason additional precision was unavailable.
An estimate above 1e-10 does not claim that density accuracy.

For neutral, nonreacting mixtures, the selected ordinary target endpoint has a
separate precision decision (issue #864, amendment r2). Eligibility requires a
square two-phase liquid/vapor feasibility cell at fixed temperature and free
pressure, one free phase amount and one exactly zero-pinned phase, more than
one component, and no caloric, critical, envelope, azeotrope, discovery, source
reference, ionic, reaction or solid equations. Association is admitted through
the shared coefficient/error owner above. Ordered starts and external coordinate
correctors retain their methods. Ordinary previous-result starts are included;
permission to search for a recovery source is a separate decision. A terminal
recovery result is included because it supplies the ordinary target endpoint.
Intermediate Newton, Ipopt, source-search and continuation candidates cannot
activate this mixture correction.

At the selected point, use its retained prepared derivatives and matching rows
and Jacobian. The shared physical error radius is the maximum of pressure error
relative to $\max(|P|,1\,\mathrm{Pa})$, log-density and log-packing errors,
absolute composition errors, and phase-amount errors relative to feed amount.
An available finite radius above `MIXTURE_FORWARD_ERROR` ($10^{-7}$), with more
long-double significand bits than double, triggers the shared correction body
once. Singular or nonfinite estimates retain their cause and cannot trigger as
finite estimates. Without extra precision, retain the admitted double endpoint,
its own estimate and actual working bits with the unavailable-precision reason.
The pure-fluid `TWO_DENSITY_FORWARD_ERROR` target above remains $10^{-10}$.

Wide coordinates remain wide through all corrections. Use full-pivot rank and
finite-row/Jacobian/direction checks independently of the free-amount Newton
predicate; require every trial to stay in the finite admitted box. The ordinary
double damping fractions remain $(1,1/2,1/4)$; only wide correction also tries
$1/8$. Above the caller's row criterion, require strict raw-row contraction,
including the existing decreasing-ratio rule during extension. Below the
criterion, the bounded full-correction exception does not establish convergence:
fresh rows, the physical-radius target, and either a sufficiently small scaled
correction or the achieved-floor correction condition must all pass. Rank loss,
nonfinite data, unavailable estimates, no representable progress, no admitted
trial or spent iteration/evaluation/deadline limits retain their specific cause
and attempted work.

Export once to double. Re-evaluate the actually exported point in wide arithmetic
and include internal-to-export coordinate displacement and output rounding in
its estimate. Repeat the ordinary raw/classified-row, finite-bound, feed/zero-pin,
phase-order/topology, local-stability and density-branch admission. Replacement
requires all these checks and the returned physical-radius target; otherwise
retain the already admitted original endpoint with its own estimate and all
attempted work and first failures. Optional observation preserves a stored wide
returned-point estimate. A rejected attempt's estimate cannot describe the
retained original endpoint.

Untriggered states retain pre-refinement arithmetic, outputs, iteration count
and thermodynamic/density work. The fixed preservation inventory is independent
of the trigger. The forward-error notebook retains numerical verification on
fixed inputs; it establishes neither physical validation nor a complete-call
speed claim.

Repeated public solves and accepted continuation/envelope correctors may reuse the
previous Ipopt endpoint's primal coordinates and equality/bound multipliers when
the ordered model, phase, support, balance-basis, row and coordinate layout is
unchanged (#505). The tuple belongs to the equilibrium driver and is separate
from any Newton-polished returned state. Correctors retain their tangent-predicted
primal point. With $a=1/N_{feed}$ and amount chart $S_i$, equality multipliers
are stored as $a\lambda$ and bound multipliers as $a u_i z_i/S_i$, where
$u_i=N_{feed}$ for phase-amount coordinates and one otherwise; the next solve
reverses these transformations using its own feed total and chart.
Ipopt 3.11.9's documented warm-start initialization uses primal and multiplier
pushes of `EQUILIBRIUM_WARM_PUSH` and monotone barrier `EQUILIBRIUM_WARM_MU`
(the manual's *Warm-Starting Capabilities via AMPL*). Nonfinite operands, changed
layouts, fixed bounds, or equality/bound multipliers above Ipopt's effective
scaled cap refuse transfer. An explicit options file keeps the existing
primal-only route. Restoration or a failed warm attempt permits one primal-only
fallback within the same iteration/time budget; all attempted work is retained.
Density-branch admission runs once at its existing caller. Rejected corrector
steps disable transfer until a new state is accepted. Fit observations and fit
trials do not opt into multiplier transfer; there is no shared solver cache.

Check structural
counts/rank, construct declared feasible starts, then solve the same intensive
rows using the required common energy. Recompute achieved residuals and phase
distinctness at the returned point; optimizer termination alone is insufficient.
See [equilibrium initialization by declaration](#equilibrium-initialization-by-declaration)
for start construction, its precedence and declaration-specific exclusions.

The shared driver invokes [neutral recovery](#neutral-recovery) only within
its declared eligibility and budgets.

Caller-declared ordered starts share one thermodynamic problem and the same driver.
A start contains physical composition, packing, amount and free-intensive guesses,
a previous result, or a complete compiled-coordinate vector. Complete vectors
skip generated initialization. Physical composition vectors use the completed
declared supports: validate their original length and finite fractions in [0,1], then apply the
same species restriction as compilation before checking surviving coordinates.
Physical guesses that supply every composition, packing, free amount and free
intensive are used as given. Otherwise the start runs the plain solve's [cold start](#equilibrium-initialization-by-declaration)
and then replaces every supplied coordinate. A supplied composition is
not speciated over and skips the finite-phase sub-solve. The problem's own guesses
act exactly as in the plain solve. Work limits bound the work but never select the start, so
an empty start under limits that are not exhausted returns the plain solve's state and work (#350).
Invalid initialization is an individual attempt cause. The first available state
meeting the native row, topology and density-branch checks is selected; subsequent
starts remain unattempted. No best-residual or Gibbs comparison selects another root.
Aggregate attempts, optimizer iterations and steady-clock elapsed time bound the
ordered search. Remaining iterations and time enter the shared Ipopt options and
callbacks, of the seed sub-solves too; polish and final admission check the same deadline.
The driver does not run once the sub-solves exhaust the iterations, and the plain
solve's retry after a failed seeded attempt is not repeated. The default attempt
limit, `EQUILIBRIUM_MAX_ITERATIONS` (200), covers all of an attempt's stages together,
whereas the plain solve gives each driver run its own 200; a default start therefore
differs from the plain solve once its stages add up to more than 200 iterations. A bounded inner
kernel may finish after its deadline: elapsed time and overrun remain visible,
and a post-deadline state is unavailable. Elapsed exhaustion is a controller cause,
not a relabelled numerical stop: the common classifier retains its
nonfinite/converged/floor/budget priority and the native optimizer termination.
Measured rows that already meet the criterion can therefore remain numerically
converged with success false and no published state after deadline rejection.
This finite search establishes neither uniqueness nor global stability. Ordinary
solve initialization retains its rules. The retained [declared-start falsifiers
and default-solve cost](../../../analyses/2026-tpd-phase-discovery/results/issue238-cost/README.md)
qualify the bounded inputs stated there.

Public solves also apply the continuation owner's density-branch admission to
each returned phase. A cold pressure closure selects its declared liquid or vapor
root and refuses a sole root that ALG-DENSITY's stored negative-slope sample places
on the other side (#423), an anchored closure establishes local uniqueness, and the corrected density
must lie in that anchor's window with positive pressure slope. Failure clears the
state and retains the optimizer and row diagnostics. Continuation uses the same
predicate with its own anchors and retains its `branch_lost` cause; its internal
corrector does not independently reselect a cold branch. This is declared-branch
admission, not global phase discovery, and inherits ALG-DENSITY's finite-search
limits.
Retain the native optimizer `solver_status` alongside the convergence stop,
residual, and failed-attempt initialization/termination details. The optimizer
stops on amount-weighted rows, which bounds an intensive trace row only by
`tau N_feed/(N_p x_i)`; after a success, acceptable-point, tiny-step, restoration-failure,
iteration-budget or time-budget stop, `EQUILIBRIUM_POLISH_STEPS` Newton steps on
the unweighted rows with the exact Jacobian close that gap. Further steps continue
only while the raw row maximum exceeds the caller's row criterion and each step's raw
contraction ratio (after/before) is strictly below the previous accepted step's,
under a separate polish cap of `EQUILIBRIUM_MAX_ITERATIONS` steps (#261,
2026-10-05). In Newton's quadratic region the ratio falls toward zero; a constant or
rising ratio means the iterate lies outside it, so polishing stops. From a
box-limited start with raw rows near 10, the nominally unloaded MEA liquid (CO₂
4.9e-13, H₃O⁺ 1.6e-12) needs four steps, and three left every ionic row near
−1.8e-9 z_i² under a 7.6e-12 classified residual. Raw rows that meet the criterion
need no further polishing. A strictly raw-residual-decreasing proposal that reaches
the caller's raw-row criterion is admitted without the decreasing-ratio test (#670).
Best-point selection, freshly classified rows and the solved-state action criterion
remain independent checks. Unfinished explicit stability trials alone may continue
with the bounded adaptive damping in ALG-STABILITY; declared solves retain the
contraction-ratio policy above. Each step must reduce the raw row maximum at a
positive finite composition of every supported species; exact absence is the
compiled support (#231). No composition cutoff applies, so the polish reaches
positive roots below the optimizer box and below `EQUILIBRIUM_BOX_FLOOR`
(#87, 2026-09-21; #266, 2026-10-07: MEA at loading 1e-6, 318.15 K, CO₂ 1.1e-21).
A trace moiety row is taken relative to its own size in the polish's progress
measure, admission, stops, Newton system and #272 forward-error bound, and in the
classified residual (#633). Its size $s_k$ is the larger of the moiety's feed
inventory $\sum_i |B_{ki}| n^{\mathrm{feed}}_i$ and the magnitude
$\sum_p N_p \sum_i |B_{ki}| x_{pi} + \sum_s |\beta_{sk}| n_s$ of the terms the row
sums at the current point, both divided by $N_{\mathrm{feed}}$. When
$0 < s_k \le$ `EQUILIBRIUM_TRACE_MOIETY_FRACTION` $=10^{-6}$, the threshold of
#272's trace-species definition, the row and its Jacobian row are divided by
$s_k$. A feed-normalized trace row (9.86e-18 at 1e-13 mol) otherwise lies below
the rounding of the O(1) rows, so the raw maximum cannot fall and an available
correction is discarded. Row equilibration leaves the Newton direction unchanged
in exact arithmetic; the optimizer's constraints and the exact derivatives of a
solved state are unchanged. The term magnitude keeps a row that cancels large
solution-side amounts, such as carbon in a loaded MEA liquid, at the size its
rounding allows.
Before evaluating a trial, a fresh undamped correction may terminate polishing
when the current point is the stored best point, its classified maximum is within
the unchanged caller criterion, and every free phase amount is positive (#272).
Let $D$ be diagonal with positive free phase amounts in mol, free $T$ in K,
free $P$ in Pa, and unity for log compositions, log packing fractions and
multipliers in RT units. The numerical-policy owner sets
$\varepsilon_x =$ `EQUILIBRIUM_POLISH_CORRECTION` $=10^{-12}$.
The scaled correction must satisfy $\|D^{-1}\delta\|_\infty\le\varepsilon_x$.
The original Newton factorization must have full numerical rank and finite
$\delta$. With $A=JD$, the transpose factorization must also have full
numerical rank, and its Eigen reciprocal-condition estimate $\widehat{rcond}_\infty(A)$
must give a finite positive denominator in

$$
\|D^{-1}\delta\|_\infty +
\frac{\|J\delta+r\|_\infty}
{\widehat{rcond}_\infty(A)\,\|A\|_\infty}
\le \varepsilon_x.
$$

Here $J\delta=-r$; the matrix norm is the induced infinity norm. Eigen's
1-norm estimate on $A^T$ supplies the infinity-norm condition estimate.
The solve-admission calculation runs only when the correction-size guards
would permit the shortcut. A nonpositive free phase amount disables it;
amounts are scaled by their own phase amount, never by feed total.
Detailed `polish_correction` capture records the scaled infinity norm, its
controlling chart coordinate, maximum feed-relative free-amount correction,
forward estimate, classified maximum, current-best and positive-amount guards,
checked Newton rank (−1 when not checked), and shortcut decision. The scaled
norm includes each free amount relative to its own phase. Payloads are built
inside `Scope::observe`; `polish_stop` names a successful correction stop.
This local estimate does not certify total state error or every raw row within
the caller criterion. Normalizing compositions and converting packing to density
can amplify a chart correction, the condition number is estimated, and later
Newton corrections are not bounded by the current correction. Paired installed-wheel
comparisons qualify state preservation only on their stated inputs.
If any guard fails, the existing polish trials and best-point selection continue.

Try fractions 1, 1/2 and 1/4 of each Newton direction, stopping at the first
admissible raw contraction; a rejected direction ends polishing. Among the
optimizer start and accepted physical iterates, return the point with the lowest
classified residual, breaking ties by lower raw residual (#66, 2026-09-22).
The classified residual may rise at an intermediate step through nonlinear
curvature; requiring it to decrease at every step can reject a direction that
resolves trace chemical equilibrium within the step budget.
Selecting the best classified point prevents degradation of the optimizer
result or loss of a stricter caller tolerance attained along the path.
For eligible calls entering `solve_compiled`, the [ADR 0009 amendments](../../../docs/scientific/decisions/0009-general-platform-one-driver.md)
apply the same exact-Jacobian Newton iteration from the existing start before
constructing Ipopt. #476 retains square feasibility cells with one liquid and
one vapor, one zero-pinned and one free amount, exactly one free intensive,
and no discovery, reaction, ion, caloric, critical or envelope condition.
#545 adds square stationarity systems at fixed positive T/P with two free
amounts, liquid/liquid or liquid/vapor phases, and identical supports equal to
the retained positive-feed support of at least two components. Every model
charge must be zero. Source-reference conversion, diagnostic options and
optimizer-state transfer exclude this added domain, as do the other #476
exclusions.

For discovery declarations, #828 extends the same free-amount attempt to one
free fluid phase on exactly the positive-feed support and to the existing
admitted two-phase support structure. Non-ideal square fixed-positive-T/P
systems alone qualify; ions, reactions, source-reference conversion, solids,
caloric, critical, envelope, azeotrope, diagnostic-options, vapor/vapor and
three-phase declarations remain excluded. One shared operation owns the
attempt and density-branch admission. A rejected attempt leaves the original
seed for unchanged Ipopt. Discovery retains all donor, objective-decrease and
final phase-set checks. Its work map counts actual Ipopt calls separately from
`declaration_newton_attempts`, retaining rejected-attempt work too.

The free-amount attempts check the finite seed against the actual physical
primal box and topology before preparing iterations, then requires finite
Jacobians/corrections and full numerical LU rank. Finite amount bounds from
`GibbsMinimizer::get_bounds_info` are multiplied by the existing amount chart;
infinite bounds retain their meaning. Candidate boxes precede preparation.
The free-amount attempts permit a trial whose raw rows already meet the caller criterion
past the falling-contraction-ratio test. #476's iteration remains unchanged.
In both domains, once the raw criterion is met, one compulsory full exact
correction uses a fresh Jacobian and solve. The corrected point must meet the
fresh classified criterion. At Newton box checks, log-composition overshoots
within four ulps of physical x=1 are snapped to zero, using four successive
`nextafter` steps above one (8.881784197001248e-16 in log coordinates), then
re-prepared. Every other bound and criterion remains unchanged.

A returned point must end before the cap, satisfy the box and unchanged
classified rows without a floor allowance, and pass topology, local checks
and density-branch admission. Free-amount points also pass the existing exact
reduced-Lagrangian-Hessian check on the balance/normalization tangent space;
chemical and pressure stationarity rows are not primal constraints. The
existing reconstructed extensive-objective multipliers below apply, including
`-N_feed lambda_j` for material rows. Rank and strict positive reduced
curvature refuse degenerate amount splits, including equal-composition
azeotropic splits; composition equality alone is not phase coincidence or a
global stability test.

Otherwise, including on admission failure or exception, unchanged Ipopt runs
from the original start and retains attempted work and the shared budget.
A #545 Newton return keeps `solver_status` `SUCCESS`; the existing
`newton_fast_path` trace and its accepted outcome identify Newton. #476 retains
`NEWTON_CONVERGED`. Both have zero optimizer iterations, optimizer Hessian
callbacks and regularization amplification.
Zero constraint multipliers of row-count length apply only to #476 feasibility
cells. #545 still evaluates exact final Gibbs curvature. Evaluations count
Newton prepared points, while trace steps separately report Newton work.
Solved-state actions consume the same coordinates and rows. Declared starts,
continuation correctors, stability trials and caloric seed sub-solves retain
Ipopt alone.

The #771 explicit charged density-temperature Vapor correction uses the log
charge row and exact log-coordinate Jacobian in EqID `tangent_plane_trial_balance`.
The same pair governs Newton correction, contraction, damping, best-point
selection, stopping and the local forward estimate. Classification retains
both the original physical criterion and the log-charge criterion. Endpoint
physical and numerical ranks are recorded separately; full numerical rank is
necessary on this path, including trial polishing. Ipopt callbacks, physical
duals and solved-state derivative actions retain their physical equations.
Only for `refine_trial && !compiled.log_charge_rows.empty()`, admit the starting
point before best-point retention and each otherwise merit-admissible proposal
through `check_density_branches` with empty anchors and `initial=true`. The
inspection copy subtracts the log of each phase's supported-coordinate sum
from its log coordinates at unchanged packing. It neither repairs charge nor
changes the raw proposal. The existing declared root, convergence, positive
slope, density window and phase-order checks own admission. A rejected
proposal retains its typed cause and enters the existing damping/recovery;
a rejected starting point returns an unresolved attempt with its rows and work.
The same work/deadline limits apply. Ordinary solves, Liquid and neutral trials,
and supplied-point evaluation skip this inspection and its root work.
This representation does not qualify relative trace derivatives.

Acceptable-point, tiny-step and budget stops use the same polish; a budget stop at an
unattainable request would otherwise return the optimizer's last iterate. A selected
point whose classified rows meet the caller criterion is accepted after any
non-error optimizer stop (success, acceptable point, tiny step, square feasible point,
restoration failure or budget); invalid-number, diverging, step-computation and infeasibility stops stay
unavailable, and a floor-limited point also needs optimizer success. A failed Ipopt
initialization, a fatal application return (invalid problem or option, exception,
memory, internal error) or an unassigned status is not a stop: it leaves no
candidate to classify and is recorded `infeasible` like an unrun structural solve
(or `left_domain` when the fallback start's reference is unavailable), even when the
fallback start's rows meet the criterion (#180). Its `solver_status`
carries both exact Ipopt names, e.g. `UNASSIGNED; Initialize returned Invalid_Option`,
so a solver failure stays distinct from a `LOCAL_INFEASIBILITY` stop; the message repeats them. One non-error rule gates the polish, row classification and acceptance. No additional step or
relaxed criterion is used. Recompute diagnostics at the selected point, and
apply the unchanged caller criterion and all final physical/branch checks.
A feasibility cell carries chemical and pressure rows as optimizer equality
constraints. Its log-composition lower bound is therefore the representability limit,
the smallest normal double (#322): the all-free optimizer box or a composition cutoff would exclude valid
trace roots before a successful optimizer stop can reach polishing (a reactive bubble point
whose trace root lies at 2.4e-24 stopped `infeasible` under the former 1e-20 bound). That bound
moves feasibility-cell values at rounding level with no work change, except conditioning-limited
near-critical solves; the measured moves are retained with #322 (PR #332). In all-free cells a warm start below the box
is pulled back into it by the optimizer's bound push and re-recovered by the polish.
The retained MEA 050 case exercises rejection of a full step and acceptance of
a half step without relaxing the requested tolerance. The returned extensive-objective constraint multipliers are
reconstructed at the returned point from the row multipliers: a material row
`(B n - b)/N_feed` carries `-N_feed lambda_j`, a charge row `sum z x` carries
`-n_p phi_p`, and a free-amount phase's normalization row is zero because the
objective and material rows obey joint extensive homogeneity in the phase's
amount, composition and volume. Together with stationarity in its free amount,
this gives zero normalization multiplier; the reduced-Hessian certificate uses these values at the returned
point, and the optimizer's own multipliers remain a trace probe (#87 review,
2026-09-21). A feasibility problem keeps the optimizer's constraint multipliers
over every row; its equality-root qualification does not assert a Gibbs minimum.
The reconstruction assumes no active
variable bound at the returned point and residual potentials that depend on
composition and molar volume only. The charge-row scaling `-n_p phi_p` is
qualified by independent extensive-objective stationarity and public ionic
chemical potentials at Khudaida Figure 2 point 1 with 1 and 10 mol feeds
(`test_trace_reliability.py`); application-specific ionic qualification remains
open. The per-phase charge row is met to the absolute row tolerance
(`EQUILIBRIUM_ROW_RESIDUAL`, 1e-10, unless the caller sets `row_tolerance`). In a phase whose
ion content lies far below that tolerance, the row does not fix how the trace ions divide
between cations and anions. Only their electroneutral combinations, which the phase-equilibrium
rows fix, are determined: for example the mean-ionic pair. Individual ion amounts and
single-ion quantities of such a phase (single-ion activity, phase pH) are then undetermined
and are not results. In the Ascani 2022 acetic acid + water + toluene reproduction, the toluene
liquid returned Ac⁻ at a mole fraction of about 1e-21 and H⁺ at about 1e-36, while the mean-ionic
balance held to 2.2e-10 (#374, #375). The Ipopt driver
scales its extensive Gibbs objective by `1 / feed_total` internally; this
preserves the minimizer while keeping objective derivatives on the same
numerical scale as the normalized material-balance rows. It accepts Ipopt's
tiny-step or restoration-produced acceptable-point termination only when the
recomputed residual meets the declared equilibrium criterion. Neither status
alone establishes success; final topology, curvature and branch admission still
apply. Ipopt unscales
final constraint multipliers before the
TNLP callback receives them, and the Engine preserves that original extensive
Gibbs-objective meaning.
The compiled preparation enforces each phase's packing domain
[`PACKING_FRACTION_FLOOR`, `PACKING_FRACTION_CEILING`], so Ipopt cuts back a step
that leaves it. A Gibbs minimization of exactly one declared phase leaves ln η
unbounded in the optimizer, because the barrier of bounds on ln η pushes a
near-spinodal liquid off its branch (#314); for such a solve the packing domain
is a hard refusal that names the phase, so a state outside it (for example a
dilute vapor below the floor) is refused rather than returned. Every multi-phase
solve and every feasibility cell keeps the ln η bounds: their barrier acts as a
trust region on liquid density, and without it mid-dome methane/ethane flashes
fall to the trivial solution or exhaust the iteration budget, and the envelope's
critical-straddle corrector lands on the trivial solution. Those solves can
still show the barrier bias near a liquid spinodal; a near-critical declared
liquid inside a multi-phase declaration is not covered. Removing the barrier
changes some single-phase outcomes: identical-isomer isomerisation started at its
symmetric stationary point (a Gibbs maximum along the reaction) now refuses as not
positive definite instead of escaping to a minimum through the barrier's
perturbation, the outcome the #159 test already designates.

<a id="equilibrium-amount-step-retention"></a>

Phase amounts enter the Gibbs objective linearly. At fixed composition and density
the amount–amount second derivative is zero, so only the amount–composition,
amount–density and constraint couplings of the full Newton system limit how much
material one step moves between phases. At the #356 MEA point (loading 0.2556,
318.15 K, 3699.08 Pa) one step from the feed start moved 140 mol into the vapor
against an 86.7 mol feed. Ipopt's own fraction-to-boundary rule ($\tau\ge0.99$;
Wächter and Biegler, 2006) left 1% of the liquid, its composition drained to an
ion-bearing state and the solve ended "phase V vanished". Every optimizer trial
point is therefore refused unless, for each free phase amount with $b-l>0$,

$$x-l\ge\kappa_N\,(b-l),\qquad \kappa_N=\texttt{EQUILIBRIUM\_AMOUNT\_STEP\_RETENTION}=0.1 .$$

Here $b$ is the amount at the base of the direction being searched (Ipopt's current
iterate, which a watchdog rollback restores; in restoration the original-variable
block of the restoration iterate), $l$ the active search problem's relaxed lower
bound (restoration's own during restoration; Ipopt relaxes the zero bound to about
$-10^{-8}$, that is $-10^{-8}S_p$ mol) and $x$ the trial amount, all in Ipopt's amount
coordinate $y_p=n_p/S_p$ defined below. The rule is unchanged by a positive per-column
scale, so it reads the same in mol. Along one shrinking
direction $d$ this is $\alpha\le(1-\kappa_N)(b-l)/(-d)$: a fraction-to-boundary rule
with $\tau_N=0.9$ on the free amounts only, the tighter of Ipopt's $\tau$ and 0.9
applying to a shrinking amount. It adds no constraint, bound or amount floor and
starts no solve. Interior stationary points are unchanged, but the basin a solve
reaches can change. A refusal is an evaluation error, not an exact step cap:
ordinary backtracking halves the step, soft restoration repeats the refused point up
to three times before full restoration, and a refusal can suppress a second-order
correction or end a watchdog attempt. Refused objective callbacks count as
evaluations; refused constraint callbacks do not. The rule keeps the distance from
the relaxed bound positive, not the amount: with $b=10^{-9}$ and $l=-10^{-8}$ (units of $S_p$),
$x=-8.9\times10^{-9}$ is admitted, and the topology, balance and density checks
still decide the returned state. Measured on the amount itself ($x\ge\kappa_N b$) the
rule refuses every step toward the barrier optimum, which lies between the relaxed
bound and zero: the neutral binary liquid at 0.8 times its dew pressure then stops
at the 200-iteration budget instead of "phase L vanished" (E4-C pin). Along accepted steps with one fixed $l$ and no rollback,
$N_m-l\ge\kappa_N^m(N_0-l)$, so reaching $q=\epsilon_N N_{feed}$ with $l<q<N_0$ takes
at least $\lceil\log[(q-l)/(N_0-l)]/\log\kappa_N\rceil$ accepted steps. That is a lower
bound on steps; the added iterations have no upper bound. Genuine disappearances keep
their classification at a cost: C3 (NaCl liquid, water-only vapor) 27 to 97
iterations in both phase orders, the neutral liquid case 42 to 50, the neutral vapor
case 120 to 113, and E4-C reaches the same 0.6075 with 2,172 to 2,181 optimizer
callbacks. Over the 58 MEA family points accepted by both starts, 12 need more of
some work component (worst 49/69/50 to 145/446/144) while the sums fall from
2778/4388/2835 to 1866/3479/1922. The qualification trace found the identified base
and active bound at all 16,161 admissions over 2,540 callbacks, including restoration
and watchdog rollback. The value 0.1 is provisional (#356). The rule reads Ipopt
3.11.9 internals (the TNLP adapter's ordering, the current iterate and the active
problem's bounds; IpBacktrackingLineSearch.cpp, IpRestoIpoptNLP.cpp). A missing
mapping or Ipopt's own variable scaling, which would put its iterate and bounds outside the
callback coordinates, stops the solve with a named error, and a compile-time
check refuses another Ipopt release, which must requalify the rule (3.14's
`get_curr_iterate` returns iterates, not these relaxed bounds).
Only `gibbs_minimizer.cpp` includes the private algorithm headers and dereferences
these objects. Its declaration keeps opaque types; destruction is defined in
the implementation so the reference-counted live iterate owner remains intact.
The #821 qualification rejected 3.14.20 on the unchanged near-critical
methane/ethane comparisons with MUMPS and MA57. The deployed dependency remains
exactly 3.11.9; confinement does not remove the private dependency or qualify
another release.

<a id="equilibrium-start-amount-scaling"></a>

The primary Ipopt attempt measures each free phase amount of a Gibbs minimization in units of its own start amount (#243): it works in $y_p=n_p/S_p$ with
$S_p=\max(n_p^0,\epsilon_N N_{\mathrm{feed}})$, $\epsilon_N=$`PHASE_DISAPPEARANCE_AMOUNT`$=10^{-8}$,
fixed for one solve after every start owner has run; every other coordinate keeps unit one,
and feasibility cells keep physical coordinates. The callbacks return $S\nabla f$, $J S$ and
$S H S$ at the physical point $n=Sy$. Objective, constraints and their multipliers, starts,
bounds, Ipopt options, polish, classification, continuation and sensitivity are unchanged.
In mol the solve was conditioned by the feed size. Ipopt starts every bound multiplier at the
constant `bound_mult_init_val` 1, so in feed-relative units $q=n/N$ the initial amount bound
multiplier is $z_q=N z_n=N$, while composition and packing multipliers do not depend on $N$.
This mismatch is demonstrated but no experiment has isolated it as the sole cause; what matched
across feeds was the normalized iteration-0 start, not the iteration-1 point. In start-amount
units the #312 feed-scaled declaration matrix passes at 1, 100 and 10,000 mol with feed-independent
optimizer iterations (PR #400). The unit was chosen empirically, as the only qualified candidate
that also adds no key-solve iterations: feed-relative $n/N$ passes the matrix but raises the
MEA flash iterations; other order-one multiples of the
start add key-solve iterations or swap the neutral phase labels; Ipopt's gradient-based
`nlp_scaling` does not scale variables; user variable scaling would conflict with the guard
above; `mu-based` bound-multiplier initialization is untested. The floor regularizes the
coordinate scale only. It does not clamp $n$, add an amount bound or classify a predictor as
vanished, and a below-cutoff predictor does not establish physical disappearance (a
continuation corrector can predict an amount at or below zero; without the floor E4-C stops
nonfinite at path 0.6075). The floor is homogeneous in the feed, so the iteration is invariant
under positive feed scaling except for the $\max(1,N_{\mathrm{feed}})$ dual-tolerance clip
below 1 mol and Ipopt's absolute `compl_inf_tol`. C3 still refuses its vanished vapor; its
convergence category is now `stalled` instead of `converged_to_floor`.

For a reacting Gibbs declaration whose selected primary attempt is rejected by `run_driver_once` for `PhaseDisappeared` after its existing non-error stop decision, the shared driver may make one additional invocation with physical-mol amount coordinates, $S=I$. `check_topology` alone owns disappearance classification. Successful attempts, nonreacting declarations, feasibility cells, solids, critical candidates, error/fatal stops and other admission refusals do not enter this recovery. Ordinary solves rebuild the original declaration's start after their existing attempts; discovery uses its conserved declaration coordinates and spends another declaration-solve slot; explicit starts and continuation use their declared or predicted coordinates. The physical-mol invocation carries no optimizer multipliers and uses the same equations, exact derivatives, options, polish and physical admission. It never restarts from a vanished endpoint and never retries itself. All caller limits remain shared and are checked again after preparation. Retain primary, preparation and recovery work and both causes on failure, including exceptions; an unavailable public result has no usable physical-state arrays. If discovery has no new raw endpoint, retain its primary point only for a budget-permitted conserved removal and keep the failed invocation's raw-point record empty. Remaining disappearance is returned to the existing caller; recovery supplies no physical-absence or global-stability conclusion.

The one Ipopt driver uses exact weighted Lagrangian curvature for every admitted
specification, including pinned/free-intensive and incipient bubble/dew cases.
Its exact Hessian callback must obey EqID: `weighted_lagrangian_hessian` at
arbitrary callback points, objective factors and optimizer multipliers. The
objective factor scales objective curvature only. Pinned nonlinear pressure and
potential constraints retain their underlying curvature when an objective is
constant; a zero multiplier makes that constraint’s weighted term zero, without
removing the required constraint Hessian. Root and limited-memory equilibrium
replacements are not admitted, and a specification that cannot be represented
by this shared exact-curvature formulation is unavailable with a named
diagnostic. Derive all required third energy derivatives on the declared regular
domain; lifting variables does not automatically remove that requirement. This
is a capability requirement, not a claim that the current implementation
satisfies it.

<a id="ipopt-diagnostic-options"></a>
The shared Ipopt driver ignores an incidental `ipopt.opt` in the working directory.
For deliberate diagnostics, set `EPCSAFT_IPOPT_OPTIONS_FILE` to an options-file path;
the driver passes that path to `Initialize` for both exact-Hessian equilibrium and
limited-memory parameter fitting. A nonempty path must name a readable regular
file; otherwise the driver reports an error naming the variable and path before
initialization. An unset or empty variable disables file loading.
Relative paths resolve from the working directory. The named file can override
the driver's solver options, including iteration limits; unset the variable after
the experiment. Iterate logging can use `output_file <log-path>` and
`file_print_level 12` in that file.

The same owner chooses Ipopt's linear solver once per loaded native library
(owner decisions 2026-10-16, 2026-10-17 and 2026-10-18). `EPCSAFT_IPOPT_LINEAR_SOLVER=mumps`
or `ma57` forces a solver; any other nonempty value, or a forced `ma57` that
cannot load, is an error at the first Ipopt solve or `_core.linear_solver()`
call, before that solve starts; a solve the free-amount Newton path finishes
never reads the variable. Otherwise Ipopt's own HSL loader
opens `EPCSAFT_HSL_LIBRARY` (unset or empty: `libhsl.so` on the runtime linker
path), and the driver uses MA57 with AMD ordering (`ma57_pivot_order 2`) when
that library provides it, and MUMPS otherwise. MA57 runs with its own scaling
(`ma57_automatic_scaling yes`, ICNTL(15)) because, unscaled, a trace moiety row
reads as dependent; an options file can override it. HSL is user-supplied under its
own licence and is never bundled; work whose results used MA57 acknowledges
"HSL, a collection of Fortran codes for large-scale scientific computation. See
https://www.hsl.rl.ac.uk/". An options file still overrides the choice. Each
equilibrium and minimax result's `linear_solver` reads back the solver Ipopt
ran (empty when the Newton path finished without Ipopt), and
`_core.linear_solver()` records the process choice that fits and predictions use.

The caller supplies one dimensionless row tolerance per solve request; the
named `EQUILIBRIUM_ROW_RESIDUAL` policy is its only default. That value alone
classifies the recomputed rows; the Ipopt `tol`, `acceptable_tol` and
`constr_viol_tol` settings are the resolved request
times `EQUILIBRIUM_INNER_TOLERANCE_FACTOR`, because the
recomputed row maximum tracks the optimizer tolerance and overshoots it by a
hair when both coincide. A feasibility cell (pinned amount, caloric or critical
condition: constant objective, followed by the exact Newton polish) stops at the criterion itself; the polish, not an
inner decade, removes the overshoot (owner decision 2026-09-23).
The shared Ipopt option owner sets the absolute unscaled `dual_inf_tol` to
that inner tolerance times $\max(1,N_{\mathrm{feed}})$, with feed total in mol
(owner decision 2026-10-06, #243). It computes this factor once from the
existing objective factor $1/N_{\mathrm{feed}}$; the scaled overall `tol`
and material-row tolerance retain their values. The parameter-fitting call's
objective factor is one, so its dual stop retains its value. Ipopt's
counter-based acceptable stop remains disabled (`acceptable_iter=0`).
Both the inner target and its feed-scaled dual target must be finite and positive;
an unrepresentable combined request is rejected before Ipopt initialization.
The result's `solver_settings` reports consumed `tol`, `acceptable_tol`,
`constr_viol_tol` and `dual_inf_tol`, after the diagnostic file has been read;
Newton-only results have no consumed Ipopt settings. No starts, bounds,
barrier settings, complementarity tolerance or iteration budget change.
A caller-valued criterion is owned
by the convergence record so copies into results, traces and exceptions keep
it alive. The result reports `requested_tolerance_met` exactly as recomputed
residual ≤ requested tolerance whenever a finite residual exists; an available
floor-limited state may miss the request and keeps its state with a warning.
A declaration's generic-point Jacobian is tested with Eigen full-pivot LU after
dividing each nonzero column by its infinity norm; zero columns remain zero.
The default relative pivot threshold is unchanged. This removes dependence on
free-variable units up to floating-point rounding, without rescaling equations.
A structurally square declaration whose normalized Jacobian loses rank
returns `NOT_RUN` with `Stop::infeasible`, an unmeasured residual, no
criterion and the row, unknown, rank and nullity counts, before any start or
optimizer work. Within one solve the driver prepares each distinct callback
point once and reuses its phase evaluation for rows, Jacobian, gradient and
curvature; the weighted constraint curvature is recombined when the optimizer
multipliers change, and the final row certification is a fresh evaluation.
Optimizer callbacks, multiplier recovery and Newton polishing prepare phase
Hessians only (ALG-ASSOCIATION); their rows, Jacobian, gradient and curvature read
only V, support-amount and free-T entries. The final certification, solved-state
sensitivities, lever Jacobians and boundary checks use full preparation, so the
certified floor and reduced-Hessian certificate keep every coordinate.
That fresh final preparation also supplies the reduced-Hessian certificate;
only the coordinate width is reduced to the primal chart, without reevaluating
the same phase energies separately for its Jacobian and Hessian.
The certificate decides on inertia, not on the raw eigenvalue (#159). In the log
chart a trace species' direction has curvature of order its amount `n x_i`, so
on a trace state the smallest eigenvalue of `Z^T H Z` falls below the Hessian's
own rounding: in the packet-3 MEA liquid with CO2 at `x ≈ 1e-15` it is about
`1e-13`, while `||H|| ≈ 2.5e3` (the volume coordinate) and `eps ||H|| ≈ 5e-13`,
so its computed sign is noise. Let `H = H_g + H_c` be the objective and
multiplier-weighted constraint curvatures, whose composition diagonals cancel at
stationarity, `S = |H_g| + |H_c|` their split magnitude and
`D = diag(S)^{-1/2}` (1 where `S_ii = 0`, which includes every free phase-amount
coordinate, so `kappa(A D)` and the band depend on the amount unit). The substitution `z = D y` maps the
tangent space `A z = 0` onto `A D y = 0`, so the reduced Hessians of `H` and
`D H D` are congruent and have the same inertia (Sylvester's law);
equilibration changes the conditioning only. A first-order rounding model
bounds each assembled entry by `|dH| <= n eps S`. By Weyl's inequality with an
orthonormal nullspace basis this moves the smallest equilibrated eigenvalue by
at most `n eps ||D S D||`. The basis from the singular value decomposition of
`A D` lies within an angle `eps kappa(A D)` of the exact one, which adds
`2 eps kappa(A D) ||D H D||`. The band is

$$
\tau=\varepsilon\left(n\,\lVert DSD\rVert_F+2\,\kappa(AD)\,\lVert DHD\rVert_F\right),
$$

where Frobenius norms bound the spectral norms from above. `lambda_min(D H D) > tau`
certifies the strict local minimum, and `< -tau` rejects it ("not positive
definite"). A value inside `|lambda| <= tau` is the typed `Stop::singular`
outcome ("is indeterminate"), with no state: at working precision the inertia
is unresolved. `reduced_hessian_minimum` stays the unscaled smallest eigenvalue,
the physical log-chart curvature (`N x_B` for a trace isomer). Its sign carries
information only above its own rounding. The band is a first-order estimate,
not a rigorous bound. It covers rounding of the assembled Hessian and nullspace,
but not cancellation inside the objective curvature sum or the returned point's
own solve error. That error generally exceeds `tau`: at a row residual near
`1e-12` the equilibrated minimum moves by order `1e-13`. It decides the outcome
only when `|lambda_min(D H D)|` is below it, which happens near a genuinely
degenerate state. Take two identical ethane pseudo-species with an isomerisation
at 200 K and 5 MPa, at the `k_ij` where the symmetric state's curvature vanishes.
Adjacent `k_ij` values give equilibrated minima of about `+5e-14` (row residual
`1.3e-12`) and `-9e-14` (row residual `9e-16`) against a band of `5e-15`. There
the decision follows the returned point, not the exact equilibrium, and the typed
indeterminate outcome is reached only when the solve error is at rounding level. Diagnostic
residual summaries reuse the certified row values. For a no-site EOS Hessian,
each unit-direction `Forward(1)`/`Reverse(2)` sweep supplies a Hessian column;
the first sweep's even entries supply the gradient, including active-parameter
coordinates. The initial value sweep and the association implicit chain remain
unchanged; a separate no-site first-order reverse sweep is unnecessary.
Standard-state terms carry exact value, first and second temperature
derivatives from the reaction target owner; no affine shortcut is assumed.
For an ordinary free-amount cell the chemical rows are not optimizer constraints: they are
recovered from the objective gradient `N x_i mu_i` through least-squares balance
multipliers, so their unscaled residual is bounded only by the dual tolerance
divided by `N x_i` and diverges where a phase or a fraction vanishes. The driver
therefore classifies the free-amount route by the KKT-scaled stationarity
residuals `N_p x_i mu_i / N_feed` and `N_p / N_feed` times the pressure row (the
first-order conditions in the solved chart; owner decision 2026-09-16), keeps
the unscaled rows in `residuals` as diagnostics, and classifies a feasibility cell by
its unscaled rows, which are its constraints. In both routes a trace moiety row is
classified relative to its size $s_k$ (polish, above). Species thermochemistry attached to a
model is compared with every target kind, including scalar `K` and affinity,
and the comparison is skipped where the temperature leaves a species interval.

A converged state carries no certificate that each phase is locally stable.
For a neutral nonreacting liquid/vapor pair with an incipient phase pinned at
zero amount, the driver applies each phase's necessary local stability
condition to a state it would otherwise accept (#352). This is the Gibbs
criterion in the Helmholtz form of the `critical_equations` sources
(Heidemann–Khalil 1980; Michelsen 1980). The bulk and incipient phases are both
checked, each at unit amount, so the pinned zero amount never weights the
matrix. With `Q = F_nn` at (T, V), the ideal `delta_ij/n_i` plus the prepared
residual Hessian, `D = diag(sqrt x)` on the phase's support, and `C` its
supported species count, the state is refused as `local phase instability`
(naming the phase, `lambda_min` and `b`) when

$$
\lambda_{\min}(DQD)<-b,\qquad b=\lVert DQD\,v-\lambda_{\min}v\rVert_2+C\,\varepsilon\,\lVert DQD\rVert_F ,
$$

with `v` the normalized eigenvector. A value inside the band is unresolved and
admitted; no positive margin is required. The band covers eigenvalue rounding
only, not error in the returned state or the equation of state. Associating
mixtures are admitted with this band alone: the 12 MEA–water bubble and dew
states checked have `lambda_min` 0.88 to 1.28, but near-critical associating
states, where `|lambda_min|` approaches the association residual floor, are
untested. The check follows the convergence
and `phases_coincide` classification on every route through the driver: plain
and explicit starts and continuation correctors. Charged and reacting cells need their
admissible perturbation subspaces first and are not checked. Positive curvature
certifies neither global stability nor another solution's absence, and a
refusal does not say that no bubble point exists. On the equimolar neutral
`bubble_at` returns of the #352 design review, near-trivial states have negative `lambda_min`
and genuine bubbles positive by orders of magnitude beyond the band (PR #363).
E4-D pins the 260 K, 2.02606 MPa guess, which returned the trivial state as a success.

PH/PS adds one extensive caloric row to the same feasibility driver (k + e > 0).
Its objective is constant; all chemical, pressure, material and caloric rows
are constraints, with exact optimizer-weighted curvature. Actual pins retain
their distinct amount meaning. H/S values use E1, Jacobians E2, curvature E3
and third solved-state actions E4. The analytic geometry differentiates ideal
calorics, amounts, raw T, packing and moving feed normalization; residual energy
contractions remain with Properties. The target controls are H [J] and S [J/K],
with partial rows −1/(N_feed RT) and −1/(N_feed R), including their higher chains.
A nonzero absent-condition control is invalid; zero is permitted for either kind.

An ordinary fixed-T seed uses the same driver and frozen C_f = Σ N_p cp,p:
T1 = T0 + (H_spec − H0)/C_f, or T1 = T0 exp((S_spec − S0)/C_f).
Only finite positive C_f and an admitted T1 change the seed. Pure two-phase
seeds instead use the pinned-zero saturation cell at P and the caloric lever
rule; an exterior target requires a caller-declared single phase. An auxiliary
failure retains its cause and the caller start; there is no outer retry. The seed
alone stops after `EQUILIBRIUM_CALORIC_SEED_ITERATIONS` (20) Ipopt iterations, within
any caller limits; one without a state by then is discarded like any failed seed.
A discarded mixture seed whose feed at $(T_\mathrm{seed}, P)$ has no liquid-side density root (two stable roots,
or a sole root denser than its negative-slope density) would start both phases near the feed at vapor
density, and the caloric row is then met by sensible cooling onto the trivial solution (#592). The start
instead is the fixed-$T$ generic start at the first $T_j=0.95^j\,T_\mathrm{seed}$ (`COOLING_RATIO`) whose
feed has a liquid side and which every caloric record admits, with no second seed (#610). No such $T_j$,
or a feed density solve that fails during the scan, keeps the start at $T_\mathrm{seed}$. Pure fluids keep
their saturation seed.
A start outside the declared phases' coexistence region otherwise spends 48–200
iterations on a seed that cannot succeed; converged seeds took 4–14 (#579).

Regular PH/PS qualification restricts the returned coordinates to a diagnostic
fixed-TP Gibbs chart, recovers its physical chemical duals and computes the
existing reduced inertia without another solve. Pure coexistence instead has
an interior lever-fixed equality root with distinct positive mechanical slopes,
pressure/chemical equality and conserved amounts. Its fixed-TP Gibbs split is
analytically flat: strict curvature and ordinary equilibrium Cp are unavailable.
Neither qualification establishes global stability. Straight continuation carries
H/S targets and controls through the same native predictor and corrector.

Declared-branch continuation follows one straight dimensionless path `s` in
`[0, 1]` between two problems with identical declarations. Fixed intensives
and resolved feed amounts may move, or a parameter-only path moves physical
EOS/Born coefficients as `theta(s) = theta0 + s (theta1 - theta0)` while T/P/feed
and all declarations remain fixed. `parameter_targets` uses the ordered active
coordinates of the supplied model; empty targets select the intensive/feed path.
Both paths use the same controller. Parameter trials rebuild the model from its
immutable base records and fixed thermochemistry; each compiled problem and
branch check uses that trial model for its complete lifetime. Accepted
`parameter_values` retain the physical coefficients. Born `c_dielectric` and
`ionic_region_relative_permittivity` are active scalars in the shared parameter
record; the latter must be positive. The tangent includes the full direct,
packing, association and source-reference parameter chain. A singular tangent
stops this straight path; it establishes neither a turning point nor absence of
an endpoint solution. Held-composition charts and arclength are outside this method. The controller predicts each step
from the implicit tangent (EqID: `implicit_sensitivity`: the first row
coefficient along the lever at the fixed state, the directional series of
ALG-ACTIONS, so source references carry their T/P chains), corrects with the same Ipopt driver seeded by the
full predicted coordinate vector (a predicted free amount below zero is
projected onto its bound, the disappearance boundary the bound-constrained
corrector then reports; a log mole fraction above zero is projected onto its
corrector bound zero; every other coordinate outside its
domain is a rejected left-domain step), and accepts a state only when the driver
accepts it, every recomputed row meets the caller `row_tolerance`, the
candidate-topology predicate reports no event, and each corrected phase density
shares the anchored locally unique root window (`DENSITY_ANCHOR_LOG_WINDOW`) of
its cold declared root at `s = 0` or of the previous accepted density
afterwards, with a positive corrected pressure slope. Finite step failures
(left domain, stalled, missing anchored root) halve the step from the last
accepted coordinates up to the caller's retry allowance and minimum step;
invalid paths, nonfinite operands, infeasible solves, a singular implicit
tangent at the accepted state (halving cannot change it), typed
disappearance/coalescence and a corrected density outside the branch window
are terminal. Passage through a critical point between consecutive accepted
states, where every nonzero component of a pair's inter-phase difference vector
(delta ln eta and each delta x_i) changes sign, is that pair's coalescence:
terminal `phases_coalesced` with the pair named and the swapped state retained
as the last result, never a relabeled continuation (owner decision 2026-09-16).
A partial prefix is never target success, a step or budget limit
is never a topology event, and every attempt is retained by the existing Core
trace scopes. Reusing only an initial seed is not continuation. The result is
a declared local branch, not a discovery or global-stability certificate.

### Polish and endpoint trace records

Detailed capture emits `polish_direction` (undamped correction infinity norm,
controlling coordinate, raw/classified defects before, extension, prior accepted
contraction ratio) and `polish_coordinate` (coordinate, phase, component, value),
whose reason names its chart quantity and units. `polish_trial` adds raw/classified
defects after, damping, acceptance, extension, prior ratio and trial index, retaining the actual rejection or preparation cause. Missing
trial quantities are NaN. `polish_stop` retains raw/classified defects, step and
the exit cause. Ordinary capture keeps the existing `polish_step` records.

`solved_endpoint_jacobian` reports defining-row count, solved-chart width, rank
and the scaled infinity-norm condition estimate after all deadlines, branch
admission and attempt selection. For the canonical defining residuals $F$ and
their exact Jacobian $J$, the endpoint calculation uses $A=J\operatorname{diag}(h)$.
The physical coordinate scale $h$ is one for log and multiplier coordinates,
the feed total for free phase and solid amounts, and $\max(|T|,1\ \mathrm{K})$ or
$\max(|P|,1\ \mathrm{Pa})$ for free temperature or pressure. `condition` is
$\|A\|_\infty\|A^{-1}\|_\infty$ and `inverse_norm` is
$\|A^{-1}\|_\infty$, not a norm of the unscaled inverse Jacobian.

The canonical row floor is
$$
\eta_i=\max\!\left(\gamma_n\sum_j|J_{ij}|h_j,
u_i+\epsilon t_i\right),\qquad
\gamma_n=\frac{n\epsilon}{1-n\epsilon},
$$
where $n$ is the solved-chart width, $\epsilon$ is the working-precision machine
epsilon, $u_i$ is the row's propagated inner error and $t_i$ its term scale.
Numerical row equilibration transports this floor by the absolute row scale.
The shared linear calculation combines the Newton correction, propagated row
floor and compensated linear defect through the exact physical-output
derivatives. Output estimates also include returned-versus-mapped displacement
and rounding. Temperature, pressure, amounts and compositions use absolute
radii; positive density and packing use log radii, exported to absolute radii
with `expm1`. Structural constants have no added rounding contribution.

Detailed capture or an explicit post-hoc request computes these estimates.
Unavailable rows, singularity or nonfinite projection retain their cause;
unresolved composition, amount, packing or phase-separation domains prevent a
conditional comparison. Source-reference chemical and first-thermal rows use
the local estimate described in [ALG-REFERENCE](#alg-reference); its returned
method annotation preserves the separate reference-replay precision and limits.
The [retained notebook](../../../analyses/2026-equilibrium-forward-error/notebook.qmd)
records the comparison coverage and limitations.

The native comparison first requires the original admission, residual, branch,
phase and inventory conditions and applies the ordinary allowance. Beyond it,
an estimated-error pass requires available, resolved estimates at both
endpoints, at least one absolute radius above that allowance, and a difference
no larger than the sum of radii in the quantity's chart. Positive density and
packing use log differences. Strict comparisons retain the ordinary allowance.
An ordinary pass never consults error estimates.

Default recorded rank may be structural or an existing solve decision; it is
not a continuously measured endpoint condition number. Failed detailed
inspection retains its cause and unavailable rank/condition. These local
linear estimates depend on coordinates and are no rigorous nonlinear error
bound or stability certificate. Their isolated EOS work is excluded from solve/search work
counts, but may update model caches/counters and adds uncapped time before return.
Detailed capture can advance finite deadlines; untimed returned-state identity
is qualified separately. Disabled capture adds no numerical work or clock reads.

Diagnostic records use ordinary construction and copies. Allocation failure
while constructing, copying or retaining a record propagates normally; it has
no incomplete-record fallback. An Ipopt log-open refusal leaves the solve
unchanged and sets `capture.trace.incomplete` with `ipopt_log_error`.
Trace-only containers and strings are constructed inside `Scope::observe`
or Core's initializer-list `probe` overload. The source rule in
`tools/check-engine-architecture.py` guards known eager-temporary and
named-local payload patterns during `check-all`.
Nested native work captures forward detailed records once with remapped
parents. Python declares capture and reads/serializes completed records only.

See the specification appendix for [cooling-onset dew](#cooling-onset-dew),
[critical-equation candidates and phase envelopes](#critical-equation-candidates-and-phase-envelopes),
and [pure stoichiometric solids](#pure-stoichiometric-solids).

The #906 saturation-loading search declares a positive increasing grid with
both interval endpoints. It solves ordinary free-pressure bubble nodes in order,
starting cold and then mapping the preceding admitted current-vector state.
Native pressure estimates resolve signs; native exact feed actions and the
augmented defining-row loading projection establish sampled local regularity.
Each strict sign bracket or pressure-uncertain node event receives one bounded
affine-loading correction through the same driver. A bracket uses its midpoint
and lower-node state; a node event uses that node and its adjacent-node span.
There is no scalar inversion, grid refinement or outer retry.

One existing recovery work envelope covers all grid, seed, fallback and
correction calls, including driver attempts. Candidate comparison uses all
loading, T/P and declared-order phase/state operands with their own native
estimates and zero comparison allowance. Every pair in a duplicate group must
agree; multiple strict brackets refuse independently of candidate deduplication.
Candidate corrections/comparisons precede the sampled slope audit. All nodes
and required events must be examined, with one admitted root group and no
refusal, before publication. “Unique within the declared search” excludes
between-node and undiscovered-branch claims. EqID `finite_addition_enthalpy`
owns the heat, moving denominator and switch derivative.

## ALG-STABILITY

The explicit native tangent-plane test and fixed-T/P discovery adopted in
[#223](https://github.com/tannerpolley/ePC-SAFT/issues/223) reuse `CompiledSystem`,
`run_driver`, density/association owners and exact CppAD value/derivative chains.
No production finite differences or Python search loop is admitted. One second minimizer is admitted (owner,
2026-10-19, #647): the native Michelsen trial minimizer below, for all-neutral, non-reacting systems only. It
changes no acceptance check. The [formulation](../equilibrium/formulation.md#discovery-and-global-claims)
owns physical admission and claim meaning; EqID: `tangent_plane_distance` owns
the reduced molar objective. The [analysis](../../../analyses/2026-tpd-phase-discovery/)
retains numerical qualification on the stated neutral and admitted
density-dependent ionic/reacting cases, explicit density-independent
charged-search refusals and unchanged declared work. Method adoption is recorded in
[ADR 0012](../../../docs/scientific/decisions/0012-native-tangent-plane-phase-discovery.md).

### Explicit calls and solved snapshots

`Discover()` dispatches natively before declared compilation, only at fixed
positive T/P. A legacy unknown-topology indicator is not permission to search.
`stability_test(...)` tests an accepted declared result without addition/removal.
Retain and verify model and complete physical declaration identity. At solved
T/P freeze actual standards, reference transfers and gauge. All-free chemical
multipliers are recovered; pinned multipliers are primal. Pinned/free-intensive
snapshots release amount pins for a fixed-T/P Gibbs-descent diagnosis, not for
following the original feasibility path. Only positive amounts supply donors.
Discovery with free intensives remains outside this method.

Declared solves, continuation, parameter fitting and ASL perform zero trial
solves. They preserve callback arithmetic/order and counted work. Returned
solved-state actions use the fixed discovered declaration without searching.
Qualification must compare the retained before/after declared numerical values
bit-for-bit, with identical iteration/evaluation/Hessian-call counts under the
same compiler, dependencies, inputs, options and thread limits. Report explicit
discovery/test wall time (median of three same-machine runs), driver roles and
trial counts separately; no universal wall-time threshold is inferred.
Explicit native searches own their work capture and pass its trace to the
counting code. The trace reports work only; no solver decision reads the ambient
trace. Declared calls retain their existing scopes and perform no trial work.

### Unfinished optimized trials

An ordinary explicit trial passes an internal refinement permission to the shared
driver. Every proposal of this polish, ordinary or recovery, limits its damping to
$\alpha \leftarrow \min(\alpha, L/\max|\Delta|)$ over the supported log mole fractions
and log packing fraction. Here $L$ = `EQUILIBRIUM_TRIAL_LOG_STEP` = $\ln 101$ is an
empirical retention policy (#811), equal to $\ln(1+f/\tau)$ for the $10^{-8}$
composition box and $10^{-10}$ row criterion; no accuracy bound follows from
that identity. At a rejected correction, recovery
starts only if the best classified point still exceeds the existing row criterion.
A row-acceptable trial, including a curvature refusal, keeps its original result.
Supplied-trial evaluation performs no optimization and grants no permission.
The scoped charged density-temperature Vapor trial also requires the normalized
declared-root inspection specified in ALG-EQUILIBRIUM before retaining the
starting point or accepting a merit-admissible correction. Rejection keeps
the raw iterate and uses the same damping and recovery limits below.

Recovery continues the same exact-row Newton loop on the same supported problem,
with fixed T/P and frozen standard/plane shift. It starts with damping scale $a=1$
and tries at most $a,a/2,a/4$ per correction. A prepared proposal must strictly
decrease the raw row maximum; the falling-ratio extension condition is suspended
only within this recovery. After accepting damping $\alpha$, set
$a=\min(1,2\alpha)$. If all three proposals fail, keep the current point and set
$a=a/8$ for the next correction. Entering recovery and unsuccessful corrections
consume slots in the original `EQUILIBRIUM_MAX_ITERATIONS` cap: at most 200
corrections and 600 candidate proposals, without a new solve or reset of work.
Existing deadline/evaluation limits, raw-row and correction stops remain; a
nonfinite direction or no representable progress stops recovery with its cause.

The shared driver keeps the best classified point and the unchanged final
row/curvature admission. Reaching a cap confers no success; detailed diagnostics
retain the cap or no-progress cause even if the best point passes final admission.
The existing positive-support charge reconstruction and reevaluated finite descent
then consume the returned raw point. Reconstruction does not revise optimizer
success. `StabilityTrial.composition` is that reconstructed composition and the
trial's state composition; `optimizer.mole_fractions` remains the raw optimizer
point, electroneutral only to the row criterion, and is not a property-state input.
This finite local refinement provides no global-stability certificate.

### Supplied trial points

`StabilityOptions.evaluate_supplied_trials` defaults to false. When true,
`stability_test(...)` requires a nonempty explicit `trials` family and keeps
its entries in caller order, including repetitions of a support/label pair.
Each full-species `composition_guess` is that entry's only seed. Require exact
structural zeros outside its support and finite positive supported entries;
apply EqID: `electroneutral_trial_reconstruction` before the log chart and
again after preparation. The returned composition is the point passed to
shared `phase_energy`, with unit trial amount and the frozen candidate plane.
No driver, composition minimization, descent classification or finite split
runs. `discover_phases` rejects this option before a declaration solve.

A supplied `packing_guess` provides eta, subject to the unchanged packing
domain; it is a coordinate replay and need not close pressure or certify a
branch. Without it, the unchanged density owner selects the Liquid/Vapor root
at the initially reconstructed seed. Both cases retain the chart's log/exp
rounding, reconstruct composition and preserve packing when computing molar
volume. No second root solve or different-root fallback follows reconstruction.
The existing quantity maps report actual `packing_fraction`,
`molar_density_mol_m3`, `pressure_row` and `input_reconstruction_l1` (the total
input-to-evaluated L1 displacement). `reconstruction_l1` still measures only
the final reconstruction. Raw residuals describe the supplied chart with its
unsolved generic multipliers. Optimizer-only dual/objective quantities are
absent. Omitted or false options preserve the ordinary search outputs.

Every available point calls the existing donor arithmetic, including
nonnegative TPD, and returns after its slope, correction and uncertainty
estimates. With a positive full-support donor its evidence is `evaluated`;
without one it is `evaluated_without_positive_donor`, with the cause
`no_positive_full_support_donor` and absent donor estimates. Packing,
reconstruction, density and energy failures retain unavailable ordered rows
and their causes, including `nonpositive/nonfinite trial composition` and
`stability point outside packing domain`. Under density-independent
permittivity the condensed-side check below applies at the supplied packing; a
refused charged point is an unavailable row with its cause. `completed` means the ordered evaluations were exhausted,
including unavailable rows; a trial budget can leave the sequence incomplete.
No stability or instability conclusion follows from supplied-point evaluation.

### Supports, seeds and trial solves

Let U be fed species plus all reaction participants and U₀ its neutral subset.
Initial discovery candidates are one free phase on U, Liquid, and also Vapor
when U is neutral. For reactions, first construct the existing positive
charge-neutral base on U, then run the unchanged ideal-solution speciation
sweeps. A zero-fed reaction participant remains on the required support; its
feed fraction is not an admissible log-chart seed. Record all initial attempts
and retain the accepted lower-G candidate; shared branch admission refuses a
candidate whose sole root lies on the other side (#423). This preparation belongs only to
discovery; declared initialization remains unchanged.
Trial defaults are (U,Liquid), (U₀,Liquid), (U₀,Vapor), deduplicating identical
pairs. Omit empty default U₀ with an explicit no-neutral-species diagnostic.
Caller support/label entries **replace** defaults and limit the tested domain;
they do not restrict the initial union candidate. An omitted caller support
resolves on U by the ALG-EQUILIBRIUM support rule, which also refuses explicit empty
supports; explicit ionic Vapor requires `DensityTemperature`. Reject empty caller lists, species outside U and one-sign
ionic supports.

Charged Liquid supports are admitted under every declared bulk-permittivity form.
`DensityTemperature` and `DipolarPerturbation` reach $\epsilon_r\to1$ as
$\rho_m\to0$ (full EOS §6.2/§6.4), and their charged trials count at every
density. Explicit charged Vapor under `DensityTemperature` must match its own
lowest stable pressure root after reconstruction, before energy/donor evaluation,
including supplied packing. The root must meet `DENSITY_ANCHOR_LOG_WINDOW`, the
existing pressure-root residual criterion and positive finite corrected slope.
A sole root established on the Liquid side is refused; a sole root without loop
evidence retains `Unspecified` admission. The generic sibling-root fallback is
unchanged for existing callers. Default trial families, discovery initial
candidates and its charged-Vapor addition veto are unchanged: no charged-Vapor
discovery is qualified. Density-independent `IonSuppression` (all choices) and
`ComponentMixing` are condensed-liquid correlations (#568). Under them a charged
trial point counts only when it coincides, within `DENSITY_ANCHOR_LOG_WINDOW`,
with a stable density root of its own composition at the test T and P, and the
density owner's side classification (`density_side`, #423) of that root is
Liquid. The seed applies the classification to its already-solved Liquid root.
The end point solves one Liquid-branch root with its sibling before the trial
energy, through the same root match that discovery addition uses; a supplied
packing is checked only there. A vapor-side root, a sole root without a loop
sample (`Unspecified`) and a point away from every stable root, including a
stalled point inside the unstable region, are unavailable with the cause
"charged trial outside the condensed domain of density-independent
permittivity". They are never witnesses or evidence of stability. Admitted
charged end points record the root's `negative_slope_density_mol_m3` when it has
one. Under density-independent forms, admission reads the trial state's side of
its own composition's mechanical loop through the density owner's
classification, not a phase label, fitted coefficient or optimizer termination.
It introduces no density cutoff, new tolerance or failed-trial filter.
Compositions with no loop at T are under-admitted, and admission establishes
condensed density only, not the correlation's fitted range. Declared solves,
continuation, fitting and ASL are unaffected.

For ordinary search, for every admitted pair enumerate a positive charge-neutral base from equal
supported weights with ionic group rescaling. A compatible candidate's strictly
positive reconstructed restriction may supply the near-pure remainder;
otherwise use the base, including species structurally absent from that candidate.
For each neutral species on a nonsingleton support make a 0.99 near-pure seed,
distributing 0.01 over the other supported species using that neutral remainder.
A singleton neutral seed is exactly one. Ionic seeds run on Liquid only. Neutral
Vapor also has an ideal-gas seed proportional to exp(π−s), shifted by its maximum
before exponentiation, without clipping or hidden floors. Include removed-phase
compositions on identical admitted support/label pairs. Reconstruct every seed
using EqID: `electroneutral_trial_reconstruction`. Density-owner Liquid/Vapor
roots supply packing starts; deduplicate identical seed/root pairs within
`DENSITY_ANCHOR_LOG_WINDOW`. Retain failed roots as unresolved trials.

For systems whose fed-plus-participant species are all neutral and that declare no reactions, minimize each
enumerated distinct start natively with Michelsen's modified tangent-plane objective (EqID:
`michelsen_modified_tangent_plane`). Take `TANGENT_PLANE_SUBSTITUTIONS` direct substitutions (Michelsen 1982 eq.
19), then second-order steps in αᵢ = 2√Wᵢ with the exact T,P Jacobian (eqs. 25–27). Take a step when its
Cholesky factor exists, α stays positive and g* does not increase; otherwise take a substitution step.
Each objective point requests the shared Helmholtz gradient. Only an attempted second-order step
requests the exact isothermal phase Hessian at that current composition and volume, without another
density solve. The first objective point reuses its seed's density root only when the full normalized
composition is exactly equal; every other point uses the trial label's unanchored density-owner root.
The root is local to that trial and is never reused at a changed composition. A trial converges when
‖ln Wᵢ + ψᵢ + sᵢ‖∞ ≤ `EQUILIBRIUM_ROW_RESIDUAL`, and its trial dual is k = −ln ΣW. A start that reaches
`EQUILIBRIUM_MAX_ITERATIONS`, throws, or leaves finite positive W runs through the shared driver unchanged. Run
every other enumerated distinct start, and every start of ionic or reacting systems, through the shared driver. Compile an
internal one-free-phase unit-amount system with reduced balance rows spanning
normalization/charge and no reactions (EqID: `tangent_plane_trial_balance`);
do not pass it through public no-reaction compilation that restores species
balances. Retain the echelon row transformation and transformed right-hand side.
Apply the frozen shift through existing standard-term callbacks, preserving
EqID: `tangent_plane_shift_log_curvature` at arbitrary points. Exact weighted
curvature, box admission, polish and local inertia remain shared. The normalized
trial affinity follows EqID: `normalized_trial_reaction_affinity`; the trial is
not required to satisfy the original zero-affinity rows.

### Reconstruction, uncertainty and finite descent

For ordinary search, retain raw coordinates only from an actual nonfatal finite, domain-valid
optimizer attempt, including a nonconverged attempt. Failed initialization has
no optimizer point. A Michelsen trial retains its converged primal chart [ln wᵢ, 1, ln η] as raw coordinates and has no raw
residuals, basis or row transformation. Keep raw residuals separate from reconstructed composition,
displacement and reevaluated energy. Normalize neutral trials; ionic rescaling
and the Liquid charge-only increase of the largest opposite-sign pivot followed
by common multiplicative normalization follow EqID:
`electroneutral_trial_reconstruction` with ties by index. Newly admitted charged
Vapor trials retain common normalization and multiplicative charge repair in
log arithmetic, without adding normalization roundoff to ionic coordinates. Recompute molar volume
from reconstructed composition and the retained packing coordinate. Reject
lost positivity/support, invalid T/P, packing or association domain. Record
actual normalization/charge defects and rounding scales; ionic defects must
meet the 1e-12 criterion. Never fabricate a species floor to admit a point.

For candidate balance/phase-charge defects let $m^0$ contain supported amounts
in mol, $J$ contain the conserved-moiety and phase-charge coefficients, and
$r=t-Jm^0$ for target $t=(b,0)$. Use one amount-weighted linear correction:

$$
D=\operatorname{diag}(m^0),\qquad a_k=\max_j|(JD)_{kj}|,\qquad
h=(SJD)^+Sr,\quad \delta=Dh,\quad m^1=\operatorname{fl}(m^0+\delta),
$$

where $S_{kk}=1/a_k$ for nonzero rows and 1 for zero rows. Divide each row
and its right-hand side directly by $a_k$; the existing complete orthogonal
decomposition handles redundant and zero charge rows. No squared inventories,
step shortening, iteration or backtracking is used. Attempt the full correction
even below tolerance; refuse nonfinite or nonpositive supported amounts and
incompatible original rows. This is local numerical reconstruction of a declared
parent, not a general constrained-feasibility solver.

Reconstruct phase amounts/compositions and check original material rows,
normalization and phase charge, including the reported log representation,
with the existing material/charge criteria. A below-tolerance residual need not
decrease strictly. Reevaluate physical energy before descent and keep this
candidate separate from the public result. Record actual stored displacement,
update rounding and physical/log representation changes using EqID:
`tangent_plane_descent_uncertainty`; every present phase and donor receives the
same final totals. For trials use normalization/charge rows and retain positive
support/domain checks. A reconstruction refusal retains its numerical cause;
it does not establish physical infeasibility.

Use the full present-phase chemical/pressure identity (EqID:
`present_phase_tangent_plane_identity`), not a chemical-residual-only threshold.
A donor must contain the entire trial support, have positive amount and finite
positive capacity (EqID: `tangent_plane_donor_descent`). Choose greatest capacity,
ties by phase index. Compute its actual slope directly from physical gradients,
recording plane value, chemical correction and pressure row separately.
Propagate EqID: `tangent_plane_descent_uncertainty` from reevaluated association
floors and site residual/conditioning, Euler energy-gradient defect, amount/volume
column scales, feasibility and arithmetic, plus supplied standard uncertainty.
A maximum native floor can conservatively replace individual column estimates;
it remains a local estimate. Final-sum rounding alone does not bound the EOS.

Negative plane and slope beyond the complete estimate require a finite split
(EqID: `conserved_finite_phase_split`). Try half the donor capacity, then at most
eight halvings. Preserve every component, hence every moiety and phase charge;
reconstruct arithmetic composition/charge, reclose the withdrawn donor at the
request pressure on its declared density branch through the shared density
owner (a sole root established on the other side fails the check, #423), recompute packing and shared energy, and verify matching conserved
totals before comparing Gibbs energies. The trial keeps its own molar volume,
and the donor's reclosed volume seeds any added declaration. A fixed donor
volume is admissible for the minimized A+PV, but its mechanical penalty
(g_vv>0 on a stable root) blocked useful finite descent within this bounded
search (#390). With exact closure on a smoothly continuing root the first-order
slope is unchanged; no smooth path is claimed across a root-selection change. A
failed donor density search is a failed finite check.
Retain each failed check and count energy/density evaluations. Admit numerical
descent only when the finite Gibbs decrease exceeds its propagated estimate.
Optimizer convergence and positive trial inertia are diagnostics, not descent
requirements. An inaccessible negative plane remains unresolved. In ordinary
negative-plane donor analysis, after reconstruction, support, domain and
positive-donor admission, a finite slope with a finite nonnegative complete
error estimate is labelled `no_resolved_descent` when $d_p\ge-E_d$; retain the
plane value, slope and error. This label means that this trial's slope does not
resolve descent against that estimate. Missing or invalid error terms remain
indeterminate. For $d_p<-E_d$, retain the negative-plane and conserved
finite-split requirements; failure to establish the finite decrease remains
indeterminate. The label neither asserts a nonnegative exact distance nor
certifies equilibrium.

The adopted output is **numerically resolved descending split with estimated
uncertainty**. Calling it a **proved instability witness** requires justified
bounds on feasibility and every error term, which native association floors do
not provide. A completed finite search reports only **no witness in the recorded
trial set**, naming unresolved trials. No global stability conclusion follows.
A finite negative plane with donor slope below minus its estimate but
unresolved finite descent makes discovery unsuccessful, keeping the candidate
for diagnosis; it is never reported as successful discovery.
Deduplicate physical trials using `PHASE_COINCIDENCE_DISTANCE` as a topology
resolution, not an energy-error threshold. Never discard actual resolved finite
descent solely because the trial is close to a present phase.

### Conserved topology changes and work limits

Rank distinct resolved descending trials by actual donor slope, more negative
first, ties by enumeration. This is deterministic, not a claim of the best
possible addition start. Seed the full added declaration with the admitted
finite split; preserve feed, reactions, standards and exact selected support.
Reuse shared density-branch admission. A unique unlabelled root takes the side
its stored negative-slope sample establishes (ALG-DENSITY, #423), keeping the
split and its Gibbs decrease: water(0.87)/1-pentanol at 150 K and 101,300 Pa with
a vapor trial family adds a second liquid. Without a sample it retains the source
declaration label without asserting a physical phase theorem; an ionic phase
established as vapor-side is refused. Require
achieved rows and an uncertainty-qualified lower Gibbs energy before accepting
an addition; otherwise try the next distinct trial. The acceptance margin sums
each present phase's molar energy estimate weighted by its amount in mol, so
the margin and extensive G/RT both have units mol. Retain this computed margin
with the addition trial; a per-mole estimate alone cannot qualify an extensive
Gibbs decrease.

The 100 mol feed-rescaling falsifier verifies this native extensive margin;
rejected cold/public declared reuse at that amount scale remains the separate
[declared-conditioning limit](../equilibrium/formulation.md#discovery-and-global-claims).
Do not infer a universal reuse guarantee or physical nonexistence from it.

Only `check_topology` produces removal events. Keep actual nonfatal finite,
domain-valid attempt coordinates internally before public `clear_state`.
Coalescence merges both component inventories and extensive volumes into the
first phase. Disappearance transfers its complete small inventory and volume
to the first remaining phase whose support contains it, ties by index. No
recipient means unresolved removal. Normalize merged amounts, reclose recipient
density on its declared branch (refusing a sole root established on the other side, #423), verify packing/association and balances/charges,
then solve the reduced declaration. Dropping a phase never drops its inventory.
The public unavailable-state meaning is unchanged.

For r=rank B, simultaneous phases cannot exceed r. The caller declaration-solve
budget defaults to **2r** and counts every initial invocation, seeded fallback,
failed addition/retry and reduced removal solve. Trial minimizations are counted
separately: for pair h with n_h neutral species and m_h removed-phase seeds,
there are at most **1+n_h+1_(neutral Vapor)+m_h** starts per round, summed over
admitted pairs before deduplication. Singleton handling can reduce this bound.
Each distinct start runs once per candidate; no hidden retry is permitted and
any shared-driver fallback counts. Caller trial-call limits further bound this
set. There are at most **2r+1 trial rounds**, including final testing; at most
**2r−I attempted additions** under the default budget, with I initial declaration
invocations and removals consuming the same remainder. A caller declaration
budget replaces 2r in that addition bound. Check limits before every driver
call. Total driver work is declaration plus trial invocations; retain roots,
reconstruction, finite-split evaluations, iterations/evaluations and Hessian calls.
Discovery option defaults belong to Equilibrium, not Core policy.

Partial enumeration is search-incomplete. Remaining resolved descent without
an admitted next declaration, exhausted budgets or descent at the balance-rank
phase limit makes discovery unsuccessful, retaining its candidate separately
for diagnosis. Rank-limit replacement and generalized multi-donor/reactive
withdrawal remain outside this slice. Completed enumeration without resolved
descent returns the declaration with the finite-search statement and unresolved
limits. Every failed/refused/indeterminate qualification follows ALG-TERMINATION's
fix-before-negative requirement; honest failure reporting alone is not completion.
Neutral qualification precedes ionic/reacting qualification. Same-wheel declared
comparisons establish integration/topology evidence; refined FeOs neutral
comparisons establish cross-implementation agreement, not physical validation.

## ALG-ACTIONS

An action request identifies observable, independent coordinate order, units,
directions, derivative order, source domain and branch/support. A successful
central solve does not imply every action is available. Return per-action
values or a specific failure reason, with rank/conditioning and linear residual
where relevant. Cache/reuse thermodynamic work for the same resolved state;
reuse must preserve all coefficient, T, density and implicit-site dependencies.

The bounded exact-inventory incipient extension in ALG-EQUILIBRIUM changes
central support admission only. Existing regular fixed-face actions retain
their support and vanished-inventory restrictions; insertion directions,
charged insertion properties and derivatives across a support change remain
outside that qualification. No new boundary enthalpy or derivative order follows
from an admitted pressure.

`CompiledSystem::lever_jacobian` assembles the first row-series coefficients
along unit T, P, feed-amount and active-parameter directions at one prepared
center. Its free T/P columns drive the corresponding state coordinate; the
row-series external control entry is unused for a free intensive. Reference
lever requests retain their existing refusal. The shared geometry includes
the ionic packing-diameter volume chain (EqID:
`packing_coefficient_parameter_derivative`), which the retired dense lever
assembly omitted. Chemical and pressure parameter columns at fixed packing
fraction therefore include that volume response. Pinned fractions and feed
scaling use the same row-series geometry as solved-state actions.

| Requested quantity | Generic energy derivatives required |
|---|---|
| pressure, chemical potential, residual h/u values | first |
| first derivatives of those quantities; fixed-composition cp/cv values | second |
| pressure Hessian on the independent T/rho chart (`pressure_state_hessian`) | third |
| cp/cv sensitivities | third |
| first equilibrium action (A1) | second |
| second equilibrium action (A2) | third |
| third equilibrium action (A3) | fourth |
| second total enthalpy action | typically third plus the equilibrium chain |

These are generic requirements, not an instruction to store dense tensors.
Directional contractions may suffice. Derivative identities apply only on a
smooth fixed branch/chart with regular defining Jacobians. DH arbitrary ionic
insertion Hessians at zero ionic strength do not satisfy that assumption.
Fixed-composition heat capacities and reacting equilibrium heat capacities
are different observables with different chains.

For the P4 native caloric action surface, each of Ideal, Residual and Total
`h`, `u`, `s`, `cp`, and `cv` admits mixed directional orders one through three
in `(T,rho,x;theta)` or `(T,P,x;theta)` on one accepted density branch.
Composition directions preserve normalization and charge balance; known fixed
parameter slots project to zero, while unknown slots are input errors. The
[accepted P4 design](https://github.com/tannerpolley/ePC-SAFT/issues/19)
specifies the qualification cases and numerical limits.

A direction entering a zero-fraction neutral component follows the exact
infinite-dilution limit (owner decision 2026-09-26, #73), with three exceptions
that report their cause and no action. The Ideal and Total entropy are singular
through `x ln x` (`support_change`). In a model with the dipolar term every result
keeps `support_change`: with every dipole absent the term's A2-only branch omits
A3 beyond second composition order (full EOS 4.1). Any direction entering a
charged component is `ionic_insertion` for every result, including the entropy,
when a neutral insertion accompanies it.

Polarization mixes each direction with the relative-`T`, relative-`V` and
unit-amount labels. Each direction is therefore divided by the power of two
nearest the largest of `|dT|/T`, `|drho|/rho`, `|dx_i|` and
`|dtheta|/max(|theta|,1)`, and each multilinear result is multiplied back exactly.
On the pressure chart the direction is scaled before its density lift, with its
pressure increment, and the lifted direction again, so an action equals the action
of the already-scaled request times the product of the scales. Unscaled, a one-pascal pressure direction in the
methane/ethane liquid at 180 K and 3 MPa (`dV/V ~ 3e-9`) misses
`sum_i x_i d ln(phi_i)/dP = 1/(rho R T) - 1/P` by 2.6e-9 relative; scaled, the
identity holds to rounding.

| Caloric family | Value energy order | First action | Second action | Third action |
|---|---:|---:|---:|---:|
| `h`, `u`, `s` | 1 | 2 | 3 | 4 |
| `cp`, `cv` | 2 | 3 | 4 | 5 |
| `ln(phi_i)` | 1 | 2 | 3 | 4 |

Component `ln(phi_i)` (EqID: `log_fugacity_coefficient`, #73) is a residual
pressure-chart request that returns one result per component, tagged by its
index; a density-chart request is an input error. Its jet is `d(A_res/RT)/dn_i`
at fixed `T`, total `V` and other amounts, the polarized path over the requested
directions plus the unit amount direction `e_i` at the unit-mole central state,
taken through the same implicit density chain as the caloric family. From it the
jet subtracts `ln Z` with `Z = P/(R T rho)` at the prescribed pressure, as
`state_at_pressure` does: the exact finite jet logarithm of the linear `P` and `T`
paths and the total density actions. In the dense methane/ethane liquid at
150 K, 0.1 MPa and `x_CH4 = 0.05` (`Z = 4.1e-3`), the density form
`ln(1 - V phi_V)` misses the `(T,P)` identity for `sum_i x_i D ln(phi_i)` by 8.3e-10
relative, against its 1e-10 criterion. With the prescribed-pressure form, every
`T` and `P` sum identity there and at 180 K, 3 MPa holds within 4.5e-14, except
the third-order `(T,P,P)` sum in the dense liquid: 5.2e-10 against its 1e-9
criterion; Gibbs-Duhem holds within 1.4e-12 of its term scale. Without
charged species the ion rows are `ionic_insertion` with no value,
as in the state. With two or three directions the batch also carries the mixed
density action `d^k rho`. Equilibrium solved-state requests do not admit this
observable.

Use the same residual-potential/mass-action tape and a single factorization of
the accepted unweighted site Jacobian for the batch. Obtain each implicit site
Taylor coefficient from the defining mass-action coefficient, then reuse each
directional path and its lower coefficients across all selected properties.
For partial extraction, auxiliary temperature labels use the fixed central
scale $T_0$: each label contributes $T_0$ K to the diagonal path. Divide the
extracted partial with $k$ temperature labels by $T_0^k$ ($k=0,1,2$).
Hold this divisor constant across the outer-direction coefficients; it is not
a temperature-dependent function to differentiate. Volume labels keep the
existing central-volume scale and $\rho_0^\ell$ restoration. A zero-label
extraction retains the central value. This changes arithmetic conditioning,
while preserving the physical increments and mixed derivatives.
The current-main density/scaling replay at `2428781b` revalidates this convention
through F6/F1, P4, ionic, Hessian/pressure/volume, F8/F9/Born and the frozen
binary128 residual-Cv check under unchanged criteria. The
[retained notebook](../../../analyses/2026-reactive-caloric-actions/notebook.qmd#integrated-density-closure-and-temperature-scaling)
records source/wheel identities, independent moved-row comparisons and the
bounded numerical-verification meaning.
For an E3 energy action with exactly two equal normalized directions, use the
sparse Taylor paths in EqID: `repeated_direction_energy_third` instead of linear
polarization. This avoids subtraction of a singular pure-composition third
derivative from a regular mixed volume response at dilute ionic references.
All-identical triples, other nonrepeated triples and E4 keep their established routes,
except the reference mixed slope below: its nonrepeated triple $D^3\psi[\nu,e_V,w]$
contains the dilute-ionic reaction direction, where the linear-polarization route
cancels catastrophically (aqueous NaCl ladder at $10^{-10},10^{-12},10^{-14}$ mol/kg:
3e-11, 6e-12, 2.5e-5 for $\partial^2C/\partial T\partial P$). It is evaluated by
balanced polarization of two repeated-direction actions,
$D^3\psi[\nu,ae_V\pm w,ae_V\pm w]$, whose difference is $4a\,D^3\psi[\nu,e_V,w]$
with $\nu$ inserted once, and $a=\max(V/T,|V_T|)$ sizes $ae_V$ to $w$ on the E3
scales (EqID: `standard_state_transfer_temperature_pressure_derivative`; the same
ladder then converges to 3.4895e-12 at all three points).
Pressure closure composes the full implicit density chain through the requested
order; a first density tangent alone is insufficient at orders two and three.
The pressure root search precedes this resolved-state action work.

Select dependencies before recording: density-chart Ideal-only work does not
consume the residual EOS, residual-only work does not consume thermal intervals,
and pressure-chart Ideal entropy retains the implicit density dependence.
At fixed density and composition, Ideal thermal fields have no explicit residual
parameter dependence. Total actions use the corresponding Ideal and Residual
actions, preserving the Ideal temperature and composition terms.

A pressure-chart request with no observables is a density-only request (#86): the
batch carries the resolved root and, per direction, the first-order derivative of that
density along the full control direction at fixed pressure closure, with a density
status that follows the support, association-residual, root and finiteness gates of
the caloric actions; the density chart has no such request.

Preserve independently available central values when an action is unavailable.
Distinguish missing thermal data, thermal interval boundaries, support changes,
singular or ill-conditioned site equations, excessive implicit residual,
unavailable density branches and nonfinite completed results. A zero or
nonfinite leading jet coefficient required for division, including the
pressure-volume derivative used by `cp`, remains a native `domain_error`; it is
outside per-result failure reporting. Branch and tape-predicate changes
are rejected by the finite-difference admission ladder before comparison; the
local action request has no perturbation size from which to report them as
runtime statuses.
Positive finite `P_rho` is required for pressure-chart actions and direct residual
or total `cp`; it is not required for direct `cv`. Malformed inputs remain
validation errors rather than scientific unavailability.

The conservative absorber's first exact residual Jacobian already needs A2
equilibrium and caloric actions. A limited-memory outer Hessian does not remove
that need. An exact outer Hessian may require A3. The older ordinary seven-state
absorber route only needing A1 does not define the new application's requirement.

Solved-state actions are evaluated by one native owner,
`solved_state_actions(compiled, central, request)`, at an accepted result; the
owner never re-solves, compiles or mutates. Admission requires the central
classified residual, the quantity the driver certifies (owner decision
2026-09-16), at or below `min(row_tolerance, EQUILIBRIUM_ROW_RESIDUAL)`; no
candidate-topology event. No composition cutoff applies (#322): in log coordinates a trace
species' column has a unit entry in its own stationarity row and O(x) entries elsewhere,
so rank, conditioning and actions do not depend on x (identical-ethane isomer actions
match the analytic values within 9e-14 for K = 1e-21 to 1e-100, condition 11.94 at every K).
Trace accuracy is bounded instead by the balance-row resolution of the state itself
(ALG-TERMINATION): a relative trace action carries the state's closure, 2.77e-6 for CO₂ at
MEA loading 1e-6 and 2.5e-3 at loading 1e-8 with row tolerance 1e-10. Admission also requires a finite square full-rank row Jacobian,
factorised once; and, for a thermal direction, a reaction target whose
temperature dependence is not frozen by a point-only supplied Gibbs vector
(`missing_thermochemistry`). A mixed action of order k is the half-sign
polarisation of k diagonal control paths q(t) = q0 + t w with the first sign
fixed positive: along each path the state y(t) = sum_j y_j t^j follows the
coefficient recurrence J y_j = -[t^j] F(y(t), q(t)) with the unknown y_j set to
zero, where the row owner composes the analytic geometry, ideal, balance, charge
and standard-state terms (Taylor in the driven or free temperature through order
three) with the model's residual gradient coefficients G_1 = H z_1,
G_2 = H z_2 + D3A[z_1, z_1, .]/2 and G_3 = H z_3 + D3A[z_1, z_2, .]
+ D4A[z_1, z_1, z_1, .]/6, the contractions coming from `residual_energy_action`
(E3 and E4 of the extensive residual Helmholtz energy, the sole higher energy
owner; it also contracts E5 for the critical rows' curvature). The central Jacobian is equilibrated once to A = R J C with exact
power-of-two diagonals (column maxima, then row maxima) and factorised; each
coefficient solves A z_k = R b_k, is refined once with the defect R b_k - A z_k
accumulated in double-double (FMA TwoProduct and TwoSum, the compensated dot
product Dot2 of Ogita, Rump and Oishi, 2005: as accurate as twice the working
precision, then rounded), and y_k = C z_k. The defect uses only double arithmetic
and the correctly rounded `fma`, so the x87 precision control cannot change it (it
assumes round-to-nearest without flush-to-zero, like all double arithmetic): AMPL's
ASL (`fpinit`, loaded by pynumero and ASL-hosted solvers) sets the x87 precision
control to 53 bits, which reduced the former `long double` defect to double and
cost the Lithium capacity-state action three digits (#173). Saving and restoring
the x87 control word instead would be x86- and glibc-specific. Each solved coefficient is
re-evaluated; its normwise backward error
eta_k = |R r_k| / (|A| |z_k| + |R b_k|) (infinity norms, Rigal–Gaches) must meet
`EQUILIBRIUM_ACTION_BACKWARD_ERROR` (`linear_solve_residual`). This certifies
only that y_k solves a nearby linear system; it catches a stale Jacobian or a
jet mismatch. The first-order forward bound
|C^-1 dy_k| / |C^-1 y_k| <= 2 cond_inf(A) eta_k (Higham, Thm 7.2), with
cond_inf(A) = |A| |A^-1| from the central factorization (the reported action
condition estimate is this equilibrated value, not the raw cond_inf(J) that
`inspect_sensitivity` reports), refuses a coefficient above
`EQUILIBRIUM_ACTION_FORWARD_ERROR` as `ill_conditioned` (owner decision
2026-09-23); the first refused coefficient names the status. Equilibration
removes diagonal unit scales such as the raw-Pa pressure column, leaving the
refusal to states that stay ill-conditioned after scaling. Refinement is not
credited by the certificate: eta is re-measured from the re-evaluated series. For k >= 2 the
bound holds given the computed lower coefficients: their errors enter b_k and
are not re-bounded. It is normwise in the equilibrated coordinates and does not
certify small components of y_k relative to themselves; it holds up to the
double-precision floor of the evaluated residual and excludes the final rounding
of each returned double. The native `e5-8/forward-bound` rows test it
against the exact bubble-point feed response of an identical m = 1 pair at
fixed T/(eps/k), the same reduced state at every energy scale: the raw
cond_inf(J) grows from 3.5e9 to 4.7e10 over the ladder while cond_inf(A) stays
22-32, every rung is admitted, and the equilibrated forward error stays within
2 cond eta + eps. Observables are
mapped along the same curved path: normalised composition, molar density 1/v,
packing fraction, N v, the phase pressure RT(1/v - a_v) and mu_res/RT from the
same gradient series; the log activity on the EOS density basis,
ln a_i = ln(x_i rho / (1 mol/m^3)) + mu_res_i/RT, is the series logarithm of
the normalised x_i rho plus that potential. It carries no standard-state term
and no ionic gauge, is rejected outside the phase support, and for an ion is
the model's unconstrained single-species value, physical only in
charge-neutral combinations. Fugacity is f_i = a_i R T (1 mol/m^3) and the full
chemical potential g_i(T)/RT + ln a_i needs the caller's thermochemical datum;
neither is a solved-state observable. Phase calorics compose the P4 density-chart actions of
each order with the lower-order curvature terms, [t^2] g = Dg[z_2]
+ D2g[z_1, z_1]/2 and [t^3] g = Dg[z_3] + D2g[z_1, z_2] + D3g[z_1, z_1, z_1]/6.
Precedence is input validation; central solve and classified rows; support,
chart and topology; finite, rank and linear residual; the selected P4
dependency; a completed nonfinite value. A zero direction returns an exact zero
after admission, and a fixed T or P keeps its identity action. Feed directions
must annihilate the preserved vanished-inventory restrictions. A direction that
changes them reports `support_change` or `ionic_insertion`, preserving central
values and independently available directions. Fixed-face EOS and reaction
columns use the unchanged full reference completion; removed amount/fraction
actions are literal zero. Continuation recompiles or reports a support event
when the inventory face changes.

`solved_state_first_actions(compiled, central, observables, directions)` evaluates
independent A1 requests through the same owner. Within that call only, it shares
the prepared center, central observable values, classified/topology admission
and defining-row Jacobian factorization. Direction-dependent reference admission,
thermal availability and propagation remain separate; a refused direction does
not invalidate siblings or erase their available central values and diagnostics.
Returned batches follow direction order. Their `state_factorizations` counts
record actual work: the batch first requiring the factorization reports one,
subsequent reuse zero. Nothing survives the call. The original request's direction
list still defines one mixed A1–A3 action, not independent directions.

Reaction controls supply affine `reaction_coefficients` increments of (a,b,c,d).
A scalar ln K uses its a increment. The native source recorder evaluates the
coefficient law at T(t), including all consumed thermal product terms; no
pre-evaluated reaction-log or local temperature-slope direction is retained.
Chemical rows use the minimum-norm stoichiometric image, EqID:
`reaction_constant_coordinate_derivative`. Caloric standards use the supplied
record anchors and the same record-less-species inverse as EqID:
`reference_species_calorics`. Supplied point Gibbs vectors cannot carry a
changing source correlation.

Referenced driven T/P caloric composition uses the anchored Gibbs derivative
columns through order four, including T/P order through three with one EOS
coordinate. Values and A1/A2 cover h, s, u_G, frozen c_p/c_v,G and v_G;
A2 permits a parameter/coefficient increment in at most one direction. Feed-only
A3 retains constant standards and frozen slopes at driven T/P. Free-intensive A1/A2 uses the same composition with its solved T/P coefficients.
The prepared standard refreshes at each current T and system-reference P,
including second row derivatives and requested T/P–parameter columns. The
existing point-local parameter cache contains both mixed columns consumed by
its free-intensive chart; driven intensive directions rebuild their selected
jet. A solved intensive cannot also have an external control. Free-intensive
F6 and apparent-Cp qualification are bounded by their separate retained #225 evidence.

For frozen c_v,G, Properties supplies the mechanical raw v_T and v_P from the
same pressure-temperature/pressure-volume jets used by its c_p calculation.
Equilibrium adds anchored g_TP/g_PP at fixed composition, forms the local
quotient in EqID: `gibbs_conjugate_molar_volume`, then follows the solved
composition path. The re-equilibrated total-volume derivative cannot enter
that quotient. Preserve mechanical volume and Gibbs-conjugate volume as
separate outputs. The local Gibbs-volume chart requires finite v_G and finite
nonzero v_G,P, in addition to the EOS/site gates. Record v_G>0 and v_G,P<0
as stable-chart diagnostics with the quotient sensitivities; their signs do not
restrict mathematical admission. Nonfinite or zero v_G,P remains PressureSingular;
independent h/s/u and volume quantities retain their own meaning when that
inverse is unavailable.

Actions use a dimensionless path variable and declared physical increments.
They retain the observable unit. Raw T/P or fitting partials divide by every
increment; Taylor factorial conversion is separate and applied once. The
[retained numerical notebook](../../../analyses/2026-reactive-caloric-actions/notebook.qmd)
owns the tested cases, tolerances, errors and interpretation limits.

## ALG-REFERENCE

Referenced solved-state third actions admit feed directions only at fixed T/P.
Thermal and mixed thermal/feed third actions retain the Engine refusal;
the solvent mole-weighting relation does not extend this derivative order.

Native fitting may supply solvent-reference and finite-activity density anchors.
Each anchor is checked for regularity before use; an anchor failure is recorded
and retried cold by the existing density owner. There is no reference-wide retry.

One ion-free solvent state per distinct source pressure supplies the analytic
reference of EqID `infinite_dilution_reference`. Pure neutral solvents and neutral
solvent mixtures are admitted for ionic activity references. A charged background,
even with zero net charge, is refused as IonicReferenceBackground. Reaction
reference declarations retain their single neutral solvent and source-basis rules.

Properties uses EqID `reference_single_insertion_taylor` for values and every
consumed derivative. The separate solvent-reference tape retains the shared Born
radial coefficient and regular contributions, with a zero DH insertion. It never
enters the ordinary-state tape cache or relaxes ordinary ionic-insertion admission.
The existing association Jacobian and implicit pressure coefficients supply all
T/P derivatives through order four and one parameter through T/P order three.
One amount insertion appears at degree r+1; extraction through degree 2r+1
cannot contain two such insertions. No ionic-concentration Hessian is required.

Requested source-reference errors for the dimensionless contraction $C$ and
its first thermal derivative $C_T$ use local higher-precision re-evaluation plus
first-order pressure/site propagation (#869). Named retained machine constants,
supplied parameters and captured path/site coefficients are fixed inputs;
computed reciprocals, products and sums are evaluated within the common scalar
expressions. The existing `long double` Taylor evaluator replays captured sites;
the selected-tape association-error owner propagates their coefficient errors
through degree five. It adds no second association or density solve for the replay.
The existing diagnostic cache owns the solvent and selected pure-reference states.
Only explicit error requests and detailed endpoint capture run this calculation.

For a returned coefficient, wider result $x_w$ and absolute extraction/assembly
scale $S_x$, the arithmetic estimate is

$$
e_x^{\mathrm{arith}}=|\widehat x-x_w|+\gamma_N^w S_x+N d_{\min}^w,
\qquad \gamma_N^w=\frac{N u_w}{1-Nu_w},\quad u_w=\epsilon_w/2.
$$

The recurrence count uses $4(d+1)^2$ operations per selected tape operation and
sign path, summed over the pressure, value and first-thermal extractions, plus
the short assembly. Differences and allowances are formed in wider arithmetic;
`nextafter` converts the final nonnegative estimate upward. The center is the
actual degree-one $C$ or normalized degree-three $C_T$ value. Pure/dilute sums
are recentered on their returned columns. With signed achieved pressure
coefficients, original target $P_s$, arithmetic/site floors $\eta_{P,k}$ and
the retained density floor converted once to Pa,

$$
b_0=|p_0-P_s|+\eta_{P,0}+\epsilon|P_s|+e_{\rho,\mathrm{Pa}},\qquad
b_1=|p_1|+\eta_{P,1},
$$
$$
e_C=e_{C,\mathrm{fixed}}+|C_{P_s}|b_0,\qquad
e_{C_T}=e_{C_T,\mathrm{fixed}}+|C_{TP_s}|b_0+|C_{P_s}|b_1/|T_1|.
$$

Raw source-pressure partials remain nonzero for a fixed source, whose public
system-pressure columns are zero. Existing species maps transport errors with
absolute coefficients, compensated equation/application defects and their
arithmetic floors. The anchored Gibbs-temperature column retains both $e_C$
and $T e_{C_T}$. Supplied caloric records receive no source-reference solve
error; higher caloric-error columns remain uncomputed. Present-solid affinity
uses its own source-reference error and the mapped species errors.

This is a local estimate, with the same comparison meaning as ALG-TERMINATION,
and can underestimate shared errors and internal cancellation. A final-scale
roundoff allowance does not bound arbitrary internal amplification. Fixed-input
arithmetic excludes supplied operand uncertainty: an exact quotient at both
precisions can still change under an uncertain denominator. The retained
underflow and reciprocal examples and the qualified stencil are in the
[forward-error notebook](../../../analyses/2026-equilibrium-forward-error/notebook.qmd).
Pressure/site propagation omits second-order sensitivity-error products and
requires the existing admitted regular branch. The result annotation states
these limits separately from failure causes and labels wider precision only as
the reference replay's precision. The endpoint/Jacobian precision and shared
`estimated_error` eligibility and judgment are unchanged. Where `long double`
is not wider than `double`, preserve ordinary values and return the explicit
cause `source-reference estimate unavailable: long double is not wider than double`.
Missing states, unavailable sensitivities and nonfinite estimates retain their
actual causes; error production never changes physical admission.

Values and columns have independent availability, actual association residuals,
site conditioning and pressure response diagnostics. Density convergence is
retained separately. No terminal sequence or terminal error estimate is measured;
those optional fields are unavailable. No failed column erases a value or sibling.
The frozen finite ladder remains qualification evidence only. Chemical rows and
record-anchored caloric completion consume the same returned columns and causes.
A model keeps its last 256 `reaction_references` results by exact inputs (temperature,
system pressure, declaration, reactions, solvent anchor, orders and jets); a repeated call,
such as each endpoint of the finite-addition rows at one (T, P), returns the stored
diagnostics unchanged (#380).

Pressure-closed fugacity evaluates ln Z from the prescribed P and resolved
density; the returned EOS pressure still records the density-root residual.
This avoids amplifying the closure residual in parameter differences of nearly
invariant reference energies. Pressure-closed columns use EqID: `log_fugacity_coefficient_parameter_derivative`.
Each reference column is the terminal-pair limit of EqID:
`reference_limit_terminal_pair`, admitted independently when s_theta
(`standard_state_transfer_parameter_derivative`) times its maximum two-ion
limit change is at most the existing reference tolerance. An unavailable
column preserves the value, molality and terminal value error. Activity and
single-ion energy actions combine these columns through EqID:
`electrolyte_activity_parameter_derivative`.

Use the explicitly declared solvent/reference path and source standard state.
The two-ion formula must contain charged components with zero stoichiometric
charge; the solvent composition must also be neutral within the composition
admission tolerance. The shared composition owner enforces this before either
reference construction or activity evaluation; neither calculation supplies an ionic gauge.
Evaluate the retained molality sequence; return the terminal change and
criterion separately. Unknown validity or failed limiting convergence remains
unavailable. Differentiate both the source correlation and EOS reference path,
including its density closure. Fixed source pressure P0 must not become trial
pressure. Reference-gauge choices that preserve reaction affinities do not
supply absolute caloric data.

For the single-solvent osmotic coefficient, reuse the pure-solvent liquid State
at the same prescribed T/P and EqID: `practical_osmotic_coefficient`. Replace
the independently terminated log-fugacity subtraction with the residual
chemical-potential difference and `log1p` density/composition ratios. Evaluate
that difference through the extensive Helmholtz Euler identity, using the
prescribed-pressure compressibility in its volume term; accumulate this scalar
contraction in extended precision before the 1/s division. The first parameter
columns remain the exact CppAD fugacity contraction. Keep the
existing density closure and ionic reference ladder. Below the provisional
`OSMOTIC_MOLALITY_FLOOR`, return `Stop::converged_to_floor` in the osmotic
quantity (precision unresolved), retaining ionic results; no admitted
single-solvent reference still means absence. The scalar remedy is applied
before this reporting restriction. Qualification in #214 is limited to the
retained aqueous chloride packets at 298.15 K and 1 bar; neither the floor nor
the observed trace error is a universal accuracy guarantee.

For an explicitly referenced neutral reaction, use
`neutral_reaction_reference_transfer`: one pure neutral solvent, charge-neutral
positive solute ratios and one declared source basis per reaction. Resolve one
ion-free solvent liquid root per distinct source pressure, requiring positive
pressure slope and regular association. The exact solvent state replaces the
former near-pure terminal state. Activity and reaction references share the same
insertion owner, with the density-standard conversion applied once.
An optional `species_reference` map assigns `InfiniteDilution`, `PureLiquid`
or `PureGas` independently. Omitted species retain solvent liquid and all
solute infinite-dilution references. Empty and full explicit default maps use
the original calculation directly, with no pure-state allocation or EOS call.
The solvent's default liquid selection reuses the exact solvent state.
For each effective pure override, resolve one neutral pure-component State per
source-pressure trace on the declared native liquid or vapor branch. A single
stable root is admitted only when the pure isotherm's ALG-DENSITY scan finds no
mechanical loop (above the EOS critical temperature it serves either declared
branch) or the root lies on the declared side of the loop. A failed or absent
root retains its species index and native cause. Retain independent availability of each requested mixed contraction and column. Every consumed derivative masks those pure species out of the
dilution insertion and adds their pure-state insertions through the existing
single-insertion Taylor owner, reusing the resolved States. Fixed-source
system-pressure columns remain analytic zeros. Nondefault mixed declarations
with `CommonMolalityInfiniteDilution` are refused with a source-basis cause;
raw mole fractions are the qualified mixed basis in issue #274.

Convert the correlation to the density target through the existing reaction/SVD
owner. Common-molality inputs
receive the solvent-mass term once; for raw mole-fraction inputs the two basis
terms cancel. Preserve each reaction's fixed or system-pressure reference.
At fixed temperature, evaluate system-pressure references at the current
nonlinear point and use EqID: `neutral_reaction_reference_pressure_jet` for the
exact first and second pressure derivatives. Reuse a reference path for
reactions sharing a source pressure. Apply the same current-point jet to row,
Jacobian, weighted Hessian and solved-state action evaluation, and refresh
returned diagnostics at the solved pressure. Fixed-source-pressure references
carry zero system-pressure derivatives; their compile-time values are reused at
every point of the compiled problem, since model, temperature, declaration and
source pressure are fixed there.
Driven T/P state and record-anchored caloric composition follow ALG-ACTIONS:
A1/A2 retain every consumed T/P/feed term and one parameter-bearing direction,
including source coefficients evaluated at the path temperature. Driven feed-only
A3 keeps the standards and frozen partials constant while composing the full
amount/composition products. The defining rows retain their minimum-norm
representative; the caloric Gibbs completion retains the supplied record anchors.
The two representatives preserve the same affinities but their T/P derivatives
are not interchangeable (EqID: `reference_species_calorics`).

The existing value and lower-slope admission remains: first thermal and EOS
slopes use the value criterion at their physical scales; the second thermal
slope uses 0.01 J/(mol K); the mixed T/EOS slope retains its declared
$Ts_\theta$ scaled criterion (EqID: `standard_state_transfer_parameter_derivative`). Higher consumed columns use ALG-REFERENCE's
three terminal molalities and both scaled consecutive changes. For frozen heat
capacities, complete both changes through the anchored Gibbs/volume quotient;
individual column admission is not a bound on an ill-conditioned quotient.

A missing column refuses only the outputs and action directions consuming it,
with its derivative, reaction, coordinate and original cause retained. Supplied
records must be consistent with each consumed reaction derivative. A source
motion must lie in the range of the record-less-species stoichiometric system;
check the combined direction, not each unit coordinate in isolation. Never
silently project an incompatible source motion into a least-squares caloric
standard. Owned interval endpoint values and isothermal actions remain valid;
only thermal actions crossing that endpoint are refused.

System-reference u_G and frozen c_v,G use the same anchored Gibbs volume,
EqID: `gibbs_conjugate_molar_volume`. Fixed-source system-pressure terms are
zero and recover mechanical u/c_v. Independently valid state, residual and
caloric outputs survive a refused Gibbs-volume chart or another missing
column. Referenced one-free-intensive pinned cells use current-point standards
and the same bounded A1/A2 chain, retaining solved-intensive motion even if a
first tangent is zero. The simultaneous-free chart remains unsupported for
its missing independent constraint. The separate F6 and apparent-Cp numerical comparisons own their qualification.

Reaction activities and $K(T)$ use one declared standard-state identity. Direct
activities are dimensionless caller inputs; mole-fraction activities use explicit
$x_i^\circ=1$; molarity and molality require explicit positive standards and the
matching density or solvent-mass inputs. Do not infer a basis from units. Ionic
contractions require either a charge-neutral reference contraction or an explicit
ionic gauge; that gauge cannot manufacture missing species calorics.

Evaluate an admitted natural-log correlation as
$a+b/T+c\ln(T/T_{ref})+dT$ with explicit positive $T_{ref}$, analytic derivative,
declared interval and standard-state identity. Convert legacy logarithm bases by
algebra only. When matching species thermochemistry is supplied, require both
the van 't Hoff enthalpy relation and heat-capacity/integral consistency; otherwise
report consistency unavailable. A mismatch is a reference inconsistency, never
a reason to alter a coefficient or anchor.

Extensive H is assembled by Equilibrium from its complete phase-property jets and
phase-amount Taylor products. Pinned-zero amount series contribute exactly zero.
`SystemEnthalpy` and `AdditionEnthalpy` in parameter fitting consume this owner.
`EquilibriumIsobaricHeatCapacity` uses EqID: `equilibrium_isobaric_heat_capacity`:
thermal H A1 supplies the value and H A2(T,w) supplies its first action, including
the feed-normalization derivative. The call shares the prepared center and
factorization across ordinary, thermal and mixed requests; nothing is reused
across calls. Driven T/P is required, and higher apparent-Cp actions are refused.
An unavailable outer action preserves an independently available Cp value.

## ALG-FIT

Fits use exact Engine derivatives under the
[native calculation rule](../../../AGENTS.md#native-calculation-rule).
Ceres remains the primary least-squares owner. An explicit `MinimaxRequest`
selects native constrained estimation; a failed least-squares check never selects
it automatically. The same native sequence owns declared starts and shared
iteration/elapsed budgets for both methods.

Use the [parameter regression formulation](../regression/science/pure-saturation-regression.md)
for observation algebra, weights, parameter coordinates and uncertainty.
Resolve each state once, evaluate its declared observables/actions, transform
residuals and Jacobians consistently, then apply the declared bounded fitting
method. Initialization, bounds, loss, covariance assumptions and budget are
scientific inputs to the fit. A failed state is a failed evaluation, never a
successful zero residual. Fitted values remain candidate parameters.

`FitResult.rank_records` publishes the existing SVD decision for the complete
weighted optimizer observation Jacobian, including admitted measured-input
columns and penalty rows. Its dimension is the column count; its absolute
cutoff is `FIT_RANK_RECIPROCAL_CONDITION` times the largest singular value.
This numerical determination is separate from statistical covariance:
unfinished, active-bound and rank-deficient candidates keep their actual rank
while withholding the normal inverse. No parameter forward-error estimate
is derived from covariance, and opaque Ceres QR or Ipopt BFGS/KKT condition
estimates remain unavailable.

`first_numerical_failure` retains a producer's first refusal through prediction,
actions, fit trials and worker replies, including successful cold replacement.
Within an evaluation, only the serial row prefix through its first refused
row is eligible; recovered causes are selected in row and producer order.
Across evaluations the fit retains the chronological first cause. A
`TrialFailure` keeps its own evaluation's cause, row and optimizer point.
Recovered causes on successful rows do not create rejected trials.
Measurement LLT refusal is infeasible, rather than necessarily rank deficient;
whitening distinguishes nonfinite solved operands from overflowing cost.
Hydration holds retain both caller criteria in detailed capture, the existing
scalar derivative cutoff and actual bracket/residual operands. Their Core
record has no invented single criterion for the two acceptance conditions.

Minimax retains `ipopt_diagnostics` only for an actual invocation. Its last
callback tuple, restoration entries and initialization/application returns
have the same observer as Gibbs minimization; a point reader creates no
optimizer attempt. Earlier concrete observation/hold causes precede generic
callback or application-status text. Final coordinates and multipliers keep
the existing constrained-estimator basis.

See the fitting appendix for [cation-mixing coordinates](#cation-mixing-coordinates),
[declared measurement uncertainty](#declared-measurement-uncertainty),
[optional fitting-density offers](#optional-fitting-density-offers) and
[guarded minimax](#guarded-minimax).
The implemented fit (R1--R3) includes pressure-closed density and total molar
isobaric heat capacity with fixed thermochemistry, driven-T/P apparent equilibrium
isobaric heat capacity, finite-addition enthalpy,
the declared equilibrium observations and direct source-basis ln K observations. It uses the cross-family affine
optimizer coordinates of EqID: `affine_parameter_coordinate_map`, a Ceres Levenberg--Marquardt trust
region with dense QR on one thread, the named `FIT_*` policy tolerances as the
only defaults, exact owner actions for every Jacobian and a covariance that is
withheld under rank deficiency or an active bound. A least-squares fit evaluates its state
groups in `FitControls.workers` processes, its own included (0, the default, means one per
hardware thread), at most one per non-electrolyte state group; minimax and the point readers
use one. Workers are exec-started copies of the Python host's interpreter on the extension
image the host loaded; each dies with its fit (a worker whose fit died before it armed that
signal exits unread), reads the request once, and the fit kills and reaps them when it
returns, fails or is interrupted. An active trace capture or a native host keeps one process.
Worker k owns every k-th non-electrolyte state group in row order for the whole fit, with that
share's warm starts and same-point centers; electrolyte rows, which share one reference cache,
stay in the fit's own process. The fit merges rows in observation order and reports the
one-process sweep's first failure, and every process publishes starts and its cache only when
the whole evaluation succeeds, so the worker count does not change the evaluation. Processes,
not threads: Ipopt 3.11.9 drives MUMPS 5.6, whose Fortran modules and C I/O state are
process-global. A Python KeyboardInterrupt stops the fit between rows only
when it runs on Python's main thread, which alone receives signals; elsewhere it runs to its end (#380). Named scalar ln K and stored
`ReactionLogPolynomial` a/b/c/d coordinates share the affine chart with EOS
coordinates. Every matching observation must declare the same reaction,
coefficients and source basis. Rebuild each trial with those values, preserving
reference transfers. Optimizer columns follow declaration order; physical rows
are independent EOS targets, held Born diameters, then reaction coordinates. Both direction families use one
solved-state first-action batch and one defining-Jacobian factorization per state.
The native point reader returns prediction columns in `optimizer_gradient` and
unit-weight residual columns in `residual_optimizer_gradient`. Difference scale a
gives dr/dv = (dy/dv)/a; natural-log LogRatio gives dr/dv = (dy/dv)/y.
The Ceres Jacobian also includes square-root row weights. Minimax consumes
prediction columns and applies its own native check-error slopes.
`reaction_enthalpy_shift` declares one source enthalpy increment in J/mol and
a required finite positive pivot in K within each supplied correlation interval.
Its form-free direction is the derived instance in
`reaction_constant_coordinate_derivative`. Every trial applies it once to its
fresh observation copy; zero shift preserves the supplied coefficients.
It starts at physical zero and admits plain coefficient coordinates on the same
reaction. Plain assignments precede shift additions, independently of declaration
order. Duplicate plain coefficients and two shifts on one reaction are refused.
The existing named and species-reference admission remains in point readers.
The fitted chart publishes the native law by reaction name, alongside the
existing physical shift row and any plain coordinate values. To replay, attach
the completed native laws and fitted EOS parameters without applying coordinates
again. The native `reaction_log10_reciprocal` helper
returns its decadic coefficients and refuses $c\ne0$.
Sequence starts and minimax explicitly refuse this coordinate; their existing
coefficient and EOS charts retain their behavior.
`LogEquilibriumConstant` uses the declared fixed positive temperature and named
source record directly (EqID: `source_log_equilibrium_constant_observation`).
It skips equilibrium compilation, reference transfer and solving; its exact
columns reuse the coefficient basis above, with zero EOS/other-reaction columns.
It admits Difference residuals only and shares named-coordinate admission.
`PhaseAmountRatio` composes component amounts by the quotient rule; any constant
caller factor belongs in the observed value and, for Difference, in its scale.
`AdditionEnthalpy` uses two declared equilibrium problems and the row's feed
state; its value and exact first actions follow EqID:
`finite_addition_enthalpy`. Both endpoint solves remain available through
residual-to-Jacobian promotion. Each endpoint uses its own solved-state action
batch and reaction-row mapping; the feed's ideal enthalpy is fixed with
thermochemistry and its residual enthalpy and EOS-coordinate action use the
Properties pressure chart. Every phase not pinned at zero contributes amount
times total enthalpy. The row forms a singleton state group, rebuilds at each
trial, and uses Difference for negative exothermic values; LogRatio retains
the existing positive-argument admission. Admission checks the material balance and identical ordered reaction
and source-reference declarations. Referenced free-temperature endpoints use
the shared native first caloric actions at each solved temperature. Failures retain the `initial:`, `final:` or `feed:` state name
with the owning solve, density, thermochemistry or action diagnostic.
Affinity targets, supplied Gibbs vectors, dependent reaction sets, constants
fixed by supplied thermochemistry are rejected by reaction name before fitting.

Admission refuses discovery options, since every fit problem declares its phases, and
validates every row against a model at the start parameters, keeping nothing:
the initial evaluation builds its trial model and compiles every state group at the resolved
optimizer point. A same-point Jacobian reuses that trial's model and compiled centers, and the
same compiled problem may serve a cold retry at that point. Failed trials still
publish neither solved starts nor a successful evaluation cache.
The two endpoint problems of an `AdditionEnthalpy` row are retained separately;
its singleton kind is always rebuilt at a new parameter point.

Electrolyte least squares admits molality-scale gamma, single-solvent Phi,
single-ion solvation and B-minus-A transfer energy (J/mol), and bulk permittivity
or a two-state ratio. Reference keys include T/P, solvent composition, ions and
stoichiometry; finite-state keys additionally include molality or full composition.
Projection metadata never selects the reference or finite state.
Dependent Born holds use EqID: `dependent_born_diameter_derivative`, safeguarded
Newton inside the declared bracket and both declared inner criteria. The same
point-local D maps all optimizer columns and conditional physical covariance.
Only positive-weight objective rows enter trial values, actions, rank and
variance. Zero-weight rows receive admission and final evaluation; their failed
final residuals retain independently available signed values and objective results.
Native final checks recompute the stored parameters without resolving a hold,
report signed row errors and max/mean groups, and use frozen source baselines
where declared. Least-squares objective results retain their accepted warm-start evaluation;
only zero-weight reports replace their final values. Unevaluated trial rows remain
explicitly unavailable. Report-only failures do not change acceptance. A failed check
is no proof of infeasibility. Singular vectors and the projected weighted
gradient describe local conditioning; dependent-target uncertainty, profiles,
parameter Hessians and physical validation are not supplied.

## ALG-TERMINATION

### Diagnose an observed refusal

The first numerical refusal and final convergence answer different questions.
A successful substitution, cold replacement or polish can retain a rejected
intermediate operation. Its presence does not establish the cause of a later
termination. Default records retain that first refusal, the latest rank per
named operation, the last actually invoked Ipopt attempt and the earliest
superseded restoration-bearing attempt. Detailed capture retains intervening
observations. Inspect `capture.trace.ipopt_log_error` and
`capture.trace.incomplete` for a failed log destination. Diagnostic retention
does not handle out-of-memory failure specially.

ASL replay uses the same real controls and requested derivative order. Its
density errors format the existing first cause, falling back to density
convergence; solved-state action errors retain rank/dimension and present
condition/backward-error operands. The condition basis is the power-of-two
column-then-row equilibrated defining Jacobian, $RJC$. Rank refusal precedes
the inverse-condition estimate and Taylor solves, so absent operands stay
absent. ASL never computes an estimate during conversion. Compare successful
native and ASL actions in the same $(T,P,w)$ coordinates, including the
normalization chain rule for $x=w/\sum w$; raw composition actions have a
different basis.

1. Retain the exact declaration, parameter/library identity, units, starts,
   requested criteria, work budgets and default result. Read final convergence
   separately from `first_numerical_failure`, `rank_records`, `last_rank`,
   `ipopt_diagnostics` and `recovered_ipopt_diagnostics`. Discovery declarations
   and trials retain their own observations; unrelated trials do not share a
   first cause.
   For a fit, inspect its rank records separately from covariance status and
   its first refusal separately from termination. For minimax, read
   `result.minimax.ipopt_diagnostics` and the existing point/multipliers.
2. With `EPCSAFT_IPOPT_OPTIONS_FILE` unset, replay those inputs through
   `with epcsaft.trace(diagnostics=True) as capture:`. Native callers use
   `core::Capture(true)`. Read the deepest association, density or callback
   refusal before attributing an outer termination to rank or restoration.
   Detailed estimates are separate observation work and can advance finite
   deadlines; untimed numerical equivalence does not promise deadline identity.
3. Request `equilibrium_numerical_estimates(model, result)` or compiled
   `numerical_estimates` at the actual admitted returned point. Retain the
   basis, cutoff, residual, floor, condition, availability and physical-quantity
   radii. A failed optimizer's point is labelled **Ipopt final point before
   polish**, not an admitted returned-state estimate. Its raw TNLP coordinates,
   lower/upper bound multipliers and constraint multipliers retain the
   originating declaration, rows, constraint-row selection, actual coordinate
   scale and objective factor. Engine coordinates are `coordinate_scale *
   tnlp_coordinates`; bound coefficients would divide by that scale. No
   transformation is reapplied to the stored arrays. Callback quantities use
   Ipopt's internal scaling basis and are not physical Gibbs-energy measures.
4. If Ipopt internals are needed, supply `ipopt_log="/absolute/unique/ipopt.log"`
   on that same capture; native callers use `core::Capture(true, destination)`.
   Keep the options-file variable unset. The destination opens lazily on the
   first actual Ipopt invocation and retains sequential warm/cold attempts
   without truncation. Nested captures inherit it unless they request another
   destination; the owning capture closes it on exit or `take`. Read every
   iteration, restoration marker, termination, factorization and line-search
   message alongside the bounded callback observations. A Newton-only solve
   opens no log. Inspect `capture.trace.ipopt_log_error` if opening fails; this
   failure leaves the calculation on its original path. J_MAIN uses level 5
   iteration summaries and other normal categories level 12. An all-category
   J_MAIN level-12 dump is not qualified on installed Ipopt 3.11.9.
5. Interpret each quantity in its own basis. Structural rank, endpoint rank,
   association pivot admission, pressure slope and Ipopt regularization describe
   different operations. Small backward error alone does not bound physical
   forward error; Ipopt's internal KKT condition is not the final defining-row
   condition. Preserve unavailable-estimate causes. For ASL retain `Errmsg`,
   stderr and real controls, then replay those exact inputs through these same
   steps. Preserve the original cause before later continuation failure;
   unavailable external functions supply no physical result.

The [options-file explanation](#ipopt-diagnostic-options) describes a deliberately
changed experiment. A nonempty `EPCSAFT_IPOPT_OPTIONS_FILE` disables warm reuse;
its presence, even with an empty value, bypasses discovery/free-amount Newton
predicates. Other eligible Newton paths can still finish without Ipopt. Preserve
that file/settings and unset the variable afterward. An incidental `ipopt.opt`
is ignored. Capture-owned logging changes no numerical option, eligibility,
retry or acceptance rule. File output and detailed observations can affect
finite deadlines; no timing identity is claimed.

Requested source-reference estimates follow [ALG-REFERENCE](#alg-reference).
Their returned method annotation supplies no new termination or admission criterion.

Core CORE1 keeps the numerical domains separate. `EQUILIBRIUM_COMPOSITION_BOX`
$\beta$ sets the lower box of the intensive `log(x_i)` coordinate inside the
optimizer subproblem: its barrier regularises the log chart's trace steps and
is not a solution domain. Species $i$ has the box in all-free cells

$$
x_i \ge \max\!\left(x_{\mathrm{cut}},\ \min\!\left(\beta,\ \varepsilon\,\frac{u_i}{N_{\mathrm{feed}}}\right)\right),
\qquad
u_i = \max\,\{\, n_i : B n = b,\ n \ge 0,\ n_j = 0 \text{ off every support} \,\},
$$

where $x_{\mathrm{cut}}$ is `EQUILIBRIUM_BOX_FLOOR`, $\varepsilon$ is
`EQUILIBRIUM_BOX_INVENTORY_FRACTION` and $u_i$ is the largest attainable total amount of $i$,
one HiGHS linear program per species, cached during compilation. With $B = I$ and every fed
fraction $\ge$ 1e-6, the feasible set is the feed point and $u_i = n_i$ is taken directly. Rows of $B$ can carry mixed signs (a hydrogen row
$n_{\mathrm{H^+}}-n_{\mathrm{LiA}}-n_{\mathrm{NaA}}$), so no single row ratio bounds $n_i$; an
unbounded or failed program keeps $\beta$. A uniform box $\beta$ is a solution domain whenever
$\sum_p \beta n_p > u_i$: it forces more of $i$ into the phases than exists and the subproblem
is infeasible, and just above that the optimum sits on the box beyond the polish's reach
(#172: 61 mol aqueous $\times$ 1e-8 = 6.1e-7 mol against 3e-7 mol Li fed; stalls up to 3e-6 mol).

`inventory_bounds(compiled, species, fraction, floor)` shares this HiGHS owner
with the optimizer box. The amount program uses $n/N_{\mathrm{feed}}$ and returns
its maximum in mol plus the separate dimensionless feed-normalized amount.
Fraction mode uses the exact transformed equations in the equilibrium formulation.
Its internal LP variable is $q=s\tau$, where $s$ is the smallest positive
$|b_k/N_{\mathrm{feed}}|$ capped at one (one when all balances vanish). This
column scaling preserves tiny nonzero inventory coefficients that HiGHS would
otherwise drop; returned witnesses and residuals use the original $\tau$ chart.
LP assembly errors or nonfinite reconstruction remain numerical failures.
It first maximizes $\tau$ to establish a finite floored inventory: a positive
feasible point, or a feasible point plus a positive-$\tau$ unbounded ray, supplies
`finite_witness`. A valid optimal $\tau=0$ means that finite set is empty.
Unresolved positivity or residuals above `EQUILIBRIUM_ROW_RESIDUAL` yield an
unavailable numerical result, with the LP cause retained. The species objective
is then restored; an optimum at $\tau=0$ is a limiting supremum only when finite
feasibility was established. `limiting` describes that maximizing witness; it
alone does not exclude a different finite maximizer (for example, $b=0$).
`witness` holds $n/N_{\mathrm{feed}}$ in amount mode
or $(y,\tau)$ in fraction mode. `amounts_mol` reconstructs the maximizing finite
point and is absent for a limiting supremum. Logarithms are available only for
positive fraction maxima. Balance, normalization and lower-bound residuals are
dimensionless. The finite witness has its own residual and termination cause;
`lp_error` retains HiGHS' reported primal infeasibility for the species program,
or for finite-feasibility termination when no species program ran.
Optimal, structural zero, infeasible inventory, unbounded amount, exhausted work
and numerical failure remain distinct. An unavailable maximum never becomes an
optimizer-box value. Structural-zero results carry zero without a logarithm.
The optimizer still keeps $\beta$ on failed/unbounded amount programs,
preserving the prior box formula and trace-inventory checks. It requests only
the compiled support union, so species outside every support add no LP run.

*Derivation of $\varepsilon$.* The boxes can hold at most $\sum_p x_{i,\mathrm{box}} n_p \le
\varepsilon (N/N_{\mathrm{feed}})\, u_i$ of species $i$, a fraction $\varepsilon$ of its attainable
inventory when the total amount $N$ is near $N_{\mathrm{feed}}$. Averaging the $S$ per-species
maximizers gives a point interior to every box jointly while $S\varepsilon N < N_{\mathrm{feed}}$
(fewer than about 100 support species at $\varepsilon$ = 1e-2). The polish must then carry the trace
composition from the box to the solution. The retained trace-Li probe bounds that reach: with the
uniform box holding 10 % of the Li inventory (6e-6 mol/L) the polish met the requested
tolerance, at 20 % it reached only the floor and at 30–60 % it stalled. $\varepsilon$ = 1e-2 sits
one decade inside the demonstrated 10 %. The minimum keeps $\beta$ wherever the uniform box
already displaces at most $\varepsilon$ of the inventory ($u_i \ge \beta N_{\mathrm{feed}}/\varepsilon$,
inventory fraction $\ge$ 1e-6), so the box is unchanged wherever $u_i/N_{\mathrm{feed}} \ge \beta/\varepsilon$;
below that it is lowered, including states the uniform box already solved (Li at 6e-6 mol/L). Lowering the box of
every minor species instead ($\varepsilon = \beta$, box $\beta u_i/N_{\mathrm{feed}}$) removed the
barrier's regularisation for a carrier at inventory fraction 5e-4 (box $\ln x$ from $-18$ to
$-26$): Lithium's four-step continuation corrector then wandered for 80 iterations and ended
locally infeasible at a Na-loaded state that the uniform box solves in 9. A globally
certified zero-inventory species is removed before box construction; a failed or
unbounded LP cannot certify absence. The box does not set accuracy: a
moiety row $(Bn-b)/N_{\mathrm{feed}}$ resolves an inventory to its achieved residual, at most the
caller's row tolerance and after the polish near rounding, $\epsilon N_{\mathrm{feed}}$ (trace
isomers fed at 1e-10, 1e-12 and 1e-14 mol in 1 mol return their amount to 4e-8, 5e-9 and 9e-3
relative; the orthonormal basis gave 8e-6, 3e-4 and a factor 4). Feasibility cells
use the smallest normal double for every species, so chemical equalities can reach trace roots.
$x_{\mathrm{cut}}$ is only the finite lower bound of the all-free box. It binds below an
attainable inventory of 1e-18 $N_{\mathrm{feed}}$, where trace values are already unresolved
(below about 1e-15 $N_{\mathrm{feed}}$); lowering it changes those values, so it stays. The
Newton polish carries a trace solution below the box to any positive root (#87,
2026-09-21; #266, 2026-10-07). No composition cutoff applies to actions, continuation or
declared seeds (#322): continuation from $K$ = 1e-16 to 9.3e-28 follows $K/(1+K)$ within 4.2e-14,
19 states below 1e-20. A declared seed needs positive finite compositions; its log
coordinates have no upper test, since the solver's own majority $\ln x$ returns up to 2.6e-13
above zero within the normalization residual. Grossly invalid declared log coordinates are
projected onto the optimizer bounds rather than refused; that this happens before any evaluation
is inferred, not tested. The lower limits are representability (a positive double) and the
solve's own acceptance. As probed: the single-phase isomer's stability decision is indeterminate by
$K$ = 1e-199 (minimum curvature 2.1e-203), and the square reactive bubble cell matches its scaled
$K$ = 1e-12 values to $K$ = 1e-150 but stops `infeasible` at 1e-200 and 1e-250 (undiagnosed). Free-phase
disappearance compares the dimensionless amount ratio
`n_phase / sum(n_feed)` with the provisional `PHASE_DISAPPEARANCE_AMOUNT`;
that topology rejection is not a stability or accuracy conclusion. Packing
admission rejects the equal-sphere face-centred-cubic geometric boundary
`eta >= CLOSE_PACKING_FRACTION`; this is a conservative PC-SAFT domain
constraint, not a mathematical singularity claim.

Numerical policy values stay in `native/core/numerical_policy.hpp`; source and
application accuracy criteria stay beside their accepted reference cases. Every
policy entry is listed in its `ENTRIES` registry, has an Engine consumer (a
lexical architecture check, not a semantic proof) and is readable as `epcsaft.tolerances()`;
changing a value means editing its entry and justification and rebuilding.
A verification check whose error exceeds its criterion but not an attainable
floor derived from reported quantities (condition estimate times rounding,
finite-difference noise over step, propagated inner residual) is *tolerance
gated*: a diagnosis recorded with the floor, its inputs and the ratio, not an
acceptable end state. The owning issue attempts the root-cause fix that lowers
the floor below the criterion (Jacobian scaling or equilibration, higher working
precision, iterative refinement, a better-resolved reference); the check stays
gated only while that fix is open or after retained evidence shows it
impossible. The floor is never a hand-picked looser tolerance, and an error
above it fails. A criterion for a near-zero multilinear action scales with
the sum of the magnitudes of its fixed chart-coordinate sub-direction actions
(temperature, pressure or density, composition), and a reaction affinity's with
the sum of its stoichiometric terms, and the reduced-Hessian stability band
with the split of the Lagrangian Hessian into objective and weighted constraint
curvature (ALG-EQUILIBRIUM), not with the near-zero result; no other split may
set the scale.
The solve counterpart is `converged_to_floor`.
Keep the achieved error, criterion, scale, iteration/evaluation budget, native
solver status and termination reason visible. Only Ipopt `MAXITER_EXCEEDED` and
`CPUTIME_EXCEEDED` count as equilibrium work-budget exhaustion;
a non-success optimizer termination with a finite non-converged residual is
unavailable, while a finite criterion-meeting result remains eligible after the
same bounded polish and physical checks. A non-success termination alone never
establishes success; floor-only results from it remain unavailable. Recomputed residual and topology criteria apply.
Distinguish converged, stagnated, exhausted,
nonfinite, left-domain, singular and infeasible outcomes. Optional traces aid
diagnosis and do not determine success. Density root producers classify all
optional diagnostics before publishing a criterion-meeting root; nonfinite
floor, contraction or amplification is unavailable and throws the existing
convergence failure.

Fix before a negative disposition (owner decision 2026-09-23): no analysis,
qualification or issue closes on a result that is refused, failed, unavailable,
reference unresolved, tolerance gated or a supported negative until its cause
is diagnosed and a root-cause fix has been attempted and shown impossible, or
the state is shown outside the stated physical domain, with evidence.
Conditioning, precision, solver and reference-resolution failures are defects
to fix. A typed refusal is a final state only for genuinely singular or
undefined mathematics, such as zero-support insertion, a Jacobian that stays
singular after equilibration or a reduced Hessian whose equilibrated minimum
stays inside its rounding band, and it states why. A state outside the stated
physical domain (e.g. packing at or above close packing) may also close,
stating the boundary. A model that still misses data after diagnosis and
attempted fixes is a physical finding that carries both. No value is ever
invented to avoid a refusal. Reviews ask whether the outcome could be fixed,
not only whether its label is honest.

Propagate inner-solve/input error into an outer accuracy estimate using the
applicable local conditioning and scale. A first-order estimate or pivot ratio
is not a rigorous error bound or a general condition number. Do not claim more
accuracy than supplied inputs support, and do not loosen accepted thresholds
to force convergence. Source rounding, numerical error, experimental uncertainty
and model discrepancy are distinct. Missing accuracy/runtime thresholds stay
owner decisions; historical microbenchmarks are not universal limits.
`Stop::nonfinite` is unavailable and has precedence over finite convergence and
floor predicates whenever the residual or an optional diagnostic is non-finite.
For equilibrium, `success` requires optimizer success or classified rows meeting
the caller criterion after a non-error stop (never after a failed initialization or fatal return), finite residuals and declared topology; an unsuccessful result retains `solver_status`, residuals,
the convergence record and failure message while clearing candidate state and
action arrays, so `has_state()` is the availability predicate. A successful
result retains the same convergence record alongside its state. The pre-solve
`NOT_RUN` record is a structural diagnostic with an unmeasured (NaN) residual
and absent criterion; it is neither a solver residual nor a physical
nonexistence claim. `requested_tolerance_met` is absent for that record.

## Verification obligations

Before a scientific behavior change is accepted, update or explicitly revalidate
the affected equation/algorithm IDs and execute the smallest independent
falsifier identified in the owning issue. Important examples are caloric
identities, association support boundaries, full weighted-Hessian directional
checks with nonzero multipliers, implicit action ladders and independent
reference-state transformations. Finite checks establish bounded evidence,
not universal symbolic correctness. Documentation readiness checks cannot
replace them. Expensive application/reproduction sweeps remain explicit.

## Appendix: Special cases and numerical qualification

### Ideal mole-fraction and molarity specifications

The ideal liquid `IdealMoleFraction` specialization uses the log-amount chart
and square direct reaction/inventory rows of EqID: `ideal_reaction_amount_coordinates`.
The molarity basis selects EqID: `ideal_molarity_reaction`, with constant reaction
Jacobian and zero reaction Hessians; exponential balance derivatives are unchanged.
Its cold start is positive on the admitted inventory support, or uses a declared
positive composition guess. The same Ipopt lifecycle and exact-Jacobian polish
classify these rows. No EOS density seed, branch admission, neutral recovery or
pre-Ipopt Newton attempt applies. Previous-result starts are refused before native
compilation. EqID: `effective_water_amounts` is a deterministic native amount map
from an available criterion-meeting result, with no further equilibrium iteration.

### Equilibrium initialization by declaration

<a id="equilibrium-initialization-by-declaration"></a>

#### Supplied-guess cap

For a consumed original declared composition guess, first apply the existing supported
positivity handling and normalization. Select species with
$0<U_i<$ `EQUILIBRIUM_TRACE_MOIETY_FRACTION`, where $U_i$ is the compiled feed
fraction when $B=I$, otherwise the positive compiled attainable amount divided by
feed total. A zero compiled bound stays uncapped. Let $R$ be the sum removed by
$\min(g_i,U_i)$ over selected species and $M$ the original normalized weight
outside that set. Only if $R>0$ and $M>0$, cap selected entries and add $R g_i/M$
to each recipient, then log the completed fractions without another normalization
(#652, R-a). Otherwise preserve the original seeding arithmetic, including when
no recipient exists. Calculated speciation, discovery and stability starts,
supplied multistart compositions, coordinates and previous states retain their
provenance and are excluded. A partially guessed fallback that discards the
guesses is excluded too. This changes only initialization: phase enrichment,
the optimizer box, stopping criteria and #407 phase amounts are unchanged.

#### Feed fugacity and support starts

For an unguessed liquid/vapor start, a sole unlabelled feed-density root cannot
supply both fugacity estimates: their self-ratio would make every $K_i=1$
and seed identical phases. In that case only the vapor fugacity estimate uses
the ideal-gas limit $\ln\phi_i^V=0$ (EqID: `log_fugacity_coefficient`,
$\mu_i^r\to0$, $Z\to1$). The resulting compositions still receive actual EOS
density closure on their declared branches; this estimate changes neither
the equations nor final admission. Two-root estimates and caller-guessed
compositions retain their existing rules. The [cold-vapor diagnosis](../../../analyses/2026-cold-vapor-branch/README.md)
retains the sole dense feed root and the recovered incipient vapor (#204).
Provisional $z/K,zK$ estimates consume fugacity coefficients only where the feed
weight $z_i$ is nonzero; an exact-zero weight contributes zero without requesting
an insertion property. This does not remove species from reactive supports (#265).
A nonreacting unguessed boundary of one liquid and one vapor, with the finite phase
F free and the incipient phase I pinned at exactly zero amount, starts F at $x_F=z$
instead of the normalized $z/K$ (liquid) or $zK$ (vapor) when
$D=\{i:z_i>0\}\setminus S_I\ne\varnothing$, with $S_I$ the compiled effective support of
I, or when F is the liquid and
$\mathcal L=[n_s>1]\lor[\rho_-\text{ available}\land\rho_L>\rho_-]$. Here $n_s$, $\rho_L$
and $\rho_-$ are the stable-root count, the largest stable root and the stored
negative-slope density of the unanchored feed-density search above at $(T_0,P_0,z)$;
$\mathcal L$ adds no EOS evaluation. At the solution the inventory rows give
$N_F=N_{feed}$ and $x_F=z$; using $z$ before the solve is a restricted initialization
policy, not that identity, and feed composition is not claimed to be the best start
elsewhere. Under ALG-DENSITY's single-loop assumption $\mathcal L$ is evidence that a
mechanically stable feed root lies on the liquid side of the scanned loop; it does not
show that the equilibrium phase is liquid. Near-critical methane/ethane feeds with one
unlabelled root and no witness keep $z/K$, which retains their solutions. The incipient
estimate, the pressure step and the density closures are unchanged (#391, #410).
Default-start results: [Cameretti 2005](../../../analyses/reproduction/2005-cameretti/README.md) and
[Gross 2002](../../../analyses/reproduction/2002-gross/README.md); `tools/check-equilibrium-corpora` pins them.
Some water/1-pentanol dews have two dew branches, one of them water-rich, so no single
liquid start is a general answer; the explicit [cooling onset](#cooling-onset-dew) selects
the intended branch.

#### Neutral Rachford–Rice start

For a neutral, non-associating, nonreactive, full-support two-phase fixed-TP request with both
amounts free, absent complete caller composition guesses, the existing feed-state
fugacity estimates define $K_i=\phi_i^L/\phi_i^V$. When finite positive $K_i$
bracket an interior root, bisect
$F(\beta)=\sum_i z_i(K_i-1)/(1-\beta+\beta K_i)$ on $[0,1]$ through the
binary precision of a double. Seed $x_i=z_i/(1-\beta+\beta K_i)$,
$y_i=K_i x_i$, $N_V=\beta N_{feed}$ and $N_L=(1-\beta)N_{feed}$, assigning
by declared phase kind. Nonfinite/nonpositive coordinates or fractions at the
existing phase-disappearance floor retain the original normalized $z/K,zK$
and default amounts below. No new phase search or equilibrium driver is introduced.

#### Reacting true-species split

A reacting fixed-TP declaration of one liquid and one vapor, identified by kind in
either order, with both amounts free and no composition guess, starts instead from a
true-species split (#350). Its predicates use the compiled effective supports, which
keep zero-feed reaction products and drop absent and unfed nonreacting species, and
compilation refuses a fed or reacting species on no support: the vapor support
lies inside the liquid support, so every potentially present species is on the
liquid support. The single-liquid speciation below, evaluated on the
two-phase compiled system's liquid chemical rows and keeping its last successful
sweep, gives the true-species fractions $z^*$ and total $N^*$. Then
$K_i=\phi_i^L(T,P,z^*)/\phi_i^V(T,P,y)$ on the vapor support and $K_i=0$ on
liquid-confined species, whose inventory makes $F(1)=-\infty$. Here $y$ only
estimates $\phi^V$, by two substitution passes $y\propto z^*\phi^L/\phi^V$ from
$\phi^V=1$; the split's vapor is $K_ix_i$. The same Rachford–Rice bisection
(Rachford and Rice, 1952; one owner for both splits) gives $\beta$,
$x_i=z_i^*/(1-\beta+\beta K_i)$, $y_i=K_ix_i$, $N_L=N^*(1-\beta)$ and $N_V=N^*\beta$.
$N^*z^*$ is conserved, so the seed carries the speciation's charge, and $Bn=b$ holds
only to the speciation's accuracy (at most 5.6e-6 mol over the qualification family):
this is an initializer with frozen reaction inventory, not a reaction and phase
equilibrium. The second $\phi^V$ pass and the split vapor's density root continue
the previous vapor root (ALG-DENSITY anchor); when the new composition moves the
root out of the anchor's window, the unanchored scan solves it. Outside this domain,
or after a failed speciation, nonfinite or nonpositive $K_i$, no bracket, a phase
below `PHASE_DISAPPEARANCE_AMOUNT` of the feed or a nonpositive fraction, the
apparent-species start above follows before the one solve; an exhausted budget
propagates and nothing retries after the solve. Supplied compositions skip the
split; supplied amounts and packings replace its coordinates. Without the [amount step rule](#equilibrium-amount-step-retention) the split loses accepted MEA states, each ending on a vanished phase; the
qualification is the [MEA flash family](../../../analyses/2026-mea-flash-family-seed/README.md) (#350, #356; PR #370).

#### Superancillary start

Only this non-associating neutral TP domain may use a bounded caller packing seed or skip
redundant density initialization for a previous successful result with compatible
positive coordinates. For a single neutral, non-associating, nonpolar component, the Bell–Deiters
PC-SAFT superancillary (2023, doi:10.1021/acs.iecr.2c02916; Zotero `BEC7SFNS`)
precedes the pressure loop below. It applies only in the same fixed-temperature,
free-pressure cell domain, with no charged-molecule declaration or segment-diameter
temperature law, $1\le m\le64$, $\tilde T\ge\tilde T_{min}$ and
$\Theta=(\tilde T-\tilde T_{min})/(\tilde T_c-\tilde T_{min})\le0.999$.
The upper start-domain limit is `SUPERANCILLARY_MAX_THETA`, beside the evaluator.
The retained NIST coefficients define piecewise degree-16 Chebyshev expansions in
$\Theta$ at 17 Chebyshev–Lobatto nodes in $w=1/m$; Clenshaw evaluation and the
barycentric interpolant on those nodes give $\tilde\rho_L$ and $\tilde\rho_V$.
Here $\tilde T=T/(\varepsilon/k)$ and
$\rho_\phi=\tilde\rho_\phi/(N_A\sigma^3)$, with molecular number density.
The start uses $P(T,\rho_V)$ from the shared state calculation, places both
densities by phase kind, and skips the initial density solves. The vapor pressure
is used because liquid stiffness amplifies density rounding at low temperature.
The unchanged equilibrium rows and admission still decide the result. Outside
that parameter or temperature domain, the pressure loop and original density
and coexistence checks below remain in effect. Meson compiles the losslessly
converted `engine/data/superancillary/pcsaft.json` into hexadecimal arrays; a
missing file or invalid table shape fails the build, and no runtime reads a file.
The source commit, file hashes, converter and NIST notice travel with the data.

#### Pressure-start loop

A neutral, nonreacting, non-associating liquid/vapor cell at fixed temperature and free
pressure, with full positive feed support and exactly one zero-pinned phase, first
attempts a Newton start in $\ln P$. It shares the driver's Newton and recovery
eligibility predicates; mixtures qualify only without a declared pressure seed,
so a seed selecting a retrograde dew branch keeps its original start. Pure fluids
qualify with or without a declared pressure seed. The finite phase starts at feed
composition $z$ and the incipient phase at its declared guess, else $z$.
Each pass solves the cold liquid and vapor density roots at its pressure and
compositions, reusing the first pass's feed states when already available.
With $K_i=\phi_i^L/\phi_i^V$, it takes $S=\sum_i z_i K_i$ for a bubble point
or $S=\sum_i z_i/K_i$ for a dew point, updates the incipient composition to
$z_i K_i/S$ or $z_i/(K_i S)$ and takes
$\Delta\ln P=\pm\ln S/(Z_V-Z_L)$ (plus for bubble, minus for dew).
This is the exact pure-fluid pressure slope; for a mixture it replaces the
feed-phase partial-molar-volume sum by its compressibility factor and is an
approximate Newton start. A sole liquid-side root uses the ideal vapor limit
($\ln\phi^V=0$, $Z_V=1$); such a pass cannot accept the start. A missing liquid
root, nonpositive $Z_V-Z_L$, failed state or nonfinite pressure rejects the loop.
Only a two-root pass with both $|\Delta\ln P|$ and the largest composition
change below `EQUILIBRIUM_START_PRESSURE_STEP` replaces the start. The cap is
`EQUILIBRIUM_START_PRESSURE_PASSES`; neither constant is an accuracy bound.
On rejection or outside that domain, the original start runs unchanged: a
multicomponent bubble declaration takes
$P_1=P_0\sum_{i\in V}x_i\phi_i^L/\phi_i^V$ at the seeded compositions.
Only when every phase is caller-guessed does it replace the vapor composition
with $y_i=x_i\phi_i^L/\phi_i^V\,P_0/P_1$; a failed state keeps $P_0$ and the
guesses. Pure fluids keep the existing mechanical-loop witness when positive.
The shared rows, Newton polish, Ipopt and admission still decide the result.

#### Reacting boundary and speciation

A cold reactive declaration at fixed temperature and free pressure with a phase pinned at
zero amount first solves its finite phases alone, at the caller's pressure seed and
with the same reactions and references; at the boundary those phases are the
solution's phases up to their pressure dependence, so their speciated compositions
replace the caller's guesses and an unguessed incipient phase takes the first
finite phase's composition on its support; the start then follows the
all-guessed rule above (#152: an unspeciated MEA liquid otherwise sends the
fugacity step to the wrong liquid's bubble pressure and the solve into feasibility
restoration). A caller start skips this seed. If the sub-solve, the seeding, the
seeded solve or its density-branch admission fails or throws, the caller's own
start runs unchanged; only the final attempt's diagnostics are returned. The
reported iterations, evaluations and Hessian calls add the work of the sub-solves
and of a failed seeded attempt, including the caloric seed's sub-solve, to the final solve's (#350).
Work that ends in an exception is not counted: a seeded attempt or a speciation sub-solve
that throws. Nor are the seed solves of the starts inside continuation, the stability
search and the no-solution fallback. A reacting
declaration of one free-amount phase, including that finite-phase sub-solve, starts
from the ideal-solution speciation of the feed balances $Bn=b$ unless a previous
successful result supplies its composition. Each chemical row's residual potential, $-\ln v$
and standard term form $c_j$, frozen at the seeded composition (the caller's guess or
the feed, absent species raised to `EQUILIBRIUM_COMPOSITION_BOX`) and its declared
density root at the seed $T$ and $P$.
The start minimizes $\sum_j n_j\,(c_j+\ln(n_j/n))$ subject to $Bn=b$ by the
Gordon–McBride element-potential Newton iteration and step control (NASA RP-1311,
1994, eqs 2.24–2.26 and 3.1–3.3; the trace-species limit applies only to a species
growing faster than the total) until every log-amount step, trace species included,
is below `EQUILIBRIUM_START_SPECIATION_STEP`; two sweeps refresh $c$ and the root at
the new composition (continuing the previous root, solved anew when speciation moves
it beyond the anchored window); the free amount starts at the speciated total and the
packing at the last sweep's root. At
an equilibrium $x^*$ the frozen terms make $x^*$ its own speciation, so a speciated
guess is kept. A failed or thrown sweep keeps the last speciation, or
the seeded start if none succeeded.
Without it, an unspeciated start (absent species at $1/|support|$, water ions at a
caller's zero) leaves under-populated species of the log-composition chart without
restoring force: their objective gradient scales with $x_i$, the Lagrangian Hessian
is indefinite there, and Ipopt cycled through restoration to the iteration budget
(#154, #157; [evidence](../../../analyses/2026-reactive-liquid-cold-start/README.md)).

#### Pure-fluid witness

Other declarations retain their original initialization, including the pure-fluid
coexistence check. A pure fluid's unguessed fixed-temperature, free-pressure
liquid/vapor start whose seed pressure $P_0$ gives one stable density root (#391)
instead takes $P_w=P(T,\rho_w)$ at the cold scan's negative-slope node $\rho_w$
(ALG-DENSITY), when $P_w$ is finite and positive, for its pressure coordinate and
both density seeds. Under the scan's single-loop shape $P_w$ lies between the
spinodal pressures, so both stable roots exist there; distinct-root, positive-slope
and branch admission are unchanged. A two-root $P_0$ is kept, and supplied
coordinates still replace generated ones. When this rule runs and finds no such
$P_w$, the start keeps $P_0$ and refuses with `no positive mechanical-loop witness
pressure`; starts the rule does not cover keep the generic coexistence refusal. A missing
witness does not prove the fluid supercritical or coexistence absent: a loop
narrower than the scan grid near the critical point has no sampled negative slope.
No pressure is invented and nothing retries; the rule was qualified on the near-critical
alcohol saturation states of the [Gross 2002 reproduction](../../../analyses/reproduction/2002-gross/README.md)
(PR #414). Final residual and stability checks are unchanged.

#### Default amounts

With more than one free phase, every free amount that no start owner estimates (no
Rachford–Rice or true-species split, speciated total, or caller or previous amount)
starts at the feed amount $N$ instead of $N/k$ for $k$ free phases (#407). This is an
infeasible default starting amount, not a material allocation: the starts hold $kN$ and
do not satisfy $Bn=b$, and Ipopt restores the balance. #243's [start-amount scale definition](#equilibrium-start-amount-scaling) is
unchanged, so these phases take $S_p=N$. Ipopt starts each bound multiplier at $z_y=1$,
and the feed-relative amount bound multiplier is $z_q=(N/S_p)z_y$: $k$ from the old
default, 1 from the new. That argument holds only for default-started free amounts of a
Gibbs minimization; in pinned, caloric or critical feasibility problems, whose
coordinates keep unit one, the change is only a different infeasible start. The
$\kappa_N$ first-step retention in the [amount step rule](#equilibrium-amount-step-retention) becomes $0.1N$ per phase instead of $0.1N/k$; it is
recomputed at every iterate, so it is no permanent floor. From $N/k$ the declared
esterification liquid pairs of the [Ascani 2023 reproduction](../../../analyses/reproduction/2023-ascani/README.md)
lost their minor liquid at 1, 100 and 10,000 mol; from $N$ they return the lower-G distinct
liquids (PR #421). This is an empirically qualified initialization convention, not a
convergence guarantee.

### Neutral recovery

<a id="neutral-recovery"></a>**Neutral recovery (#565).** An ordinary neutral, nonreacting, non-associating
declared L+V solve may perform one bounded continuation recovery only after coincidence, local phase
instability, local infeasibility or an optimizer iteration limit. A successful ordinary result returns
unchanged. Recovery uses the existing driver, start mapper, exact predictor/corrector and anchored phase
admission, with all failed initial, source and corrector work charged to the same remaining caller
budget. It cannot restart an exhausted budget, change phase identity, select another branch after a
terminal continuation event, or certify global stability. A failed recovery preserves both stopping
causes and publishes no recovered state.
The domain is a fixed-T or fixed-P boundary with one zero-pinned phase, or a fixed-T/P flash with both
amounts free, on the full positive feed support; it excludes discovery, critical, envelope, azeotrope, caloric and
source-reference declarations. Only the public solve without a caller start and cooling onset's first
solve recover; regression observations, the speciation seed, ordered starts and the recovery's own solves
do not. The triggers are a non-fatal Ipopt return at `LOCAL_INFEASIBILITY` or `MAXITER_EXCEEDED`,
coalesced phases, or a resolved negative lowest local-stability eigenvalue; a vanished phase, a
density-branch refusal and every other stop stand. Sources are ordinary solves at
$T_j=0.95^j\,T_\mathrm{seed}$ (`COOLING_RATIO`), $j=1,\dots,8$ (`EQUILIBRIUM_RECOVERY_SOURCES`), with guesses removed and an
absent free-pressure seed; the first admitted one is continued at fixed $T$ to the target, or, for a
free-$T$ boundary, first corrected at its own pressure in the free-$T$ chart (agreeing with the source
within the row tolerance) and continued at fixed $P$. Continuation steps are 0.05 to $10^{-6}$ with ten
retries; the source densities anchor its initial admission. The reached state passes the target's rows
and cold density-branch admission. Added optimizer iterations and continuation attempts are each at most
`EQUILIBRIUM_RECOVERY_WORK` (512), objective/polish evaluations four times and Hessian callbacks twice
that, within any caller remainder; an exhausted cap ends all recovery work. A refusal's message leads
with the ordinary cause; reported work is the ordinary work plus every recovery driver's.
This replaces, for this domain only, #352's rule that a no-start fallback is decided before the target
solve; the reacting and caloric no-retry rules are unchanged.

### Cooling-onset dew

`cooling_onset(compiled)` (ADR 0012 addendum) returns `(result, reports)` for a neutral,
nonreacting fixed-pressure dew: one liquid pinned at zero amount, one free vapor, full
positive supports and generated starts; other declarations are refused. The dew is the
first liquid to appear on cooling the vapor at fixed $P$ and $y=z$, the temperature where
the minimum liquid trial distance $\min\mathrm{tpd}/RT$ against the one-phase vapor reaches
zero from above. Equal chemical potentials alone also admit lower metastable roots, and
the trial distance at one temperature does not order their zero temperatures (row 558
from 375 K otherwise returns 369.425 K).

The ordinary declared solve runs once; a failure other than local instability or
infeasibility returns unchanged with an empty report list. At each sampled $T$ the
one-phase vapor reference is solved and admitted only with sufficient evidence that it is
vapor-side: the density owner's vapor-branch request at $(T, P, z)$ returns more than one
stable root, or its sole stable root lies below the scan's first negative pressure-slope
sample $\rho_-$ (the converse of the liquid-side feed test). A sole root above $\rho_-$, or
with no negative-slope sample, is unconfirmed, and the operation refuses before that
sample confirms an endpoint or sets a bracket (1-propanol(0.413)/ethylbenzene at
10,132,000 Pa from 392.05 K otherwise confirms a 227.63 K "dew" whose two phases are
both liquid-like; since #423 the ordinary solve's branch admission already refuses that vapor). The rule is sufficient, not a
necessary physical condition: without a fixed-composition mechanical loop the scan cannot
establish vapor-side identity, so valid mixture dew boundaries, for example near the
mixture critical region, can stay unconfirmed and be refused (the Gross 2002 row 173, 496
and 522 declarations raised to 6.0792, 3.039 and 10.132 MPa have dews at 563.390284,
625.552224 and 543.914200 K). #423 tracks broader reference admission. Each admitted reference then
runs the complete liquid trial family of `stability_test(...)` (at most $m+1$ trials). An ordinary success is sampled at its
own temperature and returned, with that report, when no trial resolves a descent; a
resolved descent marks a metastable root (row 722 from 475 K otherwise returns
368.061 K) and starts the search below. From the seed temperature,
heat by a factor 1.01 to the first positive minimum distance (at most 16 samples), then
cool in fixed steps of $0.005\,T_u$ to the first nonpositive one (at most 64). The
existing dew solve, started at the bracket midpoint with the lower sample's lowest-distance
liquid and $y=z$, locates the zero; its temperature must lie in the bracket within
$10^{-10}T_u$. The trial family is rerun there, and any resolved descent refuses. That is
at most 82 samples and 246 liquid trials for a binary. Each refusal appends its cause to
the original failure, or to "ordinary dew unconfirmed" after a success, and returns no state: no upper bound in 16 heating samples, no onset
in 64 downward samples, an unavailable vapor reference or liquid trial, an incomplete
trial family, an unsuccessful corrector, a corrector outside the bracket, or competing
descent at the corrected dew. All reference, density-scan, trial and corrector work is counted,
including work before a refusal.
Finite coverage: a narrow negative-distance interval between samples or a missed
composition minimum is not excluded, and a refusal means unresolved selection, not an
absent dew. PR #424 qualifies it on the water/1-pentanol dews of the
[Gross 2002 reproduction](../../../analyses/reproduction/2002-gross/README.md) from several start temperatures.

### Critical-equation candidates and phase envelopes

A critical cell (EqID: `critical_equations`, issue #224 B) declares one free
neutral phase on a positive full support, free T and P, a density seed and an
optional mode seed; reactions, ions, source references, discovery and caloric
rows are refused. Its unknowns add u between the phase block and the free
intensives: f = 2, k = 0, e = 2. It is a feasibility cell of the one driver
with constant objective. The geometry tape carries D of the moving feed, w = D u
and the ideal parts; the row owner adds d_i [F_nn(p) w]_i and F_nnn(p)[w,w,w]
along the path through Taylor order two (E2–E5). The weighted curvature of the
critical rows is the exact polarization H_ab = c2(e_a + e_b) − c2(e_a) − c2(e_b)
of their second Taylor coefficients, confined to critical cells; there is no
dense tensor or second row evaluator. The seed takes T and density from the
declaration, P from the state there and, without a mode seed, the lowest
eigenvector of M read from the mode rows' u columns.

Admission replaces the positive-slope density-root rule only for a critical
cell (a pure-fluid root has P_rho = 0; a mixture root generally does not):
the augmented Jacobian must have full rank, and
M (again the u-column block) must have one lowest eigenvalue within
|lambda| <= ||M u||/||u|| + C eps ||M|| at the solved point, above no resolved
negative one (the second eigenvalue exceeds that band). The record keeps the
spectrum, mode, cubic value, band, rank, condition estimate and a backward error
either way. The result is a critical-equation candidate, not a stability
certificate. Its first solved-state actions use E4 and its second-order actions the
order-two critical rows (E5); third order needs order three of the critical rows (E6) and is
refused as `energy_order_unavailable` with central values kept. Second actions of T, P and
density along a feed direction match Richardson-refined central differences of first actions
at re-solved critical points within 4.1e-11 relative (C3 binary, #341). A one-species candidate's
isobaric heat capacities divide by its zero P_rho and are unavailable
(`pressure_singular`); a mixture candidate's follow the ordinary pressure-slope
admission.

An envelope (EqID: `envelope_local_condition`) declares a free feed phase and a
zero-pinned incipient phase on one support, free T and P and one selected
condition: k = 1, f = 2, e = 1, raw T/P chart (Michelsen 1980, pp. 2–4). Its
continuation (`continue_equilibrium` with an envelope start problem) takes s, the
condition target, as the path parameter and first increases the declared
coordinate. At each accepted state the tangent solves J dy/ds = e (the condition
row has dF/ds = −1); the geometry owner's first Taylor coefficient along it gives
d ln T, d ln P, d ln rho_p and d ln K_i per ds. The largest of the dimensionless
|d ln T|, |d ln P| and |d ln K_i| selects the condition; a switch rescales the same
tangent and restarts s at the current coordinate value, held per accepted state
so a rejected attempt keeps its coordinate. The direction keeps the sign of the
dimensionless tangent (d ln T, d ln P, d ln K)/ds from the previous step, not the
raw chart product. Each attempt's density anchors are the Euler-predicted
densities ln rho_p + Δs d ln rho_p/ds, and the same driver corrects inside the
anchored window W = `DENSITY_ANCHOR_LOG_WINDOW` around them. The first step is
at most 0.9 W / max_p |d ln rho_p/ds| and the caller's maximum log step. After
each accepted step of size h, e = max_p |ln rho_p − ln rhô_p| is the measured
predictor error, and the next step is h min(2, sqrt(W/(2e))), capped by the
maximum log step: the Euler error grows as h², so this targets e ≈ W/2 (#571).
A rejected ordinary step halves from the attempted step down to the caller's
minimum step and retry allowance, including a corrected density outside the
anchored window (`branch_lost`) and a corrector that returns no state because
it ended on coincident phases (the near-critical approach to the straddle
threshold at cap 0.05). A singular tangent and an identity swap between
accepted states stay terminal.

A mixture passes its critical point once, by a straddle in ln K. When the selected
ln K_k is decreasing in magnitude (s times the step sign is negative), an ordinary
step lands at most at the same-sign threshold; from |s| <= threshold the straddle
requests -s, never s = 0. That attempt already expects the feed on the other
density branch: the phase kinds are swapped, and the predicted anchors carry
the densities through the crossing.
It passes only with two distinct accepted states whose pair identities swap.
Both endpoint tangents in the same s give cubic Hermite values of ln T, ln P
and both ln rho at s = 0; their mean log density and the common-volume
displacement u = D^-1 dn, dn = V*(rho_inc x_inc - rho_feed x_feed), seed the
critical cell solved by the same driver (retained in `critical_points`, seeds in
`critical_seeds`). After the passage an accepted step that brackets the initial
pressure is followed by a same-driver corrector at that pressure on the other
branch: closure. An identity swap without a requested straddle is terminal
coalescence, and a refused straddle, closing corrector, passage tangent or
candidate ends the trace as a diagnosed prefix.

A pure fluid has ln K = 0 and no straddle. Its saturation trace bounds each
predictor step by (3/8) g/|dg/ds|, g = |ln rho_L - ln rho_V|, and stops at the
first accepted state with g <= 4 W: below about 2 W the anchored window refuses
both roots. The critical cell seeded there (T, P, mean log density and the
displacement mode) ends the trace as `critical_endpoint`, with the last
distinct saturation state retained; it is not closure on another branch.
Liquid–liquid, reactive, charged, three-phase and disconnected envelopes,
critical stability and E6 remain outside this method.

An azeotrope cell (`Problem(azeotrope=True)`, EqID: `azeotrope_condition`) is a
binary bubble/dew cell on the full support of a positive feed: one free and one
zero-pinned phase of opposite liquid/vapor kinds, exactly one of T and P free,
and no other condition, ion or reaction. In the species-identity balance basis
the last row becomes the all-ones row (B' = T B, so λ'·B' = λ·B leaves the
chemical rows unchanged); only that total-amount row is emitted, and the other
C − 1 balance rows are replaced by ln x_k^inc − ln x_k^feed = 0 for the first
C − 1 support species. The layout and degrees of freedom are the bubble cell's;
the declared feed composition is only the seed, and the solved feed-phase
composition is the azeotrope. Published multipliers are in the B' basis. The
cell is regular at a simple x = y crossing (by Gibbs–Duhem the free T or P
column lies outside the composition columns' range, Δs_vap or Δv_vap ≠ 0) and
Ipopt solves it alone (ADR 0009, #571). A refusal means no azeotrope was
reached from this start, not that none exists. Every critical point and
same-root state also satisfies x = y: the topology predicate refuses
coincident phases, and local phase stability refuses the near-pseudo-critical
states it misses. Qualification covers homogeneous binary VLE azeotropes
against FeOs `binary_azeotrope` (A1); uniqueness, global stability,
heteroazeotropes, C > 2 and azeotrope-line continuation are outside it.

### Pure stoichiometric solids

A free solid's amount n_s is an Ipopt column with lower bound zero, no upper
bound and no stationarity row: its bound multiplier carries the complementarity
of EqID `solid_amount_complementarity`. One classification owner acts at the
Ipopt endpoint by complementarity: present when the feed-relative amount exceeds
the bound multiplier, n_s/N_feed > z_s (z_s = −A_s/RT with A_s/RT = λ·Bν_s −
g_s/RT; row g_s/RT − λ·Bν_s = 0), absent otherwise (n_s = 0). This replaced an
affinity test, present when A_s/RT ≥ −ε against the caller row criterion: the
barrier leaves a stable solid at −A_s/RT ≈ N_feed μ/n_s, above ε for small excess
feeds (6.29891–6.5 mol/kg halite), so the stable solid started absent and, with a
second free solid, both flipped and the round failed (#573 review r1). The
existing polish runs on that square classified declaration. Afterwards a present
solid must keep n_s > 0 and an absent one A_s/RT ≤ ε; every violating solid
flips in one round and the declaration is polished once more. Near onset
(halite feeds up to 1e-4 mol/kg above saturation) the barrier amount is still
unresolved and the stable solid takes that flip. A later violation is the typed
failure naming the solids; a singular classified Jacobian (dependent Bν_s among
present solids, or more present phases than the phase rule allows) is a singular
failure naming the present solids, as for halite with an identical free copy.
Admission does not refuse dependent stoichiometries, because polymorphs are
physical: a free copy of halite with ln K + 0.01 stays absent at SI = −0.01/ln 10.
An absent solid's row n_s = 0 is a constraint of the reduced-Hessian check. The Newton
attempts are ineligible (a solid requires a source reference). Solved-state
actions rebuild the classification from the central amounts (present exactly
when n_s > 0) and are one-sided at the onset kink; T and P directions are
admitted through order one, feed directions through the general order. The
solid's ln K^ρ is compiled once at fixed T at the source pressure, with its first
T slope where the reference admits it; a scalar K refuses T actions.

### Cation-mixing coordinates

The optional `cation_pair_mixing` and `cation_triple_mixing` scalar families
use the affine chart with target `None` and their declared empty-or-one slots.
The common Helmholtz tape owns their exact columns, including at zero
coefficients and on an absent-cation face. Native construction owns topology,
units and coefficient admission; no observation-level activity correction is
applied. EqID `cation_pair_triple_mixing_helmholtz` states the liquid domain.

### Declared measurement uncertainty

A caller-declared Gaussian measurement-error model augments the existing
Ceres least-squares fit with standardized physical input corrections (EqIDs
`errors_in_variables_objective`, `independent_feed_measurement_coordinates`,
`errors_in_variables_residual_derivative`). Output differences divided by
declared physical standard deviations and one input penalty per shared
measurement form the single objective. Exact native state/reference actions
provide columns. Zero input sigma eliminates that coordinate; absent or empty
declarations request no input action and retain ordinary fitting work.
The [formulation](../regression/science/pure-saturation-regression.md#declared-measurement-uncertainty)
owns admitted observation/input combinations, physical units, Difference-only
weighting, the per-group reference-action limit and exclusions. One optional
known covariance or correlation spans the positive-weight output errors, then
the active input corrections, across measurement groups. The measurement chart
admits it before fitting (size, finite entries, unit diagonal and symmetry
within 64 m epsilon in standardized units, positive Eigen LLT pivots) and
factors it once; it is never jittered, repaired or pseudo-inverted. After every
successful row evaluation the parent solves L z = e on the gathered compact rows
of the residual and of every complete column, scatters them back and leaves
excluded rows zero; a nonfinite result is a trial failure. Raw prediction
columns are formed from the physical actions and the coordinate map, never from
whitened rows. Requests without a matrix keep their arithmetic bitwise.
Unsupported declarations and nonunit positive relative weights are
refused. All model/input coordinates share Ceres controls, the declared box,
iteration and elapsed budget. The one admitted measurement-membership map
supplies allocation, worker partitioning, penalties and result identities;
equal states do not establish shared measurements. Same-point promotion keeps
its accepted state. Changed relevant physical inputs withhold optional prior
density seeds; input movement rebuilds trial state/reference declarations.
Failures publish no new trial cache or start and preserve their native cause.
Local covariance uses the complete standardized (whitened, with a matrix)
Jacobian and eliminates input directions; it is unavailable at an active bound,
deficient rank or unfinished fit. Known absolute sigmas or covariance supply no
estimated common variance factor.
Reconciled-input predictions and original measured-input checks retain separate
meanings. Fitted model parameters and inputs remain candidates. The frozen
synthetic comparison is retained in `analyses/2026-errors-in-variables/` when run;
the known-covariance comparison is in `analyses/2026-correlated-measurement-errors/`.

### Optional fitting-density offers

Optional fitting densities come only from the previous fully successful native
evaluator trial. Exact closure identity, the complete 10% physical-chart gate,
and regularity b <= w/100 admit an offer. Changed parameters additionally need
a retained negative-slope witness rechecked at the target model, with the whole
anchor window on the requested branch side. All thermodynamics and columns are
recomputed. A caught optional-anchor failure is recorded before exactly one cold
retry; a final cold refusal keeps its complete cause and publishes no trial seed.
The witness identifies the physical extreme in the admitted single-loop domain;
agreement with the bounded cold search is a qualification result, not a proof.
Since #275 a cold scan with a failed grid node refuses a sole root, while a
witnessed anchor can still close; such a seeded state has no cold comparator.
Qualification compares matched returned observables, densities, sites and exact
columns within 1e-10 relative, with an absolute floor of 1e-12 near zero. The
osmotic coefficient carries the ln a_w floor as max(1e-10 abs(phi), 1e-12/s),
s = M_w sum(nu) m, which grows as 1/m. The independent pressure-slope
derivative check uses 1e-6 max(1, abs(reference)).
Qualification uses independent same-point cold replays; the fitting controller
adds neither a per-trial nor an endpoint-only cold comparison.
Seeded residuals differ from cold ones within these floors and depend on which
trial seeded them. Once a fit's cost reaches that level, seeded and cold fits
can accept different steps and stop at different points; tighter stopping
tolerances do not remove this endpoint difference (#245 round 4).

### Guarded minimax

The constrained path minimizes $F=t/q$, with $t\geq0$ in decades and declared
positive $q$. Targets use unweighted $e=\log_{10}(y/y_{obs})$ and paired rows
$e-t\leq0$, $-e-t\leq0$. An enforced maximum uses the same pair with its frozen
threshold. Each mean-group member has its own nonnegative slack in the error's
unit, paired $e-s\leq0$, $-e-s\leq0$ rows, and one group-mean row. Overlapping
groups retain distinct slacks. Final means use the actual absolute errors.
Report-only rows do not enter trial evaluation; enforced rows enter regardless
of Ceres weights. Baseline thresholds use the separate stored source parameters
and the existing check owner, without solving a hydration hold.

One parameter-fitting TNLP consumes `ObservationEvaluator` values and exact raw
and complete held-chart columns. The same check error/slope producer assembles
both estimators' reports. Every inequality is divided by its declared positive
scale; zero thresholds remain zero. The complete map includes all target coefficients, numerical
scales and eliminated Born rows. The epigraph and slack columns are analytic.
Only exact first actions are supplied; Ipopt uses limited-memory BFGS, with no
finite-difference derivative checker or parameter Hessian claim. Equilibrium
retains exact Lagrangian curvature. The common checked Ipopt option/lifecycle
owner uses the [diagnostic options-file route](#ipopt-diagnostic-options) for both Hessian modes.
Required positive scaled tolerance,
iteration/CPU limits and a steady-clock elapsed callback bound each estimate.
`acceptable_iter=0` and `bound_relax_factor=0` preserve the declared limits.

At the finalized point, recompute constraint/bound feasibility, stationarity,
dual signs and complementarity from exact columns and its own multipliers.
Return scaled and physical inequalities, row labels/scales, physical epigraph
and slacks, and parameter/auxiliary bound distances separately. For scaled
$c=C/\sigma$, physical-objective multipliers are $q\lambda/\sigma$ and $qz$ for
bounds. For a regular interior eliminated diameter, reconstruct
$\nu=-(\partial L_{without\ H}/\partial d)/H_d$ from raw physical columns and
report the reconstructed held and independent full-space stationarity equations.
The held component vanishes by construction; independent multiplier qualification
uses the manufactured analytic KKT conditions. Singular
or bracket-boundary holds retain the reduced diagnostics and an explicit
unavailable full-space cause. Root-error estimates use exact held-row slopes
and the achieved diameter uncertainty $\min(|H/H_d|,\ \text{bracket width})$,
counted once; refinement qualification belongs to the retained study.

Termination, point/action availability, feasibility, final checks and
stationarity are distinct. Final checks evaluate report-only rows omitted from
the minimax trial; their raw Jacobian rows carry no evaluated columns.
Final checks use the stored fitted diameters,
without resolving holds. A fatal initialization has no candidate; local
infeasibility proves no global infeasibility. No covariance, unique optimum,
global certificate or physical validation follows from this local minimax
estimate. `evaluate_minimax` evaluates the same compiled rows at a declared
point without optimizing or publishing a stationarity certificate.
