# Scientific algorithms

Accepted requirements, reconciled against the 2026-09-06 transfer audit.
These methods consume the [equation definitions](../equations.md),
[EOS conventions](full-eos.md) and [equilibrium formulation](../equilibrium/formulation.md).
Algorithm IDs below are stable references. Current implementation and required
numerical falsifiers are mapped in the owning [GitHub issues](https://github.com/tannerpolley/ePC-SAFT-greenfield/issues).
No section claims that all required actions currently exist. The implementation
plan owns current versus deferred scope: PH/PS/UV and phase discovery are future
work. Reuse current/old numerical owners and independent cases under its
adoption gate before implementing a method. Declared-branch continuation and
mechanical/local stability checks remain current requirements.

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

## ALG-DENSITY

Inputs include T, composition, target pressure, admissible packing interval,
root-selection request, numerical policy and declared search budget.

Density closure requests only the reduced pressure gradient, a_VV and the
association floor. With association sites these come from one adjoint sweep
(the same site multipliers as the Hessian path) and one volume tangent through
the site sensitivity; the full site Jacobian is formed only for Hessian requests.

A neutral reference trace at system pressure first tries the compile-time
pure-solvent liquid density as an anchor for its pure-solvent point: the
anchored solve must find a locally unique stable root in its window. Any anchored
failure falls back to the cold scan and its branch admission; later trace points
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
per cell, discard an entirely negative-slope interval; bisect any slope sign
    transition to adjacent representable densities or an exact slope zero,
    stopping early at a stable sample whose pressure passes the target
if no stable sample brackets the target and it is within the existing
    reduced-residual tolerance of that turning-point pressure, report an
    unresolved spinodal
trim the cell to its positive-slope side and refine each strict pressure
    crossing once with the existing safeguarded Newton method
select by declared branch policy and return the stable-root count
```

The supported shape is one mechanical loop per isotherm, with one interior
minimum of `P_rho`. Once a negative-slope sample exposes that loop, its two
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
nonexistence. The root residual tolerance is unchanged. For an approximate
root, refinement also requires `(rho_upper - rho_lower) |dr/d rho| <= tolerance`
with a finite, nonzero derivative. Once the residual passes, the next sample is a
quarter radius past the Newton root toward the far bracket end, so a sign change
closes the bracket inside the radius in about two evaluations; midpoint
bisection remains the fallback when that sample leaves the bracket. Thus the represented
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
A genuine sole root remains unlabelled. An anchored request uses the same
partition and counting method within its existing local window and keeps the
requested branch only after local uniqueness is established. That window
does not establish uniqueness outside itself or global branch persistence;
a failed anchor must not silently switch to an unrelated cold root.

## ALG-EQUILIBRIUM

All admitted specifications use the one equilibrium driver required by the
[project map](../../../docs/scientific/CONTEXT.md). Compile balances, supports,
pins and free intensives from data. Check structural
counts/rank, construct declared feasible starts, then solve the same intensive
rows using the required common energy. Recompute achieved residuals and phase
distinctness at the returned point; optimizer termination alone is insufficient.
For a neutral, non-associating, nonreactive, full-support two-phase fixed-TP request with both
amounts free, absent complete caller composition guesses, the existing feed-state
fugacity estimates define $K_i=\phi_i^L/\phi_i^V$. When finite positive $K_i$
bracket an interior root, bisect
$F(\beta)=\sum_i z_i(K_i-1)/(1-\beta+\beta K_i)$ on $[0,1]$ through the
binary precision of a double. Seed $x_i=z_i/(1-\beta+\beta K_i)$,
$y_i=K_i x_i$, $N_V=\beta N_{feed}$ and $N_L=(1-\beta)N_{feed}$, assigning
by declared phase kind. Nonfinite/nonpositive coordinates or fractions at the
existing phase-disappearance floor retain the original normalized $z/K,zK$
and equal-amount seed. No new phase search or equilibrium driver is introduced.
Only this non-associating neutral TP domain may use a bounded caller packing seed or skip
redundant density initialization for a previous successful result with compatible
positive coordinates. A multicomponent bubble declaration (fixed temperature, free
pressure, one liquid and one vapor pinned at zero amount) moves the caller's
pressure seed by one fugacity step,
$P_1=P_0\sum_{i\in V}x_i\,\phi_i^L(T,P_0,x)/\phi_i^V(T,P_0,y)$, from the seeded
compositions; the pressure coordinate and every density seed use $P_1$, and a failed
state keeps $P_0$. Liquid fugacity coefficients scale as $1/P$, so $P_1$ is the
bubble pressure of the seeded liquid; its accuracy is that of the seeded $x$. Other
declarations retain their original initialization, including the pure-fluid
coexistence check. Final residual and stability checks are unchanged.
Public solves also apply the continuation owner's density-branch admission to
each returned phase. A cold pressure closure selects its declared liquid or vapor
root, an anchored closure establishes local uniqueness, and the corrected density
must lie in that anchor's window with positive pressure slope. Failure clears the
state and retains the optimizer and row diagnostics. Continuation uses the same
predicate with its own anchors and retains its `branch_lost` cause; its internal
corrector does not independently reselect a cold branch. This is declared-branch
admission, not global phase discovery, and inherits ALG-DENSITY's finite-search
limits.
Retain the native optimizer `solver_status` alongside the convergence stop,
residual, and failed-attempt initialization/termination details. The optimizer
stops on amount-weighted rows, which bounds an intensive trace row only by
`tau N_feed/(N_p x_i)`; after a success, acceptable-point or tiny-step stop, at most
`EQUILIBRIUM_POLISH_STEPS` Newton steps on the unweighted rows with the exact
Jacobian close that gap. Each step must contract the raw row maximum and keep
every composition above the structural cutoff `STABILITY_COMPOSITION_FLOOR`,
below the optimizer box when the solution lies there (#87, 2026-09-21).
Try fractions 1, 1/2 and 1/4 of each Newton direction, stopping at the first
admissible raw contraction; a rejected direction ends polishing. Among the
optimizer start and accepted physical iterates, return the point with the lowest
classified residual, breaking ties by lower raw residual (#66, 2026-09-22).
The classified residual may rise at an intermediate step through nonlinear
curvature; requiring it to decrease at every step can reject a direction that
resolves trace chemical equilibrium in the existing three-step budget.
Selecting the best classified point prevents degradation of the optimizer
result or loss of a stricter caller tolerance attained along the path.
Acceptable-point and tiny-step stops use the same polish. No additional step or
relaxed criterion is used. Recompute diagnostics at the selected point, and
apply the unchanged caller criterion and all final physical/branch checks.
A square pinned cell carries chemical and pressure rows as optimizer equality
constraints. Its log-composition bounds therefore use the
structural cutoff directly: the all-free optimizer box would exclude valid
trace roots before a successful optimizer stop can reach polishing. In all-free cells a warm start below the box
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
2026-09-21). A pinned problem keeps the optimizer's constraint multipliers
over every row; its equality-root qualification does not assert a Gibbs minimum.
The reconstruction assumes no active
variable bound at the returned point and residual potentials that depend on
composition and molar volume only. The charge-row scaling `-n_p phi_p` is
qualified by independent extensive-objective stationarity and public ionic
chemical potentials at Khudaida Figure 2 point 1 with 1 and 10 mol feeds
(`test_trace_reliability.py`); application-specific ionic qualification remains
open. The Ipopt driver
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

The caller supplies one dimensionless row tolerance per solve request; the
named `EQUILIBRIUM_ROW_RESIDUAL` policy is its only default. That value alone
classifies the recomputed rows; the Ipopt `tol`, `acceptable_tol`,
`dual_inf_tol` and `constr_viol_tol` settings are the smaller of the request
and the default, times `EQUILIBRIUM_INNER_TOLERANCE_FACTOR`, because the
recomputed row maximum tracks the optimizer tolerance and overshoots it by a
hair when both coincide. A caller-valued criterion is owned
by the convergence record so copies into results, traces and exceptions keep
it alive. The result reports `requested_tolerance_met` exactly as recomputed
residual ≤ requested tolerance whenever a finite residual exists; an available
floor-limited state may miss the request and keeps its state with a warning.
A structurally square declaration whose generic-point Jacobian loses rank
returns `NOT_RUN` with `Stop::infeasible`, an unmeasured residual, no
criterion and the row, unknown, rank and nullity counts, before any start or
optimizer work. Within one solve the driver prepares each distinct callback
point once and reuses its phase evaluation for rows, Jacobian, gradient and
curvature; the weighted constraint curvature is recombined when the optimizer
multipliers change, and the final row certification is a fresh evaluation.
That fresh final preparation also supplies the reduced-Hessian certificate;
only the coordinate width is reduced to the primal chart, without reevaluating
the same phase energies separately for its Jacobian and Hessian. Diagnostic
residual summaries reuse the certified row values. For a no-site EOS Hessian,
each unit-direction `Forward(1)`/`Reverse(2)` sweep supplies a Hessian column;
the first sweep's even entries supply the gradient, including active-parameter
coordinates. The initial value sweep and the association implicit chain remain
unchanged; a separate no-site first-order reverse sweep is unnecessary.
Standard-state terms carry exact value, first and second temperature
derivatives from the reaction target owner; no affine shortcut is assumed.
For an unpinned cell the chemical rows are not optimizer constraints: they are
recovered from the objective gradient `N x_i mu_i` through least-squares balance
multipliers, so their unscaled residual is bounded only by the dual tolerance
divided by `N x_i` and diverges where a phase or a fraction vanishes. The driver
therefore classifies the free-amount route by the KKT-scaled stationarity
residuals `N_p x_i mu_i / N_feed` and `N_p / N_feed` times the pressure row (the
first-order conditions in the solved chart; owner decision 2026-09-16), keeps
the unscaled rows in `residuals` as diagnostics, and classifies a pinned cell by
its unscaled rows, which are its constraints. Species thermochemistry attached to a
model is compared with every target kind, including scalar `K` and affinity,
and the comparison is skipped where the temperature leaves a species interval.

Declared-branch continuation follows one straight dimensionless path `q` in
`[0, 1]` between two problems with identical declarations; only fixed
intensives and resolved feed amounts move. The controller predicts each step
from the implicit tangent (EqID: `implicit_sensitivity` with the existing lever
Jacobian and linear solve), corrects with the same Ipopt driver seeded by the
full predicted coordinate vector (a predicted free amount below zero is
projected onto its bound, the disappearance boundary the bound-constrained
corrector then reports; every other coordinate outside its domain is a
rejected left-domain step), and accepts a state only when the driver
accepts it, every recomputed row meets the caller `row_tolerance`, the
candidate-topology predicate reports no event, and each corrected phase density
shares the anchored locally unique root window (`DENSITY_ANCHOR_LOG_WINDOW`) of
its cold declared root at `q = 0` or of the previous accepted density
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

## ALG-ACTIONS

An action request identifies observable, independent coordinate order, units,
directions, derivative order, source domain and branch/support. A successful
central solve does not imply every action is available. Return per-action
values or a specific failure reason, with rank/conditioning and linear residual
where relevant. Cache/reuse thermodynamic work for the same resolved state;
reuse must preserve all coefficient, T, density and implicit-site dependencies.

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
Composition directions preserve normalization, charge balance and exact support;
known fixed parameter slots project to zero, while unknown slots are input
errors. The [accepted P4 design](https://github.com/tannerpolley/ePC-SAFT-greenfield/issues/19)
specifies the qualification cases and numerical limits.

| Caloric family | Value energy order | First action | Second action | Third action |
|---|---:|---:|---:|---:|
| `h`, `u`, `s` | 1 | 2 | 3 | 4 |
| `cp`, `cv` | 2 | 3 | 4 | 5 |

Use the same residual-potential/mass-action tape and a single factorization of
the accepted unweighted site Jacobian for the batch. Obtain each implicit site
Taylor coefficient from the defining mass-action coefficient, then reuse each
directional path and its lower coefficients across all selected properties.
For an E3 energy action with exactly two equal normalized directions, use the
sparse Taylor paths in EqID: `repeated_direction_energy_third` instead of linear
polarization. This avoids subtraction of a singular pure-composition third
derivative from a regular mixed volume response at dilute ionic references.
All-identical triples, nonrepeated triples and E4 keep their established routes.
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
outside the per-result typed-status contract. Branch and tape-predicate changes
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

Solved-state actions (E5) are evaluated by one native owner,
`solved_state_actions(compiled, central, request)`, at an accepted result; the
owner never re-solves, compiles or mutates. Admission requires the central
classified residual, the quantity the driver certifies (owner decision
2026-09-16), at or below `min(row_tolerance, EQUILIBRIUM_ROW_RESIDUAL)`; no
candidate-topology event and every supported composition above
`STABILITY_COMPOSITION_FLOOR`; a finite square full-rank row Jacobian,
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
owner). Each solved coefficient is re-evaluated; its normwise backward error
eta_k = |r_k| / (|J| |y_k| + |b_k|) (infinity norms, Rigal–Gaches) must meet
`EQUILIBRIUM_ACTION_BACKWARD_ERROR` (`linear_solve_residual`). This certifies
only that y_k solves a nearby linear system; it catches a stale Jacobian or a
jet mismatch and is invariant to direction and chart scale. It does not bound
forward error: that is about cond(J) eta_k for the largest coefficient, so an
availability claim for an ill-conditioned state needs the caller's
finite-difference ladder with the reported condition estimate, including small
coefficients. Observables are
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
after admission, and a fixed T or P keeps its identity action. `support_change`
and `ionic_insertion` are unforced: a zero-feed non-participant species leaves
a rank-deficient balance row, so no admitted central state carries a removed
species.

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

## ALG-REFERENCE

Use the explicitly declared solvent/reference path and source standard state.
The two-ion formula must contain charged components with zero stoichiometric
charge; the solvent composition must also be neutral within the composition
admission tolerance. The shared composition owner enforces this before either
reference construction or activity evaluation; neither API supplies an ionic gauge.
Evaluate the retained molality sequence; return the terminal change and
criterion separately. Unknown validity or failed limiting convergence remains
unavailable. Differentiate both the source correlation and EOS reference path,
including its density closure. Fixed source pressure P0 must not become trial
pressure. Reference-gauge choices that preserve reaction affinities do not
supply absolute caloric data.

For an explicitly referenced neutral reaction, use
`neutral_reaction_reference_transfer`: one pure neutral solvent, positive
all-solute molality ratios with zero net charge, and one declared source basis
per reaction. At each resolved source pressure, evaluate the five molalities
from 1e-6 through 1e-14 mol/kg on the requested liquid branch. The density owner
admits a unique mechanically stable cold root with an unspecified branch label;
this is not a physical liquid/vapor classification. Subsequent anchored roots
retain the requested branch. Require finite
states, positive pressure slope, branch continuity and both last contraction
changes within 5e-5. Retain the terminal contraction and observed change.
Do not extrapolate this sequence: electrolyte terms can scale as sqrt(m).
Qualify the finite estimate against a second charge-neutral dilution path;
agreement is an observed domain check, not a universal limit-error certificate.

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
The implemented slice admits fixed T/P central solves and feed actions, and
square pinned free-pressure cells. Pressure/feed actions are supported through
order two; fixed-T/P feed-only order three remains available. Referenced
continuation, free temperature, raw lever derivatives, nonzero temperature or
parameter actions, and third-order actions involving pressure remain
unavailable. Ideal/Total caloric values and actions remain
unavailable even with an attached thermal record, while independently valid
residual properties remain usable. A central finite reference cannot certify
the missing reference derivatives.

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

## ALG-FIT

Use the [regression formulation](../regression/science/pure-saturation-regression.md)
for observation algebra, weights, parameter coordinates and uncertainty.
Resolve each state once, evaluate its declared observables/actions, transform
residuals and Jacobians consistently, then apply the declared bounded fitting
method. Initialization, bounds, loss, covariance assumptions and budget are
scientific inputs to the fit. A failed state is a failed evaluation, never a
successful zero residual. Fitted values remain candidate parameters.
The implemented fit (R1--R3) uses the built-in observations and affine, tied
optimizer coordinates of the formulation, a Ceres Levenberg--Marquardt trust
region with dense QR on one thread, the named `FIT_*` policy tolerances as the
only defaults, exact owner actions for every Jacobian and a covariance that is
withheld under rank deficiency or an active bound.

## ALG-TERMINATION

Core CORE1 keeps the numerical domains separate. `EQUILIBRIUM_COMPOSITION_BOX`
is the lower box of the intensive `log(x_i)` coordinate inside the all-free
optimizer subproblem: its barrier regularises the log chart's trace steps and
is not a solution domain. Square pinned equation systems use the structural
cutoff instead, so chemical equalities can reach trace roots below that box.
`STABILITY_COMPOSITION_FLOOR` is the structural-absence
cutoff for continuation admission, solved-state differentiability and the
Newton polish, which carries a trace solution below the box down to the cutoff
(#87, 2026-09-21). The delivered falsifier reaches x ≈ 1e-12; continuation
admission and action availability between 1e-12 and the cutoff are unexercised
and remain a qualification limit. Free-phase
disappearance compares the dimensionless amount ratio
`n_phase / sum(n_feed)` with the provisional `STABILITY_PHASE_AMOUNT_FLOOR`;
that topology rejection is not a stability or accuracy conclusion. Packing
admission rejects the equal-sphere face-centred-cubic geometric boundary
`eta >= CLOSE_PACKING_FRACTION`; this is a conservative PC-SAFT domain
constraint, not a mathematical singularity claim.

Numerical policy values stay in `native/core/numerical_policy.hpp`; source and
application accuracy criteria stay beside their accepted reference cases.
Keep the achieved error, criterion, scale, iteration/evaluation budget, native
solver status and termination reason visible. Only Ipopt `MAXITER_EXCEEDED` and
`CPUTIME_EXCEEDED` count as equilibrium work-budget exhaustion;
`STOP_AT_TINY_STEP` or `STOP_AT_ACCEPTABLE_POINT` with a finite non-converged
residual is unavailable, while a finite criterion-meeting result remains
eligible after the same bounded polish and physical checks. A tiny-step or
acceptable-point status alone never establishes success; floor-only results
from either remain unavailable. Recomputed residual and topology criteria apply.
Distinguish converged, stagnated, exhausted,
nonfinite, left-domain, singular and infeasible outcomes. Optional traces aid
diagnosis and do not determine success. Density root producers classify all
optional diagnostics before publishing a criterion-meeting root; nonfinite
floor, contraction or amplification is unavailable and throws the existing
convergence failure.

Propagate inner-solve/input error into an outer accuracy estimate using the
applicable local conditioning and scale. A first-order estimate or pivot ratio
is not a rigorous error bound or a general condition number. Do not claim more
accuracy than supplied inputs support, and do not loosen accepted thresholds
to force convergence. Source rounding, numerical error, experimental uncertainty
and model discrepancy are distinct. Missing accuracy/runtime thresholds stay
owner decisions; historical microbenchmarks are not universal limits.
`Stop::nonfinite` is unavailable and has precedence over finite convergence and
floor predicates whenever the residual or an optional diagnostic is non-finite.
For equilibrium, `success` requires accepted solver status, finite residuals and
declared topology; an unsuccessful result retains `solver_status`, residuals,
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
