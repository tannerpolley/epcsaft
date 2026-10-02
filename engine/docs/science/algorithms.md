# Scientific algorithms

Accepted requirements, reconciled against the 2026-09-06 transfer audit.
These methods consume the [equation definitions](../equations.md),
[EOS conventions](full-eos.md) and [equilibrium formulation](../equilibrium/formulation.md).
Algorithm IDs below are stable references. Current implementation and required
numerical falsifiers are mapped in the owning [GitHub issues](https://github.com/tannerpolley/ePC-SAFT/issues).
No section claims that all required actions currently exist. The implementation
plan owns current versus deferred scope. PH/PS/UV remain future work. The
explicit native fixed-T/P discovery and tangent-plane testing path in
ALG-STABILITY has retained numerical qualification on the stated neutral and
admitted ionic/reacting cases; independent delivered review is pending. Reuse accepted numerical owners and independent cases.
Declared-branch continuation and mechanical/local checks remain required.

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

`BornDiameter`, `DebyeHuckelDiameter` and `PackingDiameter` own positive
Angstrom values in charged component slots. `IonSuppressionCoefficient` owns
one dimensionless charged-slot value; a common C ties all declared charged
slots in the existing affine chart. Their values use `ParameterRecord<Scalar>`
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

## ALG-DENSITY

Inputs include T, composition, target pressure, admissible packing interval,
root-selection request, numerical policy and declared search budget.

Density closure requests only the reduced pressure gradient, a_VV and the
association floor. With association sites these come from one adjoint sweep
(the same site multipliers as the Hessian path) and one volume tangent through
the site sensitivity; no request forms the full site Jacobian.

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
pins and free intensives from data. The balance basis $B$ spans the left null space
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
Check structural
counts/rank, construct declared feasible starts, then solve the same intensive
rows using the required common energy. Recompute achieved residuals and phase
distinctness at the returned point; optimizer termination alone is insufficient.
For an unguessed liquid/vapor start, a sole unlabelled feed-density root cannot
supply both fugacity estimates: their self-ratio would make every $K_i=1$
and seed identical phases. In that case only the vapor fugacity estimate uses
the ideal-gas limit $\ln\phi_i^V=0$ (EqID: `log_fugacity_coefficient`,
$\mu_i^r\to0$, $Z\to1$). The resulting compositions still receive actual EOS
density closure on their declared branches; this estimate changes neither
the equations nor final admission. Two-root estimates and caller-guessed
compositions retain their existing rules. The [cold-vapor diagnosis](../../../analyses/2026-cold-vapor-branch/README.md)
retains the sole dense feed root and the recovered incipient vapor (#204).
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
compositions; the pressure coordinate and every density seed use $P_1$, the
incipient vapor seed becomes the same step's
$y_i=x_i\,\phi_i^L/\phi_i^V\,P_0/P_1$ over its support when the caller-guessed start is used
(every phase guessed; a vapor guess adds no independent information once the
liquid is given, while the default $z/K$ liquid would return $y=z$, the trivial
solution, and keeps its vapor seed), and a failed state keeps $P_0$ and the
vapor guess. Liquid fugacity coefficients scale as $1/P$, so $P_1$ is the
bubble pressure of the seeded liquid; its accuracy is that of the seeded $x$. A cold
reactive declaration at fixed temperature and free pressure with a phase pinned at
zero amount first solves its finite phases alone, at the caller's pressure seed and
with the same reactions and references; at the boundary those phases are the
solution's phases up to their pressure dependence, so their speciated compositions
replace the caller's guesses and an unguessed incipient phase takes the first
finite phase's composition on its support; the start then follows the
all-guessed rule above (#152: an unspeciated MEA liquid otherwise sends the
fugacity step to the wrong liquid's bubble pressure and the solve into feasibility
restoration). A caller start skips this seed. If the sub-solve, the seeding, the
seeded solve or its density-branch admission fails or throws, the caller's own
start runs unchanged; only the final attempt's diagnostics are returned. A reacting
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
Other declarations retain their original initialization, including the pure-fluid
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
`tau N_feed/(N_p x_i)`; after a success, acceptable-point, tiny-step, restoration-failure,
iteration-budget or time-budget stop, at most
`EQUILIBRIUM_POLISH_STEPS` Newton steps on the unweighted rows with the exact
Jacobian close that gap. Each step must reduce the raw row maximum and keep
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
hair when both coincide. A pinned square system (constant objective, followed
by the exact Newton polish) stops at the criterion itself; the polish, not an
inner decade, removes the overshoot (owner decision 2026-09-23). A caller-valued criterion is owned
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
from the implicit tangent (EqID: `implicit_sensitivity`: the first row
coefficient along the lever at the fixed state, the directional series of
ALG-ACTIONS, so source references carry their T/P chains), corrects with the same Ipopt driver seeded by the
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

## ALG-STABILITY

The explicit native tangent-plane test and fixed-T/P discovery adopted in
[#223](https://github.com/tannerpolley/ePC-SAFT/issues/223) reuse `CompiledSystem`,
`run_driver`, density/association owners and exact CppAD value/derivative chains.
No second minimizer, production finite differences or Python search loop is
admitted. The [formulation](../equilibrium/formulation.md#discovery-and-global-claims)
owns physical admission and claim meaning; EqID: `tangent_plane_distance` owns
the reduced molar objective. The [analysis](../../../analyses/2026-tpd-phase-discovery/)
retains numerical qualification on the stated neutral and admitted
density-dependent ionic/reacting cases, explicit density-independent
charged-search refusals and unchanged declared work. Method adoption is recorded in
[ADR 0012](../../../docs/scientific/adr/0012-native-tangent-plane-phase-discovery.md).

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

### Supports, seeds and trial solves

Let U be fed species plus all reaction participants and U₀ its neutral subset.
Initial discovery candidates are one free phase on U, Liquid, and also Vapor
when U is neutral. For reactions, first construct the existing positive
charge-neutral base on U, then run the unchanged ideal-solution speciation
sweeps. A zero-fed reaction participant remains on the required support; its
feed fraction is not an admissible log-chart seed. Record all initial attempts
and retain the accepted lower-G candidate. This preparation belongs only to
discovery; declared initialization remains unchanged.
Trial defaults are (U,Liquid), (U₀,Liquid), (U₀,Vapor), deduplicating identical
pairs. Omit empty default U₀ with an explicit no-neutral-species diagnostic.
Caller support/label entries **replace** defaults and limit the tested domain;
they do not restrict the initial union candidate. Reject empty caller lists,
empty supports, species outside U, one-sign ionic supports and ionic Vapor.

Before enumerating starts, one shared native support-admission predicate requires
`DensityTemperature` or `DipolarPerturbation` whenever a requested support
contains ions. These are the existing forms with $\epsilon_r\to1$ as
$\rho_m\to0$ (full EOS §6.2/§6.4). Density-independent `IonSuppression` and
`ComponentMixing` describe condensed-phase electrolyte terms and raise the typed
`ParameterError`: "permittivity model has no low-density limit; declare phases".
Reject the requested family explicitly; do not silently drop charged trials or
interpret refusal as stability. This predicate is used by the explicit test and
by discovery before its first declaration solve, and never by declared solves,
continuation, fitting or ASL. Neutral-only requested families are unaffected.
Admission reads the form, not a sampled density, fitted coefficient, phase label
or optimizer termination. It introduces no density cutoff or failed-trial filter.

For every admitted pair enumerate a positive charge-neutral base from equal
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

Run **every enumerated distinct start** through the shared driver. Compile an
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

Retain raw coordinates only from an actual nonfatal finite, domain-valid
optimizer attempt, including a nonconverged attempt. Failed initialization has
no optimizer point. Keep raw residuals separate from reconstructed composition,
displacement and reevaluated energy. Normalize neutral trials; ionic rescaling
and the largest-entry positive/negative pivot correction follow EqID:
`electroneutral_trial_reconstruction` with ties by index. Recompute molar volume
from reconstructed composition and the retained packing coordinate. Reject
lost positivity/support, invalid T/P, packing or association domain. Record
actual normalization/charge defects and rounding scales; ionic defects must
meet the 1e-12 criterion. Never fabricate a species floor to admit a point.

For candidate balance/phase-charge defects solve the local linear equality
correction in supported **amount** variables using existing native algebra.
For trials use normalization/charge rows. Reject rank loss or incompatible
corrections, enforce positive support/domain and reevaluate before descent.
Keep the nearby corrected candidate separate from the public result. Residuals,
correction norm and the energy-gradient contraction estimate feasibility error;
arithmetic equality roundoff remains explicit.

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
reconstruct arithmetic composition/charge, recompute density/packing and shared
energy, and verify matching conserved totals before comparing Gibbs energies.
Retain each failed check and count energy/density evaluations. Admit numerical
descent only when the finite Gibbs decrease exceeds its propagated estimate.
Optimizer convergence and positive trial inertia are diagnostics, not descent
requirements. An inaccessible negative plane remains unresolved; missing error
terms, near-zero margins or failed finite descent are indeterminate.

The adopted output is **numerically resolved descending split with estimated
uncertainty**. Calling it a **proved instability witness** requires justified
bounds on feasibility and every error term, which native association floors do
not provide. A completed finite search reports only **no witness in the recorded
trial set**, naming unresolved trials. No global stability conclusion follows.
Deduplicate physical trials using `PHASE_COINCIDENCE_DISTANCE` as a topology
resolution, not an energy-error threshold. Never discard actual resolved finite
descent solely because the trial is close to a present phase.

### Conserved topology changes and work limits

Rank distinct resolved descending trials by actual donor slope, more negative
first, ties by enumeration. This is deterministic, not a claim of the best
possible addition start. Seed the full added declaration with the admitted
finite split; preserve feed, reactions, standards and exact selected support.
Reuse shared density-branch admission; a unique unlabelled root retains the
source declaration label without asserting a physical phase theorem. Require
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
density on its declared branch, verify packing/association and balances/charges,
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
Composition directions preserve normalization, charge balance and exact support;
known fixed parameter slots project to zero, while unknown slots are input
errors. The [accepted P4 design](https://github.com/tannerpolley/ePC-SAFT/issues/19)
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
owner). The central Jacobian is equilibrated once to A = R J C with exact
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

A first-action direction may also supply one constant d(ln K)/dt per declared
reaction (empty means zero). The chemical-row series adds its minimum-norm
species image through the same standard-state owner, EqID:
`reaction_constant_coordinate_derivative`; residual-energy controls keep their
original width. A direction may also supply d(d ln K/dT)/dt in
`reaction_log_constant_temperature_slopes` (1/K; empty means zero). At fixed T,
these slopes contribute to referenced enthalpy actions through the same
minimum-norm species owner; they do not alter the chemical-row derivative.
Supplied standard Gibbs vectors fix K and reject nonzero reaction directions.
Nonzero reaction directions, including slope-only directions, refuse actions above
order one while retaining admitted central values. Reaction Hessians and
heat-capacity coefficient actions are not qualified.

## ALG-REFERENCE

Pressure-closed fugacity evaluates ln Z from the prescribed P and resolved
density; the returned EOS pressure still records the density-root residual.
This avoids amplifying the closure residual in parameter differences of nearly
invariant reference energies. Pressure-closed columns use EqID: `log_fugacity_coefficient_parameter_derivative`.
Each terminal reference column is admitted independently by the maximum scaled
two-ion change, max(|theta|,1) times the column difference, at the existing
reference tolerance. The terminal pair is compared first; a witness checks a
failed pair without changing the returned terminal column. An unavailable
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
order two; fixed-T/P feed-only order three remains available. A requested
temperature direction rebuilds every reference trace with E2 at its last three
points and uses EqID: `standard_state_transfer_temperature_derivative`: the
slope passes the value criterion in RT units, and its target slope
$s_r'=-d\ln K_{\rho,r}/dT$ follows the value's minimum-norm owner. Temperature
and mixed T/P actions are therefore available at order one; a refused slope
leaves the value and sibling directions intact. An order-two temperature
action adds $s_r''$ (EqID: `standard_state_transfer_temperature_second_derivative`)
and, when pressure also moves, $\partial^2 s_r/\partial T\partial P$ (EqID:
`standard_state_transfer_temperature_pressure_derivative`) to the row jet, so
temperature actions alone or mixed with feed or pressure are available at order two. A requested active EOS-parameter
direction uses the same three E2 evaluations and EqID:
`standard_state_transfer_parameter_derivative`: each requested slope passes the
value criterion at the parameter scale max(|theta|, 1), and
$\partial s_r/\partial\theta=-dC_r/d\theta$ follows the same owner, so parameter
and mixed T/P/parameter actions are available at order one. Continuation
tangents use the same order-one chains. Free temperature, raw lever derivatives, temperature actions above order two,
parameter actions above order one, heat-capacity parameter or reaction-coefficient
actions, and third-order actions involving pressure remain unavailable.
Ideal/Total enthalpy and heat capacity follow EqID: `reference_species_calorics`:
supplied records anchor their species, record-less species are completed from
the EOS-basis reaction enthalpies and heat capacities, and only a charge gauge
may remain. A heat-capacity value or enthalpy temperature action requests the
second slope (EqID: `standard_state_transfer_temperature_second_derivative`),
admitted when its terminal change is at most 0.01 J/(mol K) in heat-capacity
units; its refusal leaves enthalpy values and feed actions available. Enthalpy values
and first feed, temperature and pressure actions (a driven pressure, or the
solved pressure of a free-pressure cell, through the mixed slope's
$\partial\Delta_rh_\rho/\partial P$) and heat-capacity values and first feed
actions at fixed pressure are available; heat-capacity temperature or pressure
actions, caloric orders above one, u, s and c_v remain unavailable, while independently
valid residual properties remain usable.

For a referenced enthalpy first action at fixed temperature, request only the
consumed EOS columns of $C_{T\theta}$ (EqID:
`standard_state_transfer_temperature_parameter_derivative`) and reaction
coefficient columns of $\phi'_k(T)$ (EqID:
`reaction_constant_coordinate_derivative`). The EOS mixed derivative uses
balanced polarization of repeated-direction E3 actions and is admitted when
$T\max(|\theta|,1)$ times the observed change over the last three reference
ladder points is at most $5\times10^{-5}$. A missing column refuses only an
action that consumes it and names its reaction and parameter; central values
and sibling directions remain available. The reference enthalpy action uses
the fixed-T endpoint's own solved temperature and reaction-row mapping. It does
not admit heat-capacity parameter actions or referenced free-temperature
actions.

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

Fits use exact Engine derivatives and run natively. Python declares specifications
and data and reads results; finite differences are numerical verification only.
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
The implemented fit (R1--R3) includes pressure-closed density and total molar
isobaric heat capacity with fixed thermochemistry, finite-addition enthalpy,
the declared equilibrium observations and direct source-basis ln K observations. It uses affine, tied
optimizer coordinates of the formulation, a Ceres Levenberg--Marquardt trust
region with dense QR on one thread, the named `FIT_*` policy tolerances as the
only defaults, exact owner actions for every Jacobian and a covariance that is
withheld under rank deficiency or an active bound. Named scalar ln K and stored
`ReactionLogPolynomial` a/b/c/d coordinates share the affine chart with EOS
coordinates. Every matching observation must declare the same reaction,
coefficients and source basis. Rebuild each trial with those values, preserving
reference transfers. Optimizer columns follow declaration order; physical rows
are EOS slots first, then reaction coordinates. Both direction families use one
solved-state first-action batch and one defining-Jacobian factorization per state.
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
and source-reference declarations. Referenced free-temperature endpoints are
unavailable. Failures retain the `initial:`, `final:` or `feed:` state name
with the owning solve, density, thermochemistry or action diagnostic.
Affinity targets, supplied Gibbs vectors, dependent reaction sets, constants
fixed by supplied thermochemistry and referenced free-temperature requests are
rejected by reaction name before fitting.

Keep admission's model and compiled state groups for the initial evaluation.
Reuse them only when the bounded affine map reproduces every active EOS value
and the rewritten reaction declarations exactly; otherwise rebuild at the
resolved optimizer point. The physical-declaration comparison includes the
scalar K or all stored correlation coefficients, so a reaction-coordinate
change invalidates admission preparation even when EOS values are unchanged. The same
compiled problem may serve a cold retry at that point. Failed trials still
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
scale; zero thresholds remain zero. The complete map includes all tied affine
scales and eliminated Born rows. The epigraph and slack columns are analytic.
Only exact first actions are supplied; Ipopt uses limited-memory BFGS, with no
finite-difference derivative checker or parameter Hessian claim. Equilibrium
retains exact Lagrangian curvature. The common checked Ipopt option/lifecycle
owner preserves equilibrium's `Initialize()` and uses `Initialize("")` for
estimation, excluding file overrides. Required positive scaled tolerance,
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
and achieved diameter uncertainty/correction; refinement qualification belongs
to the retained study.

Termination, point/action availability, feasibility, final checks and
stationarity are distinct. Final checks evaluate report-only rows omitted from
the minimax trial; their raw Jacobian rows carry no evaluated columns.
Final checks use the stored fitted diameters,
without resolving holds. A fatal initialization has no candidate; local
infeasibility proves no global infeasibility. No covariance, unique optimum,
global certificate or physical validation follows from this local minimax
estimate. `evaluate_minimax` evaluates the same compiled rows at a declared
point without optimizing or publishing a stationarity certificate.

## ALG-TERMINATION

Core CORE1 keeps the numerical domains separate. `EQUILIBRIUM_COMPOSITION_BOX`
$\beta$ sets the lower box of the intensive `log(x_i)` coordinate inside the
optimizer subproblem: its barrier regularises the log chart's trace steps and
is not a solution domain. Species $i$ has the box

$$
x_i \ge \max\!\left(x_{\mathrm{cut}},\ \min\!\left(\beta,\ \varepsilon\,\frac{u_i}{N_{\mathrm{feed}}}\right)\right),
\qquad
u_i = \max\,\{\, n_i : B n = b,\ n \ge 0,\ n_j = 0 \text{ off every support} \,\},
$$

where $x_{\mathrm{cut}}$ is `STABILITY_COMPOSITION_FLOOR`, $\varepsilon$ is
`EQUILIBRIUM_BOX_INVENTORY_FRACTION` and $u_i$ is the largest attainable total amount of $i$,
one HiGHS linear program per species. Rows of $B$ can carry mixed signs (a hydrogen row
$n_{\mathrm{H^+}}-n_{\mathrm{LiA}}-n_{\mathrm{NaA}}$), so no single row ratio bounds $n_i$; an
unbounded or failed program keeps $\beta$. A uniform box $\beta$ is a solution domain whenever
$\sum_p \beta n_p > u_i$: it forces more of $i$ into the phases than exists and the subproblem
is infeasible, and just above that the optimum sits on the box beyond the polish's reach
(#172: 61 mol aqueous $\times$ 1e-8 = 6.1e-7 mol against 3e-7 mol Li fed; stalls up to 3e-6 mol).

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
locally infeasible at a Na-loaded state that the uniform box solves in 9. The cutoff keeps
the bound finite for an unattainable species ($u_i = 0$). The box does not set accuracy: a
moiety row $(Bn-b)/N_{\mathrm{feed}}$ resolves an inventory to its achieved residual, at most the
caller's row tolerance and after the polish near rounding, $\epsilon N_{\mathrm{feed}}$ (trace
isomers fed at 1e-10, 1e-12 and 1e-14 mol in 1 mol return their amount to 4e-8, 5e-9 and 9e-3
relative; the orthonormal basis gave 8e-6, 3e-4 and a factor 4). Square pinned equation systems
use the structural cutoff for every species, so chemical equalities can reach trace roots.
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
