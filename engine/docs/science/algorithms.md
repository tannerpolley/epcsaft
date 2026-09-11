# Scientific algorithms

Accepted requirements, reconciled against the 2026-09-06 transfer audit.
These methods consume the [equation definitions](../equations.md),
[EOS conventions](full-eos.md) and [equilibrium formulation](../equilibrium/formulation.md).
Algorithm IDs below are stable references. Current implementation and required
numerical falsifiers are mapped in [implementation plan](../../../docs/scientific/implementation-plan.md).
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
verify derivative linear residual and conditioning
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

```text
scan deterministically in log density within the admitted packing domain
bracket sign-changing pressure residuals
refine bracketed candidates with safeguarded steps
retain finite roots meeting pressure residual and positive-slope criteria
select only by declared root/branch policy; return root count and branch identity
```

Use gradient-only energy work during the scan when sufficient; evaluate the
correct reduced pressure slope for refinement and conditioning. Tangent roots
at a spinodal need separate treatment and are not stable nonsingular roots.
An unlabelled request with multiple stable roots is ambiguous. Tracking a
branch through a continuation step must not silently relabel a remaining root
when the prior branch disappears. Search completeness is bounded by the
specified scan; mechanical stability does not establish global stability.

## ALG-EQUILIBRIUM

All admitted specifications use the one equilibrium driver required by the
[project map](../../../docs/scientific/CONTEXT.md). Compile balances, supports,
pins and free intensives from data. Check structural
counts/rank, construct declared feasible starts, then solve the same intensive
rows using the required common energy. Recompute achieved residuals and phase
distinctness at the returned point; optimizer termination alone is insufficient.
Retain the native optimizer `solver_status` alongside the convergence stop,
residual, and failed-attempt initialization/termination details. The Ipopt driver
scales its extensive Gibbs objective by `1 / feed_total` internally; this
preserves the minimizer while keeping objective derivatives on the same
numerical scale as the normalized material-balance rows. It accepts Ipopt's
tiny-step termination only when the recomputed residual meets the declared
equilibrium criterion; a tiny-step status alone is never success. Ipopt unscales
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

Continuation requires a predictor/corrector, step admission, branch checks and
explicit stopping behavior. Reusing only an initial seed is not continuation.
Controller and step policy remain to be qualified on the retained boundary
cases; no new numerical constants are selected here.

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

The conservative absorber's first exact residual Jacobian already needs A2
equilibrium and caloric actions. A limited-memory outer Hessian does not remove
that need. An exact outer Hessian may require A3. The older ordinary seven-state
absorber route only needing A1 does not define the new application's requirement.

## ALG-REFERENCE

Use the explicitly declared solvent/reference path and source standard state.
Evaluate the retained molality sequence; return the terminal change and
criterion separately. Unknown validity or failed limiting convergence remains
unavailable. Differentiate both the source correlation and EOS reference path,
including its density closure. Fixed source pressure P0 must not become trial
pressure. Reference-gauge choices that preserve reaction affinities do not
supply absolute caloric data.

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

## ALG-TERMINATION

Core CORE1 keeps the numerical domains separate: `STABILITY_COMPOSITION_FLOOR`
is only the lower bound for the intensive `log(x_i)` coordinate. Free-phase
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
`STOP_AT_TINY_STEP` with a finite non-converged residual is stalled and
unavailable, while a finite criterion-meeting tiny-step result remains
eligible. A tiny-step status alone never establishes success; recomputed
residual and topology criteria still apply.
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
result retains the same convergence record alongside its state.

## Verification obligations

Before a scientific behavior change is accepted, update or explicitly revalidate
the affected equation/algorithm IDs and execute the smallest independent
falsifier identified in the transfer plan. Important examples are caloric
identities, association support boundaries, full weighted-Hessian directional
checks with nonzero multipliers, implicit action ladders and independent
reference-state transformations. Finite checks establish bounded evidence,
not universal symbolic correctness. Documentation readiness checks cannot
replace them. Expensive application/reproduction sweeps remain explicit.
