# Trace composition and action diagnosis

**Frozen record — not rerunnable as a current-main calculation.** Record snapshot:
commit `ca719012c31e8b41edcb98ae82a388117b58d343`. Original wheel identities, source
inputs and retained numerical evidence keep their recorded provenance. The original
probes pin a pre-repair wheel and do not apply to other wheels; #341 removed them, with
their method retained at record commit `e3b61eb9`.

Question (#87): does the amount-weighted convergence classification provide
useful intensive composition and small-action accuracy at a trace participant?
This is a bounded numerical diagnosis, not application qualification or authority
to alter equilibrium equations, charts, tolerances or state-availability policy.

The reference is the existing E5-3 identical-ethane isomerisation A ⇌ B at
350 K and 200000 Pa. Equal residual chemical potentials imply
`x_B = K/(1+K)` exactly. The retained temperature law has
`L = ln K = a + 1200/T + 0.4 ln(T/300) + 0.002 T`, with `a` set so
`K(350)` is respectively 2, 1e-4 and 1e-8. Inputs are dimensionless mole
fractions, K, Pa and mol; actions use unit kelvin directions. No thermal
property record is needed for this reaction-law composition identity.

Before execution, the controls are frozen at those three K values, each with
feed amounts `(N/2,N/2)` for `N=1,10 mol`, one declared vapor phase and
`row_tolerance=1e-10`. The six solves use the existing single driver and
non-editable wheel. The probe stops after this matrix (60 s wall-clock cap);
it does not retry with relaxed criteria or search other starts.

Independent reference derivatives follow the logistic chain. Writing `p=x_B`,
`L1=dL/dT`, `L2=d²L/dT²`, `L3=d³L/dT³`, they are
`p1=p(1-p)L1`, `p2=p(1-p)[L2+(1-2p)L1²]` and
`p3=p(1-p)[L3+3(1-2p)L1 L2+(1-6p+6p²)L1³]`.
Every reference is nonzero. The diagnostic screens are relative composition
error ≤1e-6, and relative action error ≤5e-6,5e-5,5e-4 for orders 1,2,3.
The latter retain E5's order-dependent error allowances but remove its
absolute floor for these analytically known nonzero actions: a dropped action
must fail. These screens investigate numerical capability; they are not MEA,
Lithium or film physical error limits. A failed solve/action stays a failed
case, not evidence of physical nonexistence.

For normalized material constraints `(B n-b)/N`, the extensive objective
multiplier is `lambda_constraint = -N lambda_balance`. Compare the first
material multiplier using a relative 1e-8 screen (and inspect its nonzero
magnitude). Record chemical and classified residuals, achieved stopping cause,
success, whether the requested tolerance was met, and returned quantities together.

The probe, `scripts/probe.py` at record commit `e3b61eb9`, ran on the installed
wheel `9f73aeb4…` built from source `f909bb21a93171c38465adba8777cd3b0bc6c53d`
and reused the input builder in `engine/tests/native_inputs.py`; it refuses
any other wheel, so #341 removed it from the tree. The probe records the wheel,
module, helper and probe identities in a compact provenance file; numerical rows are
retained once in CSV. Kernel source is unchanged. Fresh independent Sol CSE
review supports this bounded result: high-precision reference differentiation,
byte-identical replay of the 22 rows and installed-wheel file verification
passed. No correction round was needed; #87 owns the review disposition.

## Retained result and source diagnosis

The six controls completed in under 0.1 s on the identified wheel. The 22 rows in
[`results/trace-controls.csv`](results/trace-controls.csv) include both rejected
central states; no missing actions at a rejected state are counted as zeros.

| K at 350 K | Central outcome for 1 and 10 mol | Numerical information |
|---|---|---|
| 2 | accepted, request met | x_B relative error ≈3.2e-13; A3 relative error ≤2.7e-11 |
| 1e-4 | accepted, request met | x_B and A1–A3 relative errors ≤1.96e-8; nonzero A3 ≈−1.807e-10 K^-3 passes without an absolute floor |
| 1e-8 | both rejected despite Ipopt SUCCESS | classified residuals 2.98e-10 and 1.09e-10 exceed the 1e-10 request; unscaled maxima 0.02895 and 0.01082; stop is stalled and state arrays are cleared |

All four accepted material-dual identities satisfy the 1e-8 relative screen
(largest measured defect 1.45e-12). This supports the existing multiplier
transport on these normalized constraints at nonunit feed, not arbitrary
application scaling or all constraint kinds.

The 1e-8 case has the explicit analytic solution
`x_B=9.999999900000002e-9`; numerical rejection cannot establish its physical
absence. The current rejection correctly preserves the requested criterion.
The control is handed to #90 for initialization/recovery diagnosis, with #87
retaining ownership of intensive accuracy and result meaning. No correction
to the optimizer or acceptance chart follows from this experiment alone.

Current-source inspection confirms the following boundaries:

- `request_compile.cpp:797` classifies free chemical rows after multiplication
  by `N_phase*x_i/N_feed`. Therefore classified tolerance `tau` only bounds
  an individual intensive row by `tau*N_feed/(N_phase*x_i)`. Pinned rows are
  unscaled. This is an accepted chart choice, not automatically a defect.
- `gibbs_minimizer.cpp:768` accepts tiny-step termination only with
  `Stop::converged`, while solver SUCCESS may accompany an accepted
  `converged_to_floor` state. B1–B4 retained evidence already exercises the
  rejected tiny-step/floor branch; it does not exercise the accepted
  solver-success/floor branch. Neither this experiment nor that historical
  test changes the accepted availability rule.
- `implicit_sensitivity.cpp:574` admits actions using the classified
  residual, then solves unscaled defining-row coefficient equations. A
  classified central gate alone is not a relative trace-observable bound.
- The same action owner preserves typed unavailable Properties rows, mapped
  values and valid siblings for supported dependency failures. Its series
  inverse/logarithm throw domain errors, and there is no local containment
  around every row/path evaluation. This is a source-confirmed propagation
  risk; an admissible current-wheel triggering case is still needed before
  specifying a repair. No singular action is to be regularized or set to zero.

These controls support no code change to the already accurate mild-trace
actions or multiplier mapping. Remaining #87 work is consumer-specific
intensive error assessment, a reachable singular/domain partial-result case,
and precise documentation of the two existing floor-limited outcomes. MEA,
Lithium and film contracts remain distinct inputs; missing chemistry choices
belong to those application owners. The reviewed diagnosis establishes no
implementation-ready kernel change.

## Three-request floor-result follow-up

After the six-case result was reviewed, #87 froze three tolerances on the same
K=2, 350 K, 200000 Pa, 1 mol control: `1e-12`, `1e-14`, `1e-16`, one cold
solve each and a 10 s total cap. The question was whether this existing smooth
state reaches the accepted solver-success/unmet-tolerance branch. The probe,
`scripts/floor_probe.py` at `e3b61eb9`, ran on the same pinned wheel.

The three rows in [floor-controls.csv](results/floor-controls.csv), with
[their provenance file](results/floor-manifest.json), completed in under 0.02 s. The
`1e-12` and `1e-14` requests succeed and meet their requests; composition and
A1 relative errors are below `3.4e-15`. The `1e-16` request stops at a tiny
step and is rejected (`Stop.stalled`, residual `4.09635e-16`); central state
and A1 are unavailable. All three report no numerical-floor estimate.

This probe **does not reach the accepted floor-limited branch**. It supports
the retained rejection/availability meaning at the tight request, and does
not justify weakening that request, extrapolating a floor estimate or another
status search. The earlier six-case numerical results are unchanged. Fresh
independent Sol integrated review supports this follow-up with no material
findings and no correction rounds: the reviewer checked the analytic
composition/A1 at high precision, availability meaning against source, and
all 14 installed wheel payload files. No additional solver replay was needed.

## Repair outcome — 2026-09-21

Cause of the K = 1e-8 rejections: the analytic solution `x_B = 9.9999999e-9`
lies below the driver's `ln(x_i)` box at `ln(1e-8)`. Ipopt satisfies its own
KKT conditions with the bound multiplier `z = μ/s` absorbing the stationarity
residual, so it reports SUCCESS, while the Engine's classified row carries no
bound multiplier and reads `x_B·s`; with the intensive residual `r ≈ s`, the
retained classified values follow `x_B s² = μ` (μ ≈ 8.4e-12 and 1.2e-12 for the
two feeds). It was not an initialization or scaling problem. A second, chart
limit then appeared once the box moved away: the optimizer stops on
amount-weighted rows, so an accepted iterate bounds an intensive trace row only
by `τ N_feed/(N_p x_i)`; the K = 1e-8 state was accepted with a 4.4e-5 relative
error in `x_B`, K = 1e-12 with 350%, both below 1e-12 classified.

Repair in the driver owner (`engine/native/equilibrium/gibbs_minimizer.cpp`,
ALG-TERMINATION): the optimizer box stays at `EQUILIBRIUM_COMPOSITION_BOX =
1e-8`, because its barrier is the step regulariser of the log chart (moving it
to 1e-20 let the line search take steps of order 1e4 in `ln x`: the E5-3
corner solves coalesced and the S6 single-phase solve aborted), while the
structural cutoff `STABILITY_COMPOSITION_FLOOR = 1e-20` governs continuation
admission, action differentiability and the polish domain. After an Ipopt
success at most `EQUILIBRIUM_POLISH_STEPS = 3` Newton steps on the unweighted
rows with the exact row Jacobian close the intensive gap and carry a trace
solution below the box (the trace row is nearly linear in `ln x`, so one step
spans decades); a step stands only when it contracts both the raw row maximum
and the classified residual and keeps every composition above the cutoff, and
the polished point is classified on the caller criterion exactly as before. No tolerance, weight, objective, initialization or acceptance rule
changed; tiny-step stops are not polished.

Persistent falsifier: `engine/tests/equilibrium/test_trace_reliability.py`
(K ∈ {1e-4, 1e-8, 1e-12, 1e-16} × N ∈ {1, 10} mol, the same law, screens and logistic
references as this analysis). On the repaired wheel identified in the pull
request, all eight solves are accepted with the 1e-10 request met, classified and
unscaled residuals about 1e-16, `x_B` relative error ≤ 1.4e-14 (before: 1.96e-8
at K = 1e-4; rejected at 1e-8), A1–A3 within the retained 5e-6/5e-5/5e-4
screens and the independent extensive-objective KKT residual to 1e-12. The falsifier reaches
x ≈ 1e-16 for composition and A1–A3. Continuation admission at these trace
levels and action availability between 1e-16 and the 1e-20 cutoff remain
unexercised qualification limits, as are
pinned declared phases (their rows are optimizer equality constraints, so a
solution below the box is infeasible there and not polished), warm starts
below the box (pulled back into it by the optimizer's bound push and
re-recovered by the polish), and the polish rejection branch. The
returned extensive-objective constraint multipliers are reconstructed at the
returned point from the row multipliers (`-N_feed lambda` for material rows,
`-n_p phi` for charge rows, zero for a free phase's normalization row), and the
reduced-Hessian certificate uses them there; the owner review of the pull
request had found the certificate combining the polished primal point with the
optimizer's pre-polish multipliers (material-dual defect 1.1e-8 to 1.4e-8
relative at K = 1e-14 to 2e-20). The persistent falsifier now also checks the
balance multiplier against the species chemical potential over RT formed from
public state quantities, and the certificate against a fresh solve started at
the returned point, down to K = 1e-16. The 22 pre-repair rows above are
retained as the defect record; the probe pinned the pre-repair wheel.

The same persistent test now qualifies the charge-constraint multiplier on
Khudaida Figure 2 point 1, packet v6, at 293.15 K and 100000 Pa, with 1 and
10 mol feeds. Independent objective-gradient/Jacobian stationarity defects are
1.34e-14 and 6.76e-14 (screen 1e-12). Public species chemical potentials
recover the nonzero electrochemical offset 0.6456098173 and neutral-salt
transfer equality to 1e-12; each liquid is electroneutral. Substituting phase
fraction for phase amount at 10 mol fails stationarity by 6.26e-4. This closes
that ionic scaling falsifier on the stated fixed state; application-specific
ionic trace and action qualification remains open. Existing association
physics is consumed unchanged, and no full LLE campaign was run.

Partial-result propagation: the `logarithm` throw for a `PhaseLogActivity`
centre is unreachable by an admissible request, because the observable
validation rejects a log activity of a species outside the phase support with
a typed `invalid_argument`, and every supported composition is positive in the
log chart. No guard or falsifier was added; the risk is recorded as a stated
limit, not a demonstrated defect.
