# #90 initialization and continuation diagnosis — 2026-09-21

**Frozen record — not rerunnable as a current-main calculation.** Record snapshot:
commit `ca719012c31e8b41edcb98ae82a388117b58d343`. Original wheel identities, source
inputs and retained numerical evidence keep their recorded provenance. The two producers
require the original hash-verified wheel and packet identities; they are excluded from
current-main replay by #696.

This is the bounded Diagnosis/Design handoff for issue #90. It tests the
existing single equilibrium driver through the qualified non-editable wheel.
It does not change the solver, add a retry path, weaken a tolerance, or claim
that a failed initial guess represents physical nonexistence.

## Reproduction and owners

- Engine source: `f909bb21a93171c38465adba8777cd3b0bc6c53d` on
  `roadmap/application-preparation`.
- Installed wheel:
  `build/environment-wheel/epcsaft-0.2.0.dev0-cp313-cp313-linux_x86_64.whl`.
- Required wheel SHA-256:
  `9f73aeb466c54eff80c89a8776c90d49cb9c9c2d130f0eb15bf0311c1f76e606`.
- Neutral parameter packet:
  `data/packets/gross-2001-methane-ethane/2/parameters/parameters.json`.
- Parameter packet SHA-256:
  `7dd1899712aacfbc3395d86349768b0d93dd5dfd75b2e0f0700640be618aa9c2`.
- Executable owner: `engine/native/equilibrium/gibbs_minimizer.cpp`.
- Continuation owner: `engine/native/equilibrium/phase_boundary.cpp`.
- Public route: `engine/src/epcsaft/equilibrium/problem.py`.

The primary reproducer is the retained R1-F5 neutral methane/ethane boundary:
feed amounts `(0.5, 0.5) mol`, declared phase order `L,V` with a pinned-zero
vapor, `T*=180 K`, and `P*=1596859.1113026983 Pa`. The free-temperature problem
holds `P*` fixed and supplies a cold temperature guess. One accepted fixed-`T`
anchor at 180 K is retained in the manifest and in the comparison columns. A
target result is assigned the same local branch only when its phase-ordered
amounts, compositions and densities compare with that anchor under the existing
E4 criteria. The R1-F5 driver row establishes the above-start result from
`190 K` and the typed below-start topology event from `170 K`.

The secondary control is the existing synthetic reactive fixture from
`engine/tests/equilibrium/test_native_transport.py`: two liquid phases,
`A+B=C`, feed `(0.05, 0.475, 0.475) mol`, `P=100000 Pa`, and the accepted
EOS-basis correlation `ln K(T)=1500(1/T-1/298.15)`, evaluated at `298.15 K`.
This is a numerical initialization control, not application chemistry or a
physical prediction.

## Frozen matrix and criteria

The matrix is declared as constants in `scripts/run.py` before retained rows
are written. It contains:

| Family | Frozen controls | Expected interpretation |
|---|---|---|
| Free-temperature cold | guesses `170, 179, 180, 181, 190 K` at `P*` | An accepted state must return `T=180 K`; a typed topology event preserves a numerical basin/branch limitation and is not a physical nonexistence result |
| Warm free-temperature | accepted fixed-temperature anchors at `179` and `181 K` into the `T`-free, `P*`-fixed problem | Recover the same target branch with the same row criterion |
| Boundary continuation | `179→180 K` and reverse `181→180 K`, `initial_step=minimum_step=0.25`, `max_steps=4`, `max_retries_per_step=1` configured but unreachable at this minimum step; no retry performed | Reach the target with an accepted prefix and no topology event |
| Long continuation control | `179→181 K` and reverse `181→179 K` with the same step and `max_steps=8` | Preserve the existing bounded branch-window loss and its accepted prefix; this is a falsifier against silently skipping branch admission |
| Reactive fixed state | `T=298.15 K`, `P=100000 Pa`, the synthetic two-liquid reaction case | Accepted state with reaction and material rows under the existing driver |
| Reactive warm/continuation | warm solve and `298.15→299.15 K` continuation, step `0.5 K`, two steps | Accepted states on the same declared reactive branch |

Units are K, Pa, mol, and dimensionless residuals. The frozen caller row
tolerance is `1e-10`; accepted rows require `success=true`,
`requested_tolerance_met=true`, no topology event, and convergence residual
`<=1e-10`. For the canonical free-temperature state, the physical round-trip
criteria are `|T-180 K| <= 1e-6 K` and relative pressure error `<=1e-8`.
The anchor comparison retains phase order, amounts, molar densities and mole
fractions. It uses E4's `1e-8` absolute composition and amount criteria, the
existing `1e-8` relative pressure criterion, and the existing P2 density
envelope `|log(rho/rho_anchor)| <= 0.009950330853168083`. Intermediate
continuation states retain these arrays and are checked against the anchor only
at the 180 K endpoint; their accepted branch identity is governed by the
existing continuation density/topology checks. The analysis records both the
driver's convergence residual and the maximum absolute raw/classified residual
row. Solver iterations and evaluations are the work measures; no wall-time
claim is used as acceptance.

## Retained evidence

- [`scripts/run.py`](scripts/run.py) is the reproducible read-only driver.
- [`results/cases.csv`](results/cases.csv) retains one row per frozen direct or
  continuation case, including status, phase-ordered state arrays, anchor
  comparison, residuals, work and failure meaning.
- [`results/continuation_states.csv`](results/continuation_states.csv) retains
  each accepted continuation state and path value with its phase-ordered
  amounts, densities, compositions and endpoint comparison when applicable.
- [`manifest.json`](results/manifest.json) records the wheel, packet, fixture hashes,
  frozen controls, the fixed-180 anchor and run metadata.

The script verifies the wheel and packet hashes before importing the model and
never builds or installs anything. It imports the existing synthetic fixture
only to avoid creating a second parameter owner. No application input is
invented: MEA, absorber and Lithium cases remain pending their independent
input contracts.

## Diagnosis and proposed design boundary

The primary known cause is initialization/branch admission sensitivity in the
declared two-phase free-temperature route. The driver can converge to a
stationary state that is then rejected when the two phases coalesce; the result
retains the topology event and finite convergence record. Accepted warm states
compare with the fixed-180 anchor in phase order, and bounded continuation
recovers that endpoint while each predictor step remains inside the existing
local branch window. A longer path may lose the anchored branch; that is a
legitimate continuation limitation, not evidence that the target state does not
exist.

The supported no-code workflow is already present in the public driver: a caller
may supply an accepted compatible result to `solve_equilibrium`, and may supply
an accepted `start_result` to `continue_equilibrium`. The declaration, phase
order, support, amount kind, reaction laws and feed layout must remain
compatible. This preserves the existing topology and branch-loss semantics and
does not promise universal convergence.

A generic cold-start fallback/search remains unresolved. The retained matrix
does not establish a seed selection rule, an automatic continuation chart or a
caller-facing exhaustion contract, so no kernel design or invocation budget is
proposed. Any future capability must preserve typed topology/branch failures,
declared phases, the existing physical criteria and the single-driver route.

The matrix is sufficient to falsify the bounded design if a compatible warm
start or quarter-kelvin continuation fails on the accepted branch, if a typed
failure loses its diagnostics, or if a recovered state misses the existing
row/branch criteria. It does not qualify MEA, the absorber, Lithium, density
closure, phase discovery, or physical prediction error.

## Retained result

The installed-wheel run retained 14 case rows and 15 continuation states.

- Cold free-temperature starts at `170`, `179` and `180 K` converged in the
  optimizer but were rejected with `PhasesCoalesced`; their finite convergence
  residuals were `3.47e-13`, `3.08e-12` and `2.66e-15`, respectively. These
  are typed branch/topology failures, not physical nonexistence.
- Cold starts at `181` and `190 K` recovered `T=180 K` at the fixed pressure,
  with residuals `2.66e-15` and `2.98e-14`, respectively. The result reproduces
  the canonical R1-F5 above/below asymmetry and shows that the basin boundary
  is close to the target.
- Warm starts from accepted fixed-temperature states at `179` and `181 K`
  both recovered the target in 3 iterations and 4 evaluations, with residuals
  `1.39e-13` and `1.70e-13`. Their phase-ordered anchor errors were at most
  `7.36e-15` in composition, `6.73e-14` in density log error and zero in
  amount error.
- Forward `179→180 K` and reverse `181→180 K` continuation both reached the
  target in four accepted quarter-kelvin steps, with all state residuals below
  `4.4e-14` and `3.2e-14`, respectively. The longer `179→181 K` and
  `181→179 K` controls both stopped at the first attempted half-kelvin
  correction (`q=0.25`) with `branch_lost`/`step_too_small`; each retained its accepted
  start and the attempted-state density diagnostic. The endpoint anchor errors
  were at most `3.40e-16` in composition, `4.42e-14` in pressure and
  `4.32e-14` in density log error.
- The synthetic reactive fixed state, its warm `298.15→299.15 K` solve and
  its two-step continuation all reached accepted states with classified
  residuals between `1.69e-12` and `1.75e-12`. The largest unscaled raw row
  was about `6.7e-9`, which remains a diagnostic scaling distinction owned by
  #87 rather than a reason to alter the #90 convergence criterion.

The reproduced causal boundary is therefore initialization and local branch
admission: cold starts can reach a stationary coalesced solution, while an
accepted compatible state or sufficiently short continuation recovers the
declared branch. The evidence supports the existing caller-supplied warm-start
route as a usable workflow. It does not establish a generic automatic fallback,
fixed seed schedule or initialization chart.

The continuation coordinate is dimensionless. For `179→180 K`, a `q` increment
of `0.25` is a `0.25 K` temperature increment; for `179→181 K`, the same `q`
increment is `0.5 K`. In the retained long paths `minimum_step` equals
`initial_step`, so the single allowed retry cannot halve the step; the observed
`branch_lost`/`step_too_small` result is therefore the declared minimum-step
termination, not evidence of a tested smaller-step recovery.

The separate #87 trace/action handoff remains distinct. Its `K=1e-8` central
states retain `Ipopt SUCCESS` followed by `Stop.stalled`, classified residuals
about `2.980519e-10` and `1.093328e-10`, and raw trace-B chemical rows about
`0.02895457` and `0.01081567`; cleared states provide no claimed actions. That
diagnosis concerns trace/action availability and scaling, not this
initialization/branch result.

## Design boundary and falsifiers

There is no implementation-ready generic fallback in this correction. The
smallest remaining design question is whether current consumers need a
caller-visible contract for supplying an accepted compatible anchor, or whether
the existing `start`/`start_result` arguments and their diagnostics suffice.
That question can be answered from the settled application input contracts and
does not authorize a cold scan, a second solver, phase discovery, a topology
switch, a tolerance change or a fabricated stationary result.

The retained evidence falsifies the supported warm/continuation workflow if a
phase-ordered target result misses the fixed-180 anchor criteria, if a typed
failure loses its diagnostics, or if a recovered state changes branch identity.
It is not qualified here for MEA, the absorber or Lithium; their exact state inputs
belong to the application repositories (#82, #83 and #91 are closed). The synthetic reactive state continues to
serve only as a numerical initialization control.

## Review and correction history

The initial handoff was reviewed as a supported diagnosis with a revised design
required: branch identity needed state-vector evidence, the continuation retry
meaning needed correction, and the generic fallback contract could not use
known-answer anchors. Correction round 1 adds the fixed-180 anchor snapshot and
phase-ordered comparisons, corrects the dimensionless-`q` and minimum-step
interpretation, removes the unsupported invocation ceiling, and records the
existing caller-supplied warm-start route. A fresh independent execution review
of this final correction remains pending.

## Staged-start result — 2026-09-21

The frozen follow-up question was executed by
[`scripts/staged_start.py`](scripts/staged_start.py) against the same installed
wheel. For each caller temperature guess, the runner made one fixed-`T`,
free-`P` anchor call using that same guess and the generic `100000 Pa` pressure
guess. Every accepted anchor was passed to exactly one original free-`T`,
fixed-`P` call at `P*=1596859.1113026983 Pa`; no retry, interpolation,
continuation or known-answer seed was used. The run made all ten permitted
driver calls in `0.231 s` and retained the five rows in
[`results/staged-start.csv`](results/staged-start.csv), with the compact
provenance record in
[`results/staged-start-manifest.json`](results/staged-start-manifest.json).

| Caller guess (K) | Anchor P (Pa) | Target iterations/evaluations | Target residual | Branch match |
|---:|---:|---:|---:|:---:|
| 170 | 1175945.1075064444 | 5 / 6 | 5.33e-15 | true |
| 179 | 1551256.2800381570 | 3 / 4 | 1.39e-13 | true |
| 180 | 1596859.1113026983 | 0 / 1 | 1.94e-13 | true |
| 181 | 1643244.8286910641 | 3 / 4 | 1.70e-13 | true |
| 190 | 2095823.2840849985 | 5 / 6 | 1.24e-14 | true |

All five target rows were accepted with `topology_event=None_`,
`requested_tolerance_met=true`, `T=180 K` within `5.46e-12 K`, and pressure
equal to the target within the recorded precision. Compared with the retained
fixed-180 K anchor, the largest target composition error was `7.36e-15`, the
largest density log error was `2.84e-13`, the largest amount error was zero,
and the pressure error was zero. The anchor pressure varies with the caller
temperature as expected; it is only an intermediate accepted state and is not
the target answer.

For this methane/ethane case family, the result supports a bounded two-call
workflow as a design candidate: the existing initialization path makes one
caller-temperature fixed-`T`, free-`P` anchor, then the existing single driver
solves the original free-`T`, fixed-`P` problem using that accepted result as
`start`. This is a diagnosis/design candidate only. It establishes no
general convergence guarantee and authorizes no kernel change, continuation
search or automatic retry. Nonunit feed scales, application input contracts
and recovered-state action qualification remain untested.

Independent integrated result review supports this bounded result. The only
staged-result correction changed the documented target-solve argument from
`start_result` to `start`; numerical evidence was unchanged. This review does
not extend the result beyond the five frozen starts or establish Build readiness.
