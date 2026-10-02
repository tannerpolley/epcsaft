# Engine code budgets

The owner requested a hypothetical floor as the design goal, with soft and hard
ceilings to constrain duplication and unnecessary machinery. Here **floor
target** means an engineering estimate of a lean complete implementation. It is
not a proven minimum, an allocation to spend or a lower acceptance bound.

The single numeric table is `ARCHITECTURE.yaml:code_budgets`. Run
`python3 tools/diffstat.py --budgets` to display current counts against it.
Both normal local check commands enforce it. `check-affected.sh` provides the
documentation-only route without a wheel build; `check-all.sh` always builds
and installs after readiness checks. Numerical evidence is a separate
acceptance condition.

## Counting and scope

Count physical lines, including comments, blank lines and type stubs, in the
complete working `engine/native/` and `engine/src/epcsaft/` trees. Count tracked
and untracked source alike, regardless of diff attributes. Report native layers
and Python separately, plus their total. Reject unknown native layers, runtime
file types and symlinks rather than silently dropping them from the count.
Tests, documentation, source data, build outputs and third-party dependencies
are outside this runtime budget. They remain subject to scope/ownership review;
none may hide runtime logic to evade a ceiling.

Each total entry may be tighter than the sum of its layer entries, never looser:
a layer ceiling bounds that layer, and the total independently bounds the whole.

Physical lines are a reproducible size proxy, not a measure of scientific
complexity. Keep conventional readable formatting. Removing useful comments,
packing statements onto lines, hiding code in generated strings, relocating
runtime code, or splitting/renaming files to evade a budget is unacceptable.
Spare lines never justify a second thermodynamic or solver owner.

## Initial estimate and uncertainty

The September 6 baseline is 8,807 lines in 37 source files. The target allocation
starts with the existing Core, Contributions and admission owners. Properties allows
roughly 675 additional lines for reference thermochemistry, property chains
and compact derivative actions. Equilibrium allows roughly 1,820: historical
planning estimated about 1,400 for stability, with the balance assuming reuse
of the compiler, rows and factorizations for continuation and derivatives.
This historical equilibrium allowance includes search work now deferred by the
owning GitHub issue; it is not permission to build that work or spend the spare
lines. All numeric ceilings stay unchanged. Before each slice, inspect its
current/old owners and justify only the bounded repair, adaptation or missing
capability. Verify-only slices begin with zero new runtime. Reuse decisions and
measured equivalent-work costs, not line count alone, justify replacement.

The parameter regression target follows the retained small fitting-wrapper and
four-node observation design. Binding/Python growth assumes replacement of interface
pieces as all applications migrate; there is no compatibility layer allowance.
The allocation covers the currently documented fluid transfer, including its
required actions, not solids or unrestricted dense derivative tensors. These
estimates are not evidence that the missing science already fits or works.

## Floor, soft ceiling and hard ceiling

- Work toward the floor target. Fewer lines never prove completeness; do not
  add code just to reach it.
- Above a soft ceiling, the command warns. Review each added responsibility,
  reused owner and removable duplication; retain a short reason in the change
  description before acceptance. The soft ceiling grants room to finish and
  simplify, not automatic budget growth.
- Above any layer or total hard ceiling, readiness fails. Refactor/reduce or
  obtain an explicit owner-approved budget amendment with the scientific
  requirement, alternatives and revised numeric table. No override flag or
  automatic ratchet follows the current line count.
- A temporary over-budget working state may be simplified and repaired. Do not
  freeze uncompilable code because an intermediate write exceeded a number.
  The completed candidate must satisfy the hard ceiling before acceptance.

Never satisfy a budget by dropping required physics, domain checks, derivative
terms, failure diagnostics or numerical falsifiers. If the smallest correct
method does not fit, report that evidence and revise the design/budget openly.
Budget changes must not silently change scientific scope.

## Amendments

- 2026-10-02, #238 declared equilibrium starts, measured against integrated
  main `edf7672e`: **374 gross added, 50 removed, +324 net runtime lines**.
  Original estimate 150 gross; first draft measured 227. Lead-authorized cap
  rose to 260 gross. The lead
  authorized up to 275 for admission corrections, then the owner paused hard
  ceilings for this sprint (#248); 275 is now a target. Growth beyond that target
  supplies omitted partial/pinned seed admission, deadlines between initializer
  and admission kernels, and measured rows/work across every failure boundary.
  Previous-coordinate reuse, declaration comparison, the caller-tolerance result,
  and the existing bounded polish remain shared; optimizer iterations and polish
  evaluations are distinct work quantities.

  | Layer | Base | Gross | Removed | Net | Final | Shared hard |
  |---|---:|---:|---:|---:|---:|---:|
  | Core | 917 | 0 | 0 | 0 | 917 | 917 |
  | Contributions | 1547 | 0 | 0 | 0 | 1547 | 1547 |
  | Properties | 4509 | 0 | 0 | 0 | 4509 | 4800 |
  | Equilibrium | 6921 | 299 | 50 | +249 | 7170 | 8800 |
  | Parameter regression | 2398 | 0 | 0 | 0 | 2398 | 2478 |
  | Surface | 532 | 0 | 0 | 0 | 532 | 650 |
  | Bindings | 1886 | 31 | 0 | +31 | 1917 | 2100 |
  | Python | 2635 | 44 | 0 | +44 | 2679 | 2850 |
  | Aggregate | 21345 | 374 | 50 | +324 | 21669 | 24000 |

  One compiled problem, existing generic/boundary seeds, the shared packing
  predicate, exact Ipopt driver and density admission remain the owners.
  The fitting sequence's remaining-iteration and steady-clock accounting pattern
  is reused; the common Ipopt option owner applies positive CPU-time limits
  to exact solves as well as estimation. The driver callback carries the elapsed
  deadline; no second optimizer or acceptance rule is introduced. Deadline
  rejection clears state availability, but the common residual classifier keeps
  its nonfinite/converged/floor/budget priority. Controller causes and measured
  overruns record elapsed exhaustion separately from native termination. The
  native packing coefficient is transported for the analytic density-closed seed.
  The soft excess retains these numerical responsibilities and their work
  diagnostics. The shared #244 envelope and layer ceilings remain unchanged;
  no exclusive allowance is reserved for this task and no ceiling is edited.
  Numerical qualification
  and independent delivered review are separate from this line accounting.

- 2026-10-02, #239 density-diagnostic closeout against integrated main
  `8999d6a2`: **402 gross added, 40 removed, +362 net runtime lines**.
  Core 917, Contributions 1547, Properties 4763, Equilibrium 6921,
  Parameter regression 2398, Surface 532, Bindings 1922 and Python 2707;
  aggregate **21707**. Per-layer gross/removed/net: Properties **294/40/+254**,
  Bindings **36/0/+36**, Python including stubs **72/0/+72**; every other
  layer **0/0/0**. Properties includes the already merged #225 slice 1.
  The sprint pause from #248 applies; no numeric ceiling is changed here.
  The additions report stable and unstable pressure crossings, slope-zero
  brackets, declared-grid crossings separately from adaptive discoveries,
  finite-evaluation boundaries and unresolved brackets in the existing
  density owner. Pressure and slope evaluation, adaptive slope search,
  packing admission and ordinary positive-slope closure remain shared.
  The removed density lines are replaced scan/evaluation blocks, not removed
  physical checks. Delivered-review corrections retain already resolved roots
  by connected-partition bracket ownership, rather than re-refining and
  replacing them; exact sampled zeros survive early exits. The shared slope
  search reports its active bracket on failure. These corrections add 17 net
  lines; the measured gross is two lines above the design's 400-line gross
  ceiling.
  The owner-directed sprint pause applies rather than removing these required
  diagnostics or compressing the existing search. Review checks this growth. Bindings and Python transport explicit inputs and results;
  they perform no root, spinodal or domain iteration. Batch 12 retires its
  Python root/bisection/doubling loops and seven corresponding allowlist rows.
  Numerical tests and retained before/after tables establish qualification
  separately from these counts. The independent
  [delivered review](reviews/issue-239-delivered-review.md) returned Supported
  after the required partial-result and failure-cause corrections, with no
  remaining required findings.

- 2026-10-02, sprint pause (#220, owner decision): hard ceilings are not
  enforced until the current sprint ends. `tools/diffstat.py --budgets` still
  counts every layer and prints FAIL rows, but exits successfully while
  `HARD_CEILINGS_PAUSED` is set. Each build reports its measured gross and net
  runtime lines, and its review justifies growth and checks for duplication.
  After the sprint, a consolidation audit runs and the ceilings are reset
  from measured counts; then the pause flag is removed.

- 2026-10-02, roadmap envelope raise (#220, owner approved): with #223 and
  #225 slice 1 delivered on branches, the measured counts reach about 21,340 total,
  with Properties 4509, Bindings 1886 and Python 2635. The Ready designs still to
  build (#237–#240, #225 slices 2–5, #224) project about 23,700. New hard
  ceilings: aggregate 22000 → 24000, Properties → 4800, Equilibrium → 8800,
  Bindings → 2100, Python → 2850; the other layers are unchanged. Each build's
  closeout lowers the layers it touched to measured plus remaining planned
  allowance; consolidation audits continue.

- 2026-10-02, #223 implementation measurement against `443d8289`, including
  owner-adopted permittivity-form admission: **996 gross added, 19 removed,
  +977 net runtime physical lines**. Equilibrium adds 867, Bindings 43 and
  Python including stubs 67. Current counts are Core 917, Contributions 1547,
  Properties 4274, Equilibrium 6921, Parameter regression 2398, Surface 532,
  Bindings 1886 and Python 2635; aggregate 21110. Every existing hard ceiling
  passes, leaving Equilibrium 654, Bindings 54, Python 65 and aggregate 890
  lines. No ceiling changes or Core/Contributions growth are adopted. The
  soft excess retains one native controller, shared trial rows and exact
  derivatives, conservation and uncertainty calculations, retained failures and
  explicit work diagnostics. The charged-trial admission change is four native
  Equilibrium lines, with four further lines preparing the reacting discovery
  union before existing speciation; two further lines correct the extensive
  addition uncertainty margin and retain its computed value. Zero Binding/Python
  growth. Transport additions carry
  declarations/results rather than search. The
  [qualification notebook](../../analyses/2026-tpd-phase-discovery/notebook.qmd)
  separates admitted numerical checks from the historical solvent-only descent
  and current charged-search refusals. Final wheel qualification and all 4,374
  frozen declared/continuation/fitting/ASL rows pass, with zero implicit trials;
  the initial independent delivered review identified the dimensional margin.
  Its native correction passes the 100 mol amount-weighting falsifier and
  refreshed qualification. The same reviewer's complete affected-lens recheck
  returned Supported with no required defects for the stated cases and limits
  ([delivered review](reviews/issue-223-delivered-review.md)).

- 2026-10-01, #221 finding 5 derivative correction, measured against `bdde4003`:
  **17 gross added, 92 removed, -75 net runtime lines**, all in Equilibrium.
  Equilibrium is **6129 → 6054** and aggregate **20208 → 20133**; every other
  layer is unchanged. The [lever comparison notebook](../../analyses/2026-lever-row-series-comparison/notebook.qmd)
  retains the ionic packing-diameter stencil falsifier, other-column and
  solved-state comparisons, cost measurements and limits. The owner authorized
  the derivative correction after the initial preservation stop: the shared
  row-series geometry replaces the dense lever assembly and its unused
  analytic packing-parameter derivative. Reference refusal, coordinate order,
  shape, free T/P columns, pinned fractions and feed scaling stay with their
  existing definitions. The original 1.53–1.91 times per-call cost is accepted
  by the owner. Separate installed-build measurements under changing shared
  host load gave 1.92–2.35 times the baseline cost (about 0.1–0.4 ms per call);
  they do not isolate algorithm cost as closely as the interleaved study.
  The independent review returned Supported, with no required code correction,
  and recommends owner acknowledgment of that higher observed range. No
  speedup is claimed.

  Closeout lowers the Equilibrium hard ceiling **7650 → 7575**, measured 6054
  plus the **1521-line remaining shared roadmap allowance** from the integrated
  envelope (7650 minus the 6129-line baseline). This preserves that existing
  pool for #223/#225/#224; it is no exclusive allocation or claim that their
  changing estimates fit. Other hard ceilings and the 22000 aggregate remain
  unchanged, with the owner's measured revisit before #224 still required.
  The 54-line soft excess retains the shared equilibrium driver, exact row
  chains, continuation and solved-state actions; no scientific calculation is
  moved outside its runtime owner to meet the count. Numerical qualification
  and independent delivered review are separate conditions of this closeout.

- 2026-10-01, Bindings and Python amendment for #223 (owner approved): after
  #219 slice B the measured counts are Bindings 1843 and Python 2568, total
  20208. #223's estimate (Bindings +35–60, Python +40–75) did not fit. New hard
  ceilings: Bindings 1880 → 1940 and Python 2620 → 2700. The aggregate stays
  22000; the roadmap projects about 22,950, which is decided with measured
  counts before #224. A second consolidation audit aimed at Bindings and Python
  (far above their 650 and 1600 soft ceilings) follows under #221. Pre-envelope
  hard-ceiling changes without amendment entries are listed in the #221
  reconstruction comment (5942662477).

- 2026-10-01, #219 slice B closeout, measured against `96e9ecde`
  (runtime identical to `1fc771b0`, 19,654): **662 gross added, 108 removed,
  +554 net**, total **20,208**. The owner approved the measured sketch through
  [#230](https://github.com/tannerpolley/ePC-SAFT/pull/230), raising Parameter regression
  to 2,480 while keeping every other ceiling and the 22,000 aggregate.
  The final candidate is below the approximately +600 net stop condition.

  | Layer | Base | Gross | Removed | Net | Final | Closeout hard |
  |---|---:|---:|---:|---:|---:|---:|
  | Core | 917 | 0 | 0 | 0 | 917 | 917 |
  | Contributions | 1547 | 0 | 0 | 0 | 1547 | 1547 |
  | Properties | 4274 | 0 | 0 | 0 | 4274 | 4600 |
  | Equilibrium | 6095 | 86 | 52 | +34 | 6129 | 7650 |
  | Parameter regression | 2025 | 418 | 45 | +373 | 2398 | 2478 |
  | Surface | 532 | 0 | 0 | 0 | 532 | 650 |
  | Bindings | 1766 | 80 | 3 | +77 | 1843 | 1880 |
  | Python | 2498 | 78 | 8 | +70 | 2568 | 2620 |
  | Aggregate | 19654 | 662 | 108 | +554 | 20208 | 22000 |

  Parameter regression is three net lines above the approximately 370 estimate.
  One TNLP assembles target/max pairs and mean-slack links over the existing
  selected-row evaluator and eliminated-hold chart. Shared check arithmetic,
  admission and candidate/check publication replace their inline Ceres
  producers. The existing native sequence retains starts and shared budgets.
  Equilibrium's +34 is one line above the estimate: its common checked Ipopt
  lifecycle/options and exact return names replace the old option/name blocks;
  equilibrium retains exact curvature and file initialization. Bindings +77
  exceeds the +55 estimate because labelled constraints, frozen checks,
  physical multipliers, root-error reporting, and the same native one-point
  observation/constraint evaluation are all transported. The formula-composition
  binding reuses the existing molality owner. Python +70 stays below +80 and
  contains declarations, calls, exports and type stubs only. No thermodynamic
  equation, root owner, optimizer fallback or profiling framework is added.

  Closeout lowers Parameter regression to **2,398 measured + #225's 80 planned = 2,478**.
  Per the owner instruction, other ceilings and the aggregate remain unchanged.
  Remaining shared hard headroom is Equilibrium 1,521, Bindings 37, Python 52,
  and aggregate 1,792. These are ceilings, not claims that #223, #225 and #224
  fit; those leaves remain and must remeasure their own designs. Their scope or
  allowance is not silently changed by this closeout. Soft excess pays for the
  exact held-chart estimators, row-specific constraints and numerical reporting.
  [B1-B4 numerical verification](../../analyses/2026-electrolyte-parameter-gradients/notebook.qmd)
  identifies the wheel, operands, tolerances and local-solution limits.

- 2026-10-01, #219 slice B Parameter regression amendment (owner approved): slice B
  stopped at its budget gate. A counted sketch of the D7 native Ipopt
  constrained estimation (one Parameter regression TNLP over the slice-A evaluator) gives
  about +370 Parameter regression lines and +538 net total, against the 175-line
  allowance left after slice A. The owner chose to build it now: Parameter regression hard
  ceiling 2200 → 2480, covering slice B plus #225's planned +80. Other ceilings
  and the 22000 aggregate are unchanged. The roadmap now projects about 22,260
  against 22,000, so the aggregate is revisited with measured counts before #224
  builds. Slice B's closeout records the measured counts and lowers Parameter regression
  to measured plus the remaining #225 allowance.

- 2026-10-01, #219 slice A closeout, measured against `89c1a337` (18,932):
  **926 gross added, 204 removed, +722 net**, total **19,654**. The owner
  approved 950 gross / 770 net after the 850/770 intermediate allowance;
  the design's 610/550 figures were estimates. Counts include comments,
  blanks and type declarations; tests, documentation and retained numerical
  tables are separate.

  | Layer | Base | Gross | Removed | Net | Final | Envelope hard |
  |---|---:|---:|---:|---:|---:|---:|
  | Core | 917 | 5 | 5 | 0 | 917 | 1000 |
  | Contributions | 1593 | 22 | 68 | -46 | 1547 | 1612 |
  | Properties | 4176 | 106 | 8 | +98 | 4274 | 4600 |
  | Equilibrium | 6095 | 0 | 0 | 0 | 6095 | 7650 |
  | Parameter regression | 1607 | 533 | 115 | +418 | 2025 | 2200 |
  | Surface | 532 | 0 | 0 | 0 | 532 | 650 |
  | Bindings | 1674 | 96 | 4 | +92 | 1766 | 1880 |
  | Python | 2338 | 164 | 4 | +160 | 2498 | 2620 |
  | Aggregate | 18932 | 926 | 204 | +722 | 19654 | 22000 |

  The largest underestimate was parameter regression: 533 gross against
  the design's upper 285. Native check/baseline admission, signed metrics,
  unavailable members and stored-parameter rebuilding cost more than the
  roughly 25-line check estimate. The shared selected-row observation
  evaluator, native manufacturing and multistart/shared-budget operations,
  positive-objective SVD labels and dependent-chart failures were also larger
  than estimated. Bindings grew 96 against 55, and Python 164 against 100:
  these transport the same hold/check/chart/sequence records and preserve
  distinct result types. Properties, Core and Contributions stayed below
  their upper estimates. Removed duplicated fit preparation/reporting and
  unused contribution overrides account for the measured 204 deletions;
  those deletions do not buy extra gross allowance. Every added operation
  serves the five observations, exact held chart, final checks or slice-B
  interfaces. No new EOS physics or optimizer is included.

  **Parameter regression uses 2,025 of 2,200, leaving 175 for slice B**,
  against D7's estimate of up to 195. The owner [explicitly approved this
  measurement](https://github.com/tannerpolley/ePC-SAFT/issues/219#issuecomment-5941530617)
  rather than broader fitter changes. Delivered-review corrections preserve
  this allowance: checks rebuild only their selected rows and zero-weight
  reports, while objective rows keep the accepted evaluation; unevaluated rows
  remain unavailable. Removing unnecessary internal aliases and access
  specifiers keeps parameter regression at the approved measurement. That handoff required
  B to re-estimate before Build; its subsequent amendment below supersedes the
  175-line limit.
  The native profile driver remains deferred under D8: fixed-coordinate
  removal, declared starts, shared limits, neighboring starts and per-point
  failure/check records exceed the 35-line driver allowance in the current
  owner. SVD and covariance are delivered instead.

  Closeout lowers Core to 917 and Contributions to 1547 (no remaining
  additions planned there). Other touched hard ceilings retain shared
  remaining allowances: Properties 326, parameter regression 175, Bindings
  114, Python 122 and aggregate 2346. These are limits shared by the remaining
  roadmap work, not exclusive grants or a claim that its older estimates fit.
  Remaining designs must remeasure against the integrated tree. Soft excess
  pays for the shared EOS/exact actions, native held-coordinate fitting and
  named physical declarations; no line packing or displaced calculation was
  used to meet a ceiling. [Numerical qualification](../../analyses/2026-electrolyte-parameter-gradients/README.md)
  retains F1-F6, the acceptance-check falsifier and diagnosed misses. This
  establishes numerical verification at the declared inputs, with no #209
  parameter adoption or physical validation.

- 2026-10-01, general-engine roadmap envelope (#220, owner approved): the
  owner replaced per-issue amendments for #219, #223, #224 and #225 with one
  envelope. Measured base `4cd1ebfe` (after the #221 consolidation, −168): total
  18932; Core 917, Contributions 1593, Properties 4176, Equilibrium 6095,
  Parameter regression 1607, Surface 532, Bindings 1674, Python 2338. Projection from the
  designs' upper estimates (unmeasured; #225's broadened scope pending): about
  21,600. New hard ceilings: Properties 4258 → 4600, Equilibrium 7000 → 7650,
  Parameter regression 1619 → 2200, Bindings 1726 → 1880, Python 2374 → 2620 and aggregate
  19100 → 22000; Core 1000, Contributions 1612 and Surface 650 unchanged. The
  aggregate is deliberately tighter than the sum of the layers. The #219 D6/D7
  and #225 D4 per-issue requests are superseded. Each build's closeout records
  its measured counts and lowers every layer it touched to the measured count
  plus the allowance still planned for the remaining roadmap issues; this is an
  explicit recorded amendment at each closeout, not an automatic ratchet. A
  design or build that cannot fit the envelope returns to the owner. Floors and
  soft ceilings are unchanged; soft-ceiling reasons are still required.

- 2026-09-30, direct source-basis ln K observations (#211, owner approved):
  Parameter regression hard ceiling 1463 → 1484, Python 2355 → 2357 and aggregate
  18805 → 18830, exactly the measured counts against `95ba231e` (+25 net
  runtime lines: Parameter regression +21, Bindings +2, Python type declarations +2).
  Bindings measures 1723 within its unchanged 1725 hard ceiling. The direct
  value uses the existing reaction-property owner; exact columns reuse the
  existing reaction-coordinate basis and admission reuses named-coordinate
  rejection. No equilibrium solve or source-reference transfer is required.
  Floors, soft ceilings and other hard ceilings remain unchanged.

- 2026-09-30, total molar isobaric heat-capacity fitting (#205 O2, owner approved):
  Parameter regression hard ceiling 1440 → 1463, Python 2350 → 2355 and aggregate
  18790 → 18805, the measured candidate counts. Relative to `d47d0475`,
  runtime growth is 37 physical lines: Parameter regression 23, Bindings 2 and Python 12.
  The owner approved 65 gross added runtime lines, with 28 removed (+37 net),
  replacing the design's gross-45 ceiling while retaining all net and layer caps.
  Properties remains 4223; Bindings remains under its 1725 ceiling.
  The pressure-closed observation shares Density admission and fitting paths,
  uses the existing Properties caloric action once per parameter column, and
  keeps fixed ideal-gas thermochemistry. Python reuses its existing species
  validator. No new thermodynamic owner or parameter adoption is included.
  Floors, soft ceilings and other hard ceilings are unchanged.

- 2026-09-30, parameter regression hard ceiling 1300 → 1440, bindings 1708 →
  1725, Python 2313 → 2350 and aggregate 18556 → 18790 (#187 O1, explicitly
  approved in the owner's September 30 comment). The candidate measured against
  `fe00aa6a` after integrating #200 uses parameter regression 1440, bindings
  1719, Python 2343, equilibrium 6018 and aggregate 18765 physical runtime
  lines. Net additions against that main are 211: parameter regression 142,
  bindings 11, Python 32 and equilibrium 26.
  The new responsibilities are named source-basis reaction coordinates,
  reaction shifts in the shared solved-state first actions and component-amount
  ratio observations. The minimum-norm mapping, defining-Jacobian factorization,
  affine physical-value map and reaction-declaration equality retain one owner
  each. Validation preparation is reused only when its actual EOS values and
  resolved reaction declarations match exactly; a reaction-coordinate change
  invalidates the compiled preparation. A parallel derivative implementation and application calibrations are
  excluded. The equilibrium soft-ceiling excess is 18 lines for the reaction
  controls, admission and shared row term. Floors, soft ceilings and all other
  hard ceilings stay fixed. Tests and numerical evidence remain separate from
  these counts; no physical validation is claimed by the synthetic cases.


- 2026-09-29, Properties hard ceiling 4211 → 4223, parameter regression
  1200 → 1300 and aggregate 18365 → 18556 (#129/#193, explicitly owner-approved for the
  measured candidate). Requirement: exact request-local equilibrium-state
  sharing, residual-to-Jacobian reuse at the same parameter point, and held-set
  refresh from the accepted full weighted gradient. Preserve every observable,
  action, residual transform, weight, failure safeguard and numerical tolerance.
  Measured growth against `11ee6fbe` is +194 runtime lines: Properties +12 for
  reference-work tracing and parameter regression +182. Consolidation removed ten
  lines from the first candidate by sharing gradient and membership responsibilities.
  Alternatives: retaining repeated scalar work loses the measured initial-call
  saving; removing release, gradient or failed-trial safeguards invalidates the
  required behavior. No further safe duplicate removal was demonstrated within
  this slice; this is not proof of a minimum. Wider refactoring needs its own
  evidence. Floors, soft ceilings, other layer ceilings and scientific scope
  remain unchanged. The [fitting notebook](../../analyses/2026-fitting-state-reuse/notebook.qmd)
  retains the numerical comparisons and limits; complete installed-wheel
  verification remains a separate acceptance condition.

- 2026-09-26, aggregate hard ceiling 18353 → 18365 (#172, owner-approved
  2026-09-26 in the #172 decision comment for the measured growth, estimated at
  15–25 lines). Requirement: all-free equilibrium states whose trace-species
  inventory lies below the composition box times a phase amount (Lithium's
  five-stage lean end) were infeasible or stalled, and their actions were refused.
  Measured net growth against e5687520 is +19: Equilibrium +10 (the per-species
  box from one HiGHS linear program per species on the already linked solver,
  +26; the moiety echelon balance basis, +21; the unused runtime
  `same_subspace`, which moves into its native test, −37) and Core +9 (the
  `EQUILIBRIUM_BOX_INVENTORY_FRACTION` policy entry). Alternatives: a
  structural-floor box for every cell (removes the barrier's regularisation and
  needs full requalification), a row-ratio bound (invalid on the mixed-sign
  balance basis) and a box scaled for every minor species (weakened the
  continuation corrector). No layer ceiling or tolerance changes.

- 2026-09-24, Bindings hard ceiling 1700 → 1708, Python 2300 → 2313 and aggregate
  18145 → 18353 (#152, owner-approved 2026-09-24 in two steps: +187 as first
  measured, then +21 from review-required correctness fixes). Requirement: MEA's
  reactive refit through `epcsaft.regression.fit` with exact derivatives, no
  slower than its scipy forward-difference refit and reaching its objective.
  Measured net growth against 34d06cd2: the bounded fit (active set held in a box
  chart, a Jacobian cache at accepted points, per-iteration reporting;
  parameter regression +123 with the coordinate bounds owner), without which
  refit C stalled at 123.38 against scipy's 121.36; the speciated start of a cold
  reactive boundary with its fallbacks (Equilibrium +67 with the reference-jet
  reuse and the compiled solve), which removed four cold failures and the
  3–21 s cold tail; the partial-pressure observation and batched first actions;
  the public surface (Bindings +14: partial_pressure, the slope family,
  `differentiate`, iteration reporting; Python +19 with the stubs); the
  `FIT_ACTIVE_BOUND_DISTANCE` policy (Core +7); offset by one family→storage map
  replacing duplicated switches (Properties −22). Alternatives rejected:
  Ceres' own bounds (projected line search: 45 Jacobian evaluations on six steps;
  per-coordinate clamp: stalled), snapping onto bounds (four review rounds of
  defects), dropping the reporting item (owner kept it). Other layers stay under
  their hard ceilings; no thermodynamic method or tolerance changes.

- 2026-09-23, Contributions hard ceiling 1500 → 1612, Properties 4145 → 4150 and
  aggregate hard ceiling 17760 → 17903 (#28). The owner pre-approved at most +250 runtime lines for an
  optional Gross–Vrabec dipole–dipole term (ADR 0011) and reserved quadrupole and
  cross terms for a new approval. The measured cost is +143 lines: Contributions
  +112 (Table 1 constants, the J-integral interpolation and the Padé term), Surface +19
  (explicit `polar` family and unit-checked `dipole_moment`/`dipole_count`
  admission), Properties +5, Core +4, Bindings +4 (Python and ASL document
  transport) and Python −1. Each raised ceiling rises by exactly the measured polar count (+112,
  +5 and +143), so headroom that existed before this term, including the #84 use, is
  not consumed.

- 2026-09-22, aggregate hard ceiling 17600 → 17601 (#66). The latest merged
  main adds one Python physical line to the 17600-line accepted integration.
  The owner authorized necessary IDAES integration budget increases. No layer
  ceiling, thermodynamic method, numerical tolerance or acceptance claim changes.

- 2026-09-22, aggregate hard ceiling 17450 → 17600 (#66). The owner approved
  necessary budget increases and explicitly requested integration of the latest
  local reference/trace-recovery and Lithium surrogate-proof work. Combining
  those reviewed Engine changes with the existing live adapter uses about
  17544 physical lines. Keep the adapter Bindings ceiling 1700 and the upstream
  layer ceilings; no new thermodynamic owner or compatibility path is added.
  Removing required reference derivatives or adapter admission checks, or packing
  code to conceal the combined size, would invalidate the requested integration.
  This bounded amendment covers composition of the two lanes, not new physics.

- 2026-09-21, Properties hard ceiling 3460 → 3470 (#89, owner-approved small
  increase). Requirement: replace the density owner's fixed-width slope search
  with a bounded search using a quadratic-model estimate and an
  explicit unresolved failure and insert only loop-exposing samples, after a
  fixed-width stop lost a liquid root 1e-4 K below the critical point.
  Measured use 3468 lines with single-statement guards in `solve_density`
  (3478 with the previous three-line guard blocks; both forms are conventional
  and the file's neighbours use both), including the explicit bracket-collapse
  failure, its `<limits>` include and the endpoint rule of the refinement
  added in the correction review rounds. Alternatives costed and rejected: merging
  the anchored Newton loop into `refine` would obscure its typed left-domain,
  singular and infeasible stops; moving the search into the equilibrium layer
  would split the density owner's responsibility; dropping the golden-section
  constant's comment or the model-limit comment would remove the stated
  domain. The quadratic interpolation is not a certificate for unsampled EOS
  curvature. This amendment left two lines above its measured 3468-line use;
  any later change exceeding the ceiling needs its own justified amendment.

## Duplication constraints

One common Helmholtz definition supplies properties and their derivatives.
Association elimination, density closure, reference transformations, equilibrium
row assembly and physical-input resolution each retain their accountable owner.
Python maps inputs/outputs; applications supply chemistry and specifications as
data. All equilibrium specifications use one driver; association, density,
fitting and discovery retain their distinct numerical operations using the
shared owners. Existing architecture checks cover known
include/constant/call patterns; size and regex checks cannot establish absence
of semantic duplication. Scientific/code review must inspect the actual calls.

This budget is a local readiness rule, not CI/CD and not a numerical accuracy or
runtime-speed threshold.
