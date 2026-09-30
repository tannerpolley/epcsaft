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

- 2026-09-30, total molar isobaric heat-capacity fitting (#205 O2, owner approved):
  Regression hard ceiling 1440 → 1463, Python 2350 → 2355 and aggregate
  18790 → 18805, the measured candidate counts. Relative to `d47d0475`,
  runtime growth is 37 physical lines: Regression 23, Bindings 2 and Python 12.
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
