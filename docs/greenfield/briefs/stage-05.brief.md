# Stage 05 — Equilibrium

Branch `stage-05`, sources pinned at `greenfield@88b5a2a`; attempt 4 starts from the tag `stage-05-attempt-3` (`d02f3c9`). Machine-readable claims are in
`briefs/stage-05.brief.json`; that file governs where the two disagree.

**`EQUILIBRIUM.md` is the specification, not background.** Every row, every rejection, every
check and the mutation that must redden it, and the twenty-five traps of §10 are decided
there. This brief is an index and a set of budgets. Where the two differ, `EQUILIBRIUM.md`
wins — and that is a `source_contradicted` stop, not a judgement call.

## What this stage is

One equilibrium, one driver, one stage. It was four — a specification stage, a driver stage,
a starts stage and a boundary stage — and every split was route thinking that had already
been removed from the code. There are no rows for a second stage to write: compile emits them
once and the driver solves them for every `k`.

    rows − unknowns = k − f

with `k` pinned phases and `f` free intensives; square iff `f = k`, and rank at a generic
point is the sufficient check. The flash is `k = f = 0`, a boundary is `k = f = 1`, a
three-phase point is `k = f = 2`. **No route is named anywhere.** Bubble and dew differ only
in which phase is pinned.

## The order it is written in

`EQUILIBRIUM.md` §9 sets this, and it is the cheapest preflight the stage has:

1. `equilibrium_support.hpp`, then **C1–C3**, before any solver exists. They are green on the
   compile step alone. If C1's predictor and the compile step disagree on a cell, one of them
   misreads §2.3 — and that is found before a line of `gibbs_minimizer.cpp` is written.
2. The rest of the compile tier, C4–C9.
3. The driver, the one-component cell, starts and topology; S1–S16 except S12, which is
   Stage 7's.

If you cannot finish in one pass, stop cleanly at a numbered boundary above and report which.
A partial stage reported honestly is worth more than a complete one that guessed.

## Traps

`EQUILIBRIUM.md` §10 is twenty-five items and **the brief copies that list by reference**:
read it before writing, and again before reporting. They are the places an implementer would
otherwise guess, or where a check as first written could only be passed by fudging. Four of
them are corrections the design made to its own first draft. The ones most likely to be hit
silently: 1 (unit-amount tapes), 5 (the assembly lives in `request_compile.cpp` and the driver
never writes its own), 7 (compare against `max(tolerance, floor)` from the record, never a
bare policy number), 10 (reach the tape through the existing predicate-signature cache key),
15 (one electroneutrality row per ion-containing phase *except one of nonzero amount*), and
21 (no identifier under `native/` is named for saturation).

## What you inherit and may not change

- `CORE.md` §5.2: a solver never names a residual stop. Call `classify`. Rule 14 enforces it.
- `CORE.md` §6.1: every bound function opens exactly one `Scope`. Rule 13 enforces it.
- `nlp_scaling_method = none`, bounds on `N_p ≥ 0`, and **no outer Newton around
  `solve_equilibrium`** — three commitments S1, S15 and this stage's gate inherit.
- Python does no arithmetic. Feed-basis conversion is arithmetic, so it is native.
- Include direction: `equilibrium/` may include `model/`; nothing in `model/` may include
  `equilibrium/`.

## Budgets

The per-file `low`, estimate and ceiling in the JSON are what you are held to. They sum to
2,450 / 3,240 / 4,130 for the eight native files, 2,900 / 3,890 / 5,010 with the bindings and
Python, and 3,423 / 4,649 / 6,065 across all seventeen writes. `TRANSFER.md`'s heading now
carries the same numbers. Do not trim a file to meet a total. If a file genuinely needs more
than its ceiling, that is `budget_exceeded` — stop and report the number.

## The checks are the stage, and they are named

A first attempt passed all twelve gates and was rejected. The gates ran the meson targets;
the targets were smoke tests. `boundary_checks.cpp` was twenty-seven lines asserting that
sizes are finite and nullity is zero. **C1–C9 and S1–S16 did not exist.**

So this is now explicit. `EQUILIBRIUM.md` §4.1 names nine compile-tier claims and §4.2 names
sixteen solve-tier claims (S12 is Stage 7's, and S7 has no mutation by design). **Each one is
a function named for it** — `check_c1`, `check_c4`, `check_s10` — called from `main` in the
file §9 assigns it, asserting the claim as §4 words it, not a shape or a finiteness proxy.
You also write `tools/check-equilibrium-claims.py`, which reads the names out of
`EQUILIBRIUM.md` and fails if any is missing or uncalled. A claim with no function is a
missing check; a function that is never called is worse.

Concretely, the ones the first attempt proxied and must not be proxied again:

- **C1** walks the lattice — every cell from the four enums, thousands of specifications,
  none sampled — comparing compile's verdict to a predictor built from the specification
  alone. One loop over three feed bases is C4, not C1.
- **S10** compares `F_y` against a **centred difference** of `F` in every unknown. Asserting
  `jacobian.size() == values.size() * point.size()` is a shape, and a shape is not a
  derivative.
- **S3** solves bubble, then solves dew *seeded cold* from the result and recovers the input.
- **S9** compares `dy/dθ` against a centred difference of the whole solve.
- **S1** and **S13** read `max(tolerance, floor)` from the result's `Convergence` record.

## Attempt three: what the second attempt got wrong

The second attempt passed fourteen gates and was rejected on reading. Its assembly had no
analytic `F_y`, no `F_θ` and no exact Hessian: `CompiledSystem::jacobian` was a centred
finite difference of `evaluate`, `eval_h` called a `numerical_hessian`, and no line read
`derivative.hessian`. §2.4 and §5 item 3 say all three come from the tape. Every consequence
below follows from that one omission.

**Derivatives are the assembly's, never the test's.**

- `F_y` is `CompiledSystem::jacobian`, analytic from the tape's Hessian and the
  log-coordinate chain rule. `request_compile.cpp` must contain `derivative.hessian[`.
- `F_θ` is the assembly's `lever_jacobian(point)`, `θ ∈ {Temperature, Pressure}`: the tape's
  Hessian column in `T` for the chemical-potential and pressure rows, `−1/RT` on the pressure
  row for `P`, zero elsewhere. `implicit_sensitivity.cpp` owns `sensitivity(compiled, result,
  lever)` and must contain `-F_theta` and `lever_jacobian`. The model's active parameters are
  the same call with more columns and land with regression.
- Ipopt's `eval_h` reads `compiled_.hessian`. `gibbs_minimizer.cpp` must contain it.

**S10** centre-differences `evaluate` itself in every unknown at the boundary cell
(`boundary_checks.cpp`) and at the flash cell (`equilibrium_checks.cpp`) and compares to 1e-6;
and compares `eval_h` at the flash solution with the converged multipliers against a centred
difference of `eval_grad_f + Σ λ_j ∇c_j` in every primal unknown, to 1e-6. No
`hessian_calls` assertion. No `isfinite` sweep.

**S9** compares `dy/dθ` of the binary flash at fixed `T, P` against a centred difference of
two full solves, for each of three levers: `T` (`δ = 1e-3 K`), `P`, and one species-feed
coordinate `n_feed,i` that carries a moiety. `F_θ` for the feed column is `−B[:, i] / Σ_j n_feed,j`
on the balance rows and zero elsewhere, one extra column per species in
`implicit_sensitivity.cpp`; 1e-6 relative for every lever; and on the reactive fixture the
derivative along a reaction direction (`B d = 0`), formed from the species columns, is exactly
zero. The conserved amounts `b = B n_feed` are not a lever. The test supplies no derivative
of anything.

**Eight thin checks must assert what §4 words.** C6, S4 and S5 are adequate as one
assertion. These are not:

| check | attempt two asserted | write instead |
|---|---|---|
| C4 | composition to 1e-12; `rows.size()` | pairwise `rows[i].kind, phase, index, eqid` equal across the four bases |
| C7 | nullity and `rows.size()` equal | the permuted point built from the first compile's generic point through the layout (a second `generic_point()` is not the permutation: `η_p` depends on phase index); `residual_norm()` and sorted `|values|` equal to 1e-14 under phase reorder; then relabel components by permuting the fixture's per-component vectors, `k_ij` and the feed (helper in `equilibrium_support.hpp`) and compare the same |
| S1 | one `{L}` cell, residual rows only | every solved cell (`{L}`, `{V}`, `{V,L}` binary; ethane boundary; the reactive `{L,L}`; electrolyte `{L}`); per cell `limit = max(tolerance, floor)` from `result.convergence`; every residual, `g_p − λ·B x_p` per phase, `G − λ·B n_feed`, `Σ_i ν_i μ_{p,i}` per reaction per phase, `|z·x_p|` per ion phase, all `≤ limit` |
| S2 | `reduced_hessian_minimum > 0` on `{L}` | the same field on the `{V,L}` binary flash and the reactive `{L,L}` flash, computed from the assembly's exact Hessian with the converged multipliers, never `numerical_hessian` |
| S6 | `objective` at two generic points, slack 1e6 | solve `{L}` and `{L,L}` on the ternary without the reaction at S7's feed; `G₂ < G₁ − limit`. The `G₃` clause runs when the `{V,L,L}` cell is admitted, else is reported skipped by name |
| S7 | two scaled feeds, no solve | solve `{L,L}` once; feeds at 0.3 and 0.7 along the segment between the two phase compositions; re-solve; both compositions equal to 1e-8 |
| S11 | the affinity transform on a hand-built reaction | solve the reactive `{V}` cell with the isomerisation at `P = 10⁴, 10³, 10² Pa`; `gap(P) = |ln(x_β/x_α) + Δg°/RT|` decreasing and `gap(10²) ≤ 1e-3` |
| S16 | rank and `‖n_feed‖` on the neutral binary | solve `ascani_reactive_problem()` (it exists in the support file; attempt two never called it); three feeds on the segment between the phases' moiety vectors `B x_p`; re-solve each; the same two compositions to `limit`; `Σ_i ν_i μ_{p,i} ≤ limit` in both phases with one `λ` |

**Mutations are applied to your tree, and a mutant that does not build proves nothing.**
Seven gates in the JSON carry a sed against a line you write, each tied to a `must_contain`
token so the sed is stable. A mutant that fails to compile, or on which the gate stays green,
is `mutation_invalid` and rejects the attempt. The `boundary_checks` mutation in particular
edits `evaluate` alone: an analytic `F_y` goes red, a finite-difference `F_y` follows the
edit and stays green, and green is the rejection.

## Line counts are a signal, not a formality

Every file in the first attempt came in below its `low` — 1,826 lines against 4,085. `low`
is the design's estimate of what the file must *contain*. Coming in under it is now a
`budget_undershot` stop: report the file and the number rather than declaring the stage
complete. If you believe a file genuinely needs less than its `low`, that is a finding about
`EQUILIBRIUM.md` worth reporting — say so and stop, and the coordinator will take it to the
planner.

## Stopping

Stop and report rather than adapting if: a `must_contain` marker is absent; a write would pass
its ceiling; a gate green at base goes red for a reason other than the write it names; the
work seems to require a route name, a second driver, an outer Newton, or a change to a
tolerance or to `classify`; a write comes in under its `low`; a gate's mutation cannot be applied to what you wrote
(`mutation_invalid`); or **`EQUILIBRIUM.md` does not decide something the code must
decide** — report the cell or the row and stop. Do not choose. Do not repair your environment,
and do not accept your own work: the coordinator rebuilds from a clean directory, re-runs
every gate, and applies every mutation independently.

## Attempt four: what the third attempt got wrong

The third attempt stopped honestly with two solve-tier checks red, and the planner read
them as one defect: the coincidence refusal was guarded by `pinned_count == 0`, so every
boundary cell returned the trivial solution (`x_L = y_V`, equal densities) as converged.
Behind it, four places the design let you guess. Each is now written, and each has a gate.

- **Coordinates and bounds** (`EQUILIBRIUM.md` §2.1): `ln x`, `ln η` with the packing
  ceiling `PACKING_FRACTION_CEILING = 0.74`; a plain volume with no ceiling admitted the
  `v → ∞`, `P → 0` bubble point.
- **The rows are dimensionless** (§2.2): the pressure row is `Z − Pv/RT`, the balance rows
  are divided by the feed total.
- **The coincidence rule is a mechanism** (§3.4): `d = max(max_i |Δx_i|, |Δ ln η|)`, refused
  at `PHASE_COINCIDENCE_DISTANCE = 1e-6`, on **every** multiphase cell and every pair,
  pinned or free. S15 gains clause (c) and the gate `coincidence_every_cell`.
- **The start rule** (§5 step 4, trap 6): the caller's guess; else a one-step K-factor split
  from the two density branches for a `V`/`L` pair; else the refusal `start_symmetric` for an
  `L`/`L` pair. No other. Gate `start_splits_composition`.
- **The rows' own tolerance** (§5, tolerance and budget): `EQUILIBRIUM_ROW_RESIDUAL = 1e-10`,
  `EQUILIBRIUM_MAX_ITERATIONS = 200`, both in `numerical_policy.hpp`, which joins `writes`.
  Nothing under `equilibrium/` names `DENSITY_ROOT_RESIDUAL`. Gate `own_row_tolerance`.
- **Fixtures carry states** (§4.2, the state table): vapour–liquid cells run on the new
  `neutral_ternary()` (methane, ethane, propane at 180 K); methane–ethane keeps the `C = 2`
  compile cells and S3; MEA–water is the Cai 1996 isobaric oracle and never a synthetic
  state; the Ascani ternary without its reaction is the neutral two-liquid example. S6
  certifies every synthetic state before its dependents run; a failed certificate skips them
  by name.
- **Checks report, never abort** (§4 preamble): every solve under `core::Capture`; a failing
  check prints its trace with every Ipopt iteration. Gates `checks_do_not_abort` and
  `dew_seeded_cold` (trap 13: the dew is `solve_equilibrium(model, dew)` with no start).
- **Diagnostics** (§4.4): the brief's `diagnostics` array, which the preflight now requires.

**A budget breach is a stop, never a formatting task.** `request_compile.cpp` is budgeted
600 / 750 / 900 (low / estimate / ceiling): attempt 3 measured 716 lines of real content
with the analytic assembly in it, and this revision adds the feed column, the dimensionless
rows and the `standard_gibbs` call. If a file you write crosses its ceiling, stop
`budget_ceiling` with the number. Combining statements or removing blank lines to fit a
count is rejected on reading whatever the total, and the file is measured as first written.

Three test budgets moved after the first honest ceiling stop (333 lines in
`equilibrium_checks.cpp` with S1's four consequences per cell still unwritten):
`equilibrium_checks.cpp` 230 / 380 / 450, `specification_checks.cpp` 180 / 240 / 300,
`equilibrium_support.hpp` 150 / 250 / 320. The ceilings are 1.2× measured-plus-remaining;
the lows stay where content is guarded.

You start from the tag `stage-05-attempt-3`. Carried: the exact-tape assembly, analytic `F_y`
and `F_θ`, `stoichiometry`, `problem.hpp`, the bindings, the Python and the compile-tier
checks. Rewritten under the sections above: the driver's starts, the coincidence refusal, the
coordinates and bounds, the tolerance, the check macro and the fixture states. The
attempt-three rulings on S9, S10 and the eight thin checks stand.
