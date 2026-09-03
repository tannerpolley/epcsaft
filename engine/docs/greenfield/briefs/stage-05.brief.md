# Stage 05 — Equilibrium

Branch `stage-05`, based on `greenfield@88b5a2a`. Machine-readable claims are in
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

## A budget inconsistency, recorded rather than resolved

`TRANSFER.md` heads this stage "~3,040 lines, ceiling ~3,870", then lists files summing to
about **3,710 before tests, bindings or Python**. The headline is an older decomposition
(specification 1,000 / driver and cell 1,330 / starts and topology 710) and the file list is
newer. **The per-file ceilings in the JSON are what you are held to**; the headline is
informational and is the coordinator's to reconcile, not yours. Total estimate across all
sixteen writes is ~4,580 with a ceiling of ~5,955, and that is expected to be larger than
either headline because it counts the four test files, the bindings and the Python layer,
which the headline never did.

Do not trim a file to meet a headline. If a file genuinely needs more than its ceiling, that
is `budget_exceeded` — stop and report the number.

## Stopping

Stop and report rather than adapting if: a `must_contain` marker is absent; a write would pass
its ceiling; a gate green at base goes red for a reason other than the write it names; the
work seems to require a route name, a second driver, an outer Newton, or a change to a
tolerance or to `classify`; or **`EQUILIBRIUM.md` does not decide something the code must
decide** — report the cell or the row and stop. Do not choose. Do not repair your environment,
and do not accept your own work: the coordinator rebuilds from a clean directory, re-runs
every gate, and applies every mutation independently.
