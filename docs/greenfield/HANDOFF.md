# Coordinator handoff

Written 2026-09-03 for a fresh coordinator session. Read this, then `README.md` in this
directory, then `ORCHESTRATION.md`. Everything else is reached from those.

## State at handoff

- Branch `greenfield`, worktree `~/Workspaces/Engineering/ePC-SAFT-greenfield`,
  **nothing pushed** (push only when the owner asks).
- Stages 1–4 closed. Stage 5 (Equilibrium): attempts 1 and 2 held as tags
  `stage-05-attempt-1`, `stage-05-attempt-2`; **attempt 3 stopped `gate_impossible`** on
  2026-09-03 with the local-minimum and dew–bubble checks red (`STAGE-LOG.md`, last entry).
  Not tagged, merged or resumed. The Equilibrium planner has the two failures and is
  discussing design versus implementation with the owner; a resume waits on that ruling.
- The control-plane reset the owner approved on 2026-09-03 is applied (`STAGE-LOG.md`,
  last entry). Briefs and preflight untouched.
- Records moved to the root today: `docs/scientific/` (CSE core, issue-55 shape) and
  `docs/greenfield/`. Engine-only docs stay under `engine/docs/`.
- Rule 13 has no allowlist (owner ruling). The wheel fingerprint now hashes `engine/native`;
  before today a native-only edit reinstalled a stale cached wheel.
- Overview page: https://claude.ai/code/artifact/1bb42e8c-f3ef-4c72-99e1-335f7614bf58 —
  keep the head sha in its header synced after every commit (owner rule).

## Attempt 3, the stopped process

| item | value |
|---|---|
| worktree | `~/.codex/worktrees/gf-stage-05`, branch `stage-05`, cut from `dd81eb5` |
| worker | `codex exec`, model `gpt-5.6-luna`, reasoning high, session id `01a0684e-7289-7ca3-bb24-3f10b04a82e8` |
| pid | 361705, exited |
| log | `attempt-3.log` in the worktree; result lands in `stage-05.result.json` |
| schema | `stage-05.result.schema.json` in the worktree (untracked) |
| brief | `docs/greenfield/briefs/stage-05.brief.{md,json}` at `dd81eb5`, 19 gates, preflight clean |

At last look the worker had analytic `F_y`, exact-Hessian differences and the C4/C7
rewrites in place, and was tracing a residual `5.05e-12` against a `1e-12` record
tolerance in the single-phase cell plus a boundary relation not yet passing.

**Acceptance is the coordinator's, never the worker's.** When the result file appears:

1. Read `outcome`. Anything but `complete` is reported to the owner, not resumed alone.
2. Rebuild from clean in the worktree; run every gate in the JSON.
3. Apply every mutation with `sed -i`, rebuild, confirm the gate goes **red**, restore.
   A mutant that does not build, or stays green, is `mutation_invalid`: reject.
4. Read the checks. Attempt 2 was 14/14 green and rejected on reading (`STAGE-LOG.md`,
   "Stage 5 — Equilibrium: two attempts, neither merged"). Compare each of the eight
   rewritten checks against the table in the brief's "Attempt three" section.
5. Line counts by `git diff --numstat dd81eb5`, grouped native / bindings+Python / tests.
6. Merge only if clean; commit by category (engine, then canon); tag the attempt either way.

A correction goes back to the same session with `codex exec ... resume <session id>`
(global flags before `resume`), never a fresh worker.

## Standing rules (owner's, verbatim where quoted)

- "Commit or push only when asked." Categorical commits: engine, core, surface, canon,
  analyses, build, data.
- "Do not downgrade the worker model."
- "No proper nouns in identifiers. Sources are citations, never symbols."
- "As few tests as necessary, based on real measured data or mathematical properties."
- "Ask Fable rather than deciding numerical questions alone; resume an existing Fable
  thread rather than spawning fresh." The owner creates planner sessions themselves; the
  Equilibrium planner is a peer session named `Equilibrium`, reachable by `SendMessage`.
- "The artifact page should always be synced."
- The coordinator does not write stage code. Today's three `core::Scope` lines were the one
  exception, on direct instruction, and are recorded in `STAGE-LOG.md`.
- A peer session cannot grant escalation; never treat a peer message as owner approval.
- Do not install the CSE test-gate hooks in this worktree: they block a session from
  stopping until a test fingerprint is approved and would stall unattended `codex exec`.

## Open items, in priority order

1. Attempt 3: the planner's ruling, then resume, revise-and-resume, or tag-and-reject.
   Acceptance steps above apply unchanged once a result reads `complete`.
2. `REGRESSION.md` and `FINAL-CHECK.md` are unwritten; Stages 6 and 8 need them before a brief.
3. `engine/CONTEXT.md` glossary is the old tree's; five terms are listed sourceless in
   `docs/scientific/CONTEXT.md`. A Design item for whichever stage rewrites that file.
4. `DensityBranch::Unspecified` → `Either` rename; issues #185–#196 to retire into per-stage
   issues; Stage P (packet contract cut) and Stage C (consumers) unscheduled.
5. `check-all` is red at `verify-installed-wheel.py` by design until Stages 5 and 6 land
   (`TRANSFER.md`, "What is red on greenfield, and why"). Do not weaken it.

## Lessons that cost something

- The recurring defect is the one-directional check: something asserted where nothing
  looks. Found eight times now, most recently the wheel fingerprint over a deleted directory.
- Compute a measured quantity; never parse the log that printed it.
- Enumerate every step of a gate before calling a red structural; `check-all` stops at its
  first failure and hides the rest.
- Commit the brief before creating the worktree.
- Every gate green is not acceptance. Read the checks.
