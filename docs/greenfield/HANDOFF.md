# Coordinator handoff

Written 2026-09-03 for a fresh coordinator session. Read this, then `README.md` in this
directory, then `ORCHESTRATION.md`. Everything else is reached from those.

## State at handoff

**Reassessed, owner's ruling 2026-09-03 evening, four steps in order.** A fourth audit established
that Stage 5's acceptance rests on structural and mathematical checks plus the two inherited EOS
activity comparisons: no accepted gate runs the NIST, Cai, Khudaida, Jou or Jakobsen rows through
the new solver, and the check runner retains pass/fail only. The ruling: (1) finish and accept the
Stage 5 repair; (2) a measured-data campaign on the repaired tree before any regression work,
the planner's five oracle cells (ethane saturation against NIST; May 2015 methane–ethane bubble
points, a packet transcription task; Cai 1996 MEA–water at fixed P with T free; Khudaida tie lines
under declared {L,L}, gated on the Stage 4 variant precondition; Jou 1995 pCO2 through the
nine-species bundle with the vapour pinned), each converging with its record asserted and its
deviations reported against a stated bar, a miss being a finding about model or data unless the
trace shows a solver stop; (3) the harness first: every check writes inputs, solved state, record
and deviations to one results file per binary; (4) then Stage 6 reconciled to the repaired
interface. No Stage 6 work until (1)–(3) are done.


- Branch `greenfield` at `b40d5a7` and after, worktree `~/Workspaces/Engineering/ePC-SAFT-greenfield`,
  **nothing pushed** (push only when the owner asks).
- **Stages 1–5 closed.** Stage 5 (Equilibrium) accepted on attempt 4: tags `stage-05-attempt-1`
  to `-4`; the stage worktree `~/.codex/worktrees/gf-stage-05` still exists, branch `stage-05`
  at `d3f3392`, merged. No worker is running.
- Serena is pointed at the greenfield worktree for Claude sessions and at the stage worktree
  for codex workers (`STAGE-LOG.md`, Stage 5 closure); onboarding memories live under
  `.serena/memories/` and are not tracked.
- Overview page: https://claude.ai/code/artifact/1bb42e8c-f3ef-4c72-99e1-335f7614bf58 —
  keep the head sha in its header synced after every commit (owner rule).

## Stage 6, the live process

| item | value |
|---|---|
| worktree | `~/.codex/worktrees/gf-stage-06`, branch `stage-06`, cut from `greenfield@fb47283`, merged with `c5ef4a4` |
| worker | `codex exec`, model `gpt-5.6-sol` (owner's choice), reasoning high, session id `01a069b3-a577-7fd1-94e8-f037660b4535`, launched 2026-09-03 |
| pid | none: stopped twice on 2026-09-03, `source_contradicted` (the scope ruling was not yet in the design documents; fixed at `be3f329`, resumed) then `gate_impossible` (the assembly has no EOS-parameter `F_θ` columns, which `EQUILIBRIUM.md` §12 promised as Stage 5's). **Held until the Stage 5 repair merges**; then merge `greenfield` into `stage-06` and resume session `01a069b3` with a precondition on the new columns |
| log | `attempt-1.log` in the worktree; result lands in `stage-06.result.json`; schema `stage-06.result.schema.json` (untracked, copied from Stage 5's) |
| brief | `docs/greenfield/briefs/stage-06.brief.{md,json}` at `fb47283`, sources pinned at `6bdc05f`, 22 gates, preflight clean, no deferred mutations |
| Serena | `mcp_servers.serena` in `~/.codex/config.toml` points at this worktree; the brief's reading step reports as `serena_reading` |

Acceptance as for Stage 5 (`accept.sh` shape): wipe and rebuild, 22 gates in brief order, 22
mutations with the three static gates' mutants not required to build.

## Stage 5 repair, in flight alongside Stage 6

Eight findings from the 2026-09-03 read-only audit (`STAGE-LOG.md`, "Stage 5 repair"), briefed by
the planner as `briefs/stage-05-repair.brief.{md,json}`: the pinned-amount feed term, the standard
term's temperature derivative evaluated at the iterate, the deficit on dimensionless rows, the
convergence record into the trace, `keep_alive` on `compile_problem`, reject-by-name in the Python
mapping, the two physical defaults removed with rule 12 widened to the subpackage, and (from Stage
6's stop) the EOS-parameter `F_θ` columns. Worktree `~/.codex/worktrees/gf-stage-05-repair`, branch `stage-05-repair` at `be6ddf0`.

| item | value |
|---|---|
| worker | `codex exec`, `gpt-5.6-sol` at high (owner), session id `01a069d6-7868-7ec3-a66e-308fdd5ce3d6`, launched 2026-09-03 |
| pid | 806870, the fourteenth resume (the deficit check on the support-restricted electrolyte cell; log `attempt-1.resume14.log`). Earlier: the thirteenth resume stopped `source_contradicted` with a measurement (the subcooled ternary at 30 bar refuses as `phases_coincide`, deficits at 1e-15, so the deficit identity has no non-trivial point there); planner asked for a vanishing-by-amount cell or a "reported, not asserted" ruling for the deficit. Last open item before acceptance rerun and merge. Earlier: 769465, the thirteenth resume (one check: S15-deficit at a non-trivial point; log `attempt-1.resume13.log`). Acceptance of `4dde17e`: 18 gates green, 17 of 18 mutations red, `active_parameter_column` red on the re-anchored sed, `deficit_identity` blind at a trivial point. Earlier: the twelfth resume reported `complete`, 18 gates green; tree committed at `4dde17e` on `stage-05-repair`; acceptance (`tools/greenfield/accept-stage.sh`) running. Earlier: 662177, the twelfth resume (eleventh stopped `source_contradicted`: a relative step on a zero `k_ij` is no step; the step is now `1e-4 · max(|θ|, 1)`; log `attempt-1.resume12.log`). Earlier: 647446, the eleventh resume (active-parameter cell moved to the ternary flash per the planner; log `attempt-1.resume11.log`). Earlier: the tenth resume stopped `gate_impossible`; the S9 active-parameter cell's baseline, the binary `{V,L}` flash at 180 K and 1 bar, does not converge even with nothing active (S3, S4, S5, S10, S14 green). Planner asked to move the cell to the ternary flash and record the binary non-convergence as a campaign finding; trace `attempt-1.resume10.log`. Earlier: 556031, the tenth resume (ninth stopped `gate_impossible` on an ill-posed fixed-`T` cell, `k = 1, f = 0`; read as `P` free; log `attempt-1.resume10.log`). Earlier: 532753, the ninth resume (eighth stopped `gate_impossible`: the accessor name `active()` collided with rule 10's reserved trace token; renamed `active_parameters()`; log `attempt-1.resume9.log`). Earlier: 511602, the eighth resume (seventh stopped `source_contradicted`: the chain-rule term needs the model's active list and no model write was listed; `helmholtz_model.hpp` now carries one accessor; log `attempt-1.resume8.log`). Earlier: 492541, the seventh resume (log `attempt-1.resume7.log`; carries the corrected finding 8, the absent free-value state, the disabling mutations and the bug-fix clarification). Earlier: the sixth resume stopped `source_contradicted` at the build gate because the coordinator's prompt said to keep `boundary_checks.cpp` "exactly as written" and the file has an undeclared identifier (`reactive_compiled` for `reactive_free`); a bug fix is not compaction, and the next prompt says so. Held for the planner's corrected finding 8 (the packing-coordinate chain rule) so the next resume carries everything at once. Earlier: 458504, the sixth resume (fifth stopped `budget_ceiling` by three lines after the planner's two-part cell; re-budgeted 145/355; log `attempt-1.resume6.log`). Earlier: 428563, the fifth resume (log `attempt-1.resume5.log`; the S9 standard-term check is now a derivative comparison at a point). Earlier: the fourth resume stopped `gate_impossible` (the S9 reactive free-T cell stalls at 200 iterations; the Jacobian matches centred differences, so the state is wrong, not the assembly); waiting on the planner's cell text, then resume session `01a069d6` a fifth time. Earlier: 342270, the fourth resume (third stopped on the coordinator's own field error: two `ceiling_total` values; log `attempt-1.resume4.log`). Earlier: 327197, the third resume (second resume stopped `budget_ceiling` on `boundary_checks.cpp` 327/300; planner re-budgeted 130 per write, 340 total; log `attempt-1.resume3.log`). Earlier: the second resume (first run stopped `source_contradicted` on rule 12's scope; first resume interrupted before writing to carry the planner's empty-structure caveat). Logs `attempt-1.log`, `attempt-1.resume.log`, `attempt-1.resume2.log`. Find codex with `ps -eo pid,args \| grep "codex exec"` |
| log | `attempt-1.log`; result `stage-05-repair.result.json`; schema `stage-05-repair.result.schema.json` (untracked) |
| brief | `briefs/stage-05-repair.brief.{md,json}` at `363ff5d` (amended twice after honest stops), sources pinned at `4fca81d`, 17 gates, preflight clean, 8 mutations bite only after the write |
| Serena | `mcp_servers.serena` in `~/.codex/config.toml` points at this worktree; re-point to `gf-stage-06` before resuming Stage 6 |

Merges before Stage 6 acceptance; then `stage-06` merges `greenfield` and its session resumes with
the `active_parameter_column` precondition already in its brief.

## Next stage: 6 Regression

`REGRESSION.md` exists. Before its brief, in `EQUILIBRIUM.md`/`REGRESSION.md` and the brief:

1. An observation that names several solved targets, so finite-dose heat (two solves) is
   expressible; the `amplification` definition in the diagnostics block (Ipopt's
   `regularization_size` today); re-pin the brief's sources at the current design commit;
   the owner's ruling on `low` line floors and `budget_undershot`: **ruled 2026-09-03, drop
   both** (`CANON.md`, Budget); Stage 6's brief carries ceilings and named checks only.
2. Two Stage 5 deviations to fold in: "skipped by name" reporting when a fixture certificate
   fails; the reaction-direction zero assertion on the feed sensitivity (S9).
3. The brief tells the worker to read the tree through Serena (symbol overview, find symbol,
   references) before writing; `mcp_servers.serena` in `~/.codex/config.toml` points at the
   stage worktree, so re-point it when the Stage 6 worktree is created.
4. Model and effort per the owner: attempt 4 ran on `gpt-5.6-sol` at high; never downgrade.
5. **Scope, ruled by the owner 2026-09-03 (accepted as the planner proposed):** Stage 6 writes the regression shell and
   the `PropertyOrder` first-derivative block on `State` (`EOS.md` Properties §4, both
   coordinate sets) that a `Property` leaf's Jacobian needs. It does not write
   `Thermochemistry`, `physics/ideal.hpp`, the `State` totals, `standard_gibbs(T, order)` or
   the seventh active family: those are EOS-category writes the NIST gate never touches and
   form a **Properties stage** after this one, to be numbered by the owner. Until it lands,
   `Property` leaves name residual fields and the active parameters are the six EOS families.

**Acceptance is the coordinator's, never the worker's.** `accept.sh` (in the session
scratchpad; recreate from `STAGE-LOG.md`'s description if lost): `meson setup --wipe
<build> engine`, compile, every gate in brief order, then every mutation: sed, rebuild,
run, expect red, `git checkout` the file, rebuild. A mutant that does not build, or stays
green, is `mutation_invalid`; sort those into tooling artifact, stale sed and real finding
before reporting, and send the stale ones to the planner for re-anchoring. Then read the
checks against the design's table. Measure any file the worker compacted as first written,
from its log.

## Standing rules (owner's, verbatim where quoted)

- **Standing approval, 2026-09-03** ("I am just going to hit yes every time"): the coordinator
  commits by category, cuts worktrees, launches, interrupts and resumes workers on the model the
  owner last chose, merges a stage that passes acceptance clean, publishes the page, and makes
  the config edits already authorized, all without asking. It still asks before a **push**, before
  **deleting** anything it did not create in the session (worktrees, branches, files), and for any
  **scientific or design ruling** that is the owner's or the planner's. Categorical commits:
  engine, core, surface, canon, analyses, build, data.
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

1. Stage 6 brief (above).
2. `FINAL-CHECK.md` is unwritten; Stage 8 needs it before a brief. `REGRESSION.md` exists;
   before the Stage 6 brief it owes four items from the 2026-09-03 reviewer audit: (a) an
   observation that names several solved targets, so finite-dose heat (two solves) is
   expressible; (b) the `amplification` definition in the diagnostics block, currently Ipopt's
   `regularization_size`; (c) re-pin the brief's sources at the current design commit, since
   the preflight verifies the old pin while the worker reads the merged tree; (d) the owner's
   ruling on `low` line floors and the `budget_undershot` stop now that `CANON.md` says a
   minimum proves nothing.
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
