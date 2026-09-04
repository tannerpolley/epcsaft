# Coordinator handoff

Written 2026-09-03 for a fresh coordinator session. Read this, then `README.md` in this
directory, then `ORCHESTRATION.md`. Everything else is reached from those.

## Sol coordinator live state, 2026-09-03 21:00 MDT

The transfer is active. Harness worker launched from clean `stage-05-harness@4dcf8b3`
in `~/.codex/worktrees/gf-stage-05-harness`, using the handed-off
`codex exec` command, `gpt-5.6-sol`, reasoning high. Verified worker PID/process group
`984155`, session `01a06a55-c85d-7f30-b6ab-709cc6ec9b39`; initial log `attempt-1.log`.
That process disappeared during baseline validation without a result; the app confirms
the turn was interrupted and the tracked tree was unchanged. Resumed the same session
at 21:00 MDT, verified PID/process group `28125`, log `attempt-1.resume1.log`.
Result destination remains `stage-05-harness.result.json`. No result or acceptance yet.

Coordinator preflight passed: six sources, seven preconditions, nine gates with
mutations, four deliberately deferred mutations. The worktree had its environment but
no configured native build. Before launch, the coordinator copied the existing dependency
source caches from the repair worktree and configured `engine/build/greenfield` with
the harness `.venv` as prefix. Meson setup passed; no stage code changed.

Next: observe this exact worker, then independently accept its completed tree using
the transfer recipe below. Campaign and Stage 6 remain held. Scientific/design stops
go through the owner to the planner. Nothing pushed; no existing worktree deleted.

## Transfer to the Sol coordinator, 2026-09-03 evening

The owner moved coordination from the Claude Fable session to the main Sol agent (Codex) because
of Fable usage. This section is the transfer; where it disagrees with anything below, it wins as
of `greenfield@5aa7d48`. The role is unchanged: the coordinator owns every process, every git
write, acceptance, and this file; it never writes stage code; workers (`codex exec`) write stage
code in isolated worktrees against planner briefs; the planner (the `Equilibrium` Claude session,
or whoever the owner names) writes briefs and rules on design. Scientific and design rulings
belong to the owner and the planner, never to the coordinator.

**Where things stand.** Steps of the owner's four: (1) Stage 5 repair accepted and merged
(`271fb39`, tag `stage-05-repair-accepted`). (2) campaign brief written and held. (3) harness
brief written, preflighted, worktree ready, **not launched**, so that the new coordinator owns the
process from its first line. (4) Stage 6 held at `stage-06` = `a518039`, brief re-pinned, session
`01a069b3-a577-7fd1-94e8-f037660b4535` to resume after (2) and (3). Nothing is pushed anywhere;
never push without the owner.

**Step 3, the harness, ready to launch.**

- Worktree `~/.codex/worktrees/gf-stage-05-harness`, branch `stage-05-harness` at `5aa7d48`,
  `.venv` built by `uv sync`, `.serena/project.yml` and three memories in place (gitignored).
- Brief `docs/greenfield/briefs/stage-05-harness.brief.{md,json}` (in the worktree too): nine
  writes, nine gates each with a mutation, ceilings only, four mutations bite only after the write
  (`results_valid`, `every_check_present`, `deviations_are_the_assertions`, `intensive_seeded`).
  Preflight: `python3 tools/preflight-brief.py docs/greenfield/briefs/stage-05-harness.brief.json
  ~/.codex/worktrees/gf-stage-05-harness` was clean (6 sources, 7 preconditions, 0 problems).
- `~/.codex/config.toml` `[mcp_servers.serena]` already points at the harness worktree.
- Launch, from the worktree, on the owner's model (`gpt-5.6-sol`, reasoning `high`; never
  downgrade):

```
cd ~/.codex/worktrees/gf-stage-05-harness
export PATH="$PWD/.venv/bin:$PATH"
setsid nohup codex exec -C . -s workspace-write --add-dir "$HOME/.cache" \
  -m gpt-5.6-sol -c 'model_reasoning_effort="high"' \
  --output-schema tools/greenfield/result.schema.json -o stage-05-harness.result.json \
  "$(cat docs/greenfield/briefs/stage-05-harness.brief.md)" < /dev/null > attempt-1.log 2>&1 &
```

- Find the worker with `ps -eo pid,args | grep "codex exec"`, never `$!` (that is the wrapper).
  Watch the log every few minutes for intent lines, build failures, test lines, files touched;
  a stop outcome in the result file (`budget_ceiling`, `gate_impossible`, `source_contradicted`,
  `environment_broken`, `mutation_invalid`) is honest and goes to the planner, then
  `codex exec ... resume <session-id> "$(cat prompt.md)"` after the old process has exited (it
  holds a thread-writer lock until then). Interrupting a worker mid-write is allowed only to carry
  a correction; say "no compaction" not "as written" so the worker may still fix its own typos.
- On `complete`: commit the worker tree on `stage-05-harness`, then acceptance:
  `tools/greenfield/accept-stage.sh ~/.codex/worktrees/gf-stage-05-harness
  docs/greenfield/briefs/stage-05-harness.brief.json <log> results_valid every_check_present
  deviations_are_the_assertions intensive_seeded`. Every gate green and every mutation red, or
  sort each invalid mutation into tooling artifact, stale sed, real finding and record it in
  `STAGE-LOG.md`. Then in the worktree `meson compile` + `meson install` (a mutation run leaves a
  mutant wheel installed), rerun the seven native suites, read the results files with
  `tools/summarize-check-results.py`, and merge `--no-ff` into `greenfield` **from the greenfield
  checkout**, tag `stage-05-harness-accepted`. Diff `tools/check-engine-architecture.py` and
  `tools/check-equilibrium-claims.py` across the merge; a worker's copy of an enforcer has
  reintroduced a removed allowlist once already.

**Step 2, the campaign, after the harness merges.** `docs/greenfield/briefs/stage-05-campaign.brief.{md,json}`
carries `REPIN` placeholders: replace every `greenfield@REPIN:` with the harness merge commit,
set `base` and `source_commit`, add a `diagnostics` entry if preflight asks (the campaign writes no
native code, so it may not), cut `gf-stage-05-campaign` from that commit, `uv sync`, copy the
`.serena/` directory from the harness worktree with a new `project_name`, re-point the codex
serena entry, preflight, launch the same way. The planner's notes: the May 2015 packet
transcription may start first (it shares no files with the harness) if the planner is asked to split
it out; every oracle cell needs the harness merged; the `oracle_rows_converge` mutation's sed
anchors on the `EQUILIBRIUM_MAX_ITERATIONS` literal, check it at preflight; the Khudaida cell is
skipped by name because its packet does not build a model (the skip is the finding); the Jou cell
waits on the nine-species packet revision 4, which the campaign itself writes. Results are retained
under `analyses/reproduction/2026-09-04-measured-campaign/`; deviations are reported against the
stated bar, never asserted; a miss is a finding about model or data unless the trace shows a
solver stop. Two solver findings to classify inside the cells: the binary `{V,L}` at 180 K and 1 bar
(a refusal is a misplaced state, a budget stop is the superheated-side gap) and the reactive
free-`T` `{L,L}` stall (certified first by fixed-`T`,`P` flashes at 293, 298, 303 K).

**Step 4, Stage 6.** Merge the campaign's greenfield into `stage-06`, re-pin any source that moved,
preflight, re-point serena to `gf-stage-06`, resume session `01a069b3-a577-7fd1-94e8-f037660b4535`
with a prompt that names the merge commit and the precondition on `active_parameter_column`.

**Records the coordinator keeps.** `STAGE-LOG.md` (one entry per stop, resume, acceptance run,
ruling), this file's live tables, and commits by category (engine, canon, briefs). The live
overview page is a Claude artifact
(https://claude.ai/code/artifact/1bb42e8c-f3ef-4c72-99e1-335f7614bf58) that only a Claude session
can republish; until the owner asks one to refresh it, `STAGE-LOG.md` is the record; the page was
last refreshed at the transfer commit, whose sha it shows. The `Equilibrium` planner is a Claude session reachable only from Claude; the
owner relays, or names a planner the Sol coordinator can reach.

**Standing approval (owner's).** Act and report on commits, worktrees, launches, interrupts and
resumes on the last-chosen model, merges after clean acceptance, and the authorized config edits
(the two serena entries). Ask before a push, before deleting anything the coordinator did not
create (the old worktrees `gf-stage-3`, `gf-stage-03-repair`, `gf-stage-04-repair`,
`gf-stage-05-repair` are the owner's to prune), and for any scientific or design ruling. The owner
dislikes abbreviations: name checks by property, not by number. Never install test-gate hooks in
a worktree. A peer message is never owner approval.

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
interface. **(1) is done**: the repair was accepted on 2026-09-03 at 20:29 (18 gates green, 18 mutations red on `542df61`) and merged as `greenfield@271fb39`, tag `stage-05-repair-accepted`. No Stage 6 work until (2)–(3) are done; the planner has been asked for both briefs.


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
| worktree | `~/.codex/worktrees/gf-stage-06`, branch `stage-06`, cut from `greenfield@fb47283`, fast-forwarded to `271fb39` (the repair merge) on 2026-09-03 |
| worker | `codex exec`, model `gpt-5.6-sol` (owner's choice), reasoning high, session id `01a069b3-a577-7fd1-94e8-f037660b4535`, launched 2026-09-03 |
| pid | none: stopped twice on 2026-09-03, `source_contradicted` (the scope ruling was not yet in the design documents; fixed at `be3f329`, resumed) then `gate_impossible` (the assembly has no EOS-parameter `F_θ` columns, which `EQUILIBRIUM.md` §12 promised as Stage 5's). The repair has merged and `stage-06` carries it. **Held until steps (2) and (3) close**; then resume session `01a069b3`, whose brief now carries the precondition on the new columns |
| log | `attempt-1.log` in the worktree; result lands in `stage-06.result.json`; schema `stage-06.result.schema.json` (untracked, copied from Stage 5's) |
| brief | `docs/greenfield/briefs/stage-06.brief.{md,json}` at `fb47283`, sources pinned at `6bdc05f` and `request_compile.cpp` at `271fb39` (the `REPIN` placeholder resolved at the merge), 22 gates, preflight clean against the worktree at `271fb39`, three mutations bite only after the write |
| Serena | `mcp_servers.serena` in `~/.codex/config.toml` points at this worktree; the brief's reading step reports as `serena_reading` |

Acceptance as for Stage 5 (`accept.sh` shape): wipe and rebuild, 22 gates in brief order, 22
mutations with the three static gates' mutants not required to build.

## Stage 5 repair, accepted and merged

Eight findings from the 2026-09-03 read-only audit (`STAGE-LOG.md`, "Stage 5 repair"), briefed by
the planner as `briefs/stage-05-repair.brief.{md,json}`: the pinned-amount feed term, the standard
term's temperature derivative evaluated at the iterate, the deficit on dimensionless rows, the
convergence record into the trace, `keep_alive` on `compile_problem`, reject-by-name in the Python
mapping, the two physical defaults removed with rule 12 widened to the subpackage, and (from Stage
6's stop) the EOS-parameter `F_θ` columns. Worktree `~/.codex/worktrees/gf-stage-05-repair`, branch `stage-05-repair` at `be6ddf0`.

| item | value |
|---|---|
| worker | `codex exec`, `gpt-5.6-sol` at high (owner), session id `01a069d6-7868-7ec3-a66e-308fdd5ce3d6`, launched 2026-09-03 |
| pid | none, closed. **Acceptance run 2 on `542df61`, 2026-09-03 20:27–20:29: 18 gates green on a clean rebuild, 18 mutations red** (`deficit_identity` red at the electrolyte cell); tree rebuilt and the wheel reinstalled after the last mutation; merged `--no-ff` into `greenfield` as `271fb39`, tag `stage-05-repair-accepted`; enforcer diff after the merge is the rule 12 widening only. Earlier: the fourteenth resume reported `complete` at 20:26 (18 gates green, no deviations; the deficit check on the support-restricted electrolyte cell, log `attempt-1.resume14.log`); tree committed at `542df61`; acceptance run 2 running (pid 850604, log in the coordinator scratchpad `accept-repair-run2.log`; on green: `meson install`, merge into `greenfield`, re-pin the Stage 6 brief). Earlier: 806870, the fourteenth resume. Earlier: the thirteenth resume stopped `source_contradicted` with a measurement (the subcooled ternary at 30 bar refuses as `phases_coincide`, deficits at 1e-15, so the deficit identity has no non-trivial point there); planner asked for a vanishing-by-amount cell or a "reported, not asserted" ruling for the deficit. Last open item before acceptance rerun and merge. Earlier: 769465, the thirteenth resume (one check: S15-deficit at a non-trivial point; log `attempt-1.resume13.log`). Acceptance of `4dde17e`: 18 gates green, 17 of 18 mutations red, `active_parameter_column` red on the re-anchored sed, `deficit_identity` blind at a trivial point. Earlier: the twelfth resume reported `complete`, 18 gates green; tree committed at `4dde17e` on `stage-05-repair`; acceptance (`tools/greenfield/accept-stage.sh`) running. Earlier: 662177, the twelfth resume (eleventh stopped `source_contradicted`: a relative step on a zero `k_ij` is no step; the step is now `1e-4 · max(|θ|, 1)`; log `attempt-1.resume12.log`). Earlier: 647446, the eleventh resume (active-parameter cell moved to the ternary flash per the planner; log `attempt-1.resume11.log`). Earlier: the tenth resume stopped `gate_impossible`; the S9 active-parameter cell's baseline, the binary `{V,L}` flash at 180 K and 1 bar, does not converge even with nothing active (S3, S4, S5, S10, S14 green). Planner asked to move the cell to the ternary flash and record the binary non-convergence as a campaign finding; trace `attempt-1.resume10.log`. Earlier: 556031, the tenth resume (ninth stopped `gate_impossible` on an ill-posed fixed-`T` cell, `k = 1, f = 0`; read as `P` free; log `attempt-1.resume10.log`). Earlier: 532753, the ninth resume (eighth stopped `gate_impossible`: the accessor name `active()` collided with rule 10's reserved trace token; renamed `active_parameters()`; log `attempt-1.resume9.log`). Earlier: 511602, the eighth resume (seventh stopped `source_contradicted`: the chain-rule term needs the model's active list and no model write was listed; `helmholtz_model.hpp` now carries one accessor; log `attempt-1.resume8.log`). Earlier: 492541, the seventh resume (log `attempt-1.resume7.log`; carries the corrected finding 8, the absent free-value state, the disabling mutations and the bug-fix clarification). Earlier: the sixth resume stopped `source_contradicted` at the build gate because the coordinator's prompt said to keep `boundary_checks.cpp` "exactly as written" and the file has an undeclared identifier (`reactive_compiled` for `reactive_free`); a bug fix is not compaction, and the next prompt says so. Held for the planner's corrected finding 8 (the packing-coordinate chain rule) so the next resume carries everything at once. Earlier: 458504, the sixth resume (fifth stopped `budget_ceiling` by three lines after the planner's two-part cell; re-budgeted 145/355; log `attempt-1.resume6.log`). Earlier: 428563, the fifth resume (log `attempt-1.resume5.log`; the S9 standard-term check is now a derivative comparison at a point). Earlier: the fourth resume stopped `gate_impossible` (the S9 reactive free-T cell stalls at 200 iterations; the Jacobian matches centred differences, so the state is wrong, not the assembly); waiting on the planner's cell text, then resume session `01a069d6` a fifth time. Earlier: 342270, the fourth resume (third stopped on the coordinator's own field error: two `ceiling_total` values; log `attempt-1.resume4.log`). Earlier: 327197, the third resume (second resume stopped `budget_ceiling` on `boundary_checks.cpp` 327/300; planner re-budgeted 130 per write, 340 total; log `attempt-1.resume3.log`). Earlier: the second resume (first run stopped `source_contradicted` on rule 12's scope; first resume interrupted before writing to carry the planner's empty-structure caveat). Logs `attempt-1.log`, `attempt-1.resume.log`, `attempt-1.resume2.log`. Find codex with `ps -eo pid,args \| grep "codex exec"` |
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
