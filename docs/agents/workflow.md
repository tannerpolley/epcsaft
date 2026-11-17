# Leading a slice

One accountable lead owns each slice; one writer per file.

## Before the build

Separate unresolved scientific choices from specified implementation before
dependent coding. Scope readiness to the slice being built: its equations,
consumed inputs and falsifying checks come before implementation. Line
allowances follow [code budgets](../scientific/code-budgets.md). Size a slice to
deliver within a few rounds; split a growing issue into child issues that merge
independently. Each PR closes its issue; work left over gets a new issue before
the merge ([roadmap skill](../../.agents/skills/manage-github-roadmap/SKILL.md#handoff-and-completion)).

**Merge current main at the start of every build round.**

## Builders

A builder is the sole writer in its worktree and starts no subagents. It may
implement and qualify a clear local defect fix of about 20 lines or fewer with an
obvious numerical falsifier. It stops and reports to the lead before an
unresolved scientific choice, a new solver, derivative or method design, a
hard-ceiling excess or a scope change. It never loosens a tolerance or adds
per-case epsilons, application selectors or outer retries; never pushes, opens
pull requests, merges or edits issue bodies or labels (comments are fine); and
keeps evidence under `build/` or the analysis leaf, never `/tmp`. A promoted
study's producer source, unique inputs and run command live in its analysis
leaf; `build/` holds scratch output and logs only. The current main wheel's
key-solve comparison is listed under `numerical_evidence` in the main checkout's
`build/wheels/main-<short-commit>/release.json`.

## Review tiers

- **Small change** (about 20 runtime lines or fewer, no new equation, criterion
  or physical claim): one independent review after the build.
- **Substantive change** (new mathematics, solver or derivative behavior):
  readiness review before the build and delivered review after it.
- **Reporting-only fixes** (wording, counts, labels; no numerical change): the
  lead verifies them against the finding without another review round.

Reviewers apply `CODING_STANDARDS.md`.

## Development method

Develop with `nice -n 10 tools/native-run <program> [case]`, the native test
programs, or small C++ drivers that write CSV. Python declares inputs and reads
results. The runner prepares missing dependencies automatically; native iteration
needs no wheel or binding management. Use one-thread numerical pools only for
timing measurements.

Adapt justified existing native fitting and evaluator logic through the plan's
adoption gate, and implement only observation forms consumed by current scope;
the parameter regression formulation's broader algebra defines meaning rather
than requiring unused nodes.

Before adding a derivative, iteration, or action, name the returned quantity,
numerical step, or diagnostic that uses its result. Run the calculation only
when that consumer needs it (owner, 2026-11-13).

`CODING_STANDARDS.md` owns revalidating the mathematics with a scientific change
and treating a gated result as a diagnosis; run the owning issue's smallest
independent falsifier before acceptance (examples in ALG-TERMINATION
"Verification meaning").

## Gates

This is the only gate list. Verify internal Engine changes locally with `tools/native-run`, the native `*_checks` programs, or a small C++ driver. Run the case that exercises the change and keep the issue's numerical acceptance evidence. After several changes have accumulated and a downstream application needs a new wheel, the lead builds the Python extension and wheel once for that batch. Comparisons select prepared Python interpreters or build/environment directories, with verified source and loaded-library identities. Use the previous qualified batch environment for A and the frozen candidate environment for B. Retain an unchanged baseline environment before changing its source.

| When | Run |
|---|---|
| Every delivery state | `bash tools/check-affected.sh origin/main`: run the local source checks and affected C++ tests. The command prints the applicable checks for the lead's batch |
| Lead's downstream wheel batch, once on the frozen candidate | `agent-heavy tools/check-pr-native --pr N --head <full-sha> --base <previous-qualified-sha>` (clang-tidy on the native units changed since the previous qualified batch and every unit that includes a changed file, read from the lint build's dependency records; any `.clang-tidy`, Meson, subproject, `pyproject.toml`, `uv.lock` or `.codex/environments/` change sweeps every unit, as does `--full` or a Clang or external-header change since the last passing full sweep; plus the existing ASan/UBSan cases in one instrumented tree) |
| Lead's batch, when it contains a change to native equilibrium or a layer below it | `tools/check-equilibrium-corpora <prepared-A> <prepared-B> --out build/<dir>`: Cameretti, Gross 2002, neutral sweep, #243 matrix and replay in one table; every difference must be the claimed change |
| Before handoff of a change to native equilibrium or parameter regression code | `bash engine/tools/run-native-tests.sh 'equilibrium/Explicit/*' 'regression/Explicit/*'` (the binary128 reactive cases, the charged PH/PS case, the reacting Cp fit and the envelope traces); `check-affected.sh` selects them only when their own code changes |
| Lead's batch, when it contains a change to Engine runtime code (`engine/native` or `engine/src`) | `tools/time-key-solves <prepared-A> <prepared-B> --pairs 1 --out build/<dir>` and an A=A control (two workers from A). Exit 0 is the gate; `tools/time-key-solves.md` owns its work-count, value and CPU limits. Use `--expect-change` only when a slice in the batch declares an intended change in native work. Report the printed A and B counts to the owner alongside the intended change and its slice. Every other count difference still needs an explanation. The option preserves the existing value limits, CPU ceilings, and requirements for complete and repeatable evidence. If a batch comparison fails, narrow it down commit by commit using prepared environments for the revisions compared |
| Lead's batch, when it contains a change to fit scheduling, caching, optimizer work, observation composition or derivatives, reference computations, equilibrium solve work, or density algorithms | Add `--fit --control <A=A-summary.json>` to the same key-solve command. First run A=A with `--fit` from the previous qualified batch environment to retain the same-host control. `tools/time-key-solves.md` owns the fit limits and measurement boundaries |

Run jobs that share a Meson tree sequentially. At the downstream wheel batch, freeze the candidate and run `bash tools/check-affected.sh --python <previous-qualified-sha>` for the Python suite and the architecture checks that need linked images. Run `check-pr-native` once for the frozen candidate before delivering its wheel. It uses a separate instrumented tree under `build/pr-gates/pr-N/<head>/` and retains its logs and check records. Publish the wheel after the applicable batch checks pass.
Native tests use `MALLOC_PERTURB_=0` (owner, #428); set a nonzero byte to restore
random heap fill.

Work counts are the key-solve gate; `--pairs` above 1 is only for a timing
claim, run back to back with the A=A control and judged against its spread.
Rerun the affected native cases if the merge changed code used by the slice. If the frozen batch candidate changes, check the candidate used for wheel delivery. Reporting-only fixes get focused checks.

F1 and the reacting Cp row run once per environment in the default command. Full fits are explicit
because they repeat residual and Jacobian evaluations across optimizer steps;
input-mapping-only edits use the default command unless they change this work.
The fixed fixture, CPU ceiling and full-fit measurements are owned by the
key-solve command's manual, `tools/fixtures/mea-f1-identity.json` and `tools/fixtures/mea-cp-identity.json`.

## Acceptance

The lead checks reviewer findings against source and evidence; a build,
readiness label or reviewer opinion is not numerical acceptance. When a
capability claim changes, rewrite the issue body's current state and update
repository evidence together.

Fix before a negative disposition (owner decision 2026-09-23): no analysis,
qualification or issue closes on a result that is refused, failed, unavailable,
reference unresolved, tolerance gated or a supported negative until its cause
is diagnosed and a root-cause fix has been attempted and shown impossible, or
the state is shown outside the stated physical domain, with evidence. Reviews
ask whether the outcome could be fixed, not only whether its label is honest.
Missing accuracy/runtime thresholds stay owner decisions.

## Closeout

After each merge, the lead runs `tools/prune-worktrees WORKTREE ...` (owner,
2026-10-20) on each task worktree whose agent task is terminal with no queued
turn and whose pull request is merged or issue closed. Report dirty or active
worktrees and leave them in place. The daily `epcsaft-clean-build.timer` is
disabled. For a manual cleanup, the lead reviews `tools/clean-build --dry-run`
from the main checkout before running `tools/clean-build`. For each entry in
`build/evidence-to-promote.md`, the lead commits its exact bytes to the analysis
or deletes it and records the loss there by issue number. Once a task is
finished, the lead commits or discards the changes in its listed worktree in
`build/worktrees-to-resolve.md`, then prunes it.
