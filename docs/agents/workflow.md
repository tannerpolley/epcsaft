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
keeps evidence under `build/` or the analysis leaf, never `/tmp`. Current
key-solve work counts are in `build/wheels/main-current/README` of the main
checkout.

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
needs no wheel or binding management. Builds use explicit `-j` parallelism sized
from CPUs and available memory; `EPCSAFT_BUILD_JOBS` can select it for a measured
run. Use one-thread numerical pools only for timing measurements.

## Gates

This is the only gate list. Comparisons select prepared Python interpreters or
build/environment directories, with verified source and loaded-library identities.
Use the main checkout's prepared environment for A and the slice's for B; retain
an unchanged baseline environment before changing its source. Wheels serve
historical studies and merged-main delivery, outside these candidate checks.

| When | Run |
|---|---|
| Every delivery state | `bash tools/check-affected.sh origin/main` and `bash tools/check-all.sh` |
| Native C++ change, final delivery state only | `nice -n 10 tools/check-pr-native --pr N --head <full-sha> --base <merge-base-sha>` (clang-tidy on the native units changed since the merge base and every unit that includes a changed file, read from the lint build's dependency records; any `.clang-tidy`, Meson, subproject, `pyproject.toml`, `uv.lock` or `.codex/environments/` change sweeps every unit, as does `--full` or a Clang or external-header change since the last passing full sweep; plus the existing ASan/UBSan cases in one instrumented tree) |
| Change to the native equilibrium layer or a layer below it | `tools/check-equilibrium-corpora <main environment> <slice environment> --out build/<dir>`: Cameretti, Gross 2002, neutral sweep, #243 matrix and replay in one table; every difference must be the claimed change |
| Before handoff of a change to native equilibrium or parameter regression code | `bash engine/tools/run-native-tests.sh 'equilibrium/Explicit/*' 'regression/Explicit/*'` (the binary128 reactive cases, the reacting Cp fit and the envelope traces); `check-affected.sh` selects them only when their own code changes |
| Any Engine runtime change (`engine/native`, `engine/src`) | `tools/time-key-solves <main environment> <slice environment> --pairs 1 --out build/<dir>` and an A=A control (two workers from the same prepared main environment). Exit 0 is the gate: identical ionic LLE and MEA flash work counts, no increase in the frozen F1 and reacting Cp evaluations' reference computations/columns, density solves or equilibrium iterations, numerical values within 1e-10 of each quantity's scale, and each one's serial CPU at most 1.25 times its frozen value. A slice that intends a work change adds `--expect-change` and reports the printed counts to the owner; this never waives ceilings or missing evidence |
| Change to fit scheduling, caching, optimizer work, observation composition or derivatives, reference computations, equilibrium solve work, or density algorithms | Add `--fit --control <main-A=A-summary.json>` to the same key-solve command. First run main A=A with `--fit` to retain the same-host control. Every F1 fit must finish within 120 s wall using native default workers (32 on the reference host). The report includes the control's wall spread and C_fit = sum over evaluations of the slowest worker's compute CPU. See `tools/time-key-solves.md` for measurement boundaries and fixed values |

Run checks that share the ordinary editable Meson tree sequentially; independent
checks can use `tools/run-together.sh`. After freezing a local candidate commit,
the lead runs the heavy gate, including before publication approval with the
candidate's eventual PR number. The gate locks
`build/pr-gates/pr-N/<head>/` in the main checkout and prepares one detached source,
one locked environment and one Clang debug ASan/UBSan Meson tree. The existing
`check-native.sh --prepared <root> [--changed <base>]` and
`check-native-sanitizers.sh --prepared <root>` own their checks and never build.
The combined command owns preparation and retains their logs and test records.

Passing evidence is reused only for the same head, base, check selection,
toolchain, external headers, locked environment and build flags. Compiler,
header, environment or flag changes replace stale Meson objects; base/selection
changes rerun checks with reusable objects. Failed checks resume incomplete work.
Native tests use `MALLOC_PERTURB_=0` (owner, #428); set a nonzero byte to restore
random heap fill. Ordinary setup and the shared heavy gate both honor
`EPCSAFT_BUILD_JOBS` with explicit Ninja parallelism.

Work counts are the key-solve gate; `--pairs` above 1 is only for a timing
claim, run back to back with the A=A control and judged against its spread.
After merging main, gate again only when the merge touched code the slice
depends on. Reporting-only fixes get focused checks.

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

## Closeout

After each merge, and in a regular sweep, the lead runs
`tools/prune-worktrees WORKTREE ...` (owner, 2026-10-20), naming each task
worktree whose agent task is terminal with no queued turn and whose pull request
is merged or issue closed. The tool refuses a worktree with a process inside it
or an open pull request, archives the rest through `tools/archive-worktree`, and
deletes their local branches and any remote branch whose tip is on
`origin/main`. Unmerged branches are kept unless `--archive-unmerged` is given
and their issue is closed; they are then first tagged `archive/<branch>` on the
exact commit. `--dry-run` plans without changing anything. Gate source
checkouts under `build/pr-gates/` are left to `tools/check-pr-native`.
