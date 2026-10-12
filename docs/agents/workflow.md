# Leading a slice

One accountable lead owns each slice; one writer per file.

## Before the build

Separate unresolved scientific choices from specified implementation before
dependent coding. Scope readiness to the slice being built: its equations,
consumed inputs and falsifying checks come before implementation. Line
allowances follow [code budgets](../scientific/code-budgets.md). Size a slice to
deliver within a few rounds; split a growing issue into child issues that merge
independently.

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

## Gates

This is the only gate list. "Main wheel" is `build/wheels/main-current/` in the
main checkout; build the slice's wheel with `bash tools/build-engine.sh --wheel
build/wheels/<name>`.

| When | Run |
|---|---|
| Every delivery state | `bash tools/check-affected.sh origin/main` and `bash tools/check-all.sh` |
| Native C++ change, final delivery state only | `bash engine/tools/check-native.sh` (clang-tidy) and `bash engine/tools/check-native-sanitizers.sh`; earlier rounds may lint changed files only with `check-native.sh --changed` once #428 adds it |
| Change to the native equilibrium layer or a layer below it | `tools/check-equilibrium-corpora <main wheel> <slice wheel> --out build/<dir>`: Cameretti, Gross 2002, neutral sweep, #243 matrix and replay in one table; every difference must be the claimed change |
| Before handoff of a change to native equilibrium or regression code | `bash engine/tools/run-native-tests.sh 'equilibrium/Explicit/*' 'regression/Explicit/*'` (the binary128 reactive cases, the reacting Cp fit and the envelope traces); `check-affected.sh` selects them only when their own code changes |
| Any Engine runtime change (`engine/native`, `engine/src`) | `tools/time-key-solves <main wheel> <slice wheel> --pairs 1 --out build/<dir>` and an A=A control (the main wheel against a copy of itself). Exit 0 is the gate: identical ionic LLE and MEA flash work counts and returned states within 1e-10 of each quantity's scale. A slice that intends a work change adds `--expect-change` and reports the printed counts to the owner |

Work counts are the key-solve gate; `--pairs` above 1 is only for a timing
claim, run back to back with the A=A control and judged against its spread.
After merging main, gate again only when the merge touched code the slice
depends on. Reporting-only fixes get focused checks.

## Acceptance

The lead checks reviewer findings against source and evidence; a build,
readiness label or reviewer opinion is not numerical acceptance. When a
capability claim changes, rewrite the issue body's current state and update
repository evidence together.
