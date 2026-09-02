# How the remaining stages get built

Three models, one writer per stage, and a rule about what crosses each boundary.

## Roles

| | does | does not |
|---|---|---|
| **Claude Opus 5** | talks to the user, owns every process and every write to `docs/greenfield/*.md`, invokes `codex exec`, runs acceptance | write code in a stage worktree |
| **Fable 5.1** | plans a stage: reads the tree, writes the stage brief, solves the novel problems | execute, or supervise the executor |
| **Codex 5.6** | writes the stage's code in an isolated worktree, one stage per invocation | plan, or accept its own work |

`codex exec` is driven from Bash, verified working:

```bash
codex exec -C <worktree> -s workspace-write --output-schema <schema> -o result.json "<brief>"
```

Useful flags: `-C` working root, `-s read-only|workspace-write` sandbox, `--output-schema`
a JSON Schema the final response must satisfy, `-o` where to write it, `--json` for a JSONL
event stream, `resume`/`fork` to continue a session, and `codex exec review` for a
fresh-context review.

## Why Opus holds the Codex calls

Fable authors the plan, so it is tempting to let Fable run the implementer. It does not
work here. A subagent reports only through its final response, so if Codex were Fable's
child, Opus could not observe a running stage or steer it — a stage would be a multi-hour
opaque turn ending in a summary of a summary. That is exactly the posture that let Stage
1's wrong source reference through: a faithful implementer copying the wrong file and
reporting success. Opus→Fable by `SendMessage` and Opus→Codex by `resume` are both
parent-to-child, which is the only route the user's agent rules permit anyway.

## The stage brief

The thing that crosses Fable→Opus→Codex is not prose. Its load-bearing claims are
executable, and are executed against the tree **before** Codex starts:

- **sources** as `branch@sha:path` plus content markers that distinguish the wanted version
  from the wrong one
- **writes** as paths with a line estimate *and a ceiling*
- **gates** as runnable commands, each tagged with what it needs (`native`, `python`)
- **a named mutation per gate** — the edit that must make it fail
- **preconditions** as commands that must exit zero

A brief with a gate and no mutation is invalid. So is a gate whose command touches a path
that neither exists nor appears in `writes`.

Codex is additionally told: verify each source and gate yourself before writing; if a claim
is false, stop and say which. `budget_ceiling`, `gate_impossible` and `source_contradicted`
are first-class outcomes, not failures — an implementer that quietly absorbs a wrong plan
is the failure.

## Why this is the answer to "who reviews the planner"

Not another model reading the plan. Stage 1's plan was read carefully and was still wrong
three ways. What catches a wrong plan is running its claims:

| Stage 1 failure | caught by preflight? |
|---|---|
| source named a commit that did not contain the solver | **yes** — `git show <sha>:<path> \| grep -F <marker>` fails, *if* the marker distinguishes it |
| 950 lines budgeted, 1,349 honest | **flagged, not caught** — no check predicts line count; the ceiling and the `budget_ceiling` stop govern the outcome instead |
| gate named a Python path before bindings existed | **yes** — one path check |
| `NDEBUG` erased every assertion, suite passed vacuously | **no** — the code did not exist yet. Only the mutation at acceptance catches this class. |

Two of four, one flagged, one that plan review cannot reach. Stated plainly because a
mechanism claiming four of four would be the thing to distrust.

It has already paid for itself: applying the path check to Stage 2 as written showed its
performance gate could not run, because `measure_floor.py` imports `epcsaft` and reads
`_core._instrumentation()` while `_core` is not built until Stage 3. Corrected in
`TRANSFER.md` before Stage 2 began.

## The `orchestrate` skill

Its mechanism is not usable here: it drives `read_thread`, `wait_threads`, `list_threads`
and "Codex project tasks attached to a managed worktree", which are Codex-desktop tools
absent from this environment, and it merges to the canonical default branch while stage
PRs target `greenfield`. Three of its policies are adopted verbatim, at no cost:

1. exactly one writer and one worktree per stage at a time;
2. a worker's own summary is never acceptance;
3. a correction goes to the same worker, not a replacement.

## Running Codex: the verified invocation

Established by probing before Stage 2, not by reading documentation. Each line
below cost something to find.

```bash
cd <worktree>
codex exec -C . -s workspace-write --add-dir "$HOME/.cache" \
  --output-schema <schema.json> -o <result.json> "$(cat brief.md)" \
  < /dev/null > events.log 2>&1 &
```

| flag / step | why |
|---|---|
| `--add-dir "$HOME/.cache"` | `workspace-write` grants `[workdir, /tmp, $TMPDIR]` and **not** `$HOME`. Meson's CMake dependency lookup writes under `~/.cache`, so without this `find_package(Ceres)` returns `PACKAGE_FOUND is false` and configure dies at `meson.build:26`. Reproduced locally with a read-only HOME; granting `.cmake` alone does not fix it, granting `.cache` alone does. Verified end to end: configure, compile and test all pass inside the sandbox with it. |
| `< /dev/null` | without it `codex exec` blocks on "Reading additional input from stdin..." even when the prompt is an argument, and does nothing until killed. |
| run in background | a stage takes hours; the orchestrator's shell tool caps at ten minutes. Backgrounding also keeps the orchestrator available to the user mid-stage. |
| `-s workspace-write`, never `danger-full-access` | the narrow grant above is sufficient for this build. |

### Worktree setup, per stage

```bash
git worktree add -b stage-N ~/.codex/worktrees/gf-stage-N greenfield
cp -r <greenfield>/engine/build/subprojects <worktree>/engine/build/    # gitignored
cd <worktree> && uv sync                                                # never copy a venv
```

`engine/build/` is gitignored, so a fresh worktree has only `.wrap` files and
Meson would refetch CppAD, HiGHS, Boost and NLopt.

**Never copy a `.venv` between worktrees.** A copied venv keeps the source venv's
absolute path in every script shebang, so a wheel installed through it lands in
the *other* worktree's environment. This project has already lost a working
environment to exactly that, by symlink rather than by copy. `uv sync` builds an
independent one in seconds.

## Environment failure is a stop, not a repair

Add `environment_broken` to the stop conditions, alongside `gate_impossible` and
`source_contradicted`.

The sandbox probe showed why. A code-writing agent that hits "Ceres not found"
mid-stage does not stop -- it helps. The plausible repair from inside that
sandbox is to edit `meson.build` and drop or guard the Ceres dependency, after
which everything configures, compiles, passes the physics checks, and reports
success with the regression backend quietly amputated. That is the Stage 1
failure shape again: a faithful agent producing a confident wrong result because
its environment lied to it.

So the brief must say, and the schema must allow: **an implementer may not modify
`meson.build`'s dependency block, the virtual environment, the subprojects, or
the sandbox invocation to make a build succeed.** A build that fails for
environmental reasons is the orchestrator's problem. Stop and name it.

## Sequencing

Stage 2 runs **by hand**: Opus writes the brief and drives one `codex exec` directly, so
the tooling is built against what actually breaks — sandbox friction, latency, where a
brief was underspecified — rather than against a prediction of it. The scripts
(`preflight.sh`, `run-stage.sh`, `accept-stage.sh`, two schemas, ~350 lines under
`tools/greenfield/`) are written after Stage 2, from that evidence.

## Signals that this is not working

- a brief with no content markers, or one gate — the checks have gone vacuous
- more than two `resume` cycles on a stage — Codex is stopping on trivia
- `result.json` claiming values that acceptance does not reproduce — the script is the
  acceptance, the mismatch is the finding
- two stages in a row with zero deviations and zero preflight findings — either the
  planner became perfect or the checks stopped checking; assume the second
