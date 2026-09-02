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
