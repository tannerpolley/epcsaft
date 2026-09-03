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

Not another model reading the plan. Stage 4's plan was read carefully and was still wrong
three ways. What catches a wrong plan is running its claims:

| Stage 4 failure | caught by preflight? |
|---|---|
| source named a commit that did not contain the solver | **yes** — `git show <sha>:<path> \| grep -F <marker>` fails, *if* the marker distinguishes it |
| 950 lines budgeted, 1,349 honest | **flagged, not caught** — no check predicts line count; the ceiling and the `budget_ceiling` stop govern the outcome instead |
| gate named a Python path before bindings existed | **yes** — one path check |
| `NDEBUG` erased every assertion, suite passed vacuously | **no** — the code did not exist yet. Only the mutation at acceptance catches this class. |

Two of four, one flagged, one that plan review cannot reach. Stated plainly because a
mechanism claiming four of four would be the thing to distrust.

It has already paid for itself: applying the path check to Stage 4 as written showed its
performance gate could not run, because `measure_floor.py` imports `epcsaft` and reads
`_core._instrumentation()` while `_core` is not built until Stage 3. Corrected in
`TRANSFER.md` before Stage 4 began.

## The `orchestrate` skill

Its mechanism is not usable here: it drives `read_thread`, `wait_threads`, `list_threads`
and "Codex project tasks attached to a managed worktree", which are Codex-desktop tools
absent from this environment, and it merges to the canonical default branch while stage
PRs target `greenfield`. Three of its policies are adopted verbatim, at no cost:

1. exactly one writer and one worktree per stage at a time;
2. a worker's own summary is never acceptance;
3. a correction goes to the same worker, not a replacement.

## Running Codex: the verified invocation

Established by probing before Stage 4, not by reading documentation. Each line
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
| `PATH` must lead with the worktree's `.venv/bin` | Meson re-resolves `python3` during *regeneration*, which any `meson.build` edit triggers — and a stage always edits it to add its sources. Configure-time PATH is not remembered. Without this, regeneration finds `/usr/bin/python3`, `import pybind11` fails, and the build dies at `meson.build:33` **after** the stage has written all its files. Stage 4 hit exactly this. |

### Worktree setup, per stage

**Commit the brief pair before creating the worktree.** A worktree is branched from
`greenfield` as it stands, so a brief committed afterwards is not in the checkout. The worker
reads `briefs/stage-NN.brief.json` from its own tree and stops `environment_broken` when it is
missing — correctly, since it cannot verify a single source marker without it. Stage 03.repair
lost a run to this.

```bash
git add docs/greenfield/briefs/stage-NN.brief.{md,json} && git commit   # first
git worktree add -b stage-NN ~/.codex/worktrees/gf-stage-NN greenfield
cp -r <greenfield>/engine/build/subprojects <worktree>/engine/build/    # gitignored
cd <worktree> && uv sync                                                # never copy a venv
```

If the worktree already exists and the worker wrote nothing, `git reset --hard <brief commit>`
is enough; the venv and build directory are untracked and survive it.

`NN` is the two-digit stage number from `TRANSFER.md`'s execution table, and it is the same
`NN` as the brief pair `briefs/stage-NN.brief.{md,json}`. One number names the stage, its
brief, its branch and its worktree; nothing else numbers a stage.

`engine/build/` is gitignored, so a fresh worktree has only `.wrap` files and
Meson would refetch CppAD, HiGHS, Boost and NLopt.

**Never copy a `.venv` between worktrees.** A copied venv keeps the source venv's
absolute path in every script shebang, so a wheel installed through it lands in
the *other* worktree's environment. This project has already lost a working
environment to exactly that, by symlink rather than by copy. `uv sync` builds an
independent one in seconds.

## The stop condition earned itself on its first use

Stage 4 stopped with `environment_broken` after writing all eight of its files. It had
verified all 36 source markers, run every precondition, and confirmed the baseline built
and tested clean; then `meson compile` failed on the pybind11 lookup above. It did not
touch the dependency block.

That is the whole point. The plausible repair from inside that sandbox is to delete or
guard the pybind11 detection, after which the archive builds, the physics test passes, the
architecture checker is green, and the extension is silently broken — a stage that reports
success with the Python bindings amputated. The rule cost four lines in the brief and
bought a result that would have been wrong in a way no gate in the list would have caught.

The correction went back to the *same* session with `codex exec resume <id>`, not to a
fresh agent: the environment fix is the orchestrator's, the compile errors it then exposed
are the implementer's, and the implementer still has its own context for them.

## Environment failure is a stop, not a repair

Add `environment_broken` to the stop conditions, alongside `gate_impossible` and
`source_contradicted`.

The sandbox probe showed why. A code-writing agent that hits "Ceres not found"
mid-stage does not stop -- it helps. The plausible repair from inside that
sandbox is to edit `meson.build` and drop or guard the Ceres dependency, after
which everything configures, compiles, passes the physics checks, and reports
success with the regression backend quietly amputated. That is the Stage 4
failure shape again: a faithful agent producing a confident wrong result because
its environment lied to it.

So the brief must say, and the schema must allow: **an implementer may not modify
`meson.build`'s dependency block, the virtual environment, the subprojects, or
the sandbox invocation to make a build succeed.** A build that fails for
environmental reasons is the orchestrator's problem. Stop and name it.

## Sequencing

Stage 4 runs **by hand**: Opus writes the brief and drives one `codex exec` directly, so
the tooling is built against what actually breaks — sandbox friction, latency, where a
brief was underspecified — rather than against a prediction of it. The scripts
(`preflight.sh`, `run-stage.sh`, `accept-stage.sh`, two schemas, ~350 lines under
`tools/greenfield/`) are written after Stage 4, from that evidence.

## Signals that this is not working

- a brief with no content markers, or one gate — the checks have gone vacuous
- more than two `resume` cycles on a stage — Codex is stopping on trivia
- `result.json` claiming values that acceptance does not reproduce — the script is the
  acceptance, the mismatch is the finding
- two stages in a row with zero deviations and zero preflight findings — either the
  planner became perfect or the checks stopped checking; assume the second

## The result schema

`--output-schema` is enforced strictly: **every** property of every object must appear in
its `required` array, and every object needs `additionalProperties: false`. A schema that
omits one optional field is rejected by the API with `invalid_json_schema` before the
worker starts, which costs a launch. The stage schema is
`stage-NN.result.schema.json` in the worktree, alongside the brief:

```json
{"outcome": "complete | budget_ceiling | gate_impossible | source_contradicted | environment_broken",
 "files_written": [{"path": "...", "lines": 0}],
 "gates": [{"id": "...", "ran": true, "passed": true, "note": "..."}],
 "deviations": ["..."],
 "stop_detail": "..."}
```

## Resuming a stopped stage

A correction goes to the same session, never to a replacement worker. The global options
come **before** the subcommand; `codex exec resume -C .` is rejected with
`unexpected argument '-C' found`:

```bash
codex exec -C . -s workspace-write --add-dir "$HOME/.cache" \
  --output-schema <schema.json> -o <result.json> \
  resume <session-uuid> "$(cat resume-prompt.md)" \
  < /dev/null > resume.log 2>&1 &
```

The session id is in the launch log's first lines (`session id: <uuid>`). Update the brief
in the repository first, copy it into the worktree, and have the resume prompt say what
changed and why — the worker still holds its own context for everything else, which is the
whole reason a resume beats a fresh agent.

## One commit, one category

A commit changes one category. Not one file, not one stage — **one category**, from
`CATEGORIES.md`: core, surface, engine, canon or analyses. A commit that lands a design
document and an engine fix together makes the engine fix invisible, and a diff of the week
then reads as thousands of green lines with no way to tell which of them anyone has to
review as code.

This is enforced by discipline rather than by a hook, because the exceptions are real: a
stage's worker writes engine and canon together (the code and its checks), and a rename
touches everything. When a commit must cross a category, its message says which categories
and why.

`tools/diffstat.py` reports any range by **language** and by **category**, and crosses the
two:

```bash
tools/diffstat.py                  # since the tree was emptied — the default, and the right question
tools/diffstat.py --vs-main        # the whole rebuild, the 45,000-line deletion included
tools/diffstat.py --code           # executable code only: C++, Python, shell, build
tools/diffstat.py --files          # every file, classified
tools/diffstat.py HEAD~5..HEAD     # any range
```

It also splits the code additions into **statements, comments and blank lines**. No
off-the-shelf counter does that for a *diff* — `tokei`, `scc` and `cloc` all classify a
tree, not a range — and it is the last step between "lines added" and "code written".
Measured on this branch: **88.6 % statements, 4.4 % comments, 7.0 % blank.**

**On buying rather than building.** The language half of this is a solved problem and the
tools are good: `tokei` and `scc` are fast and accurate, `cloc` is in Debian. The category
half cannot be bought, because no tool knows that `engine/native/core/` is the core and
`analyses/` may be depended on by nothing — that mapping is this project's own, and it is
the half that makes the number mean something. Nothing off the shelf crosses the two, and
nothing does the statement/comment/blank split on a diff. Hence 200 lines here rather than
a dependency.

**The baseline matters more than the filter.** A diff against `main` is dominated by the
superseded tree Stage 1 deleted — 45,951 lines nobody wrote — so it answers "what changed"
and never "what have we built". The default baseline is `e4d964c`, the commit that laid the
layer tree and emptied the old one. Use `--vs-main` deliberately, when the 78 % cut is the
story.

Four repo-local aliases, set with `git config --local` and therefore not committed; recreate
them from here if the checkout is fresh:

| alias | what it shows |
|---|---|
| `git built` | the full language × category report |
| `git code` | the same, executable code only |
| `git src` | `--stat` over `*.cpp *.hpp *.py *.pyi` only |
| `git nodata` | `--stat` with `:(exclude)` on every data, prose and asset extension |
| `git live` | the report, redrawn whenever anything changes, **uncommitted work included** |

`--watch` is the one to leave open in a second pane. It polls a cheap fingerprint — `HEAD`,
`git status --porcelain`, and the numstat against `HEAD` — once a second and redraws only on
a change, so it costs nothing while you think. It diffs the bare baseline rather than
`base..HEAD`, which is what makes it live: edits count before they are committed. The header
carries the commit, how many files are dirty, and the delta since the last redraw.

`.gitattributes` does three things. It marks the machine-generated and measured files
`-diff`, so regenerating the equation view or landing a parameter document does not read as
thousands of hand-written additions — they stay versioned and reviewed, they stop inflating
a code diff. It sets `diff=cpp` and `diff=python`, so a hunk header names the enclosing
function rather than a stray brace. And `* text=auto` keeps line endings from ever becoming
a diff of their own.

**Measured 2026-09-02**, from the right baseline. Since the tree was emptied: **+5,950 of
executable code** — 4,654 C++, 1,247 Python, 49 build — of which **5,274 are statements**,
262 comments and 414 blank, against 10,113 lines that are prose, data or configuration. By category the code is 3,371 engine, 1,596 canon (the
native checks and the enforcers), 674 surface, 203 analysis scripts, 49 build; the core is
0, because Stage 2 is designed and not written.

Against `main` the same tree reads +15,836 / −69,022, and both numbers are true. One
answers what was built and the other what was removed, and no single figure answers both.

## A pinned source names the file as it was at that sha

A brief's `sources` are `branch@sha:path`, and the path must be the one the file had **at
that sha**, not the one it has now. Renaming a document and sweeping the new name through
the briefs breaks them silently: the preflight reports `source missing` and the sha is
perfectly valid.

This cost a repair on 2026-09-02, when four design documents were merged and three pinned
refs were rewritten to their successors. The fix is to leave the historical path and note
where the content lives now, which is what those three sources carry. Rebasing a brief to a
newer sha is the alternative, and it means re-verifying every marker.

## Preflight

`preflight.py` runs the brief's own claims against the worktree before the worker starts:
every source resolves at its sha and contains every marker; every gate has a mutation whose
file exists or is in `writes`; every precondition exits zero. It is the check that a brief
is executable rather than merely written, and it costs seconds.
