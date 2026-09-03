# Greenfield

**Overview and live status:** https://claude.ai/code/artifact/1bb42e8c-f3ef-4c72-99e1-335f7614bf58

That page is the reading surface. The documents below are what agents execute against;
where they disagree, they are canonical and the page is stale.

| document | what it owns |
|---|---|
| `TRANSFER.md` | **the map of the work**: the eight stages in category order, each one's writes, design document and gate; the oracle; the crosswalk from the old numbering |
| `ORCHESTRATION.md` | who does what, the verified `codex exec` invocation, stop conditions, environment traps |
| `CARRY-OVER.md` | what is copied, what is written fresh, what has no successor |
| `STAGE-LOG.md` | estimate vs actual per stage, and every decision taken while writing it |
| `HANDOFF.md` | what a fresh coordinator session reads first: live process, acceptance steps, standing rules, open items |
| `target-layout.md` | the finished tree, file by file, with line budgets |

Design documents. **One per stage, named for the stage**, and each is what that stage's
brief is written against:

| document | stage | what it owns |
|---|---|---|
| `CANON.md` | 1 Canon | the tree, the archive, the seven rules, and what may be tested |
| `CORE.md` | 2 Core | where the core lives, numerical policy absorbed, the declared-choice and convergence contracts, debug mode, and the `diagnostics` block every later brief must carry |
| `SURFACE.md` | 3 Surface | what a caller declares and reads: the parameter builder, the bindings, the Python mapping layer, and where the tree departs from its brief |
| `EOS.md` | 4 EOS | the layer audit, the like-ion dispersion convention, the density root, and the measured performance floor |
| `EQUILIBRIUM.md` | 5 Equilibrium | how the routes prove they are general: the count rule, the lattice, what each check may assert |
| `REGRESSION.md` | 6 Regression | *(to write)* |
| `PHASE-STABILITY.md` | 7 Phase stability | the HELD2 kernel read against Pereira 2012, finding by finding |
| `FINAL-CHECK.md` | 8 Final check | *(to write)* |

Per stage, and the data:

| | what it owns |
|---|---|
| `briefs/stage-NN.brief.{md,json}` | one stage: sources with markers, budgets with ceilings, gates with mutations. `NN` is the stage number in `TRANSFER.md` |
| `../../engine/tests/oracle/oracle.json` | the acceptance data: measured rows only |
| `../../tools/check-engine-architecture.py` | the seven rules |

## Roles

**Coordinator (Opus 5).** Owns every process, every git write, and acceptance. Runs the
preflight, drives the worker, and re-runs every gate and mutation independently. Never
writes stage code.

**Planner (Fable 5.1).** Reads the tree, decides the design, writes the brief. Solves the
novel problems. Never executes, never supervises the executor — a subagent reports only
through its final response, so a worker under the planner would be unobservable.

**Worker.** One stage per invocation, in an isolated worktree, against a brief that leaves
nothing to judgement. Verifies each source marker before writing. Stops rather than
adapting. Never plans, never accepts its own work, never repairs its environment.

A worker's summary is never acceptance. The coordinator rebuilds from a clean directory,
runs the gates, and applies every mutation to confirm each gate goes red.

## Carry rule

Carry an equation, write a decision. An equation is a fact, checkable against a closed
form; copying it is safe. A structure is a decision, and copying it imports every
assumption behind it. Every defect the audit found lives in carried structure.
