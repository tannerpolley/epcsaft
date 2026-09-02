# Greenfield

**Overview and live status:** https://claude.ai/code/artifact/1bb42e8c-f3ef-4c72-99e1-335f7614bf58

That page is the reading surface. The documents below are what agents execute against;
where they disagree, they are canonical and the page is stale.

| document | what it owns |
|---|---|
| `TRANSFER.md` | the stages: order, scope, budgets, gates |
| `ORCHESTRATION.md` | who does what, the verified `codex exec` invocation, stop conditions, environment traps |
| `TESTING.md` | what may be tested, and what may not |
| `CARRY-OVER.md` | what is copied, what is written fresh |
| `STAGE-LOG.md` | estimate vs actual per stage, and every decision taken while writing it |
| `briefs/stage-N.brief.{md,json}` | one stage: sources with markers, budgets with ceilings, gates with mutations |
| `../../tests/oracle/oracle.json` | the acceptance data: measured rows only |
| `../../../tools/check-engine-architecture.py` | the seven rules |

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
