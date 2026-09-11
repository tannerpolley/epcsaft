# ePC-SAFT agent operating agreement

Build a general fluid thermodynamics platform. The three official applications
are its first required acceptance uses; application names never select physics.

## Start with the current owners

Use the installed CSE skills. Read `docs/scientific/README.md`,
`docs/scientific/CONTEXT.md` and the nearest subtree `AGENTS.md` before selecting
the smallest applicable CSE route. Then read only the relevant equations,
algorithms, accepted decisions and `docs/scientific/implementation-plan.md`
entries. Before each milestone implementation, apply that plan's reuse map and
required adoption decision: preserve justified current/old behavior and change
only the demonstrated gap. Its deferred scope is not a completion prerequisite. The project map owns architectural detail and directory responsibilities.

Before Serena inspection, activate the exact checkout root and verify the returned
path. Use Serena for C++/Python symbol and caller navigation; normal environment
setup maintains its root compilation database. Symbol queries synchronize source
file changes. Read the current owner files for guidance; activation messages and
memories are navigation pointers, not cached authority. Ordinary source and doc
edits require no manual refresh. After changing Serena's own configuration, verify
it with a fresh server session; same-project activation can reuse cached settings.

All equilibrium specifications use one driver and shared thermodynamics.
Preserve one installed Engine wheel and explicit physical inputs. A scientific
change updates or explicitly revalidates its mathematics, algorithms, code and
numerical evidence together. Code behavior cannot redefine accepted doctrine;
unresolved scientific choices remain explicit in the implementation plan.

## Execution and collaboration

Complete authorized work and relevant verification without repeated permission.
Resolve routine implementation choices from current evidence. Ask when missing
scientific inputs or choices materially change the result; continue independent
work while waiting. Prepare a concrete, reviewable result before requesting any
additional authority needed for external or destructive actions.

Scope readiness to the scientific slice being built: establish its equations,
consumed inputs and falsifying checks before implementation, then qualify its
results before claiming completion. Resolve unrelated preparation with its own
slice. Treat downstream clones as independent current applications; investigate
historical inconsistencies only when they affect the current calculation.
Update affected authoritative docs and their live views in the same change;
remove superseded instructions instead of appending competing policy. Keep
progress in the implementation plan and delivery in linked issues. Retain
records for their scientific use, never merely to narrate task history.

One accountable lead owns each scientific slice across its required files.
Keep one writer per file. Delegate bounded independent inspections when parallel
work improves evidence or speed; substantive mathematics, solver and derivative
changes receive independent review. The lead checks reviewer findings against
source and evidence; a reviewer opinion alone is not acceptance.

Use the user-selected task model as lead. Keep model identifiers, reasoning
settings and fallback choices in task configuration. Never silently substitute
another model for an explicitly requested Astra task or delegate. Do not impose
a permanent planner/coordinator/worker hierarchy. Keep progress and completion
reports concise; give detailed scientific explanations when requested.

Work in the attached checkout and preserve unrelated work. Do not create nested
repositories, submodules, runtime sibling-source dependencies or separate Engine
installations. Superseded setup and worker instructions live only in Git history;
consult them only for a specific provenance or recovery question. Retained
scientific evidence applies to its stated inputs, domain and inspected snapshot.

## Official downstream applications

For current application migration work, use the independent clones in
`downstream/`; read `downstream/README.md`. This explicit exception permits those
nested repositories only. Revise instructions/specifications there before any
change in the real application repositories. Do not import their source at runtime.

- `~/Workspaces/Engineering/MEA-Thermodynamics`
- `~/Workspaces/Engineering/MEA-Absorption-Column`
- `~/Workspaces/Engineering/Lithium_Extraction`

When a requested Engine slice needs application implementation or evidence, the
accountable agent may inspect, edit, test and run Git commands in these named
repositories without moving the task attachment, subject to filesystem access.
Before each downstream Git mutation, verify its root, branch and working tree,
read its nearest `AGENTS.md`, and keep its changes separately reviewable.
Applications consume an explicitly identified non-editable Engine wheel and
record its SHA-256; they never import the Engine source checkout.

## Budgets and verification

Use `docs/scientific/code-budgets.md`; the sole numeric owner is
`ARCHITECTURE.yaml:code_budgets`. Work toward the floor target, explain
soft-ceiling overruns, and obtain an owner-approved amendment before accepting
hard-ceiling excess. Preserve readable code and required science when simplifying.

Use Numerical Verification-Driven Development. Persistent numerical tests must
falsify claims through values, identities, derivatives, limits or justified
reference cases. Do not add API-contract, compatibility, snapshot, prose, mock
or coverage-only tests. Retain homogeneous numerical rows once in CSV/TSV and keep JSON compact,
with diagnostics, provenance, failures and claim limits. Expensive installed-wheel
replays belong in analyses. Readiness checks do not establish numerical agreement.

The sole environment is `.codex/environments/environment.toml`. Build and install
the current non-editable wheel before Engine numerical tests or application
acceptance. Documentation-only and source readiness checks can run without
building the wheel.
Use `bash tools/check-affected.sh HEAD` for uncommitted work;
`bash tools/check-all.sh` is the full fast installed-wheel gate, with a warm-wheel
target below ten seconds. Slow evidence, analysis reproductions, full native
clang-tidy and sanitizers require an explicit owner command justified by the
changed claim. `bash tools/run-slice.sh khudaida` runs the canonical phase-equilibrium
slice. Broaden or repeat verification only for new changes, failures or unresolved
risks. This project has no CI/CD.

## Completion

Report the changed behavior, checks actually executed, numerical evidence and
its limits, unresolved work and Git state. Update the implementation plan only
when code and evidence support the new capability claim. A successful build,
readiness check or reviewer opinion is not numerical acceptance.

Run `bash "$HOME/.codex/hooks/codex-cleanup.sh" --repo-root .` after file changes
or task-owned processes/artifacts. Agents run it; never ask the user to do so.
Remove only task-created ignored artifacts or exact user-authorized candidates.

For roadmap work, update the implementation plan's authoritative milestone row at
start, completion and handoff: status, remaining/actual Astra elapsed hours and
a dated evidence record. Reopen affected and dependent qualifications when a
claim changes. Run `python3 tools/check-roadmap.py --write` after row edits; the
normal readiness gate checks the derived summary. No readiness check alone
qualifies a scientific row.
