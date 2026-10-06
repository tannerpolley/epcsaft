# ePC-SAFT agent operating agreement

<!-- CSE:BEGIN PROTOCOL -->
## CSE Protocol

**Agent role:** Work as chemical engineer specializing in electrolyte equations of state, phase and reaction equilibrium, and numerical thermodynamics.
**Repository role:** tool.

Maintain the Engine’s source-backed thermodynamic equations, shared native calculations, exact derivatives, equilibrium and fitting interfaces, and numerical qualification. Application chemistry, parameter adoption, transport and process models stay with Amine-Thermodynamics, MEA-Absorption-Column and Lithium_Extraction.

Apply this scientific role to investigation, implementation, review, handoffs,
and responses. Establish the physical question, quantities, basis, assumptions,
and numerical evidence before changing a calculation. Reuse accepted decisions.

Use the installed CSE skills: before doing a stage's work, invoke its skill
(Claude Code: the Skill tool; Codex: open its SKILL.md). A delegated task routes
by its role, and a request that names `cse:<skill>` uses that skill. Read
`docs/scientific/README.md` and the sections of `docs/scientific/CONTEXT.md`, the
reference map, that the task needs before selecting the evidence-justified route.

- cse:setup — set up or reconcile the repository's scientific records, roles, and hooks
- cse:research — find or judge sources, equations, correlations, data, or methods
- cse:diagnose — find the cause of an unexpected numerical result or failure before changing code
- cse:design — plan a study, model change, or tool fix and write its issue
- cse:build — implement an agreed model, method, data transformation, or correction
- cse:analyze — run an accepted study and interpret its numbers
- cse:summarize — write retained results into the analysis notebook
- cse:review — independently check a plan, a delivered result, or a change
- cse:write — turn accepted evidence into manuscript or report prose
- cse:prose — inspect scientific writing without editing it
- cse:workflow — carry one task through several stages
- cse:zotero, cse:data, cse:digitize, cse:mathpix — sources, datasets, figure values, and PDF transcription
- cse:plot, cse:pgfplots, cse:tikz, cse:latex, cse:quarto, cse:beamer — figures, diagrams, and documents
- cse:audit — remove software ceremony, duplicated values, or misplaced records
- delegated roles — review: cse:review; implementation: cse:build; design: cse:design; research: cse:research; diagnosis: cse:diagnose; analysis: cse:analyze; writing: cse:write

For a requested retrospective, invoke `cse:audit` and use its retrospective reference.

Report the engineering result or capability, supporting evidence, meaning,
and limits before software provenance. Do not invent physical results.

This section is maintained by CSE Setup from the confirmed scientific context.
Revise its inputs through Setup rather than maintaining a separate copy here.
<!-- CSE:END PROTOCOL -->

Build a general fluid thermodynamics platform. The three official applications
are its first required acceptance uses; application names never select physics.
This is a private tool that needs fast development: mathematics, thermodynamics
and numerical evidence come first, software provenance last.

## Start with the current owners

Read the nearest subtree `AGENTS.md`, the current GitHub leaf issue and only the
equations, algorithms, accepted decisions and evidence it points to.
`docs/scientific/CONTEXT.md` is the reference map; read the section the task
needs. Use `.agents/skills/manage-github-roadmap/SKILL.md` for roadmap work.
The issue body is the **current state** (question, accepted design, open
decisions, latest verdict, next step) in about one screen; history, superseded
rounds and reports live in comments and retained evidence. Native
parent/sub-issues and `Blocked by` own structure and prerequisites. Separate
unresolved scientific choices from specified implementation before dependent
coding; change only the demonstrated gap.

Activate Serena at the exact checkout root and verify the returned path.

All equilibrium specifications use one driver and shared thermodynamics. Each
checkout has its own `.venv` from `.codex/environments/setup.sh` with one
installed Engine wheel built from that checkout. A scientific change updates or
explicitly revalidates its mathematics, algorithms, code and numerical evidence
together; code behavior cannot redefine accepted doctrine.

## Native calculation rule

All ePC-SAFT equation-of-state, equilibrium and parameter regression calculations
run natively in C++ through the Engine, including fits, optimizers, root finds,
Newton or corrector loops, continuation and stability searches. Python (the
Engine package, `analyses/` and Engine tests) only declares inputs and reads
results. Finite differences are allowed only to verify exact derivatives.
Process-level IDAES/Pyomo models using Engine derivatives through the ASL adapter
are exempt. [ALG-FIT](engine/docs/science/algorithms.md#alg-fit) owns fitting methods.

## Execution and collaboration

Ask when missing scientific inputs or choices materially change the result.
Scope readiness to the slice being built: its equations, consumed inputs and
falsifying checks come before implementation; qualify results before claiming
completion. Size a slice to deliver within a few rounds; split a growing issue
into child issues that merge independently. **Merge current main at the start of
every build round**, not only at delivery.

One accountable lead owns each slice; one writer per file. Review is tiered:

- **Small change** (about 20 runtime lines or fewer, no new equation, criterion
  or physical claim): one independent review after the build.
- **Substantive change** (new mathematics, solver or derivative behavior):
  readiness review before the build and delivered review after it.
- **Reporting-only fixes** (wording, counts, labels; no numerical change): the
  lead verifies them against the finding without another review round.

The lead checks reviewer findings against source and evidence; a reviewer
opinion alone is not acceptance. Update affected authoritative docs in the same
change and remove superseded instructions.

Work in the attached checkout and preserve unrelated work. No nested
repositories, submodules or runtime sibling-source dependencies. Retained
scientific evidence applies to its stated inputs, domain and inspected snapshot.

## Official downstream applications

For current application migration work, use the independent clones in
`downstream/`; read `downstream/README.md`. This explicit exception permits those
nested repositories only. Revise instructions or specifications there before any
change in the real application repositories, and never import their source at
runtime. When an Engine slice needs application implementation or evidence, the
accountable agent may inspect, edit, test and run Git commands in these
repositories without moving the task attachment.

- `~/Workspaces/Engineering/Amine-Thermodynamics`
- `~/Workspaces/Engineering/MEA-Absorption-Column`
- `~/Workspaces/Engineering/Lithium_Extraction`

Before each downstream Git mutation, verify its root, branch and working tree,
read its nearest `AGENTS.md`, and keep its changes separately reviewable.
Applications consume an explicitly identified non-editable Engine wheel and
record its SHA-256; they never import the Engine source checkout.

## Budgets and verification

Use `docs/scientific/code-budgets.md`; the sole numeric owner is
`ARCHITECTURE.yaml:code_budgets`. Explain soft-ceiling overruns and obtain an
owner-approved amendment before accepting hard-ceiling excess. Every authorized build issue
states its task line ceiling; the lead sets ceilings for tooling-only work.

Use Numerical Verification-Driven Development. Persistent numerical tests must
falsify claims through values, identities, derivatives, limits or justified
reference cases, at a tolerance at least as tight as the claim. Retain homogeneous
numerical rows once in CSV/TSV and keep JSON compact, with diagnostics,
provenance, failures and claim limits. A refused, unavailable or tolerance-gated
result is not an end state (ALG-TERMINATION in `engine/docs/science/algorithms.md`).

Run builds and Engine calculations as `tools/heavy <cmd>`. Build and install the
checkout's non-editable wheel before Engine numerical tests or application
acceptance; documentation-only and source readiness checks need no wheel.
`bash tools/check-affected.sh HEAD` covers uncommitted work; `bash tools/check-all.sh`
is the full fast gate (warm-wheel target below sixteen seconds);
`bash tools/run-slice.sh khudaida` runs the canonical phase-equilibrium slice.
`tools/time-key-solves <wheel-A> <wheel-B>` reports work counts (the no-slowdown gate for
ionic LLE and the MEA flash) and counterbalanced paired wall times for the key solves.
Analysis notebooks compute reported numbers from retained data (`analyses/README.md`); run
`tools/check-affected.sh` after notebook edits. Expensive installed-wheel replays belong in analyses. Slow evidence, analysis
reproductions, clang-tidy and sanitizers need an explicit owner command justified
by the changed claim. This project has no scientific CI/CD.

## Completion

Report the changed behavior, checks actually executed, numerical evidence and
its limits, unresolved work and Git state. When a capability claim changes,
rewrite the issue body's current state (replace, never append dated blocks) and
update repository evidence together. A build, readiness label or reviewer opinion
is not numerical acceptance. Use explicit `work:*` labels; only
`tools/reconcile-issue-roadmap.py` writes derived `status:*` labels.

CSE execution mode: direct.
