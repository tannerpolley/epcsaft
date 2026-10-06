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

## Every task

Read the nearest subtree `AGENTS.md`, the current GitHub leaf issue and only the
equations, algorithms, accepted decisions and evidence it points to;
`docs/scientific/CONTEXT.md` is the reference map, read by section. Separate
unresolved scientific choices from specified implementation before dependent
coding; change only the demonstrated gap.

Activate Serena at the exact checkout root and verify the returned path. Each
checkout has its own `.venv` from `.codex/environments/setup.sh` with one
installed non-editable Engine wheel built from that checkout; run builds and
Engine calculations as `tools/heavy <cmd>`. All equilibrium specifications use
one driver and shared thermodynamics.

## Native calculation rule

All ePC-SAFT equation-of-state, equilibrium and parameter regression calculations
run natively in C++ through the Engine, including fits, optimizers, root finds,
Newton or corrector loops, continuation and stability searches. Python (the
Engine package, `analyses/` and Engine tests) only declares inputs and reads
results. Finite differences are allowed only to verify exact derivatives.
Process-level IDAES/Pyomo models using Engine derivatives through the ASL adapter
are exempt. [ALG-FIT](engine/docs/science/algorithms.md#alg-fit) owns fitting methods.

## Leading a slice

Scope readiness to the slice being built: its equations, consumed inputs and
falsifying checks come before implementation. Size a slice to deliver within a
few rounds; split a growing issue into child issues that merge independently.
**Merge current main at the start of every build round.**

One accountable lead owns each slice; one writer per file. Review is tiered:

- **Small change** (about 20 runtime lines or fewer, no new equation, criterion
  or physical claim): one independent review after the build.
- **Substantive change** (new mathematics, solver or derivative behavior):
  readiness review before the build and delivered review after it.
- **Reporting-only fixes** (wording, counts, labels; no numerical change): the
  lead verifies them against the finding without another review round.

The lead checks reviewer findings against source and evidence; a build,
readiness label or reviewer opinion is not numerical acceptance. When a
capability claim changes, rewrite the issue body's current state and update
repository evidence together.

## Where things live

- **Roadmap work:** `.agents/skills/manage-github-roadmap/SKILL.md`.
- **Review standards:** `CODING_STANDARDS.md`.
- **Checks:** `README.md` owns the commands. `tools/check-all.sh` is the full
  fast gate (warm target below sixteen seconds), and `tools/time-key-solves`
  work counts are the no-slowdown gate for ionic LLE and the MEA flash. Slow
  evidence, analysis reproductions, clang-tidy and sanitizers need an explicit
  owner command; there is no scientific CI/CD.
- **Code budgets:** `docs/scientific/code-budgets.md`. Every authorized build
  issue states its task line ceiling; the lead sets ceilings for tooling-only work.
- **Application work:** the clones in `downstream/` under `downstream/README.md`.
  An Engine slice that needs application implementation or evidence may inspect,
  edit, test and run Git in the application repositories without moving the
  task attachment.
- **Agent skills:** `docs/agents/` (issue tracker, triage labels, domain docs).

CSE execution mode: direct.
