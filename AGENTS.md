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
are its first acceptance uses; application names never select physics.
Mathematics, thermodynamics and numerical evidence come first, software
provenance last.

## Every task

Read the nearest subtree `AGENTS.md`, the current GitHub leaf issue and only the
equations, decisions and evidence it points to; change only the demonstrated
gap. For domain terms and decisions, read the list in `docs/agents/domain.md`.
The leaf issue, not `engine/docs/`, says what the installed Engine
implements. Activate Serena at the exact checkout root. `README.md` "Checks"
owns the per-checkout `.venv` with its editable Engine install.

## Native calculation rule

All ePC-SAFT equation-of-state, equilibrium and parameter regression calculations
run natively in C++ through the Engine, including fits, optimizers, root finds,
Newton or corrector loops, continuation and stability searches. Python (the
Engine package, `analyses/` and Engine tests) only declares inputs and reads
results. Finite differences are allowed only to verify exact derivatives.
Process-level IDAES/Pyomo models using Engine derivatives through the ASL adapter
are exempt. [ALG-FIT](engine/docs/science/algorithms.md#alg-fit) owns fitting methods.

## Pointers

- **Leading a slice** (readiness, review tiers, merging main, acceptance): `docs/agents/workflow.md`.
- **Roadmap work:** `.agents/skills/manage-github-roadmap/SKILL.md`.
- **Review standards:** `CODING_STANDARDS.md`.
- **Gates and the no-slowdown gate:** `docs/agents/workflow.md` "Gates"; commands and environment: `README.md` "Checks".
- **Code budgets:** `docs/scientific/code-budgets.md`.
- **Application work:** `downstream/README.md`.
- **Skill configuration** (domain docs, issue tracker, triage labels): `docs/agents/`.

CSE execution mode: direct.
