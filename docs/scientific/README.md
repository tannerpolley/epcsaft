# Scientific records

This is the CSE owner map for the ePC-SAFT monorepo. Start with
[CONTEXT.md](CONTEXT.md) for mathematical authority, methods, implementation
and evidence owners. [ADRs](adr/README.md) retain the sole accepted decision
sequence; [writing](writing/STYLE.md) owns manuscript and notebook style.

Canonical mathematics and algorithms use Markdown with LaTeX display equations
and plain pseudocode under `engine/docs/`. The equation definitions are edited
once; their YAML index is generated. Optional HTML/PDF views must render the
same source. Manuscripts may use LaTeX and analysis notebooks may use Quarto;
neither becomes a second Engine specification. This is the repository's accepted CSE Markdown layout.

The [GitHub issue roadmap](https://github.com/tannerpolley/ePC-SAFT-greenfield/milestones)
is the sole current gap and completion plan. Issue bodies own design, scope,
acceptance, reuse and handoff; native relationships own hierarchy and
prerequisites. The local [roadmap authority note](implementation-plan.md) contains
no editable status. [Dated reviews](reviews/transfer-audit-2026-09-06.md) retain
numerical evidence and inspected limitations; they do not issue current
instructions. Keep a record only while it supports a scientific result,
derivation, accepted decision or diagnosis. Update current scientific owners in
place and remove duplicate progress summaries; Git preserves edit history.

[Code budgets](code-budgets.md) define the floor target, soft ceiling and hard
ceiling; their single numeric owner is `ARCHITECTURE.yaml`.

The [modular-model proposal](proposals/modular-models/README.md) retains a
source-linked design and separate five-step roadmap for selecting complete
thermodynamic formulations with compatible parameters. It is deferred until
the current greenfield roadmap is complete and adds no current milestones or
accepted physics.

## Sources

The Zotero collection **ePC-SAFT-Favorites** (`YR8SNI5D`) is the Engine's reading
list, mirrored by the Zotero Companion project `epcsaft-equilibrium`. Start with
`cse-zotero wiki-context epcsaft-equilibrium --query "<terms>" --json` (add
`--global` when the project view is empty) and read exact passages with
`cse-zotero wiki-evidence <KEY> --query "<terms>" --json`. Equation text and
stable IDs stay in [equations.md](../../engine/docs/equations.md); Companion
`implementation_reference` records point to them. Citation keys used by the
Engine docs are listed in `engine/docs/latex/cited_keys.txt`.

## Live companion

The root Quarto project presents the current map, a link to the GitHub roadmap,
ordinary check reports, and retained thermodynamic comparisons as one read-only
local website. Start it with
`quarto preview . --host 127.0.0.1 --port 8766 --no-browser --render all --no-execute`.
Quarto owns rendering, navigation, watching, serving, and reload. Its
presentation reader only reads current owners and completed records; it never
builds, collects, tests, calculates, or mirrors issue status. GitHub issues own
roadmap state; repository evidence owns scientific qualification. The
[test-system design](test-system-design.md) defines the one-way result and freshness
and discovery behavior. Its [implementation acceptance](reviews/test-system-acceptance-2026-09-09.md)
records the retired system's bounded checks; it is historical evidence, not a
current execution contract.

## Maintaining agreement

A scientific change updates or explicitly revalidates its equation IDs,
algorithm IDs and domain assumptions in the same change as implementation and
numerical evidence. Record that disposition in the change description; do not
copy equations into another ledger. A code defect is corrected against the
accepted mathematics. A proposed change to mathematics needs its own source or
derivation and review; observed code behavior cannot silently redefine it.

Run `python3 engine/scripts/docs/sync_equation_registry.py` after editing
[definitions](../../engine/docs/equations.md), then `bash tools/check-affected.sh HEAD`
for uncommitted work. The normal checks reject malformed or duplicate equation
IDs, duplicate algorithm IDs, unresolved equation/algorithm references, and a
stale generated index. They cannot prove
scientific agreement: result-level numerical checks and review must falsify
wrong coefficients, derivative chains, bases and stopping claims. Update the
owning issue and retained evidence when a capability changes.
