# Scientific records

This is the CSE owner map for the ePC-SAFT monorepo. Start with
[CONTEXT.md](CONTEXT.md) for mathematical authority, methods, implementation
and evidence owners. [ADRs](decisions/README.md) retain the sole accepted decision
sequence; [writing](writing/STYLE.md) owns manuscript and notebook style.
Research documents and dated evidence keep their original source locators and
historical software names; they inform scientific choices but do not replace the
current specifications.

Canonical mathematics and algorithms use Markdown with LaTeX display equations
and plain pseudocode under `engine/docs/`. The equation definitions are edited
once; their YAML index is generated. Optional HTML/PDF views must render the
same source. Manuscripts may use LaTeX and analysis notebooks may use Quarto;
neither becomes a second Engine specification. This is the repository's accepted CSE Markdown layout.

[ADR 0010](decisions/0010-github-issue-roadmap-authority.md) gives GitHub issues the plan, scope and acceptance.
[Dated reviews](reviews/transfer-audit-2026-09-06.md) retain
numerical evidence and inspected limitations; they do not issue current
instructions. Keep a record only while it supports a scientific result,
derivation, accepted decision or diagnosis. Update current scientific owners in
place and remove duplicate progress summaries; Git preserves edit history.

[Code budgets](code-budgets.md) define how runtime lines are counted and how
each ceiling is raised or lowered; their single numeric owner is `ARCHITECTURE.yaml`.

[ADR 0011](decisions/0011-optional-dipole-dipole-term.md) owns the optional
Gross–Vrabec dipole–dipole admission; remaining multipoles stay retired under
[ADR 0006](decisions/0006-retire-residual-polar.md).

## Sources

The Zotero collection **ePC-SAFT-Favorites** (`YR8SNI5D`) is the Engine's reading
list, mirrored by the Zotero Companion project `epcsaft-equilibrium`. Start with
`cse-zotero wiki-context epcsaft-equilibrium --query "<terms>" --json` (add
`--global` when the project view is empty) and read exact passages with
`cse-zotero wiki-evidence <KEY> --query "<terms>" --json`. Equation text and
stable IDs stay in [equations.md](../../engine/docs/equations.md); Companion
`implementation_reference` records point to them. Citation keys used by the
Engine docs are listed in `engine/docs/latex/cited_keys.txt`.

## Status report

`tools/report` writes `build/report.html` from retained test, key-solve and analysis
records; the [test-system design](test-system-design.md) defines its one-way flow and
freshness flags.

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
