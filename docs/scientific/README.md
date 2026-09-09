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

The [implementation plan](implementation-plan.md) is the single current gap
and completion plan. [Dated reviews](reviews/transfer-audit-2026-09-06.md)
retain numerical evidence and inspected limitations; they do not issue current
instructions. Keep a record only while it supports a scientific result,
derivation, accepted decision or diagnosis. Update current owners in place and
remove superseded instructions and duplicate progress summaries; Git preserves
edit history. GitHub issues track bounded delivery and link to the plan's
scientific qualification rows rather than repeating their status.

[Code budgets](code-budgets.md) define the floor target, soft ceiling and hard
ceiling; their single numeric owner is `engine/ARCHITECTURE.yaml`.

## Live companion

The Quarto companion in `docs/companion/` presents the current implementation
plan, thermodynamic comparisons, and implemented/planned verification as one
local website with Map, Plan and Tests pages. Start the read-only viewer with
`python3 tools/serve-roadmap.py --port 8766`. Rendering pages reads current owners
and completed records without running scientific calculations or tests; results
come from the normal check commands. The implementation plan owns the
live-view behavior and limits. The
[test-system design](test-system-design.md) defines automatic recording, freshness
and discovery behavior. Its [implementation acceptance](reviews/test-system-acceptance-2026-09-09.md)
records the bounded checks. The implementation plan reopens startup persistence
and integrity-check visibility before accepting this prerequisite to Core work.

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
transfer plan only when identified evidence supports the changed capability.
