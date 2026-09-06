# Scientific records

This is the CSE owner map for the ePC-SAFT monorepo. Start with
[CONTEXT.md](CONTEXT.md) for mathematical authority, methods, implementation
and evidence owners. [ADRs](adr/README.md) retain the sole accepted decision
sequence; [writing](writing/STYLE.md) owns manuscript and notebook style.

Canonical mathematics and algorithms use Markdown with LaTeX display equations
and plain pseudocode under `engine/docs/`. The equation definitions are edited
once; their YAML index is generated. Optional HTML/PDF views must render the
same source. Manuscripts may use LaTeX and analysis notebooks may use Quarto;
neither becomes a second Engine specification. This preserves the repository's
accepted CSE Markdown layout (issue 55).

[TRANSFER.md](../greenfield/TRANSFER.md) is the single current implementation
gap plan. Dated reviews, prior briefs and stage documents retain provenance;
they neither override permanent doctrine nor authorize execution. The
[numerical-policy audit](numerical-policy-audit.md) is historical evidence.

[Code budgets](code-budgets.md) define the floor target, soft ceiling and hard
ceiling; their single numeric owner is `engine/ARCHITECTURE.yaml`.

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
