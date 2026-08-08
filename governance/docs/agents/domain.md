# Governance domain authority

This is one monorepo with one root context and directory-specific context.
Before Governance planning or authority changes, read root `CONTEXT.md`,
`governance/CONTEXT.md`, and `governance/GOVERNANCE.md`. Engine design decisions
remain under `engine/docs/adr/`; Governance does not create a parallel ADR tree.

Use these revision-13 terms:

- **Management**: the coordination role performed by the root task, not a
  repository, directory, or runtime owner.
- **Governance**: doctrine, topology, manifests, ledgers, readiness checks, and
  external-operation receipts under `governance/`.
- **Engine**: the sole `epcsaft` distribution under `engine/`, containing EOS,
  Equilibrium, and Regression.
- **Vertical slice**: one Engine-accountable task that carries a capability from
  exact Data input through production implementation to installed-wheel
  Validation evidence.
- **Evidence home**: a directory that owns durable inputs or results without
  becoming a separate implementation team.

Former repositories and names are provenance. They never select current task
residence, runtime behavior, or lifecycle state.
