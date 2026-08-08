# Data domain authority

Before Data work, read root `CONTEXT.md`, `data/CONTEXT.md`, and
`data/docs/packet-contract.md`. Engine architecture decisions remain under
`engine/docs/adr/`; Data owns scientific packets, not runtime design or a
parallel ADR tree.

Keep durable packet identity, provenance, units, domains, hashes, and use basis
in Data. Lifecycle state remains in GitHub. An Engine-accountable monorepo task
may revise the exact packet and consumers required by one vertical slice.
