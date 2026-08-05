# Data Repository Agent Contract

Repository Profile: scientific-computing

Before work, read doctrine revision 10 and the machine-readable topology from
`../ePC-SAFT-governance/GOVERNANCE.md` and
`../ePC-SAFT-governance/manifests/management.json`. Their published home is
`https://github.com/ePC-SAFT/ePC-SAFT-governance`.

This repository is the canonical store for reusable reviewed scientific
packets: parameter bundles, source-derived datasets, normalized observations,
provenance, units, domains, hashes, and license or use bases.

It must not own runtime algorithms, generated Validation results, package
artifacts, full papers, or exploratory and rejected material. Installed
packages must not import from or fetch this repository at runtime. Engine may ship
only an admitted self-contained copy generated or copied from a pinned Data
commit; Validation and Engine Regression pin exact packet commits and fingerprints.

Packets live only at `packets/<packet-id>/<version>/packet.toml` and follow
`docs/packet-contract.md`. Do not add empty speculative packet directories.

Use Numerical Verification-Driven Development from Governance section 7.
Data owns no pytest or unittest suite. Packet schema, hashes, provenance,
units, domains, and declared tolerance metadata fail closed through
`tools/validate_packets.py`. The named `Validate packets` environment action is
the sole repository readiness gate, not a scientific test suite. It may check
packet-internal normalization, unit and basis transformations, declared
arithmetic, fingerprints, and consistency. Numerical model agreement belongs
to Engine or Validation. This private repository has no CI/CD. Run that action
and the repository cleanup audit before handoff.

## Agent skills

<!-- project-truss:setup:start -->
### Project Truss

Use `project-truss:setup` once, then `project-truss:start` to begin or resume Matt-first engineering work. All six Project Truss skills are callable; Start normally invokes the relevant Matt disciplines and routes each stage until verified closeout.

### Issue tracker

Project Truss uses GitHub-native root and leaf issues. See `docs/agents/issue-tracker.md`.

### Domain docs

Read the configured domain context and ADRs. See `docs/agents/domain.md`.
<!-- project-truss:setup:end -->
