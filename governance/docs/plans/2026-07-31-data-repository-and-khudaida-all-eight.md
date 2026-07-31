# Canonical Data repository and Khudaida Figure 2 plan

Status: approved for execution

Parent: https://github.com/ePC-SAFT/ePC-SAFT-management/issues/9

## Outcome

Create `ePC-SAFT/ePC-SAFT-data` as the canonical store for reusable scientific
packets and use an immutable Khudaida Figure 2 packet to reproduce all eight
reported points through the installed EOS and Equilibrium public artifacts.

## Ownership

- Data stores reviewed source-derived bytes, metadata, provenance, units,
  domains, hashes, and license or use bases.
- EOS ships only an admitted, self-contained runtime-catalog subset copied or
  generated from a pinned Data release.
- Validation stores campaign definitions, generated results, plots, and
  receipts that pin exact Data commits and packet fingerprints.
- Lab stores exploratory, rejected, broad, or incomplete source material.
- Management stores policy, registries, access rules, operation receipts, and
  package artifact objects.

These roles are independent: `storage_owner = data`, `runtime_owner = eos`
only for an admitted catalog copy, and `evidence_owner = validation` when a
campaign exists.

## Execution sequence

1. Merge doctrine revision 6 and its machine-readable topology.
2. Create the public Data repository and record its immutable identity.
3. Publish the packet contract and migrate the reviewed Khudaida packet.
4. Pin the packet commit and fingerprint in each consuming campaign.
5. Make EOS packing-constraint callbacks evaluable at optimizer trial states
   while retaining explicit accepted-state feasibility gates.
6. Run all eight Figure 2 rows against immutable installed artifacts, diagnose
   generic remaining failures, and fix the owning package without row-specific
   branches or relaxed scientific tolerances.
7. Retain independent review and final Validation evidence for the exact
   accepted subjects.

## Acceptance boundary

Passing means every one of the eight declared observations completes the
intended equilibrium route and satisfies the source-backed numerical tolerance.
Solver termination alone is insufficient. The campaign must retain packet,
artifact, commit, tolerance, and failure-boundary fingerprints.

The campaign does not claim global phase completeness, electrolyte capability
beyond the declared domain, or runtime dependence on the Data repository.
