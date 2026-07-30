# Active SDK Channel Cleanup

Status: executed locally on 2026-07-30; publish through the Management
change review.

## Decision

The project has one live wheel per installable package. Version equality is
not enough to identify an SDK, so every active subject is bound to its source
repository commit and full wheel SHA-256 in
`manifests/artifact-store-policy.json`.

| Distribution | Owner | Source commit | Wheel SHA-256 |
| --- | --- | --- | --- |
| `epcsaft` | `ePC-SAFT/ePC-SAFT-eos` | `24ab1bdfb3a1558eb87815e9c3b5d97fbe6a025d` | `66b7ea8fb29e0a268b555cbdf401c3502517c088669a4157e8f64ab985b59ce9` |
| `epcsaft-equilibrium` | `ePC-SAFT/ePC-SAFT-equilibrium` | `13b171b0a8fbff5b2ee6d97d3c0ad21c7868d30b` | `1482778b113df890931de692262c291111b979bb4762f85a7e352dbd2d5b0266` |
| `epcsaft-regression` | `ePC-SAFT/ePC-SAFT-regression` | `fc42629253125ad327150473872994d7944545ad` | `779c9346e47e6fba109d7147f6177c42bb84b5c105b705fd1311b75f40e0d409` |

The live aliases are under `artifacts/active/<distribution>/`. Each has one
matching immutable object under `artifacts/objects/sha256/<full-sha256>/`.
Equilibrium and Regression Codex setup scripts pin the EOS object above and
reject any other installed EOS wheel.

## Historical disposition

The former 107 wheel aliases and 106 derived object copies were removed from
the live channel after their bytes and pre-cleanup indexes were captured in
`artifacts/archive/legacy-sdk-wheels-20260730.tar.zst`. Its hash and counts
are recorded in `manifests/artifact-archive.json`. The archive is evidence for
replay/provenance only; no setup script or active index may install from it.

## Ongoing rule

Replacing an active wheel is a coordinated three-package operation:

1. build in disposable staging from a named repository commit;
2. install and smoke-test the exact wheel in an isolated environment;
3. update the one active alias, its full-SHA object, and the Management policy;
4. rebuild dependent wheels against the new EOS object when their native SDK
   contract requires it; and
5. run `artifact_store.py check` and strict Management preflight before the
   change is published.

No second alias, mutable repository `dist/` wheel, or uncompressed historical
wheel may be added. If replay evidence is required, append a new compressed
archive and record its manifest rather than enlarging the active SDK set.
