# Unified ePC-SAFT Engine context

`governance_doctrine_revision: 10`

Canonical doctrine and cutover state live in
`../ePC-SAFT-governance/GOVERNANCE.md` and
`../ePC-SAFT-governance/manifests/runtime-consolidation.json`.

This source tree is the durable Unified Engine repository. It builds the sole
Python 3.13 runtime distribution, `epcsaft`, containing EOS, equilibrium, and
regression. Release `v0.2.0.dev0` at source commit
`c712fe9730faef0732341ccdabc4d799fef2d953` is the accepted package subject.

## Runtime structure

- `epcsaft` owns typed parameter records, bundle I/O, resolved models, the
  residual Helmholtz equation, density closure, state/property results,
  reference states, and the model-bound native transport.
- `epcsaft.equilibrium` owns saturation, phase equilibrium, HELD/HELD2, and
  chemical/species equilibrium results and diagnostics.
- `epcsaft.regression` owns typed observations, objectives, solver controls,
  fits, diagnostics, and immutable transfer of accepted fitted values back to
  `Parameters`.

All native modules link in the same root CMake build and ship in the same wheel.
Internal native boundaries may stay modular, but they are not independent SDK
installations or release identities.

## Scientific transfer

Reusable parameters and scientific datasets remain owned by Data. Callers
materialize an exact packet and pass its explicit parameter directory to
`Parameters.from_bundle`; Engine never searches for Data. A materialization
receipt records the Data commit, packet identity, path, version, full packet
fingerprint, destination, and file hashes.

Regression produces a typed immutable `RegressionResult`. A usable result can
create a new parameter bundle only through `RegressionResult.to_parameters`,
which verifies the source parameter fingerprint, exact parameter identities,
source provenance, validity domain, and destination bundle identity. The
resulting `Parameters` object is consumed directly by EOS and Equilibrium with
no serialization, package discovery, or filesystem handoff.

## Evidence and authority

Retained unit tests protect numerical algorithms and scientific invariants.
Cross-module integration tests protect meaningful scientific workflows.
Validation owns installed-wheel campaigns and promotion evidence. Data owns
packet validation. Lab and applications consume immutable accepted artifacts.

Historical plans, reviews, and receipts retain their original names and hashes.
They do not define current imports, setup paths, package identities, or runtime
authority. No compatibility layer is maintained for the former split-package
topology.
