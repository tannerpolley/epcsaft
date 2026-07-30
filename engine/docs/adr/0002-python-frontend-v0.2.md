# ADR 0002: Python frontend v0.2

- Status: accepted for the `0.2.0.dev0` frontend scope
- Date: 2026-07-28
- Supersedes: the pre-0.2 Python walkthrough and naming surface
- Does not supersede: ADR 0001 or any promotion receipt

## Context

The EOS implementation now has one typed parameter lifecycle and one
state/property lifecycle. The active documentation must expose that lifecycle
directly and must not preserve a second facade through aliases or duplicate
result representations. The frontend remains a Python presentation layer over
the EOS's native Helmholtz and density owners.

## Decision

The sole public Python frontend is:

```text
Parameters -> Mixture -> State
                         |
                         -> fugacity
                         -> activity for reference-created states
```

### Parameter entry points

- `Parameters.from_catalog(bundle_id, components=..., version=None)` loads an
  immutable installed selection. Omitting `version` selects the greatest
  installed integer version for the exact bundle ID.
- `Parameters.from_bundle(path, components=...)` loads and selects an exact
  strict version-1 bundle directory.
- `Parameters.from_records(...)` builds an immutable caller-owned selection
  from typed records with explicit source and validity-domain references.
- `Parameters.create_bundle_template(...)` creates a non-overwriting authoring
  skeleton. Profiles are `neutral`, `associating`, `electrolyte`, and
  `associating-electrolyte`; `correlations` is optional.

The selected component order is the order for every state, property, fugacity,
and activity result. There is no mutable catalog registry, hidden physical
default, generic parameter overlay, or compatibility alias.

### State and reference entry points

`Mixture(parameters)` resolves one immutable model. `Mixture.state` accepts
exactly one of `P` or `rho`, along with `T` and ordered `x`; an optional phase
applies only to pressure closure. `State` exposes the residual methods
`ares`, `ures`, `hres`, `gres`, `sres`, `cvres`, and `cpres`, each with
`basis="molar"` by default and `basis="reduced"` as the only alternate basis.
`State.fugacity` is the component-order fugacity projection when defined.

`Mixture.reference(T, P, solvent_x, phase="liquid")` constructs the admitted
liquid electrolyte reference. `reference.state(molality)` returns a `State`
with `State.activity` attached. This path remains limited to its source-bound
one-formula-unit reference contract.

## Consequences

The README and Python API guide can be copied as frontend documentation without
translating names or selecting a backend. Existing scientific bounds, source
provenance, catalog IDs, result ordering, structured errors, density-root
diagnostics, and residual property bases remain explicit.

The change is intentionally not a package release or a scientific promotion.
CppAD remains the sole nonlinear production derivative authority. Native C++
continues to own the hard-chain, dispersion, association, Debye--Huckel, Born
SSM+DS, and density equations. Phase equilibrium, discovery, regression, Ceres,
cross-package validation, and reaction/species equilibrium remain outside the
EOS.

The native consumer boundary remains governed by ADR 0001 and
`promotion-0015-provider-native-sdk-v1`. This decision does not broaden that
boundary or grant authority to local append-only candidates.

## Rejected alternatives

- retaining aliases for pre-0.2 class or method names;
- exposing backend-selection flags or mutable parameter dictionaries;
- adding a second serialized input or result representation; and
- silently filling missing source records or physical values.

## Evidence and references

The installed catalog table and source/evidence limits remain in
`docs/python-api.md`. Historical EOS capability, migration, and authority
statements remain in `CONTEXT.md`, the README, ADR 0001, and the accepted
promotion receipts. The full scientific equation owner is
`docs/science/full-provider-eos.md`.
