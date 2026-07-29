# Python frontend migration: 0.2.0.dev0

This record aligns the active Python documentation with the v0.2 frontend. It
is not a package release, catalog update, scientific promotion, or native
consumer authority transfer.

The canonical path is:

```text
Parameters -> Mixture -> State
```

This is a breaking replacement:

| Pre-0.2 call | 0.2 call |
| --- | --- |
| `ParameterBundle.from_catalog(...).select(components)` | `Parameters.from_catalog(..., components=components)` |
| `ParameterBundle.from_path(...).select(components)` | `Parameters.from_bundle(..., components=components)` |
| `ParameterBundle.from_records(...).select(components)` | `Parameters.from_records(..., selected_components=components)` |
| `EPCSAFT(parameters)` | `Mixture(parameters)` |
| `model.evaluate(temperature=..., pressure=... or molar_density=..., mole_fractions=...)` | `mixture.state(T=..., P=... or rho=..., x=...)` |
| `EosResult` | `State` |
| `model.reference_state(...)` | `mixture.reference(...)` |
| `reference.at_molality(value)` | `reference.state(molality=value)` |

The pre-0.2 names are not exported and have no aliases.

Resolve parameters with `Parameters.from_catalog(...)`,
`Parameters.from_bundle(...)`, `Parameters.from_records(...)`, or
`Parameters.create_bundle_template(...)`. Template profiles are `neutral`,
`associating`, `electrolyte`, and `associating-electrolyte`; `correlations` is
optional.

Evaluate one state with `Mixture(parameters).state(T=..., P=..., x=...,
phase=...)` or `Mixture(parameters).state(T=..., rho=..., x=...)`. Read
residual properties through `State.ares`, `ures`, `hres`, `gres`, `sres`,
`cvres`, and `cpres`, using the molar basis by default or
`basis="reduced"`. Read component fugacity through `State.fugacity`.

For the admitted electrolyte path, build
`Mixture.reference(T, P, solvent_x, phase="liquid")`, call
`reference.state(molality=...)`, and read `State.activity`.

The selected component order remains the order for every composition and
result. Pint quantities remain required for dimensional inputs. Supply exactly
one of `P` and `rho`; use `phase` only with pressure closure. Source,
validity-domain, topology, and parameter fingerprints remain part of scientific
input identity. Branch on structured error `code` and `context` rather than
exception wording.

No compatibility aliases, mutable registries, hidden defaults, or alternate
equation owners are introduced by this migration. The pre-0.2 facade is
historical evidence only and is not a compatibility layer in `0.2.0.dev0`.

The catalog table and its evidence domains remain in `docs/python-api.md`.
CppAD remains the sole nonlinear production derivative authority; native C++
continues to own the hard-chain, dispersion, association, Debye--Huckel, Born
SSM+DS, and density equations. Provider ownership and the boundaries on phase
equilibrium, discovery, regression, Ceres, cross-package validation, and
reaction/species equilibrium are unchanged.

Accepted receipt `promotion-0012-provider-full-candidate` still governs the
full-provider scope. Migration receipt `promotion-0015-provider-native-sdk-v1`
alone determines authority for the ADR 0001 capsule; this migration does not
transfer it. Local append-only native candidates remain bounded evidence until
their own review and promotion receipts, and the private Held water-ionization
candidate retains `authority_effect: none`.
