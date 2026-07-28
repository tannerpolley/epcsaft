# Provider Repository Agent Contract

Repository Profile: scientific-computing

Before work, read doctrine revision 4 at
`../ePC-SAFT-organization/GOVERNANCE.md`. Its future published home is
`https://github.com/ePC-SAFT/.github/blob/main/GOVERNANCE.md`. Local policy may
narrow that doctrine but must not contradict it.

This repository owns one state/property lifecycle:

```text
ParameterBundle -> ParameterSet -> EPCSAFT.evaluate -> EosResult
```

It also owns one narrow native consumer boundary:

```text
EPCSAFT -> native_sdk(model) -> one stateless pure phase block
```

The accepted provider evaluates neutral, associating, and bounded ionic SSM+DS
models at explicit molar density or specified pressure. Native C++ owns the
hard-chain, dispersion, association, Debye--Huckel, Born SSM+DS, and density
equations. CppAD is the sole nonlinear production derivative authority. Pint
is the sole Python runtime dependency.

Accepted receipt `promotion-0012-provider-full-candidate` governs this exact
full-provider scope. It transfers production authority but does not publish a
package release or broaden any capability beyond the receipt.

The provider owns parameters, EOS/state/property evaluation, density closure,
CppAD, and any admitted native consumer transport. It does not own phase
equilibrium or discovery, Ipopt, regression, Ceres, cross-package validation,
or reaction/species equilibrium. Migration receipt
`promotion-0015-provider-native-sdk-v1` is the sole authority record for the
ADR 0001 capsule; this documentation does not transfer that authority. Its
fixed prefix owns only pure nonassociating phase values and derivatives in
amount and volume. Its append-only parameterized entry owns only exact
pure-phase derivatives with respect to `m`, `sigma`, and `epsilon/k`; target
residuals, Ceres, and fit acceptance remain regression-owned. Each consumer
still requires its own cross-wheel gate. Do not add
compatibility APIs, mutable scientific dictionaries,
backend selectors, hidden physical defaults, duplicate equation owners, or
runtime dependencies on the lab, migration repository, or sibling sources.

## Agent home and Validation work

- The Provider agent remains rooted in this repository and is the sole owner of
  Provider implementation, parameter resolution, EOS/state/property behavior,
  derivatives, density closure, and the native SDK.
- When Migration assigns a bounded Provider campaign, this same agent may use a
  clean checkout or worktree of the sibling Validation repository to author and
  execute that campaign. Its task home does not move to Validation.
- Validation work must use exact immutable installed public artifacts. It must
  not import this source checkout, use private seams, copy Provider equations,
  or add implementation code to Validation.
- The Provider agent owns the campaign execution and evidence commit, then
  returns to Provider work. It does not create or wait on a resident Validation
  task, accept its own promotion, or broaden Provider authority from its own
  evidence.
- Migration serializes Validation writers and dispatches any required distinct
  review after the exact stable subject is available.

The append-only neutral methane/ethane mixture value, gradient, and Hessian
entry and its further append-only active-`k_ij` and active-`l_ij` entries are
local provider candidates pending their own review and promotion receipts.
Each active entry accepts one finite trial unlike-pair value and returns
derivatives in `(n1,n2,V,pair_parameter)` without owning fit bounds or
scaling. These candidates change
no runtime authority and admit no association, ions, density solve, phase
identity, third derivatives, equilibrium, regression targets, or fitting.

The same local candidate tail advertises three scalar pure-component
capabilities for neutral, nonassociating, constant-diameter models:
`segment_count`, `segment_diameter`, and `dispersion_energy_over_k`. One
callback projects the accepted joint pure-parameter Hessian onto
`(amount, volume, active_parameter)`. It enables one-family-at-a-time
regression without replacing the accepted joint three-parameter callback or
moving bounds, target data, Ceres, or fit acceptance into Provider.
The capability set also advertises the already-owned exact first-derivative
callbacks for current-catalog single-ion Born diameter and NaBr water
solvation factor. These direct-observable descriptors add no EOS equation,
target, residual, fit policy, or broader model admission.
The same exact five-ion water domain exposes one model-bound active
ionic-region relative permittivity through the existing solvation-Gibbs
reference sequence. The published value is fixed model input rather than a
fitted target; the callback adds derivative transport, not catalog authority
or a shared multi-model fit.
The append-only tail also exposes the existing infinite-dilution reference
sequence with exact first derivatives in the three ordered `k_ij` coordinates
for the source-backed Figiel methanol/K+/Br- and ethanol/Na+/Cl- selections.
It owns no constructed targets, residuals, bounds, Ceres, or fitted values.

Migration D-024 additionally records local append-only full-electrolyte and
molar-volume-domain candidates. The phase callback evaluates the complete
active Helmholtz block in explicit amounts and extensive volume. One domain
callback maps caller-owned packing-fraction limits to composition- and
temperature-specific molar-volume bounds; a further callback returns packing
fraction and exact CppAD gradient/Hessian in those phase coordinates. Pressure
closure and both domain callbacks share the same Provider packing owner. None
owns equilibrium, phase discovery, solver policy, fixed phase volumes,
fitting, or authority.

Project Truss Provider leaves `ePC-SAFT/ePC-SAFT#1` and `#2` record one
aqueous Figiel MIAC value/derivative callback with exactly
three active interactions in `(water--cation, water--anion, cation--anion)`
order. Its sole Provider design owner is
`docs/designs/2026-07-21-aqueous-active-kij-miac-native-sdk-design.md`.
Status is
`PROVIDER_AQUEOUS_ACTIVE_KIJ_IMPLEMENTATION_CANDIDATE_PASSED_LOCAL_GATES_INDEPENDENT_REVIEW_PENDING`.
The exact implementation remains a local candidate until bounded independent
subagent review. It does not authorize regression, target ownership, catalog
writes, organic-solvent breadth, a generic parameter overlay, runtime
authority, or promotion.

Project Truss Provider issue `#22` appends one NaBr active-water-solvation-
factor MIAC callback to support the frozen staged Figiel recovery. Its sole
Provider design owner is
`docs/designs/2026-07-22-aqueous-active-solvation-factor-native-sdk-design.md`.
The callback consumes the existing CppAD, density, fugacity, association, and
reference-sequence owners and returns only the value and exact total fixed-
pressure first derivative. Source-complete user-provided staged bundles are
narrowly admitted to the three existing family callbacks without adding a
generic parameter overlay or catalog persistence.
