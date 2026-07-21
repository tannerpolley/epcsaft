# Provider Repository Agent Contract

Repository Profile: scientific-computing

Before work, read doctrine revision 3 at
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

The Provider task remains rooted in this repository and owns execution of its
assigned black-box campaigns in the sibling Validation repository. Use only
immutable installed public artifacts there; never import this source checkout or
private seams. Do not create a resident Validation task. Package-authored
evidence cannot self-accept promotion or broaden Provider authority.

The append-only neutral methane/ethane mixture value, gradient, and Hessian
entry and its further append-only active-`k_ij` entry are local provider
candidates pending their own review and promotion receipts. The active entry
accepts one finite trial unlike-pair value and returns derivatives in
`(n1,n2,V,kij)` without owning fit bounds or scaling. These candidates change
no runtime authority and admit no association, ions, density solve, phase
identity, third derivatives, equilibrium, regression targets, or fitting.

Migration D-024 additionally records local append-only full-electrolyte and
molar-volume-domain candidates. The phase callback evaluates the complete
active Helmholtz block in explicit amounts and extensive volume. One domain
callback maps caller-owned packing-fraction limits to composition- and
temperature-specific molar-volume bounds; a further callback returns packing
fraction and exact CppAD gradient/Hessian in those phase coordinates. Pressure
closure and both domain callbacks share the same Provider packing owner. None
owns equilibrium, phase discovery, solver policy, fixed phase volumes,
fitting, or authority.
