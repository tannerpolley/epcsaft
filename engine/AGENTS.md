# Provider Repository Agent Contract

Repository Profile: scientific-computing

Before work, read doctrine revision 2 at
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

The append-only neutral methane/ethane mixture value, gradient, and Hessian
entry is a local provider candidate pending its own review and promotion
receipt. It changes no runtime authority and admits no association, ions,
density solve, phase identity, third derivatives, equilibrium, or regression.
