# Provider Repository Agent Contract

Repository Profile: scientific-computing

Before work, read doctrine revision 2 at
`../ePC-SAFT-organization/GOVERNANCE.md`. Its future published home is
`https://github.com/ePC-SAFT/.github/blob/main/GOVERNANCE.md`. Local policy may
narrow that doctrine but must not contradict it.

This repository owns one provider lifecycle:

```text
ParameterBundle -> ParameterSet -> EPCSAFT.evaluate -> EosResult
```

The local candidate evaluates neutral, associating, and bounded ionic SSM+DS
models at explicit molar density or specified pressure. Native C++ owns the
hard-chain, dispersion, association, Debye--Huckel, Born SSM+DS, and density
equations. CppAD is the sole nonlinear production derivative authority. Pint
is the sole Python runtime dependency.

Accepted receipt `promotion-0009-provider-slice-1` still governs only the
published neutral explicit-density Slice 1 baseline. Later local commits are
not promoted, released, or authoritative outside this checkout until a new
receipt is accepted.

The provider owns parameters, EOS/state/property evaluation, density closure,
CppAD, and any admitted native consumer transport. It does not own phase
equilibrium or discovery, Ipopt, regression, Ceres, cross-package validation,
or reaction/species equilibrium. ADR 0001 selects the future native transport
but adds no SDK runtime until its separate evidence gate passes. Do not add
compatibility APIs, mutable scientific dictionaries,
backend selectors, hidden physical defaults, duplicate equation owners, or
runtime dependencies on the lab, migration repository, or sibling sources.
