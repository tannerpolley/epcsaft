# Provider Repository Context

This clean repository is the production home for ePC-SAFT provider code. Its
optimization objective is minimum maintained surface subject to source-backed,
reproducible thermodynamic results.

`governance_doctrine_revision: 2`

Canonical local doctrine: `../ePC-SAFT-organization/GOVERNANCE.md`.

The accepted provider contains one typed parameter lifecycle and one scalar-
generic residual Helmholtz kernel for hard-chain, dispersion, association,
Debye--Huckel, and corrected Born SSM+DS contributions. It returns residual
chemical potentials, fugacity coefficients when defined, residual caloric
properties, pressure, and branch diagnostics for pressure-specified states.
The installed catalogs are Gross--Sadowski methane/ethane and propane records
and the bounded Figiel--Yu--Held water/methanol/ethanol alkali-halide source bundle. One reference-state
module owns x-process infinite-dilution solvation and transfer plus single-ion
and stoichiometric mean ionic activity on one explicit electroneutral formula-
unit path.

Pressure closure accepts exactly one of molar density or pressure. An
unlabelled pressure state must have one mechanically stable root. A genuine
multi-root state requires `phase="vapor"` or `phase="liquid"`; a phase hint
never renames a lone root.

Accepted authority receipt `promotion-0012-provider-full-candidate` promotes
this exact provider scope to `ePC-SAFT/ePC-SAFT`. Its immutable installed-wheel
proof includes 64 package tests and direct-data campaigns owned by
`ePC-SAFT-validation`. This is a production-authority decision, not a package
release or a claim beyond the recorded capability boundaries.

Phase equilibrium, phase discovery/global stability, reaction and species
equilibrium, multi-salt reference paths, regression, and ideal-gas caloric data
remain outside this repository's current capability. Migration receipt
`promotion-0015-provider-native-sdk-v1` alone determines authority for the ADR
0001 model-bound capsule; this documentation does not transfer it. The fixed
prefix evaluates pure, nonassociating phase Helmholtz values and CppAD
derivatives through third order in amount and volume. Its append-only entry
returns the same pure function's gradient and Hessian when `m`, `sigma`, and
`epsilon/k` are active coordinates for the first methane regression slice. It
has no phase identity, density root, target residual, Ipopt, Ceres, equilibrium
result, or fit result. Each consumer remains subject to its own installed-wheel
gate and promotion receipt.

The local Consumer Maturity Slice 2 prerequisite appends one candidate
methane/ethane neutral-mixture phase entry to the same capsule and one further
active-`k_ij` callback. The latter reuses the mixture result and returns the
same extensive ideal-plus-residual value and exact CppAD gradient/Hessian in
`(n1,n2,V,kij)`. It accepts any finite trial unlike-pair value; regression owns
bounds, scaling, targets, and fit acceptance. Neither entry is authoritative
until a separate promotion receipt is accepted, and neither adds phase,
density, equilibrium, or regression ownership.

The architecture and cost ratchet are recorded in `ARCHITECTURE.yaml`; the
scientific equation owner is `docs/science/full-provider-eos.md`.
