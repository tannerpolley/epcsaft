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
The installed catalogs are Gross--Sadowski methane/ethane and propane records,
the bounded Figiel--Yu--Held water/methanol/ethanol chloride/bromide source
bundle, a distinct aqueous Na+/K+/I- source bundle, and one distinct Khudaida
et al. water/ethanol/isobutanol/Na+/Cl- bundle at 293.15 K. One reference-state
module owns x-process
infinite-dilution solvation and transfer plus single-ion and stoichiometric mean
ionic activity on one explicit electroneutral formula-unit path. Its aqueous
single-solvent result also reports the practical molal osmotic coefficient
from the solvent-activity ratio on that pressure-state path.

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

Migration D-024 records a local append-only full-electrolyte phase-block
candidate in explicit component amounts and extensive volume. The callback
uses the sole Provider Helmholtz and association owners for every active
ideal, hard-chain, dispersion, association, Debye--Huckel, and corrected Born
SSM+DS contribution. It admits exact immutable Figiel selections and the
distinct Khudaida bundle and returns value, gradient, Hessian, pressure,
fingerprint, component order, and charges. It adds no phase discovery,
equilibrium, trace-species policy, fitting, or authority.

The D-024 volume-domain tails share one Provider packing owner with pressure
closure. The value convenience accepts the exact model fingerprint,
temperature, ordered normalized composition, and caller-owned packing-fraction
interval, returning the corresponding lower and upper molar volumes in
`m^3/mol`. The derivative tail accepts positive ordered amounts and extensive
volume and returns packing fraction plus exact CppAD gradient/Hessian in
`(n_1,...,n_C,V)`. Perdomo limits, constraint construction, search variables,
initialization, and acceptance remain Equilibrium-owned.

Migration D-023 additionally records a local implementation candidate for one
model-bound active-Born-diameter evaluation of the existing 298.15 K,
1 bar water single-ion x-process solvation reference. It is limited to exact
current-catalog Li+, Na+, K+, Cl-, and Br- model selections and returns
only the scalar solvation Gibbs energy and its first derivative. It does not
authorize targets, Ceres, fit policy, catalog persistence, other
electrolyte parameters, or broader species. The exact conditional contract and
evidence gate live in
`docs/designs/2026-07-20-active-born-solvation-native-sdk-design.md`.
The implementation moves the existing four-molality composition,
schedule, convergence, and terminal-selection algorithm behind one private
native reference-sequence owner. The public Python path delegates through
the existing `_core` module with unchanged behavior, and the callback
calls the same owner directly; neither path may duplicate the algorithm.

The architecture and cost ratchet are recorded in `ARCHITECTURE.yaml`; the
scientific equation owner is `docs/science/full-provider-eos.md`.
