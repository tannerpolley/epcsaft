# Provider Repository Context

This clean repository is the production home for ePC-SAFT provider code. Its
optimization objective is minimum maintained surface subject to source-backed,
reproducible thermodynamic results.

`governance_doctrine_revision: 4`

Canonical local doctrine: `../ePC-SAFT-organization/GOVERNANCE.md`.

Provider owns both its package evidence and execution of assigned installed-
artifact campaigns in the sibling Validation repository. The task remains based
here; Validation remains the durable black-box evidence home and has no resident
worker. Migration serializes exact campaign subjects and dispatches a separate
reviewer only when required. This execution model does not change scientific
admission or runtime authority.

The accepted provider contains one typed parameter lifecycle and one scalar-
generic residual Helmholtz kernel for hard-chain, dispersion, association,
Debye--Huckel, and corrected Born SSM+DS contributions. It returns residual
chemical potentials, fugacity coefficients when defined, residual caloric
properties, pressure, and branch diagnostics for pressure-specified states.
Callers can create `user-provided` bundles through
`ParameterBundle.from_records(...)`. `to_path(...)` and `from_path(...)`
persist and reload the strict version-1 schema without a second input
representation. Callers keep provenance and validity domains explicit. The
purpose neither promotes values into the installed catalog nor validates their
scientific quality.
An exact selected model from a persisted `user-provided` bundle can use the
generic model-bound native phase, volume-domain, and packing callbacks. This
transport admission does not admit source-specific native reference or active-
parameter callbacks and does not make caller-owned values predictive.
`ParameterBundle.catalog()` exposes the immutable installed bundles for
discovery. State and ionic-activity results repeat their component order, and
state results retain canonical temperature, composition, and parameter
fingerprint. Structured state errors preserve the rejected field, reason, and
applicable units, limits, or tolerances. The distributed Python interface is
marked with `py.typed` and exposes literal types for phases, purposes, domains,
and parameter families; runtime validation remains the scientific authority.
The installed catalogs are Gross--Sadowski methane/ethane and propane records,
the bounded Figiel--Yu--Held water/methanol/ethanol chloride/bromide source
bundle, a distinct aqueous Na+/K+/I- source bundle, one distinct Khudaida et al.
water/ethanol/isobutanol/Na+/Cl- bundle at 293.15 K, and one distinct Ascani et
al. original-Born water/1-butanol/Na+/K+/Cl- bundle at 298.15 K. One reference-state
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
ideal, hard-chain, dispersion, association, Debye--Huckel, and model-bound Born
contribution. It admits exact immutable Figiel selections, the distinct
Khudaida bundle, and the distinct Ascani original-Born bundle and returns value,
gradient, Hessian, pressure,
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

Project Truss Provider leaves `ePC-SAFT/ePC-SAFT#1` and `#2` record the
next Figiel family. One model-bound
callback evaluates `ln(gamma_pm^m)` for the exact six current-catalog aqueous
alkali-halide models while accepting only the three local trial interactions
`(water--cation, water--anion, cation--anion)`. It returns the value and exact
total fixed-pressure first derivative from the existing Helmholtz,
association, density, and reference-sequence owners. The contract and evidence
gate live in
`docs/designs/2026-07-21-aqueous-active-kij-miac-native-sdk-design.md`.
Project Truss issue `ePC-SAFT/ePC-SAFT#29` appends typed batch tails for this
callback and the NaBr active-solvation-factor callback. A batch preserves row
order and exact scalar results while constructing the active model and
derivative reference sequence once per salt and parameter vector; the exact
catalog solvation-factor coordinate retains one additional batch-level
public-value parity sequence. Provider owns the
bounded CppAD row concurrency; Regression must reject a scalar-only artifact
instead of entering an unbounded per-row fit. This throughput transport changes
no equation, derivative order, admission boundary, fitted value, or authority.
Current status is
`PROVIDER_AQUEOUS_ACTIVE_KIJ_IMPLEMENTATION_CANDIDATE_PASSED_LOCAL_GATES_INDEPENDENT_REVIEW_PENDING`;
the callback is a local implementation candidate and no retained artifact or
downstream fit exists. Bounded independent subagent review changes sequencing
only; it does not transfer runtime, catalog, regression, scientific,
promotion, or publication authority.

Project Truss Provider issue `#22` adds the one missing staged-recovery seam:
an append-only NaBr `ln(gamma_pm^m)` value plus exact total fixed-pressure
first derivative with respect to active `f_water`. The implementation tapes
that coordinate in the existing CppAD Born owner and reuses the existing
association, density, fugacity, and electrolyte reference sequence. Its exact
contract is
`docs/designs/2026-07-22-aqueous-active-solvation-factor-native-sdk-design.md`.
Installed NaBr remains fingerprint-bound; only source-complete, ordered,
corrected-SSM+DS Figiel user bundles with exact applicability may use the
staged active callbacks. This is Provider transport evidence, not a fit,
catalog change, scientific acceptance, or authority transfer.

Project Truss Provider issue `#10` adds a separate append-only installed-SDK
candidate for charge-neutral reacting standard references. For an admitted
aqueous phase model with one explicit water reference carrier, it returns a
complete neutral stoichiometric basis and gauge-invariant reference
log-fugacity contractions for any number of neutral components, cations, and
anions. The exact current Figiel selections retain their existing admission;
full admitted catalog and authenticated `user-provided` phase models can use
the same tail. Unknown caller-owned applicability transports the dimensions
but fails value evaluation until finite source temperature and pressure bounds
exist. The same reference-sequence and pressure-state owners supply the values;
no single-ion standard, reaction constant, reaction solver, chemistry selector,
or catalog admission is introduced. The result explicitly reports that `T`,
`P`, composition, and parameter derivatives are unavailable.

The architecture and cost ratchet are recorded in `ARCHITECTURE.yaml`; the
scientific equation owner is `docs/science/full-provider-eos.md`.
