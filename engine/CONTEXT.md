# EOS Repository Context

This repository is the production home for the `epcsaft` equation-of-state
package. Its
optimization objective is minimum maintained surface subject to source-backed,
reproducible thermodynamic results.

`governance_doctrine_revision: 5`

Canonical doctrine:
`../ePC-SAFT-management/GOVERNANCE.md`. During the local checkout transition,
`../ePC-SAFT-organization/GOVERNANCE.md` is the accepted Management alias.
The local Management manifest is
`../ePC-SAFT-organization/manifests/management.json` during that transition.

## Language

**Canonical progression source**:
A model-defining PC-SAFT or ePC-SAFT paper admitted by the reviewed
`ePC-SAFT Progression` source list. Membership makes its parameters eligible
for catalog review; it does not automatically admit a runtime bundle.
_Avoid_: any published paper, case-study source

**Source-complete catalog bundle**:
One immutable installed parameter bundle containing every model-table species
and every applicable pure, pair, topology, correlation, model, provenance, and
validity record needed to reproduce the source-defined selections.
_Avoid_: convenience bundle, exercised subset

**Validation campaign bundle**:
An explicit, hash-bound parameter input retained with a Validation campaign for
replay or qualification without becoming an EOS runtime default.
_Avoid_: Validation catalog, fallback catalog

EOS owns both its package evidence and execution of assigned installed-
artifact campaigns in the sibling Validation repository. The task remains based
here; Validation remains the durable black-box evidence home and has no resident
worker. Management serializes exact campaign subjects and dispatches a separate
reviewer only when required. EOS accesses Validation through a target-owned
worktree created on demand for the assigned campaign; its task-specific path
and branch remain owned by Validation. This execution model does not change
scientific admission or runtime authority.

The accepted EOS contains one typed parameter lifecycle and one scalar-
generic residual Helmholtz kernel for hard-chain, dispersion, association,
Debye--Huckel, and corrected Born SSM+DS contributions. It returns residual
chemical potentials, fugacity coefficients when defined, residual caloric
properties, pressure, and branch diagnostics for pressure-specified states.
The 0.2.0.dev0 Python frontend exposes `Parameters.from_catalog(...)`,
`Parameters.from_bundle(...)`, `Parameters.from_records(...)`, and
`Parameters.create_bundle_template(...)`; its template profiles are
`neutral`, `associating`, `electrolyte`, and `associating-electrolyte`, with
optional correlations. Callers keep provenance and validity domains explicit.
The purpose neither promotes values into the installed catalog nor validates
their scientific quality. The installed catalog is immutable package data and
its exact IDs and evidence domains are retained in `docs/python-api.md`.
An exact selected model from a caller-owned bundle can use the
generic model-bound native phase, volume-domain, and packing callbacks. This
transport admission does not admit source-specific native reference or active-
parameter callbacks and does not make caller-owned values predictive.
`Mixture(parameters).state(...)` returns immutable `State` values, and
`Mixture.reference(...).state(molality=...)` attaches an immutable `Activity`
projection to the resulting state. State and activity results repeat their
component order, and state results retain canonical temperature, composition,
and parameter fingerprint. Structured state errors preserve the rejected
field, reason, and applicable units, limits, or tolerances. The distributed
Python interface is marked with `py.typed` and exposes literal types for
phases, purposes, domains, and parameter families; runtime validation remains
the scientific authority. No compatibility aliases are retained.
The installed catalogs are Crespo et al. neutral 2B DEEA records, the
source-complete 78-species Gross--Sadowski 2001 bundle, the bounded
Figiel--Yu--Held water/methanol/ethanol chloride/bromide source bundle, a
distinct aqueous Na+/K+/I- source bundle, and one distinct Ascani et al.
original-Born water/1-butanol/Na+/K+/Cl- bundle at 298.15 K. The old
Gross--Sadowski methane/ethane and propane IDs remain only until the governed
catalog cutover removes them without aliases. One reference-state module owns x-process
infinite-dilution solvation and transfer plus single-ion and stoichiometric mean
ionic activity on one explicit electroneutral formula-unit path. Its aqueous
single-solvent result also reports the practical molal osmotic coefficient
from the solvent-activity ratio on that pressure-state path.

Historical Migration decision D-028 binds one additional private EOS candidate
for the exact
ordered water/H3O+/OH- mechanical model of Held, Cameretti, and Sadowski (2008)
at 298.15 K and 1 bar. The installed record set closes component identity,
charge, association, water--ion dispersion, excluded ion--ion dispersion,
constant-water-dielectric Debye--Huckel behavior, and zero Born contribution.
The append-only SDK tail declares the exact Helmholtz coordinate basis and
source pressure bounds and returns one neutral formula-unit reference sum plus
the pure-water reference from the existing native dilution owner. This is not
a reactive standard-state model: Equilibrium exclusively owns `K_w`, its
reference conversion, reaction/balance constraints, Ipopt, and certificates.
The candidate has `authority_effect: none`; no receipt, promotion, release, or
source-backed nonideal reaction claim follows from its presence here.

Pressure closure accepts exactly one of molar density or pressure. An
unlabelled pressure state must have one mechanically stable root. A genuine
multi-root state requires `phase="vapor"` or `phase="liquid"`; a phase hint
never renames a lone root.

Accepted authority receipt `promotion-0012-provider-full-candidate` promotes
this exact EOS scope. Its immutable installed-wheel
proof includes 64 package tests and direct-data campaigns owned by
`ePC-SAFT-validation`. This is a production-authority decision, not a package
release or a claim beyond the recorded capability boundaries.

Phase equilibrium, phase discovery/global stability, reaction and species
equilibrium, multi-salt reference paths, regression, and ideal-gas caloric data
remain outside this repository's current capability. Historical Migration
receipt
`promotion-0015-provider-native-sdk-v1` alone determines authority for the ADR
0001 model-bound capsule; this documentation does not transfer it. The fixed
prefix evaluates pure, nonassociating phase Helmholtz values and CppAD
derivatives through third order in amount and volume. Its append-only entry
returns the same pure function's gradient and Hessian when `m`, `sigma`, and
`epsilon/k` are active coordinates for the first methane regression slice. It
has no phase identity, density root, target residual, Ipopt, Ceres, equilibrium
result, or fit result. Each consumer remains subject to its own installed-wheel
gate and promotion receipt.

The local Consumer Maturity Slice 2 prerequisite appends one candidate neutral
nonassociating binary-mixture phase entry to the same capsule and further
active-`k_ij` and active-`l_ij` callbacks. Gross--Sadowski methane/ethane is
the first reference evidence, while a source-bound caller model carries its
own exact component and topology fingerprints. The same append-only tail adds
one scalar pure-parameter callback with derivative-ready descriptors for
`segment_count`, `segment_diameter`, and `dispersion_energy_over_k`. The
binary callbacks return exact CppAD derivatives in
`(n1,n2,V,pair_parameter)`; the pure callback projects the existing exact
joint pure Hessian into `(n,V,active_parameter)`. All accept finite trial
values while regression owns bounds, scaling, targets, and fit acceptance.
The scalar pure path supports independent family fits and does not replace
the accepted joint `(m,sigma,epsilon/k)` workflow. None is authoritative
until a separate promotion receipt is accepted, and neither adds phase,
density, equilibrium, or regression ownership.
For an exact neutral pure 2B model with one donor and one acceptor site, the
tail also exposes scalar association-energy/volume derivatives and two
distinct joint seven-coordinate pure-phase Hessians. The original joint
callback uses the Baygi--Pahlavanzadeh Eq. 4 temperature-dependent `d_ij^3`
basis. The appended ordinary callback and its exact typed descriptor use the
Gross--Sadowski `sigma_ij^3 kappa` basis owned by the public EOS. The
callbacks supply derivative transport for bounded Regression reconstructions
and do not admit other association topologies, targets, fitted values, or
catalog persistence.
The same finite descriptor set now binds the existing first-derivative
callbacks for current-catalog water single-ion Born diameter and NaBr water
solvation factor to their exact component identity, direct-observable
contract, units, fixed 298.15 K domain, topology, and parameter fingerprint.
This is metadata for callbacks already present in the SDK, not a new
thermodynamic implementation or a general electrolyte admission.

Historical Migration decision D-024 records a local append-only full-electrolyte phase-block
candidate in explicit component amounts and extensive volume. The callback
uses the sole EOS Helmholtz and association owners for every active
ideal, hard-chain, dispersion, association, Debye--Huckel, and model-bound Born
contribution. Its historical candidate record included exact immutable Figiel
selections, the now-retired Khudaida bundle, and the distinct Ascani original-Born
bundle and returns value,
gradient, Hessian, pressure,
fingerprint, component order, and charges. It adds no phase discovery,
equilibrium, trace-species policy, fitting, or authority.

The D-024 volume-domain tails share one EOS packing owner with pressure
closure. The value convenience accepts the exact model fingerprint,
temperature, ordered normalized composition, and caller-owned packing-fraction
interval, returning the corresponding lower and upper molar volumes in
`m^3/mol`. The derivative tail accepts positive ordered amounts and extensive
volume and returns packing fraction plus exact CppAD gradient/Hessian in
`(n_1,...,n_C,V)`. The append-only inverse geometry tail accepts a strict
`log(eta)` coordinate and ordered advertised active component parameters, then
returns positive extensive volume and its exact CppAD gradient/Hessian in
`(n_1,...,n_C,log(eta),theta_1,...,theta_P)`. It evaluates packing geometry
only: ion-fraction admission, phase Helmholtz evaluation, inner roots, limits,
constraint construction, search variables, initialization, and acceptance
remain outside this callback.

Historical Migration decision D-023 additionally records a local implementation candidate for one
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

Historical Project Truss EOS leaves `ePC-SAFT/ePC-SAFT#1` and `#2` record the
next Figiel family. One model-bound
callback evaluates `ln(gamma_pm^m)` for the exact six current-catalog aqueous
alkali-halide models while accepting only the three local trial interactions
`(water--cation, water--anion, cation--anion)`. It returns the value and exact
total fixed-pressure first derivative from the existing Helmholtz,
association, density, and reference-sequence owners. The contract and evidence
gate live in
`docs/designs/2026-07-21-aqueous-active-kij-miac-native-sdk-design.md`.
Historical Project Truss issue `ePC-SAFT/ePC-SAFT#29` appends typed batch tails
for this
callback and the NaBr active-solvation-factor callback. A batch preserves row
order and exact scalar results while constructing the active model and
derivative reference sequence once per salt and parameter vector; the exact
catalog solvation-factor coordinate retains one additional batch-level
public-value parity sequence. EOS owns the
bounded CppAD row concurrency; Regression must reject a scalar-only artifact
instead of entering an unbounded per-row fit. This throughput transport changes
no equation, derivative order, admission boundary, fitted value, or authority.
Current status is
`PROVIDER_AQUEOUS_ACTIVE_KIJ_IMPLEMENTATION_CANDIDATE_PASSED_LOCAL_GATES_INDEPENDENT_REVIEW_PENDING`;
the callback is a local implementation candidate and no retained artifact or
downstream fit exists. Bounded independent subagent review changes sequencing
only; it does not transfer runtime, catalog, regression, scientific,
promotion, or publication authority.

Historical Project Truss EOS issue `#22` adds the one missing staged-recovery seam:
an append-only NaBr `ln(gamma_pm^m)` value plus exact total fixed-pressure
first derivative with respect to active `f_water`. The implementation tapes
that coordinate in the existing CppAD Born owner and reuses the existing
association, density, fugacity, and electrolyte reference sequence. Its exact
contract is
`docs/designs/2026-07-22-aqueous-active-solvation-factor-native-sdk-design.md`.
Installed NaBr remains fingerprint-bound; only source-complete, ordered,
corrected-SSM+DS Figiel user bundles with exact applicability may use the
staged active callbacks. This is EOS transport evidence, not a fit,
catalog change, scientific acceptance, or authority transfer.

The general Regression integration tail additionally exposes one
infinite-dilution ion-solvation Gibbs result with exact total fixed-pressure
derivatives in ordered `(solvent--cation, solvent--anion, cation--anion)`
`k_ij` coordinates. Admission is limited to the catalog-bound
methanol/K+/Br- and ethanol/Na+/Cl- selections used by the retained Figiel
constructed-target ledger. The callback reuses the same active-`k_ij`
Helmholtz and electrolyte-reference owners; it does not own the equation-19
target construction, source interpretation, fitting, or catalog updates.

The same integration tail exposes the SSM+DS ionic-region relative
permittivity as one model-bound active coordinate for each admitted
current-catalog water single-ion model. It returns the existing x-process
solvation Gibbs value and its exact total fixed-pressure first derivative
through the same CppAD, density, fugacity, and reference-sequence owners.
Figiel reports the catalog value as fixed model input, so this callback proves
derivative and fitting mechanics only; it does not claim recovery of a
published fitted value, a shared cross-model parameter, or catalog authority.

The next append-only entry exposes the water component's ordinary
`relative_permittivity` as a separate active coordinate on those same five
model identities. Its solvation-Gibbs derivative includes the bulk dielectric,
Debye--Huckel, Born, fixed-pressure density/fugacity, and reference-sequence
chains. It must not be confused with the radial
`ionic_region_relative_permittivity`; neither fitted result persists to the
catalog.

Historical Project Truss EOS issue `#10` adds a separate append-only installed-SDK
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

The private Held water-ionization bundle uses that same generic neutral-reference
tail and the existing fixed-pressure density closure. Consumers obtain its
pure-water liquid molar volume from the ordinary public pure-water state; no
second standard-reference callback or density owner is retained.

The architecture and cost ratchet are recorded in `ARCHITECTURE.yaml`; the
scientific equation owner is `docs/science/full-provider-eos.md`.
