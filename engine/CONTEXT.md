# EOS Repository Context

This repository is the production home for the `epcsaft` equation-of-state
package. Its
optimization objective is minimum maintained surface subject to source-backed,
reproducible thermodynamic results.

`governance_doctrine_revision: 7`

Canonical doctrine:
`../ePC-SAFT-governance/GOVERNANCE.md`.
The Governance topology manifest is
`../ePC-SAFT-governance/manifests/management.json`.

## Parameter and evidence boundary

Reusable parameter bundles and reusable raw or organized property datasets are
owned by `ePC-SAFT-data`. Paper-, table-, and figure-specific campaign extracts
may remain in Validation when they are evidence for one campaign rather than a
reusable scientific packet. EOS ships neither category as runtime package data.

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
aggregate Gross--Vrabec DD/QQ/DQ35 point-multipole, Debye--Huckel, and corrected
Born SSM+DS contributions. Physical dipole and quadrupole moments are optional,
source-unit inputs for neutral components; the public state reports only
`ares_polar`. It returns residual
chemical potentials, fugacity coefficients when defined, residual caloric
properties, pressure, and branch diagnostics for pressure-specified states.
The 0.2.0.dev0 Python frontend exposes `Parameters.from_dictionary(...)`,
`Parameters.from_bundle(...)`, `Parameters.create_bundle_template(...)`, and
the convenience constructor `Mixture.from_dictionary(...)`.
The dictionary route uses a strict fixed-unit, versioned schema; the bundle
route loads an explicit caller-supplied path. Neither route imports, fetches, or
discovers Data. Native SDK capabilities are admitted from resolved topology,
required fields, component order, and declared validity, never paper identity.
`Mixture(parameters).state(...)` returns immutable `State` values, and
`Mixture.reference(...).state(molality=...)` attaches an immutable `Activity`
projection to the resulting state. State and activity results repeat their
component order, and state results retain canonical temperature, composition,
and parameter fingerprint. Structured state errors preserve the rejected
field, reason, and applicable units, limits, or tolerances. The distributed
Python interface is marked with `py.typed`; runtime validation remains the
scientific authority. No catalog API or compatibility alias is retained. One
reference-state module owns x-process
infinite-dilution solvation and transfer plus single-ion and stoichiometric mean
ionic activity on one explicit electroneutral formula-unit path. Its aqueous
single-solvent result also reports the practical molal osmotic coefficient
from the solvent-activity ratio on that pressure-state path.

Historical Migration decision D-028 binds one additional private native
reference candidate for the exact
ordered water/H3O+/OH- mechanical model of Held, Cameretti, and Sadowski (2008)
at 298.15 K and 1 bar. The explicit record set closes component identity,
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
The append-only fixed-topology association candidate exposes the ordinary
`sigma_ij^3 kappa` mechanics without adding topology selection. Its model-
bound descriptor publishes component and site identities, positive site
multiplicities, allowed continuous slots, parameter/topology/artifact
fingerprints, temperature domain, basis, derivative order, and tensor layout.
A caller may activate any nonempty ordered subset of the advertised `m`,
`sigma`, `epsilon/k`, site-pair association-energy, and site-pair
association-volume slots. The callback returns exact value, gradient, and
complete row-major Hessian coordinates in `(n,V,requested slots...)` with
exact caller-owned capacities. Unknown, duplicate, malformed, nonfinite,
nonpositive, or fingerprint-mismatched requests fail closed. It adds no scalar
association callback, 2B-specific path, fitted values, topology inference,
regression policy, parameter persistence, or production authority.
The same finite descriptor set now binds the existing first-derivative
callbacks for resolved water single-ion Born diameter and NaBr water
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
explicitly supplied Li+, Na+, K+, Cl-, and Br- model selections and returns
only the scalar solvation Gibbs energy and its first derivative. It does not
authorize targets, Ceres, fit policy, parameter persistence, other
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
callback evaluates `ln(gamma_pm^m)` for the exact six source-defined aqueous
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
baseline solvation-factor coordinate retains one additional batch-level
public-value parity sequence. EOS owns the
bounded CppAD row concurrency; Regression must reject a scalar-only artifact
instead of entering an unbounded per-row fit. This throughput transport changes
no equation, derivative order, admission boundary, fitted value, or authority.
Current status is
`PROVIDER_AQUEOUS_ACTIVE_KIJ_IMPLEMENTATION_CANDIDATE_PASSED_LOCAL_GATES_INDEPENDENT_REVIEW_PENDING`;
the callback is a local implementation candidate and no retained artifact or
downstream fit exists. Bounded independent subagent review changes sequencing
only; it does not transfer runtime, parameter persistence, regression, scientific,
promotion, or publication authority.

Historical Project Truss EOS issue `#22` adds the one missing staged-recovery seam:
an append-only NaBr `ln(gamma_pm^m)` value plus exact total fixed-pressure
first derivative with respect to active `f_water`. The implementation tapes
that coordinate in the existing CppAD Born owner and reuses the existing
association, density, fugacity, and electrolyte reference sequence. Its exact
contract is
`docs/designs/2026-07-22-aqueous-active-solvation-factor-native-sdk-design.md`.
NaBr callback calls remain fingerprint-bound to their resolved model; ordered,
corrected-SSM+DS inputs with exact applicability may use the staged active
callbacks. This is EOS transport evidence, not a fit, parameter mutation,
scientific acceptance, or authority transfer.

The general Regression integration tail additionally exposes one
infinite-dilution ion-solvation Gibbs result with exact total fixed-pressure
derivatives in ordered `(solvent--cation, solvent--anion, cation--anion)`
`k_ij` coordinates. Admission is limited to the explicitly supplied
methanol/K+/Br- and ethanol/Na+/Cl- selections used by the retained Figiel
constructed-target ledger. The callback reuses the same active-`k_ij`
Helmholtz and electrolyte-reference owners; it does not own the equation-19
target construction, source interpretation, fitting, or parameter updates.

The same integration tail exposes the SSM+DS ionic-region relative
permittivity as one model-bound active coordinate for each admitted
resolved water single-ion model. It returns the existing x-process
solvation Gibbs value and its exact total fixed-pressure first derivative
through the same CppAD, density, fugacity, and reference-sequence owners.
Figiel reports the baseline value as fixed model input, so this callback proves
derivative and fitting mechanics only; it does not claim recovery of a
published fitted value, a shared cross-model parameter, or parameter authority.

The next append-only entry exposes the water component's ordinary
`relative_permittivity` as a separate active coordinate on those same five
model identities. Its solvation-Gibbs derivative includes the bulk dielectric,
Debye--Huckel, Born, fixed-pressure density/fugacity, and reference-sequence
chains. It must not be confused with the radial
`ionic_region_relative_permittivity`; neither fitted result mutates the
caller-owned parameter input.

Historical Project Truss EOS issue `#10` and its generic completion leaf `#75`
own the append-only installed-SDK candidate for charge-neutral reacting
standard references. Every admitted phase model with one or more neutral
molecular components, at least one cation, at least one anion, and a complete
finite source domain uses the same callback. It returns a complete
charge-nullspace basis and gauge-invariant reference log-fugacity contractions
for arbitrary charge magnitudes. The salt-free reference composition is
equimolar over all neutral components. Basis rows and ionic anchors are ordered
by exact component identity, so equivalent component permutations reproduce
the same rows after column remapping. Unknown caller-owned applicability
transports no evaluable neutral-reference capability. The same
reference-sequence, pressure-state, and CppAD owners supply values and admitted
pressure/parameter derivatives; no single-ion standard, reaction constant,
reaction solver, chemistry selector, or source-identity admission is
introduced.

An explicit Held aqueous-electrolyte record set uses that same generic
neutral-reference tail and the existing fixed-pressure density closure.
Consumers obtain its pure-water liquid molar volume from the ordinary public
pure-water state; no second standard-reference callback or density owner is
retained.

The architecture and cost ratchet are recorded in `ARCHITECTURE.yaml`; the
scientific equation owner is `docs/science/full-provider-eos.md`.
