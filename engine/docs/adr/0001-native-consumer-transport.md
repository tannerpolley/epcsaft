# ADR 0001: Versioned model-bound capsule for native consumers

Status: accepted pure prefix authority determined by migration promotion
receipt 0015; append-only neutral-mixture, active-kij, active-lij, and scalar
pure-parameter entries are local candidates pending independent review and
their own promotion receipts; a generic pure fixed-topology association entry
is an issue-61 implementation candidate with no authority effect; a
model-bound active-Born solvation entry is design-only pending bounded
independent review, implementation, evidence, and its own authority decision;
an append-only charge-neutral reference-basis entry and its declared-solvent
source-reference transfer are Project Truss implementation candidates with no
authority effect;
consumer gates separate

Date: 2026-07-17

Parameter-catalog language in this historical transport decision is superseded
by ADR 0003. Current admission uses resolved topology, required fields,
component order, and applicability from explicit caller input.

## Context

The equilibrium package needs exact EOS derivatives inside native solver
callbacks. Python remains the orchestration boundary but is not called from an
Ipopt hot loop. Cross-repository builds may consume only an installed EOS
artifact; they may not compile EOS sources, link private C++ symbols, or
copy the EOS kernel. A public C++ ABI would add a larger compatibility promise
than the first consumer requires.

The first consumer is fixed-temperature pure, nonassociating methane, ethane,
or propane saturation. The local pressure/potential formulation directly
represents the unknown saturation boundary. It consumes third phase derivatives
because pressure and chemical potential are first derivatives of Helmholtz
energy and the admitted Ipopt route supplies their exact constraint Hessians.

Later specified-`T,P,N` flash and phase-discovery routes should instead use
direct total free-energy minimization with linear material balances where that
is the natural ensemble. Those routes normally consume only the phase Hessian.
The project does not force one formulation onto every equilibrium family merely
to standardize derivative order.

The first regression consumer fits the three pure methane segment parameters
against observed saturation states. It needs exact parameter derivatives but
does not need a second EOS kernel or density-root sensitivity path. The
regression formulation therefore lifts the liquid-like and vapor-like volumes
as row-local variables and obtains their pressure and chemical-potential
residual Jacobians from one parameterized phase Hessian.
Each row fixes `n = 1 mol`; otherwise simultaneous scaling of `n` and `V`
creates an unnecessary gauge while leaving density unchanged. The consumer
must reject coalesced phases and require finite, distinct volumes with positive
mechanical-stability slopes before treating the two local states as liquid-like
and vapor-like.

## Thermodynamic contract

The EOS evaluates one phase at a time. It receives positive finite
temperature `T [K]`, amount `n [mol]`, and volume `V [m^3]`; density is the
explicit identity `rho_m = n/V [mol m^-3]`. No density root, phase label,
continuation state, or implicit density sensitivity enters this contract.

For the first pure nonreactive contract, the defining scaled phase function is

```text
Phi(T,n,V) = A(T,n,V)/(R T n_ref)
           = n_hat * [ln(rho_m/rho_ref) - 1 + a^r(T,rho_m)]

n_ref   = 1 mol
rho_ref = 1 mol m^-3
n_hat   = n/n_ref
```

The logarithmic term is the complete fixed-temperature ideal density term and
`a^r` is the EOS's residual molar Helmholtz owner. The reference adds only
a species-linear constant, which cancels between phases and under nonreactive
material balances. It is not an absolute reactive standard-state convention
and does not add ideal-gas caloric data.

With numerical coordinates ordered as `(n_mol, V_m3)`, the EOS returns
`Phi`, its gradient, symmetric Hessian, and symmetric third-order tensor. The
same derivatives provide the diagnostics

```text
P [Pa]  = -R T n_ref * partial(Phi)/partial(V_m3)
mu/(RT) = partial(Phi)/partial(n_mol).
```

One scalar-generic defining function is recorded with CppAD. Every returned
nonlinear tensor comes from that function. The public result fixes coordinate
order, row-major tensor layout, units, reference basis, parameter fingerprint,
and explicit status/error codes.

The additive parameterized entry uses coordinates
`(n_mol, V_m3, m, sigma_angstrom, epsilon_over_k_kelvin)` and returns `Phi`, its
five-coordinate gradient, and its symmetric row-major Hessian. In particular,

```text
partial(P)/partial(q_j)       = -R T n_ref Phi_(V,q_j)
partial(mu/RT)/partial(q_j)   = Phi_(n,q_j).
```

These identities provide exact Jacobian columns for a regression-owned
residual without assigning target definitions, weights, bounds, nuisance
variables, Ceres, or fit acceptance to the EOS.
The regression problem and receipt must retain the exact trial parameter tuple,
units, transforms/scales, bounds, and immutable source fingerprint; the source
fingerprint returned by the EOS never identifies a fitted candidate as a
new source-backed parameter record.

## Native boundary

The public Python function `epcsaft.native_sdk(model)` returns a model-bound
CPython capsule named `epcsaft.native_sdk.v1`. Its pointer addresses a
size-tagged pure-C function table containing:

- ABI version, table size, and result-structure size;
- one opaque immutable EOS model context; and
- one fixed-parameter pure-phase evaluation function; and
- one additive pure-phase evaluation function with active `m`, `sigma`, and
  `epsilon/k` coordinates.

The local Slice 2 candidate appends immutable component count, mixture-result
size, and one neutral-mixture evaluation function without changing the
accepted prefix. The entry uses caller-owned amount, gradient, and Hessian
buffers. Its model-bound coordinate order is `(n_1,...,n_C,V)`, with component
order fixed by the selected parameter set, and it returns derivatives only
through order two. One further callback pointer reuses that result structure
for the exact bound-model coordinate order `(n_1,n_2,V,k_ij)`. Its single
active unlike-pair value is written symmetrically to both off-diagonal pair
positions before the canonical Helmholtz owner is taped. A final append-only
pointer applies the same contract to the distinct `l_ij` cross-diameter
coordinate.

One final scalar pure-parameter pointer reuses the mixture-sized result buffer
for exactly one active component parameter at a time. Its family code selects
`segment_count`, `segment_diameter`, or `dispersion_energy_over_k`; inactive
pure parameters retain the immutable model values. It calls the existing
five-coordinate pure-parameter owner and returns the exact projected value,
pressure, gradient, and row-major Hessian in `(n,V,active_parameter)`. A
descriptor binds the component identity, units, topology and source
fingerprints, derivative order, and source temperature domain. This projection
does not replace the joint five-coordinate callback used by accepted
three-parameter pure-saturation regression.

The issue-61 candidate appends one variable-width pure fixed-topology
association callback and descriptor. The descriptor binds exact component and
site order, site multiplicities, typed allowed slots, parameter, topology, and
installed-artifact fingerprints, temperature domain, ordinary
`sigma_ij^3 kappa` basis, per-slot open physical bounds, derivative order,
and row-major tensor layout. Its slots comprise `m`, `sigma`, `epsilon/k`,
and one association-energy and association-volume coordinate for each
canonical allowed site pair. A caller supplies any nonempty ordered subset and
exact-size gradient and Hessian buffers; output coordinates are
`(n,V,requested slots...)`. Unknown, duplicated, noncanonical, malformed,
nonfinite, nonpositive, or fingerprint-mismatched requests fail closed. EOS
retains the allowed-bond graph, inactive values, implicit site-fraction solve,
and exact first and second derivatives. No scalar association callback,
alternate basis, 2B wrapper, or compatibility path remains.

The historical D-024 candidate appended the full-electrolyte phase callback and
two model-bound volume-domain callbacks. Its unchanged phase callback recorded
exact immutable Figiel and Khudaida corrected-SSM+DS selections and the
distinct Ascani original-Born water/1-butanol/Na+/K+/Cl- selection. The
Khudaida selection is now retired and remains provenance only; it is not a
current catalog or SDK admission. The Ascani source
schema supplies explicit fixed `l_ij` cross-diameter records and precomputed
source-bound cross-association records; these are immutable model data, not an
active parameter ABI family. No callback, result structure, or table field is
added for the Ascani bundle. For an exact fingerprint, temperature, and
composition in the table's component order, the value convenience uses

```text
K_eta(T,x) = (pi/6) (N_A/10^30) sum_i x_i m_i d_i(T)^3
eta        = K_eta rho_m
v_m        = K_eta/eta
```

from the same EOS packing owner used by pressure closure. Given
`0 < eta_min < eta_max < pi/(3 sqrt(2))`, it writes
`[K_eta/eta_max, K_eta/eta_min]` in `m^3/mol`. Packing limits are explicit
caller inputs; the EOS does not own Perdomo initialization or solver
policy. The phase objective and packing derivative callback remain evaluable at
finite positive-volume constrained-optimizer trial points outside that accepted
packing interval. The molar-volume callback and consumer-owned constraints and
final-state gates remain the authoritative fail-closed physical contracts.
The callback uses a size-gated two-double caller buffer and no new
result structure. A bound-count mismatch returns `ABI_MISMATCH` without a
write; all failures are return-code authoritative because there is no second
status/error payload owner.

The adjacent derivative callback evaluates the same owner's

```text
eta(T,n,V) = (pi/6) (N_A/10^30) V^-1 sum_i n_i m_i d_i(T)^3
```

in coordinate order `(n_1,...,n_C,V)`. It returns the dimensionless value,
gradient, and row-major Hessian through exact-size caller buffers. CppAD is the
production derivative authority. Buffer-size mismatch returns `ABI_MISMATCH`;
fingerprint, model, source-temperature, positive-state, electroneutrality,
ion-fraction, and nonfinite-geometry checks fail before any output write. A
finite value above physical close packing is still returned because the
callback owns no packing limits or equilibrium constraint: Equilibrium supplies
the limits, consumes the value and derivatives in its NLP, and rejects an
infeasible final state.

D-023 identifies a local candidate electrolyte tail for the first Figiel
regression tracer. Its minimum contract is one model-bound callback for the
existing 298.15 K, 1 bar water x-process single-ion solvation free energy with
only the selected current-catalog ion Born diameter active. It returns one
value and its first total fixed-pressure derivative. The exact model allowlist,
reference path, implicit density sensitivity, proposed append-only ABI, and
evidence gate are specified in
`docs/designs/2026-07-20-active-born-solvation-native-sdk-design.md`. This ADR
records the implementation shape but no authority for that tail.

The former Python-owned four-molality algorithm is refactored behind one
private EOS-native sequence owner,
`epcsaft::detail::reference_limit_sequence` in
`src/epcsaft/_native/electrolyte_reference.hpp`. It owns the exact schedule,
convergence adjudication, and terminal selection; the electrolyte wrapper in
the same module owns formula-unit composition construction. The public
`reference_state` path delegates to that module
through the existing private `_core` module with unchanged public behavior;
the callback calls it directly. Neither path may retain a second copy
or expose this refactor as a new public abstraction.

The final local candidate tail exposes an equivalently complete charge-neutral
reference basis rather than accepting reaction rows. An admitted phase model
must contain one or more neutral molecular components, at least one cation, at
least one anion, and complete finite source temperature and pressure bounds.
The table reports `C-1` rows: one unit row for every neutral component, then
minimal formula-unit rows joining the canonical first cation to every anion and
every remaining cation to the canonical first anion. Components are ordered by
their exact identifiers for basis construction. This deterministic spanning
tree is a full-rank basis of the ordered charge-nullspace, including unequal
charge magnitudes, and is invariant to component permutation after columns are
mapped back by identity.

The callback evaluates neutral rows in the salt-free composition assigning an
equal fraction to every neutral component and evaluates every formula-unit row
through the canonical four-molality reference sequence. It returns the ordered
basis, dimensionless log-fugacity contractions, reference composition,
composition-weighted neutral molar mass, exact trial `T,P`, worst row
convergence error, source fingerprint, and typed Helmholtz basis. It
re-evaluates this reference at every admitted trial pressure; a source
standard-state pressure carried by an Equilibrium record is not substituted
for that trial pressure. An unknown caller-owned domain advertises no evaluable
neutral-reference capability. Equivalent reaction bases require only
consumer-side linear algebra; the EOS assigns no absolute single-ion
reference. The value entry retains its zero derivative bit mask, while the
append-only derivative entry supplies exact fixed-pressure and admitted active-
parameter derivatives through the existing reference-sequence and CppAD owners.

Issue #78 appends a declared-solvent operation without changing that basis or
the Provider Helmholtz gauge. Its source-reference declaration binds exact
component order, salt-free solvent mole fractions, activity convention,
standard molality, source-pressure provenance, and trial temperature and
pressure. Positive neutral fractions are finite-composition solvent
components; zero neutral fractions and every charge-neutral ionic row approach
infinite dilution on one deterministic x-process. For source solvent molar
mass M_S, the operation returns separate fugacity and activity-scale
contractions. A solvent unit row has activity-scale term ln(x_i^S), a
zero-composition neutral-solute row has ln(M_S m_standard), and an ionic
formula-unit row has the stoichiometric sum of that same term. Their sum is the
source-to-Provider transfer in the reported basis. No individual ionic
fugacity is returned.

The callback evaluates at trial pressure. Source reference pressure is retained
only in the immutable declaration receipt. Parameter, topology, component
order, source reference, domain, installed artifact, and Helmholtz-basis
fingerprints close the result identity. The first operation exposes no
reference derivatives; any nonempty derivative requirement fails with the
dedicated derivative-unavailable status before a consumer can optimize.
Success explicitly records admitted-domain and converged-limit statuses.
Unknown applicability or legacy artifacts without the new component-order and
domain identities do not advertise this appended operation; their existing
neutral-reference callback remains loadable.

The capsule owns the model context. A consumer keeps a Python reference to the
capsule for as long as it calls the table; ordinary `Py_INCREF`/`Py_DECREF`
provides retention, so a second native lifetime protocol is unnecessary. The
EOS wheel installs the declaration-only C header
`epcsaft/include/epcsaft/native_sdk_v1.h`. The header contains no equations,
templates, implementation source, CMake fragments, C++ standard-library types,
exceptions, Eigen objects, or pybind11 objects.

A single immutable EOS model is evaluated independently for every phase
state. Equilibrium may call the same capsule for liquid-like and vapor-like
`(n,V)` points, but the EOS never constructs liquid, vapor, phase, or
mixture-solver objects. The equilibrium package owns phase count, ordering,
variables, constraints, Ipopt assembly, initialization, continuation,
certification, and results.

The accepted pure entries reject associating, ionic, mixture,
out-of-source-domain, out-of-EOS-domain, and nonfinite inputs before returning
a nominal result. The parameterized entry additionally rejects nonpositive
active parameters. The local mixture candidate admits only a two-component
source-bound neutral nonassociating model in its bound component order and
rejects pure, associating, ionic, malformed-buffer, and nonpositive states.
Gross--Sadowski methane/ethane is its first reference evidence, not an
allowlist.
The active-`k_ij` and active-`l_ij` callbacks additionally reject nonfinite
trial values but do not impose fit bounds or scaling. Their returned
fingerprint identifies the
immutable source model from which the fit problem descends; candidate values
remain explicit coordinates and never masquerade as a new source record.

The active-Born candidate admits only the five exact water/ion/counterion
selections and fingerprints listed in its design. It fixes the active ion at
component index 1, requires a positive finite trial diameter, and rejects any
other solvent, component order, source fingerprint, ion, or state. H+, I-,
sulfate, transition-metal/V species, generic electrolyte parameters, targets,
and fit policy remain outside the contract.

## Compatibility

Consumers import the capsule through the public Python function and compile
only against the installed declaration header. They do not import
`epcsaft._core`, private headers, EOS source paths, or native object
layouts. A consumer accepts only the exact capsule major name and verifies the
advertised result size before calling an entry. It checks `table_size` against
the minimum prefix needed for that entry before reading appended fields; it
does not require equality with the newest whole-table size.

Backward-compatible entries may be appended within version 1. A changed field
meaning, unit, coordinate order, tensor layout, reference basis, or ownership
rule requires a new capsule major version.

## Admission gate

The EOS-side evidence requires:

- isolated installed-wheel availability with no sibling source import;
- exact capsule name, ABI version, table/result sizes, lifetime, and failure
  behavior;
- EOS artifact hash, exact admitted parameter fingerprints, and rejection
  of synthetic or structurally similar unadmitted models;
- thermodynamic agreement of `P` and `mu/(RT)` with the public EOS state;
- directional gradient, Hessian-vector, and third-tensor contraction checks
  against an independent numerical oracle;
- fixed/parameterized value agreement at the source methane parameters and
  independent five-coordinate gradient and Hessian-vector checks;
- fixed-topology 2B descriptor parity and a non-2B three-site, two-bond case
  whose value, every gradient entry, and every Hessian entry are checked
  against independently rebuilt public-model finite differences;
- arbitrary ordered subsets, nonunit multiplicity, exact ABI capacities,
  failed inner solves, and parameter/topology/artifact identity rejection;
- methane, ethane, and propane in-domain coverage; and
- absence of EOS implementation sources from the EOS wheel.

Those EOS tests pass for the fixed prefix and append-only parameterized
tail. Historical Migration receipt `promotion-0015-provider-native-sdk-v1` is the sole
authority record; this ADR does not transfer authority. Each consumer gate must
build its wheel against the installed EOS wheel, verify capsule/header
compatibility, and prove that no EOS source or duplicate EOS symbol enters
the consumer artifact. Equilibrium implementation, regression target/Ceres
implementation, cross-wheel evidence, and promotion remain separate bounded
slices and receipts.

The local mixture candidate adds EOS-owned checks for explicit-density
value and pressure agreement, centered directional gradient and Hessian
agreement, component permutation, accepted-prefix offsets, buffer ownership,
source domain, and unsupported-model rejection. These tests and an isolated
wheel are candidate evidence only; they do not amend receipt 0015 or authorize
consumer use. The further active-pair evidence checks two real
methane/ethane states, fixed-callback parity at the catalog value, the complete
four-coordinate gradient and Hessian direction, its exact appended offset,
finite-only trial admission, and binary/neutral/nonassociating scope.

The local active-Born evidence gate requires source-owner parity for all
five bound models, independent Figiel target points once Validation freezes
them, an independent high-precision corrected-Born check, centered-difference
checks of the complete first derivative, reference-value and derivative-limit
convergence, exact prefix/result offsets, no-write ABI mismatch behavior, and
installed-artifact negative space. It consumes no second derivative.
Source inspection must additionally prove that the schedule, composition
arithmetic, convergence comparison, and terminal selection occur only in the
single internal reference-sequence owner; the Python and SDK paths contain
delegation and result translation only.

## Deferred extensions

- Neutral nonassociating models with more than two components require separate
  scope and evidence.
- Binary-interaction families beyond the active `k_ij` and `l_ij` candidates,
  cross-association observations, topology selection, electrolyte, and other
  fitted-parameter families require separate contracts and evidence. The
  issue-61 entry varies continuous coordinates only inside one immutable,
  EOS-advertised neutral pure topology.
- The D-023 active-Born solvation candidate does not admit fitted parameters,
  additional ions, solvation-factor or dielectric
  parameters, other solvents/states, a generic parameter overlay, or catalog
  persistence.
- Eliminated association requires every mixed implicit sensitivity through the
  highest consumed order; alternatively a route may lift the true site state
  and its defining constraints.
- Electrolyte equilibrium requires approved electroneutral coordinates and an
  electrochemical-potential/Galvani convention.
- Reactive equilibrium requires source-backed standard-state chemical
  potentials; the nonreactive reference convention above is insufficient.

## Rejected alternatives

- Separate liquid and vapor EOS models: duplicates one EOS identity and
  assigns phase meaning to the wrong owner.
- Source/CMake SDK: creates a second compiled owner of EOS equations.
- Public shared-library or C++ ABI: broader platform and compatibility burden
  than the admitted consumer needs.
- Private pybind imports or native object layouts: unversioned and bypass the
  public ownership boundary.
- Python callbacks inside the native solver loop: retain orchestration overhead
  and do not establish a native derivative contract.
- Density-root derivative transport for the first slice: unnecessary because
  amount and volume are explicit.
- Direct free-energy minimization for the unknown pure saturation boundary:
  requires an artificial overall volume or phase split and is less direct than
  the pressure/potential boundary equations.
