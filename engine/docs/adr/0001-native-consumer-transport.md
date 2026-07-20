# ADR 0001: Versioned model-bound capsule for native consumers

Status: accepted pure prefix authority determined by migration promotion
receipt 0015; append-only neutral-mixture and active-kij entries are local
candidates pending independent review and their own promotion receipts; a
model-bound active-Born solvation entry is design-only pending permanent-lab
review, implementation, evidence, and its own authority decision; consumer
gates separate

Date: 2026-07-17

## Context

The equilibrium package needs exact provider derivatives inside native solver
callbacks. Python remains the orchestration boundary but is not called from an
Ipopt hot loop. Cross-repository builds may consume only an installed provider
artifact; they may not compile provider sources, link private C++ symbols, or
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

The provider evaluates one phase at a time. It receives positive finite
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
`a^r` is the provider's residual molar Helmholtz owner. The reference adds only
a species-linear constant, which cancels between phases and under nonreactive
material balances. It is not an absolute reactive standard-state convention
and does not add ideal-gas caloric data.

With numerical coordinates ordered as `(n_mol, V_m3)`, the provider returns
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
variables, Ceres, or fit acceptance to the provider.
The regression problem and receipt must retain the exact trial parameter tuple,
units, transforms/scales, bounds, and immutable source fingerprint; the source
fingerprint returned by the provider never identifies a fitted candidate as a
new source-backed parameter record.

## Native boundary

The public Python function `epcsaft.native_sdk(model)` returns a model-bound
CPython capsule named `epcsaft.native_sdk.v1`. Its pointer addresses a
size-tagged pure-C function table containing:

- ABI version, table size, and result-structure size;
- one opaque immutable provider model context; and
- one fixed-parameter pure-phase evaluation function; and
- one additive pure-phase evaluation function with active `m`, `sigma`, and
  `epsilon/k` coordinates.

The local Slice 2 candidate appends immutable component count, mixture-result
size, and one neutral-mixture evaluation function without changing the
accepted prefix. The entry uses caller-owned amount, gradient, and Hessian
buffers. Its model-bound coordinate order is `(n_1,...,n_C,V)`, with component
order fixed by the selected parameter set, and it returns derivatives only
through order two. One further callback pointer reuses that result structure
for the exact methane/ethane coordinate order `(n_1,n_2,V,k_ij)`. Its single
active unlike-pair value is written symmetrically to both off-diagonal pair
positions before the canonical Helmholtz owner is taped.

The D-024 candidate appends the full-electrolyte phase callback and one
model-bound molar-volume-domain callback. For an exact fingerprint,
temperature, and composition in the table's component order, the latter uses

```text
K_eta(T,x) = (pi/6) (N_A/10^30) sum_i x_i m_i d_i(T)^3
eta        = K_eta rho_m
v_m        = K_eta/eta
```

from the same Provider packing owner used by pressure closure. Given
`0 < eta_min < eta_max < pi/(3 sqrt(2))`, it writes
`[K_eta/eta_max, K_eta/eta_min]` in `m^3/mol`. Packing limits are explicit
caller inputs; the Provider does not own Perdomo initialization or solver
policy. The callback uses a size-gated two-double caller buffer and no new
result structure. A bound-count mismatch returns `ABI_MISMATCH` without a
write; all failures are return-code authoritative because there is no second
status/error payload owner.

D-023 identifies a later, design-only electrolyte tail for the first Figiel
regression tracer. Its minimum contract is one model-bound callback for the
existing 298.15 K, 1 bar water x-process single-ion solvation free energy with
only the selected current-catalog ion Born diameter active. It returns one
value and its first total fixed-pressure derivative. The exact model allowlist,
reference path, implicit density sensitivity, proposed append-only ABI, and
evidence gate are specified in
`docs/designs/2026-07-20-active-born-solvation-native-sdk-design.md`. This ADR
records neither implementation nor authority for that tail.

Before that callback is implemented, the existing Python-owned four-molality
algorithm must be refactored behind one private Provider-native owner,
`epcsaft::detail::electrolyte_reference_sequence` in
`src/epcsaft/_native/electrolyte_reference.hpp`. It alone owns electroneutral
composition construction, the exact schedule, convergence adjudication, and
terminal fugacity selection. The public `reference_state` path delegates to it
through the existing private `_core` module with unchanged public behavior;
the future callback calls it directly. Neither path may retain a second copy
or expose this refactor as a new public abstraction.

The capsule owns the model context. A consumer keeps a Python reference to the
capsule for as long as it calls the table; ordinary `Py_INCREF`/`Py_DECREF`
provides retention, so a second native lifetime protocol is unnecessary. The
provider wheel installs the declaration-only C header
`epcsaft/include/epcsaft/native_sdk_v1.h`. The header contains no equations,
templates, implementation source, CMake fragments, C++ standard-library types,
exceptions, Eigen objects, or pybind11 objects.

A single immutable provider model is evaluated independently for every phase
state. Equilibrium may call the same capsule for liquid-like and vapor-like
`(n,V)` points, but the provider never constructs liquid, vapor, phase, or
mixture-solver objects. The equilibrium package owns phase count, ordering,
variables, constraints, Ipopt assembly, initialization, continuation,
certification, and results.

The accepted pure entries reject associating, ionic, mixture,
out-of-source-domain, out-of-EOS-domain, and nonfinite inputs before returning
a nominal result. The parameterized entry additionally rejects nonpositive
active parameters. The local mixture candidate admits only the two-component
source-backed Gross--Sadowski methane/ethane model in its bound component order
and rejects pure, associating, ionic, malformed-buffer, and nonpositive states.
The active-`k_ij` callback additionally rejects nonfinite trial values but does
not impose fit bounds or scaling. Its returned fingerprint identifies the
immutable source model from which the fit problem descends; candidate values
remain explicit coordinates and never masquerade as a new source record.

The active-Born design admits only the five exact water/ion/counterion
selections and fingerprints listed in its design. It fixes the active ion at
component index 1, requires a positive finite trial diameter, and rejects any
other solvent, component order, source fingerprint, ion, or state. H+, I-,
sulfate, transition-metal/V species, generic electrolyte parameters, targets,
and fit policy remain outside the contract.

## Compatibility

Consumers import the capsule through the public Python function and compile
only against the installed declaration header. They do not import
`epcsaft._core`, private headers, provider source paths, or native object
layouts. A consumer accepts only the exact capsule major name and verifies the
advertised result size before calling an entry. It checks `table_size` against
the minimum prefix needed for that entry before reading appended fields; it
does not require equality with the newest whole-table size.

Backward-compatible entries may be appended within version 1. A changed field
meaning, unit, coordinate order, tensor layout, reference basis, or ownership
rule requires a new capsule major version.

## Admission gate

The provider-side evidence requires:

- isolated installed-wheel availability with no sibling source import;
- exact capsule name, ABI version, table/result sizes, lifetime, and failure
  behavior;
- provider artifact hash, exact admitted parameter fingerprints, and rejection
  of synthetic or structurally similar unadmitted models;
- thermodynamic agreement of `P` and `mu/(RT)` with the public provider state;
- directional gradient, Hessian-vector, and third-tensor contraction checks
  against an independent numerical oracle;
- fixed/parameterized value agreement at the source methane parameters and
  independent five-coordinate gradient and Hessian-vector checks;
- methane, ethane, and propane in-domain coverage; and
- absence of provider implementation sources from the provider wheel.

Those provider tests pass for the fixed prefix and append-only parameterized
tail. Migration receipt `promotion-0015-provider-native-sdk-v1` is the sole
authority record; this ADR does not transfer authority. Each consumer gate must
build its wheel against the installed provider wheel, verify capsule/header
compatibility, and prove that no provider source or duplicate EOS symbol enters
the consumer artifact. Equilibrium implementation, regression target/Ceres
implementation, cross-wheel evidence, and promotion remain separate bounded
slices and receipts.

The local mixture candidate adds provider-owned checks for explicit-density
value and pressure agreement, centered directional gradient and Hessian
agreement, component permutation, accepted-prefix offsets, buffer ownership,
source domain, and unsupported-model rejection. These tests and an isolated
wheel are candidate evidence only; they do not amend receipt 0015 or authorize
consumer use. The further active-`k_ij` evidence checks two real
methane/ethane states, fixed-callback parity at the catalog value, the complete
four-coordinate gradient and Hessian direction, its exact appended offset,
finite-only trial admission, and binary/neutral/nonassociating scope.

The design-only active-Born evidence gate requires source-owner parity for all
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

- Neutral nonassociating mixture families beyond the exact source-backed
  methane/ethane candidate require separate scope and evidence.
- Binary-interaction families beyond the exact local methane/ethane active-
  `k_ij` candidate, association, electrolyte, and other fitted-parameter
  families require separate parameterized phase contracts and evidence. They
  are not implied by the pure three-parameter entry or this bounded candidate.
- The D-023 active-Born solvation design does not admit runtime implementation,
  fitted parameters, additional ions, solvation-factor or dielectric
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

- Separate liquid and vapor provider models: duplicates one EOS identity and
  assigns phase meaning to the wrong owner.
- Source/CMake SDK: creates a second compiled owner of provider equations.
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
