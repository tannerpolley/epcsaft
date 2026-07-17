# ADR 0001: Versioned model-bound capsule for native consumers

Status: implemented provider candidate; consumer gate pending

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

## Native boundary

The public Python function `epcsaft.native_sdk(model)` returns a model-bound
CPython capsule named `epcsaft.native_sdk.v1`. Its pointer addresses a
size-tagged pure-C function table containing:

- ABI version, table size, and result-structure size;
- one opaque immutable provider model context; and
- one pure-phase evaluation function.

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

Version 1 rejects associating, ionic, mixture, out-of-source-domain,
out-of-EOS-domain, and nonfinite inputs before returning a nominal result.

## Compatibility

Consumers import the capsule through the public Python function and compile
only against the installed declaration header. They do not import
`epcsaft._core`, private headers, provider source paths, or native object
layouts. A consumer accepts only the exact capsule major name and verifies the
advertised table and result sizes before calling an entry.

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
- methane, ethane, and propane in-domain coverage; and
- absence of provider implementation sources from the provider wheel.

Those provider tests pass in the local candidate and `ARCHITECTURE.yaml`
records only that bounded fact. The still-pending consumer gate must build an
equilibrium wheel against the installed provider wheel, verify capsule/header
compatibility, and prove that no provider source or duplicate EOS symbol enters
the consumer artifact. Equilibrium implementation, cross-wheel evidence, and
promotion remain a separate bounded slice and receipt.

## Deferred extensions

- Nonassociating mixtures require an extensive phase function in ordered
  `(n_1,...,n_C,V)` coordinates, including ideal mixing and only the derivative
  orders consumed by their approved formulation.
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
