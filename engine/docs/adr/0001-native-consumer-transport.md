# ADR 0001: Versioned capsule transport for native consumers

Status: accepted design; runtime not yet implemented

Date: 2026-07-17

## Context

The equilibrium package needs exact provider derivatives inside native Ipopt
callbacks. Python calls are adequate for orchestration but not for that hot
loop. Cross-repository builds may consume only an installed provider artifact;
they may not compile provider sources, link private C++ symbols, or copy the
EOS kernel. A stable public C++ ABI would add a larger compatibility promise
than the first consumer requires.

The first consumer is pure, nonassociating methane, ethane, or propane
saturation at fixed temperature. It uses phase amount and volume as explicit
variables. Density is therefore `n/V`; no density root or implicit density
sensitivity belongs in this first transport contract.

## Decision

The provider will expose `epcsaft.native_sdk()` as a public Python discovery
function returning a CPython `PyCapsule` named `epcsaft.native_sdk.v1`. The capsule contains a
size-tagged C function table. It is a runtime call boundary, not a public C++
ABI or a source SDK.

Version 1 contains only:

- ABI version and function-table size fields;
- creation of an opaque, reference-counted model handle from an immutable
  provider parameter snapshot;
- retain and release operations for that handle; and
- evaluation of one fixed-temperature, pure, nonassociating phase Helmholtz
  block at positive finite amount `n` and volume `V`.

The phase-block result contains the dimensionless phase objective `F/(RT)`,
its gradient in `(n,V)`, its symmetric Hessian, and its symmetric third-order
tensor. It also contains pressure and chemical potential as diagnostics, the
parameter fingerprint, and an explicit status/error code. Array order and
units are fixed by the versioned C structure.

The provider records one scalar-generic defining phase function with CppAD and
derives every returned nonlinear tensor from that function. Iterative density
or association solver history is not taped. The first contract rejects
associating, ionic, mixture, out-of-domain, and nonfinite inputs before
evaluation.

The equilibrium package owns the NLP variables, constraints, Ipopt adapter,
and equilibrium result. It may add exact linear transformations of provider
tensors, but it does not reproduce provider Helmholtz equations or compile
provider implementation sources.

## Compatibility

Consumers discover the capsule through the public Python function; they do not
import `epcsaft._core` or private headers. A consumer accepts only the exact
major ABI name and verifies the advertised table and result-structure sizes
before calling any entry. Backward-compatible entries may be appended within
version 1. A changed field meaning, unit, variable order, or ownership rule
requires a new capsule major version.

No source paths, CMake fragments, templates, Eigen types, C++ standard-library
types, exceptions, or pybind11 objects cross the function table. The installed
provider wheel remains the sole compiled owner of provider equations.

## Admission gate

This ADR does not add a public symbol or claim an SDK capability. Implementation
requires a separate bounded slice that proves:

- isolated wheel-to-wheel build with no sibling source checkout;
- exact capsule name, size, version, ownership, and failure behavior;
- provider and consumer artifact hashes and parameter fingerprints;
- directional gradient, Hessian-vector, and third-tensor contraction checks
  against an independent high-precision oracle;
- methane, ethane, and propane in-domain coverage; and
- absence of provider source files or duplicate EOS symbols in the equilibrium
  artifact.

Only after that evidence passes may `ARCHITECTURE.yaml` list a native SDK entry
or the equilibrium repository begin its pure-saturation runtime slice.

## Rejected alternatives

- Source/CMake SDK: creates a second compiled owner of provider equations.
- Public shared-library or C++ ABI: broader platform and compatibility burden
  than the admitted consumer needs.
- Private pybind imports or native object layouts: unversioned and bypass the
  public ownership boundary.
- Python callbacks inside the Ipopt hot loop: retain orchestration overhead and
  do not establish a native derivative contract.
- Density-root derivative transport for the first slice: unnecessary because
  amount and volume are explicit NLP variables.
