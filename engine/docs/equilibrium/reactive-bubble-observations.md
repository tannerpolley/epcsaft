# Reactive bubble/VLE observations

The admitted reactive bubble family solves a fixed-topology, fixed-temperature
incipient-vapor problem. The application declares the reacting liquid problem,
ordered neutral vapor support, pressure interval and starts, observation rows,
and immutable local continuation anchor. Charged components cannot enter the
vapor support.

`ProviderNonidealVapor` is the primary model. `IdealGasVaporComparison` is an
explicit comparison type and is never selected as a fallback. A measured
partial pressure row evaluates `y_i P`; phase equilibrium continues to use the
Provider fugacity.

The owner eliminates the liquid reaction block at each pressure and solves the
scalar bubble closure. For the coexistence boundary `B(x, theta)`, it uses

```text
F_P = B_x x_P - 1
F_theta = B_theta|x + B_x x_theta
dP*/dtheta = -F_theta / F_P
```

Liquid amount derivatives are normalized by the quotient rule. Every reported
liquid, vapor, and pressure output receives its direct derivative plus its
pressure-propagation term. Source-reference pressure dependence is inside
`x_P`; direct Provider parameter effects enter `B_theta|x` exactly once.

An evaluated result binds the exact certified chemical `(T,n,V)` state to the
coexistence liquid root by component order, normalized composition, molar
volume, phase role, active trial, topology, and branch identity. It also binds
the ordered vapor support and partner Provider fingerprint. The result is a
certified local equilibrium only. Search completeness and global equilibrium
are not claimed.

The owner fails closed with the public Non-Evaluable Trial vocabulary for no
incipient vapor, unresolved branches, branch loss, support loss or active-bound
ambiguity, Provider/domain rejection, topology mismatch, source-reference
failure, and singular or ill-conditioned sensitivity. It never returns a
numeric penalty for those states.

Regression consumes this family through the existing `ParameterFit`, `fit`,
and `FitResult` path. Equilibrium owns values, exact total Jacobians, physical
branch identity, certification, and opaque receipt validation. Regression owns
targets, residual transforms, weights, affine coordinates, bounds, caller
starts, optimizer accounting, and parameter-space optimizer clusters.

Admission evidence is recorded in
`governance/receipts/cap-13-reactive-bubble-vle-admission.json`.
