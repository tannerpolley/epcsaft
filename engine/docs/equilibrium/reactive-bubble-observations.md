# Reactive bubble/VLE observations

The admitted reactive bubble family solves a fixed-topology, fixed-temperature
incipient-vapor problem. The application declares the reacting liquid problem,
ordered neutral vapor support, pressure interval and starts, observation rows,
and immutable local continuation anchor. Charged components cannot enter the
vapor support.

`EOSNonidealVapor` is the primary model. `IdealGasVaporComparison` is an
explicit comparison type and is never selected as a fallback. A measured
partial pressure row evaluates `y_i P`; phase equilibrium continues to use the
EOS fugacity.

The declared phase count and roles are inputs: one certified reacting liquid
and one incipient vapor. The bubble owner does not run phase discovery or ask a
general coexistence solver to reconstruct the liquid. At each trial pressure,
the homogeneous owner supplies the exact certified liquid `(T,n,V)` state and
its pressure/parameter sensitivities. The bubble owner then solves only the
declared neutral vapor density branch and composition.

With vapor log-ratio coordinates `q`, the owner solves `g(q,P,theta)=0` for
independent chemical-potential differences and `h(q,P,theta)=0` for one
absolute chemical-potential equality. Exact implicit differentiation gives

```text
dq/dalpha = -g_q^-1 g_alpha
H_alpha = h_alpha + h_q dq/dalpha
dP*/dtheta = -H_theta / H_P
```

Liquid amount derivatives are normalized by the quotient rule. Every reported
liquid, vapor, and pressure output receives its direct derivative plus its
pressure-propagation term. Source-reference pressure dependence is inside
the liquid pressure sensitivity; direct EOS parameter effects enter each
liquid or vapor tangent exactly once.

An evaluated result carries the exact certified chemical `(T,n,V)` state as its
liquid phase and verifies component order, normalized composition, molar
volume, phase role, active trial, topology, and branch identity. It also binds
the ordered vapor support and partner EOS fingerprint. Declared starts
search only pressure/vapor closure outcomes; they never search for a second
liquid state. The result is a certified local equilibrium only. Search
completeness and global equilibrium are not claimed.

The owner fails closed with the public Non-Evaluable Trial vocabulary for no
incipient vapor, unresolved branches, branch loss, support loss or active-bound
ambiguity, EOS/domain rejection, topology mismatch, source-reference
failure, and singular or ill-conditioned sensitivity. It never returns a
numeric penalty for those states.

Regression consumes this family through the existing `ParameterFit`, `fit`,
and `FitResult` path. Equilibrium owns values, exact total Jacobians, physical
branch identity, certification, and opaque receipt validation. Regression owns
targets, residual transforms, weights, affine coordinates, bounds, caller
starts, optimizer accounting, and parameter-space optimizer clusters.

Admission evidence is recorded in
`governance/receipts/cap-13-reactive-bubble-vle-admission.json`.
