# Phase-Equilibrium Documentation Authority

Status: canonical package overview

This document describes the implemented nonreactive phase-equilibrium owner and
its deliberately bounded scientific claim.

## Authority

The active decision is
[one general nonreactive phase-equilibrium owner](adr/0001-one-general-phase-equilibrium-owner.md).
The complete implementation sequence and acceptance criteria are in the
[General Phase-Equilibrium Plan](plans/2026-08-07-general-phase-equilibrium.md).
`CONTEXT.md` defines shared language and `ARCHITECTURE.yaml` records the
machine-readable ownership boundary.

Superseded source-specific designs and execution plans have been removed.
Primary-source research, citations, source locators, raw data, receipts,
results, and still-relevant numerical evidence remain provenance; none selects
a runtime route.

## Ownership

| Problem | Owner | Boundary |
| --- | --- | --- |
| Pure-fluid saturation | EOS `Pure.saturation` | Pure vapor--liquid boundary and exact selected-parameter sensitivities |
| Nonreactive phase equilibrium | Equilibrium `phase_equilibrium` | Discovery or constrained specification for one to three liquid/vapor phases |
| Homogeneous chemical equilibrium | Equilibrium `chemical_equilibrium` | One homogeneous reacting state over explicit species, reactions, balances, and charge |
| Reactive multiphase equilibrium | Deferred | Shared state/phase/parameter/sensitivity seam only; no runtime route |
| Solid equilibrium | Excluded | Requires a separate physical model and plan |

The destination `phase_equilibrium` owner is not a facade over author-named
neutral and electrolyte controllers. Neutral systems use the same owner and
omit charge-specific operations because their topology has no charge.

## Problem modes

`discovery` evaluates the recorded finite one-/two-/three-phase role topologies
and deterministic start set. It returns the lowest certified local candidate
found in that scope and explicitly reports that globality is not established.
Search exhaustion without a certified candidate fails closed.

`constrained` accepts one to three declared liquid/vapor roles,
solves only that topology, and certifies only that claim. A constrained success
does not prove that another topology is absent.

Hints may supply initialization, continuation, composition, density, or phase
ordering. They reorder work but cannot alter equations, tolerances,
certificates, or the scientific claim.

## Inputs and results

Inputs make model component order, feed amounts, temperature, pressure, phase
mode/roles, hints, and optional `ActiveParameterSet` explicit. Charge belongs to
the resolved model topology. No literal component name, fingerprint, paper,
campaign, or fixed position carries behavior.

Results contain ordered phase states and amounts, compositions, densities,
equilibrium quantities, balance/stationarity residuals, certificate scope,
conditioning, deterministic diagnostics, and structured failure. When active
parameters are requested, Equilibrium forms exact total solved-state
sensitivities from EOS CppAD partials by implicit linear algebra in the same
active-set order.

## Implemented numerical boundary

The native owner minimizes fixed-temperature, fixed-pressure Gibbs energy in
amount--log-packing coordinates. Charged phases use an exact electroneutral
chart: neutral amounts, total charge equivalents, and positive/negative ionic
log-ratios. This keeps every provider evaluation in the electrolyte EOS domain;
neutral phases omit the chart operations naturally. Material constraints,
objective and constraint Hessians, and selected-parameter partials are composed
into the declared KKT `F_y` and `F_p` blocks.

Model temperature and ion-fraction domains are rejected before candidate
generation. For charged multiphase states, the reported chemical-potential
residual uses the charge-reduced basis obtained by eliminating one ionic
potential; a literal cross-phase ionic chemical-potential difference would be
gauge-dependent. Accepted states are strictly interior to every numerical
chart bound, so the equality-constrained KKT system is also the sensitivity
system. An implicit Jacobian is available only when its rank, condition number,
and dimension-scaled backward-error certificate all pass.

Liquid and vapor roles are physical claims, not initialization labels. The
certificate classifies them with the result's explicit compressibility-factor
limits and leaves a fail-closed gap between the liquid and vapor ranges. Hints
stable-sort the unchanged finite starts using role, composition, and density.
They do not remove a start or alter a certificate.

The former `tp_flash`, HELD/HELD2 runtime symbols, source bounds, fingerprint
dispatch, binary result fields, and their compatibility/audit tests have been
removed. Historical research, packets, results, and source locators retain
their original names as provenance.

`chemical_equilibrium` remains separate during this work. Its arbitrary
species/reaction/balance/charge compiler and solver should be preserved and
revised only for a demonstrated defect or the later reactive-phase design.

## Evidence boundary

Engine retains compact numerical checks for balances, stationarity, stability,
phase permutation, derivative correctness, conditioning, and one-/two-/three-
phase algorithm behavior. Named publication cases, broad phase envelopes, and
application acceptance run against the installed wheel in Validation. A match
to literature cannot substitute for local mathematical certificates, and a
local optimum cannot be reported as a global phase result.
