# Phase-Equilibrium Documentation Authority

Status: canonical package overview

This document describes current ownership and the approved hard-cutover
destination. It does not claim that planned capability is already implemented
or scientifically admitted.

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

`discovery` determines whether the feed is one stable phase or splits into two
or three liquid/vapor phases and reports the stability/globality evidence
actually established.

`constrained_specification` accepts one to three declared liquid/vapor roles,
solves only that topology, and certifies only that claim. A constrained success
does not prove that another topology is absent.

Hints may supply initialization, continuation, composition, density, or phase
ordering. They may reduce work but cannot alter equations, tolerances,
certificates, or the scientific claim.

## Inputs and results

Inputs make component/species order, feed, ensemble, charge, solvent
composition, ion-pair/stoichiometry, phase mode/roles, reference convention,
and optional `ActiveParameterSet` explicit. No literal component name,
fingerprint, paper, campaign, or fixed position carries behavior.

Results contain ordered phase states and amounts, compositions, densities,
equilibrium quantities, balance/stationarity residuals, certificate scope,
conditioning, deterministic diagnostics, and structured failure. When active
parameters are requested, Equilibrium forms exact total solved-state
sensitivities from EOS CppAD partials by implicit linear algebra in the same
active-set order.

## Current implementation versus destination

The checkout may still contain `tp_flash`, HELD/HELD2 symbols, source bounds,
fingerprint dispatch, binary result fields, or family-limited provider
parameter transport while the cutover is in progress. Those are inventory,
not accepted architecture. They must be removed atomically when the general
owner and its compact numerical evidence satisfy the canonical plan.

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
