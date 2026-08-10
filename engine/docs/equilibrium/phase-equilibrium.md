# Phase-Equilibrium Documentation Authority

Status: protected HELD transition; candidate not promoted

This document describes the implemented nonreactive phase-equilibrium boundary
and the evidence required before the incumbent may be retired.

## Authority

The active transition decision is
[protect HELD until the phase-equilibrium promotion gate passes](adr/0002-protect-held-until-phase-equilibrium-promotion.md).
The destination remains
[one general nonreactive phase-equilibrium owner](adr/0001-one-general-phase-equilibrium-owner.md),
but its atomic-retirement sequence is superseded. The complete target design is
in the
[General Phase-Equilibrium Plan](plans/2026-08-07-general-phase-equilibrium.md).
`CONTEXT.md` defines shared language and `ARCHITECTURE.yaml` records the current
machine-readable boundary.

Historical research, source locators, results, and receipts keep their original
HELD/HELD2 names as provenance. During this protected transition the same names
also identify the exact incumbent implementation whose evidence remains live.

## Ownership

| Problem | Current route | Boundary |
| --- | --- | --- |
| Pure-fluid saturation | EOS `Pure.saturation` | Pure vapor--liquid boundary and exact selected-parameter sensitivities |
| Nonreactive incumbent | Equilibrium `tp_flash` | Exact retained HELD/HELD2 discovery workflow and certificates |
| Nonreactive candidate | Equilibrium `phase_equilibrium` | HELD-derived amount-feed interface and declared-role admission over the retained result |
| Homogeneous chemical equilibrium | Equilibrium `chemical_equilibrium` | One homogeneous reacting state over explicit species, reactions, balances, and charge |
| Reactive multiphase equilibrium | Deferred | Shared state/phase/parameter/sensitivity seam only; no runtime route |
| Solid equilibrium | Excluded | Requires a separate physical model and plan |

`tp_flash` and `phase_equilibrium` are two explicit interfaces over one
HELD/HELD2 numerical core. Neither is an alias or fallback, and neither silently
substitutes the other. The candidate must evolve this core; an unrelated
Gibbs-minimization solver cannot claim incumbent lineage.

## Implemented candidate modes

`discovery` runs the retained HELD or HELD2 calculation and reports its exact
bounded claim. It does not upgrade the incumbent's `not_guaranteed` globality
certificate.

`constrained` currently declares one or two liquid/vapor roles, runs the
retained calculation, classifies returned phases through the recorded
compressibility-factor limits, and accepts only an equal role multiset. A
mismatch fails closed with the retained state as the best candidate. This is a
declared-role constraint over the retained discovery result; it is not yet an
arbitrary constrained-topology solve. Three-role requests are rejected before
the incumbent runs, so the candidate cannot complete the target
one-to-three-phase owner claim.

Hints are type- and domain-validated but are not applied to HELD/HELD2 work
ordering. Diagnostics report `hints_applied=False`. Applying hints without
changing equations, tolerances, candidates, certificates, or results remains a
promotion requirement.

## Inputs and results

The incumbent retains its exact normalized-feed input, source-domain admission,
typed `PhaseState`, HELD/HELD2 diagnostics, certificates, and selected-parameter
sensitivity. These are protected scientific behavior, not a compatibility
facade.

The candidate accepts explicit feed amounts, temperature, pressure, mode,
roles, hints, and optional `ActiveParameterSet`. It rescales the normalized
incumbent result to the declared total amount, preserves phase compositions,
volumes, pressures, chemical potentials, objective, parameter fingerprint, and
sensitivity order, and adds typed roles, component amounts, packing fraction,
and transition diagnostics. A candidate failure never becomes an incumbent
success automatically.

## Numerical boundary

The restored neutral HELD and strong-electrolyte HELD2 sources are the exact
pre-CAP-09 files. Neutral and charged selection inside that incumbent remains
the retained model-topology behavior, including the audited source-domain gate
and HELD2 ten-step controller. The candidate does not change those equations,
start sets, tolerances, local optimizers, physical certificates, or exact
solved-state sensitivity route.

Liquid and vapor labels use explicit compressibility-factor bounds with a
fail-closed gap. This classification is an adapter claim only; it does not
rewrite or preselect HELD2 candidates, is not route-native phase-identity
evidence, and is insufficient for promotion. Candidate phases retain incumbent
native order so phase-indexed sensitivities cannot be desynchronized. Current
accepted HELD2 public results contain one or two phases. Three-phase constrained
solving remains unimplemented and is rejected explicitly.

The unrelated native CAP-09 scratch phase solver has been removed. The exact
`tp_flash`, HELD/HELD2 sources, diagnostic payloads, certificate audits, and
compact certificate evidence are restored and remain protected until
promotion.

## Promotion evidence

Engine retains compact numerical checks for HELD/HELD2 certificates and exact
incumbent/candidate state parity. Validation replays all eight accepted
Khudaida Figure 2 feeds through the installed wheel. It compares independent
`tp_flash` and `phase_equilibrium` solves at the incumbent's historical state,
physical, pressure, and KKT tolerances and requires the candidate's descriptive
paper error to be no worse beyond numerical tolerance.

This proves recovery and present lineage; it does not prove the missing
arbitrary constrained one-to-three-phase behavior. Retirement additionally
requires that behavior, route-native neutral and charged phase-role evidence,
applied-hint invariance, exact multiphase sensitivity/order evidence,
deterministic failure semantics, independent phase-science and
solver-architecture review, and explicit user approval bound to one accepted
receipt.

`chemical_equilibrium` remains separate. Its arbitrary
species/reaction/balance/charge compiler and solver are unchanged by this
transition.
