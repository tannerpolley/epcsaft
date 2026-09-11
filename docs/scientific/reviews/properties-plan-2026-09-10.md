# Properties milestone plan review

2026-09-10 — roadmap handoff review; no numerical qualification.

## Inspected state

The review started from clean commit
`f774b8bd6b3118a718f23ac414cb3074a9c601fd` on
`codex/contributions-qualification`. It inspected the Properties scope and
P1–P5 rows, canonical property equations, ALG-CALORICS,
ALG-ASSOCIATION, ALG-DENSITY, ALG-ACTIONS and ALG-REFERENCE, the current native
Properties owners and tests, and the mapped old thermochemistry and derivative
sources. No wheel or numerical model was run because this change only revises
the execution plan.

## Findings and disposition

The dependency order P1 → P2 → P3 → P4 → P5 is sound. The previous rows were
not sufficient implementation handoffs:

- P1 did not identify the unqualified `c_v^r` and `c_p^r` checks or the missing
  `c_p^r` conditioning basis.
- P2 combined zero-support association, density branch identity and pressure
  chains without naming the contracts that must be resolved before repair.
- P3 did not require an output/availability/reference matrix before adapting
  old thermochemistry.
- P4 did not reconcile requested action order with observable-specific energy
  derivative order, and its old derivative-header locator was incomplete.
- P5 did not state which integrated freshness, route-equivalence, availability
  and equivalent-work checks close Properties.

The implementation plan now contains one concise Resolved/Build/Check/Stop
handoff for every Properties row and marks P1 planning active. P1 is verify-first
and begins with zero runtime additions. P2–P4 remain explicitly design-gated;
P5 remains verify-only. The stale current-queue C3 state was also corrected to
match its qualified milestone and record.

Disposition: the roadmap is explicit enough to draft the P1 delivery issue and
assign its bounded design/verification work to a Luna max worker. It is not yet
Design Ready: the issue must record the identified revision and adoption
decision, resolve the `c_p^r` conditioning basis, and receive fresh independent
pre-Build review. No later Properties implementation is authorized until its
own named design gate and prerequisites are complete.

## Evidence limits and next action

The historical F1 correction remains valid preserved evidence only for its
recorded inputs and wheel. This review does not establish current-wheel caloric
agreement, settle P2–P4 scientific choices, qualify any Properties milestone,
or authorize Equilibrium work. Next, publish and independently review the P1
issue from the roadmap handoff; then run its smallest caloric falsifiers before
changing runtime code.
