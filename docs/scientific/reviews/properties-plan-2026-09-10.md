# Properties milestone plan review

2026-09-10 — roadmap handoff review; no numerical qualification.

## Inspected state

The review started from clean commit
`f774b8bd6b3118a718f23ac414cb3074a9c601fd` on
`codex/contributions-qualification`. It inspected the Properties scope and
P1–P5 rows, canonical property equations, ALG-CALORICS,
ALG-ASSOCIATION, ALG-DENSITY, ALG-ACTIONS and ALG-REFERENCE, the current native
Properties owners and tests, and the mapped old thermochemistry and derivative
sources. The final C4 wheel
`sha256:0624f7334d18d76b412888c929e5c09ddcb00129010b2bbc47a85d7c092ee3a5`
supplied a bounded design probe; this is not P1 qualification.

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

The `c_p^r` conditioning question is resolved without a new tolerance: the
central value is available only for finite positive
`D=Z+rho Z_rho=P_rho/(RT)` and a finite quotient. Numerical qualification also
requires a converged same-branch fixed-pressure enthalpy-derivative ladder;
near-zero positive `D` is retained as conditioning stress outside the claim.
This preserves the existing state behavior and avoids an unsupported universal
cutoff. General action conditioning remains P4.

The design probe used centered temperature steps `{0.1, 0.03, 0.01, 0.003} K`.
At the last step, the absolute differences between the analytic field and
`d h^r/dT|P` or `d u^r/dT|rho` were:

| State | `D` | `c_p^r` difference J/(mol K) | `c_v^r` difference J/(mol K) |
|---|---:|---:|---:|
| Ethane, 200 K, 217800.69494692178 Pa, vapor | 0.8973490716 | 9.91e-10 | 7.15e-11 |
| Ethane, 200 K, 217800.69494692178 Pa, liquid | 13.6585137441 | 6.65e-9 | 5.25e-11 |
| Nine-species MEA state, 313.15 K, 1e5 Pa, liquid | 29.4595601695 | 7.76e-9 | 3.00e-9 |
| Water/Na+/Cl- state, 298.15 K, 1e5 Pa, liquid | 32.8861097005 | 1.02e-8 | 3.92e-9 |

Both final ladder steps satisfy
`abs(model-reference) <= 1e-8 max(1, abs(reference))` for these states. At an
ethane fixed-density stress state `(280 K, 3000 mol/m3)`, `D=0.0916154` and the
finite `c_p^r=181.677 J/(mol K)` is deliberately not a general accuracy claim;
at `(280 K, 4000 mol/m3)`, `D=-0.0777101` and `c_p^r` is unavailable. These
probes justify the handoff cases and availability rule only; the worker must
implement retained comparison rows and run them from the final P1 wheel.

Disposition: [P1 delivery issue #14](https://github.com/tannerpolley/ePC-SAFT-greenfield/issues/14)
publishes the identified-revision adoption decision, exact cases, ladders,
criteria, change budget, commands, evidence and stop rules for a Luna max
worker. Fresh independent pre-Build review is the remaining Design Ready gate.
No later Properties implementation is authorized until its own named design
gate and prerequisites are complete.

## Evidence limits and next action

The historical F1 correction remains valid preserved evidence only for its
recorded inputs and wheel. This review does not establish current-wheel caloric
agreement, settle P2–P4 scientific choices, qualify any Properties milestone,
or authorize Equilibrium work. Next, publish and independently review the P1
issue from the roadmap handoff; then run its smallest caloric falsifiers before
changing runtime code.
