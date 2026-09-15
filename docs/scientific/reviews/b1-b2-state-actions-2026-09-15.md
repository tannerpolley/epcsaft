# B1/B2 bounded state-action review - 2026-09-15

## Accepted boundary

Issues #56 and #57 own the reviewed public contract. Commit `9de30498` exposes
the existing Properties-owned `density_actions` and `pressure_actions` through
the native Python extension and the `Mixture` and `Pure` facades. Directions are
SI coordinate increments along a dimensionless path parameter; returned actions
retain the requested observable's unit. One to three ordered directions select
mixed derivative order one to three.

No thermodynamic equation, contribution, derivative algorithm, or density solve
changed. The binding opens the required diagnostic scope and forwards directly
to `HelmholtzModel`. The public result preserves the native optional value,
optional action, status, observable, and energy order.

## Numerical evidence

The installed-wheel test `engine/tests/eos/test_actions.py` compares public
residual-enthalpy actions against independent centered mixed stencils at two
steps for orders one through three in both density and pressure charts. It also
checks the central value, pure endpoint identity, missing-thermochemistry partial
results, support change, empty and duplicate requests, invalid order and width,
nonfinite directions, unbalanced composition directions, and ionic charge
imbalance.

The unchanged IDAES/Pyomo campaign was replayed from clean commit `9de30498` and
promoted under `analyses/2026-idaes-asl-adapter/results/`. It passed 72 scalar
derivative ladders, 216 whole-NLP directional ladders, 27 Hessian matrices, four
process cases, 54 conditioning evaluations, and seven independent references.
The ASL library SHA-256 remains byte-identical at
`8f83088b976066163a1672c422a1be40de85720a0ab9681828b6e2ce36de6d13`.

Independent code and scientific/API reviews found no material correctness,
compatibility, security, evidence-integrity, or ownership issue. The focused
public test passed three tests, and `bash tools/check-all.sh` passed 12 Python
tests, all four native suites, Ruff, mypy, packet validation, analysis layout,
and architecture. Bindings use 692 physical lines, above the 650 soft ceiling
but below the 800 hard ceiling; the explicit public enum and result mapping
accounts for the bounded overage.

## Branch consolidation

`idaes/main` is the sole local `idaes/*` branch and is based on the current
`origin/main` commit `fa65d132`. Rejected histories remain recoverable as annotated
tags: `archive/idaes/historical-prototype` at `c995e393`,
`archive/idaes/neutral-current-engine` at `0194d4eb`, and
`archive/idaes/pre-reset-recovery` at `448e5604`.

## Limits

Formal B1 and B2 remain open because Surface S4 is unfinished. Parameter
directions, thermochemistry-backed ideal and total actions, and public diagnostic
mapping are not exposed. The IDAES campaign qualifies its existing neutral
propane ASL slice only; it does not qualify mixture, electrolyte, equilibrium,
global stability, PROMMIS, or downstream application behavior.
