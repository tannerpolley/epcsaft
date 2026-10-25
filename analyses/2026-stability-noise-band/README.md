# Final stability decision at rounding noise (#159)

Question: the final local stability check on a converged reacting liquid passed or failed on the sign of a
reduced-Hessian eigenvalue near ±1e-13. What separates a real negative curvature from rounding?

**Inputs.**
- **Case.** The packet-3 MEA nine-species liquid at 318.15 K, with the absorber's fixed density-basis
  ln K = (−605.348, −587.8, −948.2312, −25.1014, −298.6415). This is the #154/#157 worker's case.
- **Sweep.** 150 seeded states: P from 1.05 to 1.12 bar, loading from 0.03 to 0.56 on 9.7292821 mol MEA
  and 76.9751925 mol water.
- **Rejection check.** Two identical ethane pseudo-species A ⇌ B (K = 1) at 5 MPa, held at the symmetric
  stationary point.
- **Wheels.** Before: 684b213a and 4105e2a4, the issue's two wheels. Main 99e333cd gives rows identical to 4105e2a4's and is not retained; 4105e2a4's ethane-pair rows equal 684b213a's, so only `pair-before-684b213a` is kept. After: c3f6e676.
  Every JSON records its SHA-256.

```bash
OMP_NUM_THREADS=1 .venv/bin/python analyses/2026-stability-noise-band/scripts/run.py LABEL
```

## Diagnosis (`states-after.csv`; spectrum rebuilt in NumPy from the public Jacobian and Hessian)

In the log-composition chart a species' direction has curvature of order its amount n x_i.

| Quantity (150 states) | min / median / max |
|---|---|
| x_CO2 | 3.8e-17 / – / 5.4e-7 (126 states below 1e-12) |
| ‖H‖₂ = largest eigenvalue (the volume coordinate) | 2.24e3 / 2.42e3 / 2.58e3 |
| condition of the reduced Hessian | 4.8e7 / 7.8e15 / 5.5e18 |
| ε‖H‖ | 5.0e-13 / 5.4e-13 / 5.7e-13 |
| n ε ‖S‖, with S = the split magnitude \|H_objective\| + \|H_constraint\| | 8.3e-11 / 8.4e-11 / 8.4e-11 |

- **Raw sign.** At trace CO2 the smallest raw eigenvalue is about n_CO2, 1e-15 to 1e-13. That is below its own rounding.
  - The objective and multiplier-weighted constraint diagonals cancel at stationarity. On water their magnitudes sum to 3.4e4, while the diagonal entry is 23.
  - The computed sign is noise. On states that pass, raw minima run from −3.5e-13 to +7.5e-11.
- **Before.**
  - Each wheel repeats its own decisions.
  - 684b213a rejects 23 of 150 states, and 4105e2a4 rejects 29.
  - The two issue wheels disagree on 30 states.
  - Every rejected raw minimum lies in [−5.5e-13, −1e-16].

## Criterion (algorithms.md ALG-EQUILIBRIUM)

1. Equilibrate by the congruence D = diag(S)^−½. The reduced Hessians of H and DHD have the same inertia (Sylvester's law), so equilibration changes the conditioning, not the answer.
2. Compare λ_min(DHD) on the nullspace of AD with the first-order rounding band τ = ε(n‖DSD‖_F + 2κ(AD)‖DHD‖_F):

| λ_min(DHD) | Outcome |
|---|---|
| > τ | stable |
| < −τ | rejected: "not positive definite" |
| within ±τ | typed `Stop.singular`: "is indeterminate" |

`reduced_hessian_minimum` stays the raw eigenvalue.

| After (c3f6e676) | Value |
|---|---|
| outcomes | 150 stable; 20 repeated solves identical |
| equilibrated minimum | 1.36e-3 to 1.47e-3 (maximum 1.00) |
| band | 2.1e-14 to 6.4e-14 |
| minimum margin λ̃/τ | 2.2e10 |

The margin is ten orders of magnitude, so the decision does not depend on rounding differences between builds.

**Genuine negatives stay rejected** (`pair-after.csv`). The symmetric state's equilibrated minimum is
−0.047 to −0.65 for k_ij = 0.1 to 0.5 at 150–200 K, against bands of about 5e-15. Stable controls (k_ij = 0 and 0.1 at 200 K)
pass.

## Limits

- **Point error at degenerate states.** The band covers rounding of the assembled Hessian and nullspace. It does not cover the returned
  point's solve error.
  - The solve error generally exceeds the band: at a row residual near 1e-12 the equilibrated minimum moves by order 1e-13. It decides the outcome only when |λ̃| is below it.
  - That happens near a genuinely degenerate state. For the ethane pair at k_ij ≈ 0.128102889834 and 200 K, the symmetric curvature vanishes, and the 4 degenerate rows in `pair-after.csv` pass or fail on the returned point.
  - The passing rows stop at a row residual of 1.3e-12 with a raw minimum of +2.6e-13. The failing rows stop at 9e-16 with −4.3e-13 (equilibrated −9e-14), against a band of 5e-15.
  - Including propagated solve error in the band would need a Newton correction of the returned point. That is not implemented.
- **Indeterminate outcome not exercised.** No retained case reaches the band itself. The outcome can fire only when the solve error is at rounding level.
- **Amount-unit dependence.** The d = 1 fallback applies to every free phase-amount coordinate, so κ(AD) and the band depend on the amount unit. The inertia does not.
- **Rounding model.** The band is a first-order estimate, not a rigorous bound. Cancellation inside the objective curvature sum is not in the split.
