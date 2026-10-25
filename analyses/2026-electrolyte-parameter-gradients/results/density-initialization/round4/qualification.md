| Fit start | Arm | Optimizer iterations | Last accepted step: iteration, cost | First rejected step | Stopping test | Final gradient max norm or dual infeasibility |
|---|---|---:|---|---:|---|---:|
| F2 0 | cold | 24 | 18, 2.39e-19 | 12 | Parameter tolerance reached. Relative step_norm: 3.879579e-15 <= 1.000000e-14. | 1.48e-10 |
| F2 0 | seeded | 19 | 11, 6.27e-18 | 10 | Parameter tolerance reached. Relative step_norm: 1.511484e-16 <= 1.000000e-14. | 1.15e-10 |
| B4 1 | cold | 11 | 11, objective 0.227361067361 | none | Ipopt tolerance 1e-10: primal 2.12e-11, dual 4.45e-12 | 4.45e-12 |
| B4 1 | seeded | 11 | 11, objective 0.227361067359 | none | Ipopt tolerance 1e-10: primal 2.15e-11, dual 4.45e-12 | 4.45e-12 |

Cost is one half of the weighted squared residual sum. From the first rejected step onward, the F2 costs lie between 2.4e-19 and 2.1e-17 in both arms; accepted and rejected steps alternate there, and the trust region contracts until the relative step falls below the parameter tolerance. Both B4 starts in both arms end at the same minimax objective within 1.5e-12, although the two starts' coordinates differ by up to 0.086. Two starts with one objective point to an optimum that is a set rather than a point (inferred, not mapped), so the stop's position along it depends on the path.

| Fit start | First seeded evaluation | Optimizer point difference | Largest weighted residual difference | Row: kind, molality (mol/kg) | Largest relative Jacobian difference at the same point |
|---|---:|---:|---:|---|---:|
| F2 0 | 5 | 0 | 9.58e-10 | OsmoticCoefficient, 0.005 | 4.81e-11 |
| B4 0 | 6 | 0 | 1.28e-12 | row 15 | 1.37e-12 |
| B4 1 | 2 | 0 | 1.35e-12 | row 3 | 1.35e-12 |

The F2 seeded objective also depends on which previous trial seeded it: one point evaluated from two different seeds differs by 1.16e-09 in the weighted residual, while the cold objective repeats exactly. Mapped through the pseudo-inverse of the endpoint Jacobian (condition number 1027), the first seeded residual difference moves the least-squares point by up to 5.3e-11; the seeded and cold stops lie up to 1.8e-11 and 1.4e-11 from their own Gauss-Newton points.

| Controls (function, gradient, parameter or Ipopt tolerance) | Fit start | Residual / Jacobian evaluations, cold | Residual / Jacobian evaluations, seeded | Largest optimizer-endpoint difference | Difference / allowance |
|---|---|---:|---:|---:|---:|
| F2 [0, 0, 0] | F2 0 | 11 / 12 | 12 / 12 | 1.17e-10 | 11.7 |
| F2 [1e-15, 1e-14, 1e-14] | F2 0 | 25 / 15 | 20 / 12 | 8.86e-11 | 8.86 |
| F2 [1e-16, 1e-16, 1e-16] | F2 0 | 26 / 15 | 21 / 12 | 8.86e-11 | 8.86 |
| F2 [1e-20, 1e-20, 1e-20] | F2 0 | 26 / 15 | 21 / 12 | 8.86e-11 | 8.86 |
| F2 [1e-30, 1e-10, 1e-30] | F2 0 | 11 / 12 | 21 / 12 | 9.25e-11 | 9.25 |
| B4 1e-10 | B4 0 | 12 / 13 | 12 / 13 | 1.13e-12 | 0.242 |
| B4 1e-10 | B4 1 | 11 / 12 | 11 / 12 | 2.11e-11 | 3.83 |
| B4 1e-11 | B4 0 | 12 / 13 | 12 / 13 | 1.14e-12 | 0.242 |
| B4 1e-11 | B4 1 | 11 / 12 | 11 / 12 | 2.11e-11 | 3.83 |
| B4 1e-12 | B4 0 | 12 / 13 | 12 / 13 | 1.14e-12 | 0.242 |
| B4 1e-12 | B4 1 | 12 / 13 | 13 / 13 | 2.11e-11 | 3.83 |
| F2 [1e-15, 1e-14, 1e-08] | F2 0 | 10 / 11 | 10 / 11 | 1.12e-10 | 10.4 |
| F2 [1e-15, 1e-14, 1e-09] | F2 0 | 10 / 11 | 10 / 11 | 1.12e-10 | 10.4 |
| F2 [0, 0, 1e-08] | F2 0 | 10 / 11 | 10 / 11 | 1.12e-10 | 10.4 |

Allowance max(1e-12, 1e-10 abs(cold)) on each optimizer coordinate. Both arms report convergence in every probe. Outcome of the open endpoint decision (owner, 2026-10-06): this coordinate allowance is superseded. Round 5 judges F2 against per-coordinate allowances propagated from the accepted observable allowances through the endpoint linearization, where every coordinate uses at most 0.076% at the script and work-matched controls, and judges B4 by its minimax objective, where seeded and cold agree within the 1e-12 absolute floor (at most 3.0% used); B4 coordinate differences are diagnostics of the flat optimum set. The Engine defaults are 1e-12, 1e-12, 1e-10; the round-2 and round-3 F2 calls used 1e-15, 1e-14, 1e-14.

| Controls | Complete call | Cold median (s) | Seeded median (s) | Median paired saving | Paired range | Paired median absolute deviation | Evaluations cold / seeded | Per-evaluation cold / seeded (ms) | Median per-evaluation saving | Outcome |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| F2 parameter tolerance 1e-8 | F2 full | 5.804 | 4.469 | 24.6% | 17.5% to 25.1% | 0.5 percentage points | 21 / 21 | 276.4 / 212.8 | 24.6% | Pass |
| F2 parameter tolerance 1e-8 | B4 full | 7.968 | 7.326 | 7.0% | 2.0% to 15.1% | 1.3 percentage points | 48 / 48 | 166.0 / 152.6 | 7.0% | Pass |
| F2 1e-15, 1e-14, 1e-14 | F2 full | 11.528 | 7.121 | 38.2% | 37.0% to 40.0% | 0.7 percentage points | 40 / 32 | 288.2 / 222.5 | 22.8% | Pass |
| F2 1e-15, 1e-14, 1e-14 | B4 full | 7.851 | 7.514 | 5.5% | 3.5% to 8.9% | 1.3 percentage points | 48 / 48 | 163.6 / 156.5 | 5.5% | Pass |
| Engine defaults | F2 full | 5.974 | 4.984 | 17.5% | 11.7% to 18.9% | 0.9 percentage points | 23 / 24 | 259.7 / 207.7 | 21.0% | Pass |
| Engine defaults | B4 full | 7.958 | 7.484 | 5.9% | 4.8% to 9.3% | 0.5 percentage points | 48 / 48 | 165.8 / 155.9 | 5.9% | Pass |

Seven counterbalanced pairs per row, one numerical thread. An evaluation is one residual or one Jacobian evaluation; a Jacobian evaluation at an already evaluated point reuses its states. The F2 gate is at least 5% median paired saving; the B4 gate (owner amendment, 2026-10-06, replacing at least 5% for B4: no slower beyond noise) is a nonnegative median or a negative median within the paired median absolute deviation. B4 uses its declared Ipopt tolerance in every row.

| Key cold solve | Main iterations / evaluations / Hessians | Seeded iterations / evaluations / Hessians |
|---|---:|---:|
| case-6/figure-2-point-7/cold | 16 / 53 / 17 | 16 / 53 / 17 |
| case-6/figure-2-point-7/accepted | 7 / 14 / 8 | 7 / 14 / 8 |
| case-7/flash-5918.5286-Pa/cold | 43 / 57 / 44 | 43 / 57 / 44 |
| case-7/flash-5918.5286-Pa/accepted | 17 / 22 / 18 | 17 / 22 / 18 |
| khudaida/sequential-8 | 66 / 119 / 74 | 66 / 119 / 74 |
| mea-nominally-unloaded/353.15-K | 13 / 19 / 14 | 13 / 19 / 14 |
