| Quantity | Count | Maximum relative difference | Largest fraction of allowance | Criterion | Outcome |
|---|---:|---:|---:|---|---|
| density column | 1 | 0 | 0 | max(1e-12, 1e-10 max(abs(seeded), abs(cold))) | Pass |
| density | 5 | 1.71e-15 | 1.71e-05 | max(1e-12, 1e-10 max(abs(seeded), abs(cold))) | Pass |
| slope | 5 | 1.26e-14 | 0.000126 | max(1e-12, 1e-10 max(abs(seeded), abs(cold))) | Pass |
| slope derivative | 5 | 7.02e-10 | 0.000702 | 1e-6 max(1, abs(reference)) | Pass |
| sites | 8 | 6.34e-15 | 6.34e-05 | max(1e-12, 1e-10 max(abs(seeded), abs(cold))) | Pass |
| columns and observables | 28 | 5.19e-13 | 0.00519 | max(1e-12, 1e-10 max(abs(seeded), abs(cold))) | Pass |
| same-point replay chart identity | 39 | 0 | 0 | max(1e-12, 1e-10 abs(cold)) | Pass |
| observables | 204 | 1.08e-10 | 0.011 | max(1e-12 or osmotic floor, 1e-10 abs(cold)) | Pass |
| physical columns | 2496 | 1.37e-12 | 0.0137 | max(1e-12, 1e-10 abs(cold)) | Pass |

The relative difference divides by max(abs(cold), 0.01). The allowance is the row criterion; a fraction at most one passes.

| Fit | Molality (mol/kg) | Rows | Osmotic floor | Maximum absolute difference | Equivalent ln a_w difference | Largest fraction of allowance |
|---|---:|---:|---:|---:|---:|---:|
| B4-full | 0.1 | 6 | 2.78e-10 | 1.55e-12 | 5.6e-15 | 0.0056 |
| B4-full | 1 | 6 | 2.78e-11 | 9.44e-14 | 3.4e-15 | 0.000976 |
| B4-full | 4.5 | 6 | 6.17e-12 | 1.01e-13 | 1.63e-14 | 0.00106 |
| F2-full | 0.001 | 3 | 2.78e-08 | 9.13e-11 | 3.29e-15 | 0.00329 |
| F2-full | 0.002 | 3 | 1.39e-08 | 1.06e-10 | 7.63e-15 | 0.00763 |
| F2-full | 0.005 | 3 | 5.55e-09 | 6.13e-11 | 1.1e-14 | 0.011 |
| F2-full | 0.01 | 3 | 2.78e-09 | 6e-12 | 2.16e-15 | 0.00216 |
| F2-full | 0.02 | 3 | 1.39e-09 | 2.68e-12 | 1.93e-15 | 0.00193 |
| F2-full | 0.05 | 3 | 5.55e-10 | 1.97e-12 | 3.55e-15 | 0.00355 |
| F2-full | 0.1 | 3 | 2.78e-10 | 9.57e-13 | 3.45e-15 | 0.00345 |
| F2-full | 0.2 | 3 | 1.39e-10 | 4.12e-13 | 2.97e-15 | 0.00297 |
| F2-full | 0.3 | 3 | 9.25e-11 | 2.11e-13 | 2.28e-15 | 0.00228 |
| F2-full | 0.4 | 3 | 6.94e-11 | 1.37e-13 | 1.97e-15 | 0.00147 |
| F2-full | 0.5 | 3 | 5.55e-11 | 2.14e-13 | 3.86e-15 | 0.00221 |
| F2-full | 0.6 | 3 | 4.63e-11 | 2.57e-13 | 5.55e-15 | 0.00261 |
| F2-full | 0.7 | 3 | 3.96e-11 | 2.08e-13 | 5.25e-15 | 0.00209 |
| F2-full | 0.8 | 3 | 3.47e-11 | 5.91e-14 | 1.7e-15 | 0.000587 |
| F2-full | 0.9 | 3 | 3.08e-11 | 1.35e-13 | 4.39e-15 | 0.00133 |
| F2-full | 1 | 3 | 2.78e-11 | 1.69e-13 | 6.07e-15 | 0.00164 |
| F2-full | 1.2 | 3 | 2.31e-11 | 2.42e-13 | 1.05e-14 | 0.0023 |
| F2-full | 1.4 | 3 | 1.98e-11 | 4.2e-14 | 2.12e-15 | 0.00039 |
| F2-full | 1.6 | 3 | 1.73e-11 | 9.74e-14 | 5.61e-15 | 0.000989 |
| F2-full | 1.8 | 3 | 1.54e-11 | 9.3e-14 | 6.03e-15 | 0.000937 |
| F2-full | 2 | 3 | 1.39e-11 | 3.13e-14 | 2.26e-15 | 0.000313 |
| F2-full | 2.5 | 3 | 1.11e-11 | 3.75e-14 | 3.38e-15 | 0.000414 |
| F2-full | 3 | 3 | 9.25e-12 | 5e-14 | 5.4e-15 | 0.000554 |
| F2-full | 3.5 | 3 | 7.93e-12 | 6.84e-14 | 8.62e-15 | 0.000646 |
| F2-full | 4 | 3 | 6.94e-12 | 5.44e-14 | 7.84e-15 | 0.000607 |
| F2-full | 4.5 | 3 | 6.17e-12 | 9.31e-14 | 1.51e-14 | 0.00104 |

| Fit start | Residual / Jacobian evaluations, cold | Residual / Jacobian evaluations, seeded | Other work counts equal | Largest optimizer-endpoint difference / allowance | Outcome |
|---|---:|---:|---|---:|---|
| F2 0 | 25 / 15 | 20 / 12 | Yes | 8.86 | Fail |
| B4 0 | 12 / 13 | 12 / 13 | Yes | 0.242 | Pass |
| B4 1 | 11 / 12 | 11 / 12 | Yes | 3.83 | Fail |

Optimizer endpoints compare under max(1e-12, 1e-10 abs(cold)). 4 coordinates exceed it: F2 start 0 optimizer2, F2 start 0 optimizer4, F2 start 0 optimizer5, B4 start 1 optimizer5. The F2 fit is a manufactured recovery: the largest seeded distance from the manufactured coordinates is 8.2e-11, against 1.1e-11 for the cold fit, so the seeded F2 fit stopped earlier at a less refined point. Both fits report convergence under the same controls. The same-point replay chart identity holds by construction and is not evidence of an unchanged optimum.

| Complete call, base 5fec07f8 | Cold median (s) | Seeded median (s) | Median paired saving | Paired range | Paired median absolute deviation | Pairs | Outcome |
|---|---:|---:|---:|---:|---:|---:|---|
| F2 full | 13.654 | 8.823 | 40.2% | 24.5% to 48.1% | 7.9 percentage points | 7 | Pass |
| B4 full | 9.616 | 9.151 | 9.8% | -7.6% to 18.2% | 5.6 percentage points | 7 | Pass |

The paired times were measured in round 2 at base 5fec07f8 and were not rerun on the merged wheel. Round 2 checked equal dimensions, outcomes and trial failures in each pair but not evaluation counts; the round-3 endpoint comparison shows unequal F2 evaluation counts, so the F2 saving is not established for equivalent accepted work. F2 requires at least 5% median paired saving; B4 requires a nonnegative median or a negative median within the paired median absolute deviation.


| Key cold solve | Main iterations / evaluations / Hessians | Seeded iterations / evaluations / Hessians | Maximum relative state difference |
|---|---:|---:|---:|
| case-6/figure-2-point-7/cold | 16 / 53 / 17 | 16 / 53 / 17 | 0.0 |
| case-7/flash-5918.5286-Pa/cold | 43 / 57 / 44 | 43 / 57 / 44 | 0.0 |

The manufactured-recovery variance factor is 3.47653699266e-27, versus 3.47493504088e-27; the absolute difference 1.6e-30 is within 6.67e-21, twice the final-cost criterion divided by the three residual degrees of freedom.

| Invalidation or fallback | Native verdict | Anchor cause |
|---|---|---|
| Missing branch witness after a parameter change | Pass | ineligible |
| Liquid density offered for a vapor request | Pass | ineligible |
| Irregular retained closure | Pass | ineligible |
| Density outside the anchor window | Pass | left_domain |
| Empty cold domain | Pass | density domain is empty |

Missing-witness and irregular offers also produce no eligible or accepted seed event.
