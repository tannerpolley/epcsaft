| Numerical comparison | Count | Maximum error / criterion | Outcome |
|---|---:|---:|---|
| issue-245 | 39 | 0.254237637 | Pass for these comparisons |
| r2-minimax | 34 | 0 | Pass for these comparisons |
| F2-one_step | 2197 | 0 | Pass for these comparisons |
| B4-one_step | 23 | 1.1936855 | Fail |
| m5-pressure-slope | 4 | 0.000288832568 | Pass for these comparisons |

The failed dimensionless osmotic coefficient is 0.936125529328445, compared with cold 0.936125529329562. Its error 1.11743947e-12 exceeds 9.36125529e-13. All 13 physical-chart entries agree exactly. The traced capped call contains 17 accepted fitting-density anchors.

The full qualification is incomplete. These passing individual comparisons do not close the missing intermediate replay, invalidation or cause-preserving fallback cases.

| Complete-call boundary | Cold median (s) | Seeded median (s) | Median reduction | Median paired saving | Pair range | Pairs | Timing criterion |
|---|---:|---:|---:|---:|---:|---:|---|
| F2 one_step | 2.273542 | 2.263418 | 0.445% | -0.444% | -94.092% to 19.512% | 7 | Fail |
| F2 full | 14.001041 | 8.423364 | 39.838% | 37.450% | 25.669% to 46.533% | 7 | Pass for timing only |
| B4 one_step | 2.077041 | 2.112389 | -1.702% | -0.671% | -7.632% to 5.785% | 7 | Fail |
| B4 full | 10.243999 | 10.776452 | -5.198% | 4.538% | -13.993% to 23.625% | 7 | Fail |

Positive reduction means a shorter seeded call. F2 uses alternating AB/BA pairs starting with AB; B4 starts with BA. All timed pairs retain identical returned dimensions and solver outcomes. The equality failure prevents a claim of equivalent accepted work. The F2 one-step trace contains no accepted fitting-density anchor; every recorded chart gate rejects an offer.

| Complete-call boundary | AB median saving | BA median saving | Paired median absolute deviation |
|---|---:|---:|---:|
| F2 one_step | 0.000% | -0.465% | 3.556 percentage points |
| F2 full | 38.975% | 37.450% | 4.023 percentage points |
| B4 one_step | -0.671% | -2.457% | 5.706 percentage points |
| B4 full | 0.042% | 5.889% | 4.496 percentage points |

The full installed-wheel check passes the native suites, architecture, packet, type and lint checks, but fails the native/public manufactured-recovery variance-factor comparison. Its operands and unchanged relative criterion are retained with the other comparisons. The affected check against the moving main reference stops on pre-existing numeric literals in the unrelated reaction-fitting notebook. The selected density-initialization notebook has no flagged numeric literals.

| Key cold solve | A iterations / evaluations / Hessians | B iterations / evaluations / Hessians | Maximum relative state difference |
|---|---:|---:|---:|
| case-6/figure-2-point-7/cold | 16 / 53 / 17 | 16 / 53 / 17 | 0.0 |
| case-7/flash-5918.5286-Pa/cold | 43 / 57 / 44 | 43 / 57 / 44 | 0.0 |
