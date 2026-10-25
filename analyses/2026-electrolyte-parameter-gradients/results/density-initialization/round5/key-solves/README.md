Key-solve work, main wheel 1fa7ec18 (`build/wheels/main-current`, sha256 7181e0a0…5151) against the #245 round-5 seeded wheel (merged with 1fa7ec18, sha256 02d65247…a680), `tools/time-key-solves --pairs 1`.

| Case | Main iterations / evaluations / Hessians | Seeded iterations / evaluations / Hessians |
|---|---:|---:|
| case-6/figure-2-point-7/cold (ionic LLE) | 16 / 47 / 17 | 16 / 47 / 17 |
| case-6/figure-2-point-7/accepted | 7 / 11 / 8 | 7 / 11 / 8 |
| case-7/flash-5918.5286-Pa/cold (MEA flash) | 43 / 57 / 44 | 43 / 57 / 44 |
| case-7/flash-5918.5286-Pa/accepted | 17 / 22 / 18 | 17 / 22 / 18 |
| khudaida/sequential-8 | 66 / 93 / 74 | 66 / 93 / 74 |
| mea-nominally-unloaded/353.15-K | 13 / 19 / 14 | 13 / 19 / 14 |

Work is identical in all six cases and returned states agree exactly (state difference 0). The single timing pair is not a speed claim; its rows remain at record commit `21d7a112`.
