The saving is unresolved: the strict empirical envelope criterion is not met.

The six paired whole-call observations and same-image controls give:

| R | C | −ln R | Saving resolved | Median paired saving / µs |
|---|---|---|---|---|
| 0.874773 | 0.155749 | 0.133791 | False | 5.66 |

| Arm | Median / µs | Min block median / µs | Max block median / µs |
|---|---|---|---|
| A | 47.52 | 44.21 | 54.2 |
| B | 40.71 | 39.07 | 50.199 |
| FeOs | 26.81 | 26 | 31.39 |
| Clapeyron | 27.689 | 27.02 | 30.44 |

| Fresh comparison | Pairs | Median ratio | Full paired range |
|---|---|---|---|
| A-control/A | 3 | 0.990589 | 0.909225–1.1097 |
| B-control/B | 3 | 1.04351 | 1.01625–1.16853 |
| B/A | 6 | 0.874773 | 0.792601–1.01865 |
| B/Clapeyron | 6 | 1.45635 | 1.28351–1.85785 |
| B/FeOs | 6 | 1.52958 | 1.36856–1.5992 |

| Native returned quantity | Value | Unit |
|---|---|---|
| Saturation pressure | 998660.8957399128 | Pa |
| Caller full-row criterion | 1e-12 | 1 |
| Full raw residual maximum | 4.440892098500626e-15 | 1 |
| liquid density | 11100.25123320961 | mol/m³ |
| vapor density | 482.5121268150347 | mol/m³ |

All 10040 retained values from 251 regular states agree bit for bit; availability differences: 0.

| Gate | Retained observation |
|---|---|
| Native actions and fitting columns | 985 rows; maximum A/B scaled difference zero |
| Equilibrium corpora | Zero differing values; baseline admission and refusal inventories preserved |
| Explicit native cases | 19 passed, zero failures |
| Key solves and frozen work | All work counts identical; maximum state difference zero |
| F1 reference work | 295 computations; 3388 columns; 1223 density solves; 1426 equilibrium iterations |
| Reacting heat-capacity work | 3 computations; 280 columns; 11 density solves; 13 equilibrium iterations |
| F1 serial CPU / s | A 4.880778; B 4.749343; ceiling 20.349267 |
| Reacting heat-capacity serial CPU / s | A 0.633522; B 0.647181; ceiling 2.473741 |
| Full fit wall / s | A 7.858261; B 7.909250; ceiling 120 |
| Full fit critical-path CPU / s | A 6.015962; B 6.033381 |
| Same-host full-fit control / s | 9.161995–9.269454 |
| Affected and full checks | Passed |
| Native lint and sanitizers | Passed for the implementation; frozen-head records retained with PR evidence |
