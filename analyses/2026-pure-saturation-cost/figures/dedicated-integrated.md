The saving is unresolved: the strict empirical envelope criterion is not met.

The six paired whole-call observations and same-image controls give:

| R | C | −ln R | Saving resolved | Median paired saving / µs |
|---|---|---|---|---|
| 0.877719 | 0.149696 | 0.130428 | False | 5.005 |

| Arm | Median / µs | Min block median / µs | Max block median / µs |
|---|---|---|---|
| A | 40.71 | 40.13 | 42.21 |
| B | 36.805 | 34.91 | 43.3 |
| FeOs | 31.17 | 24.9 | 31.8 |
| Clapeyron | 26.93 | 26.51 | 34.55 |

| Fresh comparison | Pairs | Median ratio | Full paired range |
|---|---|---|---|
| A-control/A | 3 | 1.01199 | 1.00822–1.01688 |
| B-control/B | 3 | 0.985703 | 0.86097–1.01461 |
| B/A | 6 | 0.877719 | 0.864746–1.07899 |
| B/Clapeyron | 6 | 1.30579 | 1.0576–1.63335 |
| B/FeOs | 6 | 1.26698 | 1.11999–1.48876 |

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
| Native actions and fitting columns | 985 rows; maximum scaled difference zero |
| Regular states | 251 states; 10040 bitwise equal values; no availability differences |
| Equilibrium corpora | Zero differing values in all six inventories |
| Explicit native cases | 19 passed; zero failed |
| Key solves and frozen work | Identical states and work; no gate failures |
| F1 serial CPU / s | A 8.636359; B 5.915056; ceiling 20.349267 |
| Reacting heat-capacity serial CPU / s | A 0.729559; B 0.701602; ceiling 2.473741 |
| Full fit wall / s | A 10.027686; B 10.557869; ceiling 120 |
| Same-host full-fit control / s | [9.075746299989987, 9.558330710977316] |
| Affected and full checks | Passed; final delivery logs retained |
| Native lint and sanitizers | Passed on the frozen integrated Engine image; independent delivered review pending |
