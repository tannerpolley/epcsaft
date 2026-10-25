Wheel SHA-256: `ceb6ab8a48b671fbb48330ff1195b14e2b9372c3829edb3f4a49d268e3895707`.

| Check | Comparisons | Largest error / tolerance | Result |
|---|---:|---:|---|
| F1 | 3167 | 0.0144757 | Pass |
| F1a | 32 | 0.00529064 | Pass |
| F2 | 172 | 0.185207 | Pass |
| F3 | 175 | 0.000116415 | Pass |
| F4 | 7 | 0.607946 | Pass |
| F5 | 12 | 2.83107e-05 | Pass |
| F6 | 131 | 0.0784679 | Pass |

F1 uses 383 rows and 3,167 applicable columns. Its allowance is $2\times10^{-5}\max|J|+\delta$, with $\delta=10^{-7}$ in dimensionless units and $RT\,10^{-7}$ J/mol for energy. Unprescribed difference steps are retained without a verdict.

| Manufactured recovery | Largest error | Required limit |
|---|---:|---:|
| k coordinates | 7.34985e-11 | 1e-06 |
| held diameters (Angstrom) | 1.22524e-11 | 1e-07 |
| log10 residuals (decade) | 1.85207e-11 | 1e-10 |
| C | 2.14939e-12 | 1e-5 |

| Held ion | Stored diameter (Angstrom) | Independent root (Angstrom) | Energy error (J/mol) | Bracket width (Angstrom) |
|---|---:|---:|---:|---:|
| 1 | 2.8973221567948948 | 2.8973221567949281 | 5.2386894822120667e-09 | 7.9580786405131221e-13 |
| 2 | 3.4159107262851882 | 3.4159107262851762 | -1.280568540096283e-09 | 7.9580786405131221e-13 |
| 3 | 4.3232582622092464 | 4.3232582622095492 | 2.1420419216156006e-08 | 7.9580786405131221e-13 |

F3 requires diameter agreement within 1e-6 Angstrom and energy error within 1e-2 J/mol. Every inner hold also meets 1e-6 J/mol and 1e-10 Angstrom.

| Native check group | Maximum | Mean | Outcome |
|---|---:|---:|---|
| fail | 0.03 | 0.02 | fail |
| pass | 0.015 | 0.01 | pass |
| energy | 4 | 3 | pass |
| epsilon | 0.2 | 0.15 | fail |
| baseline | 0.03 | 0.02 | pass |
| unavailable-report | 381931 | 381931 | unavailable |
| unavailable-guard | 381931 | 381931 | unavailable |

The unavailable negative-energy logarithmic residual retains its signed J/mol difference. It blocks a guard; the report-only declaration retains a valid zero objective and passes overall acceptance. The epsilon group deliberately passes its maximum and fails its mean. Metric discrepancies are at most 2.9e-15 in their reported units, against 1e-10.

F6 maximum induced root uncertainty: 7.84679e-12 Angstrom, against 1e-10. The root/map study uses source and manufactured points, three steps and nonunit scales with tied C. Four native failure/cache comparisons require exact agreement.

| Native fit | Weighted objective | Seconds | Covariance | Active coordinates |
|---|---:|---:|---|---|
| F2 | 1.24818e-18 | 9.11792 | available | [] |
| F3 | 0.333078 | 15.5715 | active_bound | [3, 5] |

The native raw-value/column calls took 3.31903 s across F1, and the corresponding central-difference batches took 22.3529 s. These measured calls include construction, reference admission and derivative work. They do not establish an optimized finite-difference baseline or a general speedup.
