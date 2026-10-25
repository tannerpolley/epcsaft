| Controls | Coordinate | Seeded − cold | Allowance | Fraction used | Allowance without absolute floors | Osmotic-row share of allowance |
|---|---:|---:|---:|---:|---:|---:|
| Round-2 script controls (1e-15, 1e-14, 1e-14) | 0 | -3.15e-11 | 6.1e-08 | 0.052% | 4.79e-08 | 53.8% |
| Round-2 script controls (1e-15, 1e-14, 1e-14) | 1 | -6.58e-11 | 1.27e-07 | 0.052% | 9.91e-08 | 53.9% |
| Round-2 script controls (1e-15, 1e-14, 1e-14) | 2 | -3.97e-11 | 8.89e-08 | 0.045% | 6.95e-08 | 54.7% |
| Round-2 script controls (1e-15, 1e-14, 1e-14) | 3 | -8.55e-11 | 1.9e-07 | 0.045% | 1.48e-07 | 54.7% |
| Round-2 script controls (1e-15, 1e-14, 1e-14) | 4 | -3.05e-11 | 6.88e-08 | 0.044% | 5.42e-08 | 54.9% |
| Round-2 script controls (1e-15, 1e-14, 1e-14) | 5 | -8.86e-11 | 1.99e-07 | 0.045% | 1.57e-07 | 55.0% |
| Round-2 script controls (1e-15, 1e-14, 1e-14) | 6 | 1.68e-11 | 3.22e-08 | 0.052% | 2.42e-08 | 54.3% |
| Work-matched (1e-15, 1e-14, 1e-8) | 0 | -4.14e-11 | 6.1e-08 | 0.068% | 4.79e-08 | 53.8% |
| Work-matched (1e-15, 1e-14, 1e-8) | 1 | -8.67e-11 | 1.27e-07 | 0.069% | 9.91e-08 | 53.9% |
| Work-matched (1e-15, 1e-14, 1e-8) | 2 | -5.18e-11 | 8.89e-08 | 0.058% | 6.95e-08 | 54.7% |
| Work-matched (1e-15, 1e-14, 1e-8) | 3 | -1.12e-10 | 1.9e-07 | 0.059% | 1.48e-07 | 54.7% |
| Work-matched (1e-15, 1e-14, 1e-8) | 4 | -3.58e-11 | 6.88e-08 | 0.052% | 5.42e-08 | 54.9% |
| Work-matched (1e-15, 1e-14, 1e-8) | 5 | -1.04e-10 | 1.99e-07 | 0.052% | 1.57e-07 | 55.0% |
| Work-matched (1e-15, 1e-14, 1e-8) | 6 | 2.43e-11 | 3.22e-08 | 0.076% | 2.42e-08 | 54.3% |

Every F2 coordinate passes at both control sets; the largest fraction used is 0.076%. Endpoint Jacobian condition number 1027; the returned factors reproduce the weighted Jacobian to 1.0e-15 of its largest singular value. The allowance at the seeded endpoint differs from the cold one by at most 9.5e-11 relative.

| Controls | B4 start | Seeded objective | Cold objective | Relative difference | Allowance max(1e-12, 1e-10·abs(cold)) | Fraction used | Outcome |
|---|---:|---:|---:|---:|---:|---:|---|
| Round-2 script controls (1e-15, 1e-14, 1e-14) | 0 | 0.00454722134718699 | 0.00454722134719985 | -2.83e-12 | 1e-12 | 1.3% | Pass |
| Round-2 script controls (1e-15, 1e-14, 1e-14) | 1 | 0.00454722134718342 | 0.00454722134721388 | -6.70e-12 | 1e-12 | 3.0% | Pass |
| Work-matched (1e-15, 1e-14, 1e-8) | 0 | 0.00454722134718699 | 0.00454722134719985 | -2.83e-12 | 1e-12 | 1.3% | Pass |
| Work-matched (1e-15, 1e-14, 1e-8) | 1 | 0.00454722134718342 | 0.00454722134721388 | -6.70e-12 | 1e-12 | 3.0% | Pass |

B4 controls are its declared Ipopt tolerance in both sets, so both sets repeat the same calls.

| B4 start | Largest coordinate difference, seeded − cold (diagnostic) | Coordinate |
|---:|---:|---:|
| 0 | -1.13e-12 | 1 |
| 1 | -2.11e-11 | 1 |

| Key solve | Main iterations / evaluations / Hessians | Seeded iterations / evaluations / Hessians |
|---|---:|---:|
| case-6/figure-2-point-7/cold | 16 / 47 / 17 | 16 / 47 / 17 |
| case-6/figure-2-point-7/accepted | 7 / 11 / 8 | 7 / 11 / 8 |
| case-7/flash-5918.5286-Pa/cold | 43 / 57 / 44 | 43 / 57 / 44 |
| case-7/flash-5918.5286-Pa/accepted | 17 / 22 / 18 | 17 / 22 / 18 |
| khudaida/sequential-8 | 66 / 93 / 74 | 66 / 93 / 74 |
| mea-nominally-unloaded/353.15-K | 13 / 19 / 14 | 13 / 19 / 14 |
