All ionic qualification assertions passed. Each operation has three timing/work samples; the numerical report retains the first sample.

| Declared state | T (K) | P (Pa) | Charged admitted / total | Discovery |
|---|---:|---:|---:|---|
| `case-6/figure-2-point-7` | 293.15 | 100000.0 | 4 / 4 | converged |
| `case-7/bubble` | 318.15 | 7398.160760169676 | 3 / 4 | not requested (free-pressure bubble) |
| `case-7/flash-5918.5286-Pa` | 318.15 | 5918.5286 | 3 / 4 | resolved descent without an admitted lower-G declaration |
| `case-7/liquid-8877.7929-Pa` | 318.15 | 8877.7929 | 3 / 4 | converged |
| `case-6/aqueous-salt-one-liquid` | 298.15 | 100000.0 | 2 / 2 | converged |

Every charged seed matches the approved availability expectation. TPD/RT is dimensionless; its sign near zero is not a stability certificate.

| Declared state | Charged seed | Expected | Returned evidence | TPD/RT | Matches design |
|---|---|---|---|---:|---|
| `case-6/figure-2-point-7` | `positive-neutral-base` | admitted | `no_resolved_descent` | 0 | yes |
| `case-6/figure-2-point-7` | `near-pure-0` | admitted | `no_resolved_descent` | 3.5527136788e-15 | yes |
| `case-6/figure-2-point-7` | `near-pure-1` | admitted | `no_resolved_descent` | 0 | yes |
| `case-6/figure-2-point-7` | `near-pure-2` | admitted | `no_resolved_descent` | 3.77475828373e-15 | yes |
| `case-7/bubble` | `positive-neutral-base` | admitted | `no_resolved_descent` | -5.68434188608e-14 | yes |
| `case-7/bubble` | `near-pure-0` | unresolved | `unresolved` | absent | yes |
| `case-7/bubble` | `near-pure-1` | admitted | `no_resolved_descent` | -8.52651282912e-14 | yes |
| `case-7/bubble` | `near-pure-2` | admitted | `no_resolved_descent` | -5.68434188608e-14 | yes |
| `case-7/flash-5918.5286-Pa` | `positive-neutral-base` | admitted | `no_resolved_descent` | 2.84217094304e-14 | yes |
| `case-7/flash-5918.5286-Pa` | `near-pure-0` | unresolved | `unresolved` | absent | yes |
| `case-7/flash-5918.5286-Pa` | `near-pure-1` | admitted | `no_resolved_descent` | 0 | yes |
| `case-7/flash-5918.5286-Pa` | `near-pure-2` | admitted | `no_resolved_descent` | 5.68434188608e-14 | yes |
| `case-7/liquid-8877.7929-Pa` | `positive-neutral-base` | admitted | `no_resolved_descent` | 0 | yes |
| `case-7/liquid-8877.7929-Pa` | `near-pure-0` | unresolved | `unresolved` | absent | yes |
| `case-7/liquid-8877.7929-Pa` | `near-pure-1` | admitted | `no_resolved_descent` | -1.13686837722e-13 | yes |
| `case-7/liquid-8877.7929-Pa` | `near-pure-2` | admitted | `no_resolved_descent` | -8.52651282912e-14 | yes |
| `case-6/aqueous-salt-one-liquid` | `positive-neutral-base` | admitted | `no_resolved_descent` | 5.3290705182e-15 | yes |
| `case-6/aqueous-salt-one-liquid` | `near-pure-0` | admitted | `no_resolved_descent` | 7.9936057773e-15 | yes |

The five supplied-point checks use native roots and native supplied-trial evaluation.

| Supplied point | T (K) | P (Pa) | Supplied density (mol/m³) | Expected | Returned evidence | TPD/RT | All checks pass |
|---|---:|---:|---:|---|---|---:|---|
| aqueous-liquid | 298.15 | 100000.0 | 55794.7352172 | admitted | `evaluated` | 7.9936057773e-15 | yes |
| aqueous-vapor | 298.15 | 100000.0 | 43.8492327235 | unresolved | `unresolved` | absent | yes |
| salt-rich-liquid | 298.15 | 100000.0 | 37443.5674576 | admitted | `evaluated` | 5.48427036798 | yes |
| salt-rich-half-density | 298.15 | 100000.0 | 18721.7837288 | unresolved | `unresolved` | absent | yes |
| historical-loop-less-mea | 318.15 | 5918.5286 | 2.32146653766 | unresolved | `unresolved` | absent | yes |

The inherited identity and reconstruction checks retain their existing tolerances:

| Quantity | Maximum absolute explicit-trial readback | Criterion |
|---|---:|---:|
| `normalization_residual` | 2.63677968348e-16 | 1e-12 |
| `charge_residual` | 6.93889390391e-18 | 1e-12 |
| `trial_identity_defect_over_rt` | 1.94057704188e-13 | 1e-10 |

The MEA flash discovery stop is retained unchanged: `resolved descent without an admitted lower-G declaration` (#688). The ionic admission assertions pass without asserting successful MEA phase discovery.
