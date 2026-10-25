B3/B4 replay wheel SHA-256: `4e331275a44ab795d0daa18f09455acfaeaa7681372d2fd07f8e809655008bfd`.

| Analytic manufactured comparison | Largest absolute error | Limit |
|---|---:|---:|
| B1 coordinates | 6.0003391e-11 | 1e-07 |
| B1/B2 epigraph | 1.000023e-12 | 1e-09 |
| B1 mean slacks | 2.0000147e-13 | 1e-09 |
| B1 multipliers | 8.0971563e-11 | 1e-06 |

B1 uses x=0.75, z=0.5, d=1.5 and t=0.0175 decade, with both starts and the rescaled repeat. Its equality multiplier is 0.005 decade per H-unit. B2 recovers x=1.5, d=3 and t=0.015 decade; the analytic least-squares point x=0.6 has a 0.024-decade maximum, above 0.02.

| Physical-hold manufactured solve | Start 0 | Start 1 | Required limit |
|---|---:|---:|---:|
| t (decade) | 0.00454722134722 | 0.0045472213472 | 0.00500001 |
| epigraph gap (decade) | 2.03862968173e-13 | 1.88610412011e-13 | 1e-08 |
| scaled primal | 0 | 0 | 1e-08 |
| scaled stationarity | 3.89609289796e-12 | 4.45250449128e-12 | 1e-08 |
| scaled dual-sign defect | 0 | 0 | 1e-08 |
| scaled complementarity | 1.0198778726e-11 | 1.08064602358e-11 | 1e-08 |
| permittivity maximum (1) | 0.099999999875 | 0.099999999875 | 0.1 |
| permittivity mean (1) | 0.0999999995262 | 0.0999999995262 | 0.1 |
| stored hydration maximum (J/mol) | 1.22818164527e-08 | 1.7462298274e-10 | 10 |

The two native solves used 12 and 11 iterations, with a measured batch duration of 1.06122 s (300/120 s per fit; 600/240 s shared).
Maximum estimated scaled root error: 5.5042979e-13, against 1e-9. Maximum reconstructed full-space stationarity component: 8.905009e-14, against 1e-8. The held components vanish by the multiplier definition; B1's analytic KKT conditions independently qualify that multiplier.
Requiring 1e-7 J/mol and 1e-11 Angstrom roots changes each start's 45 scaled rows and 450 columns by 0, against 1e-9. The original roots already satisfy the tighter criteria, so this replay returns the same iterate and establishes no additional error bound.

| Returned parameter | Start 0 | Start 1 | Unit |
|---|---:|---:|---|
| Li-water k | -0.579152261 | -0.533913505 | 1 |
| Li-Cl k | 0.460384168 | 0.546244294 | 1 |
| Na-water k | -0.999999999 | -0.999999999 | 1 |
| Li Born diameter | 2.82593943 | 2.81519909 | Angstrom |
| K Born diameter | 4.1481693 | 4.1434935 | Angstrom |

The two starts return distinct lithium and potassium parameters with nearly identical objective values. They supply no unique parameter estimate; Na-water k is at its lower bound in both.

B3 has 4656 retained numeric comparisons at source, manufactured truth and both returned points. Its criterion remains 1e-12 relative, with 1e-12 absolute only at zero. At each fitted point both prediction evaluators replay with common cold density initialization; stored-fit minus cold predictions are retained separately as initialization diagnostics. All five kinds, nonunit affine scales, tied C, frozen thresholds and stored checks are covered. The report-only transfer check fails while the enforced rows remain identical and the fitted point passes. Zero-Ceres-weight permittivity rows enter B's constraints. B3 checks consistency across selections of the shared evaluator; independent derivative evidence remains slice A's F1 on its identified wheel.

These synthetic rows verify local solver integration and numerical stationarity, not physical validation, #209 adoption, uniqueness or global optimality. The admitted finite reference endpoint retains its existing limits.


| Row-3 density sensitivity | Value | Unit |
|---|---:|---|
| Exact density derivative | 0.0003046811946 | m³/mol |
| One density ulp | 7.275957614e-12 | mol/m³ |
| Smooth derivative times one ulp | 2.216847457e-15 | 1 |
| Observed absolute change | 9.688916336e-13 | 1 |
| Observed / smooth change | 437.0583237 | 1 |