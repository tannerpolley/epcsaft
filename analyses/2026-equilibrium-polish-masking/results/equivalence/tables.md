### Four checks against the equivalence standard

| Comparison | Main | Candidate | Former absolute limit | New absolute limit | Scaled difference | Passed |
| --- | --- | --- | --- | --- | --- | --- |
| trace/1e-08/10/kkt_residual | 1.4210854715202004e-14 | 8.000711204658728e-12 | 1e-12 | 4.9896302276469585e-09 | 1.600619682334627e-13 | True |
| solved/thermal/3 | -126.93949786517828 | -126.93949786518097 | 1.2693949786518097e-10 | 1.269394978651783e-08 | 2.115851713881684e-14 | True |
| e5-7/thermal-320/phase[0].property/TotalIsobaricHeatCapacity/A1 | -0.03905144680598874 | -0.03905144680598565 | 1e-12 | 3.905144680598874e-12 | 7.907025321185247e-14 | True |
| r4/recovery/variance_factor | 3.4749350408803575e-27 | 4.600761778464367e-27 | 4.600761778464367e-39 | 1e-24 | 1.1258267375840098e-13 | True |

### Near-zero comparisons and stationarity basis

| Quantity | Value |
| --- | --- |
| Relative equivalence limit | 1e-10 |
| Extensive objective-gradient maximum | 49.896302276469584 |
| Candidate stationarity per mole of the 10 mol feed | 8.000711204658728e-13 |
| Variance absolute floor: squared normalized-residual floor | 1e-24 |
| Variance raw relative difference (absolute floor controls) | 0.3239849736295465 |

### Maximum returned-state preservation differences

| Budget | Maximum | Limit |
| --- | --- | --- |
| observable-relative | 1.090752974497604e-14 | 1e-11 |
| trace-relative | 0.0 | 1e-11 |
| phase-relative | 1.5758037583014238e-14 | 1e-11 |

### Work and correction stop

| Case | Main | Candidate | Stop fired |
| --- | --- | --- | --- |
| masking | 6/11/7 | 6/11/7 | False |
| case-6/figure-2-point-7 | 16/53/17 | 16/47/17 | True |
| case-7/flash-5918.5286-Pa | 43/57/44 | 43/57/44 | False |
| mea-nominally-unloaded | 14/20/15 | 14/20/15 | False |
| free-T | 7/13/8 | 7/9/8 | True |
| free-P | 7/13/8 | 7/10/8 | True |
| small-positive-phase | 4/9/5 | 4/9/5 | False |
| khudaida/point-1 | 11/16/12 | 11/14/12 | True |
| khudaida/point-2 | 10/17/11 | 10/13/11 | True |
| khudaida/point-3 | 7/13/8 | 7/11/8 | True |
| khudaida/point-4 | 8/16/9 | 8/12/9 | True |
| khudaida/point-5 | 8/17/9 | 8/12/9 | True |
| khudaida/point-6 | 8/15/9 | 8/11/9 | True |
| khudaida/point-7 | 7/12/8 | 7/10/8 | True |
| khudaida/point-8 | 7/13/8 | 7/10/8 | True |

### Sequential work

| Case | Main | Candidate | Solved points |
| --- | --- | --- | --- |
| Khudaida | 66/119/74 | 66/93/74 | 8 |

### Small-phase and unloaded-MEA safeguards

| Quantity | Value |
| --- | --- |
| Small-phase amount correction / phase amount | 5.551115148669801e-09 |
| Same amount correction / feed | 1.110223024625157e-16 |
| Correction-stop threshold | 1e-12 |
| Topology phase-amount floor | 1e-08 |
| Unloaded fourth accepted undamped correction | 4.931909477903928e-06 |

### Paired solve timing

| Case | Main median (s) | Candidate median (s) | Median change (%) | Range (%) | Median absolute deviation (%) |
| --- | --- | --- | --- | --- | --- |
| case-6/figure-2-point-7/cold | 0.067785 | 0.0661269 | -7.891 | [-24.502, +3.961] | 5.082 |
| case-6/figure-2-point-7/accepted | 0.028535 | 0.0291375 | -6.788 | [-14.779, +119.673] | 6.708 |
| case-7/flash-5918.5286-Pa/cold | 0.112514 | 0.104426 | -6.754 | [-12.254, +35.713] | 5.324 |
| case-7/flash-5918.5286-Pa/accepted | 0.049778 | 0.0481629 | -3.833 | [-19.027, +59.765] | 4.662 |
| khudaida/sequential-8 | 0.293468 | 0.279005 | -1.761 | [-14.674, +21.562] | 3.227 |
| mea-nominally-unloaded/353.15-K | 0.0808708 | 0.07558 | +0.359 | [-7.987, +12.005] | 5.866 |

### Timing controls

| Quantity | Value |
| --- | --- |
| pairs | 7 |
| cpu | 0 |
| nice | 10 |
| threads | 1 |
