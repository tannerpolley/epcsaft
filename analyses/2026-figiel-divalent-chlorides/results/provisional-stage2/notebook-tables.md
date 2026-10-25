| Temperature (K) | Pressure (Pa) | Molality interval (mol/kg water) | Basis |
|---|---|---|---|
| 298.15 | 100000.0 | 0.001–4.5 | binary MgCl2/CaCl2, mol formula units/kg water, on unchanged aqueous LiCl/NaCl/KCl /1 basis |

| Salt | γ± maximum (decade) | γ± mean absolute (decade) | Φ maximum (decade) | Φ mean absolute (decade) |
|---|---|---|---|---|
| MgCl2 | 0.06407599491 | 0.02947203674 | 0.03957717574 | 0.01578270588 |
| CaCl2 | 0.04834780429 | 0.02311442551 | 0.0334417264 | 0.01328778411 |

The original every-row limit is 0.02 decade. Of 168 targets, 73 exceed it. The owner relaxed that limit on 2026-10-06; the use basis is **provisional: relaxed target**. The largest stored-row error change from the original stage-2 fit is 5.274990405e-11 decade, below the unchanged native minimax tolerance 1e-10.

| Quantity | MgCl₂ | CaCl₂ |
|---|---|---|
| sigma_angstrom | 2.717877321 | 2.839274308 |
| packing_angstrom | 2.391732043 | 2.498561391 |
| DH_angstrom | 2.391732043 | 2.498561391 |
| k_water | -1 | -1 |
| k_chloride | 1 | 1 |
| born_angstrom | 3.15860961 | 3.852587089 |

| Native search | MgCl₂ | CaCl₂ |
|---|---|---|
| Completed Ceres starts | 28 | 28 |
| Usable Ceres starts | 28 | 27 |
| Completed minimax starts | 1 | 1 |

| Numerical quantity | MgCl₂ | CaCl₂ | Requirement |
|---|---|---|---|
| Hydration error (J/mol) | -4.563480616e-08 | -9.54605639e-09 | absolute ≤ 10 |
| Born root correction (Å) | 8.491761827e-14 | 2.642643717e-14 | absolute ≤ 1e-10 |
| Maximum pressure error (Pa) | 3.659851791e-05 | 4.045838432e-05 | absolute < 0.01 |
| Minimum pressure-density slope (Pa m³/mol) | 81858.11103 | 81856.78005 | positive |
| Maximum charge-weighted mole fraction | 0 | 0 | zero |
| Verified held sigma columns | 84 | 84 | all retained rows |
| Minimax stationarity diagnostic | 6.122879981e-13 | 8.247685868e-11 | native stationary point |
| Minimax complementarity diagnostic | 1.025551376e-11 | 3.810741308e-11 | native diagnostic |

All 948 stage-1 predictions are preserved against the current-wheel reference; maximum absolute change is 0. All 45 enforced groups pass among 57 evaluated groups; G1 bromides, G2 solvent transfer and G5 water keep their native source-baseline rules.

| Guard | Metric | Stored | Source baseline | Threshold | Native result |
|---|---|---|---|---|---|
| G1 LiBr gamma | maximum | 0.05099144387 | 0.05099144388 | 0.050991443883889286 | pass |
| G1 LiBr gamma | mean | 0.02134337844 | 0.02170102655 | 0.021701026552609484 | pass |
| G1 LiBr osmotic | maximum | 0.01992325786 | 0.04238722139 | 0.04238722138799963 | pass |
| G1 LiBr osmotic | mean | 0.009646301158 | 0.01104945735 | 0.011049457352045015 | pass |
| G1 NaBr gamma | maximum | 0.02504297838 | 0.06656739032 | 0.066567390319515 | pass |
| G1 NaBr gamma | mean | 0.007007315396 | 0.03229097607 | 0.032290976068764105 | pass |
| G1 NaBr osmotic | maximum | 0.02306895657 | 0.02425387407 | 0.02425387406557926 | pass |
| G1 NaBr osmotic | mean | 0.005303186725 | 0.01177496873 | 0.01177496872822146 | pass |
| G1 KBr gamma | maximum | 0.02487646675 | 0.02614813663 | 0.02614813662600866 | pass |
| G1 KBr gamma | mean | 0.00616311168 | 0.01278146558 | 0.012781465576804976 | pass |
| G1 KBr osmotic | maximum | 0.02546399595 | 0.03511858129 | 0.03511858129021413 | pass |
| G1 KBr osmotic | mean | 0.004811163284 | 0.007790953483 | 0.007790953483218645 | pass |
| G2 Na+ water-ethanol | maximum | 4443.341548 | 4993.835888 | 4993.835888128529 | pass |
| G2 Na+ water-ethanol | mean | 2662.645236 | 3170.000691 | 3170.0006911293767 | pass |
| G2 K+ water-methanol | maximum | 1657.576423 | 1803.19083 | 1803.190830452876 | pass |
| G2 K+ water-methanol | mean | 469.748099 | 1034.487199 | 1034.487199469702 | pass |
| G5 water | maximum | 6.757713957 | 6.757713957 | 6.757713957314948 | pass |
| G5 water | mean | 1.828865076 | 2.922470749 | 2.9224707488544044 | pass |

| Re-baseline evidence | Value |
|---|---|
| Reference rows | 948 |
| Original-reference rows beyond their preservation tolerance | 40 |
| Maximum original-reference change | 2.848403735e-10 |
| Pressure closures improved | 40 |
| Individual Φ identity errors improved/worsened | 20/20 |
| Old/new maximum absolute Φ Gibbs–Duhem inconsistency | 2.990925285e-10/1.020306062e-10 |
| Preservation rule | abs(candidate-reference) <= 1e-12*abs(reference) for nonzero reference, otherwise <=1e-12. |

Source commit: `71c8a69c3feefd01c0b1aa77b00e9418180cca29`. Merged main: `ec0bafa71b5022821f1774b2b9fffb06622be87c`.

Wheel SHA-256: `fdf3fbb51244eb4991f0743a6acd377dc0296d5909b77e896d050c6e36fdf57c`.

Packet fingerprint: `d79865152abaf527b4602c13b20e3e2c46db8eb0b3f2b294b138b8bdde65a865`.

Unchanged stage-1 parameter SHA-256: `cec1bfc9c833f7196c7bc11d51d21aa8bc224cc388b04f165712acc5c4473626`.

| Check | Result |
|---|---|
| packet_validation | 0 |
| check_affected | 0 |
| check_all | 0 |
| render | 0 |
| notebook_numbers | 0 |
