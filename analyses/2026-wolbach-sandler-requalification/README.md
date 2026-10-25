# Wolbach–Sandler size factor: Figiel and MEA requalification

Issue #155 restored the Gross–Sadowski 2002 Eq. 3 size factor
`[sqrt(sigma_i sigma_j)/((sigma_i+sigma_j)/2)]^3` in the unlike association
volumes of four packets. The factor is `0.9940` for water/ethanol and
`0.9923` for water/methanol at 298.15 K, and `0.9972` for water/MEA at 313.15 K. This
leaf measures how much Figiel and MEA evidence moves between the old packet
version (v3, no factor) and the new version (v4, with the factor). The Khudaida
effect is recorded in the [settler reference](../2026-khudaida-idaes-settler/README.md).

## Figiel 2025 mixed-solvent activity coefficients

[`scripts/figiel_mixed_solvent.py`](scripts/figiel_mixed_solvent.py) evaluates
every mixed-solvent Figure 9 literature row at 298.15 K and 1 bar
([rows](results/figiel-figure-09.csv), [summary](results/figiel-figure-09.json)).

| packet family | system | rows | AARD v3 | AARD v4 | max \|γ±(v4)/γ±(v3) − 1\| |
|---|---|---:|---:|---:|---:|
| both | water–ethanol NaCl | 123 | 0.063484 | 0.063485 | 4.3e-6 |
| both | water–ethanol NaBr | 120 | 0.176383 | 0.176382 | 6.8e-6 |
| figiel-yu-held-2025 | water–methanol NaCl | 58 | 0.080928 | 0.080927 | 4.1e-6 |
| figiel-yu-held-2025 | water–methanol NaBr | 108 | 0.074429 | 0.074427 | 5.6e-6 |

The molality-scale mean ionic activity coefficient is referenced to infinite
dilution in the same mixed solvent, so the solvent–solvent cross-association
change largely cancels. Figiel mixed-solvent activity evidence is unchanged at
the reported precision. `figiel-2025-reference-electrolytes` has no
water–methanol cross-association edge in any version, so its methanol rows do
not change.

## MEA–water VLE (Cai 1996)

[`scripts/mea_water.py`](scripts/mea_water.py) solves the bubble temperature
and vapor composition at the 25 admitted Cai rows (66.66 and 101.33 kPa, 2B
MEA) with the retained objective scales of 1 K and 0.01 in vapor water mole
fraction ([rows](results/mea-water-cai.csv), [summary](results/mea-water-cai.json)).

| quantity | packet v3 | packet v4 |
|---|---:|---:|
| mean scaled square at packet `k_MEA,H2O = -0.0735275` | 4.39500 | 4.42649 |
| max \|ΔT_bubble\| vs v3 | – | 0.064 K |
| max \|Δy_water\| vs v3 | – | 0.00099 |
| one-coordinate refit of `k_MEA,H2O` | −0.0735275 | −0.0738817 |
| mean scaled square at refit | 4.39500 | 4.42423 |

The v3 refit reproduces the retained neutral-replay coefficient to `3e-9`,
which validates this current-API replay against the retained fit. The v4 best
coefficient moves by `−3.5e-4`, inside the retained 95% interval
`[-0.07792, -0.06914]`. The topology selection was not rerun: the retained grouped
score of 3B (102.7) is eleven times that of 2B (9.21), and the objective changes by
0.7%.

The [Born-form comparison](../2026-mea-parameter-estimation/results/analysis/ionic-electrostatics/README.md)
uses current packet v5. The v3-to-v4 size-factor update moved fixed-parameter
fugacities by at most 0.53%; v5 additionally converts MEA 2B κ to the Engine's
σ³ convention (`κ = 0.036287`). The converged v4/v5 refits preserve the same
ranking: SSM has the lowest χ² but is rejected on physical grounds for
implausible ion–solvent k_ij and negative packet-coordinate osmotic
coefficients; SSM+DS is weakly favoured over original Born. Original Born
stays the runtime baseline; #161 closed with the Born-form choice left to MEA
(MEA-Thermodynamics #121, decision 21).
The resolved-step refit and the decision evidence are retained in #163 and its
ionic-electrostatics analysis.

Derivative actions for the MEA selection on v4 are requalified in the
[shared derivative qualification](../2026-shared-derivative-qualification/README.md).

## Reproduce

```bash
OMP_NUM_THREADS=1 .venv/bin/python analyses/2026-wolbach-sandler-requalification/scripts/figiel_mixed_solvent.py
OMP_NUM_THREADS=1 .venv/bin/python analyses/2026-wolbach-sandler-requalification/scripts/mea_water.py
```

Claim limits: the stated literature rows and fixed parameters only; no
uncertainty interval for the v4 refit and no reactive MEA evaluation.
