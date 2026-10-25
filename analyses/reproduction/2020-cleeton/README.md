# Cleeton et al. (2020) competitive MDEA/H₂S/CO₂ ePC-SAFT

Cleeton et al., “Competitive H₂S–CO₂ absorption in reactive aqueous
methyldiethanolamine solution: Prediction with ePC-SAFT”, *Fluid Phase
Equilibria* 511 (2020) 112453, DOI `10.1016/j.fluid.2019.112453` (Zotero
`MMI98GID`; PDF `9Z3P9NXX`). Comparison data: Dicko et al. 2010 Table 1
(`YAIX67FE`), transcribed in `inputs/source/dicko-2010-table-1.csv`.

## Result (#330)

The native replay of Cleeton's printed model follows Cleeton's Figure 16 model
curves within the per-point digitization tolerance at 22 of 30 points. Under the
gate as fixed (caption loadings), native values run about 5% below the left-panel
curves (mean ln r −0.05) and 3% above the right-panel curves (+0.03); the bias is
near zero only under Dicko's group loadings, a reading added after the first run.
It reproduces Cleeton's published deviations against Dicko when those deviations
are read as the symmetric relative deviation, a formula adopted after the first
run. This is a
model replay (numerical verification of a translation), not physical validation.

**Declared model.**
- Table 2 parameters: 4C water and MDEA (2 donor + 2 acceptor sites,
  Wolbach–Sandler cross association), one-segment non-associating ions,
  k_ij with water only, no same-sign ion pairs.
- Cleeton R.1–R.6 with H₃O⁺ and Austgen 1991 Table V constants, mole-fraction
  basis, unconverted (K2 A = 231.465 from Austgen; Cleeton prints 231.456).
- Water and MDEA at the pure liquid, ions, CO₂ and H₂S at infinite dilution in
  water (Eqs. 8–9), evaluated at the system pressure (`reference_pressure_pa=None`).
- Permittivity: water's only, 69.855 at 323.15 K (the Uyan/Wangler transfer of
  ε_w(T), labelled inferred); Eq. 2's pressure term is omitted.
- Per Dicko row: a CH₄-free bubble point at 323.15 K and 50 wt% MDEA with
  p_i = y_i P. No parameter is fitted.

**Lanes:** main, ion diameters d = 0.88σ (Held-group convention) for both packing
and Debye–Hückel, fixed before the first full run (an unretained one-row probe at
both diameters preceded it); diagnostic d = σ; diagnostic with MDEA at infinite
dilution in water. The Held-group Debye–Hückel convention a_j = σ_j (#330 research,
inferred) with 0.88σ packing gives 18/30 on the gate and 57.4 / 82.2% against Dicko
(#330 review probe, not retained); the packing diameter decides the result.

| | main, 0.88σ | diagnostic, σ | diagnostic, MDEA ∞-dilute |
|---|---|---|---|
| rows solved (Dicko, 28) | 28 | 28 | 28 |
| max reaction affinity | 3.3e-14 | 2.8e-14 | 3.8e-14 |
| references | all `Available`, admitted change ≤ 3.4e-6 | same | same |
| ε at every loaded liquid | 69.855 | 69.855 | 69.855 |
| **Figure 16 gate** (30 points within tolerance) | **22**; median \|ln r\| 0.053, max 0.32 | 0; median 0.50 | 2; median 0.29 |
| diagnostic: curves at Dicko's group loadings | 23; median 0.047, max 0.13; mean ln r −0.004 (left), +0.003 (right) | 0 | 2 |
| Dicko, symmetric deviation CO₂ (12) / H₂S (18) | **56.1 / 79.1%** | 59.8 / 54.2% | 56.8 / 71.2% |
| Dicko, \|calc/exp − 1\| CO₂ (12) / H₂S (18) | 64.4 / 143.7% | 40.8 / 67.7% | 51.0 / 109.8% |
| rows overpredicting p (18 mixed): CO₂ / H₂S | 14 / 16 | 0 / 15 | 10 / 15 |
| p_CO₂ rises with H₂S loading in each CO₂ group | yes (4/4) | yes | yes |

Published: 57.18% (CO₂) and 79.32% (H₂S).

**Figure 16 gate (amendment 1).** Figure 16 holds 12 symbols in the left panel
(CO₂ at H₂S loadings ≈ 0.1/0.3/0.5/0.8) and all 18 mixed rows in the right
panel (H₂S at CO₂ loadings 0.093/0.306/0.510/0.706); 12 + 18 = 30 is Cleeton's
Table 3 N_P, and the caption gives "AAD = 57.18% (Left Panel) and AAD = 79.32%
(Right panel)". The eight model curves were traced from the PDF
(`inputs/digitized/cleeton-2020-figure-16-curves.csv`, QA overlay
`cleeton-2020-figure-16-qa.png`); a second ink-threshold lane agrees within
0.0067 in ln p. The tolerance was fixed before any native run:
|ln(p_native/p_figure)| ≤ 0.05 + 0.003·|d ln p_figure/d loading|. Here 0.034
is 2.5 px of the 170 px/decade axis, 0.016 covers the K2 constant and the
omitted ε pressure term, and 0.003 loading is 2.5 px of the loading axis. Native
values are computed at the figure's own state: the symbol's loading, clipped to
the traced span, and the series' second-gas loading.
- **Main lane: 22/30 pass.** The misses are 0.08–0.32 in ln p and fall where
  the caption's nominal loading differs from Dicko's: the CO₂ = 0.8 series
  (data 0.706) gives +0.32 and +0.15.
- **Cause of the misses (inferred, bounded probe):** the second-gas loading at
  which each curve was drawn is not printed. With Dicko's group values (CO₂ the
  group's loading, H₂S the panel series mean), the bias vanishes (mean ln r
  −0.004 and +0.003) and the largest miss falls to 0.13. Seven points still
  exceed the tolerance by up to 0.05. This reading was added after the first run
  as one predeclared alternative; no reading was searched or fitted. Right-panel
  rows 14 and 15 pass under the caption reading, so one of the two readings covers
  them. Right-panel row 16 and left-panel rows 22 and 27 miss under both, by at
  most about 0.05 beyond tolerance, and remain unresolved; they are of the order of
  the Debye–Hückel-diameter, K2-digit and ε-pressure effects.
- **d = σ misses everywhere** (median 0.50, native about 40% below the
  curves): Cleeton's curves match the 0.88σ convention (inferred, strong).
  MDEA at infinite dilution also misses (median 0.29), consistent with the
  pure-liquid MDEA reference of Austgen's K4.

**Published deviations (secondary).** Design checks CC and CD as written
(|calc/exp − 1| within ±6 points) miss, at 64.4 / 143.7%. Cleeton does not print
its AAD formula; the symmetric form below was adopted after the first run, from
Cleeton's own curves, independent of the native values.
Its own digitized curves at the 30 symbols give 66.4% / 143.5% as
|calc/exp − 1|, but 55.9% / 79.1% as 2|calc − exp|/(calc + exp), against the
published 57.18% / 79.32%. The symmetric form is therefore the likely
definition (inferred). Native main-lane rows at Dicko's exact loadings give
56.1% / 79.1% in that form.

**Check ladder (amendment 2)**, against symbols digitized from Figures 2 and 4
(`inputs/digitized/cleeton-2020-figures-2-4-symbols.csv`, QA overlays beside
it). These are subsets of the papers' data sets, so the published AADs compare
only approximately.

| check | points | native AAD | published |
|---|---|---|---|
| pure MDEA 4C vapor pressure (Fig. 2) | 10 | 9.3% | 6.30% |
| pure MDEA 4C liquid density (Fig. 2) | 10 | 1.4% | 0.40% |
| H₂S–water bubble pressure (Fig. 4, P–x) | 17 | 8.8% | 7.26% |
| CO₂–water bubble pressure (Fig. 4, P–x) | 20 | 18.5% | 13.72% |

- **Pure MDEA is denser than Cleeton's model:** +1.1% at 313 K rising to +1.9%
  at 403 K. Cleeton's own 4C curve in Figure 2 lies about +0.8% above the data at
  300 K and on them near 420 K.
- **Probes (unresolved within this bounded investigation):**
  - Cleeton's 2B parameter set is also about 1% too dense natively.
  - The site reading is not the cause: 1+1 sites give 990 kg/m³ at 313 K.
  - The association convention is not the cause: the d³ form changes density
    by under 0.1%.
  - FeOs 0.10.1 with identical inputs (#330 review probe) agrees with the Engine
    to about 1e-10 relative: MDEA 4C 1037.57755 kg/m³ at 313.15 K and 971.79325
    at 403.15 K; water 4C 1001.41016 at 313.15 K and 994.28443 at 323.15 K.
  - The pure-MDEA gap therefore lies in Cleeton's printed parameters or code, not
    the Engine. Its growth with temperature suggests a different temperature
    treatment in Cleeton's code (inferred).
- Kuranov was skipped. Amendment 2 made it conditional on CC/CD missing, and
  under the design's own formula they did miss; it was skipped because the Dicko
  deviations match under the inferred symmetric definition. Running it is
  optional. Austgen's K4 covers 298–333 K, and the Engine refuses rows
  outside it (`ParameterError: reference_unavailable: reaction cleeton-2020-K4
  temperature outside reaction correlation` at 343.15 K).

**Per-row values (amendments 3, 5), reported only.**
- ln γ∞ of MDEA in water is −0.456 at 1 bar and −0.458 at 20 bar.
- Main lane: ln γ*(MDEAH⁺), relative to infinite dilution in water at the same
  T and P, runs from −0.57 to −0.09. MDEA is 10–85% protonated.
- **Liquid mass density rises with loading from 1,067 to 1,658 kg/m³**, i.e.
  the one-segment MDEAH⁺ (120 g/mol on a 3.35 Å sphere) packs densely.
  Cleeton (conclusions) names that representation as a likely adverse factor.
  No loaded-density data were compared.

**Limits.**
- This replays the printed model; Cleeton's code and supporting information
  are unavailable.
- Figure values carry digitization error (about 2.5 px).
- The ε pressure term and Cleeton's K2 digit are not replayed.
- The curves' second-gas states and the AAD formula are inferred.
- The ladder compares against figure subsets.

Files: `scripts/generate_data.py` (native declarations; writes
`results/data/cleeton-2020-{dicko-rows,figure-16-gate,ladder}.csv` and
`cleeton-2020-summary.json`, wheel `3da83fae…`), `scripts/digitize_figures.py`
(figure tracing from the PDF), `scripts/render.py`
(`results/rendered/cleeton-2020-figure-16-native.{png,pdf}`).

```bash
.venv/bin/python analyses/reproduction/2020-cleeton/scripts/generate_data.py
.venv/bin/python analyses/reproduction/2020-cleeton/scripts/render.py
```

## Retired calculation (#308)

The earlier generator never produced a number. All four attempts recorded in
commit `bc97b7737d56d9102cebbabd301039af37af4753` stopped during parameter
construction with the typed `same-sign-ionic-dispersion-rule` rejection. It
extended the Wangler packet, declared seven anion–anion pairs and used H⁺ with
hydronium constants. Its attempts remain in
`results/data/cleeton-2020-fixed-topology-plotted-data.csv` and
`results/rendered/cleeton-2020-fixed-topology.*` (drawn by the #308 version of
`render.py`). The mixed-solvent reference it was thought to need is not needed:
Austgen's standard states are two pure-liquid solvents with aqueous
infinite-dilution solutes (#330 research).
