# Optional Gross–Vrabec dipole–dipole term: qualification and model comparison (#28)

**Question.** With the GV2006 dipole–dipole term implemented as an opt-in family
(ADR 0011), does PCP-SAFT with its own published parameters describe polar fluids
and a polar/alkane binary better than nonpolar PC-SAFT with its own parameters?
And where does it not help?

**Answer (bounded).**
- **Where it helps:**
  - **Acetone saturated-liquid density.** It is described about 2.4 times better.
  - **Acetone–n-pentane VLE.** Prediction error roughly halves with k_ij = 0, and the calibrated k_ij is about half that of the nonpolar model.
  - **DMSO (3.96 D) saturated-liquid density.** It is described about 2.2 times better on measured data.
- **Weak dipoles.** For dimethyl ether and diethyl ether the term changes little, as GV2006 reports.
- **Where it did not help:** vapor pressure of acetone and DMSO. With the published polar sets, the predeclared vapor-pressure criterion fails for both. When both models are refit on the same rows, the polar model is the better one. Neither vapor-pressure gap is larger than the uncertainty of its reference data.

## Implementation qualification (falsifiers)

| Claim | Check | Result |
|---|---|---|
| a_dd equals the GV2006 equations | Engine Δa^r against teqp 0.23.2 `PCSAFTMixture` at fixed (T, ρ, x). Cases: pure acetone (liquid, vapor and dilute), acetone/n-pentane, DME/acetone/n-pentane (`results/teqp-dipolar-rows.csv`, `scripts/teqp_reference.py`) | max relative difference 9.7e-16 (criterion 1e-12), 8 rows |
| Low-density limit | a_dd/ρ̂ → Eq. 8 at η = 0 (B_dd = −344.4 Å³, acetone at 300 K), written independently in the test | intercept error 3.5e-9; error ratio per density decade 9.998 and 9.984 (O(ρ)) |
| Derivatives | Residual-enthalpy actions of orders 1–3 in T, ρ and x vs Richardson mixed stencils at acetone/n-pentane, 372.7 K, 9000 mol/m³ | within the predeclared 1e-8, 1e-6 and 1e-4; the dipolar share exceeds each criterion by more than 100×. Z − 1 = ρa_ρ holds to < 1e-9 and u_res = −RT²a_T to < 1e-8 |
| Absent dipole | x_acetone = 0 with the family declared | a^r, P and μ_res are bitwise equal to the nonpolar model; a_dd ∝ x² near 0 |
| Opt-in | `record_native_rows.py` rerun on the rebased build | all 444 native oracle values unchanged; check-all passes |

The persistent checks are `engine/tests/eos/test_dipolar.py`.

## Complete model/parameter pairs compared

- **Polar model.** It uses `data/packets/gross-vrabec-2006/1`: GV2006 Table 2 with μ and n_μ = 1.
- **Nonpolar model.**
  - Dimethyl ether and diethyl ether use the Gross–Sadowski 2001 published sets.
  - The library holds no published nonpolar acetone or DMSO set (`inputs/sources/observations-inventory.md`). Both therefore use a documented refit of (m, σ, ε/k) on the calibration rows only.
- **Secondary rows.** Both models are refit on the same calibration rows, removing the in-sample advantage.
- **Data.** Rows alternate calibration and validation. They are **not** the KDB/DIPPR data behind GV2006 Table 2, which is not openly available, so absolute AADs are measured against the references below, not reproductions of GV2006's numbers. Each row carries the uncertainty its source states; blank means none is stated.

| Fluid | Reference (rows) | Stated uncertainty, p_sat | Stated uncertainty, ρ_liq |
|---|---|---|---|
| acetone | Lemmon and Span 2006 EOS, valid 178.5–550 K; 180–490 K | 0.25 % at 290–390 K, 0.5 % above 270 K; none below 270 K (7 rows) | 0.1 % at 280–310 K, 0.5 % below 380 K, 1 % above |
| dimethyl ether | Wu, Zhou and Lemmon 2011 EOS, valid triple point–550 K, k = 2; 200–385 K | 0.2 % above 230 K; none below (4 rows) | 0.05 % "except near the critical point"; the 385 K row (0.96 T_c) is treated as unstated |
| diethyl ether | Thol, Piazza and Span 2014 EOS, valid 270–500 K; 270–450 K | 0.4 % (deviation envelope; all rows) | 0.2 % (deviation envelope) up to 0.95 T_c; the 450 K row (0.96 T_c) is treated as unstated |
| DMSO | measured, NIST ThermoML Archive (below); p_sat 293–463 K, ρ_liq 293–473 K | 1.1–10 % (expanded) | 0.08–0.11 % (expanded) |

  - The EOS rows come through CoolProp 8.0.0 (`scripts/reference_saturation.py`, `inputs/reference-saturation.csv`). Uncertainties are quoted from the EOS papers: Lemmon and Span, *J. Chem. Eng. Data* 51 (2006) 785, acetone subsection; Wu et al., *J. Phys. Chem. Ref. Data* 40 (2011) 023104, abstract and conclusions.
  - **Diethyl-ether correction.** CoolProp's diethyl ether is Thol, Piazza and Span, *Int. J. Thermophys.* 35 (2014) 783, not Thol et al. 2012, and it is valid only from 270 K. The earlier 220–260 K rows were extrapolations below that limit; the grid now starts at 270 K. The paper (Zotero `XRMZRYJI`, PDF `LSLF8W4R`) states no separate uncertainty. Sect. 3.1 gives deviation envelopes against the fitted data: vapor pressures "less than ±0.4 %" (p. 800), and saturated liquid densities "within ±0.2 %", with deviations increasing close to the critical point (p. 801). Validity is 270–500 K (p. 799). For pressure, p. 799 gives 100 MPa, while the abstract (p. 783) and the conclusions (p. 805) give 40 MPa; the paper is inconsistent. The paper takes the average deviation to the best data set as its indication of EOS uncertainty (p. 790). The envelopes are used here as the reference bounds; this is an upper, not a k = 2, statement.
  - **DMSO rows** (`scripts/dmso_observations.py`, `inputs/dmso-saturation.csv`) come from the ThermoML Archive v2020-09-30 (doi:10.18434/mds2-2422):
    - vapor pressure from Fulem et al. 2011 (static, 293–308 K), de Pablo et al. 2018 (static, 353–433 K) and Zhang et al. 2013 (ebulliometric, 402–463 K), with replicates averaged;
    - liquid density from Paknejad et al. 2019 (vibrating tube), at the lowest measured pressure: 100 kPa up to 333 K, 500 kPa above. Compression from p_sat to that pressure is below the 1 kg/m³ reporting step.
    - The sets agree within 1–2 % with the older Douglas 1948 and Jakli–Van Hook 1972 correlations in the NIST WebBook.
- **Is a gap resolved?** A gap is the polar AAD minus the nonpolar AAD. It is resolved when it exceeds the most that a reference error within the stated uncertainty could shift it. That worst case is computed row by row (`results/summary.json:reproduction.reference_bound`).

### Pure fluids (AAD-%, all rows / validation rows)

| Fluid | Model (own parameters) | m, σ/Å, ε/k/K | p_sat | ρ_liq |
|---|---|---|---|---|
| acetone (2.88 D) | PCP-SAFT, GV2006 | 2.7447, 3.2742, 232.99 | 1.31 / 1.04 | **0.93 / 0.97** |
| | PC-SAFT, refit | 2.7274, 3.2819, 255.77 | **0.93 / 0.94** | 2.24 / 2.43 |
| | PCP-SAFT, refit (secondary) | 2.7924, 3.2482, 230.96 | 1.40 / 1.31 | 0.61 / 0.67 |
| dimethyl ether (1.3 D) | PCP-SAFT, GV2006 | 2.2634, 3.2723, 210.29 | 0.28 / 0.26 | 1.23 / 1.20 |
| | PC-SAFT, GS2001 | 2.3071, 3.2528, 211.06 | 0.33 / 0.31 | 1.41 / 1.35 |
| | refits, polar / nonpolar | — | 0.15 / 0.13 | 0.35 / 0.53 |
| diethyl ether (1.15 D) | PCP-SAFT, GV2006 | 2.9726, 3.5127, 219.53 | 0.48 / 0.49 | 0.69 / 0.70 |
| | PC-SAFT, GS2001 | 2.9686, 3.5147, 220.09 | **0.39 / 0.39** | 0.71 / 0.71 |
| | refits, polar / nonpolar | — | 0.32 / 0.33 | 0.33 / 0.35 |
| DMSO (3.96 D) | PCP-SAFT, GV2006 | 3.0243, 3.2427, 309.36 | 1.71 / 1.68 | **0.39 / 0.36** |
| | PC-SAFT, refit | 2.8119, 3.3242, 364.52 | **0.88 / 0.82** | 0.86 / 0.79 |
| | PCP-SAFT, refit (secondary) | 3.0189, 3.2464, 310.16 | 0.50 / 0.44 | 0.41 / 0.39 |

**Predeclared reproduction disposition** (`results/summary.json:reproduction`). The test is whether our PCP-minus-PC gap has GV's sign and at least half GV's magnitude.

| Fluid, property | Our gap | GV gap | Criterion | Largest shift from reference error (rows with stated uncertainty) |
|---|---|---|---|---|
| acetone ρ_liq | −1.31 | −1.63 | pass | 0.48, resolved |
| acetone p_sat | +0.38 | −1.06 | fail | −0.40 on rows above 270 K; shift up to 0.70, not resolved |
| dimethyl ether p_sat | −0.05 | −0.01 | pass | −0.055 above 230 K; shift up to 0.053, barely resolved |
| dimethyl ether ρ_liq | −0.18 | −0.27 | pass | −0.20; shift 0.00, resolved (same-sign errors larger than 0.05 %) |
| diethyl ether p_sat | +0.09 | −0.02 | fail | shift up to 0.15, not resolved |
| diethyl ether ρ_liq | −0.01 | −0.03 | fail | −0.007 on rows ≤ 0.95 T_c; shift up to 0.014, not resolved |
| DMSO p_sat | +0.83 | −0.99 | fail | shift up to 2.68, not resolved |
| DMSO ρ_liq | −0.47 | −0.25 | pass | shift 0.00, resolved |

- **Acetone p_sat.** The GV parameters give +10 % at 180 K (p = 2.9 Pa), below the PCP-SAFT fit range of 200–508 K (Kleiner and Sadowski 2007, Table 3) and where the EOS states no uncertainty. Within 200–508 K the gap is −0.21, the right sign but under half of GV's magnitude. On rows above 270 K it is −0.40, which the 0.25–0.5 % EOS uncertainty cannot resolve. This is a limit of the reference, not a demonstrated failure of the dipolar term.
- **Diethyl ether.** Earlier this passed only on rows below 270 K, where the EOS was extrapolated.
  - **Under the envelopes** (the bounds used here), both gaps fall within the Thol 2014 bounds: p_sat +0.09 against a possible reference shift of 0.15, and ρ_liq −0.007 against 0.014. GV's gaps (≤ 0.03 points, printed to 0.01) are smaller still. The weak 1.15 D dipole then has no resolvable effect, as GV reports.
  - **This conclusion depends on the envelope choice.** Thol's p. 790 rule, which takes the AAD to the best data set as the uncertainty indication (Table 7, p. 797), gives much tighter values: 0.168 % for p_sat (Ambrose et al., 0.6–0.98 T_c) and 0.041 % for ρ_liq (Hales et al.).
    - Under it, the bounds become 0.056 and 0.0008, and both gaps are resolved.
    - The p_sat gap then has the opposite sign to GV: polar is worse.
    - The ρ_liq gap has GV's sign but less than half GV's magnitude.
    - Recomputed with `bound()` on `results/pure-rows.csv`.
  - **Absolute errors.** The p_sat AADs are 0.48 % (polar) and 0.39 % (nonpolar). The polar value exceeds the ±0.4 % envelope by 1.2 times; against the 0.168 % indication, both are model error.
  - The ρ_liq AADs (0.69 %, 0.70 %) are 3.5 times the ±0.2 % envelope, so they are model error. Both refits reduce them to 0.33–0.35 %.
- **DMSO p_sat.** The GV parameters overpredict every measured vapor pressure except at the normal boiling point:
  - by +1.4 to +2.5 % against Fulem et al. (293–308 K, uncertainty 1.4–2.1 %);
  - by +1.0 to +2.5 % against Zhang et al. (402–456 K).
  
  GV2006 reports 0.46 % against KDB, so the KDB rows and these measurements probably differ by about 2 %. The nonpolar advantage (0.88 %) comes from being fit to these rows. Refit on the same rows, the polar model does better: 0.50 vs 0.88 %.
**Identifiability.** The pure refits are locally conditioned: condition number 18–63, full rank. ε/k and m are correlated at −0.98 to −1.00 for every fluid, so the individual coefficients are weakly identifiable; only their combination is well determined. For acetone and DMSO, m exceeds the cap of 2, so only the J-integral coefficients stop depending on m. The reduced moment q_i ∝ 1/m_i and η still carry it.

### Acetone–n-pentane VLE (Campbell et al. 1986: GV2006 Fig. 15 digitized, and Table I)

**Data.**
- 30 symbols, all left of the azeotrope (x_pentane ≤ 0.60), where the branches separate: filled symbols are liquid bubble points, open symbols vapor dew points. Three detections were rejected visually as curve artifacts (listed in the script). They are digitized by `scripts/digitize_fig15.py` from a 600-dpi render, with u_x = 0.004 and u_P = 10 kPa (≈ 1–2 % of P).
- **Printed paired x–y rows.** Campbell et al. 1986 (Zotero `R9D3UYNT`, PDF `2LB6R2IN`) are transcribed from Table I (p. 426) into `inputs/campbell1986-table1.csv`, checked against the rendered page: 27 binary rows plus 6 pure-component rows, with x and y in n-pentane. The Table I "calculated" P and y columns come from the paper's four-suffix Margules eq 8, fitted per isotherm to P–x only by Barker's method (Table III, p. 427). Its y is therefore a Gibbs–Duhem prediction, and its P a fitted residual.
  - Stated uncertainties: composition precision 0.005 mole fraction (p. 425), temperature 0.1 K (p. 425), pressure read to 0.1 psi (0.7 kPa; p. 424). The p. 427 error bars use 0.005 and 0.1 psi.
- **Partitions:** k_ij is calibrated at 372.7 K and validated at 397.7 and 422.6 K. n-Pentane uses Gross–Sadowski 2001 in both models.

| Model | k_ij | bubble P AAD-%, validation | dew P AAD-%, validation |
|---|---|---|---|
| PCP-SAFT, GV2006 acetone | 0 (prediction) | 10.7 | 8.8 |
| PC-SAFT, refit acetone | 0 (prediction) | 18.6 | 16.3 |
| PCP-SAFT | 0.033 (calibrated; GV: 0.024) | **1.26** | **0.50** |
| PC-SAFT | 0.059 (calibrated; GV: 0.072) | 2.06 | 0.92 |

There were no failed calculations. The calibrated difference, about 1 percentage point, is comparable to the digitization uncertainty. The k_ij = 0 difference is not.

**Paired x–y (#161).** `scripts/binary.py` solves a bubble point at each printed (T, x) with the k_ij above, which are still calibrated on the digitized 372.7 K symbols. It reports the mean |y_model − y_data| of n-pentane and the bubble-P AAD (`results/binary-xy-rows.csv`, `results/binary-xy-summary.csv`). There are 8 rows at 372.7 K and 19 at the validation isotherms, with no failed calculations.

| Model | k_ij | mean \|Δy\|, 372.7 K / validation | max \|Δy\|, validation | bubble P AAD-%, validation |
|---|---|---|---|---|
| PCP-SAFT, GV2006 acetone | 0 | 0.036 / 0.028 | 0.060 | 9.8 |
| PC-SAFT, refit acetone | 0 | 0.069 / 0.055 | 0.109 | 17.0 |
| PCP-SAFT | 0.033 | **0.0098 / 0.0102** | 0.031 | **1.65** |
| PC-SAFT | 0.059 | 0.0171 / 0.0159 | 0.038 | 2.85 |
| Campbell four-suffix Margules (eq 8), Barker P–x fit per isotherm (Table I) | — | 0.0084 / 0.0068 | 0.016 | 0.19 |

The printed rows confirm the order found on the digitized envelope. With its calibrated k_ij, PCP-SAFT predicts the validation vapor composition to 0.010 mean |Δy|: about twice the stated composition precision of 0.005, and 1.5 times the scatter floor of the paper's own per-isotherm correlation (0.007). PC-SAFT gives 0.016. The polar advantage (0.006 in |Δy|, 1.2 points in bubble P) exceeds the 0.005 composition precision only marginally in y, but clearly in P.

## Limits

- The claims hold only for these fluids, temperature ranges and data sets.
- Absolute AADs are against the named reference EOS or measured DMSO data, not the KDB/DIPPR data behind GV2006, so they are not reproductions of GV2006 Table 2.
- The diethyl-ether reference bounds are Thol 2014's deviation envelopes, not k = 2 uncertainties.
- There is no electrolyte or ion claim, and no composition derivative above first order along an absent dipolar species.
- The dipole–quadrupole cross term remains retired (ADR 0006); the quadrupole term is admitted by ADR 0014 ([2026-polar-quadrupole](../2026-polar-quadrupole/README.md)).

## Reproduce

```bash
/tmp/teqp/bin/python analyses/2026-polar-dipole/scripts/teqp_reference.py         # teqp 0.23.2, CoolProp 8.0.0 venv
/tmp/teqp/bin/python analyses/2026-polar-dipole/scripts/reference_saturation.py  # any CoolProp 8.0.0 venv
python analyses/2026-polar-dipole/scripts/dmso_observations.py /tmp/thermoml        # extracted ThermoML.v2020-09-30.tgz
python analyses/2026-polar-dipole/scripts/digitize_fig15.py /tmp/p8-08.png        # numpy/scipy/pillow
.venv/bin/python analyses/2026-polar-dipole/scripts/compare.py                    # installed wheel
.venv/bin/python analyses/2026-polar-dipole/scripts/binary.py
```

The results were produced with wheel SHA-256 `c74a91f66fb906c56a7cbad9e9448760eb7e24c79afc740c5391c061d10b220a`. `results/summary.json` records the wheel of the latest run.

## Native pure-fluid starts

The [native comparison](results/batch3-comparison.json) retains the common
optimizer chart, all starts and conditional covariance. DMSO's selected polar
start changes between numerically tied fits without changing the model ranking.
Unavailable polar starts remain unavailable. Binary-mixture evidence keeps its
historical wheel and was not recomputed by this conversion.
