# Wangler 2018: water–MDEAH⁺ fit and Table 10 prediction

Wangler et al. (2018) fitted the water–MDEAH⁺ binary to ten osmotic coefficients of
aqueous MDEA/HCl (Tables 8–9; Figure 10, 5.05%). They then predicted the 60 wt% MDEA +
CO₂ bubble pressures at 303 K (Table 10, ≈10.25%) without adjusting anything to them
(l.27, l.658). The printed water–MDEAH⁺ k_ij row (Table 6) copies water–MDEA, so #330
fits that k_ij here natively and repeats the paper's prediction.

**Result (#330):** Figure 10 is a calibration: the k_ij is fitted to it, and it is
reproduced at 2.37% / 2.16%. Table 10 is the out-of-sample prediction: 96.2% / 89.6%
ARD against the paper's ≈10.25%. It is an unresolved model-reading gap, not a fitting
artefact (see Meaning).

## Conventions (#330 research, comment 6037145321)

- MDEAH⁺ m = 3.675, with σ and u/k transferred from MDEA (inferred: "non-spherical").
- Packing diameter 0.88σ (Eq. 7) and Debye–Hückel diameter σ (Held 2008/2012/2014).
- l_hb(water–MDEAH⁺) held at the printed Table 7 values (0.345 at 273 K, 0.388 at 313 K,
  0.377 at 303 K).
- Figure 10 permittivity: water, from the one-bar Archer–Wang approximation the
  project's legacy reference uses. Wangler cites a different correlation. Using the
  water/MDEA blend instead moves Figure 10 by about 1 point (#330 research).
- Two lanes, differing only in ion self-dispersion:
  - **S0**, excluded (Held 2014, Uyan 2015): the best-supported reading and the primary lane.
  - **S1**, kept (literal l.146). For Figure 10's 1:1 salt this is `all-ion-pairs`; in
    Table 10 it adds k_ij = 1 on the same-sign cross pairs.

## Figure 10: water–MDEAH⁺ k_ij fitted per isotherm

The objective is Wangler's OF3 (Eq. 33), Σ(φ/φ_obs − 1)², minimized by native Ceres
Levenberg–Marquardt with exact Jacobians. There is one k_ij per isotherm, which is Eq. 34
truncated to two terms, with bounds −1 to 1 and two declared starts, 0 and 0.5. The
allowance δk comes from u(φ) = 0.005 by the #245 propagation:
δk = Σᵢ |J⁺|ᵢ u(φ)/φ_obs,ᵢ
([summary](figures/figure-10/output/wangler-2018-mdeah-osmotic-summary.json),
[plotted values](figures/figure-10/output/wangler-2018-mdeah-osmotic-plotted-data.csv)).

| Lane | k_ij 273 K | k_ij 313 K | δk 273 / 313 K | starts differ by | ARD (paper 5.05%) | F1 |
|---|---:|---:|---:|---:|---:|---|
| S0 | +0.3495 | +0.3259 | 0.0020 / 0.0022 | ≤ 1.2e-9 | **2.37%** | pass |
| S1 | +0.0282 | +0.0144 | 0.0021 / 0.0024 | ≤ 1.1e-8 | **2.16%** | pass |

- **F1** requires an ARD of 6.05% or less, interior endpoints with available covariance,
  and starts that agree within δk/100. Both lanes pass.
- The printed Table 6 k_ij gives 129.4% (S0) and 46.7% (S1) on the same points.
- The model's φ falls with molality at 273 K, while the data rise from 0.84 to 0.86. A
  constant k_ij cannot follow that trend, so the residual is model discrepancy, not
  measurement scatter.
- Accordingly, the standard errors from the covariance scaled by the residual variance
  (SSE/(n − 1)) are about twice δk: 0.0040 / 0.0050 (S0) and 0.0041 / 0.0045 (S1) at
  273 / 313 K. The 313 K residual is model misfit too (data 0.869 → 0.861 → 0.874). F3
  uses δk, as the design fixed.
- **Interpolated to 303 K (Eq. 34, linear):** k_ij = +0.3318 ± 0.0021 (S0) and +0.0178 ± 0.0023 (S1).

## Table 10: predicted bubble pressures at 303 K

The declaration is unchanged from #308/#330:
- Table 2 constants on the mole-fraction basis (Eq. 19), unconverted, referenced to pure
  water and to infinite dilution in water for every other species (Eqs. 13–21);
- the fixed 60 wt% blend permittivity (ε = 44.1) on every neutral species;
- a CO₂-rich vapour start.

Each fitted lane reads k_ij(303 K) and δk from the Figure 10 summary. The printed lane
keeps Table 6's k_ij and needs a declared 257-bar seed: from the observed-pressure seed,
rows 2–3 stop at `LOCAL_INFEASIBILITY`. The fitted lanes are seeded at the observed
pressure
([table](figures/figure-14/output/wangler-2018-table-10-native.csv),
[summary](figures/figure-14/output/wangler-2018-table-10-native-summary.json)).

| Lane | p (bar), against 9.2 / 14.2 / 24.5 | ARD (paper ≈10.25 ± 3) | F2 | F3: p at k ∓ δk |
|---|---|---:|---|---:|
| fitted S0 | 0.464 / 0.539 / 0.593 | **96.2%** | miss (ARD) | ±3.7% |
| fitted S1 | 0.714 / 1.566 / 3.075 | **89.6%** | miss (ARD) | ±1.9% |
| printed S0 (Table 6 k_ij = −0.190) | 196.9 / 265.9 / 322.8 | 1677% | – | – |

- **F2:** every lane solves 3/3 rows. The reaction affinity is at most 2.5e-14, and every
  reference is admitted, with changes K1–K3 ≤ 8.5e-6. Only the ARD criterion misses.
- **F3:** the δk band moves pressures by at most 3.7%, so the miss does not come from
  the Figure 10 fit's measurement allowance.

Reported per row, without gating:

| Lane | α | ln γ*(MDEAH⁺) | protonated fraction | liquid ρ (kg/m³) | k_ij(303 K) matching the row |
|---|---:|---:|---:|---:|---:|
| fitted S0 | 0.9374 | −5.56 | 0.820 | 1261 | +0.106 |
| fitted S0 | 0.9973 | −5.61 | 0.840 | 1260 | +0.069 |
| fitted S0 | 1.0448 | −5.67 | 0.856 | 1259 | +0.011 |
| fitted S1 | 0.9374 | −3.20 | 0.998 | 1504 | −0.141 |
| fitted S1 | 0.9973 | −3.25 | 0.999 | 1510 | −0.147 |
| fitted S1 | 1.0448 | −3.28 | 1.000 | 1511 | −0.164 |
| printed S0 | 0.9374–1.0448 | +10.2 to +10.5 | 0.19–0.20 | 1131–1143 | – |

- ln γ*(MDEAH⁺) is taken relative to infinite dilution in water at the same T and P.
- The matching k_ij is a native one-coordinate fit to that row's pressure alone, used as a
  diagnostic. No lane uses it.
- In the printed lane, the second phase is a dense CO₂-rich fluid (877–949 kg/m³).

**Meaning (inferred):**
- No single k_ij(303 K) serves both data sets. Figure 10 needs +0.33 (S0) or +0.018 (S1);
  Table 10 needs +0.01 to +0.11 or −0.14 to −0.16.
- In S0, the k_ij that matches each row also changes with loading.
- The transfer diagnosis (#330 comment 6045486791) ran one-at-a-time declared probes at
  row 1. Our declaration matches Wangler Tables 4–7 and Uyan 2015 Tables 3–5 pair by pair
  and site by site; no pair or site is omitted. ln γ*(MDEAH⁺) there is a near-cancellation
  of the hard-chain term (−25.0) and the water–MDEAH⁺ attraction (+27.6), so only the
  water–MDEAH⁺ k_ij carries order-8 leverage (+0.106 gives 9.23 bar at row 1). Pairs Figure
  10 does not sample reach at most 0.078× the measured pressure under source-allowed
  changes, 0.27× under source-contradicted ones and 0.56× under unsourced ones.
- The Engine's multi-segment ion terms match FeOs on a Figure 10 state (#337, one state).
- Inferred, not reproduced: Wangler's water–MDEAH⁺ binary at 303 K differs from a Figure
  10 fit under our reading, through an undisclosed k_ij(T) form (the printed row is a
  copy) or an unprinted implementation difference in the cavity–attraction balance.
- The gap is unresolved. No parameter is tuned to Table 10. The open step is asking the
  Held group for the fitted water–MDEAH⁺ k_ij(T) and l_hb(T) (owner option (a)).

## Retired and frozen records

- **#234 batch 5:** the native refit of k_ij and the association volume, using the printed
  k_ij and 0.88σ Debye–Hückel diameters, reached the association-volume bound. Its
  comparison code is retired. `figures/figure-10/output/batch-5-*` keep its values
  (summary, before/after, exact columns, derivative verification) as frozen records. They
  are not regenerated: the new conventions would trip its 1e-7 "source unchanged" guard
  by design.
- **#330 rerun (`a1238600`):** the printed lane with 0.88σ Debye–Hückel diameters gave
  1617% ARD, and the m = 1 diagnostic lane 50.3%, at an unphysical liquid density near
  2,090 kg/m³. Both lanes are superseded here and remain in Git history.
- **Historical Figure 14 record** (commit `52c165f5`, no verified wheel identity): a
  water/MDEA blend reference through a retired Python activity-scale projection. It gave
  65.17% ARD with the paper parameters and 60.82% with a binary refit. The
  `wangler-2018-table-10-grepe-*` files keep it frozen.

## Commands

From the repository root:

```bash
bash tools/build-engine.sh --install-wheel
.venv/bin/python analyses/reproduction/2018-wangler/figures/figure-10/scripts/generate_data.py
.venv/bin/python analyses/reproduction/2018-wangler/figures/figure-14/scripts/generate_data.py
.venv/bin/python analyses/reproduction/2018-wangler/figures/figure-10/scripts/render.py
.venv/bin/python analyses/reproduction/2018-wangler/figures/figure-14/scripts/render.py
```

Both data scripts record the installed wheel's SHA-256, so install the checkout's wheel first;
rerunning `bash .codex/environments/setup.sh` restores the editable install. Figure 14 reads
the Figure 10 summary, so run Figure 10 first. The render scripts plot
retained values without evaluating the Engine.
