# MEA 2B association volume: d³ versus σ³ convention

**Frozen record — not rerunnable as a current-main calculation.** Record snapshot:
commit `ca719012c31e8b41edcb98ae82a388117b58d343`. Original wheel identities, source
inputs and retained numerical evidence keep their recorded provenance. The results used
a Python solver before the native calculation rule; the pure-MEA association-volume fit
was removed, although fixed-parameter declarations remain.

Pre-rule record: the retained results were produced by a Python solver before the native calculation rule (#234, 2026-10-01); they are history, not reproducible by the current scripts. The pure-MEA association-volume fit with the absolute-relative-error sum estimator was removed from `scripts/mea_saturation.py`; fixed-parameter saturation comparisons remain.

Issue #165. Fakouri Baygi & Pahlavanzadeh 2015 (Zotero `B8EJJPCJ`, PDF `JWH69DKG`)
define the association strength in Eq. 4 (p. 791) as

`Δ^{AiBj} = d_ij^3 g_ij^hs(d_ij) κ^{AiBj} [exp(ε^{AiBj}/kT) − 1]`, with `d_ij = (d_i + d_j)/2` and
`d_i = σ_i[1 − 0.12 exp(−3ε_i/kT)]` (Eq. 5).

The Engine uses `σ_ij^3` ([`association_strength`](../../engine/docs/equations.md#association_strength)).
Their Table 2 (p. 792) MEA 2B set is `m = 3.0353`, `σ = 3.0435 Å`, `ε/k = 277.174 K`,
`κ = 0.037470`, `ε^AB/k = 2586.3 K`, with AAD 0.62 % in `p_sat` and 0.12 % in `ρ_L,sat` against
their Table 1 correlations (Avlund 2008 DIPPR, Eqs. 9–10) over 303.15–443.15 K.

## Diameter ratio

`Δ` is linear in `κ`. At a fixed temperature, the d³ model is therefore the Engine model with
`κ(T) = κ (d_ij/σ_ij)^3`. MEA–water uses the packet's Held water with its σ(T) correlation
([rows](results/association-diameter-ratio.csv)).

| T / K | MEA–MEA (d/σ)³ | MEA–water (d_ij/σ_ij)³ |
|---:|---:|---:|
| 293.15 | 0.9790 | 0.9845 |
| 313.15 | 0.9749 | 0.9811 |
| 333.15 | 0.9706 | 0.9776 |
| 353.15 | 0.9662 | 0.9739 |
| 373.15 | 0.9617 | 0.9701 |
| 393.15 | 0.9572 | 0.9662 |
| 443.15 | 0.9459 | 0.9562 |

MEA's dispersion energy is low compared with T, so `d/σ` stays near 1. The correction is 2–5 %, not the
10–15 % estimated when the issue was opened.

## Pure MEA saturation

The retained saturation comparison uses a 5 K grid from
293.15 to 443.15 K ([rows](results/mea-saturation.csv), [summary](results/mea-saturation.json)).
AAD is in percent. The p_sat range is the signed minimum and maximum deviation.

| model | κ | Baygi range 303–443 K: p_sat AAD / ρ_L AAD | packet domain 293–393 K: p_sat AAD (range) |
|---|---:|---:|---:|
| Baygi Table 2 (reported) | 0.037470 | 0.62 / 0.12 | – |
| d³, Baygi Eq. 4 | 0.037470 | 0.613 / 0.109 | 0.79 (−0.97 to +1.87) |
| σ³, packet v4 | 0.037470 | 3.451 / 0.126 | 3.09 (−4.15 to −0.15) |
| σ³, packet v5 | 0.036287 | 1.106 / 0.117 | 1.22 (−1.55 to +2.93) |

The d³ evaluation reproduces Baygi's reported deviations, which identifies d³ as the convention they
regressed with. The unconverted value in the Engine underpredicts MEA vapor pressure by 1–4 % over
Baygi's range.

Packet v5 carries the constant `κ` that best reproduces Baygi's d³ model (Eq. 8 form, p_sat and ρ_L)
over the packet domain, with `m`, `σ` and `ε` fixed. The optimum is 0.0362870 (0.9684 κ), which equals
(d/σ)³ near 343 K. Against the d³ model over 293–393 K:

- v5: p_sat AAD 0.53 %, range −0.95 to +1.04 %; ρ_L AAD 0.003 %.
- v4: p_sat AAD 2.80 %, range −3.53 to −1.98 %.

No constant conversion is exact, because (d/σ)³ falls by 2.2 % across the domain. The remaining ±1 %
is v5's error against the source model. At 293 K, v4 lies closer to the correlation (−0.15 %) than
Baygi's own model does (+1.87 %). That is a coincidence of offsetting errors, not a better model.

Under the Wolbach–Sandler rule, the MEA–water cross volume scales by `√0.9684 = 0.984`. That is within
−0.04 to +1.9 % of the d³ cross factor (d_ij/σ_ij)³ over 293–393 K. That factor pairs Baygi's MEA with the
packet's Held water, so it is not a model Baygi ran; the `k_ij` refit below bounds its effect.

## MEA–water VLE (Cai 1996)

[`scripts/mea_water.py`](scripts/mea_water.py) reuses the #155 bubble-point replay on the 25
admitted Cai rows ([rows](results/mea-water-cai.csv), [summary](results/mea-water-cai.json)).

| quantity | v4 | v5 |
|---|---:|---:|
| mean scaled square at packet `k_MEA,H2O = −0.0735275` | 4.4265 | 4.4185 |
| max \|ΔT_bubble\| vs v4 | – | 0.55 K |
| max \|Δy_water\| vs v4 | – | 0.0049 |
| one-coordinate refit of `k_MEA,H2O` | – | −0.0734191 (4.4182) |

The refit lies inside the retained 95 % interval `[−0.07792, −0.06914]`, so v5 keeps the packet
`k_ij`.

## Convention survey

Only Baygi 2015 and Chapman 1990 (Eq. 18, SAFT) write `d^3`. Nasrifar & Tafazzol 2010 (PC-SAFT
MEA, Eq. 34) writes `σ_ij^3`. Wangler 2018 (Eq. 12), Uyan 2015 (Eq. 9), Bülow 2021 (Eq. 25) and
Gross & Sadowski 2002 state only the cross-volume rule and do not write Δ. No second packet source
documents a `d^3` convention, so the Engine keeps one convention and converts the volume in the packet.

Baygi's 3B and 4C sets need the same conversion, with larger factors (their ε/k is 204 and 238 K,
against 277 K for 2B). The independent review of #165 found the same signature for both: d³
reproduces Baygi's AAD and σ³ misses it by 5–6 % in p_sat. That check is not retained here.

## Frozen comparison limits

Claim limits: pure MEA and the 2B set only; the reference is Baygi's correlations, not raw
measurements, on a 5 K grid (Baygi's point set is not stated). The converted value is a
one-coordinate fit to the d³ model's saturation values. A constant κ also changes dΔ/dT relative to
the d³ model: from the p_sat ratios, the heat of vaporization differs by roughly 0.2 kJ/mol (about 0.3 %),
an independent-review estimate that is not retained. The CO2 dispersion energy and R4 correlation calibrated by the 2026-09-02 MEA
campaign on v3 values, and every reactive result, are not re-evaluated here.
