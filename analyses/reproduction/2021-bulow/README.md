# Bülow, Ascani & Held 2021 — ePC-SAFT advanced II

Bülow, M.; Ascani, M.; Held, C. *ePC-SAFT advanced – Part II: Application to salt
solubility in ionic and organic solvents and the impact of ion pairing.* Fluid
Phase Equilibria 537 (2021) 112989, with its corrigendum (FPE 548 (2021) 113183)
and supporting information. Zotero parent `QYRA4L6M`: paper PDF `Y3P62JXY`
(SHA-256 `7d825b0c…ab9be5`), transcription `SITUIV2V`, SI `AGLEWE3W`
(`90343331…c7533d8ea`); corrigendum parent `5E8Y9396`, PDF `ST7JW4WN`
(`17b68fe4…ce805aad`).

**Result** (numbers in [`notebook.qmd`](notebook.qmd), computed from `results/`):
the outcome depends on an unresolved inorganic-ion diameter convention, and both
variants are kept. With σ, the convention Part II states (Section 2.4), the Engine
reproduces neither the corrected Fig. 4 / SI Fig. S1 curves at high salt nor any
Fig. 8b bar. With 0.88σ for the hard-sphere, Debye–Hückel and Born diameters
together (Held 2014 Eq. 21, cited in Section 3.1), the Fig. 4 / S1 advanced curves
match over their full ranges and 10 of 18 Fig. 8b bars match, including all five
ethanol predictions. Methanol transfer energies still differ, for an undetermined
reason. At the paper's predicted Fig. 7 solubilities, the native methanol−ethanol
$\sum\nu\ln x\varphi$ intervals for KCl and KBr include zero under both conventions.
The same cross-solvent test fails between the Fig. 7 alcohol and Fig. 9
[C4mim][TfO] predictions: the declared formulation does not reproduce the paired
model bars within the stated reading bounds (ionic-liquid side higher by about
+10 to +12 for all six pairs, under both conventions; the cause is unresolved).
Ion pairing and the solid–liquid solve are Engine capability gaps.

## Layout

| Path | Content |
|---|---|
| `inputs/parameters/components.csv`, `pairs.csv` | Paper Tables 1–5 and 8 (and Table 2 IL ions), with locators |
| `inputs/digitized/` | Digitized Figs. 4, S1, 7, 8b and 9, with QA overlays in `qa/` |
| `scripts/run.py` | Declares parameter documents, transfer observations and molalities; reads native results |
| `scripts/bulow2021_states.cpp` | Native driver for solvent activities (osmotic coefficient) and, at explicit compositions, per-ion $\ln\varphi$, $\sum\nu\ln x\varphi$, ε and density; `run.py` builds it with `tools/native-driver regression_checks` |
| `scripts/plot.py` | Figures from `results/` only |
| `scripts/digitize/` | Digitization scripts and helper configurations |
| `results/` | Engine values, comparisons, `run-identity.json`, `figures/` |
| `results/tables/`, `figures/*/source/` | Source table transcriptions and figure crops from the paper |

Run from the checkout root with the per-checkout `.venv` (README "Checks"):
`.venv/bin/python analyses/reproduction/2021-bulow/scripts/run.py`, then
`scripts/plot.py`. All evaluations are native Engine states at 298.15 K and
1 bar. The Python declares inputs and reads results. Its only arithmetic on model
values is the labelled Eq. (24) Born checks (Fig. 8b and the Fig. 9 identity
I1), the permittivity-mapping check and the unit-conversion candidates for
Fig. 8b; it also declares the IL-ion hard-sphere diameter at 298.15 K.

## Inputs replaced

The earlier template inputs (`inputs/parameters/pure/any-solvent.csv`,
`mixed/binary-interaction/*`, `user-options.json`) were removed. They carried
the Na⁺–water and K⁺–water $k_{ij}$ as the bare slope (−7.981, −4.012), a
methanol molar mass of 0.100 kg/mol, and Born diameters and a water solvation
factor that the paper does not use (Section 2.4 sets the Born diameter to σ).
They remain at commit `0218a597`.

The printed Table 3 slopes (−7.981·10⁻⁴ T and −4.012·10⁻⁴ T) are a factor ten
smaller than Held 2014's (the source the paper cites). Evaluated literally at
298.15 K they give $k_{ij}$ = 2.14 and 1.28, and shift the Na⁺ and K⁺ transfer
energies by about −62 and −39 kJ mol⁻¹ (`results/printed-slope-probe.csv`). The evaluated values combine the Held 2014 slopes with the paper's
intercepts, 0.00046485 and 0.1998222. This is an inferred correction: the
corrigendum does not touch Table 3.

## Digitization

Paper figures come from the PDF's embedded 300 ppi images
(`pdfimages -f 6 -l 9`: Fig. 7 = page 8, Figs. 8 and 9 = page 9). The corrected
Fig. 4 and Fig. S1 come from the corrigendum PDF (`pdfimages`, images 3 and 4).
Axes are calibrated on two major ticks per axis. Curves were traced with the CSE
`digitize_xy_plot.py` helper (configs in `scripts/digitize/*.json`), merged into
`fig4-curves.csv` (every third column; pixels within 8 px of a marker removed).
Markers come from `fig4_markers.py`. The bar values in Figs. 7, 8b and 9 are read
at the outline centre of each bar end by the column scans in `scripts/digitize/`.
Reading uncertainty is ±1.5 px: ±0.21 kJ mol⁻¹ (Fig. 8b), ±0.0002 in mole
fraction (Fig. 7), ±0.0001 (Fig. 9), and 0.0014–0.0036 in γ (Fig. 4). Bars thinner than
their outline are recorded as upper bounds. Fig. 8b error bars are read at both
caps. The upper cap follows the error line on the column darkest just above the
bar end, because the 1–2 px line sits off the bar centre. The comparison uses the
asymmetric [lower, upper] interval.

## Conventions, tolerances, gaps and limits

These are stated in [`notebook.qmd`](notebook.qmd) ("Basis and conventions",
"Capability gaps", "Limits and provenance"). In short:

- **Composition basis.** The salt axis is formula-unit and
  $\gamma_\mathrm{solvent}$ is on the species basis. This pairing was selected
  post hoc using the same Fig. 4 curves subsequently assessed.
- **Diameters.** Inorganic ions: σ and 0.88σ unified (hard sphere, Debye–Hückel,
  Born), both retained; neither is adopted. Water σ(T).
- **Ion–alcohol $k_{ij}$.** Zero for predictions.
- **Ionic-liquid states (Fig. 9).** `all-ion-pairs` dispersion on both solvent
  sides, the nearest expressible form of the paper's pair rule; formula basis.
- **Tolerances.** Digitization-limited, fixed before running.

Wall time: about 31 min for the first case build (sources, digitization, runs
and record). The revision (diameter conventions, native driver, error-bar caps,
record) took about 14 min by the system clock, from the start-of-round merge to
the passing checks. The native run takes about 13 s (19 s when the driver
recompiles); the Fig. 9 analysis round (#418) took about 40 min.
