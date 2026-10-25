Source: Bülow, Ascani and Held (2020), ePC-SAFT advanced, part I, with its corrigendum (Fluid Phase Equilib. 548 (2021) 113184); identity in `inputs/source/source-manifest.csv`.

# Engine reproduction (issue #387)

**Result.** The Engine reproduces the paper natively at 298.15 K and 1 bar with the published parameters, once its conventions are mapped. Numbers and tables are in [`notebook.qmd`](notebook.qmd).

- **Figs. 2–5, infinite dilution.** Hydration and transfer Gibbs energies of seven alkali/halide ions in water, methanol and ethanol match the paper's bars. Hydration agrees within about 1.4 kJ/mol for the advanced model and 3.8 kJ/mol for revised; transfer agrees within 0.25 kJ/mol for both. The Born part matches Fig. 3.
- **Figs. 6–7, corrigendum, activity coefficients.**
  - The revised curves and the Fig. 6 right-panel contributions match the Engine directly.
  - The advanced curves of Fig. 6 (left) and Fig. 7 equal the Engine's mole-fraction-scale γ±\* times (1 + ν m M_s). The same figure's right panel shows ln γ±\* without that factor, so the corrigendum is internally inconsistent by ln(1 + ν m M_s) within digitization precision. The cause in the authors' plotting is inferred, not reproduced.
- **Against experiment, on the basis consistent with the plotted data.**
  - The advanced model predicts NaBr in ethanol well.
  - It underpredicts the strong rise of γ±\* for LiBr, LiCl and NaBr in methanol.
  - It captures the measured osmotic coefficients of LiCl and LiBr in methanol (Zafarani-Moattar 1998) far better than revised.

Evidence type: numerical verification against the paper's digitized model values, and physical comparison with literature data using unfitted literature parameters. Not independently reviewed or promoted.

## Conventions mapped onto Engine terms

Only the hard-chain ion diameter is printed as 0.88 σ_i. Using 0.88 σ_i for Born and Debye–Hückel too is an inferred convention: no source prints that factor for those terms.

| Ion diameter | Printed in the sources | Engine | Status and evidence |
|---|---|---|---|
| Hard chain | d_i = 0.88 σ_i (Held 2014 Eq. 21) | `packing_diameter` = 0.88 σ_i | printed |
| Born radius a_i, Eq. (2) | the "ePC-SAFT ion diameters σ_i" (Bülow 2020 Secs. 2.1, 3.2) | `born_diameter` = 0.88 σ_i; Born `c_shell = c_dielectric = 0` | inferred: with d = σ_i, Li⁺ hydration moves from −555.0 to −489.3 kJ/mol (paper −555.7; tolerance 10) |
| Debye–Hückel closest approach | σ_i (Held 2014 Sec. 2.1); "unified ion diameters" (Bülow 2020 Sec. 3.4) | `debye_huckel_diameter` = 0.88 σ_i | inferred: with d = σ_i, revised-curve median errors are 3.85–10.15 px against 0.55–1.38 px |

| Paper | Engine |
|---|---|
| ε_r(x̄) linear in mole fraction, Table 4 | `component-permittivity-mixing`; water 78.09, methanol 33.05, ethanol 24.88, ions 8 |
| ePC-SAFT revised | `fully-dissociated-debye-huckel` with `solvent-only` permittivity |

Source-packet corrections:

- Fig. 6 and Fig. 7 top-left are LiBr in **methanol** per the corrigendum; the packet file `figures/figure-06/source/ethanol-libr.csv` keeps the old label.
- The Fig. 5 caption's ethanol/methanol panel labels are swapped. The Born transfer term fixes which is which, and the packet's `water-*-contributions.csv` names inherit the caption.
- `inputs/parameters/pure/any-solvent.csv`, carried over from a template, has wrong alcohol permittivities (1) and molar masses (0.1 kg/mol), and Figiel 2025 Born diameters. It is not used. `scripts/parameters.py` builds the document from Tables 1–4.

## Limits

- The paper's values are digitized, so agreement is bounded by read precision. Tolerances were fixed before running.
- Per-term ln φ (Fig. 3/5 hard-chain, dispersion and association bars) is not reported by the Engine; only the Born term is separated.
- Digitized data discs that sit on a curve are partly missed. Curve columns touching any disc ink are excluded by one rule for every curve (`scripts/digitize.py`); after it, every screened curve is within 3 px.
- Only 298.15 K and 1 bar.

## Replay

From the repository root, with the per-checkout `.venv`:

```bash
.venv/bin/python analyses/reproduction/2020-bulow/scripts/digitize.py <dir with `pdfimages -j` images of Zotero YHZIPCDI (parent L79J5R52; byte-identical to CZQJCDPD under the corrigendum parent 9XTQJ9P7)>
.venv/bin/python analyses/reproduction/2020-bulow/scripts/parameters.py
.venv/bin/python analyses/reproduction/2020-bulow/scripts/run.py
.venv/bin/python analyses/reproduction/2020-bulow/scripts/plot.py
```

`results/run-identity.json` records the Engine commit, hashes and the calculation wall time (about 15 s). The whole reproduction (source reading, digitizing, convention resolution, runs and record) took about 26 minutes of wall-clock time on 2026-10-09; the review-fix round (data-disc rule, diameter convention table, Zotero pairing) took about 10 minutes.

Earlier README text and template or copied files removed by #341 remain at record commit `21d7a112`. Figure 7 reuses `figures/figure-06/source/ethanol-libr.csv`; #341 removed its byte-identical copy.
