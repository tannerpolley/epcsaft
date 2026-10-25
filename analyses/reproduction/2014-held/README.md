# Held 2014, ePC-SAFT revised

Source: Held, Reschke, Mohammad, Luza and Sadowski (2014), *ePC-SAFT revised*,
Chem. Eng. Res. Des. 92, DOI 10.1016/j.cherd.2014.05.017 (Zotero
3DIW87Q5/TD5T732H); identity in `inputs/source/source-manifest.csv`.
Issue #387. The record is [`notebook.qmd`](notebook.qmd); every number there is
computed from the retained files below.

## Outcome

- **Reproduced with the Engine (numerical verification against the digitized
  published lines):** Figs. 1 and 2 (NaCl, KBr, KCl; strategies 1 and 2;
  298.15 and 273.15 K) for 7 of 8 series within 0.005 in Φ; the exception is
  Fig. 2 NaCl at 298.15 K (maximum 0.007361 at a data marker, 95th
  percentile 0.001526). Fig. 3b NaAc. Fig. 6a (1-butanol + NH4Cl + water
  LLE) over the traced 0.2447–12.8874 wt% NH4Cl in the aqueous phase.
  Fig. 3b KAc and LiAc exceed the strict maximum only where the traced line
  crosses data markers or the curves merge.
- **Exploratory:** Fig. 5 (benzene solubility). Its 5 % bound was set after a
  five-point probe; NaI's maximum difference from the line, 2.354 %, is an
  observation.
- **Table 5 osmotic ARD (indirect Hamer–Wu comparison; the paper used
  Lobo–Quaresma data, not retained):** with strategy 2, the Engine ARD
  against Hamer–Wu 1972 lies within 0.5 percentage point of the printed ARD
  for most of the univalent-cation salts with retained data. The outliers
  are listed in the notebook.
- **Inferred convention:** the paper prints σ for the Debye–Hückel term.
  Under the declared parameters and the supplemental permittivity, the
  Eq. (21) diameter 0.88 σ agrees substantially better with the digitized
  curves: 7 of 8 Fig. 1–2 series lie within 0.005 with 0.88 σ, against
  series maxima of 0.016–0.019 with σ. This is inferred from the numbers,
  not proven; the 2021-bulow and 2020-bulow reproductions point the same
  way. [`analyses/2026-held-2014-chlorides`](../../2026-held-2014-chlorides/README.md)
  used σ, so its divalent verdict should be re-run with 0.88 σ before #209
  relies on it. It is not edited here.
- **Water σ(T):** Held 2014 Table 1 footnote a prints +0.01146 T in the last
  exponent (σ_w = −40.3 Å at 298.15 K). The −0.01146 T used here (2.797 Å)
  is a source-backed correction from Held 2008 Table 1 and Eq. (7), not
  verbatim Held 2014 input.
- **Not reproduced:** Fig. 3a (strategy-1 acetates): the printed Held and
  Sadowski 2009 Ac⁻ parameters give Φ far above the traced lines (KAc 1.737
  against the traced 0.995 at 3.0 mol/kg). Their Engine γ± deviations are of the same order as that
  paper's own, but a γ± check does not verify a Φ curve; whether parameters
  or figure are in error is unresolved. Fig. 5 NaCl: the Engine predicts
  less salting-out of benzene (7.351 against 4.688 mmol/kg at 3.5 mol/kg);
  the two Cl⁻–benzene k_ij probes (0.45, 0.5) do not remove it, and other
  coefficients were not tested.
- **Gaps:** Fig. 4a has no native result: the Engine's native osmotic
  coefficient admits one pure solvent and one salt, not amino acid + KCl.
  Its Python contraction is kept only as a diagnostic (capability gap for
  #387). Fig. 6b (ternary diagram) was not digitized or compared. Fig. 4b and
  Eq. (6) need melting properties that the paper does not print and an SLE
  capability that the Engine lacks. Table 5 densities, the divalent-anion
  entries, and toluene/NaBr/Na2SO4 in Table 6 have no retained data.

## Reuse

The chloride analysis remains the owner of the LiCl/NaCl/MgCl2/CaCl2 γ± and Φ
row qualification. Its Hamer–Wu LiCl and NaCl rows (`inputs/evaluated.csv`)
are read for Table 5; no σ-convention chloride result is recomputed.

## Files

- `inputs/held-2014-ions.csv`, `held-2014-table-1.csv`, `-table-3.csv`,
  `-table-4.csv`: Tables 1–4 transcribed from the PDF page images. The water
  record, model families and topology come from
  `data/packets/held-2014-strategy2/1`; `scripts/model.py` assembles each
  document; with σ closest approach the KCl document reproduces that
  packet's Φ exactly.
- `inputs/held-2014-table-5-osmotic-ard.csv`: printed Table 5 osmotic ARDs.
- `inputs/hamer-wu-1972-osmotic.csv`: Hamer–Wu 1972 Φ and γ± for 37 salts,
  extracted by `scripts/transcribe_hamer_wu.py` from the Zotero Markdown child
  VYE43TLN of 28NPJ5H6 (SHA-256 recorded by the script). Its γ± values
  match all 105 overlapping rows of `data/packets/hamer-wu-1972-aqueous-alkali-halides`.
- `inputs/digitization.json`, `paper-lines.csv`, `paper-markers.csv`: the
  published lines and data markers traced by `scripts/digitize_figures.py`
  from the PDF's embedded 300 ppi images (images and QA overlays are kept
  outside the repository).
- `scripts/run.py` (Table 5, Figs. 1–3), `scripts/ternary.py` (Figs. 4a, 5,
  6a), `scripts/probes.py` (retained Fig. 3a γ± and Fig. 5 k_ij probes),
  `scripts/plot.py`; results in `results/`, figures in `figures/`.
- The older `figures/figure-0*/source/*.csv` files predate this
  reproduction, have no recorded lineage and are not used. For example, the
  Fig. 3 "digitized" NaAc value at 3.0 mol/kg (1.225) does not match the
  plotted line (about 1.19). The former `inputs/parameters/` template, which
  stored footnote slopes as Na⁺/K⁺–water k_ij values, was replaced by the
  table transcriptions.

## Replay

```bash
.venv/bin/python analyses/reproduction/2014-held/scripts/run.py
.venv/bin/python analyses/reproduction/2014-held/scripts/ternary.py
.venv/bin/python analyses/reproduction/2014-held/scripts/probes.py
.venv/bin/python analyses/reproduction/2014-held/scripts/plot.py
```

All equilibrium, density, fugacity, reference and osmotic-coefficient
calculations are native; Python declares inputs and reads results. The one
exception is the Fig. 4a diagnostic contraction above. The full replay takes under half a
minute (the exact timing is in `results/run-identity.json`). Not
independently reviewed or promoted.
