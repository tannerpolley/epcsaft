# Acetic acid + water + toluene: dissociation with liquid–liquid equilibrium (#374)

Source: Ascani, Pabsch, Klinksiek, Gajardo-Parra, Sadowski and Held, *Chem.
Commun.* 58 (2022) 8436, doi:10.1039/D2CC02943J (Zotero AGPHR9MB), and its
ESI (CC BY 3.0), sha256 `84c3c8da4592a25f2eef1bddc5c3712096c31463c5df88e4b4fb6f1df6e9ede1`,
fetched 2026-10-08 into `build/issue-374/esi/` (untracked). ESI Tables S11
and S12 are transcribed in `inputs/`.

## Tolerances (stated before running)

- **pH against their ePC-SAFT advanced values (numerical verification).**
  Agreement if |ΔpH| ≤ 0.005. The ESI reports three decimals (±0.0005), but its
  own rows at the identical x = 0.006780777 give 2.566 and 2.569, so their
  model states were not a function of the tabulated x alone.
- **Partition coefficient K = x_aq/x_org of acetic acid (numerical
  verification).** Agreement if the relative difference is ≤ 1 %; the ESI's
  duplicate-x rows differ by 0.15 %.
- **Experimental pH (physical comparison with literature parameters).** No
  pass/fail; report AARD and bias beside the ESI's 6.902 % AARD (Table S6).
- **Solution quality.** Requested row tolerance 1e-10 met; per-phase charge
  |Σ z_i n_i| ≤ 1e-12 mol; density-basis log-affinity of both reactions (acid
  dissociation and water autoprotolysis), recomputed from public state
  potentials against the solve's compiled targets, ≤ 1e-8; reduced Hessian
  minimum > 0 (local stability only); reference terminal change ≤ 5e-5 (the
  Engine's `reference_limit_terminal_pair` acceptance rule, a convergence
  diagnostic, not an error bound). The molality-basis acid law with φ∞ from a
  separate public reference path is reported only as a reference-consistency
  diagnostic.
- **Comparison at the ESI compositions.** pH and K are interpolated in
  ln x_aq on a native feed sweep. The coarse/fine interpolation differences
  (81 against 161 feeds) must be ≤ 1e-4 pH and 1e-4 relative in K; they are
  not an error bound.

## Run

From the repository root, with the checkout's `.venv` (README "Checks"):

```bash
.venv/bin/python analyses/2026-ascani-carboxylic-acid-lle/scripts/parameters.py              # ESI Tables S1-S5 -> inputs/parameters-*.json
.venv/bin/python analyses/2026-ascani-carboxylic-acid-lle/scripts/run.py 161     # native solves -> results/*.csv
.venv/bin/python analyses/2026-ascani-carboxylic-acid-lle/scripts/conventions.py # water density, no cross-association -> results/conventions.json
.venv/bin/python analyses/2026-ascani-carboxylic-acid-lle/scripts/render.py                  # figure and results/summary.csv, no Engine calls
quarto render analyses/2026-ascani-carboxylic-acid-lle/notebook.qmd                          # HTML view (untracked)
```

Results and interpretation: [notebook.qmd](notebook.qmd).
