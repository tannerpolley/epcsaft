Source: Ascani (2023), chemical and phase equilibria with an esterification reaction; identity in `inputs/source/source-manifest.csv`.

Reproduced natively in the Engine for #387; the numbers and their interpretation are in `notebook.qmd`,
computed from `results/`. Not independently reviewed.

## What is compared

- **Paper model value (exploratory consistency).** Ascani publishes model results only as three-dimensional
  renderings (Figures 8 and 9), so the one numerical model value is K_a (Table 4). It is recomputed at each
  Senina homogeneous composition from the Engine's single-liquid reaction row (`CompiledSystem.evaluate`).
  The parameter-mapping evidence is the direct check of the declared parameters against Tables 1 and 2.
- **Experiment (physical comparison, literature parameters).** Senina et al. 2022 (Ascani ref. [72]) Tables 2
  and 3: 15 homogeneous chemical-equilibrium compositions and 9 reactive tie-lines of system 1 at 318.15 K,
  in `inputs/data/`. Each homogeneous point is a feed reacted to the model equilibrium; each tie-line is fed
  at its midpoint. Every feed is one native solve, with phase discovery, of chemical and phase equilibrium together.
- **System 2 (prediction only).** Acid + alcohol feeds without products in both systems. Schmitt's 2006
  thesis (ref. [97]) holds the system-2 compositions and is not in the library.
- **Not reproduced.** Figures 2 and 3 (the paper gives no temperature or pressure); Figure 1 (ideal
  solution) and Figures 4 to 7 (algorithm schematics) carry no model values.

## Conventions

- K_a is the Eq. 10 activity product with pure-liquid references at the system T and p, constant per system.
  It is declared as a constant `ReactionLogPolynomial` on `RawMoleFractionInfiniteDilution`, with every species on
  `SpeciesReference.PureLiquid`. Since the sum of nu is 0, no basis term arises.
- k_ij(T) follows Eq. A8, evaluated at the system temperature. Water sigma(T) is the Table 1 footnote correlation.
- Cross association uses Eq. A6 (arithmetic energy) between all two-site components. The esters have zero
  self energy and associate only through cross edges. Eq. A7 prints the size bracket squared; the primary
  `printed` documents declare it as explicit edges, and the `wolbach-sandler` documents use the Engine's
  cubed rule as the variant.
- The catalyst is omitted, as in the paper.
- `inputs/parameters/pure/` and `inputs/parameters/mixed/` are earlier transcriptions not read by the run.
  `scripts/parameters.py` reads the source tables under `results/tables/`.

## Running

```
python scripts/parameters.py   # writes inputs/parameters/system-*.json
python scripts/run.py          # native solves; writes results/*.csv, timing.csv, provenance.json
python scripts/render.py       # figures/system-1-reacting-lle.png from the retained CSVs
```

## Engine findings (fixed after this record)

Near the gap edges, phase discovery disagrees with a declared two-liquid solve that has lower Gibbs energy,
in `results/discovery-cases.csv`. Both defects were fixed under #390 (PR #408) after this record, which has not
been rerun. In system 2 at acid feed fractions 0.75
to 0.78, discovery fails with "resolved descent without an admitted lower-G declaration". In system 1 at 0.66,
it returns one liquid where a small but finite second liquid has lower Gibbs energy. The sweep's one failed
grid feed is solved with declared liquids and labelled in its `route` column.

## Tolerances stated before the reacting runs (#387)

- Native solve: Engine row residual at most 1e-10. Each returned liquid is checked independently at its
  returned composition and density: reaction residual |ln K_a(x) - ln K_a| from the Engine's single-liquid
  reaction row and the cross-phase log-fugacity difference each at most 1e-8, and the pressure at the returned
  density within 1e-8 relative of 1 bar. Added in revision: the reaction row must return ln K_a within 1e-9 at
  every converged single-liquid state.
- K_a (exploratory consistency with the reported value, not a numerical-verification pass): the paper names no
  source point. The test asks whether the Table 4 value lies inside the Engine K_a band of at least one Senina
  homogeneous point when every printed mole fraction is moved by its half rounding unit, 0.0005. It was chosen
  after the K_a values at the printed compositions had been seen, so it cannot verify the parameter mapping.
- Experiment (physical comparison with literature parameters): Ascani publishes no deviation metric and no
  numerical model values, so the metric is the source data paper's (Senina 2022) standard deviation
  sigma = 100 sqrt(sum dx^2 / (N n)) over all components and points (per phase for tie-lines). Ascani's
  "quantitatively correct" claim for homogeneous chemical equilibrium is read as every point within 0.01 mole
  fraction (five times the data's u(x) = 0.002); the tie-line claim is qualitative and is reported, not gated.
