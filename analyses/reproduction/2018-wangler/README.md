# Wangler MDEA/CO₂: printed neutral model and total-pressure discrepancy

**The published model, as printed, does not reproduce its own neutral water–MDEA
curve. The Engine's evaluation of that declaration is verified against FeOs at
353 K and liquid x(MDEA) = 0.2.** At 333/353/373 K, Wangler Figure 1 gives about
0.147/0.362/0.789 bar, compared with 0.121/0.282/0.634 bar from the Engine.
These departures exceed the conservative ±0.01 bar graphical allowance.
The 313 K difference remains unresolved at that allowance.

The same-pressure reaction-standard correction gives Table 10 AARDs of
**96.248872% (S0)** and **89.550063% (S1)**. The original 96.25% and 89.55%
discrepancies reproduce; the correction does not resolve them. The
[research notebook](notebook.qmd) presents the conditions, comparisons, uncertainty,
ranked explanations and limits. It has not been independently reviewed or
promoted for manuscript use.

## What is established

- FeOs and the Engine agree on the neutral bubble pressure to about
  4.1e-10 relative and on chemical-potential terms within 8.4e-9 in μ/(RT).
  Their dimensional association volumes are matched explicitly; different
  entered dimensionless volumes account for their different prefactors.
  This is numerical verification at finite neutral states, not qualification
  of the reacting ionic model or temperature derivatives.
- Pure MDEA at 353 K gives 89.09 Pa, bracketed by Kim's Antoine/Riedel values
  80.67–96.78 Pa. The AARD against Kim's seven neat measurements is 17.938912%.
  These observations are not the parameter source's original fit grid.
- Pahlavanzadeh's Table 1 vapor-pressure coefficients labeled MDEA repeat
  Kim's MEA row. The density correlation is also inconsistent with its own
  density axis when assigned MDEA's molar mass. No corrected density
  correlation is inferred. The mislabeled correlation is not used by the
  Engine's PC-SAFT calculation, so it does not itself explain the binary error.
- No tested single printed convention recovers all four Figure 1 readouts
  within ±0.01 bar. The literal Eq. 12 radical lowers the 353 K pressure to
  0.26703 bar. Sensitivities that contradict Wangler's adopted parameters
  are labeled as excluded controls; none is adopted or fitted to pressure.
- `reference_pressure_pa=None` implements the source's same-pressure standards
  (p. 18 Eqs. 20–21). It changes neither the Table 2 constants nor the pure-water
  and infinite-dilution-in-water species standards. It cannot change Figure 1,
  whose declaration has no reactions.

The full cause of the Table 10 discrepancy remains unresolved. Ranked remaining
explanations are: (1) actual source neutral parameters or curve/table provenance
that differ from the printed records; (2) an association convention not uniquely
specified by the article; (3) a separate ionic parameter/convention transfer,
including the duplicated water–MDEA/water–MDEAH⁺ Table 6 rows, unmatched osmotic
curve shape, underdetermined temperature law and approximate solvent
permittivities; and (4) a neutral Engine evaluation defect, strongly disfavored
at the cross-checked states. One author-evaluated neutral state with its actual
interaction and association inputs is the smallest remaining discriminator.
No author communication was sent.

## Basis and retained values

S0 excludes like-ion dispersion including self. S1 retains self-dispersion and
excludes same-sign cross dispersion. Both calibrate water–MDEAH⁺ k separately
to Tables 8–9 osmotic observations at 273 and 313 K with the existing relative
objective, exact native Jacobians, declared bounds and starts, fixed Table 7
association corrections and water permittivity. Their ten-point AARDs are about
2.37% and 2.16%. Linear transfer to 303 K is an assumption, not a uniquely
identified source law. The Table 10 predictions use apparent solvent mass
fraction w(MDEA) = 0.6, mole-fraction reaction constants and fixed apparent-blend
permittivity. Pressure observations are withheld from calibration.

| Evidence | Numerical owner | Producer |
|---|---|---|
| Neutral Figure 1 comparison | [neutral-figure-1.csv](results/neutral-figure-1.csv) | `neutral_engine.py`, `compare_neutral.py` |
| Engine versus FeOs | [cross-code.csv](results/cross-code.csv), [term values](results/cross-code-terms.csv), [derivative steps](results/term-derivative-steps.csv), [literal Eq. 12](results/cross-code-literal-eq12.csv) | `neutral_engine.py`, `neutral_feos.py`, `compare_neutral.py` |
| Pure MDEA | [saturation comparison](results/pure-mdea.csv), [measurement comparison](results/pure-observations.csv) | same neutral scripts |
| One-at-a-time conventions | [convention-sweep.csv](results/convention-sweep.csv) | `neutral_models.py`, `neutral_engine.py`, `compare_neutral.py` |
| Superseded one-bar standards | [reference-pressure-baseline.csv](results/reference-pressure-baseline.csv) | `reference_pressure.py` |
| Osmotic calibration | [row values](figures/figure-10/output/wangler-2018-mdeah-osmotic-plotted-data.csv), [summary](figures/figure-10/output/wangler-2018-mdeah-osmotic-summary.json) | Figure 10 `generate_data.py` |
| Corrected Table 10 predictions | [row values](figures/figure-14/output/wangler-2018-table-10-native.csv), [summary](figures/figure-14/output/wangler-2018-table-10-native-summary.json) | Figure 14 `generate_data.py` |

[Source locators](inputs/source/source-locator.md) identify journal pages,
figure calibration, source hashes, correlations and convention meanings.
`neutral-provenance.json` identifies the replay scripts and raw-state files;
`neutral-diagnostics.json` records numerical differences and limits.
CSV files own homogeneous numerical comparisons; notebook tables and figures
are generated views of these values.

The historical printed-interaction high-pressure diagnostic no longer converges
from its declared seed, including the previously tried alternative composition.
Current CSV rows preserve those failures and omit pressure predictions for them.
The converged fitted lanes still establish the Table 10 deficit. This separate
numerical limitation does not demonstrate physical nonexistence. Old printed-lane
values must not be combined with the current fitted predictions.

## Reproduction

From the repository root, prepare the checkout's environment and install its wheel
(the existing Figure 10 generator records the wheel hash):

```bash
nice -n 10 bash tools/setup.sh
nice -n 10 bash tools/build-engine.sh --install-wheel
nice -n 10 .venv/bin/python analyses/reproduction/2018-wangler/figures/figure-10/scripts/generate_data.py
nice -n 10 .venv/bin/python analyses/reproduction/2018-wangler/figures/figure-14/scripts/generate_data.py
nice -n 10 .venv/bin/python analyses/reproduction/2018-wangler/scripts/reference_pressure.py
nice -n 10 .venv/bin/python analyses/reproduction/2018-wangler/scripts/neutral_engine.py
nice -n 10 env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 ~/.cache/epcsaft-comparators/222/feos/bin/python analyses/reproduction/2018-wangler/scripts/neutral_feos.py
nice -n 10 .venv/bin/python analyses/reproduction/2018-wangler/scripts/compare_neutral.py
nice -n 10 .venv/bin/python analyses/reproduction/2018-wangler/scripts/present.py
nice -n 10 .venv/bin/python analyses/reproduction/2018-wangler/figures/figure-10/scripts/render.py
nice -n 10 .venv/bin/python analyses/reproduction/2018-wangler/figures/figure-14/scripts/render.py
nice -n 10 quarto render analyses/reproduction/2018-wangler/notebook.qmd --to html --no-execute
```

The FeOs interpreter is the existing local comparator environment, version 0.10.1;
`neutral_feos.py` also sets one native thread. Raw states and full derivative
probes go to `build/issue-857/build-r1/`. Engine calculations and fits run in
native C++; the independent comparator runs native FeOs. Finite differences
verify only term derivatives, never solve equilibrium. Presentation and notebook
rendering read retained values without model execution.

Rerunning setup restores the per-checkout editable installation. Requested check
logs and a complete evidence copy are retained in the main checkout at
`build/issue-857/build-r1/`. No Engine runtime code changes are part of this study.

## Historical records

`figures/figure-10/output/batch-5-*` retains superseded fits with other diameter
and association conventions. `figures/figure-14/output/wangler-2018-table-10-grepe-*`
retains the historical retired activity-scale projection without a verified wheel
identity. Neither enters the current notebook. Earlier values and interpretations
remain in Git history; the current tables identify the corrected declaration.
