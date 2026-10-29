# Measured NaCl solubility in water

The declared Held packet `/2` plus Malmberg–Maryott dielectric variant predicts a
falling saturation molality while the Pinho–Macedo measurements rise with temperature.
Every prediction interval propagated from one standard error of the dissolution
constant is disjoint from its measured 95% expanded interval.
The [notebook](notebook.qmd) separates that physical discrepancy from the native
saturation, conservation and reference-state closure checks.
[Retained comparisons](results/comparison.csv) own the individual errors and intervals;
[summary values](results/summary.json) own the aggregate errors and attempted counts.

The retained temperature-trend diagnosis identifies the sodium–water interaction
law as the dominant positive thermal contribution, opposed by the net of all
other terms. The notebook preserves both decompositions, evaluated activity and
enthalpy comparisons, rounding bounds and four discriminating sensitivities.
The local frozen-interaction comparison reverses the saturation slope at the same
central states. A global freeze reverses the endpoint trend but overshoots the
measured hot solubility; it is a sensitivity, not a source-backed correction.
Replacement temperature laws remain underdetermined.
No parameter was fitted or selected using these measured solubilities.

The inputs are the eight Pinho–Macedo Table 1 NaCl–water temperatures,
298.15–353.15 K, at the NIST ThermoML pressure of 101300 Pa.
Clarke–Glew Table 34 supplies each temperature's dissolution Gibbs energy,
enthalpy, heat capacity and their separate standard errors at 101325 Pa.
The native logarithmic dissolution law uses the molality, infinite-dilution
standard, with standard molality 1 mol/kg.
Each solve declares 1 kg initial water, 12 mol NaCl feed, one liquid and halite.
The 298.15 K doubled-feed check uses 24 mol.

The [source notes](inputs/sources.md), [thermochemical rows](inputs/clarke-glew-1985-table34.csv),
[measured rows](inputs/pinho-macedo-2005.csv) and [source hashes](inputs/source-hashes.json)
retain the evidence and uncertainty conventions.
The two `inputs/held-*-malmberg-maryott.json` files are explicit derived input
snapshots, leaving the reusable packets unchanged.
Only packet `/2` is the primary prediction; packet `/1` appears as a single
298.15 K sensitivity row.
The constant-permittivity results are a separate control of the temperature-dependent
solvent property, not another fitted prediction.

All concentrations extrapolate beyond Held's up-to-5 mol/kg density/osmotic
activity-fit range.
The ion parameters extrapolate in temperature above 298.15 K; the published
sodium–water law was motivated by osmotic behavior at 273.15 and 298.15 K.
The result is conditional on the one-liquid/halite topology.
No phase discovery with solids, hydrate search or other solid comparison was performed.
Positive pressure slope establishes only mechanical regularity of the declared liquid.
Malmberg–Maryott is an explicit dielectric variant and is not proven to match
Held's fitting relation.
Ion-parameter covariance and model-form uncertainty remain unquantified.
The predicted interval is a one-standard-error thermochemical input sensitivity,
not a complete prediction interval or a 95% confidence interval.

With this checkout's own environment prepared by README "Checks", reproduce
calculation and presentation separately from the study directory:

```bash
cd analyses/2026-nacl-solubility
export EPCSAFT_BUILD_JOBS=6
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
nice -n 10 ../../.venv/bin/python scripts/run.py
nice -n 10 ../../.venv/bin/python scripts/probe.py
nice -n 10 ../../.venv/bin/python scripts/decompose.py
nice -n 10 ../../.venv/bin/python scripts/reference.py
nice -n 10 ../../.venv/bin/python scripts/verify.py
nice -n 10 ../../.venv/bin/python scripts/present.py
```

Finish builds and numerical reruns before 02:45 America/Denver or start them
after 05:00, outside the quiet timing campaign. The diagnostic producers also
refuse to start a bounded group during 03:00–05:00.

Python supplies inputs and summarizes returned values; native C++ owns equilibrium,
liquid properties, infinite dilution and exact derivative actions.
`reference.py` converts the published evaluated coefficients by direct algebra;
it does not solve or fit an equation of state.
The notebook only reads retained presentation includes and figures.
`present.py` prepares both the original solubility comparison and temperature-trend
presentation from retained results without importing the Engine.

The diagnostic scripts share `load_engine.py`, which checks that the loaded
extension belongs to this checkout. The source rows and locators are in
[diagnostic source inputs](inputs/source/diagnosis-sources.json), including
PPB Table A-8 p.53, Table A-3 p.26, Eq. (27) p.5, and Ascani–Held (2021) Table 5.
[Diagnostic results](results/temperature-trend/saturation.csv) and the
[diagnostic identity](results/temperature-trend/run-identity.json) remain separate
from the original solubility results and `results/run-identity.json`.
The diagnostic identity records the current commit dynamically, the Engine source
tree, loaded extension/library hashes and source/script hashes.

The notebook also compares the declared model with Farelo's single NaCl–KCl
co-saturation state. NaCl's sampled calculated range is disjoint from both
conditional measured intervals; KCl overlaps. This extends the existing study
with the same parameter and thermochemistry owners. The shared input is
augmented at declaration with Held's published potassium–water temperature law.
The [source inputs](inputs/co-saturation.json),
[seven native solves](results/co-saturation/native-states.json), and
[interval summary](results/co-saturation/summary.json) retain the calculation.
Reproduce from this directory with:

```bash
nice -n 10 env EPCSAFT_BUILD_JOBS=6 ../../.venv/bin/python scripts/co_saturation.py
nice -n 10 ../../.venv/bin/python scripts/present.py
```

For the 28 October quiet campaign, finish builds and numerical runs before
02:30 America/Denver or start them after 05:30. The one-state comparison has
unspecified measurement confidence, an interpolated KCl dissolution constant,
an assumed pressure, sampled rather than proved input extrema, and the known
binary NaCl temperature-trend error. It does not establish physical validation.

The separate [fixed Figiel halite comparison](results/figiel-902/comparison.md)
uses `figiel-2025-chloride-mixture-refit/2` at room temperature, without fitting.
Both pressure-specific predictions pass the three-percent engineering screen
while exceeding their source uncertainty intervals. All native attempts and
raw residuals are retained under `results/figiel-902/`; the notebook includes
separate source comparisons and thermochemical/temperature sensitivities.
Reproduce this comparison only, from the repository root:

```bash
EPCSAFT_BUILD_JOBS=6 nice -n 10 .venv/bin/python analyses/2026-nacl-solubility/scripts/figiel_halite.py
EPCSAFT_BUILD_JOBS=6 nice -n 10 .venv/bin/python analyses/2026-nacl-solubility/scripts/figiel_present.py
EPCSAFT_BUILD_JOBS=6 nice -n 10 quarto render analyses/2026-nacl-solubility/notebook.qmd
```

For this comparison, finish compiles and solves before 02:15 America/Denver;
no compiles or solves are allowed during the 02:30–08:00 host campaign.
The producer refuses to start at or after the cutoff until 08:00.
The historical commands above retain their original study scope.
