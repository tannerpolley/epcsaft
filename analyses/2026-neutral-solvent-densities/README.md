# Fixed-input neutral-solvent densities and MIBK vapor pressure

This study delivers issue #916, the bounded part of #902: fixed TOP and MIBK
liquid densities and their separate sensitivities, descriptive MIBK vapor
comparisons, TOP parameter read-back, and pure n-dodecane stand-in predictions.
Nothing is fitted. The engineering screens are one percent for density and
five percent for vapor pressure, separate from measurement uncertainty.

The [notebook](notebook.qmd) reads retained native results. The
[input notes](inputs/README.md) identify source roles, units and assumptions.
The [screens](results/screens.tsv), [density rows](results/densities.tsv),
[vapor attempts](results/vapors.tsv), [raw residuals](results/residuals.tsv),
[parameter read-back](results/readback.tsv) and [closure maxima](results/closure-maxima.tsv)
retain the numerical evidence. All source comparisons use fixed parameters;
source-scale disagreement and failed numerical attempts remain visible.

From this checkout's root, prepare its own editable environment and reproduce:

```bash
EPCSAFT_BUILD_JOBS=6 nice -n 10 bash tools/setup.sh
EPCSAFT_BUILD_JOBS=6 nice -n 10 .venv/bin/python analyses/2026-neutral-solvent-densities/scripts/run.py
EPCSAFT_BUILD_JOBS=6 nice -n 10 .venv/bin/python analyses/2026-neutral-solvent-densities/scripts/present.py
EPCSAFT_BUILD_JOBS=6 nice -n 10 quarto render analyses/2026-neutral-solvent-densities/notebook.qmd
```

The native driver calls the shared parameter builder, density observation owner
and equilibrium driver; Python declares paths, builds the program, reads files
and presents retained values. The run records this checkout's loaded extension,
shared native library, executable and source hashes. The producer refuses to
start compiles or solves from 02:15 through 08:00 America/Denver. The lead's
host timing campaign begins at 02:30. Check commands must also finish before
02:15. No Engine runtime changes or new numerical methods are included.

The figures show observations and predictions at sampled temperatures as
points, with no interpolated model curves. Density error bars retain their
separate source meanings; vapor uncertainty is unavailable. Independent MIBK
density comparison, unknown-independence Fang comparisons, possibly fitted
vapor rows, report-only Okeme estimates and assumed-carrier predictions are
kept apart. Every carrier result is dependent on the n-dodecane stand-in.
Independent delivered review and input adoption remain with the lead.
