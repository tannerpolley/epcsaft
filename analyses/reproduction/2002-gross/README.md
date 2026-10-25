The native Engine returns every retained declaration using the published
associating PC-SAFT parameters and existing native equilibrium interfaces.
From the default start only three water/1-pentanol dew declarations fail.
Explicit neighbouring starts recover them without changing Engine source, parameters,
tolerances or local checks. The notebook reports both coverage measures.

Figure 4's pressure-coordinate extraction and two Figure 3 bubble coordinates
were corrected against the PDF's solid PC-SAFT curves. Experimental markers
embedded in Figure 2 and Figure 4 source-curve files now retain their correct
source role. Selected misses receive bounded composition-uncertainty probes;
remaining misses stay explicitly undiagnosed. Figure 2 marker-overlap proxies
were checked against the solid curve: readable points were re-extracted and
hidden model rows omitted while preserving experimental observations. The
coincident pentanol-model/nonanol-marker coordinates in Figure 1 were also
checked; the obscured pentanol curve row was omitted. Figure-local correction
CSVs and calibrated QA overlays retain every old value and omission reason. This revision awaits independent
scientific review by the lead.

Issue #832 selects those final corrected rows without another calculation.
The notebook displays the old and corrected ordinates and omission reasons,
and computes experimental-marker statistics separately for each compound and
phase. The retained mixed-branch methanol entry in `comparison-summary.csv`
is historical arithmetic, not a liquid-density error. The selected liquid
subset contains six markers and has density AAD 5.365010200647692%; the vapor
marker remains a separate comparison. Four water/pentanol VLLE coordinates
describe one physical coexistence state. Supplemental data and the DECHEMA
correlation retain their separate source roles and are not independent validation.

The original pure-component fitting datasets are absent. Published vapor-pressure
and liquid-property AAD values cannot be recomputed from the retained summaries.
The available liquid-density marker subsets and mixture experimental tables
are compared explicitly; they do not replace the missing full fitting samples.
The [notebook](notebook.qmd) computes every reported value from retained CSV/JSON
and records numerical verification separately from physical comparison.

Source: Gross and Sadowski (2002), *Application of the Perturbed-Chain SAFT
Equation of State to Associating Systems*, [DOI](https://doi.org/10.1021/ie010954d).
The Zotero parent and PDF attachment identities and PDF hash are in
[the calculation record](results/data/calculation.json). Source tables remain
under `results/tables/` and source figure values under `figures/`; their
[inventory](results/tables/tables-manifest.csv) preserves the original locators.

The association model maps to one A and one B site per associating molecule,
unit multiplicities, reciprocal unlike-site edges only, and the native
Wolbach–Sandler rule for cross-association. Water and acetic acid use the same
published two-site topology. Figure 2 uses methanol/isobutane, as its caption
and cited experimental paper specify; the Table 2 isobutanol label and the
historical matrix's unrelated alcohol assignment are documented in the notebook.
No parameters were fitted, no Engine source changed, and no solver tolerance
was relaxed.

Retained results:

- [Default-start results](results/data/default-start-comparisons.csv) and
  [the exact failure declarations and recovering starts](results/data/default-start-failures.csv).
- [Reviewed results before revision](results/data/reviewed-state-comparisons.csv).
- [Bounded miss diagnosis](results/data/model-miss-diagnosis.csv) and
  [native composition-coordinate bounds](results/data/input-coordinate-bounds.csv).
- [Figure 8 experimental and correlation markers](figures/figure-08/source/experimental-points.csv)
  and [Figure 10 vapor–liquid markers](figures/figure-10/source/experimental-points.csv),
  with calibration records and inspected QA overlays beside each CSV.
- [Every declared state and pre-run tolerance](results/data/declared-states.csv).
- [Engine values, deviations and numerical diagnostics](results/data/state-comparisons.csv).
- [Comparisons by figure, components and phase specification](results/data/comparison-summary.csv).
- [Published pure AADs, available marker deviations and pure predictions](results/data/pure-component-summary.csv).
- [Initial declaration failures](results/data/initial-state-comparisons.csv) and
  [initial timing](results/data/initial-calculation.json), retained for the
  corrected liquid-feed and three-phase-start diagnosis.
- [Model and experimental comparison figure](results/rendered/model-and-experiment.png)
  and [its exact plotted values](results/data/plotted-states.csv).

Reproduce from the checkout root:

```bash
.venv/bin/python analyses/reproduction/2002-gross/scripts/calculate.py
.venv/bin/python analyses/reproduction/2002-gross/scripts/compare_coordinate_bounds.py
.venv/bin/python analyses/reproduction/2002-gross/scripts/present.py
python3 tools/check-notebook-numbers.py analyses/reproduction/2002-gross/notebook.qmd
bash tools/check-affected.sh origin/main
```

`calculate.py` declares inputs and uses all available CPU workers; every
thermodynamic calculation and root search runs natively. The calculation record retains default and recovery batch wall-clock times;
`coordinate-diagnosis.json` retains the bounded probe time, `revision.json`
the earlier revision work, and `proxy-revision.json` the marker-overlap revision
wall-clock time. Batch timing excludes setup and presentation.
`present.py` and the notebook read results without importing the Engine.
Experimental source transcription, digitization uncertainty, unresolved model-curve
misses, and supplemental liquid-only data above the paper's
temperature range limit the claims. Independent review, investigator promotion
and remote delivery remain with the lead. Earlier records removed by #341
remain at commit `21d7a112`.
