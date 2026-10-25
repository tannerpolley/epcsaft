# Gross and Sadowski (2001): partial PC-SAFT reproduction

The published parameters produce a native saturation state for every retained compound and useful comparisons with digitized pure liquid densities, binary vapor-liquid boundaries, and ternary phase compositions.
Retracing every sampled model point on the solid PC-SAFT curves removes the large source-selection errors; a few binary points remain outside the fixed figure tolerances.
The sampled carbon-dioxide ternary compositions agree within the declared composition tolerance; the lower-pressure heptane–butane–ethane comparison remains unresolved.
The original pure-component pressure and molar-volume AAD values and the binary equilibrium-ratio AAD values remain unreproduced because their underlying experimental rows are absent.

[The research notebook](notebook.qmd) computes its displayed values from retained CSV and JSON files and states the comparison basis, tolerances, numerical stopping points, uncertainty, and coverage.
It distinguishes numerical verification against paper model curves from physical comparisons against digitized experimental symbols with literature parameters.
These results have not been independently reviewed or promoted for manuscript use.

The selected solid-PC-SAFT extraction is retained in
[model-point-review.csv](inputs/digitized/model-point-review.csv), including old
and revised pixels. The current-main rerun returns 200 states and numerically
qualifies 197: 119 of 122 figure comparisons and all 78 reference-free predictions.
Cold propylene, diethyl ether and toluene remain excluded under the unchanged
relative pressure-closure criterion. The catalogue establishes calculation
breadth; it supplies no measured references at those temperatures.

The [source identity](inputs/source/source-manifest.csv) points to [the original paper](https://doi.org/10.1021/ie0003887).
The Zotero parent is `KXZN4KLQ`; the original PDF attachment is `ZAETGZK9`, read using `pdftotext -layout`.
The parent's default PDF attachment is a later reply and was not used.
The original source tables and figure images are preserved.

The Engine uses the neutral PC-SAFT formulation without association or polar terms, including for carbon dioxide, as the paper specifies.
Parameters and molar masses are taken verbatim.
Figure-specific interaction parameters override the table-based binary matrix where the caption gives another value.
All equation-of-state evaluations, root finds, phase splits, and nonlinear solves run natively; Python declares inputs, calibrates source images, computes reporting statistics, and renders retained values.
The current run uses main `51ce7b936` with the superancillary default start.
No Engine source or numerical criteria changed and no parameters were fitted.
The original run and production CSV files are retained under
[results/historical/](results/historical/run.json), along with the earlier
bounded probes and statistics. The historical run returned 200 and qualified 198;
its 120 qualified figure rows now become 119. Toluene is the additional exclusion.
[Current counts](results/rerun-counts.csv) and [changed existing fields](results/rerun-changes.csv)
compare the preserved original population with this run.

Toluene's current relative pressure mismatch is about six times the fixed
criterion. The #847 diagnosis traces its first qualification loss to #543's
tape optimization, which moved the liquid density by one represented double
step; the later #648 pressure start gives the current value. Liquid-pressure
cancellation at very low saturation pressure makes relative closure sensitive
to represented density steps. A returned state remains distinct from a
pressure-qualified result.

The [current bounded diagnostics](results/diagnosis.csv) use each row's own
pressure and phase densities for explicit starts. The default initialization
trace is retained for all three excluded rows in
[engine-defect-cases.csv](results/engine-defect-cases.csv); the obsolete
mechanical-loop seed event is no longer required.

The source-pixel selections, axis calibration, uncertainty allowances, and extraction overlays are under `inputs/digitized/` and the individual figure directories.
The extraction is sparse and manual. Every model point was checked in enlarged source crops and again in the overlays by the same investigator; the previous and revised coordinates are retained in `inputs/digitized/model-point-review.csv`. This is not a second-person independent digitization.
The ternary tie-line pairing and inferred mole-fraction basis remain limits on its interpretation. The printed Figure 11 pressure is retained; [native pressure probes](results/source-diagnosis.csv) support an unresolved source inconsistency and a possible caption error, which remains an inference without the Mehra–Thodos observations.
The native starts are declared in [the phase-start inputs](inputs/phase-starts.json); each comparison uses one declared attempt without an outer retry.
The [historical default-start results](results/default-seed-comparisons.csv) and [bounded probes](results/diagnosis.csv) preserve earlier default-start failures and current alternative-start evidence. The historical full default-start campaign was not rerun; the production comparison owns the current cold propylene default-start result.

The historical [independent-code comparison](results/independent-code-comparison.csv) checks all returned pure-liquid densities, binary bubble pressures, and hexane–chlorobutane isobar temperatures against teqp with a fixed relative tolerance declared before calculation.
It uses the earlier `a03132aa` states preserved in `results/historical/` and published parameters; [its run record](results/independent-code-run.json) records the comparator environment, hashes, methods, maxima, and wall-clock time.
That earlier comparison includes the recovered propylene density under its unchanged rule selecting every returned density, even when Engine pressure closure fails. All 200 comparisons pass the fixed 1e-10 relative tolerance: 106 densities, 85 pressures, and 9 isobar temperatures. The maxima remain 1.9761969838327786e-13, 3.530509218307998e-13, and 4.440892098500626e-16, respectively.
This is historical independent-code numerical verification. It was not rerun for the current population and does not qualify its three excluded rows.

The retained tables are [pure predictions and published AAD values](results/pure-saturation.csv), [row-level figure comparisons](results/figure-comparisons.csv), [computed comparison statistics](results/summary.csv), [ternary compositions and equilibrium ratios](results/ternary-comparisons.csv), and [ternary numerical diagnostics](results/ternary-diagnostics.csv).
[Coverage](results/coverage.csv) identifies every original table and figure and its comparison or exact missing source input. Table 3 lacks the cited ideal-gas heat-capacity coefficients and original measured rows; Figure 4 lacks individual isotherm temperatures; Figure 6 lacks experimental critical-temperature values; Figure 16 lacks pressure.
[The run record](results/run.json), [ternary timing](results/ternary-run.json), and [elapsed preparation time](results/elapsed.json) retain calculation and task timing.

Run from the checkout root after its editable Engine environment is prepared:

```bash
bash .codex/environments/setup.sh
.venv/bin/python analyses/reproduction/2001-gross/scripts/digitize.py
.venv/bin/python analyses/reproduction/2001-gross/scripts/run.py
.venv/bin/python analyses/reproduction/2001-gross/scripts/ternary.py
.venv/bin/python analyses/reproduction/2001-gross/scripts/diagnose.py
.venv/bin/python analyses/reproduction/2001-gross/scripts/source_diagnosis.py
.venv/bin/python analyses/reproduction/2001-gross/scripts/present.py
uv pip install --python ~/.cache/epcsaft-comparators/222/feos/bin/python teqp==0.23.2
~/.cache/epcsaft-comparators/222/feos/bin/python analyses/reproduction/2001-gross/scripts/independent_code.py
python3 tools/check-notebook-numbers.py analyses/reproduction/2001-gross/notebook.qmd
bash tools/check-affected.sh origin/main
bash tools/check-all.sh
```

The render requires the local notebook packages `ipykernel`, `nbclient`, and `pyyaml` in addition to the calculation environment:

```bash
QUARTO_PYTHON="$PWD/.venv/bin/python" quarto render analyses/reproduction/2001-gross/notebook.qmd --to html --output-dir ../../../build/gross-2001/rendered
```

The rendered notebook is a presentation of retained values and does not call the Engine.
Exact aggregate-AAD reproduction, unresolved figure differences, and numerical failures remain open for the lead's review.
