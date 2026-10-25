# Why the two engines disagree on Figiel activity coefficients

Current native ePC-SAFT recovers both retained historical alkali-halide activity
series when the dispersion convention is matched. At identical temperature,
density and composition, changing only the same-sign ion dispersion selector
changes dispersion alone. Exclusion also reproduces the digitized Figiel
Figure 5 LiBr model line within the fixed aggregate digitization envelope;
the all-pairs diagnostic control misses that envelope.

The [research notebook](notebook.qmd) computes all result tables, deviations,
criteria and timing from retained output. This is numerical verification of
historical implementation values and a published model line, not physical
validation. Independent delivered review and investigator promotion for
manuscript use remain outstanding.

From the checkout root, with its own editable Engine environment prepared by
[the root Checks instructions](../../../README.md#checks):

```bash
env -u EPCSAFT_HSL_LIBRARY nice -n 10 env EPCSAFT_BUILD_JOBS=6 .venv/bin/python analyses/reproduction/2025-figiel-engine-crosswalk/scripts/crosswalk.py
```

The producer reads [comparison.json](inputs/comparison.json), admits both
native declarations, evaluates the historical compositions and densities,
solves the liquid states at 298.15 K and 100000 Pa, and evaluates each
convention's own native infinite-dilution activity reference. It writes the
current [comparison table](results/data/crosswalk.csv), explicit
[state operands and convergence records](results/data/states.json),
[reference diagnostics](results/data/references.json),
[activity and model-line comparisons](results/data/activities.json),
[native failures](results/data/failures.json), and
[Engine identity, hashes, settings and timing](results/data/run.json).
Empty native failures mean no native call failed; the notebook separately
reports numerical comparison outcomes.

Both declarations consume packet
`data/packets/figiel-2025-reference-electrolytes/4/parameters/parameters.json`.
The retained [all-pairs declaration](inputs/all-pairs.json) differs only in
`ionic_dispersion`; it is a diagnostic control contrary to Figiel's convention,
not an adopted parameterization or a second independent engine.
Packet versions /3 and /4 have identical physics for these aqueous states:
the only changed rules are two water–ethanol site edges removed by component
selection. This version change cannot explain an aqueous discrepancy.
The sibling [Figure 5 digitization](../2025-figiel/figures/figure-05/source/libr-model-digitized.md)
owns the source coordinates and uncertainties and is consumed without edits.
The broader Figiel reproduction in #748 remains separate.

The raw [superseded](results/data/superseded.json) and
[historical greenfield](results/data/greenfield.json) JSON files remain
byte-for-byte from comparison commit
`eadac899f8e4ebbbe2408c252fadacbfabd30e2b`.
Neither identifies its installed binary or exact engine commit.
Their `ares` fields have different bases: J/mol for superseded and reduced
Helmholtz energy for greenfield. The old mixed-basis crosswalk is replaced by
comparison against the five-term reduced sums; the dimensional original is
preserved and normalized only for the notebook's basis check.
The earlier rounded README same-density table came from a separate probe with
a different water molar mass; the current producer uses exact archived
compositions and densities.

The [historical literature diagnosis](results/like-ion-dispersion-convention.md)
is retained as historical evidence. Retired probes are available in Git at
`21d7a112`; they are not current rerun commands. The present paired native
experiment supplies the current convention diagnosis. Positive-pressure states
satisfy the Euler identity at the accepted tolerance; off-root states with
nonpositive compressibility have unavailable log fugacity because ln(Z) is
undefined. They retain finite reduced terms and are not failed density roots.
