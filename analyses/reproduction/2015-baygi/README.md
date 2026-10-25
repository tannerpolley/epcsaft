# Baygi 2015 pure MEA saturation

The current native Engine evaluates the unmodified Table 2 MEA 2B, 3B and 4C
parameter tuples with sigma-cubed association. Figures 2–3 are compared at the
fixed drawing-resolution tolerances. The warm 3B pressure–density mismatch
remains undetermined; exact Table 2 aggregate replay also remains undetermined
because the author grid and density-correlation grouping are unknown.

Run from the repository root with this checkout's editable Engine environment:

```bash
nice -n 10 env -u EPCSAFT_HSL_LIBRARY EPCSAFT_BUILD_JOBS=6 .venv/bin/python analyses/reproduction/2015-baygi/scripts/run_pure_saturation.py
```

The [research notebook](notebook.qmd) reads the retained
[physical values and comparisons](results/current-pure-saturation.csv) and
[native diagnostics, declarations and Engine identity](results/current-pure-saturation.json).
It contains the source/model comparison figure and explains the scientific limits.
The primary source is [Fakouri Baygi and Pahlavanzadeh (2015)](https://doi.org/10.1016/j.cherd.2014.07.017),
Section 3.1, Tables 1–2 and Figures 2–3. Source-linked tables, calibrated vector
coordinates, extraction calibration and raster overlay are in `inputs/source/`.
The original Table 2 association volumes are used without Python rescaling or
substitution of the converted application tuple. There is no physical validation
against raw measurements and no claim of global phase stability.

The accepted design was independently reviewed as Ready against main commit
`69f61380`; that review is `build/issue-725/readiness-r1.md` in the main checkout.
Independent delivered review returned Supported (`build/issue-725/delivered-r1.md` in the main checkout).
Render without solving the model by copying `notebook.qmd` into
`build/issue-725/` and running Quarto there with this leaf as `--execute-dir`
and the worktree's `.venv/bin/python` as `QUARTO_PYTHON`.

The historical [saturation table](results/baygi-pure-saturation.csv),
[summary](results/baygi-pure-saturation-summary.json),
[pre-cutover diagnostics](results/pre-cutover-baseline.json) and
[boundary probes](results/feasibility-boundary-sweep.csv) retain their original meaning.
Their producer at record commit `21d7a112` used wheel SHA-256
`a0ff54417a7e88aab5cf2396804f9f5741296097b3b5328a012123db57168d63` and input SHA-256
`474748e118e421e7a1fb4d7aa63bb300cdba1945bab46fe73e30804263ebf139`.
The retired feasibility bound `2e-8` qualifies none of the current states.
The later [association-diameter study](../../2026-mea-association-diameter/README.md)
provides retained historical paired-convention expectations, not a current rerun command.
The new source-linked Table 1/2 CSVs replace the inherited `table-002` and
`table-003` extraction copies; the other source tables and figures are retained
outside the scope of this pure-MEA comparison.
