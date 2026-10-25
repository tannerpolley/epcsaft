# Propane saturation against Glos measurements

Fixed Gross–Sadowski PC-SAFT parameters are compared with every retained Glos
propane saturation-pressure and liquid/vapor-density observation. The
[notebook](notebook.qmd) reads the current [comparison rows](results/propane-saturation.csv),
[supplementary NIST comparison](results/nist-comparison.csv), and
[closure refinements](results/refinement.csv). Source extraction, native closure
and refinement are numerical verification; Glos comparisons are physical
comparisons. Differences outside supplied uncertainty intervals remain visible.

From the repository root, with this checkout's editable Engine installed:

```bash
env -u EPCSAFT_HSL_LIBRARY nice -n 10 .venv/bin/python analyses/reproduction/2004-glos/scripts/reproduce_propane_saturation.py
```

Refresh differences without native calculation or execution-identity changes:

```bash
.venv/bin/python analyses/reproduction/2004-glos/scripts/reproduce_propane_saturation.py --compare-only
```

The runner checks the existing source-value checks before calculation. Each native
call retains its own typed status, closure criterion and termination record;
[execution identity](results/execution.json) links all calls to their loaded
Engine and source inputs. No parameter fitting or Python thermodynamic solver
is used. NIST is a supplementary reference-equation comparison between different
equations of state, not independent experimental validation.

ThermoML datasets 1–3 control the Glos source locators.
Table 2 (locator inherited from retained source audit; PDF not available and table/page not independently verified in this review).
The historical source audit and fitting partitions remain unchanged in
`inputs/source/`; they do not exclude observations from this comparison.
Historical record commits `ca719012c31e8b41edcb98ae82a388117b58d343` and
`21d7a112` preserve earlier provenance and retired recipes.

For standalone rendering, copy this leaf's notebook, inputs and results byte for
byte into `build/issue-732/report/`, convert the notebook with `quarto convert`,
execute only its presentation cells with this checkout's Jupyter kernel, and
render the resulting `.ipynb` to HTML with execution disabled. The retained
render and commands are under `build/issue-732/`; the repository Sphinx build
does not render this notebook. Independent delivered review remains pending.
