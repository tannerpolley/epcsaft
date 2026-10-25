Source: Chase (1998), NIST-JANAF Thermochemical Tables, 4th ed. (`inputs/cases/`).
`results/` keeps the earlier chemical-equilibrium validation ([summary](results/chemical-equilibrium-validation.md)); its producer is not retained. Earlier README text and template or copied files removed by #341 remain at record commit `21d7a112`.

Current numerical verification covers the fourteen printed Shomate property
rows, four historical source-derived reaction references and six ideal-gas
reacting states (eighteen species amounts). All values, errors and unchanged
limits are in [comparison.csv](results/comparison.csv), with native diagnostics
in [solutions.json](results/solutions.json). The [notebook](notebook.qmd) reads
these files without running the model. Historical amount references are not
published measurements. Independent delivered review is pending.

From the repository root, using this checkout's prepared editable Engine:

```bash
env -u EPCSAFT_HSL_LIBRARY nice -n 10 .venv/bin/python analyses/reproduction/1998-chase/scripts/reproduce.py
cd analyses/reproduction/1998-chase
QUARTO_PYTHON=../../../.venv/bin/python nice -n 10 quarto render notebook.qmd
mkdir -p ../../../build/issue-789/1998-chase/render
mv notebook.html ../../../build/issue-789/1998-chase/render/
```

The producer preserves source intervals and the 100000 Pa pressure standard;
no residual EOS model, fitted reaction constant or caller equilibrium solver
is substituted. Rerun gates are `tools/check-affected.sh origin/main`,
`python3 tools/check-analysis-layout.py` and
`python3 tools/check-notebook-numbers.py analyses/reproduction/1998-chase/notebook.qmd`.
