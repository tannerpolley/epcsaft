# May methane/ethane coexistence

The current calculation compares all seventeen Table 5 states with native
Gross–Sadowski PC-SAFT using the unchanged packet /2 and an unfitted zero
methane/ethane interaction. [comparison.csv](results/comparison.csv) retains
measured/native compositions, fixed allowances and every unavailable row;
[diagnostics.json](results/diagnostics.json) retains actual native residuals
and causes. This is a physical comparison with measured compositions.
[notebook.qmd](notebook.qmd) computes its tables and plot from retained files.
Independent delivered review is pending.

Run from the repository root after preparing this checkout's editable Engine:

```bash
env -u EPCSAFT_HSL_LIBRARY nice -n 10 .venv/bin/python analyses/reproduction/2015-may/scripts/reproduce.py
cd analyses/reproduction/2015-may
QUARTO_PYTHON=../../../.venv/bin/python nice -n 10 quarto render notebook.qmd
mkdir -p ../../../build/issue-789/2015-may/render
mv notebook.html ../../../build/issue-789/2015-may/render/
```

The retained midpoint feeds are manufactured one-mole declarations, not
experimental feeds. The coalesced, unreturned state remains undetermined;
declaring two phases does not establish global stability. Source densities
are absent, so no density-validation or fixed-density pressure fitting claim
is made. No parameters, uncertainties or comparison limits are fitted.
Historical evidence retains its recorded identity; older template copies
remain at record commit `21d7a112`.
