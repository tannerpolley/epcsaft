# Delivery checks

Checked against `origin/main` at Engine source `c358de722` in this checkout's
editable `.venv`, with `EPCSAFT_HSL_LIBRARY` unset.

- Native calculation completed: both traced branches, critical equations,
  all supplied liquid trials, and retained-point assertions passed.
- `bash tools/check-affected.sh origin/main`: passed, including declaration-only construction.
- `python3 tools/check-analysis-layout.py`: passed.
- `python3 tools/check-notebook-numbers.py analyses/2026-algorithm-illustrations/notebook.qmd`: passed.
- `nice -n 10 bash tools/build-docs.sh`: passed; zero Sphinx warnings, all eleven example notebooks executed.
- `bash tools/check-all.sh`: passed.
- Quarto HTML rendering of `notebook.qmd` with execution disabled: passed.
- PDF and PNG figures: visually inspected for labels, series identity and clipping.
- Repository formatter-based executable count: `calculate.py` 112 lines plus `figures.py` 38 lines, total 150 gross additions; zero Engine lines.

Build logs remain under `build/issue-501/algorithm-*.log` and
`build/docs/`; the HTML notebook view is
`build/issue-501/algorithm-notebook/notebook.html`.
Independent scientific review and investigator promotion remain with the paper owner.
