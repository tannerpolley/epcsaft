# Algorithm illustrations

This leaf supplies traced curves for paper figures F4 (methane/ethane phase
envelope and critical candidate) and F3 (water/hexane tangent-plane distance).
Its scope is numerical illustration of the native algorithms, with retained
point comparisons. It makes no physical-validation or global-stability claim.
The [notebook](notebook.qmd) records the quantities, conditions and limitations.

Inputs are declared in [study.json](inputs/study.json). Parameter and reference
files remain with their existing owners; their hashes and the loaded native
extension identity are in [provenance.json](results/provenance.json).
The envelope uses the Gross methane/ethane packet and the equal-composition
C3/E2 case from the open-source equilibrium comparison. C1's pure-propane
critical point is outside this binary envelope's scope.
Water/hexane reuses the parameter conversion in the TPD discovery study's
`feos_compare.py`, including its two-site water association and zero binary
interaction coefficient.

From the repository root, using this checkout's prepared environment:

```bash
nice -n 10 env -u EPCSAFT_HSL_LIBRARY .venv/bin/python analyses/2026-algorithm-illustrations/scripts/calculate.py
nice -n 10 env -u EPCSAFT_HSL_LIBRARY .venv/bin/python analyses/2026-algorithm-illustrations/scripts/figures.py
```

The calculation calls native continuation, equilibrium and stability evaluation.
Python declares a composition grid, collects results and compares retained
values; it performs no thermodynamic search or numerical differentiation.
The plotting script reads retained files only and writes PDF and PNG:
[F4 envelope](figures/F4-phase-envelope.pdf) and [F3 TPD path](figures/F3-tpd-path.pdf).

`envelope.csv` keeps accepted states in native trace order, with convergence
and achieved residuals. `trace.json` keeps the terminal condition, attempt and
retry counts, and the critical candidate and spectrum. There were no rejected
steps in this retained run. `crosschecks.csv` gives relative interpolation
errors, with the corresponding starting segment index, against every retained
FeOs envelope point. Its quantity columns are errors, not state values.
`search.csv` retains the native liquid search points, descent witnesses and
numerical error estimates. `tpd.csv` retains direct liquid trial evaluations,
pressure closure, association diagnostics and tangent-plane distances.

The notebook and figures await independent scientific review and investigator
promotion for manuscript claims. This leaf does not edit the paper.
