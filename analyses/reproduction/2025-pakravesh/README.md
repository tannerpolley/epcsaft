# Pure-MDEA PC-SAFT comparison

The current producer evaluates all 120 S6 states for density and exact response
functions and ten S9 saturation temperatures. It preserves every literal
printed PC-SAFT target, including the duplicated DEA/MDEA S9 columns.
[comparison.csv](results/comparison.csv) retains values, source targets,
errors and unchanged printed-digit limits;
[diagnostics.json](results/diagnostics.json) retains closure and numerical
causes, including adjacent-density pressure probes for the lowest-temperature
states. The [notebook](notebook.qmd) reads these retained files. This is
numerical verification against printed model outputs. All saturation states
pass the independently Ready stiffness-scaled pressure and common-standard
chemical-potential criteria; every printed-model comparison still exceeds its
unchanged source limit. The cold-state pressure readout does not establish an
Engine defect, and the S6 discrepancy cause remains undetermined. The literal
S9 duplication stays visible. Independent delivered review of this rerun is
pending; the prior Not supported verdict is not superseded here.

Run from the repository root with this checkout's editable Engine:

```bash
env -u EPCSAFT_HSL_LIBRARY nice -n 10 .venv/bin/python analyses/reproduction/2025-pakravesh/scripts/reproduce.py
cd analyses/reproduction/2025-pakravesh
QUARTO_PYTHON=../../../.venv/bin/python nice -n 10 quarto render notebook.qmd
mkdir -p ../../../build/issue-789/2025-pakravesh/render
mv notebook.html ../../../build/issue-789/2025-pakravesh/render/
```

The separate PρT-SAFT-HR primary model, caloric quantities and mixtures are
excluded. Article section 3.2 specifies mixture kij=0; remaining mixture work
is verifying association/combining declarations and extracting selected rows,
not inventing missing binary parameters. Caloric comparisons require source
ideal-gas inputs. No source values or tolerances are fitted or corrected.
Earlier removed methods remain at record commit `21d7a112`.
