# Belov and Aristova eight-gas restriction

Current scope is the explicitly restricted eight-gas ideal problem at 2000 K,
100000 Pa from one mole C2O, with source standard pressure 101325 Pa. All eight
positive amounts are compared with the independent 80-decimal element-potential
reference using absolute log-amount error 2e-6. This is numerical verification,
not physical validation or the published ten-species heterogeneous answer.

Run from the repository root with this checkout's own editable Engine:

```bash
env -u EPCSAFT_HSL_LIBRARY nice -n 10 .venv/bin/python analyses/reproduction/2023-belov/scripts/reproduce.py
cd analyses/reproduction/2023-belov
QUARTO_PYTHON=../../../.venv/bin/python nice -n 10 quarto render notebook.qmd
mkdir -p ../../../build/issue-789/2023-belov/render
mv notebook.html ../../../build/issue-789/2023-belov/render/
```

The producer compiles a native declaration/reader under `build/issue-789/`,
linked to this checkout's native library. Native code derives source masses,
encodes one-state Gibbs anchors, checks all eight pressure-standard readouts,
and invokes the existing reacting equilibrium driver once. Python transports
source inputs and records identities only. The zero-Cp polynomial with a
[2000,2001] K bookkeeping interval establishes only the supplied Gibbs values
at its 2000 K anchor; no thermal derivative or caloric claim follows.

`results/comparison.csv` retains all amounts, independent references, errors,
limits and evidence labels. Amount errors are absolute logarithmic errors;
Gibbs/RT and closure rows use absolute errors, with limits 1e-10 for Gibbs,
atom inventories in mol, raw/classified rows and relative pressure. Every
residual chemical potential is exactly zero and mechanical slope must be
positive. `results/diagnostics.txt` retains native status/cause, the actual
phase state, all compositions and raw/classified row identities.
`results/run-identity.json` binds source, mass rows, driver and native library.
The notebook computes its displays from these files without rerunning the
Engine. Independent delivered review remains open.

Source: Belov and Aristova (2023), DOI 10.5541/ijot.1317496, test 1 and article
Section 6. Exact source file/commit/blob, fetched CRLF and Git-normalized LF
hashes, source Gibbs/reaction/atom declarations and restricted reference remain
under `inputs/cases/test-01/source/`. Carbon/oxygen mass declarations use the
retained Chase CO/CO2 rows, whose file hash and row names are recorded.
Graphite/diamond, the authors' heterogeneous reference and solid pressure/volume
conventions remain unavailable. No solid volumes, trace clipping, private
solver or caller retry is introduced. The older `test-01-equilibrium.json`
remains historical under its own fingerprint.
