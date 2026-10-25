# Bülow 2019 permittivity law

Current scope: four native Eq. (7) checks on the explicit three-species basis
at 298.15 K and 1000 mol/m3, each within absolute relative-permittivity error
1e-12, including the zero-water endpoint through
`State.bulk_relative_permittivity` after Engine #792. All four checks pass.
This is numerical verification of a source equation, without experimental
dielectric validation or liquid equilibrium evidence.

From the repository root, using this checkout's editable Engine:

```bash
env -u EPCSAFT_HSL_LIBRARY nice -n 10 .venv/bin/python analyses/reproduction/2019-bulow/scripts/reproduce.py
cd analyses/reproduction/2019-bulow
QUARTO_PYTHON="$OLDPWD/.venv/bin/python" nice -n 10 quarto render notebook.qmd --to html
mkdir -p ../../../build/issue-789/2019-bulow/render
mv notebook.html ../../../build/issue-789/2019-bulow/render/
```

Rendering runs inside the leaf, writes inside `build/issue-789/`, and reads
retained files only. Its HTML is a derived local view. `results/comparison.csv` retains all four declarations,
independent references, errors, unchanged tolerances and evaluator conditions;
`results/run-identity.json` identifies inputs and the loaded native library.
The complete notebook computes its displays from those files.

Source: Bülow, Ji and Held, DOI 10.1016/j.fluid.2019.03.010, Section 2.3 Eq. (7),
Tables 2–3; exact Zotero identities and hashes are in
`inputs/source/reproduction-cases.json`. Water's original table value 78.09 is
explicitly replaced by the paper's 80. Ion parameters remain as printed; zero
binary interactions select the unfitted declaration. The current claim reads
only the native permittivity value. The retained older tables and figure
images remain source/historical records. LLE Figures 3/5 need calibration and
axis mapping, Figure 2 needs its state/derivative basis, and Table 5 needs
measurements and weights. None is declared reproduced. Independent delivered
review remains open.
