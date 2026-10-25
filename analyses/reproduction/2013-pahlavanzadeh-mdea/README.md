# Pahlavanzadeh and Fakouri Baygi 2013 aqueous MDEA

The [research notebook](notebook.qmd) records the live raw-parameter comparison
against Figure 1 pure-MDEA density and Figure 4 bubble/dew temperatures. Some
local intervals intersect; other admitted states disagree, and Figure 4d's
bubble states are refused by native local liquid-stability admission. The
Engine uses sigma-cubed association with the printed tuples, while the paper
prints d-cubed. A supplementary fixed-temperature association probe isolates
that convention difference without changing the production packet or adopting
a conversion or parameter fit.

Source: *Journal of Chemical Thermodynamics* 59 (2013), 214–221,
[doi:10.1016/j.jct.2012.12.021](https://doi.org/10.1016/j.jct.2012.12.021).
The [source record](inputs/source/source-locator.md) preserves all tables,
owned model coordinates and primary Kuranov/Kamps pressure observations.
[Comparisons](results/comparisons.csv) distinguish numerical verification,
physical targets with unavailable predictions, disagreements and refusals.
[Parameter provenance](inputs/parameters/provenance.json) binds the byte-for-byte
packet materialization. Historical source-only records remain in Git.

From this checkout root, after its own `.venv` setup:

```bash
nice -n 10 env -u EPCSAFT_HSL_LIBRARY EPCSAFT_BUILD_JOBS=6 .venv/bin/python analyses/reproduction/2013-pahlavanzadeh-mdea/scripts/reproduce.py
nice -n 10 env -u EPCSAFT_HSL_LIBRARY EPCSAFT_BUILD_JOBS=6 .venv/bin/python analyses/reproduction/2013-pahlavanzadeh-mdea/scripts/diagnose.py
nice -n 10 env -u EPCSAFT_HSL_LIBRARY EPCSAFT_BUILD_JOBS=6 .venv/bin/python analyses/reproduction/2013-pahlavanzadeh-mdea/scripts/diagnose.py --association-only
```

Use `reproduce.py --compare-only` to refresh comparison rows without calculating
or relabeling prior results. Native traces go to `build/issue-736/`; tracked
results are compact and preserve execution identities. Rendering reads retained
values only; the reviewed build stages identical notebook/input/result copies
under `build/issue-736/render/` and keeps HTML and inspection evidence there.

Ideal chemistry awaits [Engine #776](https://github.com/tannerpolley/ePC-SAFT/issues/776)
and separate decisions on the printed carbon/sign conflicts and molecular-feed
mapping. Excess enthalpy, missing primary cohorts and ambiguously owned curves
remain unavailable with their precise causes in the notebook. This is a bounded
reproduction attempt, not complete-paper reproduction or physical validation.
Independent delivered review and manuscript promotion remain open.
