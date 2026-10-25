# Nasrifar and Tafazzol 2010 MDEA

Current native attempt for the neutral MDEA–water isobars. Selected dew readings
reproduce within the frozen source rectangles; bubble readings do not, with
cause undetermined after numerical refinement and bounded native checks.
Methane scoring is withheld because the source axis conflicts with its caption.
Reactive cases are unavailable pending [Engine #776](https://github.com/tannerpolley/ePC-SAFT/issues/776).
The paper's water-balance/hydronium inconsistency is recorded separately.

Read the [research notebook](notebook.qmd) for the comparisons, physical meaning,
limits and execution identities. [Source notes](inputs/source/README.md) distinguish
paper model curves from measurements and printed resolution from uncertainty.
All original cases remain in the notebook. This is molecular PC-SAFT with the
paper's MDEA 3B and water 2B parameters, not an adopted electrolyte application model.

From this checkout, after its own environment setup:

```bash
env -u EPCSAFT_HSL_LIBRARY nice -n 10 env EPCSAFT_BUILD_JOBS=6 .venv/bin/python analyses/reproduction/2010-nasrifar-mdea/scripts/reproduce.py
env -u EPCSAFT_HSL_LIBRARY nice -n 10 env EPCSAFT_BUILD_JOBS=6 .venv/bin/python analyses/reproduction/2010-nasrifar-mdea/scripts/diagnose.py
```

Refresh comparisons without solving or changing result identity:

```bash
.venv/bin/python analyses/reproduction/2010-nasrifar-mdea/scripts/reproduce.py --compare-only
```

The source targets were committed before predictions. Manual calibration, fixed
pixel bounds and second readings are retained; crops, inspected QA overlays,
full traces and a standalone notebook render are archived under `build/issue-734/`.
The previous source-only record is preserved in Git at `21d7a112`; its retired
recipes are not current rerun instructions. Delivered scientific review remains
with the lead; the notebook is not promoted for manuscript use.
