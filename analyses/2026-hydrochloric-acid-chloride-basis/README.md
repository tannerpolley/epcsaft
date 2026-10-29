# HCl on the frozen chloride basis

This analysis asks whether effective hydrogen-ion parameters can represent the
evaluated binary acid properties and predict the reserved common-chloride
mixture measurements while preserving the existing chloride parameter basis.
The approved design and its readiness review are identified by hash in
[the retained source identities](inputs/source-identities.json).

[The conditional profile](inputs/fit-declaration-issue-901.json) finds a
binary/water-guard-preserving diagnostic point meeting the Jiang I = 5 target.
The [notebook](notebook.qmd) presents its retained values, transfer limitations
and pending independent delivered review. No packet is published. Earlier
trials below remain historical records.

[The one-coordinate declaration](inputs/fit-declaration-issue-798-r2.json)
fixes the hydrogen sizes and energy and treats zero H–Cl correction as a
modelling assumption. Only the H–water correction is fitted. This trial retains every start under `results/issue-798-r2/`. The
[predeclared retrial](inputs/fit-declaration-issue-798-retrial.json) changes only
the fixed H–Cl correction to one, removing H–Cl dispersion as a second
modelling assumption, and retains all starts under `results/issue-798-retrial/`.
Both trials identify and fail the original binary accuracy limits. With the
owner's osmotic RMS limit relaxed to 0.011, the zero H–Cl endpoint passes
binary checks but fails Roy held-out mixture qualification: mean/maximum
relative errors 19.963996%/50.374497% against 5%/10%. No packet is written.
Earlier declarations and result files are preserved.

The native producer is `scripts/activity.cpp`; `scripts/run.py` declares the
parameter document and transports it to that producer. Native Ceres owns all
optimization and exact derivatives; native fitting diagnostics supply singular
values and unscaled inverse information. The [notebook](notebook.qmd) presents
identification separately from physical validation and later qualification.

Run from the repository root with this checkout's prepared environment:

```bash
env -u EPCSAFT_HSL_LIBRARY nice -n 10 agent-heavy .venv/bin/python analyses/2026-hydrochloric-acid-chloride-basis/scripts/run.py
```

The retrial uses the same producer with explicit declaration and result paths:

```bash
env -u EPCSAFT_HSL_LIBRARY nice -n 10 agent-heavy .venv/bin/python analyses/2026-hydrochloric-acid-chloride-basis/scripts/run.py analyses/2026-hydrochloric-acid-chloride-basis/inputs/fit-declaration-issue-798-retrial.json analyses/2026-hydrochloric-acid-chloride-basis/results/issue-798-retrial
```

The producer refuses to overwrite a retained fit. A changed scientific
question or fitting declaration requires owner disposition before another fit.
Measured mixture coverage does not extend to the concentrated lithium
application request. The evaluated binary training curve and corrected water
comparison are not independent physical validation.

Earlier native normalization evidence remains retained in `results/r3/`.
The fit command does not rerun those comparisons or qualification.

Prepare presentation tables and binary scores for both trials from retained values,
then render the registered notebook:

```bash
.venv/bin/python analyses/2026-hydrochloric-acid-chloride-basis/scripts/present.py
quarto render analyses/2026-hydrochloric-acid-chloride-basis/notebook.qmd
```

The evaluation-only continuation is `scripts/qualify.py`; it refuses an existing
`results/binary-water-qualification` destination. It reuses the selected endpoint
without fitting and stops at the first required failure. Its retained tables
and the notebook distinguish calibration from held-out qualification.
