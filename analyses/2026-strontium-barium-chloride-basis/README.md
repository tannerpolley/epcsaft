# Sulfate, Sr and Ba on the frozen chloride basis

The Ba binary-plus-Robinson–Bower screen rejects every frozen point because the
two weighted dispersion columns exceed the correlation limit. All observations,
columns and state/reference diagnostics are available. Fitting, identification
and binary accuracy were not reached. The [notebook](notebook.qmd) presents the
approved Hamer–Wu reference replacement, target reproduction, complete screen,
transfer disclosure and limits. El Guendouzi's fixed check rows remain unscored;
every successful fit returns to the owner before their evaluation. No packet
was written, and independent delivered review remains with the owner.

The [Ba mixture declaration](inputs/fit-ba-mixture.json) uses the existing shared
producer, screen, scorer and presentation, preserving all earlier result bytes.
Fresh authorized execution uses this checkout's own environment and refuses to
overwrite retained results:

```bash
env -u EPCSAFT_HSL_LIBRARY nice -n 10 .venv/bin/python analyses/2026-hydrochloric-acid-chloride-basis/scripts/run.py analyses/2026-strontium-barium-chloride-basis/inputs/fit-ba-mixture.json analyses/2026-strontium-barium-chloride-basis/results/ba-mixture
env -u EPCSAFT_HSL_LIBRARY nice -n 10 .venv/bin/python analyses/2026-hydrochloric-acid-chloride-basis/scripts/present.py analyses/2026-strontium-barium-chloride-basis/inputs/fit-ba-mixture.json
```

The following describes the retained sulfate and Sr attempts.

The sulfate candidate is admitted with the unchanged chloride-basis Li–Mg–Cl
mixing term. All four prescribed points fail local resolution, with complete
finite observations, positive exact-column norms and interior coordinates.
Fitting, identification and binary accuracy scoring were not reached.
The [notebook](notebook.qmd) presents the retained screen, derivative checks,
nonempty state/reference diagnostics, transfer assumptions and limitations.

The [frozen declaration](inputs/fit-sulfate.json) uses the existing shared native
producer, transport, scorer and presentation. A fresh authorized execution is:

```bash
env -u EPCSAFT_HSL_LIBRARY nice -n 10 .venv/bin/python analyses/2026-hydrochloric-acid-chloride-basis/scripts/run.py analyses/2026-strontium-barium-chloride-basis/inputs/fit-sulfate.json analyses/2026-strontium-barium-chloride-basis/results/sulfate
```

The existing result directory refuses replacement. Presentation reads retained
values only through the same shared `present.py` and declaration. The
[earlier source/admission record](inputs/sulfate-source-provenance.json) remains
unchanged as evidence of the first attempt's admission blocker. Engine issue
#814 resolved that blocker; the new result is a numerical local-resolution
failure, without a fitted sulfate coefficient. Both reserved sulfate studies
remain unopened and unhashed; no packet is written. Any later successful fit
returns to the owner before reserved checks. Independent delivered review of
this sulfate result remains with the lead.

The following describes the retained Sr and Ba attempts.

The Sr diameter-and-water-dispersion screen fails the all-points correlation
rule. All observations and exact columns are finite, and native physical-partial,
centered-difference and prior-screen comparisons pass. Fitting, identification,
binary accuracy and selected-endpoint Clegg metrics were not reached.
The [notebook](notebook.qmd) gives the screen and nonempty state/reference
diagnostics, with all earlier results retained byte-identically.
This is a lead decision: every successful fit returns to the owner before
qualification or a parameter packet. Guo and Ba mixtures were unopened at that stage.

The [declaration](inputs/fit-sr-sigma.json) fixes the tied length coordinate,
water dispersion coordinate, assumptions, screen, conditional starts, limits,
verification and record hashes. A fresh authorized execution refuses to
overwrite retained output:

```bash
env -u EPCSAFT_HSL_LIBRARY nice -n 10 .venv/bin/python analyses/2026-hydrochloric-acid-chloride-basis/scripts/run.py analyses/2026-strontium-barium-chloride-basis/inputs/fit-sr-sigma.json analyses/2026-strontium-barium-chloride-basis/results/sr-sigma
```

Presentation reads retained values only:

```bash
.venv/bin/python analyses/2026-hydrochloric-acid-chloride-basis/scripts/present.py analyses/2026-strontium-barium-chloride-basis/inputs/fit-sr-sigma.json
quarto render analyses/2026-strontium-barium-chloride-basis/notebook.qmd --no-execute
```

The following describes the historical binary-only attempts.



Both independent two-coefficient screens reject every frozen point because
their finite native columns exceed the declared correlation limit. No
two-coefficient fit or binary scoring is reached. The original single-coordinate
trials remain retained: both identified, then failed evaluated-binary accuracy. The [notebook](notebook.qmd) reads retained sensitivities, fitting
diagnostics, selected endpoints and binary errors. Binary accuracy is calibration
evidence, not physical validation. Independent delivered review remains pending.

The frozen design and readiness records are identified in
[the Sr screen declaration](inputs/fit-sr-two-k.json) and
[the Ba screen declaration](inputs/fit-ba-two-k.json).
[Source provenance](inputs/source-provenance.json) owns table locators, input
basis and exclusions. Each candidate preserves every old component, pair,
correlation, mixing declaration and materialization value of the chloride basis.
Only its own ion is appended. No Engine source changes or reserved-source work
are included; no parameter packet is created.

From the repository root, a fresh authorized execution uses the shared producer
with this checkout's environment, Sr first and then Ba:

```bash
env -u EPCSAFT_HSL_LIBRARY nice -n 10 .venv/bin/python analyses/2026-hydrochloric-acid-chloride-basis/scripts/run.py analyses/2026-strontium-barium-chloride-basis/inputs/fit-sr-two-k.json analyses/2026-strontium-barium-chloride-basis/results/two-k/Sr
env -u EPCSAFT_HSL_LIBRARY nice -n 10 .venv/bin/python analyses/2026-hydrochloric-acid-chloride-basis/scripts/run.py analyses/2026-strontium-barium-chloride-basis/inputs/fit-ba-two-k.json analyses/2026-strontium-barium-chloride-basis/results/two-k/Ba
```

The producer refuses to overwrite either retained attempt. All frozen screen points
and their full-file and per-coordinate binary64 column hashes precede any
conditional fitting. A failed point stops its ion; identities and retained hashes
are verified before the wrapper returns or permits fitting. Native C++ owns
exact derivatives, weighting, fitting and identification. The shared scorer
summarizes retained selected-endpoint predictions without evaluating the model.

Refresh presentation only, then render the notebook:

```bash
.venv/bin/python analyses/2026-hydrochloric-acid-chloride-basis/scripts/present.py analyses/2026-strontium-barium-chloride-basis/inputs/fit-sr-two-k.json
.venv/bin/python analyses/2026-hydrochloric-acid-chloride-basis/scripts/present.py analyses/2026-strontium-barium-chloride-basis/inputs/fit-ba-two-k.json
quarto render analyses/2026-strontium-barium-chloride-basis/notebook.qmd
```
