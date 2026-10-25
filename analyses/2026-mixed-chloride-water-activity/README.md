# Mixed chloride water activity

Neither complete declaration qualifies against the measured mixture water
properties: every source and ionic-strength bucket fails its frozen limits,
throughout the finite sensitivity cases. The Figiel binary-fit improvement
fails the required transfer test. The [notebook](notebook.qmd),
[verdicts](results/verdicts.csv) and [all residuals](results/residuals.csv.gz)
retain the negative result.

Stage A is Supported at commit `2aca2f65` by delivered review round 2
([#558](https://github.com/tannerpolley/ePC-SAFT/issues/558)). The accepted scope and its binding corrections are in
[Engine issue #556](https://github.com/tannerpolley/ePC-SAFT/issues/556) and the
[lead resolution](https://github.com/tannerpolley/ePC-SAFT/issues/556#issuecomment-6204697972).

The [frozen comparison specification](inputs/comparison-specification.json)
owns the conditions, measured quantities, acceptance rules, row grouping,
failure treatment and limits. The
[explicit sensitivity realizations](inputs/frozen-realizations.csv) were
generated from printed inputs and committed, with hash-bound input bytes, at
`f594304279ed61e5d8c0e13f094815a9dfc48b1e` before any comparison residual. Each realization
is applied identically to both unchanged parameter declarations. These are
finite conditional scenarios; they are neither exhaustive interval corners nor
independent measurement standard deviations.

The isopiestic sources are judged in their source-referenced practical osmotic
coefficient; the hygrometric source is judged in measured water activity. A
water-activity pass implies no osmotic-coefficient pass. Pure-salt endpoints and
reference anchors remain separate from mixture scores. Absolute qualification
and improvement over the Held declaration are separate questions, both assessed
for each source, system and ionic-strength bucket. A native failure remains in
the record and prevents a passing verdict.

The [input locks](inputs/packet-locks.json) identify complete, unchanged Held
and provisional Figiel declarations and the reusable measured-property packet.
The [file inventory](inputs/input-files.json) hashes every frozen input.
Source-only arithmetic checks and the activity-limit calibration are retained
in `inputs/`; they contain no Engine prediction or parameter-data residual.

The factual tables are owned by
[mixed-chloride-water-observations/1](../../data/packets/mixed-chloride-water-observations/1/packet.toml).
Its source metadata preserves source locators, reference conventions, printed
inconsistencies, missing uncertainties and the ambiguous dispersion statements.
No pointwise weights or confidence claims are inferred from them.

The retained MS-2 inventory and
[application basis](inputs/application-basis.json) describe a raw feed on the
specific-gravity/TDS water-mass proxy. They establish neither a conditioned
extraction inlet nor qualification of the simultaneous Na/Ca/Mg feed. The
comparison is restricted to the specified water quantities, temperature and
measured ternary compositions; individual-ion activities and extraction
chemistry remain outside its claim.

MS-2's approximate ionic strength lies inside the measured ranges. The
extrapolation is in composition: two ternary systems cannot qualify a
simultaneous Na/Ca/Mg feed.

From the repository root, reproduce the study and presentation separately:

```bash
EPCSAFT_BUILD_JOBS=6 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  nice -n 10 .venv/bin/python analyses/2026-mixed-chloride-water-activity/scripts/compare.py
.venv/bin/python analyses/2026-mixed-chloride-water-activity/scripts/check.py
.venv/bin/python analyses/2026-mixed-chloride-water-activity/scripts/present.py
```

The calculation consumes the frozen realization list without tuning or
refitting. It retains every source-specific realized row, and
`results/realization-coverage.csv` maps each unchanged source/bucket explicitly
to its nominal rows. A source-flag omission cannot remove a native failure.
The presentation script reads retained numbers without importing the Engine.
The Quarto project registers this notebook once in its render and sidebar
lists, with execution disabled. Stage B awaits independent delivered review;
it has no investigator promotion for manuscript use.
