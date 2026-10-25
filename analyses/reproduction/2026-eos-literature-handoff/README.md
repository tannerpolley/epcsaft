# Selected EOS literature comparisons

The current native Engine reproduces the selected Figiel model solvation
energies at printed precision. The selected activity coefficients, density and
Figiel calibration averages do not reproduce their printed values. The aqueous
NaCl coefficient shows limited agreement with Hamer–Wu under the adopted Figiel
parameters; the mixed-solvent activity and ethanol/NaBr density differences
remain unresolved; the calibration-average misses are the source's own retained
model-to-average residuals. See the [notebook](notebook.qmd).

The scope is the eight distinct literature values in the source-bound tests
removed at raw-parameter adoption, plus the five separately printed Figiel
model energies. This is a small cross-source comparison, not a full-paper
study. All live cases use Figiel reference packet `/4`, with the removed density
test's explicit ionic molar masses restored in a leaf-local declaration because
the packet omits them. The initial native refusal and corrected calculation are
both retained. No packet or Engine source changes are required.

Run from this checkout's root, using its own editable Engine environment:

```bash
env -u EPCSAFT_HSL_LIBRARY nice -n 10 env EPCSAFT_BUILD_JOBS=6 .venv/bin/python analyses/reproduction/2026-eos-literature-handoff/scripts/reproduce.py
```

The script retains native quantities, physical declarations, input hashes,
returned numerical diagnostics and execution identity in
[results/native-results.json](results/native-results.json). It emits all thirteen
[comparison rows](results/comparisons.csv), with signed deviations, fixed
printed-value criteria and distinct physical-comparison/numerical-verification
labels. Historical limits never determine the current outcome.
`--compare-only` refreshes comparisons from retained quantities while preserving
their original execution identity; it performs no new thermodynamic calculation.
The notebook reads retained values only. The selected HTML render and delivery
logs belong under `build/issue-730/`.

The [source selection](inputs/source/literature-targets.csv) gives table/page,
DOI and removed-test locators; [source checks](results/source-checks.json) retain
the matching rows and hashes of the existing observation owners. The original
assertions are readable at `d7738b23d^:engine/tests/`; the full recovered
assertion and raw-parameter inventories remain in `build/issue-730/` evidence.

The former positive-pressure replay did not compare paper values. Its script,
`eos_literature_handoff.py`, was removed by #341. The
[wheel and source lock](inputs/source/eos-artifact-lock.json) and
[packet locks](inputs/source/data-packet-locks.json) preserve historical identity,
including the later stale MEA fingerprint; they never select live inputs.
The former output remains in Git at
`9a7fc51a^:analyses/reproduction/2026-eos-literature-handoff/results/receipts/eos-literature-handoff.json`.
No old wheel is installed or replayed. Full Esteso, Held and Figiel campaigns
remain in their own leaves; no predictive validation or global phase-stability
claim follows from these isolated comparisons.
