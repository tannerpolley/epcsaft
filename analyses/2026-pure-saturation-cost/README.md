# Pure-fluid saturation cost

This analysis checks the two-density reduction selected in
[issue #854](https://github.com/tannerpolley/ePC-SAFT/issues/854) against the
unchanged Engine and the release study's FeOs and Clapeyron workers.
The [notebook](notebook.qmd) reads retained state, action, failure and cost
comparisons. It awaits independent delivered review and investigator promotion.

The #854 prepared-call target was missed against both FeOs and Clapeyron.
The interim [issue #865](https://github.com/tannerpolley/ePC-SAFT/issues/865)
package preserves the 251-state regular set bitwise, but its apparent saving
remains unresolved: $-\ln R=0.202224$ is below the measured control envelope
$C=0.210807$. The notebook records this separate comparison against pinned
`0480e79df`; retained operands are in `results/issue-865/`. Beating FeOs remains
with [issue #884](https://github.com/tannerpolley/ePC-SAFT/issues/884).

The calculation reuses release declarations, guarded loading, whole-call timers,
interleaving and quiet-host admission. Its workers use newly frozen images,
recorded in `results/image-A.json` and `results/image-B.json`; the historical
release records remain unchanged. Ordinary calls include compilation. Prepared
candidate calls compile the declaration once outside the timer and start each
solve afresh inside it. Neither endpoint values nor density roots carry between
samples. FeOs and Clapeyron keep their ordinary persistent-worker calls.

Source probes and engineering-gate logs are retained in the main checkout's
`build/issue-854/build-r4/`. Their numerical tables are read by the notebook from
this leaf's `results/`. The image paths identify the prepared checkouts required
by the calculation; a repeated campaign needs new image and raw-output paths.
Rendering uses only retained tables:

```bash
.venv/bin/python analyses/2026-pure-saturation-cost/scripts/reduce.py
.venv/bin/python analyses/2026-pure-saturation-cost/scripts/present.py
quarto render analyses/2026-pure-saturation-cost/notebook.qmd --no-execute
```


The dedicated prepared entry from [issue #886](https://github.com/tannerpolley/ePC-SAFT/issues/886)
is recorded separately in `results/issue-886/` and the notebook's final section;
the earlier tables retain their original boundaries. Its owner-selected
criterion is the paired saving against main, using the same-image envelope.
Fresh comparator costs remain observations rather than acceptance thresholds.
The guarded controller accepts a new output directory through
`EPCSAFT_COST_OUTPUT`; its detached-A phase-pressure bridge and compilation
source are retained with the complete evidence under the main checkout's
`build/issue-886/build-r2/`. `EPCSAFT_COST_PRESSURES` selects that compiled bridge.
Reproduction requires rebuilding the pinned A after its disposable worktree is
removed. Presentation remains a read of retained results.

The longer 30-pair prepared campaign is retained separately in
`results/dedicated-integrated-r4/`, including its rejected first attempt and
accepted complete replacement. Its notebook table and figure are regenerated
from retained values with `.venv/bin/python analyses/2026-pure-saturation-cost/scripts/present_r4.py`.
The same-image saving criterion remains unmet.
