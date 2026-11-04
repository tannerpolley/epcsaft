# Explicit phase discovery cost

This analysis compares the same
finite phase-discovery searches before and after lazy exact Hessian requests,
identical first-point seed-root reuse and removal of the repeated initial Liquid
compilation. The [notebook](notebook.qmd) reports numerical equivalence, actual
work, diagnostic stage times and the separately admitted cost measurements.
The complete 36-block quiet comparison meets the 20% saving target for both
ordinary and prepared calls on its frozen pre-integration images. The matched integrated
A/B images also meet the target in a complete quiet campaign at niceness ten.
Their numerical refresh and cost evidence remain separate from the earlier
images. Integrated operands, raw records, image identities and reduced tables
are retained in `results/integrated-cost/`; their logs and task commands were
under `build/issue-855/build-r5/`, which no longer exists. The earlier paired ranges and control
comparisons remain in `results/ratios.csv`; all ranges are observed, not
confidence intervals. The integrated result awaits independent delivered review.

`inputs/provenance.json` identifies the frozen declarations, source owners and
main-checkout evidence. `results/` owns the displayed values. The task-local
inventory records, original qualification outputs and complete logs were under
`build/issue-855/build-r1/`, and the reconstructed baseline, fresh campaign logs,
image identities and reduced tables under `build/issue-855/build-r2/`; neither
exists now. The seven other measured cases have separate wall/CPU
operands, paired ratios, same-image controls and preserved full states in
`results/other-*`. Their task adapters, logs and reduction command were under
`build/issue-855/build-r3/`; they reused the release interleaving, worker and
`build/issue-855/probe.py:cases` constructors. Those adapters, `probe.py` and
the frozen images `scripts/run.py` loads from `build/issue-855/` are lost; of
the recorded campaign inputs only the #647 ionic falsifier was recovered with
its SHA-256, as `scripts/historical/f5_ionic.py`.
**Frozen record — not replayable from this repository.** `scripts/compare.py`, `reduce.py` and `present.py` still read
the retained records. These are reported costs without
a speed target. This leaf neither overwrites
historical analyses nor replaces their numerical constructors or criteria.

`scripts/compare.py TASK OUTPUT` compares those frozen A/B operation records,
native numerical rows and the original solver parent/trial pairs. `scripts/run.py`
requires both selected prepared environments and their task-local frozen image
guards. It reuses the release workers, declares fresh ordinary/prepared calls,
and records deferred admissions. It selects the next unused measurement file
and refuses to overwrite retained cost blocks. A resumed campaign starts a new
eight-hour admission window with the same load and compiler limits.
`scripts/reduce.py` checks complete, quiet, physically matched blocks before
calculating paired ratios. Image records must match both current frozen images;
all earlier `measurement-0` blocks are retained and excluded after reconstruction
of A. Disturbed operations are repeated in full. A host-deferred campaign supplies
no cost pass.

Presentation performs no thermodynamic calculation:

```bash
nice -n 10 .venv/bin/python analyses/2026-explicit-discovery-cost/scripts/present.py
quarto render analyses/2026-explicit-discovery-cost/notebook.qmd --no-execute
```

The explicit fast-mode section of the notebook uses the separately retained
`results/fast-discovery/` operations, comparisons, work and budget-exhaustion
operands, with frozen inputs in `inputs/fast-discovery.json.gz`. It distinguishes
local two-phase admission from complete search. Full/default preservation and
fast/full state agreement are numerical verification, not physical validation.
The scoped paired-NaN diagnostic amendment and its affected rows are retained
there. Fast-mode delivered review and manuscript promotion remain open.
