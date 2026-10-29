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
are retained in `results/integrated-cost/`, with logs and task commands in the
main checkout’s `build/issue-855/build-r5/`. The earlier paired ranges and control
comparisons remain in `results/ratios.csv`; all ranges are observed, not
confidence intervals. The integrated result awaits independent delivered review.

`inputs/provenance.json` identifies the frozen declarations, source owners and
main-checkout evidence. `results/` owns the displayed values. The task-local
inventory records, original qualification outputs and complete logs remain at
`build/issue-855/build-r1/` in the main checkout. The reconstructed baseline,
fresh campaign logs, image identities and reduced tables are retained in
`build/issue-855/build-r2/`. The seven other measured cases have separate wall/CPU
operands, paired ratios, same-image controls and preserved full states in
`results/other-*`. Their task adapters, logs and reduction command are retained
in `build/issue-855/build-r3/`; they reuse the release interleaving, worker and
`build/issue-855/probe.py:cases` constructors. These are reported costs without
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
