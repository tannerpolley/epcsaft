# 2026-pure-saturation-regression

This directory preserves historical methane, ethane, and propane saturation
fit results and the current propane replay through the installed wheel.

## Retained history (July 2026)

The compact JSON summary and the two CSV owners in `results/` retain the unique
propane parameters, comparison rows, source identities, diagnostics, and
bounded failure meaning of the historical campaign. They were produced by an
earlier optimizer build and make no held-out predictive-accuracy claim.

## Current replay (R4, #49)

[`scripts/replay_propane_fit.py`](scripts/replay_propane_fit.py) refits the
three Gross--Sadowski propane parameters through `epcsaft.regression` from
[`inputs/propane-baseline-v3.fit.json`](inputs/propane-baseline-v3.fit.json) (the
versioned request: identities, origins, scales, bounds and the four
multistart points) against the three
[`inputs/source/propane-baseline.csv`](inputs/source/propane-baseline.csv)
coexistence states (Glos, Kleinrahm and Wagner 2004, Table 2): saturation
pressure and saturated-liquid mass density at 210, 270 and 330 K, log-ratio
residuals, unit weights. The packet parameter document
(`data/packets/gross-2001-propane/3`, SHA-256 `6dc4d60f…`; #341 removed the
byte-identical leaf copies of versions 2 and 3, record commit `21d7a112`) supplies the fixed
molar mass; its digest is verified before every run.

The current replay input records explicit neutral charge. Previously retained
results remain bound to the version-2 materialization and wheel written in
those result records; changing the runnable input does not revise them.
The original fit request and source manifest remain at their historical paths
with the recorded SHA-256; the current runner reads their version-3 counterparts.

Retained outputs, dated by run:

- `results/propane-replay-<date>.json` — wheel identity, source digests, every
  start's status, fitted values, the multistart spread, the conditional
  covariance (status, rank, singular values, variance factor, physical
  standard deviations scaled by the variance factor), training residuals,
  typed failures and claim limits.
- `results/propane-replay-<date>.csv` — every Glos 2004 row (roles held_out
  and stress from `analyses/reproduction/2004-glos`) with the source value,
  its expanded uncertainty, the fitted and packet predictions and their
  relative errors; failed rows keep their typed status.

Reproduction class `MODERN_REFIT`: three training states, a different optimizer
build and unit weights, so agreement with the packet or the July tuple is
descriptive. Held-out rows are comparisons against Glos 2004 with their
expanded uncertainties, not an accuracy bar; the 110 K stress row lies outside
the fitted domain and is retained whatever its error. The historical `confirmation-1`
start failed its initial evaluation on the 330 K pressure row with
`CentralSolveUnavailable: boundary seed has no coexisting density roots`.
The [native comparison](results/batch3-comparison.json) retains that failure
and the changed availability: all four starts now converge, including in the
unmodified loop on the same current wheel. The first declared start remains
the reference. Attribution to an intervening Engine change is unverified.

## Notebook data roles

The [retained-result notebook](notebook.qmd#corrected-training-and-held-out-roles)
selects the October replay. It assigns training membership from
`inputs/source/propane-baseline.csv`, keeps the inherited stress declarations,
and treats the other states as held out. In particular, the older packet's
150 K training label becomes held out for this fit. The notebook exposes both
labels and uses the corrected cohort in every metric and plot; the dated
retained CSV is unchanged.
