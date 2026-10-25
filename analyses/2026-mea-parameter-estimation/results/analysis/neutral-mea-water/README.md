# Neutral MEA-water qualification

This analysis independently replays the 25 admitted non-endpoint Cai et al.
(1996) MEA-water VLE rows at 66.66 and 101.33 kPa through the public installed
Engine wheel. It fits only `pair/monoethanolamine/water/k_ij`. Reactive and
holdout observations are excluded.

## Result

The frozen pressure-grouped one-SE rule selects 2B without ambiguity. Its
complete-row fit is
`k_MEA,H2O = -0.07352749585113796`; the conditional local-linearized 95%
interval is recorded with the native fit in
[`results/neutral-mea-water.json`](results/neutral-mea-water.json). The interval is
strictly inside the frozen `[-0.5, 0.5]` bounds, and every 2B training, held,
and complete-refit row reaches a certified phase result with optimizer
convergence.

| Topology | MEA donor/acceptor multiplicity | Grouped mean score | Grouped SE | Phase | Bounds | Uncertainty | Disposition |
|---|---:|---:|---:|---|---|---|---|
| 2B | 1 / 1 | 9.2127559013 | 1.8448718706 | pass | pass | pass | selected |
| 3B | 1 / 2 | 102.7292591340 | 2.4396974659 | pass | pass | pass | outside one-SE threshold |
| 4C | 2 / 2 | 393.0699141866 | 2.9397681134 | pass | pass | pass | outside one-SE threshold |

The one-SE threshold is `11.0576277719`; 2B is its only member. The public
Engine resolves both MEA-water donor/acceptor cross edges from the declared
arithmetic-energy/geometric-volume rule to `2506.0 K` and
`0.0411037990944876`. Packet v4 restores the cited Wolbach–Sandler size factor; a current-API replay
moves the best `k_MEA,H2O` to `-0.0738817`, inside the interval above
([requalification](../../../../2026-wolbach-sandler-requalification/README.md)). The native 4C fit now completes, replacing the historical bracket-limited
failures; its score remains outside the selection threshold. The
[comparison](results/batch1-before-after.csv) retains old values and misses.

The provisional packet ends at 393.15 K, below 12 admitted Cai temperatures.
The replay therefore declares a candidate execution domain through 460 K and
qualifies the coordinate only over the observed 362.81–431.85 K, 66.66–101.33
kPa range. Values above 393.15 K remain `candidate_extrapolation`; this replay
makes no claim outside the admitted Cai observations.

The historical topology figure is retained in `figures/output/`. The native
table above replaces its incomplete-fit label for 4C.

Exact plotted values are in
[`figures/output/topology-selection.csv`](figures/output/topology-selection.csv).
Complete row results, public Engine fingerprints, wheel hash, input hashes,
uncertainty, gate decisions, and failure diagnostics are retained in
[`results/neutral-mea-water.json`](results/neutral-mea-water.json).

## Claim boundary

This is a `MODERN_REFIT` against Cai MEA-water VLE, not an exact author-method
replay. The uncertainty is conditional on the selected topology, frozen
residual scales, and local linearization; it is not source measurement
uncertainty. The fitted binary coordinate is correlation evidence and does not
validate reactive mixtures, Xu pressure observations, or Jakobsen speciation.

Reopen the selection if another topology enters the retained one-SE set while
passing all gates, or if 2B loses phase success, interior bounds, or finite
interior uncertainty under the unchanged Cai rows and scoring.

## Native fitting and historical comparison

`../../../scripts/neutral_mea_water_run.py` declares native fits with the same observations, residual
scales, coordinate and bounds. The earlier wheel is identified in the retained
before/after record; its old selection is retained there rather than promoted
as the current result. Packet requalification remains a separate comparison.
The current topology table is the source for the figure below.

The figure renderer reads only the retained plot CSV and still runs:

```bash
MPLBACKEND=Agg .venv/bin/python \
  analyses/2026-mea-parameter-estimation/scripts/neutral_mea_water_render.py
```
