# Exact-state reuse in bounded parameter fitting

The [historical notebook](notebook.qmd) retains the bounded state-reuse,
failed-promotion and timing evidence. It retains its exact base content:
the historical native-check duration has no retained CSV/JSON source in this
leaf. The current producer repair and its results are documented here.

## Current producer rerun (#720)

The probe uses the affine target list introduced by `1a8fdec88`. It uses `ObservationEvaluator::Impl` and `ObservationEvaluation`, introduced by `398765244`, after `333595888` removed the internal `Evaluator` and `Evaluation` aliases. The executable links only its own object and loads the prepared shared native library. The included fitting owner performs every calculation.

The [current numerical rows](results/issue-720/residual-promotion-values.csv)
and [diagnostic](results/issue-720/residual-promotion-diagnostic.json) retain
identical residual values and cost after deliberate optimizer-Jacobian
overflow. Independent pressure comparisons pass at the original tolerance.
This is numerical verification of the manufactured representability failure;
it supplies no physical parameter-fit validation or reference replacement.
The uncomputed optimizer-Jacobian entry remains unavailable, not zero.

Run with this checkout's environment and matching Meson tree. The producer
writes under `build/fit-audit-20260929/`:

```bash
env -u EPCSAFT_HSL_LIBRARY EPCSAFT_BUILD_JOBS=6 nice -n 10 \
  .venv/bin/python analyses/2026-fitting-state-reuse/scripts/run-residual-promotion-probe.py
```
