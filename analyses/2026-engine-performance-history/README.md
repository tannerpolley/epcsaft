# Engine performance history

**Frozen record — not rerunnable as a current-main calculation.** Record snapshot:
commit `ca719012c31e8b41edcb98ae82a388117b58d343`. Original wheel identities, source
inputs and retained numerical evidence keep their recorded provenance. The multi-wheel
campaign requires the identified historical wheels and external timing inputs under
`build/issue-533/`; running current declarations alone cannot reproduce that comparison.

How did the retained Engine versions change the cost of equilibrium, density and
the frozen MEA F1 residual/Jacobian evaluation on this host? The
[notebook](notebook.qmd) gives the headline table, operation-group charts,
returned-value checks and limits for [issue #549](https://github.com/tannerpolley/ePC-SAFT/issues/549),
extended through October 19 by [issue #644](https://github.com/tannerpolley/ePC-SAFT/issues/644).
The evidence is numerical verification and timing, pending the lead's review.

The study reuses each wheel's commit-matched declarations from `tools/time-key-solves`
and the existing competitor, head-to-head and methane/ethane analyses. The
same-host comparator times and density sequence come from the main checkout's
`build/issue-533/`. Each retained wheel has an isolated environment; detached
worktrees supply declaration scripts only. Preparation never compiles a wheel.

The retained campaign was prepared and measured with the historical environments
named above; its executable rerun recipe is withdrawn.

The campaign refuses to overwrite a prior or interrupted `returns.jsonl.gz`.
Its records contain
three repetitions with independently bracketed current-version A=A controls,
single-thread numerical pools, a fixed CPU and per-block load. It admits quiet
conditions when available, otherwise the authorized relaxed load gate; it waits
for active compilers. Full parallel fits are outside this campaign.

The extension retains the original operations and adds separate rows for the six
methane/ethane T4 Discover points, propane saturation at 150 K, and equimolar
methane/ethane bubble/dew points at 230, 240 and 250 K. Their inputs and sources
are in `inputs/extension.json`. An older wheel that cannot declare or solve a row
records its exact diagnostic and no timing. Missing wheels are recorded without
compilation or substitution.

The current wheel is `1cd3ab15`. All ordinary workers unset
`EPCSAFT_HSL_LIBRARY`; a second current worker and its own independent control set
it to `$HOME/.cache/epcsaft-hsl/coinhsl-lib/libhsl.so`. Each return records the
process solver selection and each equilibrium row the solver used; an empty
equilibrium solver means Newton did not call Ipopt. No HSL file is copied.
MA57 measurements acknowledge: HSL, a collection of Fortran codes for large-scale
scientific computation. https://www.hsl.rl.ac.uk/

Plot and render from retained values in a presentation-only environment with
Matplotlib and Jupyter installed, without the Engine:

```bash
bash analyses/2026-engine-performance-history/scripts/render.sh
```

The render command verifies that the presentation environment contains no
Engine, records hashes of its inputs, renders a copy under `build/`, checks that
rendering leaves the retained inputs unchanged, and copies the self-contained HTML into the main checkout's
`build/issue-549/notebook.html`.
The task checkout's presentation environment is prepared with `uv venv --python
3.13 build/issue-549/presentation` and `uv pip install --python
build/issue-549/presentation/bin/python matplotlib jupyter`.

`results/headline.csv` and `results/statistics.csv` own the displayed values.
`results/returns.jsonl.gz` retains every timed native return.
`results/return-checks.csv.gz`, `results/density-checks.csv.gz` and
`results/composition-checks.csv.gz` retain the numerical checks. Sources and
tolerances are in `results/retained-values.csv.gz`; input hashes and comparator
conditions are retained beside them. The earliest key-solve references and
older F1 declarations have the limits stated in the notebook. A missing
work-count event is separate from an unavailable numerical evaluation.
`results/work-count-checks.csv` compares all key-solve points and F1 counters
between October 17 and current MUMPS, and between current MUMPS and MA57.
`results/timing-comparisons.csv` retains the absolute differences and A=A ranges;
a timing claim requires a difference beyond the range. `results/competitors-578.csv`
supplements the #533 references with matching #578 comparator medians.

Preparation and measurement never compile wheels. The declaration and
notebook-number checks run against the retained current wheel; the requested
`bash tools/check-affected.sh origin/main` gate also owns ordinary checkout setup.
The lead owns independent delivered review and merge; this task opens a ready PR.
