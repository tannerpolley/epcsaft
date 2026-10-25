# Trace solvent activity

The [research notebook](notebook.qmd) records the common-pressure solvent-activity
calculation and its numerical qualification. The current design and delivery live
in [#214](https://github.com/tannerpolley/ePC-SAFT/issues/214). Old resolved values
remain unchanged in `analyses/2026-divalent-chloride-activity/`; the Held 2008 packet
inputs (versions 3 and 4) are read from `data/packets/held-cameretti-sadowski-2008/`,
whose files are byte-identical to the leaf copies #341 removed (the notebook cites the
copies' former path).

Using this checkout's prepared `.venv`, replay the numerical falsifiers:

```bash
nice -n 10 env EPCSAFT_BUILD_JOBS=6 .venv/bin/python analyses/2026-trace-solvent-activity/scripts/verify.py
```

The replay writes resolved old/new comparison rows, trace values and typed
reporting outcomes, and wheel/input/script identities. `scripts/plot.py` reads
these retained rows without recalculating thermodynamics.

`results/diagnosis.csv` and its identity preserve the failing baseline at
`3f5ee19f`. `scripts/diagnose.py` reproduces that diagnosis only with that baseline
wheel; running it on the corrected wheel would replace the baseline evidence.
The native integral and derivative-error probes, sources and identities are
retained investigations of an unselected alternative. They were evaluated
against the baseline native archive and are not implemented Engine capabilities.

## Current producer diagnosis (#750)

[The current diagnosis](results/current-main/diagnosis.md) records the unchanged failure,
the responsible Engine changes, retained comparisons, repaired reruns and remaining
limits. New results are separate from the original retained results.
