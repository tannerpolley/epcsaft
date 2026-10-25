# Schick 2023 CO₂ + THF Table S10

Current native fixed published-tuple replay of every Table S10 row using full
binary vapor–liquid equilibrium. The [research notebook](notebook.qmd) separates
numerical verification against the paper model from physical comparison against
Knez's measurements as transcribed by Schick. Non-reproduction, numerical
qualification and source-rounding limits remain explicit in the retained tables.
Delivered-review corrections are applied; every row is numerically qualified.
Packet validation and parameter-index generation pass; paper discrepancies remain.

From this checkout, after preparing its own editable Engine environment:

```bash
nice -n 10 env -u EPCSAFT_HSL_LIBRARY EPCSAFT_BUILD_JOBS=6 .venv/bin/python analyses/reproduction/2023-schick/scripts/run.py
nice -n 10 env MPLBACKEND=Agg .venv/bin/python analyses/reproduction/2023-schick/scripts/render.py
```

The producer reads `inputs/parameters.json` and the source CSV files. The new
`schick-2023-co2-thf/1` packet owns the published parameters; its exact
fingerprint and byte-for-byte materialization are recorded in
[`inputs/source/materialization.json`](inputs/source/materialization.json).
No fitting or Python thermodynamic iteration is performed.

Rendering reads retained results only. Render without executing scientific code:

```bash
mkdir -p build/issue-727/quarto-runtime
cd analyses/reproduction/2023-schick
nice -n 10 env TMPDIR="$(git rev-parse --show-toplevel)/build/issue-727/quarto-runtime" bash "${CSE_RENDER:?Set CSE_RENDER to the installed CSE Quarto render.sh}" notebook.qmd
mv notebook.html ../../../build/issue-727/notebook.html
```

The wrapper is provided by the installed CSE plugin. The standalone HTML is
retained in task evidence; no model is rerun or site configuration changed.

The old pure-CO₂-reservoir reconstruction and diagnostic fit remain unchanged
under `results/historical/`. Their boundary differs from full binary VLE;
the historical absence of two reservoir roots is not a physical VLE claim.
