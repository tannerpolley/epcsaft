# Head-to-head comparisons with FeOs and Clapeyron

Four comparisons for the solver capability matrix, run on 2026-10-08 with the Engine main-current wheel:

1. Fresh Engine C1–C4 solve, gradient and fixed-state timings, with every timed return checked against the retained Engine values under their existing criteria.
2. Khudaida electrolyte LLE against Clapeyron's electroneutral Michelsen flash, with constant and composition-dependent permittivity, including the Clapeyron Debye–Hückel derivative defect on the `debye_length` path.
3. Exact solved-state sensitivities (C2 bubble point, C4 LLE) against Clapeyron implicit differentiation, FeOs automatic differentiation and Engine central-difference checks.
4. Born term and concentration-dependent permittivity at fixed NaCl–water states.

All are numerical verification on declared inputs, not physical validation. The [research notebook](notebook.qmd) computes every reported number from the retained results and states the limits. It has not been independently reviewed.

This analysis reuses the [competitor cross-check](../2026-competitor-cross-check/README.md): its parameter cases and constants (`inputs/cases.json`), its Engine/FeOs/Clapeyron case builders (`scripts/engine.py`, `scripts/feos.py`, `scripts/clapeyron.jl`), its comparator environments and timing protocol ([reproduction](../2026-competitor-cross-check/README.md#reproduction)), and its retained results as reference values.

## Reproduction

Install the main-current wheel into this checkout's `.venv` and verify its RECORD hashes, then run from the repository root:

```bash
A=analyses/2026-head-to-head-comparisons X=analyses/2026-competitor-cross-check C=~/.cache/epcsaft-comparators/222
# 1. Engine timing: the cross-check's engine.py collector and three timing processes, output build/h2h/Engine-main*.json
# 2-4. Engine side and the comparator inputs, then FeOs
.venv/bin/python $A/scripts/h2h_engine.py build/h2h/engine-h2h.json
$C/feos/bin/python $A/scripts/h2h_feos.py build/h2h/engine-h2h.json build/h2h/feos-h2h.json
# Clapeyron: one Julia process; each line sent to the FIFO is a task file to include, "quit" ends it
mkfifo build/h2h/julia.fifo
JULIA_DEPOT_PATH=$C/depot systemd-run --user --scope -p MemoryMax=6G \
  $C/julia-1.13.1/bin/julia --startup-file=no --project=$X/inputs $A/scripts/h2h_clapeyron_server.jl $PWD/build/h2h/julia.fifo &
echo $PWD/$A/scripts/h2h_clapeyron_tasks.jl > build/h2h/julia.fifo; echo quit > build/h2h/julia.fifo
# Read the native outputs into the retained tables
for mode in engine-timing comparators probes; do python3 $A/scripts/head_to_head.py $mode build/h2h; done
```

The Clapeyron task file resumes from finished tasks in its output JSON, so a stopped process loses only the unfinished task. `head_to_head.py` classifies each Clapeyron flash return: an equilibrium needs a finite phase fraction on the specified pressure, and invalid returns stay in the tables with their diagnostics. The point-1 term probes in `scripts/h2h-probes/` ran once from the repository root against `build/h2h/` (Python probes on the installed Engine, Julia probes through the same FIFO); their outputs are retained in `results/h2h-probes/` and the Clapeyron probe lines in `results/h2h-raw/julia-process1.log.gz`.

## Results

- `results/h2h-engine-timing-summary.csv`, `h2h-engine-timing-statistics.csv`, `h2h-engine-timed-values.csv`; raw processes in `results/h2h-engine-raw/`.
- `results/h2h-khudaida-lle.csv`, `h2h-khudaida-convergence.csv`, `h2h-khudaida-term-checks.csv`.
- `results/h2h-sensitivities.csv`, `h2h-sensitivity-costs.csv`.
- `results/h2h-born-permittivity.csv`.
- `results/h2h-session-times.csv`; raw native outputs and Julia logs in `results/h2h-raw/`.

## Named-release rerun (#680)

The [research notebook](notebook.qmd) reads the new results under
[results/main-8298c602/](results/main-8298c602/). The earlier retained files remain
under their original identities. The new provenance records the downloaded
wheel, its verified hash, the actual installed Engine identity and the declared
producer. Unavailable portions retain their stopping states and causes.
Independent delivered review and owner promotion remain open; the manuscript
is unchanged.
