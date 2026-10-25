# CppAD-backed IDAES ASL adapter

This analysis exercises the installed Engine wheel through its compiled ASL
library and IDAES's ASL-enabled Ipopt. It is bounded to the neutral propane
packet, a steady isothermal throughput unit, and the sampled pressure branch.
[Issue #10](https://github.com/tannerpolley/ePC-SAFT/issues/10) owns the
accepted design. No association, electrolyte, equilibrium, global phase-stability
or downstream application qualification follows from this demonstration.

New runs use the byte-identical physical coefficients in the version-3 propane
packet, which explicitly records zero charge. The retained September 16 results
remain bound to version 2 and their recorded wheel; this input migration does
not requalify that campaign on the current wheel.

```bash
bash tools/build-engine.sh --install-wheel
.venv/bin/python analyses/2026-idaes-asl-adapter/scripts/run.py \
  --output-dir build/idaes-adapter-candidate \
  --ipopt "${IPOPT_EXECUTABLE:?Set IPOPT_EXECUTABLE to the Ipopt executable}"
```

Finish environment preparation before running the campaign: setup can temporarily
remove the installed Engine while refreshing it. The executable above was
identified locally as Ipopt 3.14.19 / ASL 20231111; another machine must supply an
explicit ASL-enabled executable satisfying the runner's option preflight.

`results/derivatives.csv` owns native/ASL comparisons and Pa/bar expression
chains; `lagrangian.csv` owns whole-NLP curvature; `process.csv` owns the three
starts and changed-cap solve; `references.csv` owns independent equation and
dilute-limit references; `domain.csv` owns the 45 conditioning-grid and nine
derivative-state evaluations. `cost.csv` retains real cold/warm timings and observed
tape recordings. `summary.json` retains immutable wheel/library/parameter
identities, solver options, numerical aggregates, failures and claim limits.
`results/solver/` retains the solver-consumed NL and its row/column maps and logs.
Atomic CSV replacement and scalar checkpoints preserve completed rows during
bounded child execution and cancellation.

Campaigns require a new or empty output directory and never replace existing
evidence. The retained `results/` directory is updated only from a fully accepted
candidate.

Render the retained-result figures without rerunning the model:

```bash
.venv/bin/python analyses/2026-idaes-asl-adapter/scripts/plot_results.py
```

The pressure/conditioning, derivative/cost and process plots read only the
retained CSV/JSON files. The flowsheet and adapter-path diagrams describe the
implemented model and compiled call boundary.

`results/historical/` preserves superseded loading and sweep-cost receipts under
their original source and wheel identities. They are not current-main evidence
and are excluded from the notebook's current result tables.

## Qualified baseline - 2026-09-21

The `idaes/main` branch is rebased onto merged main commit `2f44dfc1`,
including the completed Bindings and parameter-regression work. Source commit
`ec5bfdba` produced the
retained campaign from a clean tree and the one installed non-editable wheel.
The campaign passed in 21.80 s with 72 scalar
derivative ladders, 216 whole-NLP directional ladders, 27 explicit Hessian
matrices, and four optimal process cases. Maximum KKT stationarity was
`1.72e-16`; maximum complementarity was `9.10e-11`.

The domain, derivative, whole-NLP and independent-reference tables are identical
to the prior accepted campaign after line-ending normalization. Concurrent
fixed-feed and interruption checks shared the host during this replay; its
timings do not establish a performance change. The full installed-wheel gate
passed before the numerical campaigns. Setup and wheel installation must run
separately from campaigns because setup temporarily uninstalls the Engine.

`domain.csv` retains 54 passing evaluations: 45 points across 385-415 K and
5-110 mol/m3 plus the nine derivative states. The minimum sampled pressure
slope divided by RT was `0.9498`. Whole-NLP recovery and Hessian checks ran on
the exact prewritten `process-0.nl` before any solve, and every solve retained
the SHA-256 of the NL file passed to Ipopt.

The runner refuses nonempty output directories, validates the materialized
packet file inventory instead of independently redefining its fingerprint, and writes
campaign identity before numerical work. Atomic reference, domain, derivative,
and whole-NLP checkpoints remain reviewable after failure or interruption.

## Notebook

`notebook.qmd` reads retained values only. It combines exact scalar tables with
plots of the retained grid, derivative ladders, timing batches and process cases;
no visual interpolates beyond the sampled model domain.
For local rendering, install the optional Quarto Python kernel support in the
same environment after the final setup/check command:

```bash
uv pip install --python .venv/bin/python nbclient ipykernel pyyaml
mkdir -p build/adapter-notebook
QUARTO_PYTHON="$PWD/.venv/bin/python" quarto render \
  analyses/2026-idaes-asl-adapter/notebook.qmd --execute --execute-daemon 0 \
  --to html --output - > build/adapter-notebook/notebook.html
```

Rendering performs no thermodynamic calculation. The working notebook is not
promoted for manuscript use.

## Verification and handoff

The retained summary records wheel SHA-256 `9f73aeb4...e606`, ASL library
SHA-256 `d4bf537f...57c6`, packet fingerprint `4593b9db...326e2`, Ipopt
3.14.19, and ASL 20231111. The full installed-wheel gate and independent reviews
must remain green before the bounded qualification is closed. This issue-owned
analysis may use `epcsaft._core` only for native-versus-ASL comparison. The
separately tested public B1/B2 state-action surface does not alter this ASL path.

No equilibrium, phase-stability, association, electrolyte, or downstream
application claim follows from this neutral propane result. Those require their
own scientific slices and evidence.

## Manual interruption

Parent SIGTERM previously left a new-session numerical child alive and its
summary marked running. The initial probe preserved nine real scalar rows
before demonstrating that failure. One shared child-process helper now owns process
cleanup for the numerical, solver and timing paths. SIGTERM raises
`KeyboardInterrupt` so it propagates through the selector; `InterruptedError`
is swallowed by Python's `EpollSelector.select` and was rejected during diagnosis.

Both SIGTERM and Ctrl+C terminate the child group, preserve the completed scalar
prefix (nine rows in each current replay), retain campaign identity, and atomically
record `stage: interrupted`.
The probe pauses a real numerical child after its checkpoint; it changes no
thermodynamic input, formula, or tolerance.
The runner's SIGTERM handler deliberately raises `KeyboardInterrupt`; CPython
therefore reports controller return code `-SIGINT` (`-2`) for both probes, while
the recorded `signal` field distinguishes the injected signal.

Run the focused check with a new disposable output directory:

```bash
.venv/bin/python analyses/2026-idaes-asl-adapter/scripts/check_interruption.py \
  --output-dir build/adapter-interruption \
  --ipopt "${IPOPT_EXECUTABLE:?Set IPOPT_EXECUTABLE to the Ipopt executable}"
```

`results/interruption-before/` retains the reproduced defect and
`interruption-after/` retains the rejected `InterruptedError` diagnostic. Current
probe output stays in the chosen disposable directory unless explicitly selected
for retention.

## Native pressure-cap reference

The [conversion comparison](results/batch7-conversion-summary.json) verifies
native propane pressure-cap references and derivative transport. Both complete
process campaigns stop when the existing cancellation probe fails to observe
the short-lived live solver; no cost stage is reached. This conversion does
not establish full process-campaign acceptance. IDAES/Pyomo process solves
retain their exemption from the native calculation rule.

## External solver rerun (#720)

The locally selected `~/.idaes/bin/ipopt` reports Ipopt 3.13.2 / ASL 20190605
and fails the required-option preflight for `acceptable_iter`, `derivative_test`,
`derivative_test_perturbation`, `derivative_test_tol` and
`point_perturbation_radius`. The explicit Ipopt 3.14.19 / ASL 20231111 executable
in the command above supports every required option and completes this campaign
on `ca719012`. The producer now also sets outer `linear_solver=mumps` explicitly.
Keep that executable and its environment outside the repository, following the
[build-time-only IDAES setup](../../docs/site/contributing.md).
HSL is never copied into the Engine or retained analysis output.

[The new summary](results/issue-720/summary.json),
[domain rows](results/domain.csv),
[independent references](results/references.csv), and
[process rows](results/issue-720/process.csv) retain the rerun.
Full derivative tables, NL files and solver logs remain under
`build/issue-720/idaes-explicit-mumps/`; its solver logs confirm MUMPS 5.8.2.
All existing numerical criteria pass, with no reference replacements.
