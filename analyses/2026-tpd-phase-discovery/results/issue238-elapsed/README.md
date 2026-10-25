# Declared-start elapsed allowance and ordinary-state comparison

At runtime source `875dbd5a` (main `c3100e94` integrated), Born case 1 uses
313.15 K, free pressure, c_dielectric=0, the retained nine-species liquid and
three-neutral vapor, and ordered water-reaction extents 0.005, 0.02, 0.05,
0.1, 0.2 and 0.3. The current wheel SHA-256 is
`82d682b28a03f7516f27a34f994d084589e168498ada7de7f111435382c55178`.
The 1 s aggregate and 0.35 s attempt allowances stop after three starts and
24 optimizer iterations, at 1.014885997 s. No state is selected. The last
three starts consume no work. The two complete recorded row vectors have
maximum absolute values 185.47568909850332 and 221.4320092906285. The third
attempt expires without a complete row vector. All attempts retain native
`INVALID_NUMBER_DETECTED` and separate elapsed-budget causes; all state
availability fields are false. Declared phase count is diagnostic topology,
not an available state. Aggregate excess is 0.014885997 s; the exposed maximum
aggregate/per-attempt overrun is 0.017649549 s.

[Attempt work](case1-attempt-work.csv), [aggregate result](case1-summary.json),
per-attempt JSON, row CSVs and iteration CSVs retain the measured diagnostics.
[Comparison](comparison.json) identifies both wheels and source scripts.
Input owners remain `analyses/2026-born-unavailable-states/inputs/`, the
retained continuation anchors and `scripts/starts.py`; no thermodynamic
calculation or restart selection runs in Python.

Ordinary ionic liquid–liquid and reacting MEA flash states and work agree
exactly with main. All 1,176 numerical rows reuse
[`issue238-cost/main-states.csv`](../issue238-cost/main-states.csv) by SHA-256;
only the six additional text [messages](messages.csv) are retained here.
The counts remain 16/49/17 and 43/57/44 (iterations/evaluations/Hessian calls).
[Before-correction differences](before-counting-fix.csv) expose one extra
reported evaluation per ordinary solve on integrated source `f07b465d`.
Main and the candidate had each incremented the same initial-polish count.
Source `875dbd5a` removes the duplicate increment, preserving every
preparation and row calculation. This exact work comparison fails before
that correction and passes afterward.

The [main samples](main-cost.csv) and [candidate samples](branch-cost.csv)
retain three elapsed measurements per case. Best times are 0.053393316 /
0.050402749 s for the ionic point and 0.091839967 / 0.087094401 s for MEA.
Differences of -5.6010% and -5.1672% do not establish causal overhead or a
speedup on a shared host. Main wheel SHA-256 is
`8640b283926540aad7ebecdce900aa331edbc5c7196bdedb4f80daf9e67a7b31`.
Both provenance files retain literal problem inputs and input-file hashes.

Reproduction holds `/tmp/epc-shared-env.lock` before `agent-heavy`, installs
and verifies the relevant non-editable wheel, and restores/verifies the
shared reference wheel on every exit. Commands under the candidate lease:

```bash
.venv/bin/python analyses/2026-born-unavailable-states/scripts/starts.py --case 1 --seconds 1.0 --attempt-seconds 0.35 --output-dir build/issue238-elapsed-budget/timed
.venv/bin/python analyses/2026-tpd-phase-discovery/scripts/declared_start_cost.py --label branch --source-commit 875dbd5a --wheel build/environment-wheel/epcsaft-0.2.0.dev0-cp313-cp313-linux_x86_64.whl --output build/issue238-elapsed-budget/cost --retain-failed
bash tools/check-all.sh
bash tools/run-slice.sh khudaida
```

The timed replay is retained analysis, not a timing-sensitive persistent test.
Its exact timings and iteration counts may vary. The single observation
supports shared allowance exhaustion and state rejection after bounded inner
work; it gives no strict real-time guarantee or universal overrun bound.
ALG-EQUILIBRIUM, ALG-DENSITY and ALG-TERMINATION remain unchanged and were
revalidated against these measured results. The full gate passes 119 Python
tests (one deselected), six native targets, ruff, mypy, architecture and
analysis layout. Khudaida accepts all eight points at 293.15 K / 100000 Pa,
maximum original row 5.573319583618286e-14 and component balance
3.608224830031759e-16 mol/s. Complete logs and slice values are retained in
`build/issue238-elapsed-budget/`.

The [research notebook](../../notebook.qmd) interprets these discrete checks.
This is numerical verification, not physical validation, exhaustive root
finding, a statement that Born equilibrium is absent, or global stability.
Uyan chemistry/reference reconstruction remains with #234. Independent
delivered review of this evidence and the counting correction remains with
the orchestrator; manuscript use has not been promoted.
