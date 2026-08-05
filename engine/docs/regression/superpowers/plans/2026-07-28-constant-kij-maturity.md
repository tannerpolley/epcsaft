# Constant `k_ij` Maturity Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Harden the existing one-active-pair constant-`k_ij` fitter against
runaway solves and contract regressions, and add a second direct-experimental
neutral VLE reference using all 22 May 2015 methane/propane rows.

**Architecture:** Keep `RegressionProblem`, `RegressionResult`,
`fit_parameters`, the existing native module, and the existing
`neutral_binary_phase_kij_v1` callback as the sole owners. Add one declared
solver-time budget and Ceres evaluation counts to the current result. Retain
the May methane/propane table as compact source evidence and exercise it
through a caller-built public Provider parameter bundle; do not add a catalog
record, fitted-parameter overlay, or reference-case runtime API.

**Tech Stack:** Python 3.13, C++17, Ceres, Eigen, pytest, CMake.

## Global Constraints

- Fit exactly one constant unordered-pair `k_ij` per problem.
- Do not add temperature-dependent `k_ij(T)`, simultaneous multi-pair fitting,
  LLE, associating VLE, reactive equilibrium, or numerical production
  derivatives.
- Provider remains the sole EOS and exact-derivative owner.
- May et al. Table 6 is direct experimental VLE evidence. The fitted
  methane/propane value is an in-sample Regression result, not a parameter
  reported by May and not Provider-catalog authority.
- Preserve existing methane/ethane, aqueous MIAC, and ion-solvation results.
- Keep solver, numerical, workflow, physical, scientific, and predictive
  statuses separate.

---

### Task 1: Bound and diagnose every general Ceres solve

**Files:**
- Modify: `src/epcsaft_regression/parameter_regression.py`
- Modify: `src/epcsaft_regression/native/general_fit.cpp`
- Modify: `tests/test_parameter_regression_contracts.py`
- Modify: `tests/test_parameter_regression.py`

**Interfaces:**
- Consumes: `RegressionProblem.maximum_solver_time_seconds: float`.
- Produces:
  `RegressionResult.residual_evaluation_count: int` and
  `RegressionResult.jacobian_evaluation_count: int`.

- [ ] **Step 1: Add failing Python contract tests**

  Require a positive finite `maximum_solver_time_seconds`, verify it occupies
  one stable field in `_native_payload`, and reject zero, negative, infinite,
  and nonnumeric values.

- [ ] **Step 2: Run the focused tests and retain the expected failure**

  Run:

  ```bash
  .venv/bin/pytest -q \
    tests/test_parameter_regression_contracts.py \
    tests/test_parameter_regression.py::test_general_kij_fit_reports_rank_confirmation_and_partition_isolation
  ```

- [ ] **Step 3: Extend the existing payload and result in place**

  Add the time budget to `RegressionProblem`, validate it with
  `_require_finite(..., positive=True)`, serialize it directly after
  `maximum_iterations`, and add the two integer Ceres counts to
  `RegressionResult`. Update every repository-owned problem constructor with
  an explicit budget; use `30.0 s` for compact general-engine tests and
  `180.0 s` for the retained Figiel aqueous reference problems.

- [ ] **Step 4: Bind the budget and counts in the existing native owner**

  Extend `Payload` and `parse_payload` by one field, assign
  `options.max_solver_time_in_seconds`, and append
  `summary.num_residual_evaluations` and
  `summary.num_jacobian_evaluations` to the existing returned tuple. Do not
  add a timer thread, callback wrapper, backend selector, or second result.

- [ ] **Step 5: Prove bounded termination and useful accounting**

  Add one native general-engine test with an extremely small declared budget
  that returns promptly without being marked numerically converged. For the
  ordinary methane/ethane case, require positive finite evaluation counts and
  preserve all existing fit diagnostics.

- [ ] **Step 6: Run the focused contract/native matrix**

  ```bash
  .venv/bin/pytest -q \
    tests/test_parameter_regression_contracts.py \
    tests/test_parameter_regression.py -m "not campaign"
  ```

- [ ] **Step 7: Commit the bounded-solver checkpoint**

  ```bash
  git add \
    src/epcsaft_regression/parameter_regression.py \
    src/epcsaft_regression/native/general_fit.cpp \
    tests/test_parameter_regression_contracts.py \
    tests/test_parameter_regression.py
  git commit -m "Bound and diagnose general parameter solves"
  ```

### Task 2: Add the May methane/propane source-backed reference

**Files:**
- Create: `evidence/may-2015-methane-propane-vle.csv`
- Modify: `tests/test_parameter_regression.py`
- Modify: `docs/science/general-parameter-regression.md`
- Modify: `ARCHITECTURE.yaml`

**Interfaces:**
- Consumes: the public Provider `ParameterBundle.from_records`, the installed
  Gross--Sadowski methane and propane pure records, and
  `neutral_binary_phase_kij_v1`.
- Produces: one campaign-only `88 x 45` constant-`k_ij` reference fit over 22
  direct experimental rows.

- [ ] **Step 1: Add the exact source table**

  Transcribe all 22 Table 6 rows from May et al., *J. Chem. Eng. Data* 2015,
  60, 3606--3620, DOI `10.1021/acs.jced.5b00610`. Retain `T/K`, `p/kPa`,
  methane liquid fraction `x1`, its reported standard/combined uncertainties,
  propane vapor fraction `y3`, and its reported standard/combined
  uncertainties. Record the local source PDF SHA-256
  `53fd1bdd55dc6807ec76cf88626438d8dfceb3ec09149d4405ea36cfbe6b842a`.

- [ ] **Step 2: Add a failing source-identity and transformation test**

  Require 22 unique rows, the exact CSV SHA-256, exact numeric transcription,
  `P_Pa = 1000 p_kPa`, and `y_methane = 1 - y_propane`. Source uncertainties
  remain observations and are not converted into fit weights or acceptance
  cutoffs.

- [ ] **Step 3: Add the public user-bundle campaign fixture**

  Build a test-only `purpose="user-provided"` methane/propane bundle from the
  two installed Gross--Sadowski pure records plus an explicit zero active-pair
  initialization record. Select component order `("methane", "propane")`.
  Do not add or alter an installed Provider catalog.

  For each VLE row use pressure scaling by observed pressure,
  dimensionless `mu/RT` scales, liquid start/origin `4.0e-5 m3/mol`, ideal-gas
  vapor start/origin `RT/P`, liquid bounds `[3.0e-5, 2.0e-4] m3/mol`, and vapor
  bounds `[5.0e-5, 1.0e-2] m3/mol`. These starts are preflight inputs, not
  fitted observations or copied EOS equations.

- [ ] **Step 4: Prove the exact derivative before fitting**

  At a nonzero trial direction, compare the native exact `88 x 45` Jacobian
  product with a centered directional residual difference. Use the existing
  derivative-check tolerances already justified for the lifted pair engine.

- [ ] **Step 5: Prove the full reference fit**

  Require Ceres convergence, numerical convergence, complete row accounting,
  full rank `45`, projected parameter rank `1`, no active bound, three-start
  confirmation, exact Provider Hessian status for all rows, and the preflight
  anchor `k_ij = 0.0038919335722629794` within a tolerance justified by a
  repeat clean installed-artifact run. Record evaluation counts and wall time
  descriptively; make no portable wall-time assertion.

- [ ] **Step 6: Document the scientific meaning**

  Add the May methane/propane result as direct-experimental in-sample transfer
  evidence. State explicitly that May supplies the VLE observations, not the
  fitted PC-SAFT `k_ij`, and that the result is neither predictive evidence
  nor Provider-catalog authority.

- [ ] **Step 7: Run and commit the second-reference checkpoint**

  ```bash
  .venv/bin/pytest -q -m campaign \
    tests/test_parameter_regression.py::test_all_may_methane_propane_rows_reproduce_the_general_kij_reference_fit
  git add \
    evidence/may-2015-methane-propane-vle.csv \
    tests/test_parameter_regression.py \
    docs/science/general-parameter-regression.md \
    ARCHITECTURE.yaml
  git commit -m "Add methane propane kij reference campaign"
  ```

### Task 3: Complete the constant-`k_ij` robustness matrix

**Files:**
- Modify: `tests/test_parameter_regression.py`
- Modify: `tests/test_parameter_regression_contracts.py`
- Modify: `docs/science/general-parameter-regression.md`
- Modify: `ARCHITECTURE.yaml`

**Interfaces:**
- Consumes: existing neutral VLE, Figiel/Hamer--Wu NaBr aqueous MIAC, and
  organic-ion solvation problem builders.
- Produces: one bounded evidence matrix for all three admitted constant-`k_ij`
  observation domains.

- [ ] **Step 1: Add neutral invariance tests**

  Fit the methane/propane rows in source order and reversed order and require
  the same fitted parameter, cost, ranks, and row accounting within
  binary64/Ceres repeatability. Replace the parameter identity with
  `PairParameterIdentity("propane", "methane")` and require the same canonical
  identity and fit.

- [ ] **Step 2: Prove direct-observable context isolation**

  For the 21-row NaBr fit, retain copies of every row's three fixed pair
  values, fit water--sodium, and require that the immutable observations and
  two inactive values are unchanged. Require the result to contain exactly
  one fitted coordinate, rank one, no active bound, complete accounting, and
  exact first-derivative status.

- [ ] **Step 3: Complete the negative controls without duplicates**

  Reuse existing unsupported fingerprint/domain and Provider-row-failure
  tests. Extend the existing monkeypatched terminal-diagnostic test to prove
  projected rank zero and over-limit conditioning both prevent numerical
  convergence. Keep one active-bound diagnostic test. Do not duplicate
  contract failures already owned by `test_parameter_regression_contracts.py`.

- [ ] **Step 4: Run the complete maturity matrix**

  ```bash
  .venv/bin/pytest -q tests/test_parameter_regression_contracts.py
  .venv/bin/pytest -q tests/test_parameter_regression.py -m "not campaign"
  .venv/bin/pytest -q -m campaign tests/test_parameter_regression.py
  ```

- [ ] **Step 5: Run installed-artifact and package checks**

  Build one wheel with the repository's existing CMake/Python packaging path,
  install it into an isolated environment with the exact Provider wheel used
  for the package tests, rerun the same three-domain matrix, and audit wheel
  contents, imports, linkage, target count, and absence of evidence/runtime
  data leakage.

- [ ] **Step 6: Reconcile status and commit**

  Record the two neutral reference shapes (`68 x 35` methane/ethane and
  `88 x 45` methane/propane), the aqueous and solvation sentinels, solver
  budgeting, evaluation accounting, and explicit exclusions. Do not change
  scientific or predictive authority.

  ```bash
  git add \
    tests/test_parameter_regression.py \
    tests/test_parameter_regression_contracts.py \
    docs/science/general-parameter-regression.md \
    ARCHITECTURE.yaml
  git commit -m "Harden constant kij regression evidence"
  ```

### Task 4: Final scientific and minimality verification

**Files:**
- Review only: all files changed by Tasks 1--3.

**Interfaces:**
- Consumes: the committed implementation and evidence.
- Produces: one reviewable clean Regression subject.

- [ ] **Step 1: Recompute evidence identities and anchors**

  Verify the source PDF/CSV hashes, canonical dataset hashes, fitted values,
  derivative checks, ranks, bounds, confirmation deltas, row counts, and
  evaluation counts from the public installed-artifact workflow.

- [ ] **Step 2: Run the isolated full Regression suite**

  ```bash
  .venv/bin/pytest -q
  .venv/bin/pytest -q -m campaign
  ```

- [ ] **Step 3: Run code-surface and cleanup review**

  Confirm there is still one general native target, one public fit function,
  one result family, no catalog mutation, no new runtime dataset registry, no
  temperature-dependent coordinate, and no retained task-owned build/process
  artifact.

- [ ] **Step 4: Run repository cleanup and inspect the final diff**

  ```bash
  bash "$HOME/.codex/hooks/codex-cleanup.sh" --repo-root .
  git diff --check
  git status --short --branch
  ```

- [ ] **Step 5: Commit only if final reconciliation changed tracked files**

  Use one bounded documentation/status commit; otherwise retain the three
  independently green checkpoints above.
