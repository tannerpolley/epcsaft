# Reaction-equilibrium constant fitting

This Engine-owned study qualifies first derivatives and noise-free recovery
of named scalar ln K and the stored a, b, c, d coefficients of
`ReactionLogPolynomial`. [Issue #187](https://github.com/tannerpolley/ePC-SAFT/issues/187)
owns the design and owner amendments. The [notebook](notebook.qmd) distinguishes the current analytic exchange
test from the retained historical comparator; it does not rerun the model. No downstream calculation is
part of this study.

The numerical inputs live with their falsifiers:

- `engine/tests/native/regression_reaction_checks.cpp`: `check_reaction_derivatives`,
  `check_reaction_recovery`, `check_reaction_rejections` and their fixture helpers.
  Recovery includes the numerical check that changed resolved ln K invalidates
  validation preparation.
- `engine/tests/regression/test_reaction_constants.py`: the two-reaction
  analytic identical-species exchange system and the native coefficient fit.
  The older notebook recovery values retain their historical manufactured
  inputs; they are not equal-number expectations for the replacement case.
- [cold_start.py](scripts/cold_start.py): the original cold-start finding at
  the specified J1 free-temperature state. It preserves Engine defaults.
  The [density diagnosis and correction](../2026-cold-vapor-branch/README.md)
  retain #204's recovered state without replacing the original input record.

Build and install the non-editable wheel in the sole environment before
running the numerical programs. Every heavy step uses the shared lock, a
time cap and one OpenMP/BLAS thread. Run these from the repository root:

```bash
env EPCSAFT_BUILD_JOBS=2 bash tools/build-engine.sh --install-wheel
```

After preparation, read `build/environment-wheel/meson-build-dir` for the
native executable directory. Run its `regression_checks R2` program with the
same lock and environment and a 600 s cap, capturing output as
`build/issue-187/native.log`. The native equilibrium and boundary programs are
also required; `bash tools/check-affected.sh fe00aa6a` selects them against the
integrated main baseline. The complete fast gate is
`bash tools/check-all.sh`, with a 900 s cap and the same lock and environment.

The exchange calculation is deliberately slow-marked and must run explicitly:

```bash
env EPCSAFT_BUILD_JOBS=2 .venv/bin/python -m pytest -q -s -o addopts='' engine/tests/regression/test_reaction_constants.py > build/issue-187/exchange.log
env EPCSAFT_BUILD_JOBS=2 .venv/bin/python analyses/2026-reaction-constant-fitting/scripts/cold_start.py
env EPCSAFT_BUILD_JOBS=2 .venv/bin/python analyses/2026-reaction-constant-fitting/scripts/retain.py build/issue-187/native.log build/issue-187/exchange.log
env EPCSAFT_BUILD_JOBS=2 .venv/bin/python analyses/2026-reaction-constant-fitting/scripts/plot.py
```

`retain.py` transcribes values and diagnostics already emitted by the tests.
`plot.py` reads only the retained CSV. The native test owns the tolerance and
chart; the figure does not calculate missing model values. The retained
EOS step probe diagnoses truncation at the original unit-scale coordinate
and explains the accepted 0.05 interaction and 0.5 Å diameter scales.
It is not an additional fitting method or acceptance requirement.

The #557 combined coefficient/enthalpy verification uses this checkout's
prepared editable Engine environment. It consumes the frozen
[Kim–Svendsen declaration](inputs/kim-svendsen-interval.json) and the existing
F1 packet, without reading an application checkout:

```bash
nice -n 10 .venv/bin/python analyses/2026-reaction-constant-fitting/scripts/verify_pair.py build/issue-557
```

The native `reaction-enthalpy-shift`, `reaction-derivatives` and `ionic-addition`
cases provide the applied-law invariants, joint recovery, duplicate refusals
and independent binary128 endpoint comparison. The notebook reads retained
values from `results/coefficient-enthalpy-pair/`; it never reruns the model.
The table retains raw steps without pass/fail and applies the issue's binding
Richardson comparison to the native columns.
