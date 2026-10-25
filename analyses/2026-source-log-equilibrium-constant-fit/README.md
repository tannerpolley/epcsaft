# Direct source-basis ln K fitting

[Notebook](notebook.qmd): exact coefficient identities, named rejection and
manufactured joint a/b recovery for Engine issue #211. Numerical verification
only; no application chemistry, physical validation or parameter adoption.
The native input owner is `engine/tests/native/regression_reaction_checks.cpp:check_log_constant`.
The installed Python check in `test_reaction_constants.py` recovers a/b and
checks absolute-uncertainty weights through the optimizer Jacobian.

Run from the checkout root, after building and installing its non-editable
wheel with the sole environment. Every heavy command uses the shared lock,
finite timeout, low priority and requested thread caps:

```bash
env EPCSAFT_BUILD_JOBS=2 bash tools/build-engine.sh --install-wheel
env EPCSAFT_BUILD_JOBS=2 bash engine/tools/run-native-tests.sh regression_fit
```

To retain the direct numerical operands, run the configured executable under
the same controls and pass its live output to `tee`:

```bash
build_dir="$(<build/environment-wheel/meson-build-dir)"
env EPCSAFT_BUILD_JOBS=2 "$build_dir/regression_checks" ln-k | tee build/issue211-ln-k.log
.venv/bin/python analyses/2026-source-log-equilibrium-constant-fit/scripts/retain.py
env EPCSAFT_BUILD_JOBS=2 bash tools/check-all.sh
```

`results/comparisons.csv` owns numerical comparison operands, tolerances and
units. `results/summary.json` owns fitted convergence diagnostics, named
rejections, asserted covariance rank and the wheel SHA-256. The rank is two
because the native covariance check requires equality to the two fitted
coordinates; it is not inferred from optimizer convergence.
