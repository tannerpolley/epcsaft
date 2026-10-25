# Native electrolyte fitting qualification

This study numerically verifies the five electrolyte observations and the exact
held-coordinate chain delivered by #219 slice A, and the explicit native
guarded minimax path in slice B. It uses the installed,
non-editable Engine wheel. It does not fit or adopt #209's evaluated experimental
data, establish a global optimum, or physically validate a parameter set.

## Inputs and conditions

Slice-A qualified wheel SHA-256: `ceb6ab8a48b671fbb48330ff1195b14e2b9372c3829edb3f4a49d268e3895707`.

All rows use 298.15 K, 100000 Pa and the liquid density branch. The parameter
source is `figiel-2025-reference-electrolytes/4`; its JSON SHA-256 is
`2504e634e555d1e28c8ddae7ea53b0fc85e208234be1156820650b2d6e1bd42d`.
`inputs/files.json` records the hashes of every consumed source table and the
separate density-temperature permittivity selection. `scripts/declarations.py`
owns the selected components, rows, coordinates and hold declarations.

F1 contains 383 rows: 156 aqueous chloride gamma/Phi rows, 156 bromide rows,
21 transfer-energy rows, three hydration energies, 29 water permittivity rows
and 18 alcohol permittivity ratios. The bromide input is #209's retained Hamer-Wu grid,
`../2026-figiel-divalent-chlorides/inputs/bromide-hamer-wu.csv`. Other
tables are read from their existing repository owners. The Hamer-Wu values
supply the molality grid; F2/F3 manufacture observations natively.

Source-defined component subsets avoid inventing an ethanol-methanol pair
absent from the packet. Each subset differentiates every applicable coordinate;
an absent component supplies no parameter column. Global F1 column order is
six cation-water/cation-Cl k values (Li, Na, K), common C, Li/Na/K Born diameters,
then Cl sigma, Debye-Huckel diameter and packing diameter. Common C ties every
charged slot. `falsifiers.csv` identifies the global observation/column and the
finite-difference step; unprescribed steps have no pass verdict.

F2's truth shifts the source k values by 0.05 (water) and 0.10 (Cl), with signs
keeping them inside [-1,1], and C by -1 inside [3,9]. Three native dependent
Born roots hold -487020, -381930 and -308890 J/mol in water, with a [1,8] Angstrom
bracket, |H| <= 1e-6 J/mol and correction/refined width <= 1e-10 Angstrom.
The 156 manufactured gamma/Phi rows use LogRatio weights
1/(0.02 ln(10))^2. Ceres receives 100 iterations, 120 seconds and stopping
tolerances 1e-15, 1e-14, 1e-14 for function, gradient and parameters. F3 adds
+0.005 decade to every manufactured NaCl gamma datum. The independent native
bisection uses the stored fitted parameters and [1,8], without the chart cache.

## Calculations and retained results

Thermodynamics, exact columns, manufacturing, held roots, Ceres fits, final
check metrics, SVD/covariance and multistart execute in C++. Finite differences
and comparison algebra appear only in verification. Python declares the
perturbations and reads native values; it drives no optimizer or Newton loop.
The reference is the Engine's admitted finite dilution-ladder endpoint. The
observed terminal changes do not prove a universal infinite-limit error bound.

- `results/summary.json`: wheel identity, inputs, fit stopping points, hold
  errors, singular values, bounds, projected gradients, timing and failures.
- `results/falsifiers.csv`: F1-F6 comparisons and the retained step studies.
- `results/reference-ladders.csv`: terminal value/column changes and molality.
- `results/fitted-observations.csv`: manufactured values and fitted predictions.
- `results/independent-roots.csv`: independently refined held roots, energy
  errors and bracket widths; `independent-reference-inputs.txt` holds their
  native input declaration, which `scripts/electrolyte_reference.cpp` reads
  (`scripts/run.py` compiles it with `tools/native-driver`).
- `results/diagnosed-failures.json`: resolved pressure-logarithm and dilute-Phi
  cancellation misses, plus the refreshed same-build native comparison rows.

The invariant transfer column required prescribed-pressure ln Z to avoid
amplifying unequal density-root residuals. Dilute Phi required the extensive
Helmholtz Euler contraction at prescribed pressure, accumulated in extended
precision. Its exact parameter columns remain the shared CppAD fugacity
contraction. Neither correction changes an equation definition or tolerance.
Native no-root trials and a manufactured Born-free H_d=0 relation refuse a
chart; the no-root trial also preserves an earlier admitted cache point.

## Reproduction

Run from the repository root, after checking branch, dirty paths and source
modification times. Keep one writer and the same prepared native archive/wheel.

```sh
env EPCSAFT_BUILD_JOBS=2 bash tools/build-engine.sh --install-wheel
.venv/bin/python engine/tests/oracle/record_native_rows.py
.venv/bin/python analyses/2026-electrolyte-parameter-gradients/scripts/run.py
bash tools/check-all.sh
bash tools/check-affected.sh HEAD
```

The notebook presents these retained numbers without running the model.
Timing compares native value/column calls with native perturbation batches;
model construction, reference admission and each owner's derivative work are
included. It is a measured implementation comparison, with no general speedup
claim. Covariance is local conditional sensitivity for independent weighted
errors and fixed hydration targets; gamma and Phi are related by Gibbs-Duhem.
The perturbed-data case reaches active bounds, so its covariance is withheld.
Profiles are deferred under D8's 35-line driver condition. Slice B uses the
[owner-approved amendment](https://github.com/tannerpolley/ePC-SAFT/pull/230)
and its measured closeout; the older 175-line allowance is superseded.



The [independent delivered review](results/review.json) returned Supported for
the corrected code and records after two result-reporting defects were fixed
and the budget authority was linked. The investigator's complete replay on the
identified wheel passes F1-F6, including the F3 check comparison at 1e-12, and
the final `check-affected.sh HEAD` passes. The reviewer did not rerun that study;
the numerical and physical limits above remain.

## Slice B: guarded minimax

The original native B1/B2/compact B3 checks identify wheel
`7e44fadae95055db943a787a4d2c624d24c8c60bd4c75c0b6db5308d2fe26890`.
The current full B3/B4 replay is identified separately in `results/minimax.json`
and the current replay section below.
The new [notebook section](notebook.qmd#guarded-minimax-with-physical-holds)
presents the retained results. `results/minimax-native.csv` owns the analytic
operands; `minimax.csv` owns full same-point and physical-hold operands;
`minimax.json` owns declarations, packet hashes, parameter/row order, budgets,
conditions and availability. The current full B3/B4 records were rerun on `8ffd627f`;
the original native analytic records retain their original wheel and review. The six-component source includes methanol at zero
fraction in aqueous target/guard rows, solely to permit the separate transfer
probe. B4 has seven independent stage-1 coordinates and all three Born holds.

The two starts attain 0.00454722134719 and 0.00454722134722 decade maxima,
with 12/11 iterations and scaled KKT norms below 1.1e-11. Permittivity max and
mean are strictly below 0.1. Stored hydration energies, exact columns and root
refinement pass their unchanged limits. Report-only transfer failure preserves
feasibility; zero-weight permittivity guards remain enforced. This is numerical
verification with manufactured observations, not a #209 refit or physical
validation. The original F1-F6 rows retain their slice-A wheel and inspected
scope; the installed-wheel preservation checks and B3 test the shared owners.
B3 is same-evaluator consistency, not independent derivative evidence. The
tighter-root replay returns the same iterate because the original roots already
meet its criteria. Full-space stationarity is reconstructed, with independent
multiplier qualification supplied by B1's analytic equations.

```sh
build/engine-meson/<prepared-build>/regression_checks R2-D7
build/engine-meson/<prepared-build>/regression_checks minimax-consistency
.venv/bin/python analyses/2026-electrolyte-parameter-gradients/scripts/minimax.py
env -u EPCSAFT_HSL_LIBRARY EPCSAFT_BUILD_JOBS=6 nice -n 10 agent-heavy .venv/bin/python analyses/2026-electrolyte-parameter-gradients/scripts/present.py
# Reproduce the density diagnosis using the retained fitted physical inputs:
env -u EPCSAFT_HSL_LIBRARY EPCSAFT_BUILD_JOBS=6 nice -n 10 agent-heavy tools/native-driver regression_checks analyses/2026-electrolyte-parameter-gradients/scripts/density_sensitivity.cpp
env -u EPCSAFT_HSL_LIBRARY EPCSAFT_BUILD_JOBS=6 LD_LIBRARY_PATH="$(realpath "$(cat build/environment-wheel/meson-build-dir)")" nice -n 10 agent-heavy build/native-drivers/density_sensitivity/driver analyses/2026-electrolyte-parameter-gradients/results/minimax-density-parameters.txt
```

`build/environment-wheel/meson-build-dir` identifies the prepared native build.
The [independent delivered review](results/minimax-review.json) returned
Supported for `310d762d` against `1fc771b0`, covering the code, final rerun values
and notebook. Claude Opus 5.5 at high effort independently reran all 261 native
comparisons and found no code defects. The subsequent documentation closeout
records the verdict and clarifies its stated evidence limits, without changing
runtime or numerical operands. The study is numerical verification only.

Density initialization uses the existing pressure closure and recomputes every
observable and exact column. The matched cold comparator uses the same source
snapshot with previous-trial density offers disabled in `fit.cpp`; it changes
no equation, root criterion, input, fitting control or observation.
`scripts/density_seed_pairs.py` runs seven alternating AB/BA complete-call
pairs for the retained full F2 fit and both B4 starts. Imports and installation
are outside the measured boundary. The independent cold same-point replay runs
separately from timing. Its declared comparison is 1e-10 relative with a 1e-12
absolute floor near zero, following the owner's equivalence rule. The osmotic
coefficient carries the floor on ln a_w as 1e-12/(M_w sum(nu) m). The
`endpoints` mode runs one complete call per workload in each environment and
compares optimizer endpoints and evaluation counts.

```sh
.venv/bin/python analyses/2026-electrolyte-parameter-gradients/scripts/density_seed_pairs.py paired <result-directory> <cold-environment> <seeded-environment>
.venv/bin/python analyses/2026-electrolyte-parameter-gradients/scripts/density_seed_pairs.py endpoints <result-directory> <cold-environment> <seeded-environment>
printf '["F2", "full", true]\n' | .venv/bin/python analyses/2026-electrolyte-parameter-gradients/scripts/density_seed_pairs.py worker <replay-directory>
```

The [density-initialization notebook section](notebook.qmd#density-initialization-in-fitting)
links the retained comparison operands, complete-call times and wheel identities.
#341 removed the round-4 diagnosis probes from `results/density-initialization/round4/diagnosis/`
and the unread raw key-solve rows; both remain at record commit `21d7a112`. The round-5
table generator is `scripts/density_round5_tables.py`.

## Current minimax replay (#757)

The full producer passes on main `8ffd627f`, using this checkout's installed
wheel SHA-256 `4e331275a44ab795d0daa18f09455acfaeaa7681372d2fd07f8e809655008bfd`.
The repaired point evaluator from #753 returns B3 values. Both B4 starts converge
in 12 and 11 iterations with the original 300 iterations / 120 s per start and
600 iterations / 240 s shared budgets. `results/minimax.csv` and
`results/minimax.json` retain the complete rerun; the native analytic B1/B2
records and their previous independent review remain on their original wheel.
This repair awaits the parent task's independent review.

The original failure compares `result.predictions[3]` from the fitted minimax
candidate with `_core.evaluate_observations(result.parameters, cs,
result.optimizer, hs, data).predictions[3]`, a fresh cold evaluation. Row 3 is
LiCl's practical osmotic coefficient at 0.1 mol/kg, 298.15 K and 100000 Pa on
the liquid branch. `compare` sets a nonzero reference's original limit to
`COMPARE_RELATIVE * abs(reference)`: `1e-12 * 0.9440422332170179` gives
`9.440422332170179e-13`. It compares two current evaluations, not a current
prediction against a stored historical prediction. The original two missed
rows, including NaCl's row 9 after row 3, are retained in
[the diagnosis](results/minimax-diagnosis.json).

The fitted prediction uses the accepted previous-trial density initialization
introduced by `27d7a2a42` and merged in `fb73b791f` (#245). Its mixture density is
55437.90350362499 mol/m³; the cold replay uses 55437.903503625 mol/m³, one ulp
apart. Both pass the unchanged pressure criterion; their pure-water density is
identical at 55340.37925215373 mol/m³. Replaying from the source parameter record
instead of the fitted record leaves the row-3 discrepancy unchanged. The fitted
record replay reconstructs identical physical parameters, ruling out a changed
Born diameter as the cause of this comparison.

A native `solvent_activity` probe, using each density as an admitted liquid
anchor with the same fitted parameters and pure-water reference, reproduces
both predictions exactly: 0.9440422332179867 and 0.9440422332170179.
[The contraction operands](results/minimax-density-operands.csv) retain the
residual Helmholtz energies, ion chemical potentials and pressures. One density
ulp changes the mixture residual Helmholtz energy by two ulps and the weighted
chemical-potential contraction by rounding-sized amounts. Accumulating in
extended precision cannot recover the information already rounded in these
binary64 operands. The resulting difference in log water activity is
3.491e-15, carried to a 9.689e-13 difference in the osmotic coefficient.

[Wheel comparisons](results/minimax-wheel-comparisons.csv) separate the changes:
the pre-#245 local wheel `1218dc78` replays this row exactly; subsequent wheels
show seeded-versus-cold rounding differences. The published `main-21e761c3`
wheel still passes the old row-3 limit; `main-cc63a135` reproduces the current
failure and endpoint. Between those snapshots, the relevant Properties
calculation change is `488451a32`, merged in `44e6b4168` (#594): ion references
use the admitted terminal extrapolation instead of the earlier ladder value.
This changes the manufactured gamma targets and fitted endpoint. It does not
change the osmotic equation. The other changed observation projections apply
only to equilibrium rows, absent here. Later sampled wheels reproduce the same
failure. The direct same-snapshot density probe identifies the immediate cause
without attributing it to a changed equation or an Engine defect.

The exact density derivative separates smooth propagation from changes in rounded
operands. With $a=A^{res}/(nRT)$, $\mu_i=\mu_i^{res}/RT$,
$D=M_w\sum_i\nu_i m$, and $x_w$ the solvent mole fraction, differentiating the
common-pressure Euler expression in
[`practical_osmotic_coefficient`](../../engine/docs/equations.md#practical_osmotic_coefficient)
at fixed temperature, composition, prescribed pressure and pure-water reference gives

$$
\frac{\partial\phi}{\partial\rho}=-\frac{1}{D}
\left[\frac{a_\rho-P/(RT\rho^2)-\sum_{i\ne w}x_i\mu_{i,\rho}}{x_w}
+\frac{1}{\rho}\right].
$$

For one total mole, the Engine's exact Helmholtz derivatives give
$a_\rho=-A_V/\rho^2$ and $\mu_{i,\rho}=-A_{n_iV}/\rho^2$, where $A$ is the
dimensionless extensive residual Helmholtz energy and $V=1/\rho$.
The retained [native sensitivity calculation](scripts/density_sensitivity.cpp)
reads [the fitted physical inputs](results/minimax-density-parameters.txt).
It gives $\partial\phi/\partial\rho=0.00030468119455730994$ m³/mol.
Multiplication by $\operatorname{ulp}(\rho)=7.275957614183426\times10^{-12}$ mol/m³
predicts a smooth change of $2.216847457437761\times10^{-15}$.
The actual binary64 change is opposite in sign and about 437 times larger;
the finite one-ulp secant is about $-0.13316$ m³/mol. Central differences
at wider steps verify the exact derivative, as retained in the diagnosis.
Thus the smooth derivative alone cannot justify loosening this comparison.
The observed discrepancy is caused by the rounded operands changing when the
root changes, followed by amplification through division by $D$.

The repair selects common cold initialization for both point evaluators.
At each fitted optimizer point, `b3_point` now calls fresh `evaluate_minimax`
and `evaluate_observations` instances with no previous density seed. Their
prediction values retain the original relative limit, with the original
absolute limit only at zero. The stored-fit minus cold predictions remain
visible in `minimax.csv` as `fitted-minus-cold` diagnostics; their different
initializations no longer carry an equality assertion. The fitted candidate's
exact columns, chart, holds and stored checks remain compared under their
original limits. Constraints, physical guards, root refinement and stationarity
checks also retain their original limits. No prediction limit is loosened.

The rerun establishes numerical verification of manufactured observations;
it supplies no experimental validation, parameter adoption, global optimum
or unique parameter estimate. Independent delivered review belongs to the
parent task.

```sh
env -u EPCSAFT_HSL_LIBRARY EPCSAFT_BUILD_JOBS=6 nice -n 10 agent-heavy .venv/bin/python analyses/2026-electrolyte-parameter-gradients/scripts/minimax.py
env -u EPCSAFT_HSL_LIBRARY EPCSAFT_BUILD_JOBS=6 nice -n 10 agent-heavy .venv/bin/python analyses/2026-electrolyte-parameter-gradients/scripts/present.py
# Reproduce the density diagnosis using the retained fitted physical inputs:
env -u EPCSAFT_HSL_LIBRARY EPCSAFT_BUILD_JOBS=6 nice -n 10 agent-heavy tools/native-driver regression_checks analyses/2026-electrolyte-parameter-gradients/scripts/density_sensitivity.cpp
env -u EPCSAFT_HSL_LIBRARY EPCSAFT_BUILD_JOBS=6 LD_LIBRARY_PATH="$(realpath "$(cat build/environment-wheel/meson-build-dir)")" nice -n 10 agent-heavy build/native-drivers/density_sensitivity/driver analyses/2026-electrolyte-parameter-gradients/results/minimax-density-parameters.txt
```
