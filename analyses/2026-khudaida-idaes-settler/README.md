# Khudaida IDAES electrolyte-LLE settler

This is the positive neutral electrolyte-LLE reference for the IDAES adapter.
It retains the selected ePC-SAFT Engine executions and packet v7 for eight
Khudaida Figure 2 feeds at 293.15 K and 100 kPa. Below, “current” means the
retained packet-v7 result relative to the packet-v6 comparison. The fixed-feed,
live-callback and multiple-record results have separate recorded wheel
identities; none is presented as a new execution on the latest checkout.

```text
1 mol/s explicit feed: H2O + ethanol + isobutanol + Na+ + Cl-
                              |
                              v
                    current ePC-SAFT Engine
                  declared two-liquid equilibrium
                              |
                              v
                  +-------------------------+
                  | fixed-feed IDAES settler|
                  +------------+------------+
                         /             \
                water-rich outlet   organic outlet
```

The Engine owns phase equilibrium, phase amounts, compositions, densities,
and numerical certificates. The square IDAES unit maps those accepted values
to ten outlet component flows through ten equations plus four outlet
temperature/pressure equations. It has zero degrees of freedom and does not
duplicate thermodynamic equations.

## Positive result

All eight current declared two-liquid states and all eight IDAES mappings pass.
Packet v7 restores the Gross–Sadowski 2002 Eq. 3 size factor in the three
unlike association volumes (issue #155). Relative to the v6 states, the phase
compositions move by up to `0.00807` in mole fraction, phase amounts by
`0.00537 mol` and `G/RT` by `0.00430`. These are the historical HELD2 states of
the [lineage diagnosis](../2026-khudaida-lineage-diagnosis/README.md). The
retained rows use wheel SHA-256
`ad12b19fcd8b3b4a87c355490386ecc911485ac1efc4b7a3f547e2393f803e46`.
The [notebook](notebook.qmd) tabulates the largest retained numerical
residuals from the [stage rows](results/fixed-feed/stage-results.csv).

For comparison with Khudaida Table 3, each charge-balanced Engine phase is
converted from explicit ions to one formula amount
`n_NaCl=(x_Na+x_Cl)/2`, then normalized with water, ethanol, and isobutanol.
The 64-coordinate grand formula-basis AAD is:

```text
current Engine (v7):   0.0110549
packet v6, no factor:  0.0128374
paper ePC-SAFT:        0.0161
```

The current result therefore passes the paper-scale thermodynamic criterion.
Against the digitized paper ePC-SAFT organic salt-free branch, its AAD is
`0.0023782` and its maximum absolute difference is `0.0052514` (v6: `0.0069835`
and `0.0119310`). Restoring the size factor moves the Engine toward both the
data and the paper's own model curve.

This is fixed-parameter paper-correlation-scale agreement, not independent
prediction and not agreement within experimental uncertainty. The retained
standard uncertainties are `0.0016` for water, `0.0007` for ethanol, `0.0015`
for isobutanol, and `0.0002` for NaCl mole fraction.

![Current Engine and experimental tie lines on the salt-free ternary diagram](figures/ternary-tie-lines.png)

## Reproduce

Build and install the current non-editable wheel, then run the complete slice:

```bash
bash tools/run-slice.sh khudaida
```

For inspectable candidates:

```bash
.venv/bin/python analyses/2026-khudaida-idaes-settler/scripts/run_fixed_settler.py \
  --output-dir /tmp/khudaida-current
.venv/bin/python analyses/2026-khudaida-idaes-settler/scripts/qualify.py \
  --results-dir /tmp/khudaida-current \
  --output-dir /tmp/khudaida-qualification
```

The solver writes identity before numerical work and atomically checkpoints
each completed point. The qualifier validates 16 unique phase rows, the exact
Figure 2 to Table 3 mapping, the 293.15 K and 5 wt% source slice, all numerical
certificates, and every retained input hash. It exits nonzero if the current
grand AAD exceeds `0.0161`.

The reference is limited to this declared local two-liquid branch. It does not
claim phase discovery, global stability, a live differentiable callback,
optimization, calorics, reactions, or Lithium/HBTA/TOPO qualification.

## Flow basis for a live consumer

The retained outlet component amounts sum to one mole of explicit species,
and the IDAES reference uses one mole per second of explicit feed. A formula
feed of one mole per second with NaCl fraction `z_s` expands to `1 + z_s`
moles per second of explicit species. The live formula-flow consumer must
scale the retained phase-component flows by this factor before comparing
extensive quantities; phase compositions and recovery fractions are invariant.
At point 4 the factor is `1.0276818280876991`. The organic ethanol flow on
the formula-feed basis is `0.0326313940856 mol/s`, giving recovery
`0.966045121390`; dividing the unscaled reference outlet by the formula-basis
inlet would incorrectly give `0.940023551051`.

## Live exact-derivative callback

Issue #832 selects the current #461 rerun in `results/live/`, rather than
the older promoted revision `4dada5de086446d3fb7cf68ad75df7c658c8251d`.
Its identity, baseline and three optimization records, nonlinear programs,
option-file hashes, logs, named multipliers, transport rows, refinement
ladders, endpoint certificates and failure probes are retained. The notebook
checks those file identities and reports this rerun's maximum stationarity
and objective. The older manuscript stationarity bound belongs to the older
revision; it must not be combined with the current rerun. Selection is for
independent delivered review, with investigator promotion still required.
No new model calculation or optimization was run for this curation.

The [live results](results/live/summary.json) use the compiled
equilibrium callback inside an exact-Hessian IDAES/Ipopt NLP, with the same
fixed temperature, pressure and declared two-liquid branch. Three starts
converged to `si = 1.060400841356`, `ss = 0.915238531399`, objective
`1.975639372755`, and ethanol recovery `0.965000000004` (packet v6:
`si = 1.074832400379`, `ss = 0.881523893241`, objective `1.956356293620`). The
maximum independently reconstructed stationarity residual is `3.67e-12`,
complementarity `9.09e-11`, and component closure `1.11e-16 mol/s`; charge
closure is zero. The 2026-10-13 rerun on main (#461) passed all seven #76
falsifiers. Its retained rows and identity use the Engine wheel below.
The [transport rows](results/live/transport.csv) contain 1,050
identical public/ASL values and derivatives. The
[ladder rows](results/live/ladders.csv) retain coarse failures and
refined passing pairs, not only successful points.

![Refining the independent stencil resolves the exact gradient and Hessian-vector actions.](figures/live-derivatives.png)

Each point is the worst of ten outputs and eight directions. The reused neutral
analysis threshold is at least as strict as the unit-floor reporting criterion;
coarse gradient differences fail visibly. Refined pairs resolve second-order
behavior or the documented numerical floor. Lines follow the ordered step
refinement, not experimental observations. No measurement uncertainty is implied.

The 256 [endpoint records](results/live/endpoints.csv) retain branch gap, requested
tolerance, rank and action availability. The NL checks cover the full Jacobian
rows and the feed-link reduction to raw `(si,ss)` row/Lagrangian Hessians; their
maximum error is `1.33e-15`. Named multipliers are retained. Twelve isolated
probes check full F/G/H cache recovery, cold/warm/reordered requests and exact
failure meanings. #461 replaced the two topology fixtures, whose earlier
conditions (`Fi=0.001`; `Fe=1e-8`) solve since solver changes #370 and #174:

- **`boundary`.** `Fe=Fi=0.001` in a cold diagnostic record gives actual phase
  disappearance (`phase_boundary_not_differentiable`). The single liquid there has
  positive TPD against salt-free trial phases and the electrolyte trials tested;
  stability against ionic trial phases is not certified, because the native
  stability tools refuse them for this permittivity model.
- **`unavailable`.** The salt-only feed (`Fw=Fe=Fi=0`) returns "rank-deficient
  declared equations" without a central state, with typed
  `CentralSolveUnavailable` (`central_solve_unavailable`). The cause, two liquids
  of pure NaCl with undetermined amounts, is inferred, not inspected.

These are failures outside the accepted optimization box, not qualified
additional operating states.

The 2026-09-25 packet-v7 replay and its independent review support this
bounded result. On that wheel the installed-wheel `bash tools/check-all.sh` gate
passed 96 Python tests, the native suites, packet validation, Ruff, mypy,
analysis layout (64 leaves), and architecture checks. The 2026-10-13 rerun
(#461) reproduces it on main's wheel with roundoff-level differences, and has not
been independently reviewed. No physical or reactive
Lithium, multistage, global-optimum or physical-validation claim follows from
this example.

The original outer solve crashed because Engine linked Ipopt 3.11.9 while the
outer executable loaded Ipopt 3.14.19; the retained
[loader trace](results/live/diagnosis-crossbindings.log) shows cross-version
C++ symbol binding. An alternate build against 3.14.19 removed the crash but
failed existing Engine checks; that diagnostic wheel has SHA-256
`5ad5c587f0acf39744685a3a9ff2a1f3ec08604828560bbfe07b8fbacae0cb17`.
Its [full Engine gate](results/live/diagnosis-ipopt314-gate.log) fails four
retained Python/native numerical comparisons and one native exact-identity
check. It is not the selected runtime, and its references/tolerances were not
changed. The selected solution uses the installed official IDAES Ipopt 3.13.2
executable, whose static Ipopt symbols are isolated from Engine's system 3.11.9.
All required options are passed through its native option file and verified as
used. An ABI preflight rejects the incompatible shared-library combination.
The 2026-10-13 rerun used non-editable Engine wheel SHA-256
`2751a0dd8f8341940db4d865086d206477435458f3b9c5ee95318e6c33e6072b` (main
`49aac619`, run at commit `ebc84547`) with installed ASL library SHA-256
`56332a7c5d9e1d8ff5b84041acd01335a352d1b8f0aaa7888dd60dd2775d1d43`. All
seven #76 falsifiers passed with the v3 mapped transport record required by
parent #66; record SHA-256 is
`32d32c0f4625dc79249affe7c82b59a5efcaf1542cda8c35099df3e03017de85`. Ipopt uses
HSL, a collection of Fortran codes for large-scale scientific computation;
see the acknowledgement retained in each solver log.

NL row names are retained once in `baseline/stage.row` for all four cases;
optimized column names are retained once in `start-0.9/stage.col` for all
three optimization starts. Their byte-identical copies are omitted.

Reproduce with the default non-editable Engine wheel and installed IDAES tools:

```bash
export PYOMO_CONFIG_DIR="${IDAES_DATA:?Set IDAES_DATA to the external IDAES directory}"
bash tools/build-engine.sh --install-wheel
.venv/bin/python analyses/2026-khudaida-idaes-settler/scripts/run_live_settler.py \
  --output-dir build/idaes-live-check
.venv/bin/python analyses/2026-khudaida-idaes-settler/scripts/plot_live.py
```

Several records in one process (#452): with each wheel installed non-editably,
run every call sequence, then compare the main and changed wheels with the
single-record Khudaida and Lithium runs (the
[notebook](notebook.qmd) holds the result):

```bash
.venv/bin/python analyses/2026-khudaida-idaes-settler/scripts/check_records.py \
  --output-dir analyses/2026-khudaida-idaes-settler/results/live/records/main  # then .../slice
.venv/bin/python analyses/2026-khudaida-idaes-settler/scripts/check_records.py \
  --report analyses/2026-khudaida-idaes-settler/results/live/records \
  --khudaida <main-live-dir> <slice-live-dir> --lithium <main-1> <slice> <main-2>
```

The rendered `results/rendered/notebook.html` is at record commit `21d7a112`; render `notebook.qmd` to view it.
