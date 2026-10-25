# Density closure diagnosis — 2026-09-21

This record preserves the bounded #89 baseline diagnosis and its subsequent
repair. Sections before **Repair outcome** describe the original wheel and
historical design/review gates. It is a root-search and declared-branch result,
not phase discovery, global-stability evidence or application qualification.

## Baseline and owners

- Engine source: `f909bb21a93171c38465adba8777cd3b0bc6c53d` on
  `roadmap/application-preparation`.
- Qualified non-editable wheel:
  `build/environment-wheel/epcsaft-0.2.0.dev0-cp313-cp313-linux_x86_64.whl`.
- Wheel SHA-256:
  `9f73aeb466c54eff80c89a8776c90d49cb9c9c2d130f0eb15bf0311c1f76e606`.
- Parameter packet: `data/packets/gross-2001-methane-ethane/2/parameters/parameters.json`,
  SHA-256 `7dd1899712aacfbc3395d86349768b0d93dd5dfd75b2e0f0700640be618aa9c2`.
- Model: pure ethane selected from the packet, composition `[1.0]`.
- Retained rows: [`results/root-controls.csv`](results/root-controls.csv).

The source equation is `F(rho) = P_EOS(T,rho,x) - P* = 0` in SI units. The
admitted density domain is the current model close-packing domain, with the
runtime scan using `rho_lower = 1e-8 P*/(R T)` and
`rho_upper = 0.999 rho_max(T,x)`. A root is locally mechanically stable only
when the native analytic slope `P_rho = (dP/drho)_(T,x)` is finite and positive.
These are the `density_residual`, `density_upper_bound`, `ideal_gas_density` and
`mechanical_stability` owners in `engine/docs/equations.md`.

The current implementation is one owner in
`engine/native/properties/density.cpp`: a 40-interval logarithmic scan brackets
sign changes, `refine` applies safeguarded Newton steps, positive-slope roots
are sorted, and one stable root is returned with branch `Unspecified`. Multiple
stable roots require an explicit vapor/liquid request; an unspecified request is
ambiguous. The bounded scan count is not a completeness certificate.

Meaningful callers are all reached through this owner:

- `HelmholtzModel::state_at_pressure` preserves the selected density,
  convergence record, branch and stable-root count.
- `HelmholtzModel::pressure_actions` resolves the root before pressure-chart
  actions and carries a typed unavailable cause when closure throws.
- `equilibrium::boundary_start` uses declared vapor/liquid roots to seed the
  single equilibrium driver.
- `equilibrium::check_density_branches` uses the same branch admission for public
  equilibria and anchored continuation;
  `gibbs_minimizer::copy_previous_phase` re-solves pinned phases at a prior
  pressure state.
- Python `Pure.state`, `Mixture.state` and pressure actions map directly to the
  same native path.

## Independent reference and frozen controls

The reference used the installed wheel's fixed-density `state(T, rho)` path,
which does not call `state_at_pressure` or `solve_density`. For each control it:

1. found the finite-state close-packing boundary by bisection and used
   `1e-3 <= rho <= 0.999999 rho_cp` for the independent reference;
2. sampled `P_EOS(T,rho)-P*` at 16,384 logarithmic intervals;
3. refined every sign-changing bracket by bisection until the log-density
   bracket was below `1e-14`; and
4. classified each refined root independently with the returned analytic
   `P_rho`.

This is a finite reference ladder, so it supports the recorded controls and
does not prove global root completeness. The shifted 256-cell falsifier uses
the production scan endpoint `0.999 rho_cp`; the CSV records both that native
endpoint and the separate `0.999999 rho_cp` reference endpoint. It also records
the exact roots, slopes, pressure residuals and the current pressure-closure
result for each control.

The controls were frozen before probing:

| Control | Inputs | Independent expectation | Purpose |
|---|---|---|---|
| Regular one-root | `T=200 K`, `P*=1.0e8 Pa`, `x=[1]` | One positive-slope root near `20042.858 mol/m^3`; a vapor or liquid request must return it unlabelled with count 1 | Legitimate sole-root semantics |
| Regular two-root | `T=200 K`, `P*=217800.69494692178 Pa`, `x=[1]` | Stable roots near `138.1005` and `17541.2513 mol/m^3`, separated by one negative-slope root; explicit branches select the low/high roots and unspecified is ambiguous | Branch identity and ambiguity |
| Near coalescence | `T=308.8 K`, `P*=5147500 Pa`, `x=[1]` | Stable roots near `5930.3996` and `6828.2950 mol/m^3`, separated by one negative-slope root near `6422.1983` | Falsify a false one-root result when crossings fall in one coarse bracket |
| Out-of-domain/no-root | `T=200 K`, `P*=1.0e10 Pa`, `x=[1]` | No finite positive-slope crossing in the admitted density domain; preserve typed `density_unavailable` | Distinguish unsupported state from a missed regular root |

All temperatures are K, pressure is Pa, density is mol m^-3, and `P_rho` is
Pa m^3 mol^-1. Numerical convergence is separate from the physical/root
reference comparison.

## Results and demonstrated cause

The regular controls support the existing P2/#86 semantics. At `200 K` and
`1e8 Pa`, all phase requests returned the same sole stable root with branch
`Unspecified` and count `1`; treating that root as vapor or liquid would be an
overclaim. At `200 K` and `217800.69494692178 Pa`, the explicit vapor and liquid
requests returned the independent low and high stable roots with count `2`,
while an unspecified request raised `density root is ambiguous (stable count 2)`.
The `1e10 Pa` control raised `no stable density root (count 0)` and the finite
reference ladder found no admitted crossing.

The near-coalescence control demonstrates a current density-owner defect. The
independent reference has three crossings:

- stable vapor-side root: `rho=5930.3996147034695`, `P_rho=8.120790200179608`;
- unstable middle root: `rho=6422.198255005569`, `P_rho=-3.726631786036806`;
- stable liquid-side root: `rho=6828.295025027459`, `P_rho=6.911968641531208`.

The current 40-point scan has only one sign-changing cell in this region,
`rho=[4210.852286093733, 7170.419984431981]`, with reduced residual signs
corresponding to dimensional residuals `[-180141.2103, 6114.1871] Pa`. That
single cell contains all three crossings. Its safeguarded refinement returns
the low stable root, then the positive-slope filter reports one stable root and
the branch is downgraded to `Unspecified`. Both explicit vapor and liquid
requests therefore return the same low root instead of either selecting a
declared branch or preserving ambiguity. This is a missed-root/aliasing result,
not legitimate one-root physics and not proof of phase nonexistence.

The cause is bounded bracket resolution followed by a one-root-per-sign-change
refinement. `refine` is valid for a bracket containing one crossing; it cannot
enumerate multiple crossings hidden inside one coarse cell. The positive-slope
filter correctly rejects the middle root once found, but it cannot recover a
stable root that was never bracketed.

## Corrected diagnosis and bounded aliasing probe

The 40-cell result is a demonstrated missed-root/aliasing defect in the
existing density owner. A constant 256-cell replacement is **not** ready for
Build. The independent grid ladder on the frozen near-coalescence state found
one sign-changing interval at 40, 64, 80, 96, 128 and 160 cells, three at 192,
one at 224, and three at 256. The retained `scan-grids.csv` records the native
admitted bounds and computed 256-cell logarithmic step; that step is larger
than both adjacent root gaps. Therefore the unshifted 256 result is an
alignment observation, not a reliable correction.

Independent integrated review supports the recovered bounded artifact. The
recovery removed the duplicated fixed-density branch field and copied grid-step
scalar; the pressure-closure branch and exact grid step remain with their
respective retained tables. All other numerical columns were preserved. This
accepts the diagnosis, not a replacement density method.

Before the added probes, four alignments were frozen: shifts `0`, `0.25`,
`0.5` and `0.75` of each 256-cell interior, always retaining the same two
endpoints. They are applied only to the four existing controls and use the
native `0.999 rho_cp` endpoint. The retained `scan-grids.csv` shows the
near-coalescence bracket counts `3, 1, 1, 3` for those shifts. The regular
one-root and two-root controls remain stable across all four alignments, and
the out-of-domain control remains root-free. This is structural grid aliasing,
so the current count cannot be treated as a branch completeness result.

The exact replay is:

```text
.venv/bin/python analyses/2026-density-closure-diagnosis/scripts/replay.py --output-dir analyses/2026-density-closure-diagnosis/results
```

The script verifies the qualified wheel SHA-256, parameter-packet SHA-256 and
the installed wheel's `direct_url.json` before running. It regenerates the
original eight control rows and writes one table per result type:

- [`root-controls.csv`](results/root-controls.csv) retains the one-root,
  two-root, near-coalescence and out-of-domain controls.
- [`scan-grids.csv`](results/scan-grids.csv) records endpoints, residuals,
  reference crossings, the native and reference upper bounds, all four shifts
  and bracket details.
- [`pressure-actions.csv`](results/pressure-actions.csv) compares
  `d rho/dP = 1/P_rho` with an independent same-branch pressure ladder for
  each of the five stable reference roots. The ladder uses a relative
  `1e-6` pressure step and the exact #86 criteria: fine error `5e-6` and
  coarse error `4 * 5e-6`, both scaled by `max(1, |reference|)`. Native action
  identity uses the #86 `atol=0`, `rtol=1e-12` reciprocal criterion evaluated
  at the density returned by that same native solve. The `expected_*`,
  `coarse_*`, `fine_*` and `native_*` action fields are in mol m^-3 Pa^-1;
  pressure slopes are in Pa m^3 mol^-1.
- [`anchors.csv`](results/anchors.csv) retains positive vapor/liquid P2 anchor
  replays plus finite positive out-of-domain and negative-anchor controls.
- [`replay-summary.json`](results/replay-summary.json) records provenance,
  frozen probe limits and claim limits without adding another numerical table.

The pressure-action ladder is informative rather than application
qualification. The independent smooth-branch ladder meets the #86 coarse and
fine limits for all five stable roots. Native action identity is evaluated
against `1/P_rho` at the density actually returned by the native action call,
not against the reciprocal slope at an independently located reference root.
It passes for every emitted native action, including the near-coalescence low
root: native `rho=5930.39961467628`, `P_rho=8.120790201562388`, and relative
identity error `3.49e-14`. The desired near-coalescence high branch is still
recorded as `missed_by_current_cold_scan`; its independent branch ladder is
retained. These are numerical closure observations, not physical error
limits.

## Baseline design gate and ownership (historical)

The remaining design question belongs to `engine/native/properties/density.cpp`
and, if a policy constant is needed, `engine/native/core/numerical_policy.hpp`:

> What finite derivative-aware bracket or aliasing test can identify a cell
> that may contain unresolved crossings, preserve explicit branch semantics,
> and fail closed within a stated evaluation and wall-time budget?

The design must retain the current safeguarded refinement, positive finite
`P_rho` admission, one-root `Unspecified` semantics, explicit multi-root
selection and typed zero-root failure. Its falsifiers must include the four
frozen controls and the four alignment results, and must distinguish a missed
root from an honestly ambiguous branch. A further constant scan-size increase,
phase discovery, global stability search or invented robust solver is not a
diagnosis result and is not proposed here. The smallest implementation-ready
batch is therefore **not yet established**. A fresh independent CSE review is
required before any kernel Build; this analysis made no kernel, policy, wheel,
application or GitHub changes.

## Diagnosis correction history

- Initial review identified the need to separate the finite reference endpoint
  from the production scan endpoint and to compare native actions at their
  returned state.
- Correction round 1 retained the four frozen shifts, the bounded action
  ladder and the anchor controls.
- Correction round 2 uses the production `0.999 rho_cp` endpoint for the
  shifted runtime-grid falsifier, renames action-unit fields, and removes the
  prior false low-root identity failure. The candidate remains awaiting its
  independent readback; no final-review claim is made here.

## Diagnosis limits and review state (historical)

This handoff establishes one demonstrated root-search gap, its current owner,
the finite aliasing falsifier and the derivative/action and anchor controls. It
does not qualify MEA, absorber or Lithium physical predictions, establish phase
discovery, certify global stability, or authorize kernel Build. The qualified
wheel and packet identities, bounded run budget and current numerical failures
are retained in the replay summary. This correction is the second automatic
correction round for the same acceptance boundary; work stops before a third
round if the independent readback finds another material defect.

## Repair outcome — 2026-09-21

The recovery keeps pressure closure in `engine/native/properties/density.cpp`
(ALG-DENSITY). The supported isotherm has one mechanical loop. After the scan
exposes its negative-slope region, pressure-slope sign transitions divide the
cells into mechanically stable and unstable intervals. A stable interval is
monotone in pressure and has at most one root. Its strict pressure crossing
is refined once; an exactly zero stable scan node is owned once while its
neighboring cells remain connected.

This replaces the competing endpoint rules in the earlier PR #93 candidates:
double-ended Newton refinement, coincidence-radius merging, recursive
pressure-extremum splitting and the forced first step. Those rules first
counted one near-tangent vapor root twice, then dropped an exact liquid grid
root when the forced step was introduced. The exact-grid defect was verified
on head `46ffb2ca`, wheel `38c30256…`: at 308 K and
5057915.461281232 Pa all requests returned 5124.770479 mol/m³ with count one.
Independent fixed-density bisection also found the stable liquid root
7160.554343 mol/m³ and an unstable root at 6950.006505 mol/m³. One pressure ULP
on either side returned both stable roots in that candidate.

The persistent numerical checks retain the original nine loop states and add
that exact-grid pressure and its two immediate floating-point neighbors,
plus a near-tangent state where residual-only stopping exceeded the returned
local-sensitivity radius (0.0831 mol/m³ error against a 0.0776 mol/m³ radius).
Refinement now requires both the reduced-residual criterion and a remaining
bracket no wider than that local radius, unless it hits an exact numerical
zero. A residual already inside tolerance triggers bisection to contract the
remaining bracket; exhausting the existing iteration budget is a failure.
They independently bracket and bisect fixed-density pressure crossings,
classify their analytic pressure slopes, compare the two declared branch
densities and require an unspecified request to be ambiguous. The existing
regular sole-root, no-root, anchor and pressure-action cases remain in P2/P4.
The installed Python check transports both the original 308.8 K loop and the
308 K grid-coincidence case. The latter failed on the pre-recovery wheel
before the replacement was built.

A separate P2 check locates the vapor-side pressure maximum at 308 K by
independent bisection of the analytic fixed-density slope. A target within
the existing reduced pressure-residual tolerance of a fold is explicitly
unresolved. It is not returned as a distinct stable branch or a unique root.
This preserves the distinction between numerical resolution and physical
existence. Bisection of a fold has the existing 60-step budget and ends at
adjacent representable densities or an exact slope zero.

The golden-section search for a hidden negative slope remains a bounded
search: only a genuine sampled local-minimum bracket is used for its
quadratic interpolation. A positive lower bound for that interpolating
quadratic is a model estimate, not a bound on unsampled EOS curvature.
A dip narrower than the sample spacing with greater curvature is not
excluded. In particular, the critical cubic pressure expansion is in density,
so it does not make a quadratic slope model in log density exact.
No global root-completeness or phase-stability claim follows.

Historical sweep qualification on wheel `38c30256…` comprised 228 states:
171 near-critical states from 57 temperatures, 48 near-tangent states and the
nine retained cases. Three proposed temperatures were skipped because the
loop pressure span was below 0.001 Pa. Its success domain reached approximately
`Tc - 1.42e-5 K` and tangency distances of at least 0.001 of the loop span.
The comparison radius `tol |d rho/d r|` is first order; it was exceeded by
nonlinear error in additional, closer-to-tangency probes. That history does
not establish the accuracy or cost of a later wheel. Current numerical
qualification and the installed-wheel identity are reported in
[PR #93](https://github.com/tannerpolley/ePC-SAFT/pull/93); the native comparison
rows retain the actual values, slopes, criteria and root counts.

The anchor uses the same local partition and counting method within its
existing window. Local uniqueness does not establish uniqueness outside that
window or continuation through a disappearing branch. More than one
mechanical loop and exact critical-point solving remain outside this
qualification. Application physical predictions remain with #61/#63.
