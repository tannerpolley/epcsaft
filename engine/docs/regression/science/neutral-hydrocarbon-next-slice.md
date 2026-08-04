# Neutral-Hydrocarbon Next-Slice Contract

Status: D-023 defers propane and binary as authority-neutral, non-production
provenance; neither is an active migration blocker. Authority effect none.

This record preserves the smallest regression-owned contract for the two
historical checkpoints. It adds no runtime fit, public export, provider
equation, provider-catalog parameter, promotion, or predictive claim.

## Deferred checkpoint provenance and successor

Checkpoint A was defined to extend the existing pure-saturation engine to
propane.
Validation's direct-experimental source packet is accepted at
`7e51590757f1cb85f51df98e9fe1f88cd4255a88`. The exact candidate is retained
as blocked evidence at Regression
`33228253383ab8755384170b3caac7b38733e157`; the pressure-resolution design
below changes neither that subject nor its runtime.

Checkpoint B was defined to fit one constant methane/ethane `kij`. The manager-
retained provider subject is commit
`45d5764f61729d387100348a8ff91792f6e0a395`, tree
`271d4848faf73afd4cf0683efe5c855053df7d01`, wheel SHA-256
`95d2292b052ab74657931f2dec97c3ea4160d9b17812956515c7195a853e6c5b`,
and installed public-header SHA-256
`6f3a186bf5359f32449a31c544e8b8525be6c594804151d3e74a74e411ded8f4`.
The exact installed artifact passes the derivative and numerical feasibility
checks below, but runtime work is `BLOCKED_FROZEN_PRESSURE_CLOSURE`: the three
declared converged equal-weight starts did not meet the `1e-8` physical
pressure-closure gate. This result does not establish global infeasibility or
model/data incompatibility. Clean Regression HEAD retains no binary executable
target, package API, or runtime implementation.

Migration D-023 supersedes this checkpoint order without rewriting the frozen
evidence. The approved Figiel source and Provider design prerequisites and the
resulting five-ion Ceres contract are now bound only by
`docs/science/figiel-born-diameter-tracer.md`. Migration checkpoint `3a4ef0a`
preserves that complete design as a parallel, deferred Provider-dependent
track behind HELD2; runtime remains unauthorized. This historical neutral-
hydrocarbon record admits no Figiel residual, callback, fit, persistence
contract, or implementation surface.

## Frozen source and model identity

The binary target uses all 17 rows in validation's audited May 2015 packet:

- validation source commit:
  `73a37f5935e919a34d1e4fa3af285951d6fac8e7`;
- source CSV SHA-256:
  `5cd1e74925a3c6504f5106dcf911f2cae2d6e99a5133fccc20454d8991bdbc7f`;
- source packet SHA-256:
  `d43433e93b354e01f96d330c760818a24b775026461ce795e45774cfb11ac94e`;
- component order: methane, ethane; and
- expected base-provider fingerprint:
  `sha256:307fcb28d535b94782f3e3caf4012c0c8c0dc87ee4239d6c316de56553543286`.

The clean Gross--Sadowski pair record fixes the source value `kij = 0` and
identifies the mixing rule used by the provider as

```text
epsilon_ij = sqrt(epsilon_i * epsilon_j) * (1 - kij).
```

The pair-record file SHA-256 is
`747e8281c7a1e4240ee4badbc0bedd047521fb303726699b38c99fccf7f74c2a`.
The regression start is therefore exactly zero. The approved engineering
search interval is `[-0.15, 0.10]`; it is a bounded neighborhood of the source
value and keeps the cross-dispersion energy positive. It is not attributed to
the paper. The transformed parameter is

```text
kij = 0 + 0.01 * z_kij.
```

The fit is rejected if `kij` is nonfinite or within `1e-8` of either bound.
Bounds cannot be widened after observing the fit.

## Variables and scaling

For row `r`, measured `T_r`, `P_r`, liquid methane fraction `x_r`, and vapor
methane fraction `y_r` remain fixed. Each phase contains one mole:

```text
n_L = (x_r, 1 - x_r) mol
n_V = (y_r, 1 - y_r) mol
```

Each row lifts two positive volumes:

```text
V_L = 6.5e-5 m3 * exp(u_L)
V_V = (R * T_r / P_r) * exp(u_V)
```

The liquid scale lies inside the `5.54e-5` to `8.56e-5 m3` range retained by
the current clean, authority-neutral May campaign for its 16 evaluated local
states. It is initialization evidence only. The vapor scale is the
dimensionally exact ideal-gas one-mole volume. Frozen bounds are
`[2e-5, 1e-4] m3` for liquid and `[1e-4, 1e-2] m3` for vapor. They contain the
current source-state initialization evidence and enforce phase ordering
without density roots.

The complete variable order is one `z_kij`, followed by `(u_L, u_V)` in the
17 immutable source-row order. The problem therefore has 35 dimensionless
variables.

## Residuals

The provider owns

```text
Phi = A / (R T n_ref),  n_ref = 1 mol
P = -R T n_ref * Phi_V
mu_i / (R T) = Phi_ni
```

For each row, regression assembles four residuals with common numerical weight
`w = 0.25`:

```text
r_LP  = w * (P_L - P_r) / P_r
r_VP  = w * (P_V - P_r) / P_r
r_mu1 = w * (Phi_n1_L - Phi_n1_V)
r_mu2 = w * (Phi_n2_L - Phi_n2_V)
```

All 68 residuals are dimensionless. The two pressure residuals use the one
observed coexistence pressure to define each lifted phase state; they are not
treated as two independent measurements. The weights make no uncertainty
claim. Measured compositions are inputs, not lifted data-reconciliation
variables.

## Exact Jacobian contract

The provider callback must return `Phi`, its gradient, and its row-major
symmetric Hessian in exact coordinate order

```text
(n_methane_mol, n_ethane_mol, V_m3, kij).
```

For phase `a` and component `i`, the consumed derivatives are

```text
d(r_aP)/d(u_a)   = w/P_r * (-R*T*n_ref*Phi_VV_a) * V_a
d(r_aP)/d(z_kij) = w/P_r * (-R*T*n_ref*Phi_Vk_a) * 0.01

d(r_mui)/d(u_L)   =  w * Phi_niV_L * V_L
d(r_mui)/d(u_V)   = -w * Phi_niV_V * V_V
d(r_mui)/d(z_kij) =  w * (Phi_nik_L - Phi_nik_V) * 0.01
```

All other row-to-volume entries are zero. This produces the exact `68 x 35`
Ceres residual Jacobian from a second-order CppAD provider tensor. No third
derivative, density root, implicit branch derivative, equilibrium dependency,
or numerical production derivative is consumed.

## Future result and acceptance diagnostics

These are the documented requirements for the later binary implementation;
no forward result or diagnostic types are packaged by the current preparation.

The private immutable result contract retains exact source/provider hashes;
separate solver, numerical, and physical statuses; termination and usability;
initial/final costs; iteration count; start, final, movement, bounds, and
active-bound state for `kij`; every source row and its raw/scaled residuals;
phase volumes, pressures, chemical potentials, stability slopes, and rejection
reasons; full and projected Jacobian diagnostics; confirmation outcomes; and
the predictive status
`NOT_ADJUDICATED_NO_APPROVED_HELD_OUT_CUTOFF`.

Solver acceptance requires finite complete Jacobian columns, Ceres
`CONVERGENCE`, a usable finite solution, cost reduction from the source start,
full scaled-Jacobian rank 35, projected-parameter rank 1, finite nonzero
parameter movement, and a finite non-bound fitted `kij`. Rank uses

```text
sigma_max * max(68, 35) * epsilon_double.
```

The projected parameter column is

```text
j_k_perp = (I - J_V * pinv(J_V)) * j_k.
```

Its one nonzero singular value and rank are retained; no covariance,
uncertainty, or global-identifiability interpretation is allowed.

Numerical confirmation perturbs all liquid starts by `1.01`, all vapor starts
by `0.98`, and repeats from `kij = -0.05` and `+0.05`. It requires maximum
scaled-`kij` agreement `<= 1e-5` and symmetric relative-cost agreement
`<= 1e-8`.

Physical status requires exact source and provider identity, positive amounts
and volumes, provider-domain success, ordered separated phases, stable phase
slopes, and scaled pressure closure `<= 1e-8`. The fitted chemical-potential
residuals remain reported calibration errors because no accuracy cutoff was
approved.

## Installed-artifact preflight result

The preflight consumed the retained installed wheel/header and audited CSV,
with no sibling provider source. Provider derivative checks, the exact
`68 x 35` Jacobian check, rank 35, projected parameter rank 1, a finite
non-bound fitted `kij`, and confirmation from starts `-0.05` and `+0.05`
passed. The three declared converged equal-weight starts did not meet the
`1e-8` pressure-closure gate. The primary start reached
`0.022691483584976503` at row `may2015-ch4-c2h6-002`; no weight, tolerance,
bounds, variable, row, provider binding, or formulation was changed.

The compact retained receipt is `evidence/binary-kij-preflight.json`. The
rejected executable alternative and its original receipt subject remain
immutable provenance at Regression
`47a2a8d9579a01347df7ccaa977337ad7d4047af`, subject
`d51c9f0713b6a7355be719b6843e4459f41d46d16d973668d694715d36b63676`,
and Migration `4fd0e6eff5fd17399573132430d13c7e34626b98`. Clean HEAD retains no
replay harness. The corrected retained receipt subject and file SHA-256 are
`f9033e6d1b00b0bfc11b421e6ef1e388db714d95431618b981d294d11f728d62`
and `c734bb9f56768f89a710575d58c45b460543c0828e913c5ea5372603c768ff78`.

## Cross-checkpoint pressure-resolution evidence protocol

The canonical formula, candidate-derivation comparison, pure-component replay
matrix, and status semantics are in
`docs/science/pure-saturation-regression.md`. This record adds only the binary
and staged-maturity constraints. No numeric `atol_resolution` is selected.
Under D-023 this entire protocol is retained provenance, not an active
Regression, Provider, or Validation prerequisite.

The frozen candidate formula is

```text
abs(P_phase - P_report)
    <= atol_resolution + 1e-8 * abs(P_observed).
```

The `1e-8` relative term, residual weights, rows, variables, bounds, starts,
provider callbacks, and solver controls remain unchanged. Observed pressure is
a magnitude scale, not experimental uncertainty. The pure protocol derives one
raw-Pa candidate from its predeclared calibration states using the certified
high-precision discrepancy plus input-lattice bound. Propane 120 K and every
binary row are locked challenges: neither blocked closure value can derive or
widen the candidate.

### Frozen artifacts and commands

This section is immutable historical provenance, not an executable guide for
the current `main` runner. Replaying it requires the exact Regression subject
and its contemporaneous `tools/run_candidate.py`; the active Provider `0.2`
runner must not be combined with these Provider `0.1` arguments or used to
rewrite the retained receipts.

The immutable pure subject is:

| Item | Identity |
|---|---|
| Regression blocked evidence | commit `33228253383ab8755384170b3caac7b38733e157`, tree `23f9528701123165a90b4887dbebf7e245c2b053` |
| Regression implementation | commit `aab87ebd4a40cb29f21486e06687c10eb1e44624`, tree `462cda550cf5461e4b13cc7e40630708dd336bd7` |
| Regression wheel | SHA-256 `32b815fc00241516f13574594af3cf631f8fe30629bbec45ac3828d9357d705f` |
| Propane receipt | subject `c059e5381a3fbf0cf00a43aaa0cc28b67c074e325f1a33066b6937c16cc0a761`, file SHA-256 `daaa93a2f5d0e0dbf81fc83562bb809c687a36e66578071c68ea66b1988ccaec` |
| Pure provider | commit `4b10cb899c94687cae734980285badb224dc95e6`, wheel SHA-256 `f92f79c8d6f614660e5c201b7061c9b02b5cd1a25a4ed8c8fee0b59adaabf2bf`, installed header SHA-256 `414c257d28322d6be41809a8dc6b98023859dd156202a5205079d674b28b4070` |

The immutable binary subject is:

| Item | Identity |
|---|---|
| Rejected executable provenance | Regression commit `47a2a8d9579a01347df7ccaa977337ad7d4047af`, subject `d51c9f0713b6a7355be719b6843e4459f41d46d16d973668d694715d36b63676` |
| Compact blocked evidence | subject `f9033e6d1b00b0bfc11b421e6ef1e388db714d95431618b981d294d11f728d62`, file SHA-256 `c734bb9f56768f89a710575d58c45b460543c0828e913c5ea5372603c768ff78` |
| Active-`kij` provider | commit `45d5764f61729d387100348a8ff91792f6e0a395`, tree `271d4848faf73afd4cf0683efe5c855053df7d01`, wheel SHA-256 `95d2292b052ab74657931f2dec97c3ea4160d9b17812956515c7195a853e6c5b`, installed header SHA-256 `6f3a186bf5359f32449a31c544e8b8525be6c594804151d3e74a74e411ded8f4` |
| May source | Validation commit `73a37f5935e919a34d1e4fa3af285951d6fac8e7`, CSV SHA-256 `5cd1e74925a3c6504f5106dcf911f2cae2d6e99a5133fccc20454d8991bdbc7f` |

Before any campaign, these frozen inputs are checked with:

```bash
git rev-parse HEAD HEAD^{tree}
sha256sum \
  ../artifacts/regression-pure-saturation-v1/aab87eb/epcsaft_regression-0.1.0.dev0-cp313-cp313-linux_x86_64.whl \
  ../artifacts/provider-native-sdk-v1/4b10cb8/epcsaft-0.1.0.dev0-cp313-cp313-linux_x86_64.whl \
  evidence/propane-candidate-fit-receipt.json \
  evidence/binary-kij-preflight.json \
  evidence/may-2015-methane-ethane-vle.csv
```

The frozen pure artifact replay command is:

```bash
python3.13 tools/run_candidate.py \
  --component propane \
  --provider-wheel ../artifacts/provider-native-sdk-v1/4b10cb8/epcsaft-0.1.0.dev0-cp313-cp313-linux_x86_64.whl \
  --provider-test-receipt ../artifacts/provider-native-sdk-v1/4b10cb8/provider-tests.xml \
  --regression-wheel ../artifacts/regression-pure-saturation-v1/aab87eb/epcsaft_regression-0.1.0.dev0-cp313-cp313-linux_x86_64.whl \
  --output /tmp/propane-resolution-replay.json
```

Methane and ethane use the same command with only `--component` changed. The
output is temporary comparison evidence and must not replace the retained
receipt. The clean Regression tree does not regain the historical binary
replay harness. Binary resolution evidence must be produced in a temporary
detached worktree at exact commit `47a2a8d9579a01347df7ccaa977337ad7d4047af`,
consume the retained active-`kij` wheel above, write only outside the worktree,
and be deleted after its hashes and result are transferred to the independent
review receipt.

Provider's smallest prerequisite is two non-installed test-only oracle cases
in its existing `tests/test_native_sdk.py` owner: one finite pure-parameter
pressure matrix and one active-`kij` pressure matrix, each evaluated with an
independent directed-rounding high-precision transcription. Their review
command is frozen as:

```bash
python3.13 -m pytest -q \
  tests/test_native_sdk.py::test_native_sdk_pure_parameterized_pressure_resolution_oracle \
  tests/test_native_sdk.py::test_native_sdk_active_kij_pressure_resolution_oracle
```

Those test names do not exist in the retained provider subjects; their absence
is an explicit prerequisite, not permission to fabricate output. Provider
must return an exact commit, tree, wheel/header identities, oracle-data hash,
and command result before Regression can propose a number.

Validation then extends its existing installed-artifact pure-saturation
campaign rather than creating a generic tolerance framework. Its source and
campaign checks remain:

```bash
python3.13 -m pytest -q tests/test_pure_saturation_regression.py
python3.13 campaigns/pure_saturation_regression.py \
  --provider-wheel ../artifacts/provider-native-sdk-v1/4b10cb8/epcsaft-0.1.0.dev0-cp313-cp313-linux_x86_64.whl \
  --regression-wheel ../artifacts/regression-pure-saturation-v1/aab87eb/epcsaft_regression-0.1.0.dev0-cp313-cp313-linux_x86_64.whl \
  --output-dir /tmp/propane-pressure-resolution-campaign
```

The current Validation campaign does not yet admit the propane artifact or the
resolution oracle. The command freezes its exact input subject and output
location; it becomes executable evidence only after the bounded Validation
extension is reviewed. A later corrected Regression artifact requires a new
hash-bound command and cannot silently replace this blocked input.

### Binary replay matrix and gates

The binary matrix includes all 17 May rows for each of the three frozen starts
`kij = 0`, `-0.05`, and `+0.05`. For every converged solution it records the
center and the one-axis `nextafter` minus/plus probes for liquid volume, vapor
volume, and the observed common pressure used in residual subtraction. It
retains callback/oracle pressure, `dP/dV`, both component chemical potentials,
phase stability, topology, all four raw/scaled residuals, ranks, fitted `kij`,
cost, and confirmation deltas. The common pressure remains fixed source data;
its probes measure binary64 subtraction only and never reconcile data.

The one candidate `atol_resolution` may be applied to binary evidence only if
every active-`kij` oracle bound lies within the pure-derived candidate. A
failure leaves Checkpoint B blocked and cannot widen the number. The row-002
closure is then evaluated, never used to enlarge the candidate. Given the
retained `2.2691483584976503e-2` relative miss, Checkpoint B remains `NOT_READY`
unless the independently derived raw-Pa term passes the mixed criterion
without violating negative controls. No current evidence makes that claim.
The retained binary contract
also has no approved chemical-potential equality cutoff; its maximum absolute
component residual is `0.06950599355260767`. This pressure-resolution protocol
does not create that missing physical-admission criterion, so a mixed pressure
pass alone cannot make Checkpoint B ready.

Accepted methane/ethane decisions and raw numerical outputs must remain
unchanged. Ethane 100 K and propane 110 K remain failed. Both signs of a
pressure perturbation

```text
max(100 * atol_resolution, 1e-6 * abs(P_observed))
```

must reject for pure and binary representative rows. Solver convergence,
derivative, rank, bounds, topology, stability, chemical-potential, source, and
provider gates are unchanged. Predictive status remains
`NOT_ADJUDICATED_NO_APPROVED_HELD_OUT_CUTOFF`.

## Staged parameter-family readiness

`READY` requires all six columns below. A missing column makes the family
`NOT_READY`; similarity to an admitted family is not evidence.

| Parameter family | Source-backed target | Exact provider seam | Rank-sufficient formulation | Bounds/scaling/confirmation | Installed-artifact plan | Approved scope | Readiness |
|---|---|---|---|---|---|---|---|
| Pure `m`, `sigma`, `epsilon/k` | yes for methane, ethane, propane | yes, five-coordinate Hessian | yes, rank 3 and full lifted system | yes | retained propane campaign only | yes, component-specific | methane/ethane accepted; propane `NOT_READY`, deferred non-production under D-023 |
| Constant methane/ethane `kij` | yes, all 17 May rows | yes, active-`kij` value/gradient/Hessian | yes, rank 35 and projected rank 1 | yes | frozen preflight only; no runtime artifact | historical D-022 evidence only | `NOT_READY`, deferred non-production under D-023 |
| Association parameters | no admitted target | no admitted complete parameter/association-state derivative seam | no | no | no | no | `NOT_READY` |
| Polar parameters | no admitted target | no admitted polar parameter derivative seam | no | no | no | no | `NOT_READY` |
| Current-catalog ion Born diameter | approved Validation packet `8944d34f` | approved Provider design `da966048`; no runtime | frozen `5 x 5` diagonal design | frozen in sole design owner | frozen; no artifact exists | D-023 parallel/deferred design only | design complete; runtime `NOT_READY` pending a Migration-bound Provider artifact |
| Reactive or temperature-dependent parameters | no admitted target | no admitted exact derivative seam | no | no | no | no | `NOT_READY` |

Uncertainty, global identifiability, provider-catalog persistence, and a
generic target/parameter registry are not parameter families admitted by this
matrix. They remain excluded rather than predesigned.

## Explicit exclusions

This preparation does not admit predicted `x/y`, held-out binary evidence,
uncertainty, covariance, temperature-dependent interactions, other binary
pairs, association, electrolytes, equilibrium, density closure, global
identifiability, parameter persistence, a generic target registry, promotion,
publication, or authority transfer.
