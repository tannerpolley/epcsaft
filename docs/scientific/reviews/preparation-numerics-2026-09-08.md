# Numerical preparation review — 2026-09-08

Status: preparation evidence only. This record retains two runnable numerical
falsifiers before production corrections; it does not qualify a scientific
sublayer. The branch snapshot is `greenfield@4239e49b3a0f7837518148b0bd41181c35e0e8cf`.

## Question and scope

Can the current Engine expose the two demonstrated numerical defects as
independent, reproducible checks while keeping the reference, units, tolerance
rationale and evidence age explicit? The scope is the residual caloric identity
(F1) and the Ipopt weighted Hessian callback (F5). No producer correction,
roadmap edit or algorithm choice is made here.

## Evidence state

The historical caloric probe used the Gross–Sadowski methane/ethane packet,
revision 2, at `T=300 K`, `rho=10000 mol/m^3`, pure ethane. Its centered
temperature rows are retained once in
[`transfer-audit-2026-09-06-caloric.csv`](transfer-audit-2026-09-06-caloric.csv).
The audit recorded a large mismatch between the returned residual internal
energy and both the enthalpy/compressibility relation and the temperature
derivative identity. This is historical wheel evidence, not a current
qualification result.

The current source still contains the two inspected F1/F5 failure paths:

- `engine/native/model/helmholtz_model.cpp` computes the residual internal
  energy without the second temperature factor required by
  `u^r = -R T^2 (partial a^r / partial T)`.
- `engine/native/equilibrium/request_compile.cpp` returns a zero Hessian for a
  pinned problem, while `engine/native/equilibrium/gibbs_minimizer.cpp` scales
  the complete Hessian by the objective factor. The callback relation is
  `nabla^2 L = alpha nabla^2 f + sum(lambda_j nabla^2 g_j)`; callback
  multipliers are separate from chemical balance multipliers.

The baseline owner verified the current noneditable wheel before this focused
native run (`40b5507982ae8b15067afd860a9379e166a9491e641783b8fd3b1d04de1bdd67`).
The focused native build used
`build/engine-meson/cpython-313-linux-x86_64-x86_64-ad7b04df3dc8cbafbb10255b1cf307664e9b8e7b6d4114458df39608fb7901be`;
the executable hashes were `model_checks`
`aa8ee1d37b6295c95528d04ed16181e750f37b1c1c6cb7791425981de6f342e1` and
`equilibrium_checks`
`cfab2b5daa9fed9650bc0dae7643ec35cbe87f245a3f5994d4dd06f06d6bf87d`.
The current focused outputs are retained in
[`preparation-numerics-2026-09-08.csv`](preparation-numerics-2026-09-08.csv).
This slice does not claim a full gate or application replay.

## Retained falsifiers

| ID | Runnable check and input | Independent expectation | Local comparison | Status at preparation |
|---|---|---|---|---|
| F1 | `model_checks F1` function `residual_internal_energy_identity_ladder`; Gross revision-2 ethane, `T=300 K`, `rho=10000 mol/m^3`, `x=(1)` | Fixed-`rho,x` `u^r=-R T^2 a^r_T`; independently `u^r=h^r-RT(Z-1)`; centered `a^r(T+Delta T)-a^r(T-Delta T)` ladder | `Delta T=(0.1,0.03,0.01,0.003) K`; provisional `2e-7` relative ladder allowance is above the largest retained historical refinement error with a small margin; provisional `1e-10` direct-identity allowance is tied to same-state floating-point evaluation | Current run exited `134` at the state identity; four centered rows passed their local ladder check and the state mismatch is retained in the CSV |
| F5 | `equilibrium_checks F5` path in `S10`; regular two-phase neutral ternary solved point and pinned/free-intensive neutral binary generic point | For every full primal coordinate, centered derivative of `alpha grad f + sum(lambda_j grad g_j)` equals the callback lower triangle | `alpha=(0,0.5,1)`; all callback multipliers nonzero; step factors `(10,3,1)*sqrt(epsilon)*max(1,|y|)`; regular case preserves existing `1e-6` scaled comparison, pinned case uses provisional `1e-5` on every Hessian entry in its native coordinate chart | Current run exited `1`; all 18 factor/step summaries and their largest expected/actual/error rows are retained in the CSV |

The F1 allowances and the per-case F5 allowances are detection choices for
these finite-difference tests. They are not universal numerical bars, model
adequacy limits or cross-model allowances. A passing rerun must retain its
source revision, wheel hash, exact input, branch/support, steps and signed
errors before either check is treated as evidence for a layer.

In the retained CSV, `expected` is the independent identity or centered
weighted-gradient value and `actual` is the state field or callback Hessian;
`error` is the absolute difference. For F5, `row` and `column` are zero-based
lower-triangle indices of the callback Hessian entry that produced the retained
maximum for that case, objective factor and step factor. The mixed
native-coordinate units for F5 are intentional because the callback includes
logarithmic, amount, packing, intensive and embedded-multiplier coordinates.

## Missing numerical claim map

| Layer | Missing claim | Independent reference and tolerance rationale | Current disposition |
|---|---|---|---|
| Core | Defining constants, unit factors, policies and coordinate transforms | Source-qualified defining values and dimensional round trips; derive precision from representation and scale, with boundary cases around each policy | Candidate only; no new run |
| Contributions | Ideal and residual contribution formulas, association equations and ionic/Born limits | Closed-form terms, limiting slopes and independently evaluated source equations; use operation-specific floating-point/refinement plateaus | Candidate only; no new run |
| Properties | Ideal/Residual/Total fields, association elimination, fixed-pressure implicit derivatives and zero-support reductions | F1 identities and centered ladders; support-reduced/full-state agreement only on regular domains; tolerance follows conditioning and step refinement | Initial audit failure; [focused F1 correction now passes](caloric-correction-2026-09-08.md); other claims not qualified |
| Equilibrium | Weighted curvature, phase branches, equal-fugacity roots, scaling and stability | F5 full weighted-gradient derivative; independent scalar/root or conservation relations per case; retain residual, conditioning and search scope separately | Historical source diagnosis; focused current F5 run failure; remaining cases candidates |
| Regression | Parameter recovery, observation Jacobian and weight handling | Manufactured identifiable equations with analytic predictions and independent Jacobian finite differences; limits follow parameter scale and conditioning | Capability/input work pending; no run |
| Surface and bindings | Unit/specification mapping and preservation of values/failures | Round-trip physical inputs and explicit typed failures at the boundary; these checks do not replace Properties or Equilibrium numerical evidence | Readiness evidence only |
| Applications | MEA, absorber and lithium observables and actions | Source-backed observations, declared basis/units and independent holdouts; retain source uncertainty separately from solver residuals and model discrepancy | Inputs are present in named records, but conversion and acceptance studies remain open |

This map does not promote the candidate rows in
[`test-candidates.md`](../../../engine/docs/equilibrium/test-candidates.md) to
accepted references. Numerical residual limits, source uncertainty and model
discrepancy remain separate quantities.

## Convention and provenance trace

The F1 relation is registered in `engine/docs/equations.md` and its residual
versus total-property boundary is described in `engine/docs/science/full-eos.md`.
The historical CSV retains the audit ladder; the preparation CSV retains the
new focused run. The native check prints fresh values when run.

The Khudaida Figure-2 packet has association inputs in revisions 4 and 6.
Revision 4 records isobutanol self-association and cross terms derived with the
Gross–Sadowski 2002 combining rules at `k_hb=0` for its stated temperature.
Availability does not establish source transcription, conversion or model
applicability; those comparisons remain a separate candidate study.

The draft absorber `reaction-system.json` contains MEA R1–R3. Their source
reaction directions, standard states, activity basis, reference temperature and
molality conversion still need a numerical convention check. The input record
is therefore present but not an accepted Engine reference.

## Execution boundary and next evidence

No production numerical fix was made. The focused commands run against the
owner-verified wheel/build were:

```bash
.venv/bin/meson compile -C "$(cat build/environment-wheel/meson-build-dir)" model_checks equilibrium_checks
timeout 60 "$(cat build/environment-wheel/meson-build-dir)/model_checks" F1
timeout 120 "$(cat build/environment-wheel/meson-build-dir)/equilibrium_checks" F5
```

The first command rebuilds only the two native check executables. The focused
selectors avoid the broader equilibrium suite. Both current failures are
expected until their producers are corrected. After a passing rerun, retain
the printed values/errors with the wheel hash and re-review each local
tolerance before updating any scientific status.

## Independent review

A separate read-only reviewer reproduced both focused executable outcomes and
matched every retained CSV summary, including the executable hashes. The lead
corrected the initial expected/actual labelling, required a multi-step F5 ladder,
kept the regular comparison bar, and clarified that the pinned detection bar
applies to all entries. This accepts the record as evidence of the existing
failures, not the provisional limits as future capability acceptance criteria.
