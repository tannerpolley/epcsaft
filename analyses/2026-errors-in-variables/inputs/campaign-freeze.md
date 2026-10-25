## 9. Synthetic MEA/ionic bias experiment

Freeze the following experiment in the accepted plan before Build. These are
design inputs, not measured performance or evidence of a successful fit.
No parameter/noise/grid/seed tuning after viewing a fitted comparison.

### MEA case

- Use the parameter and reaction declarations in
  `analyses/2026-fitting-state-reuse/inputs/mea-initial.json`, originally selected
  at `1ea79ab29d1c5d66e0e13113ab2acd4651169a0a` and verified byte-identical at
  main `ae6a8f11` (**R4**), and its first 16 states in stored order
  (`vle_obs_0119` through `vle_obs_0134`). They specify T=313.15 K and a declared liquid plus an
  incipient neutral vapor. Replace all retained measured targets with
  native-generated CO2 partial pressures at the known truth. Preserve each
  state's phase support, reaction/source references, feed total, guesses and
  solve criterion. Free pressure is an output, not a noisy pressure input.
- Fit only the first stored model coordinate: symmetric Kij for water and
  carbamate, whose stored origin is -0.06655083320547248, scale 0.01 and
  bounds [-0.5,0.5]. Its physical truth is -0.06655083320547248, verified in
  `parameters.pairs` at `pair/carbamate-anion/water/k_ij`. Freeze the complete
  supplied parameter bundle as truth. Other stored fit coordinates stay at their supplied
  physical values, including the stored active-bound values. Their historical
  fitted status is not physical validation of this synthetic truth.
- Independent noisy input coordinates are T (sigma=1 K) and feed CO2 mole
  fraction (sigma=0.002 mol/mol). The latter replaces water at fixed total,
  leaving all other species, including the charge-balanced ionic feed,
  unchanged. Use the physically complete one-dimensional neutral exchange
  bounds 0 <= x_CO2 <= x_CO2,true+x_water,true; temperature bounds
  [303.15,323.15] K. These are declared qualification domains, not sigma
  clipping. The measured noisy coordinates may lie outside the physical box;
  record that as an invalid experiment, never redraw or clip.
- sigma_y,i is 0.01 times the absolute native true partial pressure, fixed
  before generating measurement errors; it is recorded in Pa. Refuse an
  unavailable or zero true signal. Start the model coordinate at truth+0.01,
  within its declared bounds; every latent input starts at its measurement.

### Ionic case

- Use the existing three-species native synthetic aqueous-electrolyte fixture
  `properties_support::aqueous_electrolyte`, with the same explicit changes as
  `coupled_lengths(3.2)` in `regression_f1_j1_checks.cpp`: cation sigma=3.2
  angstrom, packing and Debye-Hueckel diameters=0.88 times sigma, Born
  diameter=3.5 angstrom. These are declared synthetic parameters.
- Fit that single tied length coordinate using `coupled_coordinate`: origin
  3.0 angstrom, scale 0.2 angstrom, bounds [2.4,4.4] angstrom. All other
  parameters remain fixed; no hydration hold is introduced in this experiment.
  Start at 3.0 angstrom.
- Sixteen pure-solvent salt states have T=298.15 K, P=1e5 Pa and
  m_j=0.1*(45)^(j/15) mol/kg for j=0,...,15. Generate both molality-scale
  mean ionic activity coefficient and osmotic coefficient natively at each
  true state. Share one measured m/T/P group between those two outputs.
- Independent input sigma values are 0.5 K, 2000 Pa and 0.1 m_j mol/kg.
  Declared domains are [288.15,308.15] K, [5e4,1.5e5] Pa and [0.001,8]
  mol/kg. Each output sigma is 0.005 times its absolute true value, frozen
  before noise generation. Keep the existing reference definition and witness.

### Repeated measurements, controls and decision

For each case generate R=128 independent replicate datasets, with independent
unit-normal draws for every input measurement and output measurement. Use
native `mt19937_64`, seeds 51001 (MEA) and 51002 (ionic). Retain the actual
standard-normal draws and noisy physical values, since a library's normal
distribution transform is not a portable identity. Generate once; both
estimators consume those same recorded datasets and sigma values. No clipping,
truncation, sorting, redrawing or selective replicate deletion is allowed.

At each replicate run:

1. The current fixed-input Difference estimator using measured inputs and
   w_i=(scale_i/sigma_y,i)^2.
2. The augmented estimator using those same measurements and output sigma,
   plus the declared input sigma.

Use identical parameter starts, bounds, Ceres tolerances and controls:
100 iterations and 600 elapsed seconds per fit, the existing native-check
budget. **R3 — resolved controls:** explicitly set main function, gradient and
parameter tolerances to 1e-12, 1e-12 and 1e-10. These equal main defaults at
`core/numerical_policy.hpp:276–289`; zero at `fit.cpp:925–927` selects a policy
default and must not be multiplied to form a replay. Do not loosen stopping,
retry with alternate starts or borrow unused budget from the competing estimator.
Record all state/action evaluations, trial refusals, iterations, elapsed times,
achieved optimality and bound/rank
diagnostics. State-evaluation work may increase; this issue claims less bias,
not speed. No-input ordinary key-solve work must remain unchanged.

For method A define the signed optimizer-scaled parameter error
e_A,r=(theta_A,r-theta_true)/s_theta, using s_theta=0.01 for MEA and
0.2 angstrom for ionic. Estimate signed bias b_A=mean_r(e_A,r), its Monte
Carlo standard error sd(e_A)/sqrt(R), absolute bias |b_A| and paired absolute
bias difference Delta=|b_EIV|-|b_fixed|. Also report parameter RMSE separately;
it is not the bias acceptance metric.

**R1 — simultaneous signed-mean confidence set:** freeze this rule before any
fitted comparison. Keep 2000 paired bootstrap resamples of replicate indices
(native seed 51003, separate recorded index streams per case); both methods use
the same indices within a resample. Bootstrap signed means b_A, never Delta or
an absolute mean. For each of four signed means (two methods in two cases),
form a two-sided percentile interval I_A=[L_A,H_A] at quantiles 0.00625 and
0.99375: nominal error 0.05/4=0.0125 per interval. An empirical quantile is the
sorted resampled mean number ceil(2000 p), with one-based indexing. Their
Cartesian product is the simultaneous signed-mean confidence set before taking
absolute values. The union bound allocates 5% total error if the four marginal
intervals attain nominal coverage; method independence is not required.

For each case compute

$$
U_\Delta=\max(|L_E|,|H_E|)-d_F,\qquad
 d_F=\begin{cases}0,&L_F\le0\le H_F,\\
 \min(|L_F|,|H_F|),&\text{otherwise}.\end{cases}
$$

Require U_Delta < -epsilon_num in **each** case. Report signed intervals, point
Delta, U_Delta, numerical allowance and raw paired errors even when it fails.
An inconclusive set prevents acceptance; do not change R, resample count or
interval rules after viewing comparisons. Signed means are smooth even at zero
bias; the ordinary percentile bootstrap is applied before the nonsmooth
absolute-value transformation. Coverage remains approximate for R=128:
skewness, heavy tails and only 2000 resamples can affect tail quantiles. No
measured finite-sample coverage or power is claimed, and the union bound does
not turn approximate marginal coverage into an exact 95% guarantee. This
resolves `readiness-r1.md` R1 (Fang and Santos, Example 2.1/section 3.2) without
changing the estimator, datasets or noise model.

**R3 — scaled-parameter numerical allowance:** include two contributions;
tighter Ceres stopping alone is insufficient.

1. Replay the unchanged indices **0, 32, 64 and 96** for both methods and cases
   with explicit function/gradient/parameter tolerances **1e-13, 1e-13, 1e-11**,
   the same starts, bounds, data, 100 iterations and 600 seconds per fit. Let
   t_A=max |theta_replay-theta_main|/s_theta over those four replays. Retain
   cold same-point endpoint evaluations prescribed by ALG-FIT; tighter stopping
   cannot remove seeded/cold differences.
2. At every main fit endpoint and true-data generation state, retain physical
   prediction-error estimates nu_i and raw action-error estimates from existing
   state/reference/action diagnostics. Include density-closure/association error,
   equilibrium defining-row residuals propagated through its checked defining
   Jacobian, and terminal reference rung/witness changes. Contract ionic
   log-fugacity errors with stoichiometric weights and section 5's activity
   formulas; propagate solvent ln(a_s) error to osmotic coefficient by division
   by s=nu M_s m. Retain parameter and raw T/P action differences, then apply D,
   sigma_q and 1/sigma_y exactly as for fitted columns. Report state,
   reference-value and reference/action contributions separately. A threshold,
   missing diagnostic or successful solve is not an achieved error estimate;
   unresolved consumed prediction or action error prevents qualification.

At a regular interior endpoint let r and J be the complete standardized
residual and augmented Jacobian (parameter columns only for the fixed method).
Let eta_i=(nu_i,fit+nu_i,true)/sigma_y,i on output rows, zero on exact input-prior
rows, holding the recorded noisy data and frozen sigma fixed. Let E_ia estimate
the absolute numerical error of each standardized J_ia from the raw state/action
diagnostics; prior identity rows have E=0. Propagate to parameter component a
using the existing SVD and the local least-squares sensitivity:

$$
H=J^TJ,\qquad
h_a=\left|e_a^T H^{-1}\right|
       \left(|J|^T\eta+E^T|r|\right).
$$

Absolute values are componentwise and products are matrix products. Solve
through the SVD, without an explicit inverse or second covariance owner.
In this one-parameter campaign v_a already uses s_theta, so h_a has the same
scaled units as signed errors. Use the Gauss–Newton approximation to the
stationarity perturbation H delta z=-(J^T delta r+delta J^T r); this is a local
error estimate,
not a nonlinear/global bound. Report each contribution to h_a, including the
terminal-reference component. On the four replay endpoints include cold versus
seeded prediction/column discrepancies in eta/E by taking the larger of each
observed discrepancy and its diagnostic estimate. A contradiction of those
estimates needs diagnosis before qualification, not an assumption that
unexamined replicates are safe.

Let n_A=mean_r h_A,r over all 128 main endpoints of a case. Set

$$
\epsilon_{\rm num}=2\left(n_E+n_F+t_E+t_F\right).
$$

The sum covers both methods' signed-mean sensitivity before absolute values;
the factor two is a declared planning multiplier, not a proved bound. Report
n_A and t_A separately and apply this case-specific allowance to U_Delta.
Endpoint propagation is retained for all replicates; the four tighter-fit/cold
replays remain an observed sensitivity check, not a bound on every replicate.
The admitted terminal reference remains the synthetic target: rung differences
estimate numerical/reference sensitivity, never physical-model error or a reason
to switch to an extrapolated reference. Unfinished fits, unresolved error,
changed branches, active fitted bounds or rank deficiency prevent qualification.
Freeze the numerical allowance rule before the fitted comparison.

Retain every failed replicate and its cause. The predeclared primary comparison
requires all 128 paired fits per case to converge regularly and remain in the
declared domain; a failure prevents a passing claim. A complete-pair descriptive
table may still be shown with its changed sample count and selection limit,
but cannot replace that primary comparison. Numerical failures go to diagnosis
before changing a solver, tolerance, reference, noise level or experiment.

The finding, if successful, is numerical verification under a known synthetic
noise model in these two domains. It does not establish reduced bias for all
datasets, measured thermodynamic validity or general errors-in-variables
statistical consistency at finite sample size.

## R3 error-estimate prescription (lead, 2026-10-15)

This answers the builder's question on how to build nu_i and E_ia in §9. The existing owners expose residuals, not achieved errors. Freeze this together with the R1 and R3 rules before any fitted comparison. It adds no solver, no tolerance setting, no right-hand-side type and no derivative owner. Terms are formed on raw values, then scaled by D, sigma_q and 1/sigma_y exactly as for the fitted columns. They are summed, which is conservative.

At every main fit endpoint and every true-data generation state, for output row i and column a:

$$
\nu_i=\nu_i^{P}+\nu_i^{\rm ref}+\nu_i^{\rm obs}+4\varepsilon|f_i|,\qquad
E_{ia}=E_{ia}^{\rm lin}+E_{ia}^{\rm ref}+E_{ia}^{\rm obs}+4\varepsilon|J_{ia}|,\qquad \varepsilon=2^{-52}.
$$

- **nu^P, the density closure.** |df_i/dP|_{T,n} times the achieved absolute pressure residual of the density root. Take the derivative from the existing pressure-direction property action. Take the residual from `PropertyActionDiagnostics::pressure_residual`; confirm in source whether it is relative and multiply by P if so. This is the first-order Newton correction of the root.
- **nu^ref.** The reference rung and witness term already specified in §9.
- **nu^obs, E^obs: observed same-point discrepancy.** Take |seeded − cold| for each prediction and column. "Seeded" is the fit's own evaluation; "cold" is the same-point evaluation ALG-FIT already prescribes, now taken at **every** endpoint and true-data state, not only the four replays. This is the achieved estimate for every inner solve whose start differs between the two paths: density, association and equilibrium states. The equilibrium defining-row residual is covered here. Do not add a J^-1 r propagation.
- **E^lin: the linear-solve forward bound.** Use 2·kappa·eta·|J_ia| for each linear action solve on the column's chain (the bound named at `implicit_sensitivity.cpp:404`):
  - solved-state actions: kappa = `condition_estimate`, eta = `maximum_linear_residual`;
  - pressure-closed property actions: eta = that path's pressure-chain residual (relative, `helmholtz_model.cpp:1328`) and kappa = 1, because it is a scalar solve;
  - association site actions: eta = `association_residual` and kappa = `site_pivot_ratio` (maximum over minimum pivot, `helmholtz_model.cpp:121`).

  Sum the bounds over the solves on the chain.
- **E^ref: the reference-action term (lead, answering the builder 2026-10-15).** The raw terminal-to-preceding-rung, or witness, differences of the parameter and T/P reference actions that §5 and §9 already require. It is a separate retained term, scaled by D, sigma_q and 1/sigma_y like every column, and reported separately. It is not folded into E^lin.
- **An inner solve on a row's path without its diagnostic leaves that row unresolved.** §9 qualification then fails as written. A term is zero only when that solve does not occur on the row's path.

Retain every term, by row and endpoint, and report the components of n_A in the notebook. This is a local first-order estimate, not a bound. The factor two in epsilon_num stays a declared planning multiplier.

**§8.1 exact-answer check stopping (lead, 2026-10-15).** Only the exact-answer linear-reference verification may use the tighter explicit controls 1e-14/1e-14/1e-12. Run a bounded stopping-refinement probe first: main defaults, then 1e-13/1e-13/1e-11, then 1e-14/1e-14/1e-12. Retain every step, showing the error in each checked quantity falling as stopping tightens. That shows the default-control miss (q0 off by 3.36e-7 K against a 3.0e-8 K allowance) is stopping error, not a formulation error. The 1e-10 max(1,|reference|) criterion is unchanged. If the tightest controls still miss it, report the achieved gradient norm and the Newton-step estimate H^-1 g; don't loosen the criterion. The campaign's frozen controls and the replay tolerances (1e-13/1e-13/1e-11) are untouched.

**§8.1 shared-output plateau (lead, 2026-10-15).** The builder reported the following, retained in `build/issue-510/linear-refinement-complete.log`:
- Independent groups pass at 1e-13 and 1e-14.
- The shared-output case stops at q = 320.00000027 K at all three settings. That is 2.72e-7 K off, against a 3.2e-8 K allowance. The achieved gradient is |J^T r|_inf = 2.66e-6 and the Gauss–Newton step is |(J^T J)^-1 J^T r|_inf = 1.36e-7.

Identical endpoints at three stopping settings mean these tolerances are not what stops it. For a linear least-squares problem, one Gauss–Newton step from any point lands on the optimum, yet this step is exactly half the error. So diagnose (cse:diagnose) before changing anything:
1. Record Ceres's termination type and message and its per-iteration log (accepted and rejected steps, actual and predicted cost change) at each setting. Run one extra analytic-only probe at function tolerance 1e-16 with the other tight controls.
2. At the stalled endpoint, compare the code's residual vector and Jacobian entry by entry with the analytic residuals and Jacobian of the shared-output problem. Also compare them with central differences of the code's own residuals (verification only).
3. Compute the q component of the code's Gauss–Newton step and of the analytic one from that endpoint.
   - If they differ, for example by a factor of 2, it is a build defect in shared-output handling. Candidates: a prior counted per row instead of once per measurement, a latent column assembled twice, or a membership map that disagrees. Fix it in its owner, then rerun §8.1 at the original controls.
   - If the code's residuals and Jacobian match the analytic ones and the code's own Gauss–Newton step lands within the 1e-10 allowance, the stall is cost-roundoff stagnation. The check then passes on that one exact step from the endpoint, with the Ceres endpoint and the stall retained.
4. Don't loosen the 1e-10 criterion and don't change the analytic reference without showing it is wrong. Campaign fits stay paused until §8.1 passes.

**Analysis allowance revised (owner, 2026-10-15).** Analysis scripts may grow by up to 700 net lines, up from 250–350. This covers the native campaign driver, its runner and its result reader in the one analysis folder. Engine runtime allowances, method, data, seeds, controls and criteria are unchanged.

**R3 diagnostic gaps (lead, 2026-10-15).** The builder's source audit found three gaps:
- the pressure-chain linear residual isn't exposed on `pressure_actions`;
- the MEA first-actions chain exposes no association or source-reference action diagnostics;
- Equilibrium has zero growth allowance.

These substitutes replace "unresolved" for those rows. They use only existing owners and retained verification evidence: no new diagnostic definitions, no Equilibrium change, no allowance change.
1. **Density root.** nu^P = |df/dP|_{T,n} times the exposed absolute root residual in Pa. Do not multiply by P.
2. **Pressure-chain column solve.** In `pressure_actions` it is one scalar division by dp/drho, so its forward error is rounding. The 4 eps |J| floor covers it, and E^lin is 0 for that step. This is not a missing diagnostic.
3. **Columns on chains without kappa and eta for an inner solve** (association and reference actions on the MEA first-actions chain): E_ia += rho_f |J_ia|. rho_f is the largest relative error estimate for that column family from the exact-against-central-difference verification:
   - For each column, take the smallest |J_exact − J_FD(h)| over the frozen steps (0.02/0.01/0.005), plus |J_FD(h) − J_FD(h/2)| as the difference error at that step.
   - This is an observed relative estimate and generally an overestimate. Report it as such.
   - If the frozen derivative check doesn't cover a family, extend that verification run to cover it before any fitted comparison. This is a verification use of differences.
4. **Reference values on the MEA chain.** If the chain consumes a laddered reference, the native driver calls the existing reference owner directly at the same reference states to get rung and witness differences for nu^ref and E^ref. The MEA chain itself is unchanged. If the chain uses only the closed-form source law for ln K, then nu^ref and E^ref are 0 beyond the floor. State which applies.
5. **What remains unresolved.** Only a row whose term is still unavailable after 1–4 stays unresolved.

Re-freeze the allowance rule with these items before any fitted comparison.

**rho_f normalization (lead, 2026-10-15).** Use the column max norm, as the builder proposed:
- For each standardized residual column c and frozen step h in {0.02, 0.01, 0.005}: d_c(h) = max_i(|J_ic − FD_ic(h)| + |FD_ic(h) − FD_ic(h/2)|). The h/2 evaluations are for verification only.
- rho_c = min_h d_c(h) / max_i |J_ic|.
- rho_f is the maximum of rho_c over the columns in a family: model parameters, T, P, feed composition, molality.
- An exactly zero column with a nonzero discrepancy is unresolved. An exactly zero column with zero discrepancy contributes zero.

Apply it column-wise at each endpoint: E_ia += rho_f · max_i |J_ia|. This is the column norm at that endpoint, not rho_f·|J_ia|, so small entries aren't credited with smaller errors than the column's own scale supports.

**MEA reference propagation (lead, 2026-10-15).** The builder found that all five MEA reactions consume laddered neutral reaction references. Confirmed propagation:
- Call `reaction_references` at each cold endpoint and take its terminal-limit changes delta C_k. These must be in ln K units, the same as the source-law intercept A_k; contract per-species reference values stoichiometrically if the ladder reports them per species.
- Take S_ik = df_i/dA_k from the existing `solved_state_first_actions`, with one intercept direction per reaction.
- nu_ref,i = sum_k |S_ik| delta C_k.
- E_ref,ia = sum_k |S_ik| delta(dC_k/dz_a), scaled by D or sigma_q and by 1/sigma_y.
- Retain every raw terminal change and every propagated contribution.

The solved-state linear bound and the frozen rho_f cover the remaining chain action error. The cross term (d^2 f/dz dC)·delta C is not estimated. State that in the notebook as a limit of this first-order allowance; it is covered only nominally by the declared factor of two. This is analysis code only; no derivative owner is added.

**Density-root terms (lead, 2026-10-15).**
- **Ionic nu^P.** Sum over every contributing root of |that root's pressure action| times its exposed absolute residual: finite and terminal-reference roots for gamma; finite and pure-solvent roots for osmotic.
- **MEA.** The equilibrium endpoint's volume closure is covered by nu_obs. The hidden reference-trace roots are represented by the retained terminal-limit changes in nu_ref plus nu_obs, with no separate density term. The notebook states that these root residuals are not separately estimated.

With this, the R3 rule is settled; freeze it. Gaps found later take the closest rule here and are recorded for delivered review.

## Requalification on main (#624, lead, 2026-10-17)

Frozen before any requalification fit. The 2026-10-14 signed-bias summaries have been inspected, so nothing here changes a rule that decides the comparison. Seeds, domains, sigma rules, starts, controls, replicate count, bootstrap rule, the numerical-allowance rule and the decision rule are unchanged.

**Superseded clauses.** Two frozen clauses are superseded by the Engine, not by this analysis:
- §9 Ionic "Keep the existing reference definition and witness": #594 replaced the raw 1e-12 mol/kg reference with the √m limit (10 f(1e-14) − f(1e-12))/9 of the 1e-12/1e-14 pair. The ionic γ± truth falls by 1.180943e-6 relative at every molality; φ and the MEA pressures move only at rounding level.
- §9 "Retain the actual standard-normal draws and noisy physical values": the noisy outputs are regenerated from the new truth with the identical retained standard-normal draws, admission flags and bootstrap index streams, which `prepare` reproduces bit for bit on main.

The 2026-10-14 notebook calls such a regeneration "a different experiment". This requalification is that experiment, run under the unchanged rules above with the comparison rule below frozen first. #623 corrected the binary128 test oracle to the same limit, so `prepare` verifies both cases on main.

**Kept reference-value term.** ν^ref, E^ref and the parameter reference term keep the frozen raw 1e-10→1e-12 rung pair and the frozen witness selection; rung differences remain sensitivity estimates, never a reason to switch reference. The freeze forbids redefining a sensitivity estimate after the comparisons were viewed. No allowance has been defined or measured for the √m-limit reference; any such allowance would be a new comparison outside this requalification.

**Principle for the driver corrections.** The density-root term ν^P is an achieved error of the value the Engine returns. It therefore charges the roots that returned value is computed from, with weights equal to that value's dependence on each root: the "every contributing root" rule above applied to main's reference. A literal reading on main of the reference terms, the 1e-12→1e-14 pair, would give ionic ε_num 0.0113 instead of 0.1128; both leave the ionic decision inconclusive against U_Δ = −0.00309. The three corrections in `scripts/campaign.cpp`:
1. γ± ν^P: the terminal-reference root term charges the 1e-14 root at the 1e-14 composition with weight 10/9 and the 1e-12 root at the 1e-12 composition with weight 1/9. On main the line paired the 1e-14 composition with the 1e-12 root.
2. ν^ref selection: the witness is used only when the raw 1e-10→1e-12 change exceeds the reference convergence tolerance, as at freeze. On main `convergence_error` is the limit's change.
3. E^ref input-action selection: the witness is used only when σ times the raw 1e-10→1e-12 action change exceeds that tolerance, as at freeze. On main the witness is always evaluated, so the line charged the 1e-12→1e-14 change.

**Comparison rule.** Old values are the 2026-10-14 results, read from `git show 3fbccd954:analyses/2026-errors-in-variables/results/<case>/…`; τ = 1e-4 scaled units.

| # | Quantity | Unchanged when |
|---|---|---|
| C1 | every `paired-errors.csv` signed error, 2 cases × 2 methods × 128 | max_r abs(e_new − e_old) ≤ τ; every regular flag 1 |
| C2 | `signed_mean`, `signed_lower`, `signed_upper`, `rmse` per method and case | abs(Δ) ≤ τ |
| C3 | `U_Delta`, `delta` per case | abs(Δ) ≤ 2τ |
| C4 | `mcse` per method and case | abs(Δ) ≤ τ |
| C5 | `epsilon_num` per case | ratio new/old in [0.8, 1.25] for MEA, [0.99, 1.01] for ionic |
| C6 | `accepted`, `failed_main_fits` | MEA 1, ionic 0; 0 failed main fits per method |
| C7 | ionic `numerical-components.csv`, per method | `reference_value` ratio new/old in [0.99, 1.01]; `reference_action` and `pressure` ratio new/old in [0.5, 2] |

C2–C4 follow from C1 because the bootstrap indices are identical: each resampled mean moves by at most max_r abs(Δe), order statistics are 1-Lipschitz in that norm, and U_Δ combines two interval ends. `results/requalification.csv` records old, new and Δ for every quantity.

**Stop conditions.** C1, C5 or C7 failing, any irregular fit (including a fit stopped by its 600 s limit) or any unavailable numerical term stops the work for diagnosis before anything is reported. No replicate is dropped, redrawn or re-seeded. A decision flip is reported as a finding, not as a requalification.

**Fifth driver correction and the refused attempt (lead, 2026-10-17).** The first round-3 `run` (2026-10-17T10:35:57Z, frozen inputs `0cb5a4ce…` MEA draws) fitted nothing: all 528 fit records, 2 cases × 2 methods × 128 main fits plus the replays, were refused as `invalid_setup` before any iteration ("maximum_iterations and maximum_elapsed_time_seconds must be set and positive"). The driver's `controls()` assigned `request.controls` by position. At the freeze Engine `6b879f27` `FitControls` was {function, gradient, parameter tolerance, maximum_iterations, maximum_elapsed_time_seconds, workers}; #486 replaced the three tolerances with the six-entry `solver_options` array (association, density, row, function, gradient, parameter), so the same braces put 1e-12, 1e-12, 1e-10, 100, 600, 1 into those six tolerances and left both limits at 0. The correction assigns the controls by name: function and gradient tolerance 1e-12 (replay 1e-13), parameter tolerance 1e-10 (replay 1e-11), 100 iterations, 600 s, 1 worker. That restores the frozen controls exactly; the association, density and row tolerances are unset and take their policy defaults through `SolverOptions::value`, as they did at freeze. The misplaced tolerances also reached `prepare`'s MEA observation solves, so the MEA inputs were regenerated with the corrected driver: ionic files byte-identical; MEA truth, sigma and measured outputs move by at most 5.86e-14 relative and verification differences by at most 3.4e-9; normal draws, admission flags, input rows, ρ and bootstrap streams identical. F1b and F2 hold against the 2026-10-14 inputs (worst step change 0.78 of its bound; MEA model error/allowance 1.248e-3). The refused attempt's records are archived in the main checkout's `build/issue-624/build-r1/failed-run-results.tar.zst`; its statistics are all NaN and nothing in it was inspected as a result.
