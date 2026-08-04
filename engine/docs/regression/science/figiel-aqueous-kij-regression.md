# Figiel Staged Aqueous Current-Catalog Recovery

Status:
`CONDITIONAL_STAGE_C_NUMERICALLY_CONFIRMED_PRINTED_TUPLE_NOT_RECOVERED`.
Stage B remains unchanged. Stage C now has one source-bound conditional Ceres
owner and one retained campaign result. Every conditional coordinate converged
from both declared starts, but nine of eleven fitted values missed the
user-approved `0.05` printed-value comparison. This is a valid negative
scientific result, not catalog admission, prediction, or authority transfer.

This document is the sole design and science owner for the staged Figiel
current-catalog recovery. It replaces the falsified assumption that one
fixed-family `164 x 11` fit should recover every printed Table 4/5 interaction.
The rejected fit remains immutable evidence; it is not proof that the
source-described staged procedure is impossible.

Active installed-artifact replay binds Provider frontend `0.2.0.dev0` to
commit `14fa3745264db66b8e59c12268737d694c706f2f`, tree
`eb04f10f445957cc768bad1ef4f330038c69a293`, wheel SHA-256
`48f3a75c9fc16ba71616aa703b526f41c2dcf89a7e00eebe23f75fcb8fa24594`, and
installed-header SHA-256
`2cd2b73b83c65936dff21155fd800a87b56e81cce977df7b8491ccfb2bf4c50b`.
The final Provider JUnit receipt artifact is
`provider-python-frontend-v0.2/14fa374/provider-tests.xml`,
SHA-256
`56d1d3e9fcf47bb700df4223fa2d9a20444dc97aa22b5c39525d446a296ba3cf`.
The artifact name is relative to an executor-supplied immutable artifact root.
The retained Stage B/C results below remain bound to their historical Provider
subjects; they are not rewritten as final-artifact evidence.

## Bounded claim and sequence

The smallest source-backed sequence is:

1. **Stage A — Born diameters (`5 x 5`, rank 5).** Reuse the existing first
   tracer: five active Born diameters for Li+, Na+, K+, Cl-, and Br- against
   the five SI Table S5 reported-average water-solvation targets.
2. **Stage B — water solvation factor (`21 x 1`, rank 1).** Fit one
   ion-independent `f_water` to all 21 audited Hamer--Wu NaBr MIAC rows.
3. **Stage C — conditional aqueous interactions (eleven scalar fits).** For
   each Table 4/5 coordinate, hold the other ten interactions at the printed
   tuple and fit the active coordinate to every Hamer--Wu row whose salt
   contains that interaction. Assemble the eleven fitted values and evaluate
   all 164 rows with the exact `164 x 11` Jacobian.
4. **Confirmation.** Repeat every scalar fit from the displaced start and
   require maximum parameter disagreement at most `1e-5`.

The final candidate either reproduces the printed Table 4/5 tuple within the
frozen `0.05` maximum parameter difference and passes the in-sample observable
gates, or retains
`SOURCE_DESCRIBED_STAGED_RECOVERY_DID_NOT_REPRODUCE_PRINTED_TUPLE`. Failure
does not authorize changing rows, weights, bounds, starts, or tolerances.

The eleven immutable Stage-C coordinates are:

| Column | Interaction | Published `k_ij` | Source cell |
| ---: | --- | ---: | --- |
| 0 | water--Li+ | -0.4 | Table 5, Li+ row, water column |
| 1 | water--Na+ | -0.3 | Table 5, Na+ row, water column |
| 2 | water--K+ | -0.1 | Table 5, K+ row, water column |
| 3 | water--Cl- | -0.3 | Table 5, Cl- row, water column |
| 4 | water--Br- | -0.3 | Table 5, Br- row, water column |
| 5 | Li+--Cl- | 0.8 | Table 4, Cl- row, Li+ column |
| 6 | Na+--Cl- | 0.8 | Table 4, Cl- row, Na+ column |
| 7 | K+--Cl- | 0 | Table 4, Cl- row, K+ column |
| 8 | Li+--Br- | 0.5 | Table 4, Br- row, Li+ column |
| 9 | Na+--Br- | 0.65 | Table 4, Br- row, Na+ column |
| 10 | K+--Br- | -0.35 | Table 4, Br- row, K+ column |

The explicit fitted zero is a target. Blank Table 4/5 cells are not zeros.
Inherited water--methanol and water--ethanol interactions are excluded.

## Source statements and Regression choices

The primary paper is Figiel, Yu, and Held, *Industrial & Engineering Chemistry
Research* 64 (2025) 9406--9418, DOI `10.1021/acs.iecr.5c00475`. The retained
71,826-byte Markdown has SHA-256
`ce80533925a91bc59d8d0d8056113c40611ca26c2edf04aced76986d50bd4bae`.
Durable line locators are 277, 279, and 281:

- line 277 says Born diameters were adjusted to water solvation Gibbs energies
  and solvent-specific, ion-independent `f_k` values to NaBr MIAC data;
- line 279 says aqueous ion--water and ion--ion `k_ij` values were adjusted to
  experimental aqueous-salt MIAC literature data; and
- line 281 declares the order Born diameter, `f_k`, then `k_ij`, followed by an
  iteration whose parameter changes were small.

The paper does **not** disclose an exact MIAC objective, weights, row subset,
bounds, starts, optimizer, or cycle termination. The primary-source audit in
`docs/research/figiel-table-4-5-fitting-method.md` records those unknowns.
Every numerical choice below is Regression-owned and is not attributed to the
authors.

Held et al. (2014), retained Markdown SHA-256
`b8b1e46bf870224de5de68b5989f9cb377d17445d87109a5462a94f1efaafbda`,
provides lineage evidence rather than a Figiel-specific prescription. Lines
255--268 define a squared relative-deviation objective for osmotic coefficients,
including Eq. 20, `sum(1 - phi_calc/phi_exp)^2`, and describe sequential
parameter fitting. Regression therefore uses the analogous relative MIAC
residual, while explicitly labeling that choice.

The approved Validation packet is commit
`8944d34f7002cda1bb8760e606cc1f11696f58cd`, tree
`6c8fd350dcd6bfdd7be1918f73fd33a23e2070dd`, with:

- target-ledger SHA-256
  `f405a3e48d21cd979a8dd480d5f8cb3be40754f5d6babf368b505b5f305607f0`;
- parameter-packet SHA-256
  `932e8baa90fcefbaa8c3a8730cdeadd83a4c01f0a3b109f4e4cd0319aee9312b`;
- metadata SHA-256
  `8ea06c6ca5452d01448a03f9a76cf7d0c35bb99c9abe23ccb1729d56c71d468f`;
- SI extraction SHA-256
  `85bd39f727158d5a9d6eea6828c1673f73850e783a655b09660cc9b66d84321a`;
  and
- Hamer--Wu CSV SHA-256
  `2f63e13f06a5b0f4e8bca2980b6a8d9d7fb0f839153c43e3a71952daf9796595`.

Stage C consumes only the hash-bound Hamer--Wu CSV rows; the ledger,
parameter-packet, metadata, and SI identities remain packet-level provenance
and support rather than additional residual inputs.

All states are aqueous molality-scale observations at `298.15 K` and
`100000 Pa`. Stage A uses exactly the five SI Table S5 reported averages
already frozen by the Born-tracer contract. Stage B uses every one of the 21
audited NaBr rows in the source packet (`0.001` through `6 mol/kg`). Stage C
uses every one of the 164 LiCl, NaCl, KCl, LiBr, NaBr, and KBr rows. Each row
appears once. The `<=5 mol/kg` subset is sensitivity evidence only. All rows are
training data; there is no pointwise uncertainty or approved held-out cutoff.

## Frozen residuals, variables, and solver

Stage A retains its existing residual, units, bounds, scales, starts, exact
Provider Jacobian, and closed result contract without reinterpretation.

For Stages B and C, positive observed molality-scale MIAC
`gamma_q_observed` defines the dimensionless equal-weight residual

```text
r_q(theta) = 1 - gamma_q_model(theta) / gamma_q_observed.
```

No observed value is treated as uncertainty. Equal weights are a predeclared
Regression choice because the packet has no rowwise uncertainty.

Stage B has one dimensionless variable `f_water = z`, scale `1`, bounds
`[1,2]`, and starts `1.2` and `1.8`. Published `1.5` is comparison-only and is
not a seed. Stage C carries forward the exact fitted Stage-A diameters and the
primary Stage-B value `f_water = 1.5590515389548207`; their evidence subjects
are part of the immutable Stage-C specification. Each Stage-C scalar problem
has one dimensionless variable `k_j = z_j`, scale `1`, bounds `[-1,1]`, and
starts `0` and `0.25`. The other ten interactions are fixed to the published
tuple as the explicit conditional context. The active published value is
comparison-only and is not a seed, target, prior, or regularizer.

Every stage uses the existing Ceres owner with `DENSE_QR`, one thread, silent
logging. Stage C uses at most 50 iterations and function, gradient, and parameter
tolerances `1e-10`. The `180 s` operational deadline applies independently to
each scalar Ceres solve, not to the eleven-coordinate schedule as a whole.
Stage C has one closed result contract in the existing workflow and reuses the
sole Ceres engine, native module, and CMake target.

## Exact Provider derivative contracts

Stage A consumes the existing model-bound active-Born value/first-total-
derivative callback.

Stage B requires one appended model-bound callback for an ordered aqueous
`(water, Na+, Br-)` model. At fixed `T`, fixed `P`, and formula-unit molality,
it accepts one finite trial `f_water` and returns
`ln(gamma_pm^m)`, exact total fixed-pressure
`d ln(gamma_pm^m)/d f_water`, reference convergence diagnostics, parameter
fingerprint, component order, and structured status. Provider retains EOS,
Born/association/electrolyte reference sequence, density closure, and CppAD
ownership. Regression derives

```text
gamma_model = exp(log_gamma_model)
dr/df_water = -(gamma_model/gamma_observed) * dlog_gamma_model/df_water.
```

The retained Stage-C result consumed the bounded historical Provider callback
from commit
`8af6e5467cba70ade226cb28f93688ce88048c29`, tree
`2ec243429600de75ccaeb51955044bf2e7557fbe`. The retained wheel SHA-256 is
`1d77a46aca369269fd97b32a722e8b548617a4ae840e07645040a138f63863db`;
the installed-header SHA-256 is
`555ebde7aa2dcfc7a41a1a8f14af8a7ede678caafd24aa93cd165d291325fb4d`.
For each ordered `(water,cation,anion)` row it returns bounded
`ln(gamma_pm^m)` and exact total fixed-pressure derivatives with respect to
`(k_water_cation,k_water_anion,k_cation_anion)`. Regression maps those three
entries into the global eleven-column Jacobian and computes

```text
dr/dk_j = -(gamma_model/gamma_observed) * dlog_gamma_model/dk_j.
```

All other row entries are structural zero. Production numerical derivatives,
copied EOS/reference equations, independent density closure, and Equilibrium
dependencies are forbidden.

## Derivative, rank, and confirmation gates

Installed-artifact derivative checks use centered callback-value differences at
steps `1e-4` and `5e-5`. Each active column must meet

```text
abs(J_exact - J_h/2)
  <= max(1e-8, 20*abs(J_h - J_h/2), 2e-8*abs(J_exact)).
```

Finite differences are evidence only. They are not a runtime backend.

For each assembled result, SVD rank uses

```text
s_max * max(residual_count, parameter_count) * epsilon_binary64 * 100.
```

The assembled Stage-C `164 x 11` Jacobian must have rank 11, and each scalar
problem must have a finite nonzero active column over its affected rows.
Results retain singular values, condition number, complete/nonzero columns,
active bounds, and the least-sensitive direction. The fixed published context
defines the approved conditional question; it is not used to manufacture the
assembled rank.

Both declared Stage-C starts must converge to the same eleven-coordinate
assembled result within `1e-5` in infinity norm. There is no coordinate sweep,
outer A-B-C cycle, or published-value initialization.

## Retained rejected alternatives

`evidence/figiel-aqueous-kij-published-tuple-preflight.json`, SHA-256
`5bd86e332b94781112eeee0ca06765a0f084020a30af76169861bbc610d5743d`,
retains the exact isolated logarithmic-residual result. Its installed
`164 x 11` Jacobian has rank 11; all 492 declared derivative checks pass; the
three starts agree within `6.740175084729572e-11`; fitted cost is
`0.25101017331848846` versus published cost `0.4200114112320464`; and the
maximum published-parameter difference is `1.8` with five lower bounds active.

A read-only discriminator using the frozen relative residual also selected an
incompatible boundary optimum: all-164 published AARD `5.13246%` versus fitted
`4.19795%`; `<=5 mol/kg` published AARD `4.68197%` versus fitted `3.46284%`;
maximum parameter difference `1.8`. These results justify the staged source
contract. Neither establishes global impossibility or author error.

## Standalone Stage-B result

The two declared starts converged in seven Ceres iterations to
`1.5590515389548207` and `1.559051762650834`, a maximum absolute difference
of `2.2369601326843735e-7`, below the frozen `1e-5` start-agreement gate. The
`21x1` exact Jacobian has rank 1, no bound is active, all 21 rows are evaluated,
and the primary MIAC RMSE is `0.04582795014810974`. Solver, numerical,
physical, and workflow statuses pass independently.

The Table 5 value `1.5` is printed to one decimal and was never a hidden seed.
The fitted difference `0.05905153895482074` is sufficiently close for the
user's engineering comparison. This is descriptive adjudication, not a newly
invented acceptance tolerance and not experimental uncertainty.

## Conditional Stage-C result

The retained bounded Provider subject is commit
`8af6e5467cba70ade226cb28f93688ce88048c29`, tree
`2ec243429600de75ccaeb51955044bf2e7557fbe`. Its wheel and installed-header
SHA-256 values are
`1d77a46aca369269fd97b32a722e8b548617a4ae840e07645040a138f63863db`
and `555ebde7aa2dcfc7a41a1a8f14af8a7ede678caafd24aa93cd165d291325fb4d`.
The primary fitted tuple is

```text
(-0.5364822745, -0.3292520491, -0.1579732178,
 -0.4055029787, -0.3743902294,  1.0000000000,
  0.9394025030,  0.1342006022,  0.7607943632,
  0.6520078999, -0.1743688955)
```

All eleven primary and confirmation scalar solves terminated with usable Ceres
convergence. The maximum start disagreement is
`1.3809200687386891e-6`. Two printed parameters are within `0.05`; nine are
not. The largest difference is Li+--Br- (`0.2607944`), and Li+--Cl- reaches
the declared upper bound at `1.0`. The assembled Jacobian has rank 11 and
condition number `40784.83209599819`; all 1,804 entries pass the retained
two-step centered derivative check.
Consequently the retained status is
`SOURCE_DESCRIBED_STAGED_RECOVERY_DID_NOT_REPRODUCE_PRINTED_TUPLE`.

## Result and status semantics

The standalone Stage-B result retains source/artifact identities, input and
evaluated row IDs, both starts, the fitted coordinate, row residuals,
objectives, rank, active bounds, Provider diagnostics, and ordered failure
reasons. The Stage-C result remains separately reviewable and makes no
catalog-authority claim.

Statuses remain independent:

1. `solver_converged` is reported per Ceres stage/start.
2. For standalone Stage B, `numerically_converged` requires finite diagnostics,
   rank 1, and declared-start agreement within `1e-5`; the installed-artifact
   test separately checks the exact Jacobian at both declared centered-
   difference steps. Stage C separately requires its two-start conditional
   agreement gate.
3. `physically_valid` requires every fixed-pressure state to be Provider-usable
   with positive finite MIAC. Born-tracer physical gates remain owned by Stage
   A and are not silently re-run by the standalone Stage-B result.
4. `workflow_valid` additionally requires exact source hashes, row membership,
   units/basis, state, salt/component order, Provider artifacts/fingerprints,
   fixed-family inputs, and complete input/evaluated/failed accounting.
5. The standalone Stage-B comparison is descriptive and user-adjudicated.
   Stage-C `scientifically_valid` additionally requires every
   frozen in-sample observable gate and
   `max_j abs(k_j_fit-k_j_published) <= 0.05`.
6. `predictive_status` remains
   `NOT_ADJUDICATED_NO_APPROVED_HELD_OUT_CUTOFF`.

The `0.05` comparison is a user-approved engineering recovery criterion, not
source uncertainty or a uniqueness statement. Existing in-sample forward
checks remain pooled MIAC RMSE `<=0.17`, per-salt RMSE `<=0.35`, per-salt
maximum absolute MIAC error `<=1.25`, and first predicted MIAC for each salt
`<0.98`. A failure remains a valid scientific result but not a valid recovered
catalog candidate.

## Ownership, next gate, and negative space

Regression owns targets, residuals, Ceres execution, diagnostics, and the
Stage-B and conditional Stage-C results. Provider owns values, exact
derivatives, model records, reference sequences, density closure, and bounded
callback cancellation. Provider implementation `8af6e5467cba70ade226cb28f93688ce88048c29`
supplies the installed bounded batch seam used by Stage C.

The package checkpoint retains the failed comparison without tuning. A later
installed-artifact Validation campaign may replay the exact conditional result;
Regression evidence cannot self-promote. Recovering the historical printed
tuple from experimental rows remains unresolved because the paper omits the
objective, weights, exact rows, optimizer, bounds, starts, and sub-staging.
Artifact validation is intentionally a matrix: the Stage-C wheel with Provider
`8af6e54` passes seven routine Stage-C tests plus the explicit public-fit
campaign, while the accepted legacy Provider artifact preserves the 49
accepted pure-workflow tests. The public campaign is replayed with
`python -m pytest -o addopts='' -q
tests/test_figiel_aqueous_kij.py::test_public_aqueous_kij_fit_replays_retained_negative_result`;
it is excluded from the routine suite because it executes all 22 scalar Ceres
solves. Running the entire historical suite against `8af6e54` is diagnostic
only; its three known pure-workflow numerical differences are not Stage-C
acceptance failures.

Excluded: generic registries, mutable parameter overlays, Provider catalog
writes, compatibility shims, simultaneous all-table solves, organic-solvent or
expanded-ion scope, density/osmotic targets, association/polar/reactive/MEA
scope, uncertainty, prediction, global identifiability, an alternate solver,
result, module, or target, and any runtime dependency on Zotero, Validation source trees,
Migration, or lab code.
