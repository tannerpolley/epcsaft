# MEA reactive refit through the Engine fitter: derivatives, Jacobian cost, slow states

**Frozen record — not rerunnable as a current-main calculation.** Record snapshot:
commit `ca719012c31e8b41edcb98ae82a388117b58d343`. Original wheel identities, source
inputs and retained numerical evidence keep their recorded provenance. The comparison
requires the pinned MEA application checkout, its calibration records and the original
before/after wheels; it is not a portable current-main replay.

Question (#152): can MEA's refit C (four ion k_ij plus the MEAH⁺–water k_ij 1/T
slope, R4 fixed) run through `epcsaft.regression.fit`, with Ceres and exact
implicit derivatives, and why did MEA's 40 °C high-loading states take 20–60× the
median solve?

Inputs are read from the MEA checkout `calibration-misfit` at `e522dc4` (read
only): the state packet, the selected reactions, and MEA's production request
mapping `shared_evaluation._problem_from_request`. The parameter record is that
checkout's committed record (SHA-256 `868a5018…`), the pre-refit record refit C
started from; its start, bounds and scipy result come from
`calibration-misfit/refit-C-converged.json`. The Engine is the installed wheel.
Every result JSON records the wheel SHA-256, the MEA commit and the packet hash.

## Protocol

- **Timing** (`timing-{before,after}.csv`), all 123 packet states in process:
  - a cold solve at the pre-refit record (`base_cold`);
  - the refit-C start trial, seeded from that solve's anchor as MEA's `probe.py` does (`trial_anchor`);
  - a cold trial solve when the anchored one fails (`trial_protocol`, the sum of both);
  - the anchored trial through MEA's own forked child (`child_s`: fork, solve, snapshot with its action batch, JSON pipe).

  "Before" is MEA's pinned wheel `b66c7b96` (Engine 443a9da4). "After", and every other table here, is this branch's final wheel `684b213a`.
- **Derivatives** (`derivatives.csv`, `step-study.csv`).
  - Every refit-C coordinate against the `regression.differentiate` row.
  - Reference: centered re-solves of rebuilt models at steps h and h/2, combined by Richardson (h = 1e-4 for k_ij, 0.5 K for the slope).
  - Pass criterion: relative error ≤ 1e-6.
  - `separate` rebuilds the same row from single-direction `solved_state_actions` requests.
- **Jacobian cost** (`jacobian-{before,after}.csv`), per calibration state:
  - the pre-#152 regression path (compile, then one action request per parameter);
  - one first-action batch on the solve's compiled problem.

  "Before" is 34d06cd2: wheel `b66c7b96` has no parameter actions on referenced reactions.
- **Refit C** (`refit-c-residuals.{csv,json}`), exactly as `refit.py` scores it:
  - pCO2 rows: `partial_pressure`, `log_ratio`, weight 1/0.3²;
  - species rows: `mole_fraction`, `difference` with scale 0.1·obs + 0.001;
  - the NMR "HCO3⁻" row is the HCO3⁻ + CO3²⁻ pool;
  - every non-80 °C state the pre-refit record solves, declared from its anchor;
  - coordinates are scaled by refit C's sweep steps, inside its bounds, and capped at 6 iterations.

## Results

**Slow states: diagnosed and removed.**
- On `b66c7b96`, every state above 15 s was an anchored attempt that ran Ipopt to its 200-iteration limit. It was followed by a 1–3 s cold solve. This is visible in the 1419 MEA cache records and reproduced here.
- The anchored stall was already removed by #131's incipient-vapor seed.
- The remaining cold-start cost was Ipopt feasibility restoration: about 80 % of native samples, split between exact constraint Hessians and re-running the reference ladder at every trial pressure.
  - Its cause is the packet's unspeciated liquid guess, with molecular CO2 about 100× its speciated value.
  - That guess sends the fugacity pressure step to the wrong liquid's bubble pressure.
- The Engine now speciates the finite phases alone at the pressure seed before the boundary start (algorithms.md).
  - On 34d06cd2, vle_obs_0128, 0129, 0144 and 0187 failed cold, after 3.4, 21.5, 3.7 and 5.1 s. They now converge in 3 iterations.
  - On `b66c7b96`, vle_obs_0128, 0145 and 0227 failed cold.
- All 76 pCO2 states solved on both wheels agree to |ΔP|/P ≤ 2.9e-13, so they are on the same branch.

| Per-state solve (s), 123 states | median | p90 | max | failures |
|---|---|---|---|---|
| base cold, before | 0.35 | 2.09 | 7.7 | 6 |
| base cold, after | 0.26 | 0.36 | 2.2 | 3 |
| MEA trial protocol, before | 0.33 | 4.74 | 42.6 | 7 anchored, 1 overall |
| MEA trial protocol, after | 0.26 | 0.31 | 0.37 | 0 |
| MEA forked child, before | 0.45 | 4.43 | 41.7 | 7 |
| MEA forked child, after | 0.38 | 0.42 | 0.45 | 0 |

- Before, the 12 slowest states took 65 % of the trial wall time.
- MEA-side cost: the fork and snapshot path adds a median 0.11 s per state to the in-process solve. Its Ipopt banner prints once per child process.

**Jacobian.**

| Jacobian per state (s) | median | p90 | max |
|---|---|---|---|
| before: 4 k_ij, one request each (+ compile) | 0.338 | 0.353 | 0.406 |
| before: 4 k_ij, one batch | 0.225 | 0.234 | 0.292 |
| after: 4 k_ij, one batch | 0.038 | 0.040 | 0.041 |
| after: 5 coordinates (with the slope), one batch | 0.038 | 0.040 | 0.044 |

- Before, 99 % of the Jacobian time was the reference ladder. Each direction rebuilt it twice, with jets for all parameters.
- Now one ladder per state carries every direction. This is exact because its jets enter the rows only through Σ_a ∂s/∂θ_a (θ_a − θ_a,c).
- The batched and separate rows agree bitwise (0 relative) in all 40 derivative rows. This is center-reuse evidence only: separate requests use the same shared term. The re-solves below are what falsify the shared term.

**Derivatives.** All 40 rows pass ≤ 1e-6 relative.

| Observation / family | Rows | Worst relative error |
|---|---|---|
| Species mole fractions (Böttinger 050, 059) | 30 | 2.4e-8 |
| Partial pressure (vle_obs_0131, 0224) | 10 | 8.7e-7 |
| k_ij 1/T slope, off T_ref (059 at 333 K, 0224 at 393 K) | 4 | 4.8e-10 |

- The 8.7e-7 row is vle_obs_0224's weak carbamate–MEAH⁺ sensitivity. The step study shows it is re-solve noise: the reference converges onto the exact value as h grows (8.7e-7, 6.1e-8, 3.3e-8, 1.6e-8 at h = 1e-4 … 2e-3).
- vle_obs_0131 and Böttinger 050 lie at T_ref = 313.15 K. There dk_ij/db = 1/T − 1/T_ref = 0, so their slope rows are exactly zero on both sides and cannot falsify the slope.

**Refit C.** Refit C through `epcsaft.regression.fit`:
- converges, cost **121.35928** over the same 180 scaled residuals, against scipy's 121.35941;
- takes 6 iterations (5 steps, 6 residual and 6 Jacobian evaluations): 298 s fit plus 28 s setup, against about 2.5 h for the scipy forward-difference run;
- has no trial failures;
- ends with the bicarbonate–water k_ij on its upper bound, and carbamate–MEAH⁺ held 2e-8 inside its lower bound where scipy left it (both reported as active bounds; covariance withheld as `active_bound`);
- fitted values: k_ij 0.10691, −0.24168, 0.3, −0.3; slope −141.02 K.

Two earlier Engine runs, recorded here as diagnosis only, showed why bounds are an active set:
- Ceres' projected line search spent 45 Jacobian evaluations on 6 steps (1512 s).
- Per-coordinate clamping rejected five steps against the −0.3 bound it wanted to cross, stopping at 123.38.

## Limits

- Timing is single-threaded on a 2-CPU host with no MEA process running. It is equivalent-work evidence, not a solver certificate.
- Before/after timings use different wheels by design. Iteration counts are performance diagnostics only.
- Three liquid-only speciation states (Böttinger 042 and 046, Matin 016) exhaust Ipopt's 200-iteration budget cold at the pre-refit record on both wheels.
  - They converge in 16–24 iterations at refit C's start.
  - They lie outside both residual sets, so the comparison is unaffected.
  - Flooring the packet's zero H3O⁺/OH⁻ guesses does not change them.
  - #152 closed without resolving their liquid-speciation convergence from the packet guess.
- The osmotic / solvent-activity observation and cross-association l_ij^hb were designed in the issue, not built.
