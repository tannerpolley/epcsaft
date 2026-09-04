# Stage 05 campaign — the first measured comparisons through the new solver

Branch `stage-05-campaign`, based on the greenfield commit that merges the harness (the
coordinator pins). Machine-readable claims are in `briefs/stage-05-campaign.brief.json`.
This is the owner's step 2: before any regression work, the accepted equilibrium is run
against measured phase-equilibrium data, every number retained.

**Depends on the harness landing first**: every cell writes a results file in the harness's
schema, and acceptance reads it. The one write that does not depend on it is the packet
transcription, which may start now.

## The rule

Each oracle cell is a Python test under `engine/tests/oracle/` in the pattern of
`tests/eos/test_activity.py`: rows from `engine/tests/oracle/oracle.json`, a model from a
packet by `Parameters.from_json`, a `Problem` through `epcsaft.equilibrium`, one
`solve_equilibrium` per row. Python arithmetic is allowed in tests, never in
`src/epcsaft/`. Each row asserts two things and reports the rest: the solve converged
(`stop` is `converged` or `converged_to_floor`) and the record's `residual ≤ max(criterion,
floor)`. Deviations from the measured values are computed, compared to the cell's stated
bar, written to the results file with `within`, and **never asserted**: a miss is a finding
about the model or the data, recorded with its trace; a solver stop on a measured row is
the defect this campaign exists to find and is a failed row by name.

Numbers are retained under `analyses/reproduction/2026-09-04-measured-campaign/`: the
results files, the commit and wheel hash, and the table `tools/summarize-check-results.py
--table` prints, committed by the coordinator at acceptance.

## The five cells

| cell | rows | fixture and parameters | problem | reported against the bar |
|---|---|---|---|---|
| ethane saturation | `light_alkane_saturation`, 10, NIST SRD 69, 100–280 K | `gross-sadowski-2001` ethane | `C = 1`, `{L, V}` with `V` pinned at zero, `T` fixed, `P` free | `P` and liquid density; bar: Gross and Sadowski 2001's reported deviations for the light alkanes, 2 % in `P`, 2 % in density, named as the paper's, not ours |
| methane–ethane bubble points | `light_alkane_binary_vle`, 17, May et al. 2015 (10.1021/acs.jced.5b00610), 243–271 K, **new case, transcribed from the superseded tree's `analyses/reproduction/2015-may/inputs/source/may-2015-methane-ethane-vle.csv` into a greenfield packet with the CSV's provenance and uncertainties** | `gross-2001-methane-ethane` (`k_ij = 0` by Gross's statement) | `{L, V}` with `V` pinned, `T` and liquid `x` fixed, `P` free | `P` against `pressure_standard_uncertainty_Pa` and `y` against `y_comparison_allowance`, both from the rows; bar 3 standard uncertainties |
| MEA–water isobaric | `aqueous_amine_binary_vle`, 6, Cai 1996, 101.33 kPa | `associating_binary()`'s parameters as a packet row set (Baygi 2015 with the fixture's `k_ij = −0.052`), stated honestly; the MEA repository's fitted `−0.0735` is recorded beside it, not used | `{L, V}` with `V` pinned, `P` and `x` fixed, `T` free | `T` against `temperature_resolution_k` and `y` against `composition_uncertainty_mole_fraction`; bar 1 K and 0.02 in `y`, the published model's own deviation domain |
| Khudaida tie lines | `aqueous_alcohol_salt_lle`, 8, 293.15 K, 5 wt % NaCl | `khudaida-2026-figure-2-electrolyte-lle` revision 6 | `{L, L}` declared, `T`, `P` fixed, feed the tie-line midpoint, composition guesses the two measured endpoints | endpoint compositions; bar the paper's Table 9 AAD per component and phase (grand 0.016). **Gated**: today `Parameters.from_json` on the packet fails with `missing association input association/isobutanol/a/isobutanol/b/energy_over_k` (*measured*, 2026-09-04), so the cell is `skipped` by that name until the packet or the builder is corrected; the skip is the finding |
| Jou 1995 partial pressures | `aqueous_amine_carbon_dioxide_reactive_vle`, 10, 313–393 K, 30 wt % MEA | `mea-co2-h2o-nine-species-estimation` revision 3, which builds a nine-species model today (*measured*); its `reaction_correlations` carry R4 and R5 only, so **R1–R3's `ln K = a + b/T + c ln T + dT` coefficients are transcribed into revision 4 from the MEA repository's source contract (Austgen et al. 1989)** as a packet write | `{L, V}` with `V` restricted to `{CO2, MEA, water}` and pinned, `T` and loading fixed, `P` free; `K(T)` per reaction evaluated in the test from the packet's declared kinds at the row's `T` and passed as `standard_gibbs_j_per_mol` through the registered recipe | `pCO2 = y_CO2 · P` against `co2_partial_pressure_pa`; bar the MEA repository's current candidate, 0.24 in `log10`, recorded as theirs |

Each cell also records the two solver findings it can classify:

- **The binary `{V,L}` flash at 180 K and 1 bar with nothing active** (repair, tenth resume):
  reproduce once, retain the record. By the EOS surface's own numbers the feed is above its
  dew pressure (1.58 bar), so a `declared_phase_absent` refusal is a misplaced fixture state
  and closes the finding; a budget stop is "a vanished phase on the superheated side stalls
  to the budget", recorded for a fix brief with the trace.
- **The reactive `{L,L}` cell with `L₂` pinned at 1 % and `T` free** (repair, resume 3):
  certify the state first with fixed-`T`, `P` flashes at 293, 298 and 303 K and read whether
  `L₂`'s amount crosses 1 % in that window; if it does, the free-`T` stall is a driver
  finding for its own brief; if not, the state has no solution and the fixture row is
  corrected.

## Traps

- No tolerance is added to `numerical_policy.hpp`; bars live in the test with their source.
- No fixture is added to `model_support.hpp`; the cells build models from packets, which is
  what a user does. The MEA–water cell states its parameter provenance beside its numbers.
- The trivial solution is refused, never counted: a row whose solve returns
  `phases_coincide` is a failed row by name, not a near-miss.
- Read the tree through Serena before writing: `find_symbol` for `solve_equilibrium`,
  `problem_from_mapping`, `Parameters.from_json`; `get_symbols_overview` on
  `src/epcsaft/equilibrium/problem.py` and `tests/eos/test_activity.py`.

## Checks, named

| id | asserts | mutation (must go red) |
|---|---|---|
| oracle-rows-converge | every row of every runnable cell converges with `residual ≤ max(criterion, floor)`; skipped cells are skipped by name | `EQUILIBRIUM_MAX_ITERATIONS` read as one |
| oracle-results-retained | the campaign's results file validates and carries one cell per row with `deviations` and `bar` | the writer omits `bar` |
| packet-provenance | the two new packet revisions validate and every transcribed row names its source locator | one locator blanked |

## Stopping

`budget_ceiling`, `gate_impossible`, `source_contradicted`, `environment_broken`,
`mutation_invalid`. A cell that cannot build its model is `skipped` by name, not a stop.
