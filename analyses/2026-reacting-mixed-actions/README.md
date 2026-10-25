# Referenced mixed-temperature second actions and pressure-moving enthalpy (#147)

**Frozen record — not rerunnable as a current-main calculation.** Record snapshot:
commit `ca719012c31e8b41edcb98ae82a388117b58d343`. Original wheel identities, source
inputs and retained numerical evidence keep their recorded provenance. The retained
mixed-action comparison requires the original MEA and absorber application loaders and
wheel; those frozen loaders are excluded from current declaration checks.

Question: on reacting liquids with source-referenced reactions, do the new exact
actions equal re-solves? Families: `D²O[e_T, e]` for e = CO2 feed, P and T
(rows carry s_TT and s_TP, EqIDs `standard_state_transfer_temperature_second_derivative`
and `standard_state_transfer_temperature_pressure_derivative`) and the extensive liquid
enthalpy action `d(n_L h_L)` along P (fixed-P states) or along T and CO2 feed with the
solved pressure moving (bubble states), which uses `∂Δ_r h_ρ/∂P = RT² ∂²C/∂T∂P`.

Inputs: MEA-Thermodynamics `04a9328` (production request mapping, adopted parameters,
physical ideal-gas records of `verify_reference_calorics.py`) for Böttinger 2008 states
050 and 059 (liquid, fixed T/P) and Hilliard `vle_obs_0131`, `vle_obs_0224` (bubble,
solved P); absorber `4f2d12e9` exploratory record `868a5018…` with its reactions and
ideal-gas records at the case 3C bottom-node and S2-node liquid inputs of
`analyses/greenfield_node_qualification`. Both checkouts are read-only inputs.

## Protocol

- Observables: liquid ln a_i of all nine species and n_L; bubble states add P and y_CO2.
- Reference: the first action `D O[e]` (or re-solved `n_L h_L`) differenced along the
  second direction over warm-started re-solves, centred steps h and h/2 (T 0.2/0.1 K,
  P 2e4/1e4 Pa, CO2 2e-4/1e-4 mol), Richardson R = (4D(h/2) − D(h))/3. `vle_obs_0224`
  is at R2's correlation bound (393.15 K), so its T differences use the second-order
  backward formula of the same h² order. `resolution` = |D(h/2) − D(h)|/3 is the h²
  term Richardson removes.
- Pass: |a − R| ≤ 1e-6 |R|.

## Result (`results/actions.csv`, wheel `48a639e7…`, 25 s)

176/176 rows pass. Worst relative error per family:

| Family | MEA states | Absorber nodes |
|---|---|---|
| D²O[e_T, e_CO2] (N4 e_T) | 9.9e-8 (0224, n_L) | 1.5e-9 (3C, ln a_CO2) |
| D²O[e_T, e_P] | 1.7e-9 (050, ln a_H3O+) | 2.9e-7 (3C, ln a_CO3²⁻, value 1.7e-14/Pa) |
| D²O[e_T, e_T] | 8.6e-8 (0224, n_L) | 1.1e-10 (3C, n_L) |
| d(n_L h_L)/dP | 3.0e-10 (059) | 6.2e-9 (S2) |
| d(n_L h_L)/dT, solved P | 2.1e-11 (0224) | — |
| d(n_L h_L)/dCO2, solved P | 9.5e-10 (0224) | — |

In 80 rows the h² term exceeds 1e-6|R|; Richardson removes it and the action still
agrees to ≤ 2.9e-7, so those rows show the O(h²) approach, not an unresolved reference.

## Claim limits

Numerical verification on the six listed states only, not physical validation. The
absorber record is exploratory (not adopted). Heat-capacity T/P actions, caloric
second order and parameter second order remain typed refusals (not tested here).
