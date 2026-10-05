# Nominally unloaded MEA liquid reference

These frozen declarations reproduce the trace-ion stationarity failure in
[Engine #261](https://github.com/tannerpolley/ePC-SAFT/issues/261): 7 mol MEA/kg
water at CO₂ loading 1e−5 mol/mol (a numerical proxy for zero), 353.15 K and
101325 Pa, one nine-species liquid. They qualify the Engine's Newton polish and
the equilibrium heat-capacity action at a trace-ion state. They do not adopt
these application parameters or validate calorimetry.

Both files were copied without modification from Amine-Thermodynamics commit
`c6e3f872e928c5858f756133de43b39c4992c0ad` (#124 assessment):
`parameters.json` from
`analyses/mea_parameter_bundle/results/density-cp-160/D1-evidence/F1-double-prime-candidate-parameters.json`
and `request.json` from
`analyses/mea_parameter_bundle/results/density-cp-160/cp-124-assessment/solution-request-0.0.json`.
Component order is CO₂, MEA, water, MEAH⁺, MEACOO⁻, HCO₃⁻, CO₃²⁻, H₃O⁺, OH⁻.

| Declaration | SHA-256 |
|---|---|
| `parameters.json` | `756fec502d3a1433d538abb073c4caeabf91902013ab7e27d5459ce70396543b` |
| `request.json` | `bc8a03e13efab9a9675a63121c12538dbb9f42a18dbe4ddc7fe6742b6ab6da61` |

The ideal-gas calorics are the identical record already retained as
[`../mea-zero-inventory/ideal-calorics.json`](../mea-zero-inventory/ideal-calorics.json)
(`c94d504e…`). The test
`test_nominally_unloaded_mea_resolves_trace_ion_stationarity` declares the
application's runtime choices explicitly: Born shell and dielectric scales of
one, the request's five reaction laws and references, its feed and composition
guess, and the neutral reference (water solvent; neutral solutes 1 and summed
source ion pairs, scaled by 1/8).
