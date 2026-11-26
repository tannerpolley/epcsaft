# A7 MEA parameter reference

This frozen parameter document is the MEA–water–CO₂ set, labeled A7, that
Amine-Thermodynamics selects for its current results. It uses the
`solvent-density-temperature` permittivity rule. The Engine keeps it for
below-bubble MEA liquid–vapor cases. At 353.15 K and 41 kPa, with feed
0.358 mol CO₂, 1 mol MEA and 7.93 mol water, the liquid water mole fraction is
0.8416. The copy qualifies Engine behavior only. It does not adopt these
application parameters.

`parameters.json` was copied without modification from Amine-Thermodynamics
commit `046e7462f67c5d36b0a8d96c2d70ffcd4acce91e`:
`analyses/mea_parameter_bundle/results/selected-current-best-parameters.json`
(SHA-256 `ae17803e1bcb2ca0a1970e05c7cc235d615c4d8a3a349b5afef0b0401dedefed`).
Component order is CO₂, MEA, water, MEAH⁺, MEACOO⁻, HCO₃⁻, CO₃²⁻, H₃O⁺, OH⁻.

The reaction stoichiometry, the R1 and R3 laws and the neutral reference come
from the existing MEA declaration in `tests/native_inputs.py`. The document's
`reaction_correlations` replace the R2, R4 and R5 coefficients using Amine's
`shared_evaluation.py::_engine_reaction_records` conversions.
