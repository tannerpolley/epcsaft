At 298.15 K, the complete Lithium #188 selection (water, Li⁺, Na⁺, K⁺,
Mg²⁺, Ca²⁺, Cl⁻, diluent, BuPhen, NaA, LiB₂A, MgA₂ and CaA₂) constructs
from the unchanged `/2` aqueous packet and read-only Lithium organic records.
The [input hashes](results/cation-mixing/dry-input-hashes.json) identify the
records and declaration sources; organic parameters remain source-owned priors.

The [native evaluations](results/cation-mixing/dry-support.csv) contain
14 states: pure diluent, five solutes at mole fraction $10^{-20}$ in
diluent, and a mixed neutral composition, each at fixed density 10 mol/m³ and
at the 101325 Pa liquid root. Available neutral chemical potentials are finite.
The mixing contribution and all retained term-on/term-off value, dense
Hessian and supported third-through-fifth energy-action differences are exactly
zero. Pure-diluent liquid density is 4327.6987500593159 mol/m³;
the largest absolute liquid pressure residual is
2.0981970010325313e-06 Pa.
The fixed-density pressure is an evaluated quantity, not a pressure closure.

With all organic amounts zero, the [aqueous comparison](results/cation-mixing/dry-aqueous-equality.csv)
at 100000 and 101325 Pa has exactly equal densities, pressures, mixing values,
common gradients, Hessians and chemical potentials. The
[positive-solvent qualification comparison](results/cation-mixing/dry-qualification-comparison.json)
compares 9,015 numeric fields, including the frozen calibration
predictions and exact parameter columns, held-out and rounding values, trace
checks and chloride transfer checks: 0 bit differences.

These are numerical verifications of a fixed-support definition. They do not
supply physical validation of organic parameters, ionic solubility, extraction,
or a full-composition dry Hessian. Transverse zero entries are defined storage.
Dry two-cation states and absent-ion insertion actions remain refused; the
missing partner's scalar chemical potential and fugacity are also refused on
a dry one-cation support. Independent delivered review and investigator
promotion remain pending; no new extraction prediction is made.
