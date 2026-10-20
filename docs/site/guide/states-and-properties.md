# States and properties

A state is one homogeneous phase at a specified temperature and composition.
Supply either molar density or pressure. At pressure, declare a density branch
when it matters. The [water notebook](../examples/parameters-and-water-state.ipynb)
checks pressure closure; the [methane caloric notebook](../examples/methane-caloric-properties.ipynb)
adds published ideal thermal records. See
[density closure](../theory/algorithms.md#alg-density) and
[calorics](../theory/algorithms.md#alg-calorics) for their definitions.

## Choose the independent quantities

`model.state(T, rho=..., x=...)` evaluates the T–density chart.
`model.state(T, P=..., x=..., phase="liquid")` solves density at pressure.
For `Pure`, omit `x`. For `Mixture`, provide fractions in `component_ids`
order or a mapping from component names. Omitted mapping entries are zero.
Fractions are dimensionless; temperature is K, pressure Pa and density mol/m³.
Physical electrolyte compositions must preserve charge balance.

A pressure solve may accept an available state as `anchor`. This is a numerical
start on the requested branch, not an instruction to accept an inconsistent
density. Check `density_convergence`, the achieved pressure and
`pressure_density_derivative`. A positive slope is a mechanical condition;
it is not a complete phase-stability test.

The [model reference](../api/eos.rst) owns these call signatures and the
[native state fields](../api/native/eos-state-reference.rst) describe the
returned values and causes.

## Distinguish Ideal, Residual and Total

| Returned quantity | Meaning and unit |
|---|---|
| `residual_helmholtz` and individual contributions | residual A/(nRT), dimensionless |
| `residual_enthalpy`, `residual_internal_energy` | departure from the ideal gas at the same T and density, J/mol |
| Entropy and heat capacities | J/(mol K) |
| `residual_chemical_potential_over_rt` | residual chemical potential divided by RT, component order |
| `log_fugacity_coefficient` | dimensionless ln φ at the state's T and pressure |
| `molar_density` | mol/m³; multiply by mixture molar mass in kg/mol for kg/m³ |

Total calorics add the ideal and residual quantities at the same state.
A residual enthalpy is not an absolute total enthalpy. Heat capacities in the
state record hold composition fixed; a reacting equilibrium heat capacity
also includes composition response and is a different observable.

## Supply ideal thermal records

The Engine reads polynomial, Shomate, NASA7, NASA9 and DIPPR107 ideal thermal
forms through declared intervals and reference constants. A
`ThermochemistryRecord` follows the EOS component order and supplies its
pressure standard. `ThermochemistryPacket` reads the supported NASA9 packet
and selects named species into that order. The native reader converts the
source gas constant to the Engine convention.

Keep source interval endpoints and anchors. Gaps and ambiguous overlaps are
refused; a piecewise endpoint with an available value need not have an
available derivative across the join. Do not extend or smooth a source
correlation to obtain a requested quantity. The methane notebook checks
published NASA9 values and the exact enthalpy temperature action against Cp;
it does not establish agreement with calorimetry.

## Interpret unavailable outputs

Without thermal inputs, residual properties can remain available while Ideal
and Total calorics are `None`. Inspect `thermochemistry_status` for
`MissingRecord`, `OutsideInterval` or another named cause. Do not substitute
zero or the residual part for an unavailable total property.

`derived_properties=True` requests additional fixed-composition properties
such as compressibility and expansivity. Their availability depends on
mechanical regularity and, for caloric consumers, the required thermal inputs.
A component chemical potential may be unavailable for ionic insertion on an
ion-free support even when the central state is available.

## Electrolyte reference quantities

For activities, construct a reference at the same T, P, parameters and solvent
basis. An ionic reference declares the ion pair and formula stoichiometry;
`Reference.activity` takes formula-unit molality in mol/kg. Solvent activity
and the practical osmotic coefficient use their declared solvent-mass basis.
The [NaCl notebook](../examples/aqueous-nacl-activity.ipynb),
[reference definitions](../theory/algorithms.md#alg-reference) and
[reference calls](../api/reference.rst) show these operations. They do not
supply an arbitrary single-ion convention or a global stability proof.
