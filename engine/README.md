# ePC-SAFT Provider

This package provides source-traceable parameters and a neutral,
nonassociating, nonionic PC-SAFT evaluator at explicit molar density. Its root
surface is exactly `ParameterBundle`, `ParameterSet`, `ParameterError`,
`unit_registry`, `EPCSAFT`, and `EosResult`.

```python
from epcsaft import EPCSAFT, ParameterBundle, unit_registry

parameters = ParameterBundle.from_catalog(
    "gross-2001-methane-ethane", version=1
).select(("methane", "ethane"))

result = EPCSAFT(parameters).evaluate(
    temperature=250 * unit_registry.kelvin,
    molar_density=8000 * unit_registry.mole / unit_registry.meter**3,
    mole_fractions=(0.4, 0.6),
)
```

`EosResult` contains the dimensionless `hard_chain`, `dispersion`,
`residual_helmholtz`, and `compressibility_factor` values, plus `pressure` as
a package-registry Pint quantity in pascals. Temperature and molar density must
be Pint quantities; component and composition order follows the selected
`ParameterSet` exactly.

The installed catalog contains only the Gross and Sadowski 2001 methane/ethane
parameter bundle. The compact Figiel/Yu/Held 2025 and Ascani/Held 2021 data in
`tests/data` is a disconnected parser fixture: it is not installed, cataloged,
or a physically complete mixture.

This boundary does not provide specified-pressure density closure,
association, ions, fugacity or caloric properties, phase equilibrium,
regression, a native SDK, or a solver. Specified-pressure density closure over
this pressure evaluator is the next functional slice.
