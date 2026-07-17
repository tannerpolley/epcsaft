# ePC-SAFT Provider

`epcsaft` is a compact, source-traceable thermodynamic provider. One typed
parameter lifecycle feeds one native residual Helmholtz kernel covering
hard-chain, dispersion, association, Debye--Huckel, and corrected Born SSM+DS
contributions.

```python
from epcsaft import EPCSAFT, ParameterBundle, unit_registry as u

parameters = ParameterBundle.from_catalog(
    "gross-2001-methane-ethane", version=1
).select(("methane", "ethane"))

model = EPCSAFT(parameters)
explicit = model.evaluate(
    temperature=250 * u.kelvin,
    molar_density=8000 * u.mole / u.meter**3,
    mole_fractions=(0.4, 0.6),
)
pressure_specified = model.evaluate(
    temperature=250 * u.kelvin,
    pressure=10 * u.bar,
    mole_fractions=(0.4, 0.6),
)
```

For a bounded Figiel electrolyte, construct one explicit reference state and
reuse it across concentrations:

```python
electrolyte = EPCSAFT(
    ParameterBundle.from_catalog(
        "figiel-2025-reference-electrolytes", version=1
    ).select(("water", "sodium-cation", "chloride-anion"))
)
reference = electrolyte.reference_state(
    temperature=298.15 * u.kelvin,
    pressure=1 * u.bar,
    solvent_mole_fractions=(1.0, 0.0, 0.0),
    phase="liquid",
)
gamma_pm = reference.at_molality(
    1 * u.mole / u.kilogram
).mean_ionic_activity_coefficient_molality
```

For a native equilibrium consumer, `native_sdk(model)` returns the model-bound
capsule `epcsaft.native_sdk.v1`. Version 1 accepts only a pure nonassociating
methane, ethane, or propane model and evaluates one fixed-temperature phase at
explicit amount and volume. It returns the complete isothermal phase Helmholtz
value on the documented reference basis plus its gradient, Hessian, and third
tensor. The declaration-only C header is installed at
`epcsaft/include/epcsaft/native_sdk_v1.h`.

The capsule evaluates one phase state at a time. A consumer uses the same
immutable model for every candidate phase; the provider does not construct
liquid/vapor objects or own phase identity, Ipopt, equilibrium constraints, or
certification. See `docs/adr/0001-native-consumer-transport.md`.

`EosResult` contains contribution and total residual Helmholtz energies,
compressibility, pressure, molar density, component residual chemical
potentials, fugacity coefficients when defined, and residual entropy, internal
energy, enthalpy, `C_V`, and `C_P` when defined. Pressure-specified states also
return immutable root diagnostics with the bracket, pressure residual,
`dP/dρ`, iteration count, mechanical-stability status, and branch.

State quantities must be Pint quantities. Component and composition order is
the selected `ParameterSet` order. Supply exactly one of `molar_density` or
`pressure`. If more than one mechanically stable pressure root exists, select
`phase="vapor"` or `phase="liquid"`; otherwise the solve rejects the ambiguous
state. A phase-specific request rejects if only one stable root remains;
omit `phase` to evaluate that state with the explicit `branch="single"`
diagnostic.

The installed catalogs are:

- Gross and Sadowski (2001) methane/ethane neutral PC-SAFT parameters; and
- Gross and Sadowski (2001) pure-propane neutral PC-SAFT parameters needed by
  the first scoped saturation consumer; plus
- the bounded Figiel, Yu, and Held (2025) water/methanol/ethanol alkali-halide
  records needed for aqueous and mixed-solvent MIAC systems and explicit
  transfer reference states.

This repository does not perform phase equilibrium or phase discovery,
reaction/species equilibrium, regression, parameter fitting, ideal-gas
caloric calculations, or multi-salt reference paths. Solver
convergence is not a claim of global phase stability or literature validation.

Accepted promotion receipt `promotion-0012-provider-full-candidate` makes this
repository the sole production owner of the stated provider scope. The receipt
does not claim a release, equilibrium/global-stability capability, generic
electrolyte coverage, or promotion of the later local native-SDK candidate.
