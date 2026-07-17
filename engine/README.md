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
state. A phase hint on a one-root state reports `branch="single"` rather than
inventing a phase identity.

The installed catalogs are:

- Gross and Sadowski (2001) methane/ethane neutral PC-SAFT parameters; and
- the bounded Figiel, Yu, and Held (2025) water/methanol/ethanol alkali-halide
  records needed for aqueous and mixed-solvent MIAC systems and explicit
  transfer reference states.

This repository does not perform phase equilibrium or phase discovery,
reaction/species equilibrium, regression, parameter fitting, ideal-gas
caloric calculations, or multi-salt reference paths. Solver
convergence is not a claim of global phase stability or literature validation.

The full provider implementation is currently a local candidate. The accepted
promotion receipt still covers only the earlier neutral explicit-density
Slice 1 baseline; publication and promotion require separate validation and
approval.
