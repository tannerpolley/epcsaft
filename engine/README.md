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
activity = reference.at_molality(1 * u.mole / u.kilogram)
gamma_pm = activity.mean_ionic_activity_coefficient_molality
osmotic_coefficient = activity.practical_osmotic_coefficient_molality
```

The practical molal osmotic coefficient is available for the single-water
reference path. Mixed-solvent references continue to return MIAC values and
report this scalar as `None` because no mixed-solvent osmotic standard state is
admitted.

For a native equilibrium consumer, `native_sdk(model)` returns the model-bound
capsule `epcsaft.native_sdk.v1`. Its accepted prefix evaluates a pure
nonassociating methane, ethane, or propane phase at explicit amount and volume,
returning the complete isothermal phase Helmholtz value on the documented
reference basis plus its gradient, Hessian, and third tensor. The
declaration-only C header is installed at
`epcsaft/include/epcsaft/native_sdk_v1.h`.

The same table has an additive pure-parameter entry for the first methane
regression slice. It returns the phase value, gradient, and Hessian in
`(n,V,m,sigma,epsilon/k)` coordinates so a regression consumer can assemble
exact residual Jacobians without copying EOS equations. Regression targets,
Ceres, bounds, diagnostics, and fitted records remain outside this package;
binary and associating parameter fits are not admitted by this entry.

A local append-only candidate entry evaluates the source-backed neutral
methane/ethane phase in ordered component amounts and volume. It returns the
extensive ideal-plus-hard-chain-plus-dispersion value, exact CppAD gradient and
row-major Hessian, and pressure through caller-owned buffers. This candidate is
not authoritative until its own promotion receipt is accepted; it adds no
phase identity, density solve, third derivative, equilibrium, or regression
contract. One further callback reuses the same result for coordinates
`(n1,n2,V,kij)`, treating the unlike interaction as an explicit finite trial
coordinate while leaving bounds, scaling, targets, and fit acceptance to the
regression owner. See `docs/science/neutral-mixture-native-sdk.md`.

The append-only electrolyte candidate callback evaluates the complete
ideal-plus-residual phase Helmholtz block in ordered component amounts and
volume, including association, Debye--Huckel, and corrected Born SSM+DS terms.
It admits only immutable source-backed model bundles. The distinct Khudaida et
al. (2026) bundle adds the exact water/ethanol/isobutanol/Na+/Cl- model at
293.15 K without changing the Figiel bundle. Phase bounds, complementarity,
equilibrium, and fit policy remain consumer-owned.

One further append-only callback converts caller-owned packing-fraction limits
to model-bound molar-volume limits at an ordered composition and temperature.
The Provider retains the temperature-dependent segment-diameter and packing
geometry; consumers retain their search limits and solver policy. The callback
adds no second result structure or fixed phase-volume constants.

An adjacent append-only callback evaluates the same model packing fraction in
explicit ordered component amounts and extensive volume. It returns the scalar
value and exact CppAD gradient/Hessian in `(n_1,...,n_C,V)` through caller-owned
buffers, so a consumer can impose composition-dependent packing constraints
without copying Provider geometry. It adds no equilibrium or density owner.

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
- the bounded Figiel, Yu, and Held (2025) water/methanol/ethanol
  chloride/bromide records needed for aqueous and mixed-solvent MIAC systems
  and explicit transfer reference states; plus
- a distinct Figiel aqueous iodide bundle for source-complete NaI and KI;
  LiI remains unsupported because the published Li+/I- interaction is blank;
  plus
- the distinct Khudaida et al. (2026) water/ethanol/isobutanol/Na+/Cl- records
  for the source-bound 293.15 K electrolyte-LLE challenge.

This repository does not perform phase equilibrium or phase discovery,
reaction/species equilibrium, regression, parameter fitting, ideal-gas
caloric calculations, or multi-salt reference paths. Solver
convergence is not a claim of global phase stability or literature validation.

Accepted promotion receipt `promotion-0012-provider-full-candidate` makes this
repository the sole production owner of the stated provider scope. The receipt
does not claim a release, equilibrium/global-stability capability, generic
electrolyte coverage, or native-SDK authority. The separate authority state for
the fixed SDK v1 prefix and its append-only pure-parameter tail is recorded only
by migration receipt `promotion-0015-provider-native-sdk-v1`; this
documentation does not transfer it, and consumer promotions remain separate.
