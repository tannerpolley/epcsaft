# ePC-SAFT Provider

`epcsaft` is a compact, source-traceable thermodynamic provider. One typed
parameter lifecycle feeds one native residual Helmholtz kernel covering
hard-chain, dispersion, association, Debye--Huckel, and corrected Born SSM+DS
contributions.

## Installation and Python use

`epcsaft` requires Python 3.11 or newer and is not published to a package
index. Install a supplied wheel with `python -m pip install /path/to/wheel`, or
build the current checkout with `python -m pip install .`; source builds require
a C++17 compiler, CMake 3.20 or newer, and Ninja.

The [Python API guide](docs/python-api.md) is the canonical walkthrough for
installed bundle IDs, component discovery and ordering, Pint state inputs,
density and pressure routes, root selection, result interpretation,
electrolyte activities, and stable errors.

The [checked equation mirror](docs/equations.md) maps every retained
mathematical relation to its C++ or Python implementation owner and to the
[normative scientific contract](docs/science/full-provider-eos.md). The
implementation remains executable authority; the generated Markdown and
LaTeX views are synchronized documentation mirrors.

`ParameterBundle.catalog()` returns the installed bundles as an immutable tuple.
You can inspect exact IDs, versions, components, sources, and domains without a
mutable registry.

The same guide includes a complete [user-provided bundle
workflow](docs/python-api.md#user-provided-parameter-bundles). Create typed
caller records with `ParameterBundle.from_records(...)`, persist the bundle
with `bundle.to_path(...)`, and reload it with
`ParameterBundle.from_path(...)`. Keep source provenance, validity domains,
units, and required pair parameters explicit.

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
assert pressure_specified.component_ids == model.component_ids
assert pressure_specified.mole_fractions == (0.4, 0.6)
```

## Direct C++ debug harness

`export_native_model(model, path)` writes the Provider's resolved immutable
model format for native consumers. Installed wheels export the CMake target
`epcsaft::native_sdk` from `epcsaft/cmake`; its public C headers are under
`epcsaft/include`. A consumer loads the resolved model with
`epcsaft_native_model_load_v1`, keeps that handle alive while it uses the
returned `epcsaft_native_sdk_v1` table, and releases it with
`epcsaft_native_model_destroy_v1`. This is the supported installed transport;
it exposes no Provider C++ implementation types or EOS equations.
`epcsaft_native_model_fingerprint_v1` returns the handle-owned model identity
for downstream result and callback validation.

For local Provider debugging, configure `-DEPCSAFT_ENABLE_CPP_DEBUG=ON`, build
`epcsaft_cpp_debug`, and run its `epcsaft_cpp_debug_smoke` CTest. The executable
uses that same installed-consumer model and SDK contract. The debug executable
itself remains opt-in and is not installed.

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

For an in-process Python extension consumer, `native_sdk(model)` returns the
same model-bound table in capsule `epcsaft.native_sdk.v1`. Its accepted prefix evaluates a pure
nonassociating methane, ethane, or propane phase at explicit amount and volume,
returning the complete isothermal phase Helmholtz value on the documented
reference basis plus its gradient, Hessian, and third tensor. The
declaration-only C header is installed at
`epcsaft/include/epcsaft/native_sdk_v1.h`.

The append-only neutral-reference candidate returns a complete charge-neutral
stoichiometric basis and its dimensionless reference log-fugacity contractions
for an admitted aqueous phase model containing one explicit water reference
carrier, any neutral components, and one or more cations and anions. Neutral
unit rows use their infinite-dilution fugacity in the reported pure-water
composition; minimal cation/anion formula-unit rows use the canonical
four-molality limits. Component order, charges, source fingerprint,
composition, solvent molar mass, source-domain state, convergence, and the
Helmholtz basis are explicit. Authenticated user bundles with unknown
applicability expose the topology but fail value evaluation until finite source
temperature and pressure bounds are supplied. The callback returns no
derivatives or absolute single-ion standards; reaction constants, source
standard-state transformation, and equilibrium remain consumer-owned.

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
volume, including association, Debye--Huckel, and the bundle-selected Born
formulation. It admits immutable source-backed catalog models and exact
caller-owned models selected from persisted `user-provided` bundles. The latter
admission supplies generic phase, volume-domain, and packing transport only; it
does not validate the caller's parameter provenance or admit source-specific
reference and parameter-derivative callbacks. The distinct
Khudaida et al. (2026) bundle adds its corrected SSM+DS
water/ethanol/isobutanol/Na+/Cl- model at 293.15 K. The distinct Ascani et al.
(2022) bundle adds the original-Born water/1-butanol/Na+/K+/Cl- model at
298.15 K, with its combined dielectric rule and explicit `l_ij` records.
Neither mutates a Figiel bundle. Phase bounds, complementarity, equilibrium,
and fit policy remain consumer-owned.

A private, non-authoritative Held--Cameretti--Sadowski candidate adds the exact
ordered `(water, hydronium-cation, hydroxide-anion)` mechanical phase model at
`298.15 K` and `1 bar`. Its formulation uses the source-bound water and ion
parameters, water--ion dispersion, no ion--ion dispersion, association, and
Debye--Huckel contribution; Born is absent. The existing generic
charge-neutral reference tail returns its gauge-invariant water and
cation-plus-anion basis contractions, convergence evidence, fingerprint, and
Helmholtz-basis identity. The ordinary public pure-water state supplies the
fixed-pressure liquid molar volume. It does not own `K_w`, reaction
stoichiometry, reactive standard-state conversion, species equilibrium, or any
solver, and has `authority_effect: none`.

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
Each result also carries its canonical temperature, ordered composition,
component IDs, and parameter fingerprint. Distributed wheels include `py.typed`;
the public interface provides Pint quantity, phase, bundle-purpose, domain, and
parameter-family annotations for static tooling.
The model-neutral `born` field carries the Born contribution selected by the
immutable parameter bundle, whether corrected SSM+DS or original Born.
`EPCSAFT.permittivity_model` identifies the source-named dielectric topology,
and `EosResult.bulk_relative_permittivity` exposes the exact value shared by
Debye--Huckel and Born. In addition to the installed Figiel and Ascani
topologies, complete source-bearing record sets can resolve Schreckenberg
density, Zuber ion-specific, or Rueben perturbation topologies. These are
parameter-record capabilities, not caller-selected phase modes or new catalog
defaults; unsupported or mixed records fail closed.

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
  at 293.15 K; and
- the distinct Ascani et al. (2022) water/1-butanol/Na+/K+/Cl- original-Born
  records at 298.15 K; and
- the private Held et al. (2008) water/H3O+/OH- candidate records at
  298.15 K and 1 bar, with no authority or reactive-equilibrium claim.

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
