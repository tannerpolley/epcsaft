# Python API guide

This guide covers the supported public Python path:

```text
ParameterBundle -> ParameterSet -> EPCSAFT.evaluate -> EosResult
                                      |
                                      -> reference_state -> IonicActivity
```

It describes one model and one thermodynamic state at a time. Phase-equilibrium
discovery, global stability, reactions, regression, and parameter fitting are
separate capabilities and are not supplied by this package.

## Install and import

`epcsaft` requires Python 3.11 or newer. The project has not been published to a
package index. Install either a supplied wheel or a source checkout:

```console
python -m venv .venv
. .venv/bin/activate
python -m pip install /absolute/path/to/epcsaft-0.1.0.dev0-<platform>.whl
```

To build from a checkout, a C++17 compiler, CMake 3.20 or newer, and Ninja are
required. The build retrieves the pinned CppAD source declared by CMake.

```console
git clone https://github.com/ePC-SAFT/ePC-SAFT.git
cd ePC-SAFT
python -m venv .venv
. .venv/bin/activate
python -m pip install .
```

The public import is `epcsaft`; do not import `epcsaft._core` or other private
modules.

## Installed parameter bundles

There is deliberately no mutable catalog registry or list API. Use an exact ID
and version from this table, then inspect the immutable public `components`,
`domains`, and `sources` records on the loaded bundle.

| Bundle ID | Version | Exact component IDs | Source evidence domain and topology limits |
| --- | ---: | --- | --- |
| `gross-2001-methane-ethane` | 1 | `methane`, `ethane` | Neutral, nonassociating. Methane: 97–300 K; ethane: 90–305 K. A mixed selection uses the intersection. No pressure evidence bound is reported. |
| `gross-2001-propane` | 1 | `propane` | Pure, neutral, nonassociating; 85–523 K. No pressure evidence bound is reported. |
| `figiel-2025-reference-electrolytes` | 1 | `water`, `methanol`, `ethanol`, `lithium-cation`, `sodium-cation`, `potassium-cation`, `chloride-anion`, `bromide-anion` | Corrected Born SSM+DS at 298.15 K and 1 bar; total ion mole fraction at most 0.38. Reference states require one source-complete cation/anion formula unit. |
| `figiel-2025-aqueous-iodides` | 1 | `water`, `sodium-cation`, `potassium-cation`, `iodide-anion` | Corrected Born SSM+DS at 298.15 K and 1 bar; total ion mole fraction at most 0.38. Only the source-complete aqueous NaI and KI selections are supported; LiI is absent. |
| `khudaida-2026-figure-2-electrolyte-lle` | 1 | `water`, `ethanol`, `isobutanol`, `sodium-cation`, `chloride-anion` | Exact five-component corrected Born SSM+DS bundle at 293.15 K and 0.1 MPa. The Provider evaluates states but does not solve the reported LLE. |
| `ascani-2022-case-study-2-electrolyte-lle` | 1 | `water`, `1-butanol`, `sodium-cation`, `potassium-cation`, `chloride-anion` | Exact five-component original-Born bundle at 298.15 K and 1 bar. The Provider evaluates states but does not solve the reported multiphase problem. |

The table records scientific evidence bounds, not a promise of extrapolation.
Resolved temperature intersections and electrolyte ion-fraction caps are
checked during evaluation. Pressure bounds remain source metadata and are not
currently a runtime rejection check, so callers must keep pressure within the
reported evidence domain. A component appearing in a bundle does not imply
that every possible subset is parameter-complete; `EPCSAFT(...)` fails closed
when a selected topology lacks required source records.

`ParameterBundle.components` exposes the bundle's canonical component order.
`select(...)` accepts exact component IDs, not names or aliases, and its tuple
defines all later component and composition ordering. The bundle fingerprint
identifies its source records; the selected-set fingerprint additionally binds
the requested component order.

## Neutral state example

This example inspects the bundle, selects an explicit order, and evaluates the
same composition by direct density and by pressure closure.

```python
import math

from epcsaft import EPCSAFT, ParameterBundle, unit_registry as u

bundle = ParameterBundle.from_catalog("gross-2001-methane-ethane", version=1)
assert tuple(component.component_id for component in bundle.components) == (
    "ethane",
    "methane",
)

parameters = bundle.select(("methane", "ethane"))
model = EPCSAFT(parameters)
assert model.component_ids == ("methane", "ethane")
assert model.parameter_fingerprint == parameters.fingerprint

direct = model.evaluate(
    temperature=250 * u.kelvin,
    molar_density=8000 * u.mole / u.meter**3,
    mole_fractions=(0.4, 0.6),
)
assert direct.density_diagnostics is None
assert math.isfinite(direct.residual_helmholtz)
assert direct.log_fugacity_coefficient is not None

closed = model.evaluate(
    temperature=250 * u.kelvin,
    pressure=10 * u.bar,
    mole_fractions=(0.4, 0.6),
)
assert abs(closed.pressure.to("bar").magnitude - 10.0) < 1e-8
assert closed.density_diagnostics is not None
assert closed.density_diagnostics.branch == "single"
```

`residual_helmholtz` and each contribution are dimensionless molar residual
Helmholtz quantities on the `A_res/(nRT)` basis. `compressibility_factor` is
`Z`; `residual_chemical_potential_over_rt[i]` is `mu_i_res/(RT)`; and
`log_fugacity_coefficient[i]` is `ln(phi_i)`, all in `model.component_ids`
order. The Provider returns fugacity coefficients, not component fugacities.
Pressure, molar density, and residual caloric properties are Pint quantities.
See [the EOS science contract](science/full-provider-eos.md) for the equations,
units, and reference bases.

## Aqueous electrolyte example

Ionic mole fractions must be electroneutral. For a formula-unit molality path,
build the pressure reference once and reuse it across concentrations.

```python
import math

from epcsaft import EPCSAFT, ParameterBundle, unit_registry as u

bundle = ParameterBundle.from_catalog("figiel-2025-reference-electrolytes", version=1)
parameters = bundle.select(("water", "sodium-cation", "chloride-anion"))
electrolyte = EPCSAFT(parameters)

state = electrolyte.evaluate(
    temperature=298.15 * u.kelvin,
    pressure=1 * u.bar,
    mole_fractions=(0.90, 0.05, 0.05),
    phase="liquid",
)
assert state.density_diagnostics is not None
assert state.density_diagnostics.branch == "liquid"
assert state.debye_huckel != 0.0
assert electrolyte.permittivity_model == "figiel-composition"
assert state.bulk_relative_permittivity > 1.0

reference = electrolyte.reference_state(
    temperature=298.15 * u.kelvin,
    pressure=1 * u.bar,
    solvent_mole_fractions=(1.0, 0.0, 0.0),
    phase="liquid",
)
activity = reference.at_molality(1 * u.mole / u.kilogram)
assert activity.ion_component_ids == ("sodium-cation", "chloride-anion")
assert math.isfinite(activity.mean_ionic_activity_coefficient_molality)
assert activity.practical_osmotic_coefficient_molality is not None
```

`IonicActivity` reports formula-unit molality, the generated composition,
single-ion coefficients, stoichiometric mean coefficients on mole-fraction and
molality scales, and the practical molal osmotic coefficient for a
single-active-water reference. The osmotic field is `None` for mixed-solvent
references because no mixed-solvent osmotic standard state is admitted.

`EPCSAFT.permittivity_model` identifies the immutable, source-named topology
resolved from the parameter records. Current values are
`figiel-composition`, `ascani-mole-fraction`, `schreckenberg-density`,
`zuber-ion-specific`, and `rueben-perturbation` for ionic models. Callers do
not select this value during state evaluation. The latter three topologies are
available through complete source-bearing `ParameterBundle.from_records(...)`
inputs; they are not alternate modes silently applied to an installed catalog
bundle. Mixed or incomplete topology records fail during construction.

`EosResult.bulk_relative_permittivity` is the dimensionless state value used
by both Debye--Huckel and the bundle-selected Born formulation. Density-based
Schreckenberg and Rueben models approach one in the zero-density limit. The
Zuber implementation currently requires exactly one neutral solvent because
its source mixed-solvent equation requires solvent-volume data not represented
by the admitted parameter payload. Neutral models report
`permittivity_model="none"` and `bulk_relative_permittivity=None` rather than a
synthetic dielectric value.

## State and root contract

- `temperature`, `molar_density`, and `pressure` are Pint quantities.
  Dimensionless Python numbers are rejected for dimensioned inputs.
- `mole_fractions` and `solvent_mole_fractions` are finite, nonnegative,
  dimensionless numbers in the selected component order. Their sum must differ
  from one by no more than `1e-12`.
- Ionic states must satisfy `abs(sum(x_i * z_i)) <= 1e-12`. The Figiel bundles
  also require total ionic mole fraction at most 0.38 (with the same `1e-12`
  numerical allowance).
- Supply exactly one of `molar_density` or `pressure`. `phase` is invalid with
  direct density.
- A direct-density call evaluates the requested mathematical state and returns
  `density_diagnostics=None`; it does not certify mechanical stability or
  global phase stability.
- For pressure closure, `phase=None` succeeds only when one mechanically stable
  root is available. A multi-root state requires `phase="vapor"` or
  `phase="liquid"`. Conversely, a phase-specific request is rejected if that
  branch is unavailable; a phase label never relabels a lone root.
- Pressure results include the selected bracket, pressure residual, `dP/drho`,
  iterations, local mechanical-stability flag, and branch. These diagnostics
  certify the selected density root, not phase completeness or global
  stability.
- A pressure-specific `phase` request is a branch-identity assertion, not a
  liquid/vapor hint. If a density-dependent permittivity topology leaves only
  one stable root, omit `phase` and inspect `branch="single"`.

`log_fugacity_coefficient` is `None` when `Z` is nonfinite or nonpositive.
`residual_isobaric_heat_capacity` is `None` when the local pressure slope does
not support the constant-pressure transformation. Always check nullable fields
before using them.

## Stable errors

Parameter loading and selection raise `ParameterError`; state, density, and
reference evaluation raise `EosError`. Both derive from `ValueError` and expose
immutable `code` and `context` properties. Branch on those properties rather
than parsing the human-readable exception string.

```python
from epcsaft import EPCSAFT, EosError, ParameterBundle, ParameterError, unit_registry as u

try:
    ParameterBundle.from_catalog("not-installed", version=1)
except ParameterError as error:
    assert error.code == "missing-record"
    assert error.context == (("field", "catalog"), ("value", "not-installed/1"))

model = EPCSAFT(
    ParameterBundle.from_catalog("gross-2001-propane", version=1).select(("propane",))
)
try:
    model.evaluate(
        temperature=250,  # missing Pint units
        pressure=10 * u.bar,
        mole_fractions=(1.0,),
    )
except EosError as error:
    assert error.code == "invalid-state"
    assert error.context == (("field", "temperature"),)
```

Common `ParameterError.code` values include `missing-record`,
`invalid-selection`, `invalid-schema`, `invalid-unit`, and `invalid-domain`.
Common `EosError.code` values include `invalid-state`, `evaluation-failed`,
`invalid-reference-state`, and `reference-state-failed`. Domain violations use
`invalid-state` with a `field` such as `temperature`, `mole_fractions`, or
`electroneutrality`. Native SDK errors are a separate consumer contract and
are not part of normal Python state evaluation.

Solver convergence means only that the requested numerical state was found.
It does not establish numerical sensitivity, literature agreement, phase
equilibrium, or global stability. The native consumer boundary is documented
separately in [ADR 0001](adr/0001-native-consumer-transport.md).
