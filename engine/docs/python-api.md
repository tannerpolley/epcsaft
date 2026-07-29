# Python API guide

This guide documents the sole public Python frontend for version
`0.2.0.dev0`:

```text
Parameters -> Mixture -> State
                         |
                         -> State.fugacity
                         -> State.activity (reference-created electrolyte states)
```

The frontend evaluates one resolved model and one thermodynamic state at a
time. Phase equilibrium, phase discovery, global stability, reactions,
regression, and parameter fitting are separate capabilities and are not
supplied here. There are no compatibility aliases for an earlier frontend.

## Install and import

`epcsaft` requires Python 3.11 or newer and is not published to a package
index. Install a supplied wheel or build the checkout. Source builds require a
C++17 compiler, CMake 3.20 or newer, and Ninja. Pint is the sole Python runtime
dependency.

```console
python -m venv .venv
. .venv/bin/activate
python -m pip install /absolute/path/to/epcsaft-0.2.0.dev0-<platform>.whl
```

The supported import is from `epcsaft`. Private modules are not part of this
contract.

## Parameter resolution

`Parameters` is immutable. Its selected `component_ids` tuple establishes the
order for every composition, property, fugacity, and activity result.
Parameter values remain tied to bundle fingerprints, source records, and
validity domains; no scientific value is inferred or silently defaulted.

The four parameter entry points are:

| Entry point | Purpose |
| --- | --- |
| `Parameters.from_catalog(bundle_id, components=..., version=None)` | Load one installed catalog bundle and select an explicit component order. Omitted `version` selects the greatest installed integer version for the exact ID. |
| `Parameters.from_bundle(path, components=...)` | Load one strict version-1 bundle directory and select an explicit component order. |
| `Parameters.from_records(...)` | Build one immutable caller-owned bundle from typed records. Use `purpose="user-provided"` for caller-owned values. |
| `Parameters.create_bundle_template(path, components=..., profile=..., correlations=False)` | Create a non-overwriting authoring skeleton for a `neutral`, `associating`, `electrolyte`, or `associating-electrolyte` profile. |

### Installed catalog

The catalog is immutable package data, not a mutable registry. The following
table is the installed source-backed catalog for this frontend. Its temperature,
pressure, topology, and ion-fraction entries are evidence domains, not promises
of extrapolation.

| Bundle ID | Version | Exact component IDs | Source evidence domain and topology limits |
| --- | ---: | --- | --- |
| `gross-2001-methane-ethane` | 1 | `methane`, `ethane` | Neutral, nonassociating. Methane: 97–300 K; ethane: 90–305 K. A mixed selection uses the intersection. No pressure evidence bound is reported. |
| `gross-2001-propane` | 1 | `propane` | Pure, neutral, nonassociating; 85–523 K. No pressure evidence bound is reported. |
| `figiel-2025-reference-electrolytes` | 1 | `water`, `methanol`, `ethanol`, `lithium-cation`, `sodium-cation`, `potassium-cation`, `chloride-anion`, `bromide-anion` | Corrected Born SSM+DS at 298.15 K and 1 bar; total ion mole fraction at most 0.38. Reference states require one source-complete cation/anion formula unit. |
| `figiel-2025-aqueous-iodides` | 1 | `water`, `sodium-cation`, `potassium-cation`, `iodide-anion` | Corrected Born SSM+DS at 298.15 K and 1 bar; total ion mole fraction at most 0.38. Only the source-complete aqueous NaI and KI selections are supported; LiI is absent. |
| `held-2008-water-self-ionization` | 1 | `water`, `hydronium-cation`, `hydroxide-anion` | Debye–Hückel electrolyte model without a Born term at 298.15 K and 1 bar. The Provider supplies the mechanical state and charge-neutral reference basis, not reaction equilibrium or `K_w`. |
| `khudaida-2026-figure-2-electrolyte-lle` | 1 | `water`, `ethanol`, `isobutanol`, `sodium-cation`, `chloride-anion` | Exact five-component corrected Born SSM+DS bundle at 293.15 K and 0.1 MPa. The Provider evaluates states but does not solve the reported LLE. |
| `ascani-2022-case-study-2-electrolyte-lle` | 1 | `water`, `1-butanol`, `sodium-cation`, `potassium-cation`, `chloride-anion` | Exact five-component original-Born bundle at 298.15 K and 1 bar. The Provider evaluates states but does not solve the reported multiphase problem. |

Load an exact installed selection as follows:

```python
from epcsaft import Mixture, Parameters, unit_registry as u

parameters = Parameters.from_catalog(
    "gross-2001-methane-ethane",
    components=("methane", "ethane"),
    version=1,
)
mixture = Mixture(parameters)
assert mixture.component_ids == ("methane", "ethane")
```

The `components` argument is required for catalog and bundle loading. IDs are
exact, lowercase identifiers; component names and aliases are not alternate
selection routes. A selected subset is accepted only when its complete source
records define a valid topology. Missing unlike-pair, site, charge,
permittivity, Born, or association records fail closed during construction.

### Caller-owned bundles

`Parameters.from_records(...)` accepts typed records from `epcsaft.records`.
Every scientific record must point to a source and a validity domain present in
the same bundle. Use neutral, associating, or electrolyte records according to
the complete selected topology; do not rely on a hidden physical default.

This compact pure-component example shows the required provenance shape without
using a catalog value:

```python
from decimal import Decimal

from epcsaft import Parameters, unit_registry as u
from epcsaft.records import (
    ComponentRecord,
    SingleParameterRecord,
    SourceRecord,
    ValidityDomain,
)

source_id = "source-1"
domain_id = "domain-1"
provenance = {
    "source_id": source_id,
    "locator": "table-1",
    "domain_id": domain_id,
}

parameters = Parameters.from_records(
    bundle_id="user-neutral",
    bundle_version=1,
    purpose="user-provided",
    sources=(SourceRecord(source_id, "source citation", "model input"),),
    domains=(
        ValidityDomain(
            domain_id,
            "fit-range",
            temperature_min=250 * u.kelvin,
            temperature_max=350 * u.kelvin,
        ),
    ),
    components=(ComponentRecord("sample"),),
    singles=(
        SingleParameterRecord(
            "sample-segment-count", "sample", "segment_count", Decimal("1.5"), **provenance
        ),
        SingleParameterRecord(
            "sample-segment-diameter",
            "sample",
            "segment_diameter",
            Decimal("3.5") * u.angstrom,
            **provenance,
        ),
        SingleParameterRecord(
            "sample-dispersion-energy",
            "sample",
            "dispersion_energy_over_k",
            Decimal("200") * u.kelvin,
            **provenance,
        ),
    ),
)
```

`Parameters.create_bundle_template(...)` creates the strict directory shape
without overwriting an existing path:

```python
from epcsaft import Parameters

Parameters.create_bundle_template(
    "user-bundle",
    components=("sample",),
    profile="neutral",
)
# Populate the generated files with source-complete records, then load them:
parameters = Parameters.from_bundle("user-bundle", components=("sample",))
```

`from_bundle` validates the manifest, exact file inventory, record types, units,
provenance references, topology completeness, and canonical fingerprint. A
caller-owned bundle remains caller-owned; loading it does not mutate the
installed catalog or turn its values into a validated literature model.

## Mixture and state evaluation

Construct a mixture only from resolved parameters:

```python
mixture = Mixture(parameters)
```

The state entry point has two routes:

```python
direct = mixture.state(
    T=250 * u.kelvin,
    rho=8_000 * u.mole / u.meter**3,
    x=(0.4, 0.6),
)

closed = mixture.state(
    T=250 * u.kelvin,
    P=10 * u.bar,
    x=(0.4, 0.6),
    phase="vapor",
)
```

`T`, `P`, and `rho` must be Pint quantities with the corresponding dimensions.
Supply exactly one of `P` and `rho`; `phase` is invalid on the direct-density
route and is optional on pressure closure. `x` is a finite, nonnegative,
normalized sequence in `mixture.component_ids` order. The normalization
tolerance is `1e-12`.

For ionic models, `x` must also satisfy electroneutrality within `1e-12` and
the selected bundle's ion-fraction limit. The resolved model's generic bundle-defined
relative-permittivity formulation is a parameter capability, not a caller
selected mode.

An explicit density call evaluates the requested mathematical state. It does
not certify mechanical stability or global phase stability. A pressure call
solves density closure and returns diagnostics containing the selected bracket,
pressure residual, `dP/drho`, iteration count, local mechanical-stability flag,
and branch. With `phase=None`, pressure closure succeeds only when one
mechanically stable root is available. A genuine multi-root state requires
`phase="vapor"` or `phase="liquid"`; a phase label never relabels a lone root.

Every `State` retains the canonical temperature, ordered composition,
component IDs, parameter fingerprint, pressure, molar density, compressibility
factor, contribution values, and nullable diagnostics. The scientific equations,
reference bases, and derivative authority are defined in
[the full-provider EOS contract](science/full-provider-eos.md).

## Residual properties and fugacity

The residual methods are immutable projections of one evaluated state:

| Method | Default `basis="molar"` | `basis="reduced"` |
| --- | --- | --- |
| `State.ares()` | residual molar Helmholtz energy, J/mol | `Ares/(nRT)` |
| `State.ures()` | residual molar internal energy, J/mol | `Ures/(nRT)` |
| `State.hres()` | residual molar enthalpy, J/mol | `Hres/(nRT)` |
| `State.gres()` | residual molar Gibbs energy, J/mol | `Gres/(nRT)` |
| `State.sres()` | residual molar entropy, J/(mol K) | `Sres/(nR)` |
| `State.cvres()` | residual molar `C_V`, J/(mol K) | `C_V^res/R` |
| `State.cpres()` | residual molar `C_P`, J/(mol K), or `None` | `C_P^res/R`, or `None` |

The reduced basis returns plain floats. The molar basis returns Pint quantities.
`cpres()` is nullable when the local pressure slope does not support the
constant-pressure transformation.

```python
ares = closed.ares()
reduced_ares = closed.ares(basis="reduced")
entropy = closed.sres()
assert ares.to("joule / mole").magnitude == reduced_ares * 8.31446261815324 * 250
assert entropy.dimensionality == (u.joule / u.mole / u.kelvin).dimensionality
```

`State.fugacity` is either `None` or an immutable `Fugacity` projection with:

- `component_ids`: the mixture component order;
- `ln_coefficient`: `ln(phi_i)` in that order;
- `coefficient`: `phi_i` in that order; and
- `value`: component fugacities `x_i phi_i P` as Pint pressure quantities.

```python
fugacity = closed.fugacity
if fugacity is not None:
    assert fugacity.component_ids == closed.component_ids
    assert len(fugacity.value) == len(closed.component_ids)
```

## Electrolyte reference and activity

`Mixture.reference(T, P, solvent_x, phase="liquid")` builds the explicit
liquid reference sequence for one source-complete cation/anion formula unit.
`solvent_x` is normalized in the full mixture component order and is zero for
ionic components. Only the liquid phase is admitted by this reference path.

```python
from epcsaft import Mixture, Parameters, unit_registry as u

electrolyte = Mixture(
    Parameters.from_catalog(
        "figiel-2025-reference-electrolytes",
        components=("water", "sodium-cation", "chloride-anion"),
    )
)
reference = electrolyte.reference(
    T=298.15 * u.kelvin,
    P=1 * u.bar,
    solvent_x=(1.0, 0.0, 0.0),
    phase="liquid",
)
state = reference.state(molality=0.1 * u.mole / u.kilogram)

assert state.activity is not None
assert state.activity.ion_component_ids == (
    "sodium-cation",
    "chloride-anion",
)
```

`reference.state(molality=...)` evaluates the generated electroneutral
composition at the reference temperature, pressure, and phase. The returned
`State.activity` contains:

- `formula_unit_molality` and generated `mole_fractions`;
- `single_ion` coefficients in `ion_component_ids` order;
- `mean_mole_fraction` and `mean_molality`; and
- `osmotic_coefficient`, which is available on the single-water path and is
  `None` when no admitted mixed-solvent osmotic standard exists.

The reference object also retains the source-bound infinite-dilution sequence,
`reference_molality`, convergence error, and solvation free-energy projection.
These are reference-state results, not reaction constants, species-equilibrium
results, fitted parameters, or catalog mutations.

## Stable errors and typing

Parameter loading raises `ParameterError`; state, density, and reference
evaluation raise `EosError`. Both derive from `ValueError` and expose immutable
`code` and `context` properties. Branch on those machine-readable properties
instead of parsing exception text.

```python
from epcsaft import ParameterError, Parameters

try:
    Parameters.from_catalog("not-installed", components=("sample",))
except ParameterError:
    pass
```

The wheel includes the PEP 561 `py.typed` marker. Runtime validation remains
the scientific authority for dimensions, normalization, domains, source
completeness, density-root selection, and reference-state applicability.

## Scientific authority and retained history

CppAD is the sole nonlinear production derivative authority. Native C++ owns
the hard-chain, dispersion, association, Debye--Huckel, Born SSM+DS, and
density equations; Python validates inputs and presents typed results. The
provider owns state/property evaluation and density closure, not phase
equilibrium, discovery, regression, Ceres, cross-package validation, or
reaction/species equilibrium.

Accepted receipt `promotion-0012-provider-full-candidate` governs the exact
full-provider scope. It transfers production authority but does not publish a
package release or broaden capability. Migration receipt
`promotion-0015-provider-native-sdk-v1` is the sole authority record for the
ADR 0001 native capsule; this guide does not transfer that authority. Each
native consumer requires its own cross-wheel gate.

The append-only neutral-mixture, active-pair, scalar pure-parameter,
association, full-electrolyte, molar-volume-domain, packing, active-solvation,
active-permittivity, and active-reference entries remain bounded provider
candidates unless their own receipts say otherwise. They own no targets,
bounds, Ceres, fit acceptance, phase identity, equilibrium, or catalog writes.
Gross--Sadowski methane/ethane is first reference evidence, not a mixture
allowlist. The private Held water-ionization candidate retains
`authority_effect: none`; Equilibrium owns `K_w`, reaction conversion,
constraints, Ipopt, and certificates.

The [native consumer ADR](adr/0001-native-consumer-transport.md), [0.2
frontend ADR](adr/0002-python-frontend-v0.2.md), and [0.2 migration
record](migration-0.2.md) preserve the transport and migration history.
