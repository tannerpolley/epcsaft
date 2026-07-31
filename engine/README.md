# ePC-SAFT EOS

`epcsaft` is a compact, source-traceable equation-of-state package. The
0.2.0.dev0 Python frontend has one public state/property path:

```text
Parameters -> Mixture -> State
                         |
                         -> Fugacity
                         -> Activity (for a reference-created electrolyte state)
```

One resolved parameter selection feeds one native residual Helmholtz kernel for
hard-chain, dispersion, association, Debye--Huckel, and corrected Born SSM+DS
contributions. The frontend is immutable and has no compatibility aliases,
backend selectors, mutable scientific dictionaries, hidden physical defaults,
or duplicate equation owners.

## Installation and Python use

`epcsaft` requires Python 3.11 or newer and is not published to a package
index. Install a supplied wheel with `python -m pip install /path/to/wheel`, or
build the current checkout with `python -m pip install .`; source builds require
a C++17 compiler, CMake 3.20 or newer, and Ninja. Pint is the sole Python
runtime dependency.

The [Python API guide](docs/python-api.md) is the canonical walkthrough for
catalog IDs, parameter bundles, Pint state inputs, density and pressure
closure, residual property bases, fugacity, electrolyte references, activities,
and stable errors. The guide's catalog table records the currently installed
source-backed selections and their evidence domains.

```python
from epcsaft import Mixture, Parameters, unit_registry as u

parameters = Parameters.from_catalog(
    "gross-2001-methane-ethane",
    components=("methane", "ethane"),
)
mixture = Mixture(parameters)

state = mixture.state(
    T=250 * u.kelvin,
    P=10 * u.bar,
    x=(0.4, 0.6),
)

assert state.component_ids == ("methane", "ethane")
assert state.ares().check("[energy] / [substance]")
assert state.fugacity is not None
```

`Parameters.from_catalog(bundle_id, components=..., version=None)` selects the
requested component order from an installed bundle. Omitting `version` chooses
the greatest installed integer version for that bundle ID; pass a version to
pin a reproducible catalog identity. `Parameters.from_bundle(path,
components=...)` loads the strict version-1 directory schema. Use
`Parameters.from_records(...)` for immutable caller-owned typed records, and
`Parameters.create_bundle_template(...)` to create a non-overwriting authoring
skeleton. Caller-owned values retain their provenance and validity domains and
are not promoted into the installed catalog.

`Mixture.state(T=..., P=..., x=..., phase=...)` and
`Mixture.state(T=..., rho=..., x=...)` accept exactly one of pressure or molar
density. `phase` is optional and applies only to pressure closure. Residual
methods on `State` return molar Pint quantities by default or reduced scalar
values with `basis="reduced"`:

| Method | Molar result | Reduced result |
| --- | --- | --- |
| `ares()` | J/mol | `Ares/(nRT)` |
| `ures()` and `hres()` | J/mol | `Ures/(nRT)`, `Hres/(nRT)` |
| `gres()` | J/mol | `Gres/(nRT)` |
| `sres()` | J/(mol K) | `Sres/(nR)` |
| `cvres()` and `cpres()` | J/(mol K), when defined | value divided by `R`, when defined |

`State.fugacity` is an immutable component-order projection containing
`ln_coefficient`, `coefficient`, and component fugacity `value` quantities when
the fugacity projection is defined. `Mixture.reference(T=..., P=...,
solvent_x=..., phase="liquid")` constructs the admitted liquid electrolyte
reference. Calling `reference.state(molality=...)` returns a `State` whose
`activity` contains single-ion, mean mole-fraction, mean molality, and—on the
single-water path—practical osmotic-coefficient projections.

## Scientific and authority boundaries

The [checked equation mirror](docs/equations.md) maps retained mathematical
relations to their C++ or Python implementation owner and to the [normative
scientific contract](docs/science/full-provider-eos.md). The implementation
remains executable authority; generated Markdown and LaTeX views are
synchronized documentation mirrors.

EOS evaluates one state at a time. It does not perform phase
equilibrium or phase discovery, reaction/species equilibrium, regression,
parameter fitting, ideal-gas caloric calculations, or multi-salt reference
paths. Solver convergence is not a claim of global phase stability or
literature validation.

Pressure closure accepts exactly one of molar density or specified pressure. An
unlabelled pressure state must have one mechanically stable root. A genuine
multi-root state requires `phase="vapor"` or `phase="liquid"`; a phase hint
never renames a lone root.

The installed catalog evidence is:

- Crespo et al. neutral 2B PC-SAFT records for pure DEEA, with the separate
  PIL volume correction excluded;
- Gross--Sadowski methane/ethane and pure-propane neutral PC-SAFT records;
- the bounded Figiel--Yu--Held water/methanol/ethanol chloride/bromide source
  bundle and a distinct aqueous Na+/K+/I- source bundle;
- the distinct Ascani et al. original-Born water/1-butanol/Na+/K+/Cl- bundle at
  298.15 K.

The Held et al. source-complete aqueous-electrolyte bundle records water and
all 29 Table 2 ions at the reported 298.15 K and 1 bar reference conditions.
It remains non-authoritative for reaction equilibrium and carries no `K_w`,
ion-pairing, or Born-term claim.

The accepted authority receipt `promotion-0012-provider-full-candidate`
promotes the stated EOS scope. The receipt identifier is retained as
historical evidence. It is a
production-authority decision, not a package release and not a claim beyond
the recorded capability boundaries.

## Separate native consumer boundary

The installed native SDK is a separate, narrow consumer boundary. Its accepted
transport is documented in [ADR 0001](docs/adr/0001-native-consumer-transport.md)
and its authority is determined only by the historical Migration receipt
`promotion-0015-provider-native-sdk-v1`; this Python frontend documentation
does not transfer that authority. The capsule evaluates one stateless phase
block at a time and exposes no EOS implementation types or equations. Consumers
retain phase identity, equilibrium constraints, solver
policy, targets, bounds, scaling, and fit acceptance.

The fixed SDK prefix owns pure, nonassociating phase values and CppAD
derivatives through third order in amount and volume. Its append-only
parameterized entry owns exact pure-phase derivatives with respect to `m`,
`sigma`, and `epsilon/k` for the first methane regression slice. Each consumer
still requires its own installed-wheel gate and promotion receipt.

The local append-only neutral binary-mixture, active-`k_ij`, active-`l_ij`,
scalar pure-parameter, association, full-electrolyte, volume-domain, packing,
active-Born-diameter, active-permittivity, active-solvation-factor, and
active-reference candidates remain bounded evidence or historical Migration
records. They
do not add Python aliases, catalog authority, equilibrium, regression target
ownership, fitting policy, or broader model admission. Gross--Sadowski
methane/ethane is first reference evidence, not a mixture allowlist.

`authority_effect: none` remains the status of the Held aqueous reference
candidate. Equilibrium exclusively owns `K_w`, reactive standard-state
conversion, reaction and balance constraints, Ipopt, and certificates.

The architecture and cost ratchet are recorded in `ARCHITECTURE.yaml`; the
scientific equation owner is `docs/science/full-provider-eos.md`. The complete
0.2 frontend decision and migration record are [ADR 0002](docs/adr/0002-python-frontend-v0.2.md)
and [migration-0.2](docs/migration-0.2.md).
