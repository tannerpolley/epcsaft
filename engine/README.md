# ePC-SAFT EOS

`epcsaft` is the typed equation-of-state engine for the ePC-SAFT project. Its
public lifecycle is:

```text
Parameters -> Mixture -> State
                         |
                         -> Fugacity
                         -> Activity (electrolyte reference states)
```

EOS owns parameter resolution, residual Helmholtz evaluation, derivatives,
density closure, properties, and the model-bound native SDK. It does not own
parameter discovery, fitting, phase equilibrium, or literature validation.

## Install

Python 3.11+, CMake 3.20+, Ninja, and a C++17 compiler are required for a source
build. Pint is the sole Python runtime dependency.

```console
python -m pip install .
```

## Construct and evaluate a model

Tests and small applications can pass compact fixed-unit dictionaries directly:

```python
from epcsaft import Mixture, Parameters, unit_registry as u

raw = {
    "schema": "epcsaft.parameters",
    "schema_version": 1,
    "components": ("methane", "ethane"),
    "parameters": {
        "mw": (0.016043, 0.03007),
        "m": (1.0, 1.6069),
        "s": (3.7039, 3.5206),
        "e": (150.03, 191.42),
        "k_ij": ((0.0, 0.0), (0.0, 0.0)),
    },
    "options": {"permittivity_model": "none"},
    "validity": {"kind": "unknown"},
}

mixture = Mixture.from_dictionary(raw)
state = mixture.state(
    T=250 * u.kelvin,
    rho=1000 * u.mole / u.meter**3,
    x=(0.4, 0.6),
)

assert state.component_ids == ("methane", "ethane")
assert state.ares().check("[energy] / [substance]")
```

`Parameters.from_bundle(path, components=...)` loads the same strict external
bundle format from an explicit path. Omitting `components` preserves the bundle
component order. `Parameters.create_bundle_template(...)` creates an authoring
skeleton. EOS has no installed parameter catalog and never locates or imports
`ePC-SAFT-data` at runtime.

The external loader accepts the historical purpose value retained by immutable
Data migration packets only as inert manifest metadata; it does not restore a
catalog route, lookup, identity, or admission rule.

Reusable parameter/property packets and provenance belong in
[`ePC-SAFT-data`](https://github.com/ePC-SAFT/ePC-SAFT-data). Paper-, table-, or
figure-specific campaign extracts may remain in Validation when they are local
campaign evidence. Applications materialize the desired Data packet themselves
and pass its parameter directory to `from_bundle`.

See [the Python API guide](docs/python-api.md) for the complete dictionary
schema, units, aliases, bundle route, errors, and property API. The SDK headers
are installed under `epcsaft/include/epcsaft`; native capability admission uses
resolved topology, required fields, component order, and validity—not source or
bundle fingerprints.

## Scientific boundary

The native kernel evaluates hard-chain, dispersion, association,
Debye--Hückel, and Born contributions. State pressure, fugacity, residual
caloric properties, and density diagnostics all derive from that one owner.
The open polar equation work will expose one aggregate `ares_polar`
contribution; it will not create separate Helmholtz result terms for each polar
interaction type.

The repository test suite is fully self-contained. Compact raw cases live in
`tests/_raw_cases.py`; external bundle files are created only in temporary
directories for I/O tests.
