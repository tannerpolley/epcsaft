# Unified ePC-SAFT Engine

`epcsaft` is one installable thermodynamic engine for parameter handling, EOS
evaluation, phase/species equilibrium, and parameter regression:

```text
Data packet -> Parameters -> Regression -> fitted Parameters
                          -> EOS State -> Equilibrium
```

The distribution has one version and three public module surfaces:

```python
import epcsaft
from epcsaft import equilibrium, regression
```

There are no separate equilibrium or regression packages and no compatibility
imports.

## Install

Python 3.13, CMake, Ninja, a C++17 compiler, Ipopt, Ceres, and their native
dependencies are required for a source build.

```console
python -m pip install .
```

The one built wheel contains the Python modules, three compiled extension
modules, and the public native headers. Downstream repositories must install an
immutable wheel and hash; they must not install this checkout in editable mode.

## End-to-end use

```python
from epcsaft import Mixture, Parameters, SourceRecord, ValidityDomain
from epcsaft import equilibrium, regression, unit_registry as u

parameters = Parameters.from_bundle(
    "/materialized/data-packet/parameters",
    components=("methane", "ethane"),
)
model = Mixture(parameters)
state = model.state(
    T=250 * u.kelvin,
    rho=1000 * u.mole / u.meter**3,
    x=(0.4, 0.6),
)
```

Regression workflows return immutable typed results. A valid result transfers
fitted values back to the EOS only through
`RegressionResult.to_parameters(...)`, with explicit source, validity, bundle
identity, and locator metadata. Equilibrium then consumes the resulting
`Mixture` directly.

## Data and evidence

Reusable scientific packets remain in the Data repository. Engine never
discovers a sibling checkout; callers materialize exact packet bytes and pass
an explicit path. Test fixtures include receipts that bind the Data commit,
packet fingerprint, and every file hash.

Validation runs black-box campaigns against one immutable Engine wheel. Lab and
applications consume the same accepted wheel. Historical plans and receipts
are preserved as provenance but are not live installation instructions.
