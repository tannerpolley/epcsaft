# ePC-SAFT

Compute fluid properties, declared phase and reaction equilibrium, exact
directional derivatives and parameter estimates with the native C++ Engine.
Python supplies the physical declarations and reads values and diagnostics.
Start with [installation](installation.md) and the [quickstart](quickstart.md),
then use the [example gallery](examples/index.md) for executed calculations.
The [equations](theory/equations.md) and
[algorithms](theory/algorithms.md) define the thermodynamics and methods.

For a first liquid-water state at 298.15 K and 100,000 Pa, run from the
repository root in the prepared environment:

```python
from epcsaft import Parameters, Pure

water = Pure(Parameters.from_packets(
    "data/packets", [("held-cameretti-sadowski-2008/4", ["water"])]))
state = water.state(298.15, P=100_000.0, phase="liquid")
print(state.molar_density)  # mol/m³
```

The [water notebook](examples/parameters-and-water-state.ipynb) identifies the
published parameter source and checks pressure closure. The
[density algorithm](theory/algorithms.md#alg-density) defines the branch solve;
the [model reference](api/eos.rst) documents the call. This is a
source-parameter prediction with numerical closure verification, not a
comparison with measured water density.

Choose the next guide by the quantity you need. Parameters and reference
states remain explicit; unavailable outputs carry causes. Read
[limits](limits.md) before extending a case to another source domain or
interpreting a finite stability search as global evidence.

````{grid} 1 2 2 3
:gutter: 3

```{grid-item-card} User guide
:link: guide/index
:link-type: doc

Declare parameters, states and equilibrium specifications.
```
```{grid-item-card} Python reference
:link: api/index
:link-type: doc

Inspect Python declarations and the native interface.
```
```{grid-item-card} Theory
:link: theory/index
:link-type: doc

Read the source-owned equations, conventions and algorithms.
```
````

```{toctree}
:maxdepth: 2
:hidden:

installation
quickstart
guide/index
examples/index
api/index
theory/index
limits
contributing
```
