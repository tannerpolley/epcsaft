# ePC-SAFT

Compute fluid properties, declared phase and reaction equilibrium, exact
directional derivatives and parameter estimates with the native C++ Engine.
Python supplies the physical declarations and reads values and diagnostics.
Start with [installation](installation.md) and the [quickstart](quickstart.md),
then use the [example gallery](examples/index.md) for executed calculations.
The [equations](theory/equations.md) and
[algorithms](theory/algorithms.md) define the thermodynamics and methods.

Follow the [quickstart](quickstart.md) for a first state declaration and the
[water notebook](examples/parameters-and-water-state.ipynb) for its source,
pressure-closure check and prediction limits.

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
Special-case appendix <examples/index>
api/index
theory/index
limits
contributing
```
