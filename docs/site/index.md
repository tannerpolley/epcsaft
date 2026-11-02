# ePC-SAFT

A native C++ implementation of the PC-SAFT and electrolyte PC-SAFT (ePC-SAFT)
equations of state. It computes fluid properties, declared phase and reaction
equilibrium, phase stability tests, exact derivatives of solved states and
parameter estimates. Python supplies the physical declarations and reads
values and diagnostics; every iterative calculation runs in the Engine.

The residual Helmholtz energy follows PC-SAFT (Gross and Sadowski, 2001) with
association and optional polar terms. For electrolytes, ePC-SAFT adds a
Debye–Hückel term (Cameretti, Sadowski and Mollerup, 2005; Held, Cameretti and
Sadowski, 2008) with cation–anion dispersion (Held et al., 2014). Optional
later extensions add composition-dependent permittivity (Bülow, Ji and Held,
2019) and a Born term (Bülow, Ascani and Held, 2021; Figiel, Yu and Held,
2025). The [EOS specification](theory/eos.md) lists every source and
the [equation registry](theory/equations.md) gives each implemented equation
with its source equation number.

Start with [installation](installation.md) and the [quickstart](quickstart.md),
then the [examples](examples/index.md), which are executed when this site is
built. Parameters and reference states remain explicit; unavailable outputs
carry causes. Read the [limits](limits.md) before applying a case outside its
source domain or reading a finite stability search as global evidence.

````{grid} 1 2 2 2
:gutter: 3

```{grid-item-card} User guide
:link: guide/index
:link-type: doc

Declare parameters, states and equilibrium specifications.
```
```{grid-item-card} Examples
:link: examples/index
:link-type: doc

Executed notebooks with their sources, checks and limits.
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
Examples <examples/index>
api/index
theory/index
limits
contributing
```
