# ePC-SAFT

A native C++ implementation of the PC-SAFT and electrolyte PC-SAFT (ePC-SAFT)
equations of state. It computes fluid properties, declared phase and reaction
equilibrium, phase stability tests, exact derivatives of solved states and
parameter estimates. Python supplies the physical declarations and reads
values and diagnostics; every iterative calculation runs in the Engine.

The residual Helmholtz energy follows PC-SAFT (Gross and Sadowski, 2001) with
association, optional polar terms and, for electrolytes, the Debye–Hückel and
Born contributions of ePC-SAFT (Cameretti, Sadowski and Mollerup, 2005; Held,
Cameretti and Sadowski, 2008; Held et al., 2014) and later permittivity and
Born revisions. The [EOS specification](theory/eos.md) lists every source and
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
