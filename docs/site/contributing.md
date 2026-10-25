# Contributing

Keep each calculation's equations, source inputs, units, numerical comparisons
and limitations visible. Implement thermodynamics in the native Engine;
Python declares inputs and reads results. Follow the
[scientific ownership map](https://github.com/tannerpolley/ePC-SAFT/blob/main/docs/scientific/CONTEXT.md)
and [coding standards](https://github.com/tannerpolley/ePC-SAFT/blob/main/CODING_STANDARDS.md)
when proposing a change.

## Prepare and check a checkout

Use [installation](installation.md) for prerequisites and source setup. Each
checkout has its own `.venv` with an editable Engine. Run commands from the
repository root. The
[README Checks](https://github.com/tannerpolley/ePC-SAFT/blob/main/README.md#checks)
describes the commands; the
[scientific gates](https://github.com/tannerpolley/ePC-SAFT/blob/main/docs/agents/workflow.md#gates)
own which checks a change requires.

```bash
nice -n 10 bash tools/check-affected.sh origin/main
nice -n 10 bash tools/check-all.sh
nice -n 10 bash tools/build-docs.sh
```

The documentation command installs the locked optional documentation packages,
uses the checkout's Engine with MUMPS, executes every gallery notebook afresh,
and builds HTML and cheap doctests with warnings treated as errors. Logs,
rendered pages and timing are under `build/docs/`; the HTML home is
`build/docs/html/index.html`. Notebook execution errors fail the build.

## Document the source-owned meaning

Write NumPy-style docstrings beside Python definitions and native declaration
prose in `_core.pyi`. Explain each input's physical quantity, unit, order and
basis; explain returned availability and failure causes. Use one description
owner per field. Re-exported names point to their definition instead of
repeating the prose. The [reference pages](api/index.rst) render those sources.

Use a small doctest for a declaration or inexpensive identity. Put complete
state, equilibrium, discovery and fitting calculations in the gallery.
Every displayed code fragment must run with the installed Engine. Do not
hide a numerical solve in a test helper, a notebook utility or Python loops.
The [quickstart](quickstart.md) shows the short declaration style.

## Write a reproducible notebook

State the physical question, source locators, conditions and units before
calling the Engine. Resolve the checkout root and read canonical parameter
packets or observation inputs; the wheel does not contain those data.
Retain the source's use restrictions. Use original model-output plots and
small factual comparisons only where their source basis permits that use.
A code licence does not grant redistribution rights over every source table.

Declare requests in Python, then read native values. The Engine owns density
roots, association, phase/reaction equilibrium, fitting, continuation and
stability searches. Finite differences may verify exact derivatives.
Pyomo/IDAES process equations consuming Engine derivatives through ASL are the
explicit exception described in [process models](guide/process-models.md).

Check balances, residuals and an independent identity or reference comparison
with stated numerical tolerances. A termination word is not a numerical
result. Distinguish numerical verification, parameter calibration and physical
validation; keep uncertainty and source extrapolation explicit. The
[water](examples/parameters-and-water-state.ipynb),
[propane](examples/propane-parameter-estimation.ipynb) and
[ammonia](examples/ammonia-synthesis-equilibrium.ipynb) examples illustrate
those different evidence claims.

Keep theory in its canonical scientific document. The
[theory pages](theory/index.md) include those documents; do not copy equations
into a second independently maintained manual. Capability prose must agree
with the installed implementation and its demonstrated limits.

## Optional process-example tools

The [IDAES settler](examples/idaes-two-liquid-settler.ipynb) needs optional
external solver tools in addition to the locked Python packages. Obtain IDAES
extensions 3.4.2 solver and library archives directly from
[upstream IDAES](https://github.com/IDAES/idaes-ext/releases/tag/3.4.2).
The documentation workflow verifies the archives and extracts them into an
external temporary directory. For local execution, point `IDAES_DATA` and
`PYOMO_CONFIG_DIR` to an equivalent prepared directory outside the checkout,
with the tools under `bin/`. Preserve upstream notices.

The Engine documentation build selects MUMPS with `EPCSAFT_HSL_LIBRARY` and
`EPCSAFT_IPOPT_OPTIONS_FILE` unset. The optional outer Ipopt also selects
MUMPS. The upstream Ipopt executable contains HSL even under that selection;
solver selection does not remove the embedded code.

These tools are build-time only. Keep every solver executable, library,
archive and environment containing HSL outside the repository, Engine builds,
Engine wheels, documentation outputs, downloads, uploaded site payload and
shared caches. Upload only rendered pages, original notebooks, figures,
numerical records and text logs. Computational outputs may be displayed
without redistributing the solver. Readers obtain their optional tools
directly from upstream. A separately built HSL-free outer solver would need
its own qualification.

HSL, a collection of Fortran codes for large-scale scientific computation. See https://www.hsl.rl.ac.uk/.
