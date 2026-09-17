# Unified ePC-SAFT Engine

One installed Python distribution owns the Engine. Current state/equilibrium imports:

```python
import epcsaft
from epcsaft import equilibrium
```

Consumers import `epcsaft` from their environment. They do not run Meson or
load modules from this source tree.

## Canonical build and installation

From the repository root, build and install the current source non-editably:

```bash
bash .codex/environments/setup.sh
```

Setup keeps the checkout-local Ninja build warm, uses `ccache` when available,
and reuses a verified machine-local wheel keyed by source, ABI, platform, and
toolchain and dependency contents. It builds and installs before any test
starts. The root [`README.md`](../README.md) owns routine, full, native,
sanitizer, and scientific validation commands.

## Direct Meson development

Run root setup once so the locked Python build tools and native prerequisites
are available. A direct Meson build is useful for compile diagnostics:

```bash
export PATH="$PWD/.venv/bin:$PATH"
meson setup build/engine-dev engine
meson compile -C build/engine-dev
```

Source changes need only `meson compile`. Re-read changed build definitions or
dependency metadata with:

```bash
meson setup --reconfigure build/engine-dev engine
meson compile -C build/engine-dev
```

Use a clean configure after changing the compiler or fundamental Meson options:

```bash
meson setup --wipe build/engine-dev engine
meson compile -C build/engine-dev
```

The Meson graph includes native scientific checks. Run installed-wheel Python
checks through the root commands; expensive native evidence remains explicit. Do not use `meson
install`; the wheel is the installation boundary.

The supported direct native test route records one coherent outcome through the
shared recorder and defaults to the fast suites:

```bash
bash engine/tools/run-native-tests.sh
```

It runs `contributions`, `properties`, `diagnostics`, `specification` and
`surface_example`; the equilibrium, boundary, and Properties-floor suites remain
explicit owner commands.

## C++ usage and debugging

`tests/native/surface_example.cpp` is an ordinary Python-free caller of the
native archive. It includes only `surface/parameter_builder.hpp` and
`equilibrium/problem.hpp` and shows, in order: a `ParameterDocument` with
explicit units (species index = position in `component_ids`), one reusable
`HelmholtzModel`, a full `State` in the T-P and T-rho charts, selected
observables with first- and second-order actions, a bubble point and a
two-phase flash, exact solved-state actions with the caller-formed fugacity
`a_i R T (1 mol/m^3)`, a partial result without thermochemistry, and two
inspectable failures. Its mixture is manufactured. It returns nonzero when one
of its identities fails, and it links `epcsaft_native`; there is no installed
SDK or ABI.

```bash
meson compile -C build/engine-dev surface_example
meson test -C build/engine-dev surface_example --print-errorlogs
build/engine-dev/surface_example            # prints every value with its unit
```

For a debugger, configure an unstripped debug tree once and run the same target:

```bash
meson setup build/engine-debug engine -Dbuildtype=debug -Dstrip=false
meson compile -C build/engine-debug surface_example
gdb --args build/engine-debug/surface_example
```

Create and install wheels only through root setup. When the root environment is
already synced, `bash tools/ensure-engine-wheel.sh` owns the canonical cached
wheel build and installation.

## Native dependency graph

The package owns one Meson graph and no installed cross-module SDK. Data
packet materializations are explicit and hash-bound. Installed-wheel analyses
live in `../analyses/`. See `AGENTS.md`, `../ARCHITECTURE.yaml`, and
`../docs/scientific/CONTEXT.md`.

Meson discovers threads natively; Ipopt 3.11.9 or newer and Eigen 3.4 or newer
through pkg-config. Exact archive wraps own Boost 1.88.0, CppAD
20250000.3, HiGHS 1.15.1, and NLopt 2.11.0. Meson's CMake module translates
the unmodified NLopt project. CMake remains a third-party translation tool for
NLopt; the Engine has no CMake build route.

The HiGHS tag archive replaces its same-version source release asset because
that asset embeds Git metadata. All archive URLs and SHA-256 values are pinned.
Downloads and extracted trees stay ignored under `engine/build/subprojects/`
and never enter the wheel.

## Scientific specification and availability

The [scientific owner map](../docs/scientific/CONTEXT.md) identifies the
permanent equations and methods. The
[GitHub issue roadmap](https://github.com/tannerpolley/ePC-SAFT-greenfield/milestones)
records implementation gaps, including the absent regression module. Do not infer
availability from historical usage examples or a registered equation.
