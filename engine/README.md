# Unified ePC-SAFT Engine

One installed Python distribution provides EOS, equilibrium, and regression:

```python
import epcsaft
from epcsaft import equilibrium, regression
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

The Meson graph has no standalone test targets. Run Python tests against the
installed wheel through the root check commands. Do not use `meson
install`; the wheel is the installation boundary.

Create and install wheels only through root setup. When the root environment is
already synced, `bash tools/ensure-engine-wheel.sh` owns the canonical cached
wheel build and installation.

## Native dependency graph

The package owns one Meson graph and no installed cross-module SDK. Data
packet materializations are explicit and hash-bound. Installed-wheel analyses
live in `../analyses/`. See `AGENTS.md`, `ARCHITECTURE.yaml`, and
`docs/scientific/CONTEXT.md`.

Meson discovers threads natively; Ipopt 3.11.9 or newer and Eigen 3.4 or newer
through pkg-config. Exact archive wraps own Boost 1.88.0, CppAD
20250000.3, HiGHS 1.15.1, and NLopt 2.11.0. Meson's CMake module translates
the unmodified NLopt project. CMake remains a third-party translation tool for
NLopt; the Engine has no CMake build route.

The HiGHS tag archive replaces its same-version source release asset because
that asset embeds Git metadata. All archive URLs and SHA-256 values are pinned.
Downloads and extracted trees stay ignored under `engine/build/subprojects/`
and never enter the wheel.

## Regression example

Regression exports owner-neutral native `least_squares` and `curve_fit`:

```python
from epcsaft import regression

result = regression.least_squares(owners, observations, initial_parameters)
print(result.parameters, result.physical_jacobian)
```

EOS and Equilibrium compile their own immutable owner handles and sidecars.
Regression owns only ordered observations, caller-defined discrepancy weights,
affine coordinates, bounds, and solver accounting. The result retains both
physical owner Jacobians and weighted optimizer-coordinate Jacobians; it never
promotes or adopts parameters. Installed-wheel literature evidence remains
owned by `../analyses/reproduction/`.

The canonical capability cutover and its EOS, Regression, and Equilibrium owner
plans begin at
[`docs/plans/2026-08-07-general-scientific-capability-hard-cutover.md`](docs/plans/2026-08-07-general-scientific-capability-hard-cutover.md).
