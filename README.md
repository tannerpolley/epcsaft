# ePC-SAFT

Private scientific monorepo for the unified ePC-SAFT engine and its governed
data, numerical evidence, validation campaigns, and research oracles.

## Prerequisites

The supported development environment is Linux with Git, Python 3.13, `uv`, a
C and C++ compiler, Ninja, CMake, pkg-config, and `flock`. The locked root
environment supplies Meson, meson-python, pybind11, and `patchelf`. The native
graph also needs the environment-provided Ipopt and Eigen development
packages. The first build needs network access for the four
checksum-pinned Meson subproject archives. The native owner checks additionally
require `clang-18`, `clang++-18`, `clang-tidy-18`, `run-clang-tidy-18`, `file`,
and `nm`. `ccache` is optional.

## Setup and installation

Run commands from the repository root:

```bash
bash .codex/environments/setup.sh
bash tools/status.sh
```

Setup syncs the locked root `.venv`, builds `engine/` through meson-python,
writes one verified wheel under `build/environment-wheel/`, and installs that
wheel non-editably. Re-running setup reuses checkout-local Ninja state,
`ccache`, and the content-addressed machine-local wheel cache when their full
build fingerprints match.

Python consumers use the installed package and do not invoke Meson:

```python
import epcsaft
from epcsaft import equilibrium, regression
```

Do not use an editable install or add `engine/src` to `PYTHONPATH`. Downstream
applications install an explicitly identified wheel and record its SHA-256.

See [`engine/README.md`](engine/README.md) for direct Meson compile diagnostics
and canonical wheel build and installation guidance.

## Validation

Run the affected-area gate during normal development:

```bash
bash tools/check-affected.sh
```

Run the parallel fast installed-wheel gate with:

```bash
bash tools/check-all.sh
```

`check-affected.sh` selects checks from the changed paths. Data, governance, or
documentation-only changes can stop at readiness checks without building or
installing the wheel; selected Engine or Validation checks run setup first.
`check-all.sh` is the full fast installed-wheel gate: it always installs the
current wheel, then covers the default EOS, equilibrium, regression, static,
and readiness checks. It excludes Validation campaigns and slow or expensive
numerical evidence.

Native C++ changes use the explicit owner commands:

```bash
bash engine/tools/check-native.sh
bash engine/tools/check-native-sanitizers.sh
```

The first command builds both private extensions with Clang and runs the full
clang-tidy sweep. The second builds a separate debug ASan/UBSan wheel and runs
its compact numerical tests. Run slow Engine tests or a Validation campaign
only when the changed scientific claim requires that evidence.

Run a canonical end-to-end workflow with:

```bash
bash tools/run-slice.sh khudaida
```

See `AGENTS.md` for agent routing and `governance/GOVERNANCE.md` for the
scientific and architectural doctrine.
