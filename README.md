# ePC-SAFT

Private scientific monorepo for the unified ePC-SAFT engine, canonical data,
compact numerical evidence, and retained scientific analyses.

The independent repository is
[tannerpolley/ePC-SAFT-greenfield](https://github.com/tannerpolley/ePC-SAFT-greenfield).
Its `main` baseline preserves the accepted Core history from `cea695f`.
The original `tannerpolley/ePC-SAFT` repository remains available for historical
issues and reviews; current development and delivery use this repository.

Start with the [project map](docs/scientific/CONTEXT.md) for the purpose,
scientific owners and one equilibrium driver. The
[implementation plan](docs/scientific/implementation-plan.md) distinguishes
required capability from implemented and numerically qualified behavior.

For the live Map, Plan and Results website, run
`quarto preview . --host 127.0.0.1 --port 8766 --no-browser --render all --no-execute`
and open the
[Quarto companion](http://127.0.0.1:8766/). The companion is a read-only
viewer: Quarto renders the project map, plan, ordinary check reports, and
retained numerical results. Starting or refreshing it never builds the Engine,
collects tests, or runs a scientific calculation.

To keep the companion running for this checkout, install the user service and
enable it:

```bash
install -Dm644 tools/systemd/epcsaft-companion-greenfield.service \
  "$HOME/.config/systemd/user/epcsaft-companion-greenfield.service"
loginctl enable-linger "$USER"  # start the user service at boot without login
systemctl --user daemon-reload
systemctl --user enable --now epcsaft-companion-greenfield.service
```

Check it with `systemctl --user status epcsaft-companion-greenfield.service`.
The unit serves from this repository root and restarts after a service failure.

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
from epcsaft import equilibrium
```

Do not use an editable install or add `engine/src` to `PYTHONPATH`. Downstream
applications install an explicitly identified wheel and record its SHA-256.

See [`engine/README.md`](engine/README.md) for direct Meson compile diagnostics
and canonical wheel build and installation guidance.

## Checks

Run the affected-area gate during normal development:

```bash
bash tools/check-affected.sh
```

Run the fast installed-wheel gate with:

```bash
bash tools/check-all.sh
```

`check-affected.sh` selects checks from the changed paths. Data, analysis, or
documentation-only changes can stop at readiness checks without building or
installing the wheel; selected Engine checks run setup first.
`check-all.sh` is the full fast installed-wheel gate: it always installs the
current wheel, then covers the default EOS, equilibrium, static,
and readiness checks. It excludes analysis reproductions and slow or expensive
numerical evidence. Regression checks join the gate when that capability exists.

Native C++ changes use the explicit owner commands:

```bash
bash engine/tools/check-native.sh
bash engine/tools/check-native-sanitizers.sh
```

The first command builds the native extension with Clang and runs the full
clang-tidy sweep. The second builds a separate debug ASan/UBSan wheel and runs
its compact numerical tests. Run slow Engine tests or an analysis reproduction
only when the changed scientific claim requires that evidence.

Run a canonical end-to-end workflow with:

```bash
bash tools/run-slice.sh khudaida
```

See `AGENTS.md`, `ARCHITECTURE.yaml`, and
`docs/scientific/CONTEXT.md` for execution, architecture, and scientific
ownership.
