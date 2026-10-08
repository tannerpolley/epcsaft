# ePC-SAFT

Private scientific monorepo for the unified ePC-SAFT engine, canonical data,
compact numerical evidence, and retained scientific analyses.

The independent repository is
[tannerpolley/ePC-SAFT](https://github.com/tannerpolley/ePC-SAFT).
Its `main` baseline preserves the accepted Core history from `cea695f`.
The earlier repository is archived as
[tannerpolley/ePC-SAFT-legacy-2](https://github.com/tannerpolley/ePC-SAFT-legacy-2)
for its historical issues and reviews.

Start with the [project map](docs/scientific/CONTEXT.md) for the purpose,
scientific owners and one equilibrium driver. The
[GitHub issue roadmap](https://github.com/tannerpolley/ePC-SAFT/milestones)
owns implementation plans, dependencies and delivery status. Each leaf issue
records its reuse decision, acceptance and evidence links. Repository scientific
records distinguish required capability from implemented and numerically
qualified behavior.

For the live Map, Roadmap and Results website, run
`quarto preview . --host 127.0.0.1 --port 8766 --no-browser --render all --no-execute`
and open the
[Quarto companion](http://127.0.0.1:8766/). The companion links to the GitHub
roadmap and renders the local project map, ordinary check reports, and retained
numerical results. Starting or refreshing it never builds the Engine, collects
tests, or runs a scientific calculation.

To keep the companion running for the main checkout (the unit's working directory), install the user service and
enable it:

```bash
install -Dm644 tools/systemd/epcsaft-companion.service \
  "$HOME/.config/systemd/user/epcsaft-companion.service"
loginctl enable-linger "$USER"  # start the user service at boot without login
systemctl --user daemon-reload
systemctl --user enable --now epcsaft-companion.service
```

Check it with `systemctl --user status epcsaft-companion.service`.
The unit serves from the main checkout root and restarts after a service failure.

## Prerequisites

The supported development environment is Linux with Git, Python 3.13, `uv`, a
C and C++ compiler, Ninja, CMake, pkg-config, and `flock`. The locked root
environment supplies Meson, meson-python, pybind11, and `patchelf`. The native
graph also needs the environment-provided Ipopt and Eigen development
packages. The first build needs network access for the two
checksum-pinned Meson subproject archives (CppAD and HiGHS). The native owner checks additionally
require `clang-18`, `clang++-18`, `clang-tidy-18`, `run-clang-tidy-18`, `file`,
and `nm`. `ccache` is optional.

## Setup and installation

Run commands from the repository root:

```bash
tools/heavy bash .codex/environments/setup.sh
bash tools/status.sh
```

Setup syncs the locked root `.venv`, builds `engine/` through meson-python,
writes one verified wheel under `build/environment-wheel/`, and installs that
wheel non-editably. Re-running setup reuses checkout-local Ninja state,
`ccache`, and the content-addressed machine-local wheel cache when the wheel
fingerprint (wheel inputs and toolchain, not native test sources) matches.

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
tools/heavy bash tools/check-affected.sh [base]
```

Run the fast installed-wheel gate with:

```bash
tools/heavy bash tools/check-all.sh
```

`check-affected.sh` selects checks from the paths changed since `base` (default: the
merge base with `origin/main`, printed first) plus uncommitted files. Data, analysis, or
documentation-only changes can stop at readiness checks without building or
installing the wheel; selected Engine checks run setup first.
`check-all.sh` is the full fast installed-wheel gate: it always installs the
current wheel, then covers the default EOS, equilibrium, static,
and readiness checks. It excludes analysis reproductions and slow or expensive
numerical evidence.
`check-affected.sh` runs the 19 Explicit cases (16 binary128 reactive-caloric
and reactive-cp cases, 2 envelope traces and the reacting Cp fit) only when their code
changes; request them with
`tools/heavy bash engine/tools/run-native-tests.sh 'equilibrium/Explicit/*' 'regression/Explicit/*'`.

Develop and debug on native programs alone: `tools/native-run` rebuilds one program
in the checkout's Meson build directory and runs it, whole or for one case, without
building or installing the wheel. New studies can be small C++ drivers that write
CSV rows; plots read the CSV. Run `check-all.sh` and `check-affected.sh`, which
cover the installed wheel and the Python interface, before review.

```bash
tools/heavy tools/native-run equilibrium_checks --list
tools/heavy tools/native-run equilibrium_checks caloric-admission
tools/heavy tools/native-run equilibrium/Automatic_fast/e5-neutral
```

Find a claim's case with `grep -n '"<claim>/' engine/tests/native/*.cpp`, or with
`--list`, which prints each case's test id, arguments and alias; `native-run` also
takes the Meson test name that `run-native-tests.sh` prints.

Native C++ changes use the explicit owner commands:

```bash
tools/heavy bash engine/tools/check-native.sh
tools/heavy bash engine/tools/check-native-sanitizers.sh
```

The first command builds the native extension with Clang and runs the full
clang-tidy sweep. The second builds a separate debug ASan/UBSan wheel and runs
its compact numerical tests. Run slow Engine tests or an analysis reproduction
only when the changed scientific claim requires that evidence.

Each checkout has its own `.venv` from setup with one non-editable Engine wheel.
Run builds and Engine calculations as `tools/heavy <cmd>`, which admits one
numerical thread when memory allows. A warm `check-all.sh` takes about 20–40 s
(24.0 s measured in #319). `tools/heavy tools/time-key-solves <wheel-A> <wheel-B>` reports work counts,
the no-slowdown gate for ionic LLE and the MEA flash, with paired wall times.
After merges, `tools/heavy tools/refresh-main-wheel` restores or rebuilds main's
wheel in the main checkout, compares its key-solve work with the previous
`build/wheels/main-current/`, and replaces that directory and its README.
There is no scientific CI/CD.

Run a canonical end-to-end workflow with:

```bash
tools/heavy bash tools/run-slice.sh khudaida
```

See `AGENTS.md`, `ARCHITECTURE.yaml`, and
`docs/scientific/CONTEXT.md` for execution, architecture, and scientific
ownership.
