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

`build/report.html` shows the latest test results, key-solve work and wall-time trend per
main commit, and the analysis index, each with the commit and time of its data. `check-all`
writes it, `tools/report` rewrites it on demand, and `xdg-open build/report.html` opens it.

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

Setup syncs the locked root `.venv` and installs `engine/` editable through
meson-python, building into `build/engine-meson/shared` with `ccache` and the
shared HiGHS.
`import epcsaft` then rebuilds whatever native sources changed (Ninja's incremental
build, printed to stderr when there is work); a Python edit needs no rebuild and no setup.
Re-run setup after a lock change. Setup rebuilds the tree from scratch after a compiler
or HiGHS change; after another system library upgrade, delete `build/engine-meson/shared`
and re-run setup.

```python
import epcsaft
from epcsaft import equilibrium
```

`tools/build-engine.sh` owns every Engine build. Wheels exist only where an exact
wheel file and its SHA-256 are recorded: `tools/build-engine.sh --wheel <dir>` for application pins
and A/B timing, and `tools/build-engine.sh --install-wheel` for an analysis that
records an installed wheel (re-run setup afterwards to return to the editable
install). `epcsaft.provenance()` identifies the running Engine: the commit and
whether `engine/` is clean for the editable install, the wheel path and SHA-256
for a wheel. Downstream applications install an explicitly identified wheel and
record its SHA-256.

## Checks

Run the affected-area gate during normal development:

```bash
tools/heavy bash tools/check-affected.sh [base]
```

Run the fast full gate with:

```bash
tools/heavy bash tools/check-all.sh
```

`check-affected.sh` selects checks from the paths changed since `base` (default: the
merge base with `origin/main`, printed first) plus uncommitted files. Data, analysis, or
documentation-only changes can stop at readiness checks without building the
Engine; selected Engine checks run setup first.
`check-all.sh` is the full fast gate: it runs setup, which rebuilds what changed, then covers the default EOS, equilibrium, static,
and readiness checks. It excludes analysis reproductions and slow or expensive
numerical evidence. Both gates run strict mypy over `engine/src` and `engine/tests`, so
pytest alone misses type errors in a test edit.
`check-affected.sh` runs the 19 Explicit cases (16 binary128 reactive-caloric
and reactive-cp cases, 2 envelope traces and the reacting Cp fit) only when their code
changes; request them with
`tools/heavy bash engine/tools/run-native-tests.sh 'equilibrium/Explicit/*' 'regression/Explicit/*'`.

Develop and debug on native programs alone: `tools/native-run` rebuilds one program
in the checkout's Meson build directory and runs it, whole or for one case, without
Python. New studies can be small C++ drivers that write
CSV rows; plots read the CSV. Run `check-all.sh` and `check-affected.sh`, which
cover the Python interface too, before review.

```bash
tools/heavy tools/native-run equilibrium_checks --list
tools/heavy tools/native-run equilibrium_checks caloric-admission
tools/heavy tools/native-run equilibrium/Automatic_fast/e5-neutral
```

Find a claim's case by its prefix, e.g. `grep -n '"e5-9/' engine/tests/native/*.cpp`, or with
`--list`, which prints each case's test id, arguments and alias; `native-run` also
takes the Meson test name that `run-native-tests.sh` prints. `equilibrium_checks`,
`properties_checks` and `regression_checks` are each built from one unit per case family
(e.g. `equilibrium_e5_checks.cpp`) with shared helpers in `<program>_shared.hpp`;
`<program>.cpp` holds only the case registry, so an edit to one case recompiles one unit.

`tools/native-isolated` builds one program in a temporary detached worktree with its own
Meson build directory, leaving this checkout's source and build directory free:

```bash
tools/heavy tools/native-isolated falsify equilibrium_checks E2-2 plant.diff  # passes only if the planted case fails
tools/heavy tools/native-isolated compare regression_checks ionic-addition origin/main  # rows at 1e-10 against a commit
```

Native C++ changes use the explicit owner commands:

```bash
tools/heavy bash engine/tools/check-native.sh
tools/heavy bash engine/tools/check-native-sanitizers.sh
```

The first command builds the native extension with Clang and runs the full
clang-tidy sweep. The second builds a separate debug ASan/UBSan wheel and runs
its compact numerical tests. Run slow Engine tests or an analysis reproduction
only when the changed scientific claim requires that evidence.

Each checkout has its own `.venv` from setup with an editable Engine install.
Run builds and Engine calculations as `tools/heavy <cmd>`, which admits one
numerical thread when memory allows. It stops the whole command after `HEAVY_TIMEOUT`
(default 30m); set `HEAVY_TIMEOUT=2h` for a long chain such as a series of falsify runs.
A warm `check-all.sh` takes about 20–40 s
(24.0 s measured in #319). `tools/heavy tools/time-key-solves <wheel-A> <wheel-B>` reports work counts,
the no-slowdown gate for ionic LLE and the MEA flash, with paired wall times; build the
checkout's side with `tools/heavy bash tools/build-engine.sh --wheel build/wheels/<name>`.
After merges, `tools/heavy tools/refresh-main-wheel` builds main's
wheel in the main checkout, compares its key-solve work with the previous
`build/wheels/main-current/`, replaces that directory and its README, and keeps
`build/wheels/main-<commit>/`, the applications' stable pin.
There is no scientific CI/CD.

Run a canonical end-to-end workflow with:

```bash
tools/heavy bash tools/run-slice.sh khudaida
```

See `AGENTS.md`, `ARCHITECTURE.yaml`, and
`docs/scientific/CONTEXT.md` for execution, architecture, and scientific
ownership.
