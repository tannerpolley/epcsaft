# ePC-SAFT

Scientific source distribution of the ePC-SAFT Engine, canonical data,
compact numerical evidence, and retained scientific analyses.

Public source releases are available at
[tannerpolley/epcsaft](https://github.com/tannerpolley/epcsaft/releases).
Start with the [project map](docs/scientific/CONTEXT.md) for the purpose,
scientific owners and equilibrium formulation. Retained analyses document
numerical verification and its limits; they do not establish experimental
validation for other systems or conditions.

## Prerequisites

The supported development environment is Linux with Git, Python 3.13, `uv`, a
C and C++ compiler, Ninja, CMake, pkg-config, and `flock`. The locked root
environment supplies Meson, meson-python, nanobind, and `patchelf`. The native
graph also needs the environment-provided Ipopt, Eigen and Ceres development
packages, at the versions in [engine/README.md](engine/README.md#native-dependency-graph). The first build needs network access for the two
checksum-pinned Meson subproject archives (CppAD and HiGHS). The native owner checks additionally
require `clang-18`, `clang++-18`, `clang-tidy-18`, `run-clang-tidy-18`, `file`,
and `nm`. `ccache` is optional. The line budgets in both check commands also need
`clang-format-18`, and count Python with the locked `.venv` Ruff, so run setup
before them.

## Public installation from source

Public releases are source releases from the public release repository. Download
and unpack the release source archive, then run these commands from its root.
The archive includes the Engine, documentation, retained analyses and paper
source; it contains no Git history. The conda-forge recipe is in `packaging/conda-forge/`; its submission and
binary availability are separate from this source release.

Use Linux, Python 3.13, `uv`, a C/C++ compiler, Ninja, CMake and pkg-config.
Install the native development packages listed in
[engine/README.md](engine/README.md#native-dependency-graph): Ipopt 3.11.9,
Eigen >=3.4, Ceres >=2.2.0 and <2.3.0, and their solver dependencies.
The first build downloads the checksum-pinned CppAD and HiGHS archives.

```bash
uv sync --python 3.13 --locked --all-packages --no-install-workspace
uv pip install --python .venv/bin/python --no-build-isolation ./engine
.venv/bin/python -c 'import epcsaft; print(epcsaft.__file__)'
```

The code uses GPLv3 with an additional permission under section 7 for linking
with EPL-licensed Ipopt; see [LICENSE](LICENSE). Retained third-party data have
separate terms, including non-commercial restrictions; see [NOTICE](NOTICE).
`tools/refresh-main-wheel` distributes private development wheels and requires
a private repository. Public source releases do not use that tool.

Maintainers prepare a local snapshot with `python3 tools/release-snapshot` after
committing the candidate. It writes a new tree and a JSON manifest under
`build/release/`, excludes ignored files and `build/`, and replaces files larger
than 50 MiB with SHA-256 notes. The manifest records the HSL and journal-PDF
path checks and personal-path counts. Complete the secret scan before any
publication; snapshot generation does not publish or tag a release.

## Setup and installation

Run commands from the repository root:

```bash
bash .codex/environments/setup.sh
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

`tools/build-engine.sh` owns ordinary editable and wheel builds;
`tools/check-pr-native` owns frozen instrumented preparation. Wheels exist only where an exact
wheel file and its SHA-256 are recorded: `tools/build-engine.sh --wheel <dir>` for application pins
and merged-main delivery, and `tools/build-engine.sh --install-wheel` for an analysis that
records an installed wheel (re-run setup afterwards to return to the editable
install). `epcsaft.provenance()` identifies the running Engine: the commit and
whether `engine/` is clean for the editable install, the wheel path and SHA-256
for a wheel. Downstream applications install an explicitly identified wheel and
record its SHA-256.

## Checks

Which checks a delivery runs, and when, is listed only in
[`docs/agents/workflow.md` "Gates"](docs/agents/workflow.md#gates); this section
describes the commands. Run the affected-area gate during normal development:

```bash
bash tools/check-affected.sh [base]
```

Run the fast full gate with:

```bash
bash tools/check-all.sh
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
`bash engine/tools/run-native-tests.sh 'equilibrium/Explicit/*' 'regression/Explicit/*'`.

Develop and debug on native programs alone: `tools/native-run` rebuilds one program
in the checkout's Meson build directory and runs it, whole or for one case, without
Python. New studies can be small C++ drivers that write
CSV rows; plots read the CSV. `check-all.sh` and `check-affected.sh` cover the
Python interface too.

```bash
tools/native-run equilibrium_checks --list
tools/native-run equilibrium_checks caloric-admission
tools/native-run equilibrium/Automatic_fast/e5-neutral
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
tools/native-isolated falsify equilibrium_checks E2-2 plant.diff  # passes only if the planted case fails
tools/native-isolated compare regression_checks ionic-addition origin/main  # rows at 1e-10 against a commit
```

Freeze a candidate commit, then run
`nice -n 10 tools/check-pr-native --pr N --head <full-sha> --base <merge-base-sha>`.
It shares one locked Clang debug ASan/UBSan Meson tree and environment under the
main checkout's `build/pr-gates/pr-N/<head>/`, including GCC's binary128 headers.
`engine/tools/check-native.sh --prepared <root> --changed <base>` owns lint selection;
`engine/tools/check-native-sanitizers.sh --prepared <root>` owns instrumented checks.
Both consume the prepared tree. Identical inputs reuse passing evidence without
compiling; toolchain/header/environment changes replace stale objects, while a
changed base or check selection reruns the checks. Logs and test records stay
beside that frozen source. See the [development method](docs/agents/workflow.md#development-method).

Each checkout has its own `.venv` from setup with an editable Engine install.
`tools/check-equilibrium-corpora <A> <B>` and `tools/time-key-solves <A> <B>` select
prepared Python interpreters or build/environment directories, verify loaded
Engine identities before calculations, and retain the existing corpus rows,
work counts, states and paired wall times. Use the same prepared environment for
both sides of an A=A control. The TPD baseline/replay, IDAES-ASL, density benchmark
and shared-derivative qualification readers use the same prepared-image identity.
Historical study readers remain wheel-only: the Khudaida `run_fixed_settler.py`
campaign and TPD `feos_compare.py` and `asl_native_trace.py`. These are outside the
current corpus command's case selection.
The [key-solve command](tools/time-key-solves.md) includes the frozen 149-row MEA
F1 residual-plus-Jacobian evaluation: no increased reference computations,
reference columns, density solves or equilibrium iterations, and serial CPU
at most 1.25 times its frozen value. Add `--fit` for the 120 s full-fit wall
ceiling and slowest-worker critical-path CPU; an A/B full fit needs the same-host
A=A `--fit` summary via `--control <summary.json>`. Full fits use native default
workers (32 on the reference host); the default command runs one serial F1
evaluation per environment. The Gates table owns when the full fit is required.
After merges, `tools/refresh-main-wheel` compares retained previous/candidate
prepared environments before packaging main. Select the retained baseline with
`EPCSAFT_PREVIOUS_PREPARED`; the candidate defaults to main’s prepared build tree.
Set `EPCSAFT_NUMERICAL_EVIDENCE` to space-separated existing evidence paths.
Wheels are `auditwheel`-repaired `manylinux_2_39` wheels carrying Ipopt, MUMPS and
Ceres; before retaining one it runs `tools/check-wheel-container.sh`, which refuses
any HSL member and runs `verify-installed-wheel.py` and the Python suite in a
digest-pinned `python:3.13-slim` container without those host libraries.
It retains the wheel in `build/wheels/main-<short-commit>/`, verifies a draft's
uploaded bytes and full commit tag, then publishes the private prerelease
`main-<full-commit>` with `latest=false`. Retries reuse the retained wheel and
resume missing assets in the same draft. An interruption before an upload can
resume; a failed transfer may leave an incomplete `starter` asset. Incomplete,
unequal or unreadable existing assets stop for owner diagnosis; the tool never
deletes or replaces them. A different identity also stops without replacement.
`build/wheels/main-current/` is updated after confirmed publication. The
[downstream guide](downstream/README.md) owns authenticated download and file pins.
The documentation workflow builds the site and its examples; it does not run
the scientific delivery gates. `bash tools/check-all-container.sh` runs
`check-all.sh` on HEAD from a cold build in a digest-pinned `ubuntu:24.04`
container, the system a hosted runner would use, at no GitHub minutes.

Build the documentation with its optional documentation dependencies and the
checkout's editable Engine:

```bash
nice -n 10 bash tools/build-docs.sh
```

The command reuses `tools/build-engine.sh`, forces every example notebook to
execute in this checkout's `.venv`, and runs Sphinx HTML and cheap doctests with
warnings treated as errors. It selects MUMPS with `EPCSAFT_HSL_LIBRARY` and
`EPCSAFT_IPOPT_OPTIONS_FILE` unset. Outputs, logs and timing are under
`build/docs/`; open `build/docs/html/index.html`. The gallery is initially an
index; the notebook part of #695 supplies its eleven examples. No notebook
execution cache substitutes for this gate.

For process examples, the root environment already pins IDAES and Pyomo. Obtain
the optional IDAES extensions 3.4.2 solver/library archives directly from
[upstream](https://github.com/IDAES/idaes-ext/releases/tag/3.4.2).
The documentation workflow checks their hashes and extracts them under the
runner's temporary directory. Local process-example builds can set `IDAES_DATA`
and `PYOMO_CONFIG_DIR` to an equivalent directory outside the checkout with the
tools in its `bin/` subdirectory. The upstream Ipopt includes HSL, even when
MUMPS is selected: retain upstream notices and keep executables, libraries,
archives and their environment outside the repository, Engine builds, docs
outputs, uploads and shared caches. Only rendered pages, original notebooks,
figures, numerical records and text logs are uploaded. See the
[acknowledgement](docs/site/contributing.md#optional-process-example-tools).

Pages deployment requires a public repository, main and
`DOCS_PAGES_ENABLED=true`; leave that variable unset until launch. The build
and selected record uploads run while deployment is disabled.

Run a canonical end-to-end workflow with:

```bash
bash tools/run-slice.sh khudaida
```

See `AGENTS.md`, `ARCHITECTURE.yaml`, and
`docs/scientific/CONTEXT.md` for execution, architecture, and scientific
ownership.
