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

The supported development environment is Linux with Python 3.13, `uv`, a
C and C++ compiler, Ninja, CMake, pkg-config, and `flock`. The locked root
environment supplies Meson, meson-python, nanobind, and `patchelf`. The native
graph also needs the environment-provided Ipopt, Eigen and Ceres development
packages, at the versions in [engine/README.md](engine/README.md#native-dependency-graph). The first build needs network access for the two
checksum-pinned Meson subproject archives (CppAD and HiGHS). `ccache` is optional.

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
.venv/bin/python -c 'import epcsaft; from importlib.metadata import version; print(version("epcsaft")); print(epcsaft._core.__file__)'
```

The code uses GPLv3 with an additional permission under section 7 for linking
with EPL-licensed Ipopt; see [LICENSE](LICENSE). Retained third-party data have
separate terms, including non-commercial restrictions; see [NOTICE](NOTICE).

## Setup and installation

For an editable development installation, run from the extracted source root after
installing the prerequisites above:

```bash
uv sync --python 3.13 --locked --all-packages --no-install-workspace
nice -n 10 bash tools/build-engine.sh
.venv/bin/python -c 'import epcsaft; from importlib.metadata import version; print(version("epcsaft")); print(epcsaft._core.__file__)'
```

`tools/build-engine.sh` installs `engine/` editable through meson-python into
`build/engine-meson/shared`. Imports rebuild changed native sources with Ninja;
Python edits need no rebuild. Repeat the commands after a lock change.

```python
import epcsaft
from epcsaft import equilibrium
```
