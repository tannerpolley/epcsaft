# Installation

Use a repository checkout for the examples: their source parameters and
observation files live under `data/` and `analyses/`, outside the Python wheel.
Run the commands below from the repository root. The Engine performs the
thermodynamic calculations in C++; Python declares inputs and reads results.
The [water-state notebook](examples/parameters-and-water-state.ipynb) checks
this installation with a native pressure closure.

## Source installation

The supported development environment is Linux with Git, Python 3.13, `uv`,
a C and C++ compiler, Ninja, CMake, pkg-config and `flock`. Install the native
development libraries before setup: Ipopt 3.11.9, Eigen 3.4 or newer and
Ceres 2.2. The locked Python environment supplies Meson, meson-python,
nanobind and `patchelf`. The first build downloads the pinned CppAD and HiGHS
archives, so it needs network access.

```bash
git clone https://github.com/tannerpolley/ePC-SAFT.git
cd ePC-SAFT
nice -n 10 bash .codex/environments/setup.sh
.venv/bin/python -c 'import epcsaft; print(epcsaft.__file__)'
```

Setup creates the checkout's `.venv` and installs the Engine editable. Activate
it with `source .venv/bin/activate` before interactive Python or Jupyter work.
An import rebuilds changed native sources through Ninja; Python edits need no
rebuild. Repeat setup after changing the lock file. The
[repository setup instructions](https://github.com/tannerpolley/ePC-SAFT/blob/main/README.md#setup-and-installation)
own these commands and describe rebuilds after toolchain changes.

## Installing a wheel supplied for a calculation

If you have an identified Engine wheel, install that exact file into the
calculation's environment. Record its origin with your inputs. There is no
public package-index installation promised here. A wheel supplies the Engine,
not the repository's example data. Use the source checkout above for the
complete downloadable notebooks.

The repository's `tools/build-engine.sh --wheel <directory>` builds a wheel;
`--install-wheel` replaces the editable installation for a wheel-based run.
Run setup afterwards to restore the editable environment. See the
[Engine installation reference](https://github.com/tannerpolley/ePC-SAFT/blob/main/engine/README.md#build-and-installation)
for the installation boundary.

## Optional process calculations

Ordinary state and equilibrium calculations need neither IDAES nor Pyomo.
The repository environment includes the pinned Python packages used by the
[process-model examples](guide/process-models.md). The
[methane ASL notebook](examples/methane-asl-properties.ipynb) calls the installed
Engine property library without an external Ipopt executable. The
[IDAES settler notebook](examples/idaes-two-liquid-settler.ipynb) additionally
needs upstream IDAES solver tools. Obtain them directly from upstream and
follow the [external-tool boundary and acknowledgement](contributing.md#optional-process-example-tools).

## Check the calculation you installed

Continue with the [quickstart](quickstart.md). Check units, source selection,
branch and achieved residuals before interpreting a value. The
[parameter rules](theory/algorithms.md#alg-parameters) and
[density algorithm](theory/algorithms.md#alg-density) define these operations;
the [Python model reference](api/eos.rst) gives the call signatures.
