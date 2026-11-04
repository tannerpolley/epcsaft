# Installation

The Engine performs the thermodynamic calculations in C++; Python declares
inputs and reads results. The examples read parameter packets and observation
files under `data/` and `analyses/` of the source release. Run the examples
from the root of the unpacked source archive of the same
[release](https://github.com/tannerpolley/epcsaft/releases).
The [water-state notebook](examples/parameters-and-water-state.ipynb) checks
an installation with a native pressure closure.

## Conda package

Install release 1.3.0 for Linux (linux-64) and Windows (win-64) from the
`tannerwpolley` conda channel:

```bash
conda create -n epcsaft -c conda-forge -c tannerwpolley python=3.13 epcsaft=1.3.0
conda activate epcsaft
python -c 'import epcsaft; print(epcsaft.__file__)'
```

The package contains the Python interface, the native extension and its
property library, a private Ipopt build with MUMPS, and the third-party
notices. The conda package also includes the FeOs Gross–Sadowski 2001 and
2002 packets. Wheels do not include the packet tree. For other packets,
download the matching source archive and pass its `data/packets` directory to
`Parameters.from_packets`; see [packet selection](guide/parameters-and-models.md#installed-and-downloaded-packets).
Each packet states its own use terms.

## Build from source

Building from source is supported on Linux with Python 3.13, `uv`, a C and C++
compiler, Ninja, CMake and pkg-config. Install the native development
libraries first: Ipopt 3.11.9, Eigen 3.4 and Ceres 2.2. The locked Python
environment supplies Meson, meson-python, nanobind and `patchelf`. The first
build downloads the checksum-pinned CppAD and HiGHS archives, so it needs
network access. From the root of the unpacked source archive:

```bash
uv sync --python 3.13 --locked --all-packages --no-install-workspace
uv pip install --python .venv/bin/python --no-build-isolation ./engine
.venv/bin/python -c 'import epcsaft; print(epcsaft.__file__)'
```

Activate the environment with `source .venv/bin/activate` before interactive
Python or Jupyter work. For an editable development installation, which
rebuilds changed native sources on import, follow the
[Engine build instructions](https://github.com/tannerpolley/epcsaft/blob/main/engine/README.md#build-and-installation).

## Optional process calculations

Ordinary state and equilibrium calculations need neither IDAES nor Pyomo.
The locked source environment includes the pinned Python packages used by the
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
