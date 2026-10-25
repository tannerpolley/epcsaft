#!/usr/bin/env bash
# Build the site with the checkout's installed Engine and fresh notebook execution.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
unset EPCSAFT_HSL_LIBRARY EPCSAFT_IPOPT_OPTIONS_FILE
export EPCSAFT_IPOPT_LINEAR_SOLVER=mumps MPLBACKEND=Agg
started="$EPOCHREALTIME"
mkdir -p build/docs
rm -f build/docs/timing.json
rm -rf build/docs/html build/docs/jupyter_execute
uv sync --python 3.13 --group docs --all-packages --no-install-workspace --locked --inexact
bash tools/build-engine.sh 2>&1 | tee build/docs/native-build.log
.venv/bin/python - <<'PY'
import epcsaft
from epcsaft import _core
assert _core.linear_solver() == "mumps", _core.linear_solver()
print(f"Documentation Engine: {epcsaft.__file__}; linear solver: {_core.linear_solver()}")
PY
# When the process-example tools are supplied, verify them before any notebook.
if [[ -n "${IDAES_DATA:-}" ]]; then
    .venv/bin/python - <<'PY'
import os
from pathlib import Path
import idaes
import pyomo.environ as pyo
from pyomo.contrib.pynumero.asl import AmplInterface
root = Path(os.environ['IDAES_DATA']).resolve()
assert not root.is_relative_to(Path.cwd())
assert Path(idaes.bin_directory).resolve() == root / 'bin'
assert AmplInterface.available()
model = pyo.ConcreteModel()
model.x = pyo.Var(initialize=2, bounds=(0, None))
model.objective = pyo.Objective(expr=(model.x - 1)**2)
solver = pyo.SolverFactory('ipopt', executable=str(root / 'bin/ipopt'))
result = solver.solve(model, options={'linear_solver': 'mumps', 'print_user_options': 'yes'},
                      logfile=str(root / 'outer-mumps.log'))
log = (root / 'outer-mumps.log').read_text()
assert 'running with linear solver mumps.' in log, log
assert abs(pyo.value(model.x) - 1) < 1e-6
(Path('build/docs') / 'outer-mumps.log').write_text(log)
print(result.solver.termination_condition)
PY
fi
# Register this checkout's interpreter; a system python3 kernel can load a different Engine.
export JUPYTER_PATH="$PWD/build/docs/jupyter/share/jupyter${JUPYTER_PATH:+:$JUPYTER_PATH}"
.venv/bin/python -m ipykernel install --prefix "$PWD/build/docs/jupyter" --name python3 --display-name 'ePC-SAFT docs'
html_started="$EPOCHREALTIME"
.venv/bin/python -m sphinx -E -W --keep-going -b html \
    -D nb_execution_mode=force -d build/docs/doctrees -w build/docs/html-warnings.log \
    docs/site build/docs/html 2>&1 | tee build/docs/html-build.log
html_finished="$EPOCHREALTIME"
.venv/bin/python -m sphinx -E -W --keep-going -b doctest \
    -D nb_execution_mode=off -d build/docs/doctest-doctrees -w build/docs/doctest-warnings.log \
    docs/site build/docs/doctest 2>&1 | tee build/docs/doctest-build.log
.venv/bin/python - "$started" "$html_started" "$html_finished" "$EPOCHREALTIME" <<'PY'
import json, subprocess, sys
from pathlib import Path
start, html_start, html_end, end = map(float, sys.argv[1:])
record = {"commit": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
          "build_seconds": end - start, "html_seconds": html_end - html_start,
          "doctest_seconds": end - html_end, "sphinx_warnings": 0,
          "notebooks": len(list(Path("docs/site/examples").glob("*.ipynb"))),
          "notebook_execution": "force", "engine_linear_solver": "mumps"}
Path("build/docs/timing.json").write_text(json.dumps(record, indent=2) + "\n")
print(json.dumps(record, indent=2))
PY
