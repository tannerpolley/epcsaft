#!/usr/bin/env bash
set -euo pipefail
verification_dir=analyses/reproduction/2015-uyan/scripts/species-reference
results_dir=analyses/reproduction/2015-uyan/results/data/native/species-reference
mkdir -p build/issue-274
.venv/bin/python - <<'PY'
from pathlib import Path
import sys
sys.path.insert(0, 'engine')
from epcsaft import Pure
from epcsaft.equilibrium import solve_equilibrium
from tests.native_inputs import ETHANE, parameters
model = Pure(parameters((ETHANE,)))
state = solve_equilibrium(model, model.saturation(180.0))
assert state.success, state.message
Path('build/issue-274/ethane-saturation-pressure.txt').write_text(str(state.pressure))
PY
.venv/bin/python "$verification_dir/compile_native.py" "$verification_dir/reference_qualification.cpp"
build/issue-274/reference_qualification "$results_dir/qualification.tsv" "$(cat build/issue-274/ethane-saturation-pressure.txt)"
.venv/bin/python "$verification_dir/compile_native.py" "$verification_dir/reference_snapshot.cpp" -DMIXED_REFERENCE
build/issue-274/reference_snapshot > build/issue-274/candidate-native-values.csv
cmp "$results_dir/default-native-values.csv" build/issue-274/candidate-native-values.csv
.venv/bin/python "$verification_dir/wheel_snapshot.py"
