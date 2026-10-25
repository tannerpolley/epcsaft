#!/usr/bin/env bash
# #572: run from the repository root. Builds the native campaign driver, runs the frozen comparison (existing
# inputs/draws.csv and inputs/bootstrap-indices.csv are read, never regenerated), records input hashes and
# retains the joint analytic reference comparisons of the native measurement-correlation case and the multistart table.
set -euo pipefail
leaf=analyses/2026-correlated-measurement-errors
driver=$(nice -n 10 tools/native-driver regression_checks "$leaf/scripts/campaign.cpp")
export LD_LIBRARY_PATH="$(<build/environment-wheel/meson-build-dir)"
"$driver" "$leaf" prepare
(cd "$leaf/inputs" && if [[ -f hashes.sha256 ]]; then sha256sum -c hashes.sha256; else sha256sum draws.csv bootstrap-indices.csv > hashes.sha256; fi)
MALLOC_PERTURB_=0 nice -n 10 "$driver" "$leaf" | tee "$leaf/results/campaign.log"
nice -n 10 tools/native-run regression_checks measurement-correlation \
    | sed -n 's/^EPCSAFT_COMPARE_JSON //p' \
    | jq -r '[.comparison_id, .quantity, .unit, .model, .reference, .atol, .verdict] | @csv' \
    | sed '1i comparison,quantity,unit,model,reference,allowance,verdict' > "$leaf/results/joint-reference.csv"
.venv/bin/python "$leaf/scripts/start_evidence.py"
