#!/usr/bin/env bash
# Replay the application checks invalidated by Engine changes since each pin, on the installed
# final wheel, one step at a time and single-threaded (#64).
# Usage: replay.sh MEA_CHECKOUT LITHIUM_CHECKOUT SCRATCH [STEP...]
#   MEA_CHECKOUT      detached MEA-Thermodynamics worktree at 799b201 (record 868a5018)
#   LITHIUM_CHECKOUT  detached Lithium_Extraction worktree at 53f1714 whose .venv has the final wheel
#   SCRATCH           output directory; compare.py reads it
# Steps: mea-actions mea-sweep mea-heat mea-scores li-fit li-holdout li-volume li-live li-v3 (default: all).
# Each step waits for the 1-minute load average to fall below MAX_LOAD (default 5).
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYOMO_CONFIG_DIR="$HOME/.idaes"
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
engine="$(cd "$here/../../.." && pwd)"
mea="$(realpath "$1")" li="$(realpath "$2")" out="$(realpath -m "$3")"
shift 3
steps=("$@")
((${#steps[@]})) || steps=(mea-actions mea-sweep mea-heat mea-scores li-fit li-holdout li-volume li-live li-v3)
mkdir -p "$out"

installed() {  # the wheel file an environment's epcsaft was installed from
    "$1" -c 'import json,importlib.metadata as m;from urllib.parse import unquote,urlparse;print(unquote(urlparse(json.loads(m.distribution("epcsaft").read_text("direct_url.json"))["url"]).path))'
}
li_py="$li/.venv/bin/python"
wheel="$(installed "$li_py")"
sha="$(sha256sum "$wheel" | cut -d ' ' -f 1)"
[[ "$(sha256sum "$(installed "$engine/.venv/bin/python")" | cut -d ' ' -f 1)" == "$sha" ]] \
    || { echo "Engine and Lithium environments hold different wheels" >&2; exit 1; }
commit="$(git -C "$engine" rev-parse HEAD)"
flags=(--engine-wheel "$wheel" --engine-sha256 "$sha" --engine-commit "$commit")
printf 'wheel=%s\nsha256=%s\nengine_commit=%s\n' "$wheel" "$sha" "$commit" >"$out/identity.txt"

wait_load() {
    while awk -v max="${MAX_LOAD:-5}" '{exit !($1 >= max)}' /proc/loadavg; do sleep 30; done
}

for step in "${steps[@]}"; do
    wait_load
    echo "[$(date -Is)] start $step (load $(cut -d ' ' -f 1 /proc/loadavg))"
    start=$SECONDS
    case "$step" in
        mea-actions | mea-sweep | mea-heat)
            # run.py writes beside itself; a scratch copy keeps the #61 record intact.
            mkdir -p "$out/mea-run/scripts"
            cp "$engine/analyses/2026-mea-recovered-state-actions/scripts/run.py" "$out/mea-run/scripts/"
            "$engine/.venv/bin/python" "$out/mea-run/scripts/run.py" "$mea" "${step#mea-}" ;;
        mea-scores)
            "$engine/.venv/bin/python" "$here/mea_scores.py" "$mea" ;;
        li-fit)
            (cd "$li" && "$li_py" -c 'import json, sys
from analyses.hbta_topo_calibrated_surrogate.scripts.generate import _fit_diluent
json.dump(_fit_diluent("dodecane"), open(sys.argv[1], "w"), indent=1, default=str)' "$out/li-fit.json") ;;
        li-holdout)
            (cd "$li" && "$li_py" -m analyses.hbta_topo_configured_reactive_lle.scripts.zhang_isotherm_holdout "${flags[@]}") ;;
        li-volume)
            # Nonzero exit when the observed 3e5 limit is rejected; compare.py applies the required scales.
            (cd "$li" && "$li_py" -m analyses.hbta_topo_configured_reactive_lle.scripts.verify_volume_scaling "${flags[@]}") || true ;;
        li-live)
            (cd "$li" && "$li_py" -m analyses.hbta_topo_ms2_idaes_optimization.scripts.run_live_adapter \
                --output "$out/li-live" "${flags[@]}") ;;
        li-v3)
            # compare.py pins each variant's wheel; the disposable checkout names the final wheel instead.
            sed -i -e "s/18c7812feff19751f10ee3d60991d782a44a0fe7b4baceccd2a8ae9ef4a1bfbc/$sha/" \
                -e "s/\"engine_commit\": \"1303c119e4a21ba31e46596fe62ee3fdfe4cc253\"/\"engine_commit\": \"$commit\"/" \
                "$li/analyses/hbta_topo_ms2_engine_comparison/scripts/compare.py"
            (cd "$li" && PYTHONPATH="$li" "$li_py" analyses/hbta_topo_ms2_engine_comparison/scripts/compare.py generate \
                --variant greenfield_current --wheel "$wheel" --method v3 --protocol first_strict --workers 1 \
                --output "$out/li-v3") ;;
        *) echo "unknown step $step" >&2; exit 2 ;;
    esac
    echo "[$(date -Is)] done $step in $((SECONDS - start)) s"
done
