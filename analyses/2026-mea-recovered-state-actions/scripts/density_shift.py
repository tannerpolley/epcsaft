"""#331 falsifier: the phase-pressure row of the 0131 -> 0188 feed continuation under liquid-density shifts.

Usage: density_shift.py MEA_CHECKOUT (the sweep inputs, 04a9328). The solved liquid density is shifted by
+-1e-11 and +-1e-12 relative and the state is rechecked with run.py's criteria; the metric predicts
|row| ~ |shift| rho dP/drho / (rho R T), about 3e-10 (fails ROW_TOL) and 3e-11 (passes).
"""
from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

sys.argv = [sys.argv[0], sys.argv[1], "density-shift"]
spec = importlib.util.spec_from_file_location("run", Path(__file__).resolve().parent / "run.py")
run = importlib.util.module_from_spec(spec)
spec.loader.exec_module(run)


class Shifted:
    """The solved result with its liquid density replaced; every other value is the Engine's."""

    def __init__(self, result, densities):
        self._result, self.molar_densities = result, densities

    def __getattr__(self, name):
        return getattr(self._result, name)


ids = ("vle_obs_0131", "vle_obs_0188")
run.REQUESTS.update({i: run.request(i) for i in ids})
outcome = run.continuation(run.cold(run.REQUESTS[ids[0]]), run.REQUESTS[ids[1]])
result, problem, req = outcome["result"], outcome["problem"], run.REQUESTS[ids[1]]
liquid = run.vt.phase_index(problem, "liquid")
state = run.MODEL.state(result.temperature, rho=result.molar_densities[liquid],
                        x=list(run.np.asarray(result.mole_fractions).reshape(result.phase_count, run.WIDTH)[liquid]))
rows = []
for shift in (0.0, 1e-11, -1e-11, 1e-12, -1e-12):
    densities = list(result.molar_densities)
    densities[liquid] *= 1.0 + shift
    out = run.criteria(problem, Shifted(result, densities), req)
    predicted = abs(shift) * state.pressure_density_derivative / (run.R * result.temperature)
    rows.append({"relative_density_shift": shift, "reduced_pressure_residual_liquid": out["reduced_pressure_residual_liquid"],
                 "predicted_shift_contribution": predicted, "fugacity_max": out["fugacity_max"],
                 "criteria_ok": out["criteria_ok"]})
    print(rows[-1], flush=True)
record = {"transition": "vle_obs_0131 -> vle_obs_0188, frozen continuation", "ok": outcome["ok"], "work": outcome["work"],
          "pressure_pa": result.pressure, "liquid_molar_density": result.molar_densities[liquid],
          "row_tolerance": run.ROW_TOL, "rows": rows,
          "provenance": {"wheel": run.WHEEL_NAME, "wheel_sha256": run.shared.ENGINE_WHEEL_SHA256,
                         "mea_commit": run.subprocess.check_output(["git", "-C", str(run.MEA), "rev-parse", "HEAD"], text=True).strip()}}
(run.RESULTS / "issue-331-density-shift.json").write_text(json.dumps(record, indent=1) + "\n", encoding="utf-8")
