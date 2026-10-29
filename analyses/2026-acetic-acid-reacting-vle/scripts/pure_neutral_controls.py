"""Recheck pure and neutral saturation controls with this checkout's installed Engine."""

import hashlib
import json
from pathlib import Path
import epcsaft
from epcsaft import equilibrium as eq

leaf = Path(__file__).resolve().parents[1]
root = leaf.parents[1]
doc = json.loads((leaf / "inputs/parameters.json").read_text())
T = 323.15
assert root in Path(epcsaft._core.__file__).resolve().parents, (
    "Use this checkout's prepared .venv"
)
out = {
    "loaded_package": epcsaft.__file__,
    "loaded_extension": epcsaft._core.__file__,
    "solves": [],
}
for name in ("water", "acetic-acid"):
    mix = epcsaft.Pure(epcsaft.Parameters.from_mapping(doc, components=[name]))
    r = eq.solve_equilibrium(mix, mix.saturation(T))
    assert r.success and r.requested_tolerance_met
    out["solves"].append(
        dict(
            name=name,
            success=r.success,
            pressure=r.pressure,
            row_residual=r.convergence.residual,
            densities=list(r.molar_densities),
            termination=str(r.solver_status),
        )
    )
mix = epcsaft.Mixture(
    epcsaft.Parameters.from_mapping(doc, components=["water", "acetic-acid"])
)
for xw in (0.5781, 0.6673, 0.7683):
    q = eq.Problem(
        T=T,
        P=eq.Free(11000),
        feed=eq.Amounts({"water": xw, "acetic-acid": 1 - xw}),
        phases=[
            eq.Phase("L", kind="liquid"),
            eq.Phase("V", kind="vapor", amount=eq.Pinned(0)),
        ],
    )
    r = eq.solve_equilibrium(mix, q)
    assert r.success and r.requested_tolerance_met
    out["solves"].append(
        dict(
            x_water=xw,
            success=r.success,
            pressure=r.pressure,
            row_residual=r.convergence.residual,
            densities=list(r.molar_densities),
            termination=str(r.solver_status),
        )
    )
paths = sorted(
    {
        line.split()[-1]
        for line in Path("/proc/self/maps").read_text().splitlines()
        if "epcsaft" in line and ".so" in line
    }
)
out["loaded_binary_sha256"] = {
    p: hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in paths
}
(leaf / "results/pressure-diagnosis/pure-neutral-controls.json").write_text(
    json.dumps(out, indent=2) + "\n"
)
print(json.dumps(out, indent=2))
