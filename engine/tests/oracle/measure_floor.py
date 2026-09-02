"""Measure the installed Python EOS entry points; timings are reported only."""
from __future__ import annotations

import json
import os
import subprocess
import sys
import time
from pathlib import Path
from statistics import median
from typing import Callable

assert os.environ.get("OMP_NUM_THREADS") == "1"

from epcsaft import Mixture, Parameters

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from tests.packets import parameters_json

FLOOR_COMPOSITION = (0.80, 0.10, 0.05, 0.02, 0.02, 0.005, 0.0, 0.005, 0.0)
COUNTER_NAMES = ("forward", "reverse", "hessian_calls", "recordings")


def counter_values(model: Mixture) -> dict[str, int]:
    counters = model.counters()
    return {name: int(getattr(counters, name)) for name in COUNTER_NAMES}


def measure(model: Mixture, operation: Callable[[], object]) -> dict[str, object]:
    for _ in range(3):
        operation()
    model.reset_counters()
    before = counter_values(model)
    times: list[int] = []
    for _ in range(20):
        start = time.perf_counter_ns()
        operation()
        times.append(time.perf_counter_ns() - start)
    after = counter_values(model)
    return {"wall_us": median(times) / 1000.0, "counters_per_call": {
        name: (after[name] - before[name]) / 20 for name in COUNTER_NAMES}}


def environment() -> dict[str, str]:
    commit = subprocess.run(
        ["git", "rev-parse", "HEAD"], check=True, capture_output=True, text=True
    ).stdout.strip()
    model_name = next(
        line.split(":", 1)[1].strip()
        for line in Path("/proc/cpuinfo").read_text().splitlines()
        if line.startswith("model name")
    )
    return {"git": commit, "cpu": model_name, "python": sys.version,
            "OMP_NUM_THREADS": os.environ["OMP_NUM_THREADS"]}


def main() -> None:
    methane_ethane = Mixture(Parameters.from_json(
        parameters_json("gross-2001-methane-ethane", 2),
        components=("methane", "ethane")))
    mea = Mixture(Parameters.from_json(
        parameters_json("mea-co2-h2o-nine-species-estimation", 2),
        components=("water", "monoethanolamine", "carbon-dioxide",
                    "protonated-monoethanolamine", "carbamate-anion",
                    "bicarbonate-anion", "carbonate-anion", "hydronium-cation",
                    "hydroxide-anion")))
    figiel = Mixture(Parameters.from_json(
        parameters_json("figiel-2025-reference-electrolytes", 2),
        components=("water", "sodium-cation", "chloride-anion")))
    case_c = lambda: mea.state(313.15, P=1e5, x=FLOOR_COMPOSITION, phase="liquid")
    state_c = case_c()
    reference = figiel.reference(
        298.15, 1e5, solvent_fractions=(1.0, 0.0, 0.0),
        ions=("sodium-cation", "chloride-anion"), stoichiometry=(1, 1))
    rows = [
        {"case": "A", **measure(methane_ethane, lambda: methane_ethane.state(300.0, rho=1e4, x=(0.4, 0.6)))},
        {"case": "B", **measure(mea, lambda: mea.state(313.15, rho=3e4, x=FLOOR_COMPOSITION))},
        {"case": "C", **measure(mea, case_c)},
        {"case": "D", **measure(mea, lambda: mea.state(313.15, P=1e5, x=FLOOR_COMPOSITION, phase="liquid", anchor=state_c))},
        {"case": "E", **measure(figiel, lambda: figiel.reference(298.15, 1e5, solvent_fractions=(1.0, 0.0, 0.0), ions=("sodium-cation", "chloride-anion"), stoichiometry=(1, 1)))},
        {"case": "F", **measure(figiel, lambda: reference.activity(1.0))},
    ]
    metadata = environment()
    for name, value in metadata.items():
        print(f"{name}={value}")
    for row in rows:
        print(json.dumps(row, sort_keys=True))
    if len(sys.argv) > 1:
        Path(sys.argv[1]).write_text(json.dumps({"environment": metadata, "rows": rows}, indent=2) + "\n")


if __name__ == "__main__":
    main()
