#!/usr/bin/env python3
"""Collect a binary methane/ethane bubble point and warm runtime statistics."""

from __future__ import annotations

import csv
import json
import math
import os
import statistics
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
ANALYSIS = Path(__file__).resolve().parents[1]
RESULTS = ANALYSIS / "results"
FORKS = Path("/home/tnnrpolley21/Workspaces/Engineering/Forks-Clones")
TEMPERATURE_K = 180.0
FEED = (0.5, 0.5)
COMPONENTS = ("methane", "ethane")
WARMUPS = int(os.environ.get("RUNTIME_WARMUPS", "5"))
SAMPLES = int(os.environ.get("RUNTIME_SAMPLES", "25"))


def git_info(path: Path) -> dict[str, object]:
    def git(*args: str) -> str:
        return subprocess.run(
            ["git", "-C", str(path), *args],
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip()

    return {"path": str(path), "commit": git("rev-parse", "HEAD"), "dirty": bool(git("status", "--porcelain"))}


def scalar(value: object) -> float:
    try:
        if len(value) == 1:  # type: ignore[arg-type]
            value = value[0]  # type: ignore[index]
    except TypeError:
        pass
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"non-finite result {value!r}")
    return result


def vector(value: object) -> list[float]:
    return [scalar(item) for item in value]  # type: ignore[union-attr]


def timed(call):  # type: ignore[no-untyped-def]
    for _ in range(WARMUPS):
        call()
    samples_ms: list[float] = []
    result = None
    for _ in range(SAMPLES):
        started = time.perf_counter_ns()
        result = call()
        samples_ms.append((time.perf_counter_ns() - started) / 1.0e6)
    samples_sorted = sorted(samples_ms)
    return result, {
        "warmups": WARMUPS,
        "samples": SAMPLES,
        "solve_mean_ms": statistics.fmean(samples_ms),
        "solve_median_ms": statistics.median(samples_ms),
        "solve_p10_ms": float(np.percentile(samples_sorted, 10)),
        "solve_p90_ms": float(np.percentile(samples_sorted, 90)),
        "solve_min_ms": samples_sorted[0],
        "solve_max_ms": samples_sorted[-1],
        "solve_samples_ms": samples_ms,
    }


def result_row(
    name: str,
    source: Path,
    solver: str,
    setup_ms: float,
    runtime: dict[str, object],
    values: dict[str, object],
) -> dict[str, object]:
    return {
        "implementation": name,
        "solver": solver,
        "temperature_K": TEMPERATURE_K,
        "feed_methane": FEED[0],
        "feed_ethane": FEED[1],
        "pressure_Pa": scalar(values["pressure_Pa"]),
        "liquid_x_methane": scalar(values["liquid_x_methane"]),
        "liquid_x_ethane": scalar(values["liquid_x_ethane"]),
        "vapor_y_methane": scalar(values["vapor_y_methane"]),
        "vapor_y_ethane": scalar(values["vapor_y_ethane"]),
        "rho_liquid_mol_m3": scalar(values["rho_liquid_mol_m3"]),
        "rho_vapor_mol_m3": scalar(values["rho_vapor_mol_m3"]),
        "phase_pressure_rel_error": values.get("phase_pressure_rel_error"),
        "setup_ms": setup_ms,
        **runtime,
        "source": git_info(source),
    }


def collect_engine(*, declare_only=False):
    from epcsaft import Mixture, Parameters
    from epcsaft.equilibrium import MoleFractions, Phase, bubble_point, solve_equilibrium

    started = time.perf_counter_ns()
    packet = ROOT / "data/packets/gross-2001-methane-ethane/3/parameters/parameters.json"
    model = Mixture(Parameters.from_json(packet, components=list(COMPONENTS)))
    problem = bubble_point(
        liquid=Phase("L", kind="liquid", composition_guess=FEED),
        T=TEMPERATURE_K,
        feed=MoleFractions(dict(zip(COMPONENTS, FEED))),
    )
    if declare_only:
        return model, problem
    setup_ms = (time.perf_counter_ns() - started) / 1.0e6
    final, runtime = timed(lambda: solve_equilibrium(model, problem))
    assert final.success, final.convergence
    mole_fractions = vector(final.mole_fractions)
    densities = vector(final.molar_densities)
    residual = max(abs(scalar(value)) for value in final.residuals)
    return result_row(
        "ePC-SAFT",
        ROOT,
        f"Engine bubble_point + solve_equilibrium (max_abs_residual={residual:.3e})",
        setup_ms,
        runtime,
        {
            "pressure_Pa": final.pressure,
            "liquid_x_methane": mole_fractions[0],
            "liquid_x_ethane": mole_fractions[1],
            "vapor_y_methane": mole_fractions[2],
            "vapor_y_ethane": mole_fractions[3],
            "rho_liquid_mol_m3": densities[0],
            "rho_vapor_mol_m3": densities[1],
            "phase_pressure_rel_error": None,
        },
    )


def collect_external() -> list[dict[str, object]]:
    import si_units as si
    from feos import EquationOfState, Parameters, PhaseEquilibrium
    from thermopack import pcsaft

    feos_root = FORKS / "feos"
    started = time.perf_counter_ns()
    params = Parameters.from_json(list(COMPONENTS), str(feos_root / "parameters/pcsaft/gross2001.json"))
    feos_model = EquationOfState.pcsaft(params)
    feos_x = np.asarray(FEED)
    feos_setup_ms = (time.perf_counter_ns() - started) / 1.0e6
    feos_final, feos_runtime = timed(
        lambda: PhaseEquilibrium.bubble_point(feos_model, TEMPERATURE_K * si.KELVIN, feos_x)
    )
    feos_liquid = feos_final.liquid
    feos_vapor = feos_final.vapor
    feos_pressure_l = scalar(feos_liquid.pressure() / si.PASCAL)
    feos_pressure_v = scalar(feos_vapor.pressure() / si.PASCAL)
    density_unit = si.MOL / si.METER**3
    feos_row = result_row(
        "FeOs",
        feos_root,
        "PhaseEquilibrium.bubble_point",
        feos_setup_ms,
        feos_runtime,
        {
            "pressure_Pa": feos_pressure_l,
            "liquid_x_methane": feos_liquid.molefracs[0],
            "liquid_x_ethane": feos_liquid.molefracs[1],
            "vapor_y_methane": feos_vapor.molefracs[0],
            "vapor_y_ethane": feos_vapor.molefracs[1],
            "rho_liquid_mol_m3": feos_liquid.density / density_unit,
            "rho_vapor_mol_m3": feos_vapor.density / density_unit,
            "phase_pressure_rel_error": abs(feos_pressure_l - feos_pressure_v) / feos_pressure_l,
        },
    )

    started = time.perf_counter_ns()
    thermo_model = pcsaft.pcsaft("C1,C2")
    thermo_setup_ms = (time.perf_counter_ns() - started) / 1.0e6
    thermo_final, thermo_runtime = timed(
        lambda: thermo_model.bubble_pressure(TEMPERATURE_K, np.asarray(FEED))
    )
    thermo_pressure, thermo_y = thermo_final
    thermo_liquid_volume = thermo_model.specific_volume(
        TEMPERATURE_K, thermo_pressure, np.asarray(FEED), thermo_model.LIQPH
    )
    thermo_vapor_volume = thermo_model.specific_volume(
        TEMPERATURE_K, thermo_pressure, thermo_y, thermo_model.VAPPH
    )
    thermo_row = result_row(
        "ThermoPack",
        FORKS / "thermopack",
        "pcsaft.bubble_pressure",
        thermo_setup_ms,
        thermo_runtime,
        {
            "pressure_Pa": thermo_pressure,
            "liquid_x_methane": FEED[0],
            "liquid_x_ethane": FEED[1],
            "vapor_y_methane": thermo_y[0],
            "vapor_y_ethane": thermo_y[1],
            "rho_liquid_mol_m3": 1.0 / thermo_liquid_volume,
            "rho_vapor_mol_m3": 1.0 / thermo_vapor_volume,
            "phase_pressure_rel_error": None,
        },
    )
    return [feos_row, thermo_row]


def run_mode(mode: str) -> object:
    if mode == "engine":
        return collect_engine()
    if mode == "external":
        return collect_external()
    raise ValueError(mode)


def main() -> None:
    if len(sys.argv) > 1:
        print(json.dumps(run_mode(sys.argv[1])))
        return

    comparison_python = os.environ.get("COMPARISON_PYTHON", sys.executable)
    commands = [
        [str(ROOT / ".venv/bin/python"), str(Path(__file__)), "engine"],
        [comparison_python, str(Path(__file__)), "external"],
    ]
    records: list[dict[str, object]] = []
    command_log: list[list[str]] = []
    for command in commands:
        environment = os.environ.copy()
        completed = subprocess.run(command, check=True, capture_output=True, text=True, env=environment)
        value = json.loads(completed.stdout.strip().splitlines()[-1])
        records.extend(value if isinstance(value, list) else [value])
        command_log.append(command)

    RESULTS.mkdir(parents=True, exist_ok=True)
    fields = [
        "implementation", "solver", "temperature_K", "feed_methane", "feed_ethane",
        "pressure_Pa", "liquid_x_methane", "liquid_x_ethane", "vapor_y_methane",
        "vapor_y_ethane", "rho_liquid_mol_m3", "rho_vapor_mol_m3",
        "phase_pressure_rel_error", "setup_ms", "warmups", "samples", "solve_mean_ms",
        "solve_median_ms", "solve_p10_ms", "solve_p90_ms", "solve_min_ms", "solve_max_ms",
    ]
    with (RESULTS / "values.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows({key: value for key, value in record.items() if key in fields} for record in records)
    (RESULTS / "provenance.json").write_text(
        json.dumps(
            {
                "case": {
                    "components": COMPONENTS,
                    "temperature_K": TEMPERATURE_K,
                    "feed_mole_fractions": FEED,
                    "parameter_source": "Gross-Sadowski 2001, local methane/ethane packet, k_ij=0",
                    "runtime_protocol": {"warmups": WARMUPS, "samples": SAMPLES, "clock": "perf_counter_ns", "FEOS_MAX_THREADS": os.environ.get("FEOS_MAX_THREADS"), "OMP_NUM_THREADS": os.environ.get("OMP_NUM_THREADS")},
                },
                "commands": command_log,
                "records": records,
            },
            indent=2,
        )
        + "\n"
    )
    print(f"wrote {RESULTS / 'values.csv'}")


def declarations():
    return collect_engine(declare_only=True)


if __name__ == "__main__":
    main()
