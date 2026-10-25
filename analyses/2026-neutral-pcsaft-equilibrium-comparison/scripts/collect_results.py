#!/usr/bin/env python3
"""Collect one neutral PC-SAFT saturation point from local package forks."""

from __future__ import annotations

import csv
import json
import math
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
ANALYSIS = Path(__file__).resolve().parents[1]
RESULTS = Path(os.environ.get("EPCSAFT_COMPARISON_OUTPUT", ANALYSIS / "results"))
FORKS = Path("/home/tnnrpolley21/Workspaces/Engineering/Forks-Clones")
T_K = 300.0
PARAMETERS = {
    "component": "propane",
    "temperature_K": T_K,
    "m": 2.0020,
    "sigma_Angstrom": 3.6184,
    "epsilon_over_k_K": 208.11,
    "association": False,
    "polarity": False,
    "electrolyte": False,
}


def git_info(path: Path) -> dict[str, object]:
    def git(*args: str) -> str:
        return subprocess.run(
            ["git", "-C", str(path), *args],
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip()

    return {
        "path": str(path),
        "commit": git("rev-parse", "HEAD"),
        "dirty": bool(git("status", "--porcelain")),
    }


def float_value(value: object) -> float:
    try:
        if len(value) == 1:  # type: ignore[arg-type]
            value = value[0]  # type: ignore[index]
    except TypeError:
        pass
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"non-finite result: {value!r}")
    return result


def row(name: str, values: dict[str, object], source: Path, solver: str) -> dict[str, object]:
    return {
        "implementation": name,
        "solver": solver,
        "temperature_K": T_K,
        "pressure_Pa": float_value(values["pressure_Pa"]),
        "rho_liquid_mol_m3": float_value(values["rho_liquid_mol_m3"]),
        "rho_vapor_mol_m3": float_value(values["rho_vapor_mol_m3"]),
        "phase_pressure_rel_error": values.get("phase_pressure_rel_error"),
        "source": git_info(source),
    }


def collect_engine(*, declare_only=False):
    from epcsaft import Parameters, Pure
    from epcsaft.equilibrium import solve_equilibrium

    packet = ROOT / "data/packets/gross-2001-propane/3/parameters/parameters.json"
    model = Pure(Parameters.from_json(packet, components=["propane"]))
    problem = model.saturation(T_K)
    if declare_only:
        return model, problem
    result = solve_equilibrium(model, problem)
    assert result.success, result.convergence
    densities = [float_value(x) for x in result.molar_densities]
    residual = max(abs(float_value(x)) for x in result.residuals)
    return row(
        "ePC-SAFT",
        {
            "pressure_Pa": result.pressure,
            "rho_liquid_mol_m3": densities[0],
            "rho_vapor_mol_m3": densities[1],
            "phase_pressure_rel_error": None,
        },
        ROOT,
        f"Engine saturation + solve_equilibrium (max_abs_residual={residual:.3e})",
    )


def collect_python_packages() -> list[dict[str, object]]:
    import si_units
    from feos import EquationOfState, Parameters, PhaseEquilibrium
    from thermopack import pcsaft

    feos_root = FORKS / "feos"
    feos_params = feos_root / "parameters/pcsaft/gross2001.json"
    feos_eos = EquationOfState.pcsaft(
        Parameters.from_json(["propane"], str(feos_params))
    )
    feos_equilibrium = PhaseEquilibrium.pure(feos_eos, T_K * si_units.KELVIN)
    feos_liquid = feos_equilibrium.liquid
    feos_vapor = feos_equilibrium.vapor
    feos_pressure_l = float_value(feos_liquid.pressure() / si_units.PASCAL)
    feos_pressure_v = float_value(feos_vapor.pressure() / si_units.PASCAL)
    density_unit = si_units.MOL / si_units.METER**3

    thermopack_model = pcsaft.pcsaft("C3")
    temperatures, pressures, volumes_l, volumes_v = (
        thermopack_model.get_pure_fluid_saturation_curve(
            initial_pressure=None, initial_temperature=T_K, nmax=100
        )
    )
    index = min(range(len(temperatures)), key=lambda i: abs(float(temperatures[i]) - T_K))
    if abs(float(temperatures[index]) - T_K) > 1.0e-10:
        raise ValueError(f"ThermoPack did not return T=300 K: {temperatures[index]}")

    return [
        row(
            "FeOs",
            {
                "pressure_Pa": feos_pressure_l,
                "rho_liquid_mol_m3": feos_liquid.density / density_unit,
                "rho_vapor_mol_m3": feos_vapor.density / density_unit,
                "phase_pressure_rel_error": abs(feos_pressure_l - feos_pressure_v)
                / feos_pressure_l,
            },
            feos_root,
            "PhaseEquilibrium.pure",
        ),
        row(
            "ThermoPack",
            {
                "pressure_Pa": pressures[index],
                "rho_liquid_mol_m3": 1.0 / volumes_l[index],
                "rho_vapor_mol_m3": 1.0 / volumes_v[index],
                "phase_pressure_rel_error": None,
            },
            FORKS / "thermopack",
            "pcsaft.get_pure_fluid_saturation_curve",
        ),
    ]


def collect_teqp() -> dict[str, object]:
    import numpy as np
    import teqp

    coeff = teqp.SAFTCoeffs()
    coeff.m = PARAMETERS["m"]
    coeff.sigma_Angstrom = PARAMETERS["sigma_Angstrom"]
    coeff.epsilon_over_k = PARAMETERS["epsilon_over_k_K"]
    teqp_model = teqp.PCSAFTEOS([coeff])
    rho_l, rho_v = teqp_model.pure_VLE_T(T_K, 11100.0, 482.0, 100)
    R = 8.31446261815324

    def teqp_pressure(rho: float) -> float:
        return float_value(
            rho * R * T_K * (1.0 + teqp_model.get_Ar01(T_K, rho, np.array([1.0])))
        )

    teqp_pressure_l = teqp_pressure(float(rho_l))
    teqp_pressure_v = teqp_pressure(float(rho_v))

    return row(
        "teqp",
        {
            "pressure_Pa": teqp_pressure_l,
            "rho_liquid_mol_m3": rho_l,
            "rho_vapor_mol_m3": rho_v,
            "phase_pressure_rel_error": abs(teqp_pressure_l - teqp_pressure_v)
            / teqp_pressure_l,
        },
        FORKS / "teqp",
        "PCSAFTEOS.pure_VLE_T",
    )


def collect_clapeyron() -> dict[str, object]:
    retained = os.environ.get("CLAPEYRON_RESULT")
    if retained:
        values = json.loads(Path(retained).read_text())
        return row("Clapeyron", values, FORKS / "Clapeyron.jl", "saturation_pressure(PCSAFT)")
    julia = os.environ.get(
        "CLAPEYRON_JULIA",
        "julia",
    )
    script = Path(__file__).with_name("collect_clapeyron.jl")
    completed = subprocess.run(
        [julia, f"--project={FORKS / 'Clapeyron.jl'}", str(script), str(T_K)],
        check=True,
        capture_output=True,
        text=True,
    )
    values = json.loads(completed.stdout.strip().splitlines()[-1])
    return row("Clapeyron", values, FORKS / "Clapeyron.jl", "saturation_pressure(PCSAFT)")


def run_mode(mode: str) -> object:
    if mode == "engine":
        return collect_engine()
    if mode == "python":
        return collect_python_packages()
    if mode == "teqp":
        return collect_teqp()
    raise ValueError(mode)


def main() -> None:
    if len(sys.argv) > 1:
        print(json.dumps(run_mode(sys.argv[1])))
        return

    comparison_python = os.environ.get("COMPARISON_PYTHON", sys.executable)
    commands = [
        [sys.executable, str(Path(__file__)), "engine"],
        [comparison_python, str(Path(__file__)), "python"],
        [os.environ.get("TEQP_PYTHON", comparison_python), str(Path(__file__)), "teqp"],
    ]
    records: list[dict[str, object]] = []
    command_log: list[list[str]] = []
    for command in commands:
        environment = os.environ.copy()
        if command[-1] == "teqp":
            teqp_python = os.environ.get("TEQP_PYTHONPATH")
            if teqp_python:
                environment["PYTHONPATH"] = teqp_python + os.pathsep + environment.get("PYTHONPATH", "")
        completed = subprocess.run(
            command, check=True, capture_output=True, text=True, env=environment
        )
        value = json.loads(completed.stdout.strip().splitlines()[-1])
        records.extend(value if isinstance(value, list) else [value])
        command_log.append(command)
    records.append(collect_clapeyron())

    RESULTS.mkdir(parents=True, exist_ok=True)
    fields = [
        "implementation",
        "solver",
        "temperature_K",
        "pressure_Pa",
        "rho_liquid_mol_m3",
        "rho_vapor_mol_m3",
        "phase_pressure_rel_error",
    ]
    with (RESULTS / "values.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(
            {key: value for key, value in record.items() if key in fields}
            for record in records
        )
    (RESULTS / "provenance.json").write_text(
        json.dumps({"case": PARAMETERS, "commands": command_log, "records": records}, indent=2)
        + "\n"
    )
    print(f"wrote {RESULTS / 'values.csv'}")


def declarations():
    return collect_engine(declare_only=True)


if __name__ == "__main__":
    main()
