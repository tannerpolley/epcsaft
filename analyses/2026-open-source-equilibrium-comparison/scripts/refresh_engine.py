"""Refresh only Engine values; external rows remain owned by their original studies."""
import csv
import hashlib
from importlib.metadata import distribution
import json
import math
import os
from pathlib import Path
import resource
import runpy
import statistics
from urllib.parse import unquote, urlparse

HERE = Path(__file__).resolve().parents[1]
ANALYSES = HERE.parent
CASES = {
    "propane": "2026-neutral-pcsaft-equilibrium-comparison",
    "bubble": "2026-binary-methane-ethane-bubble-point",
    "flash": "2026-binary-methane-ethane-tp-flash",
}


def table(path):
    with path.open() as handle:
        return list(csv.DictReader(handle))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    installed = json.loads(distribution("epcsaft").read_text("direct_url.json"))
    wheel = Path(unquote(urlparse(installed["url"]).path))
    wheel_sha = sha(wheel)
    flash_dir = ANALYSES / CASES["flash"] / "results"
    flash_meta = flash_dir / "performance-current.json"
    flash_samples = flash_dir / "performance-current.csv"
    assert json.loads(flash_meta.read_text())["wheel_sha256"] == wheel_sha, \
        "Rerun performance.py for the installed wheel before reusing its flash samples"
    samples = [r for r in table(flash_samples) if r["case"] == "0" and r["mode"] == "cold"]
    assert len(samples) >= 5 and all(r["success"] == "True" for r in samples)
    first = samples[0]
    rows, timing_rows = [], []
    sources = [flash_meta, flash_samples]
    for case, folder in CASES.items():
        source = ANALYSES / folder
        if case == "flash":
            deciles = statistics.quantiles([float(r["wall_ms"]) for r in samples],
                                          n=10, method="inclusive")
            row = dict(case=case, temperature_K=first["T_K"], pressure_Pa=first["P_Pa"],
                rho_liquid_mol_m3=first["rho_L"], rho_vapor_mol_m3=first["rho_V"],
                x_methane=first["x_methane"], y_methane=first["y_methane"],
                beta_vapor=first["beta"],
                median_ms=statistics.median(float(r["wall_ms"]) for r in samples),
                p10_ms=deciles[0], p90_ms=deciles[-1], samples=len(samples), warmups=2)
        else:
            collector = source / "scripts/collect_results.py"
            record = runpy.run_path(str(collector))["collect_engine"]()
            row = {key: record[key] for key in (
                "temperature_K", "pressure_Pa", "rho_liquid_mol_m3", "rho_vapor_mol_m3")}
            row["case"] = case
            if case == "bubble":
                row.update(x_methane=record["liquid_x_methane"],
                           y_methane=record["vapor_y_methane"],
                           median_ms=record["solve_median_ms"],
                           p10_ms=record["solve_p10_ms"], p90_ms=record["solve_p90_ms"],
                           samples=record["samples"], warmups=record["warmups"])
                timing_rows = [dict(sample=i, wall_ms=v)
                               for i, v in enumerate(record["solve_samples_ms"])]
            sources.append(collector)
        oracle = source / "results/values.csv"
        for reference in table(oracle):
            for key in ("pressure_Pa", "rho_liquid_mol_m3", "rho_vapor_mol_m3",
                        "x_methane", "y_methane", "beta_vapor"):
                external_key = {"x_methane": "liquid_x_methane",
                                "y_methane": "vapor_y_methane"}.get(key, key) if case == "bubble" else key
                if key in row and reference.get(external_key):
                    value, expected = float(row[key]), float(reference[external_key])
                    tolerance = 1e-7 if key in ("pressure_Pa", "rho_liquid_mol_m3", "rho_vapor_mol_m3") else 1e-8
                    assert math.isfinite(value) and abs(value-expected) <= tolerance*max(1., abs(expected)), \
                        (case, reference["implementation"], key, value, expected)
        rows.append(row)
        sources.extend([oracle, source / "results/provenance.json"])
    output = HERE / "results"
    output.mkdir(exist_ok=True)
    for name, data, fields in (
        ("current-engine.csv", rows, ["case", "temperature_K", "pressure_Pa",
            "rho_liquid_mol_m3", "rho_vapor_mol_m3", "x_methane", "y_methane",
            "beta_vapor", "median_ms", "p10_ms", "p90_ms", "samples", "warmups"]),
        ("bubble-timings.csv", timing_rows, ["sample", "wall_ms"])):
        with (output / name).open("w") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
            writer.writeheader()
            writer.writerows(data)
    metadata = dict(wheel_sha256=wheel_sha, source_commit=record["source"]["commit"],
        source_dirty=record["source"]["dirty"],
        sources={str(p.relative_to(ANALYSES)): sha(p) for p in sources},
        parameter_sha256={name: sha(ANALYSES.parent / "data/packets" / name / "2/parameters/parameters.json")
                          for name in ("gross-2001-propane", "gross-2001-methane-ethane")},
        peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        affinity=sorted(os.sched_getaffinity(0)), nice=os.getpriority(os.PRIO_PROCESS, 0),
        limits="Engine-only refresh. Flash reuses same-wheel central cold samples, not the nine-case pooled median. External measurements are historical, not rerun. No pure-fluid runtime measurement.")
    (output / "provenance.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(json.dumps(rows, indent=2))


if __name__ == "__main__":
    main()
