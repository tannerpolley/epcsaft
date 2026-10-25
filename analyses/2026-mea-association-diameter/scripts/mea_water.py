"""MEA-water bubble points against Cai (1996): packet v4 against v5 and the v5 k_ij refit (#165).

v5 converts the MEA 2B association volume to the Engine sigma^3 convention. The bubble solve,
objective and scales are the #155 replay (analyses/2026-wolbach-sandler-requalification), reused.
"""

import csv
import importlib.metadata
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
ROOT = HERE.parents[1]
spec = importlib.util.spec_from_file_location("replay", ROOT / "analyses/2026-wolbach-sandler-requalification/scripts/mea_water.py")
replay = importlib.util.module_from_spec(spec)
spec.loader.exec_module(replay)
PACKETS = {v: ROOT / f"data/packets/mea-co2-h2o-nine-species-estimation/{v}/parameters/parameters.json" for v in (4, 5)}


def main() -> None:
    rows = [r for r in csv.DictReader(replay.CAI.open(encoding="utf-8-sig")) if r["role"] in replay.ADMITTED]
    documents = {v: json.loads(p.read_text()) for v, p in PACKETS.items()}
    kij = next(c["value"]["magnitude"] for p in documents[5]["pairs"] if {p["component_id_a"], p["component_id_b"]} == set(replay.COMPONENTS)
               for c in p["coefficients"] if c["family"] == "k_ij")
    table = []
    for row in rows:
        entry = {"row_id": row["row_id"], "pressure_pa": float(row["pressure_pa"]), "liquid_water": float(row["water_liquid_mole_fraction"]),
                 "temperature_cai_k": float(row["temperature_k"]), "vapor_water_cai": float(row["water_vapor_mole_fraction"])}
        for version, document in documents.items():
            entry[f"temperature_v{version}_k"], entry[f"vapor_water_v{version}"] = replay.bubble(replay.with_kij(document, kij), row)
        table.append(entry)
    budget = replay.regression.FitControls(maximum_iterations=100, maximum_elapsed_time_seconds=600.0)
    best = replay.fit_rows(replay.with_kij(documents[5], kij).parameters, rows, kij, (kij - 0.03, kij + 0.03), budget)
    if not best.usable:
        raise RuntimeError(best.message)
    with (HERE / "results/mea-water-cai.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(table[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(table)
    summary = {
        "packet_k_ij": kij,
        "at_packet_k_ij": {f"v{v}": {"mean_scaled_square": replay.objective(replay.with_kij(d, kij), rows),
                                     "max_abs_temperature_change_vs_v4_k": max(abs(e[f"temperature_v{v}_k"] - e["temperature_v4_k"]) for e in table),
                                     "max_abs_vapor_water_change_vs_v4": max(abs(e[f"vapor_water_v{v}"] - e["vapor_water_v4"]) for e in table)}
                           for v, d in documents.items()},
        "v5_refit": {**replay.fit_numbers(best), "retained_95_percent_interval": [-0.07792, -0.06914]},
        "identity": {"wheel_sha256": replay.sha(Path(json.loads(importlib.metadata.distribution("epcsaft").read_text("direct_url.json"))["url"].removeprefix("file://"))),
                     "packets": {f"v{v}": replay.sha(p) for v, p in PACKETS.items()}, "cai": replay.sha(replay.CAI)},
        "claim_limits": "25 admitted Cai 1996 rows, 66.66 and 101.33 kPa, 2B MEA; one-coordinate least squares with the retained "
                        "scales; no uncertainty interval, topology reselection or reactive evaluation.",
    }
    replay.compare_cai(HERE / "results", rows, documents, {5: best}, summary)
    (HERE / "results/mea-water-cai.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps(summary, indent=1))


def declarations():
    from epcsaft import Parameters

    return [
        (
            Parameters.from_json(path, components=replay.COMPONENTS),
            replay.observations(
                Parameters.from_json(path, components=replay.COMPONENTS),
                list(csv.DictReader(replay.CAI.open(encoding='utf-8-sig'))),
            ),
        )
        for path in PACKETS.values()
    ]


if __name__ == "__main__":
    main()
