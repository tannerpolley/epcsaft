"""Write the explicit water/Na+/Cl- selections compared in #27.

Each selection is the Figiel 2025 base (figiel-yu-held-2025/3) restricted to
water, Na+ and Cl-, with only its permittivity formulation exchanged. The
density-temperature and dipolar-perturbation records come from
data/packets/dielectric-response-records/1. Records the exchanged form does not
consume (Figiel's water permittivity 78.09 and C_epsilon = 7.01) are removed,
so the admitted builder sees one complete formulation per document.
"""

from __future__ import annotations

import copy
import csv
import json
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / "data/packets/figiel-yu-held-2025/3/parameters/parameters.json"
RECORDS = ROOT / "data/packets/dielectric-response-records/1/parameters/dielectric-records.csv"
OUT = Path(__file__).resolve().parents[1] / "inputs/selections"
COMPONENTS = ("water", "sodium-cation", "chloride-anion")
FORMS = {
    "density-temperature": ("solvent-density-temperature", {"dielectric_volume", "dielectric_temperature"}),
    "dipolar-perturbation": ("dipolar-perturbation", {"dipole_scaling", "polarizability_scaling", "correlation_integral_parameter"}),
}


def restricted_base() -> dict:
    document = json.loads(BASE.read_text())
    keep = set(COMPONENTS)
    document["components"] = [c for c in document["components"] if c["component_id"] in keep]
    document["pairs"] = [p for p in document["pairs"] if {p["component_id_a"], p["component_id_b"]} <= keep]
    document["correlations"] = [c for c in document["correlations"] if c["component_id"] in keep]
    topology = document["topology"]
    topology["sites"] = [s for s in topology["sites"] if s["component_id"] in keep]
    topology["edges"] = [
        e for e in topology["edges"]
        if e["endpoint_a"]["component_id"] in keep and e["endpoint_b"]["component_id"] in keep
    ]
    return document


def exchanged(base: dict, choice: str, families: set[str]) -> dict:
    document = copy.deepcopy(base)
    record = next(f for f in document["model_families"] if f["kind"] == "permittivity")
    record["choice"] = choice
    record["provenance"] = {
        "source_id": "dielectric-response-records",
        "locator": f"#27 hybrid: Figiel 2025 base with the {choice} relation (full EOS 6.2/6.4)",
        "domain_id": "aqueous-298-15-k",
    }
    document["model_coefficients"] = [
        m for m in document["model_coefficients"] if m["family"] != "ion_fraction_suppression_coefficient"
    ]
    document["correlations"] = [c for c in document["correlations"] if c["family"] != "relative_permittivity"]
    with RECORDS.open(newline="") as stream:
        rows = [r for r in csv.DictReader(stream) if r["family"] in families]
    by_id = {c["component_id"]: c for c in document["components"]}
    for row in rows:
        by_id[row["component_id"]]["coefficients"].append({
            "identity": f"component/{row['component_id']}/{row['family']}",
            "family": row["family"],
            "value": {"magnitude": float(row["magnitude"]), "unit": row["unit"]},
            "provenance": {"source_id": row["source_id"], "locator": row["source_locator"], "domain_id": "aqueous-298-15-k"},
        })
    packet = tomllib.loads((RECORDS.parents[1] / "packet.toml").read_text())
    used = {row["source_id"] for row in rows}
    document["sources"] += [dict(s) for s in packet["sources"] if s["source_id"] in used]
    document["domains"] += [dict(d) for d in packet["domains"]]
    document["document_id"] = f"issue-27-figiel-{choice}"
    return document


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    base = restricted_base()
    selections = {"ion-fraction-suppression": base}
    selections.update({name: exchanged(base, choice, families) for name, (choice, families) in FORMS.items()})
    for name, document in selections.items():
        (OUT / f"{name}.json").write_text(json.dumps(document, indent=1) + "\n")


if __name__ == "__main__":
    main()
