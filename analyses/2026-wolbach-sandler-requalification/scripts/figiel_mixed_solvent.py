"""Figiel 2025 Figure 9 mixed-solvent mean ionic activity coefficients, packet v3 against v4.

For both Figiel packets, v3 resolves the water/alcohol cross-association volume without
the Gross2002 Eq. 3 size factor; v4 restores it. The reference-electrolytes packet has
no water/methanol cross edge, so only figiel-yu-held-2025 tests that pair. Each
literature row is evaluated at 298.15 K and 1 bar on its salt-free solvent composition.
"""

import csv
import hashlib
import importlib.metadata
import json
import math
from pathlib import Path

import epcsaft

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "analyses/reproduction/2025-figiel/figures/figure-09/source"
FAMILIES = ("figiel-2025-reference-electrolytes", "figiel-yu-held-2025")
PACKETS = {(f, v): ROOT / f"data/packets/{f}/{v}/parameters/parameters.json" for f in FAMILIES for v in (3, 4)}
SALTS = {"nacl": ("sodium-cation", "chloride-anion"), "nabr": ("sodium-cation", "bromide-anion")}
T, P = 298.15, 1.0e5


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def declared(family: str, version: int) -> dict:
    """The packet, with reference-electrolytes' historical water/methanol omission declared (#569): no cross
    edge, not a physical-absence claim."""
    mapping = json.loads(PACKETS[family, version].read_text())
    if family == "figiel-2025-reference-electrolytes":
        mapping["topology"]["no_cross_association"] = [{
            "component_id_a": "water", "component_id_b": "methanol",
            "provenance": {"source_id": mapping["document_id"], "locator": f"packet {version}: no water/methanol edge"}}]
    return mapping


def main() -> None:
    rows, summary = [], {}
    for family, path in ((f, p) for f in FAMILIES for p in sorted(SOURCE.glob("water-*.csv"))):
        _, solvent, salt = path.stem.split("-")
        components = ("water", solvent, *SALTS[salt])
        models = {v: epcsaft.Mixture(epcsaft.Parameters.from_mapping(declared(family, v), components=components)) for v in (3, 4)}
        references = {}
        errors = {3: [], 4: []}
        for record in csv.DictReader(path.open()):
            x_solvent = float(record[f"x_{solvent.capitalize()}"])
            molality, measured = float(record["molality"]), float(record["miac_m"])
            if x_solvent in (0.0, 1.0) or molality <= 0.0:
                continue  # single-solvent rows carry no cross-association
            row = {"packet": family, "system": path.stem, "x_solvent_salt_free": x_solvent, "molality_mol_kg": molality,
                   "miac_literature": measured}
            for version, model in models.items():
                key = (version, x_solvent)
                if key not in references:
                    references[key] = model.reference(T, P, solvent_fractions=(1.0 - x_solvent, x_solvent, 0.0, 0.0),
                                                      ions=SALTS[salt], stoichiometry=(1, 1))
                value = references[key].activity(molality).mean_ionic_molality_scale
                row[f"miac_v{version}"] = value
                errors[version].append(abs(value / measured - 1.0))
            rows.append(row)
        summary[f"{family}/{path.stem}"] = {"rows": len(errors[3]), **{f"aard_v{v}": math.fsum(e) / len(e) for v, e in errors.items()},
                              "max_relative_change_v4_vs_v3": max(abs(r["miac_v4"] / r["miac_v3"] - 1.0)
                                                                   for r in rows if r["system"] == path.stem and r["packet"] == family)}
    output = HERE / "results"
    output.mkdir(exist_ok=True)
    with (output / "figiel-figure-09.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    wheel = Path(json.loads(importlib.metadata.distribution("epcsaft").read_text("direct_url.json"))["url"].removeprefix("file://"))
    (output / "figiel-figure-09.json").write_text(json.dumps({
        "systems": summary,
        "identity": {"wheel_sha256": sha(wheel), "packets": {f"{f}/{v}": sha(p) for (f, v), p in PACKETS.items()},
                     "sources": {p.name: sha(p) for p in sorted(SOURCE.glob("water-*.csv"))}},
        "claim_limits": "298.15 K, 1 bar, Figure 9 literature rows only; relative deviation of the "
                        "molality-scale mean ionic activity coefficient; no fitted parameter changes.",
    }, indent=1) + "\n")
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
