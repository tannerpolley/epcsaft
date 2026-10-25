#!/usr/bin/env python3
"""Copy the downstream HDA campaign (#73 S6) byte for byte into results/campaign/ with provenance.json, and the
ideal-package control's summary (diagnostic, not acceptance) as ideal_control.csv. Table logic only (downstream
decides acceptance): auxiliary_domains.csv (A10), near_active_rows.csv, measure_reference.csv. Arguments: the
downstream campaign directory, the 3b-measure s4a directory and the control's summary.csv."""

import csv
import hashlib
import json
import re
import shutil
import sys
import tomllib
from pathlib import Path

from domain_matrix import ANALYSIS, PAIRS, classify

ROOT, OUT = ANALYSIS.parents[1], ANALYSIS / "results/campaign"
COMPONENTS = ("hydrogen", "methane", "benzene", "toluene")  # downstream callback order
PARAMETERS, NASA9 = "data/packets/hda-pc-saft/1", "data/packets/nasa-glenn-2002-hda-species/1"
sha = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()  # noqa: E731
REAL_ROWS = {"fs.overhead_loss <= ub": "overhead loss", "fs.product_flow >= lb": "benzene product",
             "fs.H101.control_volume.properties_out[0.0].temperature >= lb": "H101 temperature bound",
             "fs.F102.control_volume.properties_out[0.0].pressure >= lb": "F102 pressure bound"}


def write(name, rows):
    with (OUT / name).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(dict.fromkeys(k for r in rows for k in r)), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def auxiliary_domains():
    """Every audited row against the T and P ranges of each pair, pure fit and NASA9 record of its species."""
    pure = tomllib.loads((ROOT / PARAMETERS / "sources/hda-pc-saft-records.toml").read_text())["pure_records"]
    nasa9 = json.loads((ROOT / NASA9 / "parameters/nasa9-records.json").read_text())["species"]
    refs = ([(f"pair {p}", p.split("/"), t, pr) for p, _, t, pr in PAIRS]
            + [(f"pure fit {r['component_id']}", [r["component_id"]],
                tuple(map(float, re.search(r"([\d.]+)-([\d.]+) K", r["domain"]).groups())), None) for r in pure]
            + [(f"NASA9 {s['species_id']}", [s["species_id"]], (s["intervals"][0]["lower_K"], s["intervals"][-1]["upper_K"]),
                None) for s in nasa9])
    rows = []
    for r in csv.DictReader((OUT / "audit.csv").open()):
        key = {k: r[k] for k in ("factor", "case", "block", "point")}
        if not r.get("T_point_K"):  # a refused root keeps a row
            rows.append({**key, "reference": "", "T_class": "unavailable"})
            continue
        T, P = float(r["T_point_K"]), float(r["P_Pa"])
        present = {c for c, a, b in zip(COMPONENTS, json.loads(r["x"]), json.loads(r["y"])) if a > 0 or b > 0}
        rows += [{**key, "T_point_K": T, "P_Pa": P, "reference": name, "T_class": classify((T, T), t),
                  "P_class": classify((P, P), p)} for name, species, t, p in refs if set(species) <= present]
    return rows


def near_active(campaign):
    """Every near-active row of every optimization; identity is the real-row label or the bounded quantity."""
    return [{"factor": o["factor"], "start": "/".join(map(str, o["start"])), "name": r["name"],
             "identity": REAL_ROWS.get(r["name"]) or re.sub(r".*\.(\w+\[[^\]]*\]) .*", r"\1", r["name"]),
             "counted": not r["excluded"], **{k: r[k] for k in ("slack", "row_norm", "projection_norm", "projection_error",
                                                                 "multiplier", "projection")}}
            for o in campaign["optimizations"] if "kkt" in o for r in o["kkt"]["near_active"]]


def measure_reference(s4a):
    """The four real rows and the common-optimum numbers of every 3b-measure optimization."""
    rows = []
    for path in sorted(s4a.glob("width_*_starts_*.json")):
        for o in json.loads(path.read_text())["optimizations"]:
            near = {r["name"]: r for r in o["kkt"]["near_active"]}
            rows += [{"factor": o["factor"], "start": "/".join(map(str, o["start"])), "termination": o["termination"],
                      **{k: o[k] for k in ("objective_USD_per_yr", "C101_work_W", "purity", "recovery")},
                      **{f"duty_{u}_W": v for u, v in o["duties_W"].items()},
                      **{f"vapor_fraction_{u}": v for u, v in o["vapor_fraction"].items()},
                      "row": label, "multiplier": near[name]["multiplier"], "projection": near[name]["projection"]}
                     for name, label in REAL_ROWS.items() if name in near]
    return rows


def main(source, s4a, control):
    shutil.rmtree(OUT, ignore_errors=True)
    shutil.copytree(source, OUT)
    files = {str(p.relative_to(OUT)): sha(p) for p in sorted(OUT.rglob("*")) if p.is_file()}
    assert all(sha(source / name) == digest for name, digest in files.items()), "a copy differs from its source"
    campaign = json.loads((OUT / "campaign.json").read_text())
    locks = {p["packet_path"]: p for p in json.loads((ROOT / "data/packet-locks.json").read_text())["packets"]}
    packets = {key: {**locks[path], "downstream_copy_matches": sha(ROOT / path / "parameters" / name)
                     == campaign["provenance"][f"{key}_sha256"]}
               for key, path, name in (("parameters", PARAMETERS, "parameters.json"), ("nasa9", NASA9, "nasa9-records.json"))}
    write("auxiliary_domains.csv", auxiliary_domains())
    write("near_active_rows.csv", near_active(campaign))
    write("measure_reference.csv", measure_reference(s4a))
    shutil.copyfile(control, OUT / "ideal_control.csv")
    (OUT / "provenance.json").write_text(json.dumps({
        "files": files, "downstream": campaign["provenance"], "packets": packets, "ideal_control": {str(control): sha(control)},
        "measure_reference_sources": {p.name: sha(p) for p in sorted(s4a.glob("width_*_starts_*.json"))}}, indent=1) + "\n")


if __name__ == "__main__":
    main(Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]))
