"""PHREEQC pitzer.dat (Reaktoro) water activity on every #556 nominal row: a reference model, not data, no limit.

Run with Lithium_Extraction's read-only environment:
~/.local/share/mamba/envs/lithium-extraction-reaktoro/bin/python analyses/2026-chloride-brine-refit/scripts/pitzer.py
"""

import csv
import json
import math
from collections import defaultdict
from pathlib import Path

import reaktoro as rkt

HERE = Path(__file__).resolve().parents[1]
ROWS = HERE.parents[0] / "2026-mixed-chloride-water-activity/results/nominal-residuals.csv"

aqueous = rkt.AqueousPhase("H2O H+ OH- Na+ Mg+2 Ca+2 Cl-")
aqueous.set(rkt.ActivityModelPitzer())
system = rkt.ChemicalSystem(rkt.PhreeqcDatabase("pitzer.dat"), aqueous)
water_mass = system.species().get("H2O").molarMass()
rows, groups = [], defaultdict(list)
for r in csv.DictReader(ROWS.open()):
    if r["declaration"] != "Held" or (r["role"] == "mixture" and r["included"] != "True"):
        continue
    m_na, m_d = float(r["m_NaCl"]), float(r["m_divalent"])
    state = rkt.ChemicalState(system)
    state.temperature(298.15, "K")
    state.pressure(1, "bar")
    for species, amount in (("H2O", 1 / water_mass), ("H+", 1e-7), ("OH-", 1e-7), ("Na+", m_na or 1e-16),
                            ("Mg+2" if r["salt"] == "MgCl2" else "Ca+2", m_d or 1e-16), ("Cl-", m_na + 2 * m_d)):
        state.set(species, amount, "mol")
    ln_aw = float(rkt.ChemicalProps(state).speciesActivityLn("H2O"))
    model = math.exp(ln_aw) if r["quantity"] == "aw" else -ln_aw / (water_mass * (2 * m_na + 3 * m_d))
    bucket = r["I_bucket"] if r["role"] == "mixture" else ""
    rows.append({"source": r["source"], "row_id": r["row_id"], "role": r["role"], "salt": r["salt"], "quantity": r["quantity"],
                 "I_bucket": bucket, "m_NaCl": m_na, "m_divalent": m_d, "observed": float(r["observed"]), "pitzer": model,
                 "residual": model - float(r["observed"])})
    groups[r["source"], r["role"], r["quantity"], bucket].append(model - float(r["observed"]))
with (HERE / "results/pitzer-rows.csv").open("w", newline="") as stream:
    writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
buckets = [{"source": s, "role": role, "quantity": q, "I_bucket": b, "count": len(e), "bias": sum(e) / len(e),
            "RMS": math.sqrt(sum(x * x for x in e) / len(e)), "maximum_absolute_error": max(map(abs, e))}
           for (s, role, q, b), e in sorted(groups.items())]
(HERE / "results/pitzer.json").write_text(json.dumps({"reaktoro": rkt.__version__, "database": "PHREEQC pitzer.dat",
                                                      "water_molar_mass_kg_mol": water_mass, "buckets": buckets}, indent=2) + "\n")
for b in buckets:
    print(b)
