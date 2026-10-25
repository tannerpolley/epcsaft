#!/usr/bin/env python3
"""Write data/docs/parameter-index.md: every packet in data/packet-locks.json, its current
(highest locked) version and intended use, and a component -> packet index. Run after a lock
change: python3 data/tools/parameter_index.py > data/docs/parameter-index.md"""
from __future__ import annotations

import csv
import json
from collections import defaultdict
from pathlib import Path

DATA = Path(__file__).resolve().parents[1]
USES = {
    "pahlavanzadeh-2013-mdea": "Pahlavanzadeh 2013 raw molecular MDEA/water association alternatives; Engine sigma-cubed differs from printed d-cubed; Figure 4 k_ij=0; no ideal chemistry.",
    "held-2012-alcohols": "Held 2012 solvent-specific ethanol and methanol salts at 298.15 K; printed screening sigma, 0.88 sigma packing; no Born or ion-ion dispersion.",
    "nasrifar-2010-mdea": "Nasrifar and Tafazzol 2010 molecular PC-SAFT: MDEA 3B and water 2B with neutral gases; no reactive or electrolyte calculation.",
    "el-guendouzi-2004-licl-nacl-observations": "El Guendouzi 2004 LiCl-NaCl observations read from Lassin 2015 Fig. 22; inverted legend reading retained, original table not checked.",
    "robinson-covington-1968-kcl-cacl2-observations": "Robinson and Covington 1968 Table 1 KCl-CaCl2 isopiestic observations; starred measurements and pure-salt endpoints retained.",
    "yao-1992-licl-mgcl2-observations": "Yao 1992 Table 1 LiCl-MgCl2 isopiestic observations; printed typos and the illegible row remain flagged.",
    "ascani-2022-case-study-2-electrolyte-lle": "Ascani 2022 case study 2: water/1-butanol with NaCl and KCl, electrolyte LLE at 298.15 K and 1 bar.",
    "crespo-2021-deea": "Pure diethylethanolamine (neutral PC-SAFT, Crespo 2021).",
    "dielectric-response-records": "Density-dependent relative-permittivity relations (Schreckenberg, Rueben) for water, Na+ and Cl-; formulation records, not a complete parameter set.",
    "feos-esper-2023": "Esper 2023 PCP-SAFT pure parameters for 1842 substances (FeOs); non-commercial use (CC BY-NC 4.0); no mixture parameters.",
    "feos-gross-sadowski-2001": "Gross and Sadowski 2001 Table 2 nonassociating pure fluids as distributed by FeOs; values equal gross-sadowski-2001.",
    "feos-gross-sadowski-2002": "Gross and Sadowski 2002 Table 1 associating pure fluids (FeOs) with explicit Wolbach-Sandler rules; supply k_ij for a mixture.",
    "feos-rehner-2023": "Rehner 2023 binary k_ij and induced-association parameters with their Esper 2023 pure records (FeOs); non-commercial use (CC BY-NC 4.0).",
    "figiel-2025-aqueous-iodides": "Figiel 2025 aqueous NaI and KI at 298.15 K and 1 bar.",
    "figiel-2025-chloride-mixture-refit": "Aqueous chloride brines at 298.15 K: /2 adds the issue-658 Li calibration and opt-in Li-Mg term under owner-amended limits; Yao is calibration, covariance is withheld at the Born bound. The unchanged /1 Na/Mg/Ca basis retains its r1408 exception. See each packet's domain and claim limits.",
    "figiel-2025-chloride-activity-refit": "Provisional (relaxed target): aqueous MgCl2 and CaCl2 refit on the Figiel 2025 set; the packet description states its measured errors.",
    "figiel-2025-reference-electrolytes": "Figiel 2025 alkali halides in water, methanol and ethanol (modified Born term).",
    "figiel-2025-validation-parameters": "Figiel 2025 validation-side transcription, including the iodide; for validation, not a separate model.",
    "figiel-yu-held-2025": "Figiel, Yu and Held 2025 full solvent and ion set, including H+, sulfate and vanadium ions.",
    "gross-2001-methane-ethane": "Methane-ethane binary (Gross and Sadowski 2001, k_ij = 0).",
    "gross-2005": "Gross 2005 quadrupolar CO2, N2 and benzene; valid only with the quadrupole-quadrupole family, not as nonpolar PC-SAFT parameters.",
    "gross-2001-propane": "Pure propane (Gross and Sadowski 2001).",
    "gross-sadowski-2001": "Gross and Sadowski 2001 Table 2: 78 nonassociating pure fluids (alkanes, alkenes, aromatics, gases, esters, ethers, halocarbons).",
    "gross-vrabec-2006": "Gross and Vrabec 2006 dipolar fluids; valid only with the dipole-dipole family, not as nonpolar PC-SAFT parameters.",
    "hamer-wu-1972-aqueous-alkali-halides": "Observations, not parameters: evaluated molality-scale mean ionic activity coefficients of six aqueous 1:1 salts at 298.15 K.",
    "held-2014-chloride-mixture-refit": "Unqualified candidate: Na/Mg/Ca refit on chloride binaries and NaCl-MgCl2, NaCl-KCl-CaCl2 mixtures at 298.15 K (Held 2014 form); the packet description states its held-out errors.",
    "held-2014-strategy2": "Held 2014 strategy-2 water and chloride ions at 298.15 K; /2 uses the Eq. (21) 0.88 sigma closest approach that reproduces the paper figures, /1 the text's sigma; transcription, not a qualified fit.",
    "hda-pc-saft": "Toluene hydrodealkylation H2/CH4/benzene/toluene: nonassociating PC-SAFT pure records and all six k_ij (two fitted in #73).",
    "held-cameretti-sadowski-2008": "Held 2008 water and 29 ions, fully dissociated aqueous electrolytes at 298.15 K and 1 bar.",
    "hydrogen-binary-vle": "Observations, not parameters: H2/benzene bubble points and H2/CH4 coexisting compositions, cryogenic to 423 K.",
    "khudaida-2026-figure-2-electrolyte-lle": "Khudaida 2026 Figure 2: water/ethanol/isobutanol/NaCl electrolyte LLE.",
    "mea-co2-h2o-nine-species-estimation": "CO2-MEA-water nine-species model under development; carries only reactions R4 and R5 and Born constants 0 (see engine/README.md, reaction sets).",
    "schick-2023-co2-thf": "Schick 2023 dry CO2-THF binary: published pure records and k_ij for full binary VLE replay of Table S10 (no water, salts or reactions).",
    "saad-padova-marcus-1980-nacl-kcl-cacl2-observations": "Observations, not parameters: isopiestic NaCl-KCl-CaCl2 osmotic coefficients at 298.15 K (Saad, Padova and Marcus 1980).",
    "mixed-chloride-water-observations": "Source observations and reference conversions for NaCl-MgCl2 and NaCl-CaCl2 water properties at 298.15 K; no model qualification or pointwise weights.",
    "nasa-glenn-2002-hda-species": "NASA9 ideal-gas records for H2, CH4, benzene and toluene at 1 bar; thermochemistry, not PC-SAFT parameters.",
    "zuber-2014-electrolyte-permittivity": "Zuber 2014 ion-specific dielectric suppression coefficients in water, methanol and ethanol.",
}


def components(packet: Path) -> list[str]:
    found: list[str] = []
    for path in sorted((packet / "parameters").glob("*")):
        if path.suffix == ".json":
            found += [item["component_id"] for item in json.loads(path.read_text()).get("components", [])]
        elif path.suffix == ".csv":
            found += [row["component_id"] for row in csv.DictReader(path.open())]
    return list(dict.fromkeys(found))


locked: dict[str, list[str]] = defaultdict(list)
for lock in json.loads((DATA / "packet-locks.json").read_text())["packets"]:
    locked[lock["packet_id"]].append(lock["packet_version"])
missing = set(locked) - set(USES)
if missing:
    raise SystemExit(f"state the intended use of {sorted(missing)} in {__file__}")

print("# Parameter index\n")
print("Generated by `data/tools/parameter_index.py` from `data/packet-locks.json`; do not edit by hand.")
print("The current version of a packet is its highest locked version; the locks list every admitted")
print("version. Read one with `Parameters.from_json(\"data/packets/<packet>/<version>/parameters/parameters.json\",")
print("components=[...])`. Sets that share components give different results; compare them against")
print("observations (for example `hamer-wu-1972-aqueous-alkali-halides`) for the intended use.\n")
print("| packet | current | locked | intended use |\n|---|---|---|---|")
holders: dict[str, list[str]] = defaultdict(list)
for packet_id in sorted(locked):
    versions = sorted(locked[packet_id], key=int)
    print(f"| `{packet_id}` | {versions[-1]} | {', '.join(versions)} | {USES[packet_id]} |")
    for component in components(DATA / "packets" / packet_id / versions[-1]):
        holders[component].append(f"`{packet_id}/{versions[-1]}`")
print("\n## Components\n\nCurrent packet versions whose parameter files hold the component.\n")
print("| component | packets |\n|---|---|")
for component in sorted(holders):
    print(f"| `{component}` | {', '.join(holders[component])} |")
