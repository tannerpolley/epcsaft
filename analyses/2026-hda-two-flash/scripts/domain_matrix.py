#!/usr/bin/env python3
"""Pair-by-state evidence domains (#73): every HDA state class against each binary pair's source
temperature and pressure evidence, classified below / within / above (or unstated)."""

from __future__ import annotations

import csv
from pathlib import Path

ANALYSIS = Path(__file__).resolve().parents[1]
# (pair, source, T range K, P range Pa or None when the source states none)
PAIRS = [
    ("hydrogen/methane", "Benham-Katz 1957, Sagara 1972 (fit); Hong-Kobayashi 1981 (held out)", (103.15, 183.12), (75200.0, 28406000.0)),
    ("hydrogen/benzene", "Park 1996 Table 1 (323.2/423.2 K fit, 373.2 K held out)", (323.2, 423.2), (2550000.0, 15730000.0)),
    ("hydrogen/toluene", "Saajanlehto 2014 Table 6 (Simnick 1978, Laugier 1980 data)", (502.15, 575.15), None),
    ("methane/benzene", "Gross-Sadowski 2001 Table 4", (421.15, 501.15), None),
    ("methane/toluene", "Gross-Sadowski 2001 Table 4", (422.15, 543.15), None),
    ("benzene/toluene", "Yu 2026 Table S5 (DECHEMA infinite dilution)", (213.15, 576.15), None),
]
STATES = [  # (state, T range K, P range Pa)
    ("feed", (303.2, 303.2), (350000.0, 350000.0)),
    ("heater/reactor", (500.0, 800.0), (350000.0, 350000.0)),
    ("F101", (298.0, 450.0), (350000.0, 350000.0)),
    ("F102", (298.0, 450.0), (105000.0, 150000.0)),
]


def classify(state: tuple[float, float], evidence: tuple[float, float] | None) -> str:
    if evidence is None:
        return "unstated"
    parts = []
    if state[0] < evidence[0]:
        parts.append("below")
    if state[1] >= evidence[0] and state[0] <= evidence[1]:
        parts.append("within")
    if state[1] > evidence[1]:
        parts.append("above")
    return "+".join(parts)


def main() -> None:
    rows = [{"pair": pair, "evidence": source, "T_min_K": t[0], "T_max_K": t[1],
             "P_min_Pa": p[0] if p else "", "P_max_Pa": p[1] if p else "", "state": state,
             "state_T_K": f"{st[0]}-{st[1]}", "state_P_Pa": f"{sp[0]}-{sp[1]}",
             "temperature": classify(st, t), "pressure": classify(sp, p)}
            for pair, source, t, p in PAIRS for state, st, sp in STATES]
    with (ANALYSIS / "results/domain_matrix.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    for row in rows:
        print(row["pair"], row["state"], row["temperature"], row["pressure"])


if __name__ == "__main__":
    main()
