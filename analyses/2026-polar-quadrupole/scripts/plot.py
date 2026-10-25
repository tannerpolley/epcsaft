"""Saturation-pressure deviation of each model from its reference, read from results/engine-rows.csv (no model run).

    .venv/bin/python analyses/2026-polar-quadrupole/scripts/plot.py
Points are Engine predictions at the retained grid temperatures; the zero line is the reference
(Span-Wagner CO2 EOS, DIPPR-101 benzene correlation).
"""
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parents[1]
rows = list(csv.DictReader((HERE / "results/engine-rows.csv").open()))
figure, axes = plt.subplots(1, 2, figsize=(9.0, 3.6), constrained_layout=True)
for axis, species, reference in zip(axes, ("carbon-dioxide", "benzene"), ("Span–Wagner EOS", "DIPPR-101 correlation")):
    for kind, label, marker, color in (("p_sat", "PCP-SAFT with quadrupole, Gross 2005 set", "o", "#0072B2"),
                                       ("p_sat_nonpolar_gs2001", "PC-SAFT, Gross–Sadowski 2001 set", "x", "#D55E00")):
        group = [r for r in rows if r["kind"] == kind and r["case"].startswith(species + " ")]
        T = [float(r["case"].split()[1]) for r in group]
        deviation = [100.0 * (float(r["engine"]) / float(r["reference"]) - 1.0) for r in group]
        axis.plot(T, deviation, marker, color=color, label=label, markerfacecolor="none" if marker == "o" else None)
    axis.axhline(0.0, color="0.3", linewidth=0.8)
    axis.set_xlabel("Temperature (K)")
    axis.set_ylabel(f"$p^{{sat}}$ deviation from {reference} (%)")
    axis.set_title("CO$_2$" if species == "carbon-dioxide" else "Benzene")
axes[0].legend(frameon=False, fontsize=8)
figure.savefig(HERE / "figures/saturation-deviation.png", dpi=150)
