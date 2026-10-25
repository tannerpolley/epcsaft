"""Evaluated-reference saturation rows for the pure-fluid comparison (#28).

Pure-fluid reference equations of state through CoolProp 8.0.0 (not an Engine dependency):
    python -m venv /tmp/coolprop && /tmp/coolprop/bin/pip install CoolProp==8.0.0
    /tmp/coolprop/bin/python analyses/2026-polar-dipole/scripts/reference_saturation.py
The temperature grid spans each fluid's GV2006 Table 2 range where stated (acetone: triple point
upward), clipped to the reference EOS validity range and capped at 0.97 T_c. Rows alternate
calibration/validation for the documented nonpolar refit. Each row carries the relative uncertainty the
EOS paper states for p_sat and saturated-liquid density at that temperature; blank means none is stated.
"""
import csv
from pathlib import Path

import CoolProp
from CoolProp.CoolProp import PropsSI

OUT = Path(__file__).resolve().parents[1] / "inputs/reference-saturation.csv"
FLUIDS = {  # component id: (CoolProp name, reference EOS, T_min/K, T_max/K)
    "acetone": ("Acetone", "Lemmon and Span, J. Chem. Eng. Data 51 (2006) 785", 180.0, 490.0),
    "dimethyl-ether": ("DimethylEther", "Wu, Zhou and Lemmon, J. Phys. Chem. Ref. Data 40 (2011) 023104", 200.0, 385.0),
    # CoolProp's DiethylEther is Thol-IJT-2014, valid 270-500 K; GV2006's 220-466 K range starts below it.
    "diethyl-ether": ("DiethylEther", "Thol, Piazza and Span, Int. J. Thermophys. 35 (2014) 783", 270.0, 450.0),
}


def uncertainty(component: str, T: float, Tc: float) -> tuple[float | None, float | None]:
    """Stated relative uncertainty of (p_sat, saturated-liquid density); None where the paper states none."""
    if component == "acetone":  # Lemmon and Span 2006, acetone subsection, closing uncertainty statement
        p = 0.0025 if 290.0 <= T <= 390.0 else 0.005 if T > 270.0 else None
        rho = 0.001 if 280.0 <= T <= 310.0 else 0.005 if T < 380.0 else 0.01
        return p, rho
    if component == "dimethyl-ether":  # Wu et al. 2011, abstract and Sec. 5 (k = 2)
        # "0.05 %, except for very near the critical point": the 0.95 T_c cutoff is this analysis's choice.
        return (0.002 if T > 230.0 else None), (0.0005 if T <= 0.95 * Tc else None)
    # Thol et al. 2014, Sect. 3.1 pp. 800-801 (Zotero XRMZRYJI): vapor pressures deviate "less than +-0.4 %";
    # saturated liquid densities "within +-0.2 %", increasing "close to the critical point" (same 0.95 T_c cutoff).
    # These are deviation envelopes, not a k = 2 uncertainty. The paper's p. 790 rule instead takes the AAD to the best
    # data set as the uncertainty indication (Table 7: 0.168 % p_sat, 0.041 % rho_liq); the README reports that alternative.
    return 0.004, (0.002 if T <= 0.95 * Tc else None)


rows = []
for component, (name, source, low, high) in FLUIDS.items():
    Tc = PropsSI("Tcrit", name)
    high = min(high, 0.97 * Tc)
    count = 24
    for k in range(count):
        T = round(low + (high - low) * k / (count - 1), 2)
        u_p, u_rho = uncertainty(component, T, Tc)
        rows.append({"component_id": component, "T_K": T, "p_sat_Pa": repr(PropsSI("P", "T", T, "Q", 0, name)),
                     "rho_liq_mol_m3": repr(PropsSI("Dmolar", "T", T, "Q", 0, name)),
                     "role": "calibration" if k % 2 == 0 else "validation", "reference_eos": source,
                     "u_p_sat_rel": "" if u_p is None else u_p, "u_rho_liq_rel": "" if u_rho is None else u_rho})
with OUT.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
print(f"CoolProp {CoolProp.__version__}: {len(rows)} rows -> {OUT}")
