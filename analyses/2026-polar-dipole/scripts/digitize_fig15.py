"""Digitize the Campbell et al. (1986) acetone-n-pentane symbols from GV2006 Fig. 15 (p. 1201).

The original tables (J. Chem. Eng. Data 31 (1986) 424, doi:10.1021/je00046a016) are not in the owner's
Zotero library; these are figure estimates, not printed numbers. Run with numpy/scipy/pillow available:
    pdftoppm -f 8 -l 8 -r 600 -gray -png "<Zotero DWISZT5F PDF>" /tmp/p8
    python analyses/2026-polar-dipole/scripts/digitize_fig15.py /tmp/p8-08.png
Calibration: the plot frame (x_n-pentane 0..1, P 0..2 MPa) located at 600 dpi. Filled symbols (bubble
branch, liquid x) survive a 15-px disk opening that removes the ~5-px curves; open symbols (dew branch,
vapor y) are small enclosed white interiors. Kept: x_n-pentane <= 0.60, where the branches are separated
by more than the symbol size; right of that the branches merge within 0.02 MPa and the curve bundle makes
symbol identity ambiguous. Legend and label text are outside the kept pressure bands. The pure-component
endpoints on the frame are excluded (they test only pure vapor pressures).
"""
import csv
import sys
from pathlib import Path

import numpy as np
from PIL import Image
from scipy import ndimage as nd

X0, X1, Y0, Y2 = 2964.5, 4273.5, 5465.5, 4299.5  # frame pixels: x = 0, 1; P = 0, 2 MPa
BANDS = ((372.7, 0.40, 0.70), (397.7, 0.72, 1.16), (422.6, 1.30, 1.82))  # T/K, P band / MPa
# Visually rejected detections (#28 review round 2): at 397.7 K two "filled" blobs are the curve bundle
# overlapping the open diamond at x_n-pentane 0.56-0.575; at 372.7 K one "open" interior is a white
# pocket between curves at x_n-pentane 0.559, not the open triangle.
REJECTED = {(397.7, "liquid", 0.4405), (397.7, "liquid", 0.4251), (372.7, "vapor", 0.441)}
OUT = Path(__file__).resolve().parents[1] / "inputs/gv2006-fig15-acetone-pentane-digitized.csv"


def main(page: str) -> None:
    image = np.array(Image.open(page).convert("L"))
    x_lo, y_lo = 2972, 4306
    dark = image[y_lo:5460, x_lo:4268] < 140
    yy, xx = np.mgrid[-7:8, -7:8]
    points = []
    filled, _ = nd.label(nd.binary_opening(dark, structure=(xx**2 + yy**2) <= 49))
    for i in range(1, filled.max() + 1):
        points.append(("liquid",) + nd.center_of_mass(filled == i))
    white, count = nd.label(~dark)
    for i, (size, box) in enumerate(zip(nd.sum(~dark, white, range(1, count + 1)), nd.find_objects(white)), 1):
        h, w = box[0].stop - box[0].start, box[1].stop - box[1].start
        if 40 <= size <= 400 and max(h, w) <= 30 and min(h, w) >= 5:
            points.append(("vapor",) + nd.center_of_mass(white == i))
    rows = []
    for phase, cy, cx in points:
        pentane = (cx + x_lo - X0) / (X1 - X0)
        pressure = (Y0 - (cy + y_lo)) / (Y0 - Y2) * 2.0
        band = next((T for T, low, high in BANDS if low <= pressure <= high), None)
        if band is None or pentane > 0.60 or (band, phase, round(1.0 - pentane, 4)) in REJECTED:
            continue
        rows.append({"T_K": band, "phase": phase, "x_or_y_acetone": round(1.0 - pentane, 4), "P_Pa": round(pressure * 1e6, -2),
                     "u_x_or_y_acetone": 0.004, "u_P_Pa": 10000.0})
    rows.sort(key=lambda r: (r["T_K"], r["phase"], -r["x_or_y_acetone"]))
    with OUT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f"{len(rows)} points -> {OUT}")


if __name__ == "__main__":
    main(sys.argv[1])
