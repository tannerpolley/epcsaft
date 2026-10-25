"""Digitize Cleeton 2020 Figures 2, 4 and 16 from the published PDF (Zotero 9Z3P9NXX, 300 ppi page images).

Figure 16: the eight ePC-SAFT model curves, traced column by column from a seed column where the four
curves of a panel are separate; each track takes the dark run (width <= 4 px) nearest its quadratic
prediction within GATE px and bridges only under Dicko's symbols (located from Table 1 through the
calibration). Series: at a fixed column a curve's rank from the top follows the second gas's loading
(0.8, 0.5, 0.3, 0.1), because each gas's partial pressure rises with the other's loading.
Figures 2 and 4: data symbols located by normalized cross-correlation with their legend glyphs.
Axes are calibrated on tick marks found in the images. Run from the repository root:
    python analyses/reproduction/2020-cleeton/scripts/digitize_figures.py <Cleeton PDF>
"""

import csv
import json
import math
import subprocess
import sys
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw
from scipy.ndimage import maximum_filter
from scipy.signal import fftconvolve

ROOT = Path(__file__).resolve().parents[4]
OUT = Path(__file__).resolve().parents[1] / "inputs" / "digitized"
HERE = ROOT / "build" / "cleeton-2020-figures"  # extracted page images (untracked)
IMAGE = HERE / "p12-000.jpg"
SERIES = (0.8, 0.5, 0.3, 0.1)
# pixel -> value: x linear between the 0 and the last-labelled tick, y log10 between 1e3 and 1e7.
CALIBRATION = {
    "left": {"x_px": (148.5, 830.5), "x_val": (0.0, 0.8), "y_px": (701.5, 18.5), "y_log": (3.0, 7.0),
             "columns": (151, 829), "rows": (21, 699), "legend": (553, 815, 357, 690),
             "seed": 281, "seed_rows": (312.0, 396.0, 487.0, 584.0)},
    "right": {"x_px": (1095.5, 1695.5), "x_val": (0.0, 0.8), "y_px": (694.5, 18.5), "y_log": (3.0, 7.0),
              "columns": (1099, 1769), "rows": (21, 692), "legend": (1494, 1756, 361, 680),
              "seed": 1379, "seed_rows": (213.5, 249.5, 313.5, 401.0)},
}
GATE = 5.0
MAX_GAP = 40
MARKER_RADIUS = 14.0  # markers are about 20-24 px across
DICKO = OUT.parent / "source" / "dicko-2010-table-1.csv"


def marker_pixels(cal, panel):
    """Pixel centres of Dicko's plotted symbols, from Table 1 through the calibration."""
    (xp0, xp1), (xv0, xv1) = cal["x_px"], cal["x_val"]
    (yp0, yp1), (yl0, yl1) = cal["y_px"], cal["y_log"]
    gas, series = ("co2", "figure_16_left_series_h2s_loading") if panel == "left" else \
        ("h2s", "figure_16_right_series_co2_loading")
    with open(DICKO, newline="") as stream:
        rows = [r for r in csv.DictReader(stream) if r[series]]
    return [(xp0 + (float(r[f"{gas}_loading"]) - xv0) / (xv1 - xv0) * (xp1 - xp0),
             yp0 + (math.log10(float(r[f"{gas}_pressure_mpa"]) * 1e6) - yl0) / (yl1 - yl0) * (yp1 - yp0))
            for r in rows]


def runs(column, y0, y1):
    idx = np.flatnonzero(column[y0:y1]) + y0
    out, start, prev = [], None, None
    for i in idx:
        if start is None:
            start = prev = i
        elif i == prev + 1:
            prev = i
        else:
            out.append(((start + prev) / 2.0, prev - start + 1))
            start = prev = i
    if start is not None:
        out.append(((start + prev) / 2.0, prev - start + 1))
    return out


def track(dark, cal, seed_y, step, markers):
    c0, c1 = cal["columns"]
    lx0, lx1, ly0, ly1 = cal["legend"]
    points = [(cal["seed"], seed_y)]
    col, gap = cal["seed"], 0
    while c0 <= col + step <= c1 and gap <= MAX_GAP:
        col += step
        # Quadratic through the last 30 traced columns carries the prediction under a symbol.
        recent = points[-30:]
        xb, yb = recent[-1]
        if len(recent) >= 10:
            coeff = np.polyfit([p[0] for p in recent], [p[1] for p in recent], 2)
            predicted = float(np.polyval(coeff, col))
        else:
            predicted = yb
        candidates = [y for y, w in runs(dark[:, col], *cal["rows"])
                      if w <= 4 and abs(y - predicted) <= GATE
                      and not (lx0 <= col <= lx1 and ly0 <= y <= ly1)
                      and all(math.hypot(col - mx, y - my) > MARKER_RADIUS for mx, my in markers)]
        hidden = any(math.hypot(col - mx, predicted - my) <= MARKER_RADIUS + 2 for mx, my in markers)
        if candidates:
            points.append((col, min(candidates, key=lambda y: abs(y - predicted))))
            gap = 0
        else:
            gap += 1 if hidden else MAX_GAP // 6  # bridge only under a symbol; elsewhere stop after ~6 columns
    return points if step > 0 else points[::-1]


def to_data(cal, px, py):
    (xp0, xp1), (xv0, xv1) = cal["x_px"], cal["x_val"]
    (yp0, yp1), (yl0, yl1) = cal["y_px"], cal["y_log"]
    x = xv0 + (px - xp0) / (xp1 - xp0) * (xv1 - xv0)
    y = 10 ** (yl0 + (py - yp0) / (yp1 - yp0) * (yl1 - yl0))
    return x, y


def figure_16(threshold=150, suffix=""):
    gray = np.asarray(Image.open(IMAGE).convert("L")).astype(int)
    dark = gray < threshold
    overlay = Image.open(IMAGE).convert("RGB")
    draw = ImageDraw.Draw(overlay)
    rows = []
    colors = ["#d62728", "#1f77b4", "#2ca02c", "#9467bd"]
    for panel, cal in CALIBRATION.items():
        markers = marker_pixels(cal, panel)
        for mx, my in markers:
            draw.ellipse((mx - MARKER_RADIUS, my - MARKER_RADIUS, mx + MARKER_RADIUS, my + MARKER_RADIUS),
                         outline="#ff7f0e")
        for rank, seed_y in enumerate(cal["seed_rows"]):  # rank 0 = highest curve at the seed column
            left = track(dark, cal, seed_y, -1, markers)
            right = track(dark, cal, seed_y, +1, markers)
            pts = left[:-1] + right
            draw.line(pts, fill=colors[rank], width=2)
            for px, py in pts:
                x, y = to_data(cal, px, py)
                rows.append({"panel": panel, "series_loading": SERIES[rank], "pixel_x": px, "pixel_y": py,
                             "loading": round(x, 5), "pressure_pa": round(y, 2)})
    with open(OUT / f"cleeton-2020-figure-16-curves{suffix}.csv", "w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    if not suffix:
        overlay.save(OUT / "cleeton-2020-figure-16-qa.png")
    summary = {}
    for panel in CALIBRATION:
        for rank in range(4):
            sel = [r for r in rows if r["panel"] == panel and r["series_loading"] == SERIES[rank]]
            summary[f"{panel}-{rank}"] = (sel[0]["loading"], sel[0]["pressure_pa"], sel[-1]["loading"],
                                          sel[-1]["pressure_pa"], len(sel))
    print(json.dumps(summary, indent=1))

HALF = 13
MARGIN = 8  # half-clipped symbols on the frame and frame corners are not read
# Symbols hidden under an overlapping symbol, read by eye from a 3x zoom (bounding-box centres).
MANUAL = {"fig4-co2-water": [("373", 1940, 43), ("398", 1600, 265), ("373", 1610, 266)]}
PANELS = {
    "fig2-vapor-pressure": dict(image="fig2.jpg", box=(134, 671, 21, 546), x=((355.75, 571.5), (350.0, 400.0), "linear"),
                                y=((195.5, 371.5), (2.0, 0.0), "log10"), legend=(330, 600, 380, 520),
                                glyphs={"exp": (400, 485)}),
    "fig2-density": dict(image="fig2.jpg", box=(871, 1408, 21, 546), x=((1093.5, 1308.5), (350.0, 400.0), "linear"),
                         y=((195.5, 371.5), (1000.0, 950.0), "linear"), legend=(905, 1190, 375, 515),
                         glyphs={"exp": (960, 485)}),
    "fig4-h2s-water": dict(image="fig4.jpg", box=(862, 1336, 9, 409), x=((928.5, 1200.5), (0.01, 0.03), "linear"),
                           y=((196.0, 390.5), (6.0, 5.0), "log10"), legend=(1085, 1336, 225, 392),
                           glyphs={"311": (1137, 251), "344": (1138, 286), "377": (1138, 320)}),
    "fig4-co2-water": dict(image="fig4.jpg", box=(1577, 2063, 2, 409), x=((1575.5, 1966.5), (0.0, 0.02), "linear"),
                           y=((90.0, 276.5), (7.0, 6.0), "log10"), legend=(1810, 2062, 205, 400),
                           glyphs={"298": (1865, 234), "323": (1865, 267), "373": (1865, 301), "398": (1865, 336)}),
}


def refine(dark, cx, cy):
    """Centre of the dark pixels' bounding box in a window around a legend glyph."""
    win = dark[cy - HALF:cy + HALF + 1, cx - HALF:cx + HALF + 1]
    ys, xs = np.nonzero(win)
    return int(round(cx - HALF + (xs.min() + xs.max()) / 2)), int(round(cy - HALF + (ys.min() + ys.max()) / 2))


def ncc(image, template):
    t = template - template.mean()
    num = fftconvolve(image, t[::-1, ::-1], mode="same")
    ones = np.ones_like(template)
    s1 = fftconvolve(image, ones, mode="same")
    s2 = fftconvolve(image ** 2, ones, mode="same")
    var = np.maximum(s2 - s1 ** 2 / ones.size, 1e-9)
    return num / np.sqrt(var * (t ** 2).sum())


def value(axis, pixel):
    (p0, p1), (v0, v1), scale = axis
    t = (pixel - p0) / (p1 - p0)
    return 10 ** (v0 + t * (v1 - v0)) if scale == "log10" else v0 + t * (v1 - v0)


def figures_2_4(threshold=0.55):
    rows = []
    for name, panel in PANELS.items():
        gray = np.asarray(Image.open(HERE / panel["image"]).convert("L")).astype(float)
        dark = gray < 150
        ink = 255.0 - gray
        l, r, t, b = panel["box"]
        scores = {}
        for series, (gx, gy) in panel["glyphs"].items():
            cx, cy = refine(dark, gx, gy)
            score = ncc(ink, ink[cy - HALF:cy + HALF + 1, cx - HALF:cx + HALF + 1])
            scores[series] = score
        best = np.max(np.stack(list(scores.values())), axis=0)
        peaks = (best == maximum_filter(best, size=15)) & (best > threshold)
        lx0, lx1, ly0, ly1 = panel["legend"]
        overlay = Image.open(HERE / panel["image"]).convert("RGB")
        draw = ImageDraw.Draw(overlay)
        for py, px in zip(*np.nonzero(peaks)):
            if not (l + MARGIN < px < r - MARGIN and t + MARGIN < py < b - MARGIN) or (lx0 <= px <= lx1 and ly0 <= py <= ly1):
                continue
            series = max(scores, key=lambda s: scores[s][py, px])
            rows.append({"panel": name, "series": series, "pixel_x": int(px), "pixel_y": int(py),
                         "score": round(float(best[py, px]), 3),
                         "x": value(panel["x"], px), "y": value(panel["y"], py)})
            draw.ellipse((px - 9, py - 9, px + 9, py + 9), outline="#d62728", width=2)
            draw.text((px + 10, py - 8), series, fill="#1f77b4")
        for series, px, py in MANUAL.get(name, []):
            rows.append({"panel": name, "series": series, "pixel_x": px, "pixel_y": py, "score": "manual",
                         "x": value(panel["x"], px), "y": value(panel["y"], py)})
            draw.ellipse((px - 9, py - 9, px + 9, py + 9), outline="#2ca02c", width=2)
        overlay.crop((l - 40, max(t - 10, 0), r + 5, b + 40)).save(OUT / f"cleeton-2020-{name}-qa.png")
    with open(OUT / "cleeton-2020-figures-2-4-symbols.csv", "w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    for name in PANELS:
        print(name, sum(r["panel"] == name for r in rows))


def main(pdf: str) -> None:
    HERE.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)
    for page, prefix in ((4, "p4"), (6, "p6"), (12, "p12")):
        subprocess.run(["pdfimages", "-f", str(page), "-l", str(page), "-j", pdf, str(HERE / prefix)], check=True)
    (HERE / "p4-001.jpg").replace(HERE / "fig2.jpg")
    (HERE / "p6-000.jpg").replace(HERE / "fig4.jpg")
    figure_16()
    figure_16(110, "-threshold-110")  # second lane: a darker ink threshold
    figures_2_4()


if __name__ == "__main__":
    main(sys.argv[1])
