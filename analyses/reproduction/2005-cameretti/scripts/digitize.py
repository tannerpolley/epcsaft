"""Recover the model lines and experimental crosses of Cameretti et al. (2005) Figures 1-4.

Calibration: each MATLAB axes frame spans exactly its first and last tick labels, so the
frame edges give two known points per axis. Crosses are found by matching an X template;
lines are traced column by column after the crosses are masked. Outputs:
results/digitized-lines.csv, results/digitized-points.csv and figures/figure-0N/qa-overlay.png.
"""
import csv
from pathlib import Path

import numpy as np

from water_reference import saturation_pressure
from PIL import Image, ImageDraw
from scipy.signal import correlate2d

HERE = Path(__file__).resolve().parents[1]
NACL_T = [278.15 + 5 * i for i in range(10)]
# figure, salt, quantity, frame pixels (left, right, top, bottom), axis limits (x0, x1, y0, y1),
# series temperatures in K (vapor pressure: ascending p is ascending T; density: the upper line is the colder one).
P2 = [(89.5, 440.5), (553.5, 904.5), (1016.5, 1367.5)]
R2 = [(22.5, 299.0), (388.5, 665.0), (754.0, 1030.5)]
P3 = [(100.5, 448.5), (560.5, 908.5), (1020.0, 1367.5)]
R3 = [(22.5, 296.5), (384.5, 659.0), (747.0, 1021.5)]
R4 = [(24.5, 298.5), (386.5, 659.5), (746.5, 1020.5)]
PANELS = [
    (2, "NaCl", "p", P2[0] + R2[0], (0, .3, 0, 15), NACL_T),
    (2, "LiCl", "p", P2[1] + R2[0], (0, .3, 0, 40), [303.15, 313.15, 323.15, 343.15]),
    (2, "KCl", "p", P2[2] + R2[0], (0, .3, 2, 8), [298.15, 303.15, 313.15]),
    (2, "NaBr", "p", P2[0] + R2[1], (0, .3, 0, 15), [303.15, 313.15, 323.15]),
    (2, "LiBr", "p", P2[1] + R2[1], (0, .3, 0, 15), [303.15, 313.15, 323.15]),
    (2, "KBr", "p", P2[2] + R2[1], (0, .3, 0, 8), [291.15, 303.15, 313.15]),
    (2, "NaI", "p", P2[0] + R2[2], (0, .3, 0, 20), [303.15, 313.15, 323.15, 333.15]),
    (2, "LiI", "p", P2[1] + R2[2], (0, .3, 0, 40), [303.15, 313.15, 323.15, 333.15, 343.15]),
    (2, "KI", "p", P2[2] + R2[2], (0, .3, 2, 8), [298.15, 303.15, 313.15]),
    (3, "NaCl", "rho", P3[0] + R3[0], (0, .1, 1, 1.2), [313.15, 293.15]),
    (3, "LiCl", "rho", P3[1] + R3[0], (0, .4, 1, 1.4), [313.15, 293.15]),
    (3, "KCl", "rho", P3[2] + R3[0], (0, .1, 1, 1.2), [313.15, 293.15]),
    (3, "NaBr", "rho", P3[0] + R3[1], (0, .2, 1, 1.8), [313.15, 293.15]),
    (3, "LiBr", "rho", P3[1] + R3[1], (0, .4, 1, 2), [313.15, 293.15]),
    (3, "KBr", "rho", P3[2] + R3[1], (0, .1, 1, 1.4), [313.15, 293.15]),
    (3, "NaI", "rho", P3[0] + R3[2], (0, .2, 1, 2), [313.15, 293.15]),
    (3, "LiI", "rho", P3[1] + R3[2], (0, .2, 1, 2), [313.15, 293.15]),
    (3, "KI", "rho", P3[2] + R3[2], (0, .2, 1, 1.8), [313.15, 293.15]),
    (4, "Li2SO4", "rho", (827.5, 1360.5) + R4[0], (0, .06, .9, 1.3), [333.15, 293.15]),
    (4, "Na2SO4", "rho", (827.5, 1360.5) + R4[1], (0, .05, 1, 1.4), [293.15]),
    (4, "K2SO4", "rho", (827.5, 1360.5) + R4[2], (0, .015, 1, 1.15), [293.15]),
]
CROSS = 11  # template size in pixels
MASK = 7  # pixels masked around a cross before tracing lines


def template():
    t = -0.15 * np.ones((CROSS, CROSS))
    for i in range(CROSS):
        t[i, i] = t[i, CROSS - 1 - i] = 1.0
    return t / np.abs(t).sum()


def crosses(dark):
    score = correlate2d(dark.astype(float), template(), mode="same")
    peak = (score > 0.6 * score.max())
    found = []
    for y, x in sorted(zip(*np.nonzero(peak)), key=lambda p: -score[p]):
        if all((y - b) ** 2 + (x - a) ** 2 > 16 for a, b in found):
            found.append((x, y))
    return found


def centres(column):
    rows = np.nonzero(column)[0]
    return [g.mean() for g in np.split(rows, np.nonzero(np.diff(rows) > 1)[0] + 1) if len(g) and len(g) <= 6]


def seeded(dark, guides):
    """Follow each curve rightward. Until 12 line points are found the guide row (the ideal-solution line through the
    pure-water point) predicts the curve; afterwards its own local slope does, bridging columns hidden by crosses."""
    curves = []
    for guide in guides:
        pts, slope = [], 0.0
        for col in range(MASK + 3, dark.shape[1]):
            gap = col - pts[-1][0] if pts else 0
            if gap > 80:
                break
            local = len(pts) > 12
            guess = pts[-1][1] + slope * gap if local else guide(col)
            tolerance = 3 + gap / 10 if local else min(25, 6 + gap / 4)
            near = [c for c in centres(dark[:, col]) if abs(c - guess) <= tolerance]
            if not near:
                continue
            pts.append((col, min(near, key=lambda c: abs(c - guess))))
            if len(pts) > 12:
                slope = (pts[-1][1] - pts[-12][1]) / (pts[-1][0] - pts[-12][0])
        curves.append(pts)
    return curves


def ordered(dark, nseries):
    """Paired lines spanning the same x range: assign by vertical order where all are resolved, else to the nearest."""
    columns = []
    for col in range(dark.shape[1]):
        rows = np.nonzero(dark[:, col])[0]
        columns.append([g.mean() for g in np.split(rows, np.nonzero(np.diff(rows) > 1)[0] + 1) if len(g) and len(g) <= 6])
    resolved = {col: sorted(c, reverse=True) for col, c in enumerate(columns) if len(c) == nseries}
    if not resolved:
        return [[] for _ in range(nseries)]
    known = np.array(sorted(resolved))
    curves = [[(col, resolved[col][k]) for col in known] for k in range(nseries)]
    for col, centres in enumerate(columns):
        if col in resolved or not centres or col < known[0] or col > known[-1]:
            continue
        guess = [np.interp(col, known, [resolved[c][k] for c in known]) for k in range(nseries)]
        for y in centres:
            k = int(np.argmin([abs(g - y) for g in guess]))
            if abs(guess[k] - y) <= 3:
                curves[k].append((col, y))
    return [sorted(c) for c in curves]


def smooth(pts):
    """Drop pixels more than 3 px from an iterated cubic fit of the curve (letters, crossings, jumps to a neighbour)."""
    pts = np.array(pts, dtype=float)
    for _ in range(5):
        if len(pts) < 8:
            return []
        fit = np.polyval(np.polyfit(pts[:, 0], pts[:, 1], 3), pts[:, 0])
        keep = np.abs(pts[:, 1] - fit) <= 3
        if keep.all():
            break
        pts = pts[keep]
    return [tuple(p) for p in pts]


def run():
    lines, points = [], []
    images = {}
    for figure, salt, quantity, (left, right, top, bottom), (x0, x1, y0, y1), temps in PANELS:
        path = HERE / f"figures/figure-0{figure}/source/paper-source-01-cameretti-2005-figure-00{figure}.png"
        image = images.setdefault(figure, Image.open(path).convert("RGB"))
        grey = np.array(image.convert("L"))
        box = (int(left) + 3, int(top) + 3, int(right) - 2, int(bottom) - 2)
        dark = grey[box[1]:box[3], box[0]:box[2]] < 140
        to_x = lambda px: x0 + (px - left) / (right - left) * (x1 - x0)
        to_y = lambda py: y0 + (bottom - py) / (bottom - top) * (y1 - y0)
        per_pixel = (y1 - y0) / (bottom - top)
        marks = crosses(dark)
        clean = dark.copy()
        for x, y in marks:
            clean[max(0, y - MASK):y + MASK + 1, max(0, x - MASK):x + MASK + 1] = False
        if quantity == "p":  # every vapor-pressure line starts at the pure-water saturation pressure
            def guide(T):
                return lambda col: (bottom - box[1] - (saturation_pressure(T) / 1e3 * (1 - to_x(col + box[0])) / (1 + to_x(col + box[0]))
                                                        - y0) / (y1 - y0) * (bottom - top))
            curves = seeded(clean, [guide(T) for T in temps])
        else:  # the upper density line is the colder one
            curves = ordered(clean, len(temps))
        curves = [smooth(c) for c in curves]
        draw = ImageDraw.Draw(image)
        for T, pts in zip(temps, curves):
            for col, y in pts:
                px, py = col + box[0], y + box[1]
                lines.append(dict(figure=figure, salt=salt, quantity=quantity, temperature_K=T,
                                  x_salt=round(to_x(px), 6), value=round(to_y(py), 6), value_per_pixel=round(per_pixel, 6),
                                  pixel_x=px, pixel_y=round(py, 2)))
                draw.point((px, py), fill=(220, 0, 0))
        for x, y in marks:
            px, py = x + box[0], y + box[1]
            near = [(abs(y - np.interp(x, [c for c, _ in pts], [v for _, v in pts])), T) for T, pts in zip(temps, curves)
                    if pts and pts[0][0] - 2 * MASK <= x <= pts[-1][0] + 2 * MASK]
            if not near:
                continue
            distance, T = min(near)
            if distance > 8:  # tick labels and letters, not data crosses
                continue
            points.append(dict(figure=figure, salt=salt, quantity=quantity, temperature_K=T, nearest_line_px=round(distance, 2),
                               x_salt=round(to_x(px), 6), value=round(to_y(py), 6), value_per_pixel=round(per_pixel, 6),
                               pixel_x=px, pixel_y=py))
            draw.rectangle((px - 6, py - 6, px + 6, py + 6), outline=(0, 90, 255))
    for figure, image in images.items():
        image.save(HERE / f"figures/figure-0{figure}/qa-overlay.png")
    for name, rows in (("digitized-lines.csv", lines), ("digitized-points.csv", points)):
        with (HERE / "results" / name).open("w") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)


if __name__ == "__main__":
    run()
