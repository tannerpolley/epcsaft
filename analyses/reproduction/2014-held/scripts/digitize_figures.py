"""Trace Held 2014 model lines and filled markers from the PDF's embedded raster images.

Usage: digitize_figures.py IMAGE_DIR QA_DIR. IMAGE_DIR holds `pdfimages -png` output of the
Zotero PDF (TD5T732H); images are not retained. Calibration, seeds and image hashes are in
inputs/digitization.json. Each line is tracked column by column from a seed (x, y): the dark run nearest the row predicted
from the recent slope is accepted within max_step; runs outside [min_run, max_run] (markers,
crossings, lines of the other weight) are skipped, and the tracker coasts across up to max_gap columns.
"""
import csv
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw

HERE = Path(__file__).resolve().parents[1]


def runs(column):
    out, start = [], None
    for i, dark in enumerate(list(column) + [False]):
        if dark and start is None:
            start = i
        elif not dark and start is not None:
            out.append((start, i - 1))
            start = None
    return out


def main(image_dir, qa_dir):
    config = json.loads((HERE / "inputs/digitization.json").read_text())
    lines, markers = [], []
    for panel in config["panels"]:
        path = Path(image_dir) / panel["image"]
        assert hashlib.sha256(path.read_bytes()).hexdigest() == panel["image_sha256"]
        gray = np.asarray(Image.open(path).convert("L")).astype(int)
        (c0, x0), (c1, x1) = panel["x_calibration"]
        (r0, y0), (r1, y1) = panel["y_calibration"]
        col = lambda x: c0 + (x - x0) * (c1 - c0) / (x1 - x0)
        val_x = lambda c: x0 + (c - c0) * (x1 - x0) / (c1 - c0)
        row = lambda y: r0 + (y - y0) * (r1 - r0) / (y1 - y0)
        val_y = lambda r: y0 + (r - r0) * (y1 - y0) / (r1 - r0)
        top, bottom = panel["rows"]
        overlay = Image.open(path).convert("RGB")
        draw = ImageDraw.Draw(overlay)
        for series in panel.get("lines", []):
            seed_c, seed_r = round(col(series["seed"][0])), row(series["seed"][1])
            traced = {}
            for direction in (1, -1):
                # Predict from the recent slope; coast across markers and dotted gaps.
                history, r_pred, c, gap = [(seed_c, seed_r)], seed_r, seed_c, 0
                lo, hi = round(col(series["x_range"][0])), round(col(series["x_range"][1]))
                while lo <= c <= hi and gap <= series.get("max_gap", 40):
                    candidates = [(top + (s + e) / 2) for s, e in runs(gray[top:bottom, c] < series["threshold"])
                                  if series.get("min_run", 1) <= e - s + 1 <= series["max_run"]]
                    near = [r for r in candidates if abs(r - r_pred) <= series["max_step"]]
                    if near:
                        r_best = min(near, key=lambda r: abs(r - r_pred))
                        traced[c] = r_best
                        history.append((c, r_best))
                        gap = 0
                    else:
                        gap += 1
                    (ca, ra), (cb, rb) = history[max(0, len(history) - 12)], history[-1]
                    slope = (rb - ra) / (cb - ca) if cb != ca else 0.0
                    c += direction
                    r_pred = rb + slope * (c - cb)
            for c in sorted(traced)[::panel["stride"]]:
                lines.append(dict(figure=panel["figure"], series=series["name"], x=round(val_x(c), 5),
                                  y=round(val_y(traced[c]), 5), pixel_x=c, pixel_y=traced[c]))
                draw.ellipse([c - 2, traced[c] - 2, c + 2, traced[c] + 2], outline=tuple(series["qa_color"]))
        if "markers" in panel:
            spec = panel["markers"]
            dark = gray < spec["threshold"]
            # Erosion by a square keeps blobs thicker than the bold lines (filled markers).
            k = spec["erode"]
            core = np.ones_like(dark)
            for dr in range(-k, k + 1):
                for dc in range(-k, k + 1):
                    core &= np.roll(np.roll(dark, dr, 0), dc, 1)
            seen = np.zeros_like(core)
            for r, c in zip(*np.nonzero(core)):
                left, right = spec.get("cols", (0, gray.shape[1]))
                if seen[r, c] or not top <= r < bottom or not left <= c < right:
                    continue
                stack, pixels = [(r, c)], []
                seen[r, c] = True
                while stack:
                    pr, pc = stack.pop()
                    pixels.append((pr, pc))
                    for nr, nc in ((pr + 1, pc), (pr - 1, pc), (pr, pc + 1), (pr, pc - 1)):
                        if core[nr, nc] and not seen[nr, nc]:
                            seen[nr, nc] = True
                            stack.append((nr, nc))
                mr, mc = np.mean(pixels, axis=0)
                markers.append(dict(figure=panel["figure"], x=round(val_x(mc), 5), y=round(val_y(mr), 5),
                                    pixel_x=round(mc, 1), pixel_y=round(mr, 1), core_pixels=len(pixels),
                                    core_mean_gray=round(float(np.mean([gray[q] for q in pixels])), 1)))
                draw.rectangle([mc - 9, mr - 9, mc + 9, mr + 9], outline=(255, 0, 0))
        Path(qa_dir).mkdir(parents=True, exist_ok=True)
        overlay.save(Path(qa_dir) / f"{panel['figure']}-overlay.png")
    for name, rows in (("paper-lines.csv", lines), ("paper-markers.csv", markers)):
        if rows:
            with (HERE / "inputs" / name).open("w") as stream:
                writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
                writer.writeheader()
                writer.writerows(rows)
    print(len(lines), "line points;", len(markers), "markers")


if __name__ == "__main__":
    main(*sys.argv[1:])
