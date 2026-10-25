"""Digitize the corrigendum's corrected Figs. 6-7 (Fluid Phase Equilib. 548 (2021) 113184).

Input: the two embedded 300-ppi JPEGs extracted with `pdfimages -j` from Zotero attachment YHZIPCDI, filed under the
paper parent L79J5R52 and byte-identical to CZQJCDPD under the corrigendum parent 9XTQJ9P7 (PDF SHA-256 in
`figures/digitized-calibration.json`); page images are not retained.
Usage: python scripts/digitize.py <dir holding corrimg-003.jpg and corrimg-004.jpg>
Axes are calibrated on the plot frame, whose edges are the axis end ticks; frame pixels were located
as the dark rows/columns spanning each panel. Curves: per-column median of the series colour.
Data points: dark connected blobs of disc size. Writes CSV and a QA overlay per panel.
Data-disc rule, the same for every curve: curves are 2-4 px wide and the filled data discs 12-16 px, and the
coloured curves are drawn over the discs, which can split a disc's dark pixels into undetected pieces. Ink of
any colour at least 7 px thick is therefore disc (or curve-on-disc overlap); a curve column with any series
pixel within 2 px of it is excluded and kept in the CSV as series "<name> inside data disc".
"""
import csv
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw
from scipy import ndimage

HERE = Path(__file__).resolve().parents[1]
# panel: image, frame (left px, right px, bottom px, top px), axis (x_left, x_right, y_bottom, y_top), series
PANELS = {
    "figure-06/libr-methanol-left": ("corrimg-003.jpg", (88.5, 621.0, 418.5, 17.5), (0, 0.2, 0, 4),
                                     ("advanced", "born-constant-permittivity", "points")),
    "figure-06/libr-methanol-right": ("corrimg-003.jpg", (894.5, 1427.5, 418.5, 17.5), (0, 0.2, -4, 4),
                                      ("advanced", "born-contribution", "debye-huckel-contribution", "points")),
    "figure-07/libr-methanol": ("corrimg-004.jpg", (109.5, 606.5, 389.0, 15.5), (0, 0.2, 0, 4),
                                ("advanced", "revised", "points")),
    "figure-07/nabr-ethanol": ("corrimg-004.jpg", (912.5, 1409.5, 389.0, 15.5), (0, 0.010, 0, 1),
                               ("advanced", "revised", "points")),
    "figure-07/licl-methanol": ("corrimg-004.jpg", (109.5, 606.5, 964.5, 590.5), (0, 0.15, 0, 2),
                                ("advanced", "revised", "points")),
    "figure-07/nabr-methanol": ("corrimg-004.jpg", (912.5, 1409.5, 961.5, 590.0), (0, 0.05, 0, 1),
                                ("advanced", "revised", "points")),
}
INSET = 6  # px kept clear of the frame and its inward ticks


def masks(rgb):
    r, g, b = (rgb[..., k].astype(int) for k in range(3))
    gray = (r + g + b) / 3
    return {
        "advanced": (g - r > 40) & (g - b > 40),                       # green
        "born-constant-permittivity": (r > 190) & (g > 90) & (g < 190) & (b < 130),  # orange
        "born-contribution": (r > 190) & (g > 90) & (g < 190) & (b < 130),
        "debye-huckel-contribution": (b > 170) & (r < 90) & (g < 90),  # blue
        "dark": (gray < 90) & (np.max(rgb, 2).astype(int) - np.min(rgb, 2) < 50),
    }


def main(image_dir):
    calibration = {"source_pdf_sha256": hashlib.sha256(next(Path("/home/tnnrpolley21/Zotero/storage/YHZIPCDI").glob("*.pdf")).read_bytes()).hexdigest(),
                   "zotero": "corrigendum PDF YHZIPCDI, filed under the paper parent L79J5R52; the corrigendum parent 9XTQJ9P7 holds the byte-identical PDF CZQJCDPD", "inset_px": INSET,
                   "data_disc_rule": "exclude curve columns with a pixel within 2 px of ink (mean RGB < 170) that survives opening with a 7-px disc", "panels": {}}
    for panel, (image, (left, right, bottom, top), (x0, x1, y0, y1), series) in PANELS.items():
        path = Path(image_dir) / image
        rgb = np.asarray(Image.open(path).convert("RGB"))
        m = masks(rgb)
        roi = np.zeros(rgb.shape[:2], bool)
        roi[int(top) + INSET:int(bottom) - INSET + 1, int(left) + INSET:int(right) - INSET + 1] = True
        to_x = lambda px: x0 + (px - left) * (x1 - x0) / (right - left)
        to_y = lambda py: y0 + (py - bottom) * (y1 - y0) / (top - bottom)
        # Filled data discs are ~12-16 px across and curves 2-4 px wide: an opening with a 7-px disc keeps only discs.
        yy, xx = np.mgrid[-3:4, -3:4]
        discs = ndimage.binary_opening(m["dark"] & roi, structure=xx**2 + yy**2 <= 9)
        labels, n = ndimage.label(discs)
        blobs = ndimage.find_objects(labels)
        disc = [i for i, sl in enumerate(blobs) if (sl[0].stop - sl[0].start) >= 8 and (sl[1].stop - sl[1].start) >= 8]
        curve_labels, _ = ndimage.label(ndimage.binary_dilation(discs, iterations=2))
        ink = rgb.astype(int).mean(2) < 170  # half-intensity edge of any colour
        on_disc = ndimage.binary_dilation(ndimage.binary_opening(ink & roi, structure=xx**2 + yy**2 <= 9), iterations=2)
        rows = []
        for name in series:
            if name == "points":
                for i in disc:
                    cy, cx = ndimage.center_of_mass(labels == i + 1)
                    w = blobs[i][1].stop - blobs[i][1].start
                    rows.append(dict(series="points", x=to_x(cx), y=to_y(cy), pixel_x=cx, pixel_y=cy,
                                     note="merged discs" if w > 18 else ""))
                continue
            if name == "revised":
                mask = m["dark"] & roi & (curve_labels == 0)
            else:
                mask = m[name] & roi
            for cx in range(mask.shape[1]):
                ys = np.flatnonzero(mask[:, cx])
                if len(ys):
                    if name == "revised" and np.ptp(ys) > 8:
                        continue  # column also crossed by green/points; skip ambiguous column
                    cy = float(np.median(ys))
                    label = f"{name} inside data disc" if on_disc[ys, cx].any() else name
                    rows.append(dict(series=label, x=to_x(cx), y=to_y(cy), pixel_x=cx, pixel_y=cy, note=""))
        out = HERE / "figures" / f"{panel}.csv"
        out.parent.mkdir(parents=True, exist_ok=True)
        with out.open("w") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n")
            w.writeheader()
            w.writerows(rows)
        overlay = Image.fromarray(rgb).convert("RGB")
        draw = ImageDraw.Draw(overlay)
        draw.rectangle((left, top, right, bottom), outline=(255, 0, 255))
        for r in rows:
            c = {"points": (255, 0, 0), "revised": (0, 200, 255)}.get(r["series"], (150, 150, 150) if "inside" in r["series"] else (255, 0, 255))
            draw.point((r["pixel_x"], r["pixel_y"]), fill=c)
            if r["series"] == "points":
                draw.ellipse((r["pixel_x"] - 3, r["pixel_y"] - 3, r["pixel_x"] + 3, r["pixel_y"] + 3), outline=c)
        overlay.crop((int(left) - 40, int(top) - 10, int(right) + 10, int(bottom) + 40)).save(out.with_suffix(".qa.png"))
        calibration["panels"][panel] = dict(image=image, image_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                                            frame_px=dict(left=left, right=right, bottom=bottom, top=top),
                                            axis=dict(x_left=x0, x_right=x1, y_bottom=y0, y_top=y1),
                                            px_per_x=(right - left) / (x1 - x0), px_per_y=(bottom - top) / (y1 - y0),
                                            counts={s: sum(r["series"] == s for r in rows) for s in sorted({r["series"] for r in rows})})
        print(panel, calibration["panels"][panel]["counts"])
    (HERE / "figures/digitized-calibration.json").write_text(json.dumps(calibration, indent=1) + "\n")


if __name__ == "__main__":
    main(sys.argv[1])
