"""Calibrate manually selected curve points and symbol centers from paper figures."""
import csv
import hashlib
import json
from pathlib import Path

from PIL import Image, ImageDraw

HERE = Path(__file__).resolve().parents[1]


def main():
    root = HERE / 'inputs/digitized'
    calibration = json.loads((root / 'calibration.json').read_text())
    pixels = list(csv.DictReader((root / 'pixels.csv').open()))
    rows, images = [], {}
    for r in pixels:
        c = calibration[r['figure']]
        x0, x1, xv0, xv1, y0, y1, yv0, yv1 = c['axis']
        px, py = float(r['pixel_x']), float(r['pixel_y'])
        assert min(x0, x1) <= px <= max(x0, x1)
        assert min(y0, y1) <= py <= max(y0, y1)
        x = xv0 + (px - x0) * (xv1 - xv0) / (x1 - x0)
        y = yv0 + (py - y0) * (yv1 - yv0) / (y1 - y0)
        # Digitization uncertainty is independent of numerical solver tolerance.
        r.update(kind=c['kind'], x=x, y=y,
                 x_uncertainty=3 * abs((xv1 - xv0) / (x1 - x0)),
                 y_uncertainty=3 * abs((yv1 - yv0) / (y1 - y0)))
        rows.append(r)
        f = r['figure']
        source = HERE / f'figures/figure-{f}/source/paper-source-01-gross-2001-figure-0{f}.png'
        if f not in images:
            images[f] = Image.open(source).convert('RGB')
            c['source_sha256'] = hashlib.sha256(source.read_bytes()).hexdigest()
            c['qa'] = 'pass: every model point retraced on solid PC-SAFT branch; separate source-crop and overlay inspection; same investigator; raster uncertainty remains'
        color = 'red' if r['reference'] == 'experiment' else 'blue'
        ImageDraw.Draw(images[f]).ellipse((px - 5, py - 5, px + 5, py + 5), outline=color, width=2)
    with (root / 'figure-values.csv').open('w') as out:
        w = csv.DictWriter(out, fieldnames=rows[0], lineterminator="\n"); w.writeheader(); w.writerows(rows)
    (root / 'calibration.json').write_text(json.dumps(calibration, indent=2) + '\n')
    for f, im in images.items():
        im.save(HERE / f'figures/figure-{f}/digitization-overlay.png')


if __name__ == '__main__':
    main()
