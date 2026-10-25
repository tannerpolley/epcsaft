"""Map manually selected Nasrifar figure centers; inspect QA before native runs."""
import csv
import hashlib
import json
from pathlib import Path
import sys

from PIL import Image, ImageDraw

LEAF = Path(__file__).resolve().parents[1]
QA = Path(sys.argv[1]).resolve()
SOURCE = LEAF / 'inputs/source'
calibration = json.loads((SOURCE / 'figure-calibration.json').read_text())
half = calibration['coordinate_half_width_px']
rows = []
for panel, spec in calibration['panels'].items():
    path = QA / spec['image']
    assert hashlib.sha256(path.read_bytes()).hexdigest() == spec['crop_sha256']
    image = Image.open(path).convert('RGB')
    draw = ImageDraw.Draw(image)
    x0, x1, a, b = spec['x']
    y0, y1, c, d = spec['y']
    ux, uy = half * abs((b - a) / (x1 - x0)), half * abs((d - c) / (y1 - y0))
    for series, picks in spec['series'].items():
        assert max(abs(i - j) for i, j in zip(picks['picks'][0], picks['second_pick'])) <= half
        for ordinal, (px, py) in enumerate(picks['picks'], 1):
            assert x0 + half < px < x1 - half and min(y0, y1) + half < py < max(y0, y1) - half
            x, y = a + (px - x0) * (b - a) / (x1 - x0), c + (py - y0) * (d - c) / (y1 - y0)
            row_id = f'{panel}/{series}/{ordinal}'
            rows.append(dict(row_id=row_id, panel=panel, series=series, evidence=picks['evidence'],
                             x=x, y=y, x_half_width=ux, y_half_width=uy, y_scale=spec['y_scale'],
                             pixel_x=px, pixel_y=py, source_pdf_page=spec['source_page'],
                             pressure_Pa=spec.get('pressure_Pa', ''),
                             condition_note='measurement pressure 66.7 kPa; figure label 67.7 kPa; unresolved' if panel == '6c' and series == 'measured' else ''))
            draw.rectangle((px-half, py-half, px+half, py+half), outline='#D55E00', width=1)
            draw.text((px+7, py-12), row_id, fill='#0072B2')
    for px, value in [(x0, a), (x1, b)]:
        draw.line((px, min(y0,y1), px, max(y0,y1)), fill='#009E73', width=1)
        draw.text((px+4, max(y0,y1)-20), str(value), fill='#009E73')
    for py, value in [(y0, c), (y1, d)]:
        draw.line((x0, py, x1, py), fill='#009E73', width=1)
        draw.text((x0+12, py+6), str(value), fill='#009E73')
    image.save(QA / f'figure-{panel}-qa.png')
with (SOURCE / 'figure-targets.csv').open('w') as handle:
    writer = csv.DictWriter(handle, fieldnames=rows[0], lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
print(f'{len(rows)} calibrated targets; QA overlays in {QA}')
