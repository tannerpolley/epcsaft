"""Map manually read original-source pixels to values; no Engine import."""
import csv
import hashlib
import json
import math
from pathlib import Path

from PIL import Image, ImageDraw

LEAF = Path(__file__).resolve().parents[1]
ROOT = LEAF.parents[2]
SOURCE = LEAF / 'inputs/source'
calibration = json.loads((SOURCE / 'figure-calibration.json').read_text())
traces = json.loads((SOURCE / 'figure-pixels.json').read_text())
corrections = {(int(r['figure']), r['panel'], r['salt'], int(r['pixel_x'])): r
               for r in csv.DictReader((SOURCE/'guide-corrections.csv').open())}
def corrected_guide(panel, px, guide):
    row = corrections.get((panel['figure'], panel['panel'], panel['salt'], px))
    return float(row['corrected_guide_y']) if row and row['disposition'] == 'recentered' else guide
rows, counts, gaps, images = [], [], [], {}
for panel, trace in zip(calibration['panels'], traces, strict=True):
    figure = panel['figure']
    source = ROOT / panel['source']
    assert hashlib.sha256(source.read_bytes()).hexdigest() == panel['sha256']
    if figure not in images:
        images[figure] = Image.open(source).convert('RGB')
    draw = ImageDraw.Draw(images[figure])
    left, top, right, bottom = panel['box']
    xmax = panel['x_range'][1]
    ymin, ymax = panel['y_range']
    draw.rectangle(panel['box'], outline='green')
    # Mid-axis ticks were independently checked against the original crop.
    for px, py in trace['line']:
        assert left - 1 <= px <= right + 1 and top - 1 <= py <= bottom + 1
    gray = Image.open(source).convert('L')
    neighbours = [t for p, t in zip(calibration['panels'], traces, strict=True)
                  if p['figure'] == figure and p['box'] == panel['box']]
    pixels = []
    draw.line([tuple(v) for v in trace['line']], fill='cyan', width=1)
    # Guides identify series; dark ink alone cannot establish ownership.
    for a, b in zip(trace['line'], trace['line'][1:]):
        for px in range(max(left+3, math.ceil(a[0])), min(right-2, math.ceil(b[0])), 3):
            guide = a[1] + (b[1]-a[1])*(px-a[0])/(b[0]-a[0])
            guide = corrected_guide(panel, px, guide)
            if (figure, panel['panel'], panel['salt'], px) in corrections:
                draw.line((px-2, round(guide), px+2, round(guide)), fill='yellow')
            corridor = 9 if figure in (1, 2, 3, 5, 8) else 5
            candidates = [y for y in range(max(top+2, round(guide)-corridor), min(bottom-1, round(guide)+corridor+1))
                          if gray.getpixel((px, y)) < 150
                          and not any((px-mx)**2+(y-my)**2 < 64 for t in neighbours for mx, my in t['markers'])
                          and not any(x0 <= px <= x1 and y0 <= y <= y1 for x0, y0, x1, y1 in trace['exclude_boxes'])]
            others = [(t, corrected_guide(p, px, u[1]+(v[1]-u[1])*(px-u[0])/(v[0]-u[0])))
                      for p, t in zip(calibration['panels'], traces, strict=True) if t in neighbours and t is not trace
                      for u, v in zip(t['line'], t['line'][1:]) if u[0] <= px < v[0]]
            owned = [y for y in candidates if abs(y-guide) <= 3 and all(abs(y-g) > abs(y-guide)+2 for _, g in others)]
            cause = '' if owned else 'owned stroke masked or absent; no source ink' if not candidates else 'stroke ownership unresolved against competing guides or outside three-pixel guide'
            reviewed = corrections.get((figure, panel['panel'], panel['salt'], px))
            if reviewed and reviewed['disposition'] == 'review':
                cause = reviewed['reason']
            if candidates:
                py = min(owned or candidates, key=lambda y: abs(y-guide))
                pixels.append((px, py, cause))
            else:
                gaps.append(dict(figure=figure, panel=panel['panel'], salt=panel['salt'], pixel_x=px, guide_y=guide, cause=cause))
                draw.point((px, round(guide)), fill='orange')
    for kind, color, points in [('line', 'red', pixels), ('markers', 'blue', trace['markers'])]:
        for order, point in enumerate(points):
            px, py = point[:2]
            cause = point[2] if kind == 'line' else 'marker centre overlaps ordinate axis; molality unresolved' if px <= left + 2 else ''
            row = {k: panel[k] for k in ['figure', 'panel', 'salt', 'temperature_K', 'quantity']}
            row.update(source_kind=kind, order=order, pixel_x=px, pixel_y=py,
                       molality=xmax * (px-left)/(right-left), value=ymin+(ymax-ymin)*(bottom-py)/(bottom-top),
                       molality_half_width=3*xmax/(right-left), value_half_width=3*(ymax-ymin)/(bottom-top),
                       label='numerical verification' if kind == 'line' else 'physical comparison',
                       qa='review' if cause else 'pass', cause=cause,
                       citation=trace['citation'] if kind == 'markers' else 'Held 2008 ePC-SAFT line',
                       source=panel['source'], source_sha256=panel['sha256'],
                       printed_page=92 if figure < 3 else 93 if figure < 6 else 94)
            rows.append(row)
            draw.ellipse((px-2, py-2, px+2, py+2), outline='magenta' if cause else color)
    counts.append(dict(figure=figure, panel=panel['panel'], salt=panel['salt'], temperature_K=panel['temperature_K'],
                       line_points=len(pixels), marker_centres=len(trace['markers']),
                       review_line_points=sum(bool(p[2]) for p in pixels),
                       gap_columns=sum(g['figure']==figure and g['panel']==panel['panel'] and g['salt']==panel['salt'] for g in gaps),
                       obscured=trace['obscured'], qa='source overlay inspected; unresolved strokes excluded'))
for name, data in [('figure-targets.csv', rows), ('figure-extraction.csv', counts), ('figure-gaps.csv', gaps)]:
    with (SOURCE/name).open('w') as handle:
        writer = csv.DictWriter(handle, fieldnames=data[0].keys(), lineterminator='\n')
        writer.writeheader()
        writer.writerows(data)
for figure, im in images.items():
    im.save(LEAF/f'figures/digitization/figure-{figure:02}-qa.png')
print(len(counts), 'series;', len(rows), 'source-selected coordinates')
