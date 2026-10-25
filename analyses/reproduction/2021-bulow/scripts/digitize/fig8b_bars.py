"""Column scan of Fig. 8b (paper PDF page 9, embedded image 1500x457 px, right panel).

Value axis: major ticks at rows 36 (40 kJ/mol) and 385.5 (-10 kJ/mol); zero line row 315.5.
Each bar end is the centre of its 2-px black outline, found by walking from the zero line
through contiguous non-white pixels one column inside the bar's left outline. Gray bars
carry error bars: the upper cap is the top of the error line above the bar, followed on the
column (within the bar width) that is darkest just above the bar end, since the 1-2 px line
sits off the bar centre; the lower cap is the first dark row inside the gray fill.
"""
import sys, csv, numpy as np
from PIL import Image, ImageDraw
img = Image.open(sys.argv[1]).convert('RGB'); a = np.asarray(img).astype(int); L = a.sum(2) / 3
Y0, PX = 315.5, (385.5 - 36) / 50.0
val = lambda y: round((Y0 - y) / PX, 2)
white = lambda y, x: L[y, x] > 228
groups = [("methanol", ["Na+", "K+", "Cl-", "Br-"], 990), ("ethanol", ["Na+", "K+", "Cl-", "Br-", "I-"], 1218)]
rows, draw = [], ImageDraw.Draw(img)
for solvent, ions, first in groups:
    for k, ion in enumerate(ions):
        left = first + 57 * k                       # category pitch from detected gray-bar left edges
        for series, dx in (("experiment", 0), ("prediction", 16.5), ("correlation", 33)):
            xc = int(round(left + dx + 2))
            up = int(Y0) - 2
            while not white(up - 1, xc): up -= 1    # contiguous bar above the zero line
            down = int(Y0) + 2
            while not white(down + 1, xc): down += 1
            if Y0 - up > 4: end = up + 0.5
            elif down - Y0 > 4: end = down - 0.5
            else: end = Y0                           # bar no taller than the zero line and outline
            row = dict(solvent=solvent, ion=ion, series=series, value_kJ_mol=val(end), pixel_x=xc, pixel_end_y=end,
                       upper_cap_kJ_mol="", lower_cap_kJ_mol="")
            if series == "experiment":
                y = int(end) - 3
                xm = min(range(left + 2, left + 11), key=lambda x: L[y, x])
                while not white(y - 1, xm): y -= 1
                row["upper_cap_kJ_mol"] = val(y + 0.5)
                inside = [yy for yy in range(int(end) + 3, int(Y0) - 2) if L[yy, xc] < 175]
                row["lower_cap_kJ_mol"] = val(np.mean(inside[:3])) if inside else ""
                draw.line([xm - 6, y, xm + 6, y], fill=(0, 200, 0))
            rows.append(row)
            draw.line([xc - 6, end, xc + 10, end], fill=(0, 140, 255), width=1)
with open(sys.argv[2], 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n"); w.writeheader(); w.writerows(rows)
img.save(sys.argv[3])
