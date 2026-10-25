"""Fig. 7 (paper PDF page 8, embedded image 700x451 px): salt mole fraction bars at 298 K.
Value axis: ticks at rows 414.5 (0), 225.5 (0.03), 35.5 (0.06). Bars are found as fill-colour runs on row 412
(gray = experiment, green = advanced, orange = advanced + Bjerrum); each bar end is the centre of its
2-px black outline, reached by walking up the fill-centre column through non-white pixels. A bar whose
walk ends within the baseline outline (<= 3 px) is recorded with value_upper_bound only."""
import csv, sys, numpy as np
from PIL import Image, ImageDraw
img = Image.open(sys.argv[1]).convert('RGB'); a = np.asarray(img).astype(int)
Y0, PX = 414.5, (414.5 - 35.5) / 0.06
def kind(c):
    r, g, b = c
    if r > 228 and g > 228 and b > 228: return '.'
    if abs(r - g) < 15 and abs(g - b) < 15 and 170 <= r <= 228: return 'experiment'
    if g > 150 and r < 120 and b < 120: return 'advanced'
    if r > 190 and 90 < g < 170 and b < 80: return 'advanced-bjerrum'
    return 'k'
runs, x = [], 192
while x < 695:
    k = kind(a[412, x])
    if k in ('experiment', 'advanced', 'advanced-bjerrum'):
        x0 = x
        while kind(a[412, x]) == k: x += 1
        if x - x0 >= 6: runs.append((k, x0, x - 1))
    x += 1
# slots: gray left edge L, green L+15, orange L+31 (from the detected runs); pitch 56 px
slots = [("methanol", n, 198 + 56 * k) for k, n in enumerate(("NaCl", "NaBr", "KCl", "KBr"))] + \
        [("ethanol", n, 422 + 56 * k) for k, n in enumerate(("NaCl", "NaBr", "KCl", "KBr", "KI"))]
rows, draw = [], ImageDraw.Draw(img)
for solvent, salt, L in slots:
    for series, dx in (("experiment", 0), ("advanced", 15), ("advanced-bjerrum", 31)):
        xc = L + dx + 5
        y = int(Y0) - 1
        while a[y - 1, xc].sum() / 3 <= 228: y -= 1
        height_px = Y0 - (y + 0.5)
        readable = height_px > 3
        rows.append(dict(solvent=solvent, salt=salt, series=series,
                         x_salt=round(height_px / PX, 5) if readable else "",
                         x_salt_upper_bound="" if readable else round(3 / PX, 5),
                         height_px=height_px, pixel_x=xc))
        draw.line([xc - 5, y + 0.5, xc + 5, y + 0.5], fill=(255, 0, 255))
with open(sys.argv[2], 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n"); w.writeheader(); w.writerows(rows)
img.save(sys.argv[3])
for r in rows: print(r)
