"""Fig. 9 (paper PDF page 9, embedded image 700x475 px): salt solubility in [C4mim][TfO] at 298 K.
Value axis: ticks at rows 397.5 (0), 207.5 (0.0125), 16.5 (0.025). Bars are fill-colour runs on row 395
(gray = experiment, green = advanced without Bjerrum); each bar end is the centre of its outline, reached by
walking up the fill-centre column through non-white pixels; ends within 3 px of the baseline are upper bounds."""
import csv, sys, numpy as np
from PIL import Image, ImageDraw
img = Image.open(sys.argv[1]).convert('RGB'); a = np.asarray(img).astype(int)
Y0, PX = 397.5, (397.5 - 16.5) / 0.025
def kind(c):
    r, g, b = c
    if r > 228 and g > 228 and b > 228: return '.'
    if abs(r - g) < 15 and abs(g - b) < 15 and 170 <= r <= 228: return 'experiment'
    if g > 140 and r < 140 and b < 120: return 'advanced'
    return 'k'
runs, x = [], 189
while x < 694:
    k = kind(a[395, x])
    if k in ('experiment', 'advanced'):
        x0 = x
        while kind(a[395, x]) == k: x += 1
        if x - x0 >= 6: runs.append((k, x0, x - 1))
    x += 1
print(runs)
salts = ("NaF", "NaCl", "NaBr", "NaI", "KF", "KCl", "KBr")
gray = [r for r in runs if r[0] == 'experiment']; green = [r for r in runs if r[0] == 'advanced']
# NaF's green bar sits inside the baseline outline: slots come from the six detected green runs
pitch = (green[-1][1] - green[0][1]) / (len(green) - 1)
first_green = green[-1][1] - pitch * (len(salts) - 1)
rows, draw = [], ImageDraw.Draw(img)
for k, salt in enumerate(salts):
    gl = first_green + pitch * k
    for series, left in (("advanced", gl), ("experiment", gl - (green[-1][1] - gray[-1][1]))):
        xc = int(round(left + 5))
        y = int(Y0) - 1
        while a[y - 1, xc].sum() / 3 <= 228: y -= 1
        h = Y0 - (y + 0.5)
        rows.append(dict(solvent="C4mim-TfO", salt=salt, series=series, x_salt=round(h / PX, 6) if h > 3 else "",
                         x_salt_upper_bound="" if h > 3 else round(3 / PX, 6), height_px=h, pixel_x=xc))
        draw.line([xc - 5, y + 0.5, xc + 5, y + 0.5], fill=(255, 0, 255))
with open(sys.argv[2], 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n"); w.writeheader(); w.writerows(rows)
img.save(sys.argv[3])
for r in rows: print(r)
