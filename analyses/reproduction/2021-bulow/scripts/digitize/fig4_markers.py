"""Black filled markers in corrigendum Fig. 4 / SI Fig. S1: discs of dark ink, possibly crossed by coloured curves.
Ink = dark (V<120) or saturated colour; opened with a 7-px disc; kept where it touches dark ink; opened again
with a 5-px disc. Blobs under 40 px are fragments; centroids closer than 10 px are one marker (merged=True)."""
import cv2, csv, json, sys, numpy as np
out = []
for name in ("fig4-licl", "fig4-libr", "figS1-licl-ethanol"):
    cfg = json.load(open(f"{name}.json")); ax = cfg["axis"]
    f = "cimg-004.png" if name.startswith("figS1") else "cimg-003.png"
    im = cv2.imread(f); hsv = cv2.cvtColor(im, cv2.COLOR_BGR2HSV)
    dark = hsv[..., 2] < 120; ink = ((dark | (hsv[..., 1] > 120)) * 255).astype(np.uint8)
    l, t, r, b = cfg["roi"]; sub = np.zeros_like(ink); sub[t:b, l:r] = ink[t:b, l:r]
    o = cv2.morphologyEx(sub, cv2.MORPH_OPEN, cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (7, 7)))
    o &= cv2.dilate((dark * 255).astype(np.uint8), np.ones((3, 3), np.uint8))
    o = cv2.morphologyEx(o, cv2.MORPH_OPEN, cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (5, 5)))
    n, lab, st, cen = cv2.connectedComponentsWithStats(o)
    blobs = [(c[0], c[1], s[4]) for c, s in zip(cen[1:], st[1:]) if s[4] >= 40]
    groups = []
    for x, y, area in sorted(blobs):
        if groups and abs(groups[-1][-1][0] - x) < 10 and abs(groups[-1][-1][1] - y) < 10: groups[-1].append((x, y, area))
        else: groups.append([(x, y, area)])
    sx = (ax["x2_value"] - ax["x1_value"]) / (ax["x2_pixel"][0] - ax["x1_pixel"][0])
    sy = (ax["y2_value"] - ax["y1_value"]) / (ax["y2_pixel"][1] - ax["y1_pixel"][1])
    for g in groups:
        w = sum(a for *_, a in g); px = sum(x * a for x, _, a in g) / w; py = sum(y * a for _, y, a in g) / w
        out.append(dict(panel=name, series="experiment", x=ax["x1_value"] + (px - ax["x1_pixel"][0]) * sx,
                        y=ax["y1_value"] + (py - ax["y1_pixel"][1]) * sy, pixel_x=round(px, 1), pixel_y=round(py, 1),
                        area_px=int(w), merged=len(g) > 1 or w > 170))
        cv2.circle(im, (int(px), int(py)), 7, (255, 0, 255), 1)
    cv2.imwrite(f"{name}-markers-qa.png", im)
w = csv.DictWriter(open("fig4-markers.csv", "w", newline=""), fieldnames=list(out[0]), lineterminator="\n"); w.writeheader(); w.writerows(out)
for r in out: print(r["panel"], round(r["x"], 4), round(r["y"], 4), r["area_px"], r["merged"])
