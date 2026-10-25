# LiBr model curve digitized from Figiel, Yu & Held (2025), Figure 5 (right panel)

**This file records a digitization of a published plot. It is not a measurement, it is not an
experimental data set, and it is not an oracle row.** Every number in
`libr-model-digitized.csv` was read off a raster image of a printed figure;
`engine/tests/eos/test_activity.py` compares the model against it as a curve reproduction. Treat the values as a curve traced from ink, with the stated per-point uncertainty, and never as a
published tabulation.

## What was digitized

| Item | Value |
| --- | --- |
| Paper | Figiel, P.; Yu, G.; Held, C. *Predicting Thermodynamic Properties of Ions in Single Solvents and in Mixed Solvents Using a Modified Born Term within the ePC-SAFT Framework.* Ind. Eng. Chem. Res. 2025, 64, 9406-9418 |
| Figure | Figure 5, page 7 of the article PDF |
| Panel | **Right panel** (Br- salts). The left panel is the Cl- salts and was not used |
| Series | **Orange line = LiBr model curve**, "Calculated using ePC-SAFT with the modeling approach introduced in this work using parameters from Tables 1-5" (Figure 5 caption) |
| System | LiBr in water, 298.15 K, 1 bar |
| Source PDF | `~/Zotero/storage/WE6I8RN3/Figiel et al. - 2025 - ... Modifie.pdf`, sha256 `a3c940895c530f72f47f22f8c0f4796b5ad37918a4edd274f85497c6dd81ad1f` |
| Figure raster | Embedded image object 6 on PDF page 7, extracted losslessly with `pdfimages -f 7 -l 7 -png`; 1500 x 609 px, 300 ppi, sha256 `85633c445343501360289fb1aa7084f5d5fbe2f0ace662f80a50159e08784c38` |
| Transcription cross-check | `~/Zotero/storage/2M74BKJG/Figiel, Yu, Held - 2025 - ....md` (Figure 5 caption and the "MIACS IN WATER" section) |

The caption assigns colors: purple down-triangles = H+ salts, orange up-triangles = Li+ salts,
green squares = Na+ salts, gray circles = K+ salts. Symbols are literature data (Hamer & Wu 1972,
ref 71); lines are the ePC-SAFT model of this work. Only the **line** was digitized.

## Which activity-coefficient scale

**Molality scale.** Both the caption ("Molality-based MIACs of different salts in water at
298.15 K and 1 bar") and the ordinate label (gamma with superscript `*,m` and subscript `+-`)
say so, and the abscissa is `m~ / mol*kg^-1`. The `*` denotes the infinite-dilution
(unsymmetric) reference state, which is the usual convention for ions and is what eq 20-21 of the
paper define.

So the CSV column `gamma_pm_m` is the **molality-scale** mean ionic activity coefficient. To get
the mole-fraction-scale value, multiply by `1 + M_water * m * sum(nu)` with `M_water = 0.018015
kg/mol` and `sum(nu) = 2` for LiBr. This is the same convention as the sibling experimental file
`water-libr.csv` in this directory, whose `miac_m` column is the molality scale and whose `miac`
column is the mole-fraction scale; the two differ there by exactly this factor (checked at
m = 6: 3.915 * (1 + 0.018015 * 6 * 2) = 4.761, matching that file's `miac` value). Confusing the
two scales would corrupt every comparison, so any consumer of this CSV must use the molality
scale.

## Axes

Both axes are **linear**. There is no logarithmic axis anywhere in Figure 5.

| Axis | Range shown | Labelled ticks |
| --- | --- | --- |
| x = m~ / mol*kg^-1 | 0 to 6 | 0, 3, 6 (minor ticks at 1, 2, 4, 5) |
| y = gamma_+-^(*,m) / - | 0 to 4.507 at the top frame | 0, 1, 2, 3, 4 |

## Pixel anchor points (right panel, in the extracted 1500 x 609 raster)

Tick positions were located automatically as the dark runs on the panel frame, not by eye.

| Data value | Pixel coordinate (center of the tick/frame run) |
| --- | --- |
| m = 0 | column 901.0 (left frame of right panel; run 900-902) |
| m = 3 | column 1196.0 (run 1195-1197) |
| m = 6 | column 1490.0 (right frame; run 1489-1491) |
| gamma = 4 | row 95.5 (run 94-97) |
| gamma = 3 | row 195.0 (run 194-196) |
| gamma = 2 | row 294.5 (run 293-296) |
| gamma = 0 | row 493.5 (bottom frame; run 492-495) |

Resulting calibration:

```
m     = (column - 901.0) / 98.1667     [98.1667 px per mol/kg = 589 px / 6 mol/kg]
gamma = (493.5 - row)   / 99.5         [99.5 px per gamma unit]
```

The gamma = 1 tick is overprinted by the marker pile and was not used; it is not needed, because a
straight-line fit through the four usable y anchors is exact to better than 0.02 px (spacings
99.5, 99.5, 199.0 px). The three x anchors are consistent to 0.5 px (spacings 295 and 294 px).
The left panel uses the identical scale (columns 125 / 420 / 714) and was checked but not used.

## Method

1. Extract the figure raster losslessly from the PDF (`pdfimages -png`). The publisher's embedded
   image is the highest resolution available; rendering the page at a higher DPI would only
   upsample it.
2. Locate the panel frames and ticks from columns/rows of dark pixels; derive the linear
   calibration above.
3. Build a color mask for the LiBr series, `|RGB - (255,132,0)| (L1) < 120`. The figure uses a
   small, saturated palette - purple (132,0,255), green (0,132,0), orange (255,132,0), gray
   (198,198,198) - so the series separate cleanly with essentially no color ambiguity.
4. For every image column, collect vertical runs of orange pixels. A column is accepted only if it
   contains **exactly one** run of height <= 6 px (the drawn line is 2-4 px tall; the LiBr
   triangles are 18-23 px tall) **and** that thin run is at least 3 px clear of any thick run in
   the same column. Columns where the line runs into or behind a marker therefore drop out
   automatically rather than being guessed at.
5. Reject residual false positives (chiefly the thin apexes of the orange triangles, which mimic a
   line) with an iterative degree-6 polynomial fit over the accepted centers, discarding
   points more than 2 px from the fit. The fit was used **only** as an outlier filter; every value
   written to the CSV is a measured pixel position, not a fitted value.
6. Report a point as the median of the run centers over a +-3 column window (7 columns), converted
   through the calibration above.
7. Verify by drawing the 27 written points back onto the original raster and inspecting the
   overlay; all 27 crosses sit on the orange line.

## Independent calibration check

The orange triangles are Hamer & Wu literature points at round molalities, and this repository
already holds those values in `water-libr.csv`. Comparing the measured marker positions to the
tabulated values tests the whole calibration chain end to end:

| m | tabulated gamma_+-^(m) | marker bounding-box center, this calibration | difference |
| --- | --- | --- | --- |
| 2.5 | 1.166 | 1.206 | +0.040 |
| 3.0 | 1.359 | 1.387 | +0.028 |
| 3.5 | 1.598 | 1.638 | +0.040 |
| 4.0 | 1.891 | 1.930 | +0.039 |
| 5.0 | 2.696 | 2.729 | +0.033 |
| 6.0 | 3.915 | 3.945 | +0.030 |

The offset is a constant +0.035 +- 0.006, i.e. +3.5 +- 0.6 px, which is exactly `h/6` for markers
of height h ~ 21 px: the triangles are drawn with their **area centroid** on the data point, so
their bounding-box center sits h/6 above it. After allowing for that marker geometry the residual
scatter is +-0.6 px, which bounds the calibration error well below 1 px. (The m = 2.0 marker is
partly hidden behind a purple triangle and its blob is clipped, so it is excluded from this table.)

## Uncertainty

`digitization_uncertainty` is an absolute, 1-sigma-equivalent uncertainty in gamma, computed
**per point** from that point's own pixel evidence:

```
sigma_gamma = sqrt(  (w/4 / 99.5)^2                 line-center from the local stroke height w
                   + (0.25 / sqrt(7) / 99.5)^2      averaged column-to-column trace scatter
                   + (0.5 / 99.5)^2                 y-axis anchor (tick center) location
                   + (0.0013 * gamma)^2             y-axis scale error, 0.5 px over 398 px
                   + (dgamma/dm * 0.5 / 98.1667)^2  x-axis anchor error mapped through the slope
                  )
```

Terms, and where each number comes from:

- **Line thickness.** `w` is the median height of the orange run in the 7-column window, measured
  per point: 2 px at one point, 3 px over most of the curve, 4 px at the steep high-molality end
  (the vertical run of a sloped stroke is longer than its perpendicular width). The center of a
  thresholded run estimates the stroke center with a residual bias bounded by about a quarter of
  the run height, hence `w/4`: 0.5 to 1.0 px.
- **Trace scatter.** The measured column-to-column RMS deviation of the accepted centers about a
  local quadratic is 0.267 px (n = 271 columns), damped by sqrt(7) through the windowed median.
- **Axis anchors.** Tick and frame runs are 3-4 px wide, so their centers are good to about
  +-0.5 px. The scale term follows from the same 0.5 px over the 398 px between the gamma = 0 and
  gamma = 4 anchors. The x-anchor term converts a 0.5 px horizontal error into a gamma error
  through the local slope (0.28 to 0.58 per mol/kg over the digitized range).

Result: **0.008 to 0.012 in gamma**, i.e. 0.8 % relative at the low-molality end and 0.45 % at the
high end. One point (m = 3.973) sits on a 2 px stroke and gets 0.0078. The uncertainty is
deliberately not a constant, because the stroke height, the local slope and the scale term all
vary along the curve.

Molality uncertainty is not tabulated: the abscissa of each point is a chosen pixel column, so it
is exact to +-0.5 px = +-0.005 mol/kg, and that error is already folded into the gamma budget
through the slope term.

## Coverage, and what was deliberately skipped

- **Written range: m = 2.160 to 5.847 mol/kg, 27 points.** Inside a clean window the points are
  0.081 to 0.143 mol/kg apart; consecutive windows are separated by 0.12 to 0.33 mol/kg where a
  marker covers the line. The spacing is set by those gaps, not chosen freely.
- **m < 2.1 was skipped entirely.** Below about 2 mol/kg the LiBr line is completely buried under
  the overlapping pile of purple H+ triangles and orange Li+ triangles at 0.5 mol/kg intervals (see
  the low-molality region of the right panel). A few isolated orange fragments survive the
  automatic filter down there, but visual inspection shows they are triangle apexes and slivers of
  partly-occluded markers, not the line. They were discarded. **The minimum of the LiBr curve, and
  its entire initial decrease, are therefore not in this file.**
- **The last 3 columns before the right frame (m > 5.88) were dropped** so that the black frame
  cannot contaminate the run. The curve visually reaches about gamma = 2.75 at m = 6, but that
  value is not written here.
- Gaps at every 0.5 mol/kg (each orange marker) and where the line passes behind a purple triangle
  are simply absent from the CSV. Nothing was interpolated across them.

## Judgement calls

1. **Which curve.** The single cleanest and widest curve in the whole figure is the purple **HBr**
   line in the right panel: it climbs clear of the marker pile at about m = 0.9 and is traceable
   over 0.88 to 6.0 mol/kg (416 of 590 clean columns). It was rejected because H+ is not in this
   repository's Figiel parameter bundle (`inputs/parameters/mixed/binary-interaction/k-ij.csv`
   carries Li+, Na+, K+, Cl-, Br-, I- only, and H+ is flagged `later_capability_gate` in the
   regression-target ledger), so an HBr curve could not be compared against anything today.
   Among the salts the engine can actually evaluate, LiBr was chosen over NaBr and the two chloride
   curves because (a) its line runs through open white space between the purple and green curves,
   whereas the NaBr line is tangent to the tops of the green squares along its whole length, and
   (b) it spans gamma = 1.16 to 2.69, a far larger dynamic range than NaBr's 0.85 to 1.31
  over the same usable molality window, so a
   comparison against it is a sharper test. LiCl (left panel, orange) has the same shape but fewer
   clean columns (280 vs 324).
2. **Marker-geometry offset.** The +0.035 gamma offset in the check table above is attributed to
   triangle centroid-vs-bounding-box placement, not to a calibration error. This attribution is an
   inference; it is supported by its size matching `h/6` and by being constant across six markers,
   but it is an inference. **No such offset was applied to the digitized line** - it concerns only
   how the markers are drawn, and the line is not a marker.
3. **Local slope of accepted runs.** Where the line is steep (m > 5.4) the vertical run grows to
   4 px. The run center is still the best estimate of the stroke center for a symmetric stroke, and
   the widened run is charged to the uncertainty through the `w/4` term rather than corrected.
4. **The isolated clean window at m = 3.92 to 4.02** is only 11 columns wide because a purple
   triangle crosses the LiBr line there. One point (m = 3.973) was taken from its middle. It is the
   point most exposed to contamination from the neighbouring purple marker, though its narrow 2 px
   run and its agreement with its neighbours argue that it is clean.
5. **Color threshold.** An L1 distance of 120 to the pure series color was used. The palette is
   saturated and well separated, so the result is insensitive to this choice; tightening it only
   thins the anti-aliased edges of the same runs.
6. The reported values are medians over 7 columns, so a written point is a local average of the
   line over about 0.07 mol/kg, not a single-pixel reading.

## Reproducing

Extract image object 6 from page 7 of the PDF, apply the calibration above, and follow the
thin-run/thick-run rule in the Method section. No digitizing GUI was used and no point was placed
by hand or by eye.
