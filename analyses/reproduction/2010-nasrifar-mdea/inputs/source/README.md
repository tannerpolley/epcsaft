# Source values and figure readings

Nasrifar and Tafazzol (2010), DOI 10.1021/ie901181n, Zotero CBRPAXQG,
PDF 3G4FGGY4, is the molecular PC-SAFT source. Tables retain printed decimal
strings and verified PDF/printed-page locators. The parameter packet owns
the numerical declarations; these CSVs retain source transcription only.
Table 2 footnote b reports saturated liquid **volume** deviations for water,
CO2 and methane, and density deviations for MDEA and H2S.

`figure-calibration.json` records manual picks on 600 dpi crops, tick centers,
linear or log10 transforms, source/crop hashes and independent second readings.
`figure-targets.csv` contains the fixed rectangles: three pixels on each axis,
with selected segments at least six pixels from another candidate series.
Bounds, marker centers, curve ownership and calibration were visually checked
on `build/issue-734/figure-*-qa.png` before the target-freeze commit and native
predictions. The archived build directory retains the original crops and QA.
No automatic curve tracer or interpolated target was used. The XY helper's
threshold extraction was unsuitable for these dotted curves;
`scripts/digitize.py` maps the retained manual centers directly.

Figure 6 lower branches are bubbles and upper branches dews. The model uses
the printed pressures. Figure 6c's measured markers retain the separate
66.7 kPa primary-source condition and 67.7 kPa figure label; their pressure
ownership is unresolved without Voutsas et al. (2004), DOI
10.1016/j.fluid.2004.05.012. Figure 7 retains its printed `1000*x(MDEA)` axis
and methane-solubility caption as conflicting information. No solubility
scoring or guessed methane coordinate is admitted without Jou and Mather
(2006), DOI 10.1021/je060118g.

Figure 15 isolated water, carbonate, hydroxide, hydronium, free CO2 and MDEA
curves are numerical-verification targets, with log10 mole fraction as the
stored vertical coordinate. The overlapping bicarbonate/protonated-MDEA
curve, crossings and endpoints are **not digitized: ambiguous ownership**.
Figure 16 and Figure 17 model curves are **not digitized: ambiguous ownership**;
same-style curves lack independent condition labels. No curve is assigned
by proximity to a measurement or an Engine prediction.

Kamps et al. (2001), DOI 10.1021/ie000441r, Zotero 5R8DLM23 / WA2BMMU6,
supplies Tables 2–3. Its library filename says 2000; publication is 2001.
Rho et al. (1997), DOI 10.1021/je970097d, Zotero 9PPV4M2Z / GAN4RCX5,
supplies the selected Tables 1 and 3 rows. Every numeric tuple was checked
against the exact PDF page during design. These are physical-comparison
targets, not the paper's predicted curves. Rho's measured 20.5 mass % at
373.15 K stays separate from Figure 17b's 20% model label.

Pressure agreement uses half the last printed digit, not experimental
uncertainty. Kamps pressure column `10p/MPa` converts by 100000 Pa per unit;
Rho pressure converts from kPa. Kamps reports pressure uncertainty 0.1% of
transducer full scale (0.6, 4 or 10 MPa; row identity unspecified), temperature
±0.1 K, gas mass ±0.008 g, solvent mass at most 0.7% and MDEA molality ±0.1%.
Printed Kamps temperature resolution is ±0.005 K. Rho reports ±1.5% GC
response and ±5% reproducibility, not a complete row pressure uncertainty.

Reactive predictions are unavailable pending Engine #776. The printed
hydronium reactions and Eq. (12) water balance are inconsistent; this is a
source inconsistency, distinct from the missing ideal chemistry and native
effective-water conversion. No Python chemistry or corrected source model
is substituted. Pure-MDEA Figures 2–3 and Table 2 aggregate errors remain
source-only outside the original case inventory; their complete fitting
grids and membership are unavailable.
