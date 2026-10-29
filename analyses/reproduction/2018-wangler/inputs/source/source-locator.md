# Wangler et al. 2018 source locator

- Source: `/home/tnnrpolley21/Zotero/storage/RIUUBS2B/Wangler et al.md`
- SHA-256: `f60c6ae046c71ec1d5c853d85f80382078cdf40e5d4adc609840d799ac6c1da0`
- DOI: `10.1016/j.fluid.2017.12.033`
- Parameter basis: Tables 2 and 4–7; Eqs. 3, 8, 13–21, and 34–35.
- Independent aqueous MDEA/HCl fit: Tables 8–9 and Figure 10.
- Reactive MDEA/CO2 comparison: Table 10 and Figure 14.

The campaign retains the three exact printed Table 10 observations and their
reported uncertainties in `wangler-2018-table-10.csv`. The vapor phase and its
neutral species support are declared before the solve; the calculation does
not search for a different phase topology.

## Neutral diagnosis sources

Printed journal pages are used throughout. Wangler's article PDF is Zotero
`G7ZTAJZS`; its Figure 1 is on p. 20 and Eq. 12 on p. 17. PDF and supporting
source hashes are retained in `diagnosis-sources.json`.

`wangler-2018-figure-1.csv` contains four source-curve readouts at x(MDEA) = 0.2,
not experimental observations. The embedded raster was calibrated from x = 0
at pixel 171 to x = 0.4 at 955, and p = 0 at y = 438 to p = 1 bar at 44.
Median ink coordinates at grayscale thresholds 110 and 170 support a
conservative ±0.01 bar allowance. The source raster remains in the local
Zotero PDF; the tracked analysis retains the digitized values and calibration.
No complete figure agreement is inferred from this single composition.

Pahlavanzadeh and Fakouri Baygi, *J. Chem. Thermodynamics* 59 (2013), 214–221,
DOI `10.1016/j.jct.2012.12.021`, Zotero PDF `YAUNIZGX`: p. 216 Table 2 supplies
the original MDEA parameter sets; p. 215 Eq. 4 uses an association d³ prefactor.
The inspected p. 216 Table 1 vapor-pressure row labeled MDEA duplicates Kim's
MEA row; the coefficients are retained in `diagnosis-sources.json` and the
source page remains in the local Zotero PDF. Its reported pressure/density fit
deviations cannot be recomputed without the original fit-point grid. Its
pure p–density plot has no per-point temperatures and cannot supply an
independent temperature-matched density measurement.

Kim, Svendsen and Børresen, *J. Chem. Eng. Data* 53 (2008), 2521–2531,
DOI `10.1021/je800290k`, Zotero PDF `4SMCK4SC`: p. 2522 Table 3 supplies the seven
neat-MDEA measurements in `kim-2008-neat-mdea.csv`; p. 2523 Table 5 supplies the
Riedel and Antoine coefficients retained in `diagnosis-sources.json`.
Riedel is evaluated as ln(p/Pa) = A + B/T + C ln T + D T^E, explicitly correcting
the table footnote's omitted reciprocal. Antoine is log10(p/mmHg) = A − B/(t/°C+C),
with 1 mmHg = 133.32236842105263 Pa. The source's p. 2526 uncertainties are
±0.05 K and ±0.3 kPa for temperature and the pressure controller. The two
correlations have different data bases; low-temperature evaluations extrapolate
beyond these seven new measurements.

Gross and Sadowski (2002), Eqs. 2–3, Zotero `2AJYEZS6`, transcription `FZVTDKXU`,
is the checked owner of the usual Wolbach–Sandler size-corrected association
rule. Both implementations include both unlike donor–acceptor directions,
exclude donor–donor and acceptor–acceptor bonds, and use the same neutral
packing-diameter law and physical constants. FeOs 0.10.1 multiplies its
association volume by (sigma_i sigma_j)^(3/2); the Engine uses sigma_ij³.
The dimensionless volumes are converted so their dimensional volumes coincide.

## Convention table meanings

`baseline` uses the named usual rule and polynomial temperature powers;
`eq12-literal` puts q³(1−l) inside the radical as printed on Wangler p. 17;
`literal-k`/`literal-l` multiply each coefficient by ΔT instead of its
subscript-indicated power (p. 19 Eqs. 34–35). `mdea-origin-digits` uses the
extra digits in Pahlavanzadeh Table 2; `association-diameter-d` rescales every
association edge by (d_ij/sigma_ij)³ per its Eq. 4. These are the printed readings
or parameter-source comparisons; none is adopted as a neutral correction.

Excluded diagnostic controls are geometric cross volume without size correction,
placing only (1−l) under the radical, geometric Mayer-strength cross energy,
constant water diameter, raw-temperature or truncated/sign-reversed k(T),
changed MDEA site counts and induced-only association. In `mdea-4:2` the counts
are four acceptors and two donors; `mdea-2:3` means two donors and three acceptors.
`water-pahlav-2/4` and `mdea-source-2site/6site` are complete alternate sets from
Pahlavanzadeh Table 2; `water-gross-2` is the Gross–Sadowski water set. They are
not Wangler's adopted sets. `induced-retain-cross` removes MDEA self edges only;
`induced-zero-self-energy` also reduces the cross energy. Neither is authorized
by Wangler's statement about protonated MDEAH⁺. The CSV retains failures rather
than assigning a pressure to an unstable state. No simultaneous fitting of
conventions or parameters was performed.
