# Gross–Sadowski 2001 catalog review

## Source

- Gross, J.; Sadowski, G. “PC-SAFT: An Equation of State Based on a
  Perturbation Theory for Chain Molecules.” *Industrial & Engineering
  Chemistry Research* **2001**, 40, 1244–1260.
- DOI: `10.1021/ie0003887`
- Zotero collection: `ePC-SAFT Progression`
- Inspected PDF SHA-256:
  `373b1c6224407933ee8e123bdd78b186b8bf49cabe1f17dffc9bd631d660bc72`

## Admitted records

The source audit found 78 model-parameter rows in Table 2. Every row is
represented by one component identity, one source fit domain, and the four
printed inputs: molar mass, segment count, segment diameter, and dispersion
energy over Boltzmann’s constant. Values are stored at the paper’s printed
precision; molar masses are converted exactly from grams per mole to kilograms
per mole.

The bundle also includes:

- all 24 PC-SAFT binary interaction parameters in Table 4, including the
  printed temperature ranges converted from degrees Celsius to kelvin;
- the explicit methane–ethane zero and carbon-dioxide–ethane value from the
  ternary-system discussion;
- the carbon-dioxide–decane and carbon-dioxide–cyclohexane values from Figures
  12 and 13; and
- the 1-chlorobutane–hexane value from Figure 15.

Blank table cells and unreported pairs are absent. Table 1’s universal
dispersion coefficients remain equation constants rather than parameter-bundle
records.

## Appendix B disposition

Appendix B reports three n-alkane correlations parameterized by molar mass.
The current typed correlation contract accepts temperature as its independent
variable and does not represent the Appendix B outputs
`sigma(M)`, `m(M) / M`, or `epsilon(M) / k`. Storing those coefficients would
therefore require a special-purpose schema and resolver path. They are not
runtime records in this bundle; the 20 printed n-alkane rows from Table 2 are
included directly.

## License, use basis, and tolerance

The inspected article is publisher-copyrighted source material. The package
retains factual numeric parameter values, bibliographic attribution, and short
table/figure locators; it does not redistribute article prose, tables, figures,
or the source PDF. The use basis is exact, source-traceable implementation of
the paper’s model parameterization.

Source-record acceptance has zero numeric transcription tolerance: decimal
values and explicit zeroes must equal the printed digits exactly. The only unit
transformation is molar mass from grams per mole to kilograms per mole, applied
as exact decimal division by 1000. Table 4 temperature endpoints are converted
from degrees Celsius to kelvin by exact decimal addition of 273.15. Bundle and
selection fingerprints are exact. The component-order EOS comparison uses
relative tolerance `1e-13` and absolute tolerance `1e-12 J/mol`, sufficient for
floating-point evaluation-order noise while materially smaller than the
reported parameter precision.

## Verification contract

The package test independently transcribes all 78 Table 2 rows and checks every
component identity, fit domain, unit, and pure-parameter value. It also fixes
all 29 admitted binary values, the source DOI, the bundle fingerprint, three
selection fingerprints, a newly admitted eicosane pure evaluation, and a
carbon-dioxide/cyclohexane component-order permutation through the public EOS
path.
