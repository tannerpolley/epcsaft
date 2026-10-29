# Sources and input transformations

- Pinho and Macedo (2005), *Solubility of NaCl, NaBr, and KCl in Water, Methanol,
  Ethanol, and Their Mixed Solvents*, Table 1, p.30, DOI
  [10.1021/je049922y](https://doi.org/10.1021/je049922y).
  The eight NaCl–water rows and mass-fraction uncertainties come from
  [NIST ThermoML dataset 1](https://trc.nist.gov/ThermoML/10.1021/je049922y.html).
  NIST identifies the uncertainty evaluator as the data-file compiler and the
  confidence level as 95%. The authors' temperature precision is ±0.1 K;
  their replicate-scatter statement is not the compiler's expanded uncertainty.
- Clarke and Glew (1985), *Evaluation of the Thermodynamic Functions for Aqueous
  Sodium Chloride from Equilibrium and Calorimetric Measurements below 154 °C*,
  [Table 34, p.608](https://doi.org/10.1063/1.555730).
  Each retained row contains ΔsolG°, ΔsolH° and ΔsolCp° and the printed standard
  error of each quantity. Sections 2.1 and 2.8, pp.493–495, define the molal
  standard and dissolution relation. No covariance matrix was recovered.
  Table 17A is part of the same joint evaluation and is not independent validation.
- Held et al. (2014), *ePC-SAFT revised*,
  [Tables 1–3 and Table 2 footnote a](https://doi.org/10.1016/j.cherd.2014.05.017).
  `data/packets/held-2014-strategy2/2` is the selected screening convention:
  closest approach 0.88σ, inferred by the retained source-curve reproduction.
  `/1` uses σ, following the printed screening statement. Packing uses 0.88σ
  in both. Eq. (21) does not by itself establish the executable screening choice.
  Sections 2.4 and 3.1.1 describe the density/osmotic calibration, generally at
  298.15 K up to 5 mol/kg, and the sodium–water temperature law motivated at
  273.15 and 298.15 K. That law is restored natively as
  `k_Na-water = 2.37999 - 0.007981 T/K`.
  The retained water diameter law uses the negative final exponent from
  Held (2008), Eq. (7), correcting the later table's printed sign.
- Malmberg and Maryott (1956), *Dielectric Constant of Water from 0° to 100° C*,
  [Eq. (5), p.5](https://doi.org/10.6028/jres.056.001):
  `ε = 87.740 - 0.40008 t + 0.0009398 t² - 0.000001410 t³`, `t = T/K - 273.15`.
  This native temperature correlation is supplied without a fitted offset.
  It is not identified as the original Held dielectric fit law.
- [NIST sodium chloride identity](https://webbook.nist.gov/cgi/cbook.cgi?ID=C7647145):
  molar mass 58.44277 g/mol; the molality conversion uses 0.05844277 kg/mol.
  The packet's NIST water identity gives 0.01801528 kg/mol.
- Halite volume 27.1 cm³/mol follows PHREEQC `pitzer.dat`, as cited by the
  existing `engine/tests/equilibrium/test_pure_solid.py` reference case.
  The source pressure is 101325 Pa; the study uses 101300 Pa.

The exact source-file hashes are in `source-hashes.json`.
The raw source PDFs and NIST JSON remain in the main checkout's `build/issue-858/`;
this analysis retains their selected numerical rows without duplicating the PDFs.
Each derived parameter snapshot changes only sodium–water temperature dependence
and water permittivity from its corresponding packet.
Its source, identity and locator records describe those changes explicitly.

## NaCl–KCl–water co-saturation

`co-saturation.json` owns the selected Farelo row, raw thermochemical inputs,
solid volumes, restored potassium–water law, and source-file hashes.
Farelo (2006), *Solid-Liquid Equilibria in the Quaternary System
NaCl–KCl–NH4Cl–H2O*, Polish J. Chem. 80, 161–172:
[primary institutional PDF](https://scholar.tecnico.ulisboa.pt/api/records/1abe35b8-c5cc-4b86-b1ca-7da06ebb0493/file/1aaf85bd434bce3eddbd18a819f6ab4b20f551895ae506cbf74881c2f24b38af.pdf),
Table 1 p.166, left block, first 293.2 K row, zero NH4Cl, gives 5.3 mol/kg
NaCl and 1.75 mol/kg KCl. Experimental p.165 gives ±5.2% and ±1% salt-content
uncertainties, respectively, at unspecified confidence. No numerical pressure
is given. The declared 101325 Pa and separate 100000 Pa sensitivity are assumptions.
The retained source inspection was checked against printed pp.165–166.

Clarke–Glew Table 34 p.608, 20 °C, supplies NaCl thermochemistry at 101325 Pa.
The local constant-Cp law maps it from 293.15 to 293.2 K. Its Gibbs-energy input
half-width is the sum of mapped individual standard errors and 0.051 J/mol
rounding. This does not infer joint probabilistic coverage.

Archer (1999), *Thermodynamic Properties of the KCl+H2O System*,
[DOI 10.1063/1.556034](https://doi.org/10.1063/1.556034), provides:
Eqs.(1)–(3) pp.2–3 and Table 1 p.3 for the crystal spline; Table 2 p.4 for
crystal entropy; Eq.(17) p.5 for the dissolution convention; Eq.(18) p.7 and
Tables 4 and 6 p.8 for the thermal mapping at 100000 Pa. The quintic aqueous
Hermite interpolant uses printed Gibbs increments, entropy and Cp. Its first
and second derivatives are −S and −Cp/T. The crystal Gibbs increment integrates
the published spline with lower derivative −0.25 and natural upper end.
The fitted Table 4 Gibbs anchor −5181.81 J/mol is used consistently rather than
the rounded p.15 summary −5184.1 J/mol. Table 4's 95% Gibbs half-width is
10.5 J/mol; p.15 explains the 0.042 J/(mol K) entropy uncertainty. The empirical
2 J/mol interpolation allowance covers the observed 0.41022 J/mol spread of
the inspected alternative mappings; it is not a proven error bound.

Held Table 2 footnote b p.2889 supplies k(K,water)=1.3959−0.004012 T/K.
The existing shared parameter input supplies all other quantities and retains
the sodium–water law. No reusable packet is changed. PHREEQC
[`pitzer.dat`](https://github.com/usgs-coupled/phreeqc3/blob/master/database/pitzer.dat)
Halite lines 336–340 and Sylvite lines 433–437 supply 27.1 and 37.5 cm³/mol.
Native reference transfer retains both source pressures and applies solid-volume
pressure work exactly once. These fixed volumes have no retained uncertainties.
The source PDFs, inspected images, interpolation alternatives and downloaded
database remain in the main checkout's `build/issue-878/`; exact selected source
hashes are retained in the input file.

The declared seven-solve checks require: classified residual ≤1e−10;
absolute component balance and charge error ≤1e−10 mol; absolute saturation
index ≤1e−9 decade; independent native-property affinity contraction error
<1e−8 natural-log unit; doubled-feed molality change <1e−8 relative.
These verify the declared branch separately from comparison with measurements.
