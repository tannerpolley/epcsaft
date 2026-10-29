# Frozen pure-solvent inputs

The accepted design and readiness review are retained by the lead at
`build/issue-902/design-p3-r2.md` and `build/issue-902/readiness-p3-r2.md` in the
main checkout. Issue #916 owns this bounded study. All 19 source hashes in the
supplied `p3-r2-source-hashes.json` were verified before input transcription.
`input-hashes.json` binds the local calculation snapshot and copied source text.

`parameters.tsv` has no header: component label, CAS, molar mass (kg/mol),
segment count, diameter (Å), dispersion energy/k (K), association energy/k (K),
association volume, dipole moment (D), and dipole count. The source tuple is
fixed; `MIBK-mass` changes only molar mass. `esper-parameters.json` retains the
two original component records and the complete packet's SHA-256, including its
source/license records. It is a source excerpt, not a materialized shared packet.
The executable consumes the explicit TSV declarations, not the full packet.

TOP is assigned Yu's tuple (Table 2, p.4) to CAS 78-42-2, with the NIST mass
0.4346331 kg/mol. This assignment is a modeling assumption, supported by the
branching sentence in Yu §4.1 p.5, not recovered batch identification. Effective
2B means two sites of multiplicity one, with only unlike association. No
residual dipole, quadrupole, ionic contribution or dielectric model is active;
Yu's relative permittivity 11 remains inactive source information. This is not
proof of molecular donor self-association. TOP parameter-fit observations are
unknown, so Fang's fitting independence is unknown.

MIBK uses Esper component 875 with its inherited mass 0.100089 kg/mol and the
Gross–Vrabec dipole term. The inherited inactive acceptor site is deliberately
removed in this isolated pure-fluid declaration: no association edges means
zero association in either declaration. No cross-association is inferred.
The alternate NIST mass 0.1001589 kg/mol is a separate sensitivity, not a packet
correction. Dodecane is Esper component 809 with its inherited mass and zero
polar contribution. Pure dodecane is an assumed neutral carrier composition;
every dependent prediction is **dependent on the n-dodecane stand-in**.
It does not represent characterized technical or sulfonated kerosene.

`densities.tsv` has no header: row, component, realization, T/K, P/Pa,
observed density/(kg/m³), comparison half-width/(kg/m³). Zero observations
mark stand-in predictions, for which error and measurement agreement are
unavailable. F1–F3 copy Fang Table 1 p.2718, DOI 10.1021/je8006138. Primary
nominal pressure is 100000 Pa; ±0.1 kg/m³ is a conservative upper bound on
unspecified-coverage uncertainty. Bath control is ±0.01 K, not a calibrated
scale uncertainty. ThermoML's 101000 Pa and expanded ±0.4 kg/m³ are separately
compared. Actual barometric pressure is unknown. M1 copies only Asadzadeh's
“This work” Table 2 p.3567, DOI 10.1021/acs.jced.3c00309: 796 kg/m³ at
298.15 K and 101000 Pa. ±2 kg/m³ is a chosen twice-standard-uncertainty
comparison, not an author-reported 95% interval. Separate temperature and
pressure sensitivities are ±0.05 K and ±10000 Pa; covariance is unknown.

`vapors.tsv` has no header: row, source, source temperature, source pressure.
FV rows copy Fuge Table III p.1014, DOI 10.1021/j150500a022, measured pressure
column only, in °C and mmHg. Native conversion uses +273.15 K and
133.322387415 Pa/mmHg. AV rows copy Ambrose Table I p.768, DOI
10.1016/0021-9614(88)90030-4, in IPTS-68 kelvin and kPa, including the visually
corrected 309.676 K and 7.638 kPa entries. No temperature-scale correction is
invented. Neither source has a recovered row-level uncertainty budget.
Membership in Esper's fitted vapor rows is unresolved for both sources.
The copied fitting summary preserves the disagreement between SI table and CSV
about density classification. Published fitting discrepancies are not source
measurement uncertainties. No row is held out on the available evidence.

Okeme's DOI 10.1021/acs.jced.9b01126 remains report-only: calibrated
log10(p/Pa) = −4.52 ±0.19 (standard prediction uncertainty) at nominal 298 K,
uncalibrated −4.53, and vaporization enthalpy 106 ±3.7 kJ/mol from GC retention
at 383–413 K. ThermoML's 298.15 K, 3.02e−5 ±2.81e−5 Pa expanded estimate and
106 ±7 kJ/mol expanded rendering describe the same estimates. They are not
independent direct coexistence observations, and they do not qualify TOP
volatility at a five-percent tolerance.
