# NaCl dilution enthalpy: source notes

The measurements originate with Messikomer and Wood (1975), *J. Chem. Thermodyn.* 7, 119–130, DOI [10.1016/0021-9614(75)90259-1](https://doi.org/10.1016/0021-9614(75)90259-1).
The original numerical tables and full uncertainty documentation were unavailable.
The calculation uses their measured values reproduced by Partanen and Partanen (2023), *J. Solution Chem.* 52, 1352–1385, DOI [10.1007/s10953-023-01322-y](https://doi.org/10.1007/s10953-023-01322-y), supplementary Appendix C, Tables S19–S21, printed pages 21–23.
The [retained source pages](sources/partanen2023-tables-s19-s21.txt) contain the measured initial molality, final molality and dilution enthalpy in the first three columns.
Adjacent correlation errors and adiabatic temperature calculations are excluded.
The [CSV](measurement-candidates.csv) retains 60 rows, their original source lines, table and row numbers, inclusion and half-last-printed-digit radii.
The numeric CSV fields omit some trailing zeros; the unchanged `source_line` retains their printed strings, which the producer independently checks against these source pages.
The primary population contains 54 rows and 24 temperature-specific chords.
The six rows beginning at 5.008 mol/kg remain retained but excluded because they exceed Held's up-to-5 mol/kg fit domain.
No selection uses calculated errors.

Pitzer, Peiper and Busey (1984), *JPCRD* 13, 1–102, DOI [10.1063/1.555709](https://doi.org/10.1063/1.555709), Sec. 2, Eq. (17), printed page 4, gives final minus initial apparent molar enthalpy for dilution.
[Retained pages](sources/pitzer1984-pages4-6.txt) include that relation, Sec. 3's precision discussion and Table 1, reference 38.
The tabulated dilution precision coefficient is 0.003, while reproduction by their fitted equation is separately 0.004.
The published concentration-dependent precision prescription does not establish which concentration controls a finite chord or the extension below 0.1 mol/kg here.
This study therefore chooses the uniform benchmark 0.003RT, without assigning it measurement-uncertainty coverage.
Printed enthalpy rounding is added separately.
Missing temperature, pressure, preparation and systematic uncertainties remain unknown.
Replicates share inputs and possibly systematic errors; the benchmark is never divided by their square root.

The parameter owner is `data/packets/held-2014-strategy2/2/parameters/parameters.json`.
The [derived input snapshot](parameters.json) restores the published sodium–water temperature law and declares Malmberg–Maryott solvent permittivity.
Held et al. (2014), DOI [10.1016/j.cherd.2014.05.017](https://doi.org/10.1016/j.cherd.2014.05.017), Sec. 2.4, Eqs. (19)–(20), pages 2887–2889 and Tables 1–3, reports density-difference and osmotic objectives.
Its Fig. 2 and discussion, pages 2890–2891, uses osmotic behavior at 273.15 and 298.15 K to determine the sodium–water law.
No direct caloric rows are documented in those objectives.
[Held evidence](sources/held2014-evidence.json) retains exact source passages.

The water diameter follows Held, Cameretti and Sadowski (2008), DOI [10.1016/j.fluid.2008.06.010](https://doi.org/10.1016/j.fluid.2008.06.010), Eq. (7) and Tables 1–2, pages 89–92.
Use the negative final exponent, rather than the positive-exponent misprint in Held (2014).
Its ion fit uses density and mean ionic activity coefficients.
[Held water evidence](sources/held2008-evidence.json) and [Cameretti lineage evidence](sources/cameretti2005-evidence.json), DOI [10.1021/ie0488142](https://doi.org/10.1021/ie0488142), Sec. 3, pages 3358–3359, document the water vapor-pressure/density lineage.
The original Fuchs estimation ledger and any calorimetric ancestry of compiled osmotic inputs were not fully recovered.
Thus caloric overlap is only partly audited: the comparison is withheld relative to published objectives, without a complete historical independence claim.

The packet's packing and screening diameters are 0.88σ for both ions.
This screening choice is a repository numerical inference from reproduction of source curves, not a recovered author executable.
Debye–Hückel is included; Born and same-sign ionic dispersion are excluded.
The unlike-ion dispersion and water–ion interactions remain declared.

Malmberg and Maryott (1956), DOI [10.6028/jres.056.001](https://doi.org/10.6028/jres.056.001), Eq. (5), printed page 5, supplies the salt-independent cubic permittivity law.
The [retained transcription](sources/malmberg-maryott-1956.txt) includes its accuracy discussion on page 6.
This law is consistent with the salt-independent solvent assumption, but its identity with Held's fitting law is unverified.
No coefficient or derivative covariance is supplied.
The higher study temperatures extrapolate beyond the two temperatures documented for the sodium–water law.

The thermal records are deliberately manufactured mathematical inputs.
Their isothermal differences cancel by conserved ion inventories and the added pure-water balance.
No gas-ion standard-state accuracy, absolute apparent heat capacity, salt crystal or solution enthalpy is inferred.
The numerical source checks use no new caloric fit, offset, empirical correction or parameter selection.
Optional screening/permittivity sensitivities and direct measured heat capacities are not run.
Tanner–Lamb original heat-capacity tables were unavailable; an abstract or evaluated correlation is not substituted for measurements.
The cumulative-titration ThermoML candidate is excluded because its normalization was not verified.

[Source records](source-records.json) identify retained file hashes and their full-source parents.
The supplement was acquired through its [open publisher supplement](https://media.springernature.com/original/springer-static/esm/art%3A10.1007%2Fs10953-023-01322-y/MediaObjects/10953_2023_1322_MOESM1_ESM.pdf).
The original PPB PDF is the [NIST reprint](https://srd.nist.gov/jpcrdreprint/1.555709.pdf).
All source evidence was copied from the accepted investigation; no browser or new source acquisition was used for this calculation.
