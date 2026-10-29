# Source declarations and applicability

The fixed input declarations answer the owner-selected pressure question in [#856](https://github.com/tannerpolley/ePC-SAFT/issues/856).
The accepted amendment permits unresolved historical calibration membership and requires the corresponding conservative label.
There is no local parameter fit or measurement replacement.

## Acid dissociation, E1

Harned, H. S.; Ehlers, R. W. *The Dissociation Constant of Acetic Acid from 0 to 60° Centigrade*. Journal of the American Chemical Society **1933**, *55* (2), 652–656. [DOI](https://doi.org/10.1021/ja01329a027).
Zotero parent `2BZC7C6J` was added to collection `YR8SNI5D` through Companion's local `acquire-metadata` operation after inspecting its preview.
The owner's explicit metadata-only instruction superseded the ordinary acquisition command for this operation.
No system browser, Connector save, PDF attachment, conversion or bibliography refresh was used.
The preview and applied response are retained in the main checkout's `build/issue-856/build-r1/`.

The [indexed primary PDF](https://www.electronicsandbooks.com/edt/manual/Magazine/J/Journal%20of%20the%20American%20Chemical%20Society%20US/1933%20%20%28vol%20055%29/02%20%20%28441-864%29/652-656.pdf), p. 652, Eq. (1), identifies molalities and activity coefficients and extrapolates to zero ionic strength.
Its final Table II, p. 655, gives the observed and calculated values at 50 °C as `K × 10^5 = 1.633` in both columns.
The selected row is transcribed in `reactions.json`, including the other table columns to make the row identity checkable.
The web index returned the primary-paper text, but direct PDF retrieval failed with `No route to host` and web fetch timeouts.
Visual verification of the original table page remains unavailable; indexed PDF text is the retained verification basis.

The declaration uses dimensionless molality activities normalized by 1 mol kg⁻¹, infinite dilution in pure water, and formal H⁺.
The source value does not establish measured liquid pressure dependence.
The fixed 101325 Pa source pressure is an adopted convention; a hydrogen electrode gas standard would not establish the liquid's pressure dependence.
The Engine performs the reference transfer to the equilibrium state.

## Historical calibration membership, E2

Ascani et al. **2023**, [DOI](https://doi.org/10.3390/molecules28041768), Table 2, water/acetic-acid row, and reference 107 trace `k = −0.1247` to Pabsch et al., *Influence of solvent and salt on kinetics and equilibrium of esterification reactions*, Chemical Engineering Science **2022**, *263*, 118046, [DOI](https://doi.org/10.1016/j.ces.2022.118046).
Ascani's retained source has Zotero parent/child `YHTZTEAR`/`VGBSV5AU`.
It identifies binary vapor–liquid fitting but does not identify the fitted rows.
The Pabsch abstract concerns esterification measurements at 323 K; it does not enumerate the binary calibration observations.
At the initial comparison, the full paper, supplement and original binary observation inventory could not be inspected through the available nonbrowser source path; publisher access returned 403, and no retained CES parent or companion was found.
The configured Elsevier acquisition path reports an unresolved Article Retrieval key restriction; no browser fallback was used.
The separate Pabsch JCED salt-solubility paper, `4VQ3N88M`/`IZTGW7GI`, is not evidence of membership in the CES binary fit.

Bernatová et al. **2006**, [DOI](https://doi.org/10.1016/j.fluid.2006.06.005), ThermoML groups 13/14, rows 9–11, therefore cannot be shown to be historically withheld.
All tables use **fixed-parameter comparison; historical calibration overlap unresolved**.
The comparison is neither an independently withheld prediction nor a calibration replay of the original model.
Gross–Sadowski pure-component saturation and density fitting, water and ion parameter adoption, and the Pabsch binary fit remain inherited calibration.
No calibration statistics are computed for those unavailable historical datasets.

## Water ionization, E3

[IAPWS R11-24](https://iapws.org/documents/release/Ionization.download), **2024**, p. 2 nomenclature, p. 3 Eq. (1), and p. 6 Table 3 supply the activity basis and `pKw = 13.26` at 50 °C.
Table 3 footnote b uses liquid water at 0.1 MPa below 100 °C.
The retained PDF was inspected, including the reaction page.
The source writes `2 H2O ⇌ H3O+ + OH−`, with molal ionic standards and a molecular-water standard.
The adopted formal convention `H+ := H3O+ − H2O`, or `a_H3O+ = a_H+ a_water`, gives `H2O ⇌ H+ + OH−` with the same scalar constant and one water activity.
This is an explicit convention conversion, not a substitution of hydronium parameters for hydrogen-ion parameters.
The Engine includes the water activity and performs both reaction-reference conversions.
The tabulated constant is retained at its published rounding; its pressure and rounding uncertainties are not assumed zero.
Water ionization is included in the primary calculation, so no omission is needed even though its pressure effect is below the numerical resolution of this comparison.

## Acid permittivity, E4

Wohlfarth, C., *Permittivity (Dielectric Constant) of Liquids*, in the **82nd edition** of the *CRC Handbook of Chemistry and Physics*, Lide, D. R., Ed. (2001), pp. 6-139 and 6-142, supplies the acetic-acid temperature correlation.
The introduction specifies static relative permittivity at nominal atmospheric pressure below the normal boiling point and the law `epsilon_r = a + b T + c T^2 + d T^3` with `T` in kelvin.
The C₂H₄O₂ acetic-acid row gives `a = −15.731`, `b = 0.12662 K−1`, `c = −0.00017738 K−2`, and `d = 0`, over 293–363 K.
The retained PDF page was visually inspected, including signs, powers and the row's temperature interval.
The primary calculation adopts its value at the measurement temperature; it does not reuse a room-temperature constant.
The original water-only dielectric declaration is retained as a sensitivity control, not selected after inspecting pressure accuracy.
Neither the pure-acid correlation nor its adoption validates the water/acid mixture permittivity rule or ion speciation.
The CRC source points to Wohlfarth's Landolt–Börnstein compilation for the original measurements; those measurements were not independently re-evaluated here.

## EOS family and measurements

Ascani et al. **2022**, [DOI](https://doi.org/10.1039/D2CC02943J), ESI Tables 1–3 and 5, pp. 10–12, supplies the adopted acid/ion family; parent `AGPHR9MB`, main PDF `WQ3RRAD4`.
The ESI is not a retained Zotero child; its [author-posted text](https://www.researchgate.net/publication/361665138_Prediction_of_pH_in_multiphase_multicomponent_systems_with_ePC-SAFT_advanced) and the accepted source comparison support the table locators.
Gross–Sadowski **2002**, [DOI](https://doi.org/10.1021/ie010954d), Table 1, independently supplies the acid's 211.59 K dispersion energy; Eqs. (2)–(3) supply the selected cross-association rule; parent/child `2AJYEZS6`/`FZVTDKXU`.
Ascani **2023** Table 1's water-diameter footnote supplies the inherited temperature law; its inconsistent acid energy is not adopted.
Held et al. **2014**, [DOI](https://doi.org/10.1016/j.cherd.2014.05.017), Table 2 and §2.3, supplies ion antecedents and the like-sign dispersion exclusion; parent/child `3DIW87Q5`/`TZJJENHW`.
Ascani–Held **2021**, Eq. (10) and Eq. (13), supplies the uncorrected Born and component-permittivity interpretations; parent/child `KLBK4XDV`/`6WHXVDBI`.
The cross-association choice remains a stated interpretation of the adopted family.

The [NIST ThermoML JSON](https://trc.nist.gov/ThermoML/10.1016/j.fluid.2006.06.005.json) is retained unchanged.
`extract_measurements.py` joins pressure group 13 and liquid-composition group 14 on their shared gas-water fraction.
The measured gas fraction is a join key, not the analytical liquid composition or a molecular vapor partial-pressure datum.
All joined rows are retained with selection or exclusion reasons.
The compiler's combined expanded uncertainties are not standard deviations or independently verified author estimates.
No temperature uncertainty or covariance is supplied in the selected constraint.

## Executable and numerical reference, E5

The exact accepted main revision, wheel file and hash, installed module origins, native library mappings, build-source record and input hashes are in `results/executable.json`.
The main checkout's `build/issue-856/build-r1/` retains the wheel, installed-file verification, source hashes, compiler, build log and check logs.
The study uses an installed wheel without an editable-import override.
Every thermodynamic state, density root, reaction-reference conversion, bubble solve and derivative action is native.
The finite dilution ladder is an independent solution-verification diagnostic of EqID `reference_limit_terminal_pair`, not the reference used for the primary results.
Its successive extrapolations are compared with the accepted main revision's analytic infinite-dilution insertion; their differences are declared as reaction-log perturbations in otherwise identical native bubble re-solves.
This measures the pressure effect of reference refinement without fitting a constant to experiment.
The separate native local error estimate includes higher-precision source-reference replay and first-order pressure/site propagation.
Neither finite refinement differences nor local forward estimates are rigorous global error bounds.

## Diagnostic sources and locators

The retained [source checks](source/source-checks.json) record the recovered
Pabsch supplement, independent acid vapour-pressure correlation and association
source disagreement with exact locators.
The [access record](source/source-access.json) distinguishes recovered content
from failed requests; [Elsevier's response](source/pabsch-retrieval.xml) contains
metadata only, not the main article.

Pabsch et al. (2022), DOI 10.1016/j.ces.2022.118046:
[original supplement](source/pabsch-2022-si.docx),
[Word XML text extract](source/pabsch-2022-si.txt), Table S1 water/HAc coefficient,
Table S4 water/HAc fitted-temperature row, reference 9, and
[Figure S1](source/pabsch-figure-S1.png).
The supplement SHA-256 and publisher URL are in `source-checks.json`.
The supplement was recovered after the initial comparison; its binary fit source
and printed temperature range supersede the earlier access limitation in E2.
Its reference is Ito and Yoshida (1963), DOI 10.1021/je60018a012.
The main article's complete parameter and association declaration, actual fitted
rows, and original Ito–Yoshida paper remain unavailable.
The atmospheric comparison in Figure S1 does not establish every fitted pressure.

NIST Chemistry WebBook, acetic-acid Antoine table, McDonald, Shrader and Stull
(1959), DOI 10.1021/je60004a009: the checked row uses base-10 logarithms,
pressure in bar and temperature in kelvin, with A=4.68206, B=1642.54 K,
C=−39.764 K and range 290.26–391.01 K.
The retained evaluation at 323.15 K is in `source-checks.json`; no correlation
uncertainty is supplied.
Bernatová ThermoML group 13's pure gas-water endpoints supply independent
central endpoint measurements and compiler-expanded uncertainties.

Gross–Sadowski (2002), Table 1: acid pure-pressure/liquid-density fit over
302–592 K, mean absolute relative pressure deviation 2.12%; the dispersion
energy divided by Boltzmann's constant is 211.59 K.
Eq. (3) uses a cube of the segment-size ratio in the cross-association volume.
Ascani (2023), Table 1, instead prints acid dispersion energy 311.59 K and
water association volume 0.04509; Appendix A, Eq. (A7), prints a square of the
size ratio. These are bounded transcription/rounding probes, without adopting
new parameters or claiming which convention Pabsch's program used.

Janeček and Paricaud (2013), DOI 10.1021/jp4012125,
[primary abstract](https://pubmed.ncbi.nlm.nih.gov/23806104/): doubly bonded
cyclic dimers compete with cross association in acid mixtures with water and
alcohols, with improvements particularly in dew curves relative to classical
2B. Only the abstract was inspected; full equations and parameters were not
acquired or tested. Ascani (2023), §3.5, discusses acid dimerization and cross
association as possible error sources in related mixtures. Neither establishes
the microscopic cause or a correction for the present isotherm.
