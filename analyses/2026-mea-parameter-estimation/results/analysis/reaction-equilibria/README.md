# Reaction-equilibrium audit

This slice resolves Issue #118 for reactions R1--R5 over 293.15--393.15 K. The
machine-readable authority is `reaction-audit.json`; run
`python analyses/2026-mea-parameter-estimation/scripts/check_reaction_audit.py` to verify the balances, correlations,
standard-state conversions, source domains, and fit decisions.

## Verified decisions

| Reaction | Selected source relation | Source-qualified target range (K) | Target extrapolation (K) | Fit-eligible coordinates |
|---|---|---:|---:|---|
| R1 | Austgen et al. (1991), Table V: `B = -13445.9 K` | 293.15--393.15 | none | none |
| R2 | Austgen et al. (1991), Table V | 293.15--393.15 | none | none |
| R3 | Austgen et al. (1991), Table V | 293.15--393.15 | none | none |
| R4 | Tong et al. (2012), Table 5: `ln K = 2.151 - 1545.3 K/T` | 293.15--323.00 | (323.00, 393.15] | none |
| R5 | Bates and Pinching (1951), abstract and Eq. (10): `ln K = -ln(10)(2677.91 K/T + 0.3869 + 0.0004277 T/K)` | 293.15--323.15 | (323.15, 393.15] | none |

R1's negative `1/T` coefficient is confirmed directly in Austgen Table V. The
positive transcription found in historical campaign tables is rejected: at
298.15 K it gives `ln K = 215.8678`, whereas the primary relation gives
`ln K = -40.2654` before conversion.

R1--R3 were published on an unsymmetric aqueous mole-fraction basis. They are
converted to the campaign molality basis with
`ln K_m = ln K_x - delta_nu_solute ln(M_water m_standard)`, using
`M_water = 0.01801528 kg mol-1` and `m_standard = 1 mol kg-1`. R4 and R5 are
already molality-based thermodynamic relations. The campaign's 100000 Pa
reference pressure is a provider convention, not a pressure stated by these
primary correlations. Aroua's R4 precursor is the formation-direction,
concentration-basis relation `log10 K_form = -0.934 + 671/T`, with
`K_form` in `dm3 mol-1`. Inverting it changes the direction and units to
`mol dm-3`; Tong then reports the rounded, natural-log relation selected here
on its all-`K_i` molality basis. Tong's source-authored basis conversion, not
an assumed pressure conversion, is authoritative for the selected R4 values.
The Aroua measurement temperatures (298--328 K) and Tong's declared converted-
correlation range (293--323 K) are retained separately; the latter governs the
selected relation's conservative extrapolation flag.

Pressure-only VLE cannot identify reaction coordinates. The admitted Matin
speciation data resolve R4/R5 species at only one temperature, while the
multi-temperature Böttinger data aggregate free and protonated MEA. Therefore
no raw correlation coefficient has independent multi-temperature evidence or
a source-backed prior/bound sufficient for the joint fit. R4 and R5 remain
fixed source baselines, with their above-support evaluations explicitly marked
as extrapolation rather than fitted away.

## Primary locators

- Austgen et al. (1991), *Ind. Eng. Chem. Res.* 30, 543--555, DOI
  `10.1021/ie00051a018`: printed pp. 545--547, Eq. (7), reactions (1a)--(5a),
  and Table V.
- Aroua et al. (1999), *J. Chem. Eng. Data* 44, 887--891, DOI
  `10.1021/je980290n`: Eqs. (1), (2), and (14); measurements at 298, 308, 318,
  and 328 K; `log10 K_form = -0.934 + 671/T`.
- Tong et al. (2012), *Int. J. Greenhouse Gas Control* 6, 37--47, DOI
  `10.1016/j.ijggc.2011.11.005`: printed pp. 43--44, Eqs. (5)--(12), and
  Table 5.
- Bates and Pinching (1951), *J. Res. NBS* 46, 349--352, DOI
  `10.6028/jres.046.039`: printed p. 349 abstract and p. 351 Eq. (10)/Table 2.

---

Evidence labels: the equations, coefficients, conventions, locators, and
temperature domains above are verified from the cited primary records. The
fit-ineligibility conclusion is an inference from those facts and the admitted
campaign observables; it must be revisited if independent multi-temperature
true-species or thermochemical evidence is admitted.
