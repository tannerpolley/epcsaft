# Ascani, Sadowski & Held 2022, Case Study 2: water/1-butanol/NaCl/KCl LLE

Source: Ascani et al. (2022), J. Chem. Eng. Data 67, 1972–1984,
doi:10.1021/acs.jced.1c00866 (Zotero `48XP953U`), Case Study 2, Tables 4–5, SI
Tables S1–S4; the paper's Gibbs values are in
`inputs/source/ascani-2022-case-study-2-gibbs.csv`. The [notebook](notebook.qmd) holds every result, computed from
`results/summary.json`.

## Question

Does the Engine, with the paper's ePC-SAFT advanced parameters, reproduce the
two-liquid split of the paper's feed with Na⁺, K⁺ and Cl⁻ permitted in both
liquids? The paper reports the phase compositions and the Table 5 fugacities at
298.15 K and 1 bar.

## Model and conventions

- Parameters: [packet v3](../../../data/packets/ascani-2022-case-study-2-electrolyte-lle/3/parameters/parameters.json)
  (SI Tables S1–S4; Ascani & Held 2021 Eq. 10 Born term with diameter σ and
  Eq. 13 combined permittivity rule, as Engine `component-permittivity-mixing`;
  ion dispersion between cation and anion only, `exclude-same-sign-ion-pairs`;
  water/1-butanol cross association from SI Eq. S6, k_hb = 0.026; l_ij in the
  SI Eq. S4 (1 − l_ij) form, which the Engine uses).
- Water/1-butanol k_ij: the script overrides the packet's −0.102 with Nann et
  al. 2013 Table 3, k_ij = 2.94e-4 T/K − 0.102, evaluated at 298.15 K. SI Table
  S3 prints −0.102 as k_ij(298.15 K) and cites Nann 2013; −0.102 is that
  expression's intercept. This is an inferred resolution of conflicting source
  conventions (lead decision, 2026-10-09), not a confirmed author error. With
  −0.102 a liquid disappears from the one binary start and both quaternary
  starts tried. Discovery and stability testing are unavailable for this
  permittivity form, so this does not show that no gap exists or what the
  authors ran. That run is retained.
- l_ij sign: Nann 2013 Eq. 2 uses (1 + l_ij), SI Eq. S4 uses (1 − l_ij); the
  primary run follows the SI. A diagnostic variant with l_ij = +0.0044 in the
  (1 − l_ij) form is retained and not adopted.
- Species: explicit ions. The paper's formula basis has NaCl = Na⁺ and
  KCl = K⁺, renormalized without Cl⁻. The feed is the paper's formula mass
  fractions converted to explicit mole fractions (`inputs/source/`).
- Fugacity: ln(f_i/bar) = ln x_i + ln φ_i + ln(P/bar). Mean ionic:
  ln f± = (ln f₊ + ln f₋)/2.
- Solve: two declared liquids, all five species in each, from a generic start
  and from the paper's printed phases. The Engine refuses phase discovery and
  the stability test for this permittivity form ("no low-density limit"), so
  no stability certificate exists.

## Tolerances (declared before the retained run)

- Engine solve: requested tolerance met, classified residuals ≤ 1e-10
  (row tolerance), positive reduced Hessian.
- Independent checks from returned states: material balance ≤ 1e-10,
  per-phase charge |Σ z_i x_i| ≤ 1e-10, pressure ≤ 1e-8 relative, neutral and
  mean-ionic ln f equal between phases ≤ 1e-8.
- Paper model values (numerical verification): agreement means within half a
  unit of the last printed digit. Larger differences are reported as
  deviations, not tuned.
- Physical comparison with literature parameters: mass-fraction average
  absolute deviation over all components in both phases, the printed Ascani &
  Held 2021 Eq. 15 (its preceding prose swaps Eqs. 14 and 15).

## Status and limits

- The solve converges and passes every numerical check. It does not reproduce
  the paper's printed phases. The printed phases are inconsistent with the
  printed feed: each formula-component balance alone gives a different organic
  share (0.033 to 0.819). For the paper's Table 5 ion fugacities, four
  separately inferred factors (per salt and phase) suggest a Born term 14–15 %
  stronger than Eq. 10 as printed; they are not a single factor. Neither cause
  is established; see the notebook.
- The l_ij-sign variant overshoots the paper's organic split (organic
  x_water 0.326 against the paper's 0.443; primary 0.516), so it does not fit
  clearly better.
- Physical comparison with literature parameters pending: the Pedraza et al. 2004 tie lines (Fluid
  Phase Equilib. 216, 27–31, doi:10.1016/j.fluid.2003.09.003; Zotero
  metadata-only item `UFBTEEBX`) are not in Zotero or the NIST ThermoML archive.
  Figure 4 is a 3-D perspective plot and cannot be digitized into three mass
  fractions.
- Single feed, 298.15 K, 1 bar. No Engine source change.

## Run

```bash
.venv/bin/python analyses/reproduction/2022-ascani/scripts/solve_case_study_2.py
python3 analyses/reproduction/2022-ascani/scripts/check_ascani_2022_source_contract.py
```
