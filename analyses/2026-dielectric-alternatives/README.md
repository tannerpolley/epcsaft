# Admitted permittivity alternatives on common NaCl(aq) data (#27)

**Question.** The Engine already implements two density-dependent
relative-permittivity forms: Schreckenberg density–temperature (full EOS §6.2)
and Rueben dipolar perturbation (§6.4). Once admitted with complete records,
do they describe NaCl(aq) MIAC, dielectric decrement, ion solvation and density
better than Figiel's ion-fraction suppression, at 298.15 K and 1 bar?

**Answer (bounded).** No, not as validated replacements.
- **Fixed-parameter swap.** Changing ε_r alone on the Figiel parameters roughly
  doubles the MIAC error at m ≤ 2 mol/kg and gives about 7× the error above
  2 mol/kg.
- **Shared-k refit.** With one water–ion k_ij refitted per form, the
  density–temperature form extrapolates MIAC above 2 mol/kg better than an
  equally refitted Figiel model. It needs k = −0.62 to do so, misses most of
  the measured dielectric decrement, and shifts solvation and density.
- **Status.** Both forms stay research choices. The admission itself is
  qualified: complete records, units, the ion-free limit, and the ε_r(T, ρ)
  chains.

## Admission (public input → native)

| Packet choice | Native form | Required records (unit) | Status |
|---|---|---|---|
| `none` | none | — | admitted |
| `solvent-only` | IonSuppression, α = 0 | solvent `relative_permittivity` | admitted, one equation |
| `ion-fraction-suppression` | IonSuppression, α_i = C_ε on ions | model `ion_fraction_suppression_coefficient` | admitted, one equation |
| `ion-specific-suppression` | IonSuppression, α_i per ion | ion `ion_specific_suppression_coefficient` | admitted, one equation |
| `component-permittivity-mixing` | ComponentMixing | component `relative_permittivity` | admitted |
| `solvent-density-temperature` | DensityTemperature | every neutral solvent: `dielectric_volume` (liter / mole), `dielectric_temperature` (kelvin); none on ions | **admitted in #27** |
| `dipolar-perturbation` | DipolarPerturbation | every component: `dipole_scaling` (debye ** 2), `polarizability_scaling` (angstrom ** 3), `correlation_integral_parameter` | **admitted in #27** |

Three labels reduce to the single IonSuppression equation, and no equation
was added. A choice rejects every permittivity record it does not consume;
for example, a constant solvent ε_r is rejected under the density forms.
Derivative dependencies (∂ε_r/∂T, ∂ε_r/∂ρ and ∂ε_r/∂x) are carried on the
common tape, as qualified natively in #26. Here they are checked again
through the admitted packet path.

`engine/tests/eos/test_permittivity_forms.py` runs on the installed wheel, for
both forms:
- the ion-free ε_r equals the §6.2/§6.4 relation, written independently from
  the packet records, to 1e-13 at liquid and vapor densities;
- a_DH, a_Born and a_res equal those of a constant-ε model set to the same
  ε_r, to 1e-12;
- Z − 1 = ρ ∂a/∂ρ is met to < 1e-9 (achieved ≤ 5.5e-12);
- u_res = −RT² ∂a/∂T is met to < 1e-9 relative (achieved ≤ 3.5e-13).
  Both identities use Richardson (h, h/2) same-branch differences.
- The ε_r chain contributes ~0.03 to Z and ~75 J/mol to u_res. That is more
  than 10⁵ times each criterion, so a dropped chain fails.

**Named unavailable cases.**
- The Yang ePPC-SAFT source pair needs polar, MSA and NAHS terms (deferred to
  #28), as does Yang's salt-corrected RSP-3/4 relation (Eq. 3, no native
  equation).
- The Rueben source pair (Held 2014 base) is not retained. Its water σ and its
  T-dependent water–Na⁺ k_ij differ from Figiel's, and FeOs's Born convention
  is unverified.
- The ε_r records are fixed inputs, not ParameterFamily coordinates, so they
  have no parameter derivatives.
- No temperature other than 298.15 K is assessed.

## Parameter compatibility and sources

The records are in [`dielectric-response-records/1`](../../data/packets/dielectric-response-records/1/packet.toml).
- **Water d_V and d_T.** Water d_V = 0.3777 dm³/mol and d_T = 1403 K come from
  Yang2024 Table 6 (preset from Schreckenberg). They are a pure-water ε_r
  correlation.
- **Rueben PT2 ion records.** Na⁺ (5.407 D², 0.215) and Cl⁻ (2.563 D², 0.138)
  come from Rueben2024 Table 1. Polarizability is 0 by the Section 2.3 scheme.
  These records were fitted to 298.15 K MIAC with Held 2014 ePC-SAFT. That base
  shares Figiel's ion σ and u/k and water m and u/k, but not the water σ, k_ij
  or Born treatment.
- **Water PT record.** Water (5.199 D², 0, 0.1276) comes from the FeOs source
  authors' fixture, a secondary source.

On the Figiel base, both forms are therefore **hybrids**:
[`inputs/selections/`](inputs/selections) holds the three explicit
water/Na⁺/Cl⁻ documents, written by `scripts/selections.py`. They differ only in
their permittivity family.

## Comparison (predeclared on #27 before running)

The six variants are:
- **F0, D0, P0:** the published Figiel parameters, with the ε_r form exchanged
  (fixed-parameter sensitivity).
- **F1, D1, P1:** one shared water–Na⁺ = water–Cl⁻ k_ij refitted per form.
  Refits minimize Rueben's relative MIAC loss on Hamer–Wu NaCl m ≤ 2 mol/kg
  (21 rows); the 8 rows at m > 2 validate.

Native fits use the same tied water–ion coordinate and bounds [-0.8, 0.4]. σ_k comes
from the Gauss–Newton curvature. "σ" is σ_k propagated to first order through
the compared metric itself, the validation MARD, in percentage points. Failed
calculations would be retained as rows with their cause and excluded from the
metrics. This run had none.

The [retained-result notebook](notebook.qmd#activity-and-dielectric-comparison)
computes the comparison table from the current retained CSVs, including the
corrected one-decimal F0 high-molality rounding.

**Data overlap:** Figiel fitted C_ε to these dielectric rows,
its Born diameters to the S5 solvation averages, and its k_ij to all Hamer–Wu
NaCl rows. Pure water is identical for all variants: −7.4e-5 relative to
IAPWS-95 (NIST SRD 69) at 298.15 K and 0.1 MPa.

**Verdict by property and domain** (NaCl(aq), 298.15 K, 1 bar only):
- **MIAC.**
  - The fixed-parameter swaps are not supported.
  - Under the same one-parameter refit, D1 beats F1 above 2 mol/kg by 4.4
    points, against a combined k uncertainty of 2.8 points.
  - P1's 1.5-point lead is within its combined 3.3 points, so it is not
    resolved.
  - D1 is 3× worse on the calibration rows and still above 13 %. Every m ≤ 2
    refit extrapolates worse than Figiel's all-molality k.
  - **Inconclusive; research choice.**
- **Dielectric decrement.**
  - At x_ion = 0.155 the data give 44.4. F0 gives 37.4, D0 70.8 and P0 75.2.
  - The D/P deficit exceeds the Nortemann–Buchner scatter (~3 units at the
    dilute end) by an order of magnitude.
  - Density–temperature carries ions only through solvent dilution, and the
    Rueben ion records were fitted to MIAC, not to ε_r.
  - **Not supported.**
- **Solvation.**
  - At fixed parameters all forms agree within 0.5 kJ/mol, because pure-water
    ε_r is 77.4–78.4 in all of them.
  - The refitted k moves D1/P1 by −10 kJ/mol. That is within the ±30 kJ/mol
    literature spread, so it does not discriminate.
  - **Equivalent / not discriminating.**
- **Density.**
  - Pure water: unaffected.
  - Solution density: no retained NaCl(aq) reference, so no claim. The strong
    k of D1/P1 raises the 6 mol/kg density 2.2–2.4 % above F1, which is an
    untested risk for any refit that uses it.

## Limits

- **Hybrids.** Only one k_ij is refitted; ion, Born and water parameters stay
  Figiel's, and these were fitted with ion-fraction suppression.
- **Uncertainty.** Hamer–Wu uncertainty is not retained, and the propagated σ covers
  k only. The dielectric rows are digitizations without pointwise uncertainty,
  and they may include kinetic depolarization (Rueben 2024).
- **Domain.** One salt, one solvent, 298.15 K and 1 bar. No temperature,
  mixed-solvent or multivalent claim is made.

Reproduce (≈6 s, installed wheel):
`.venv/bin/python analyses/2026-dielectric-alternatives/scripts/selections.py && .venv/bin/python analyses/2026-dielectric-alternatives/scripts/run.py`.
Rows: `results/comparison.csv` (one row per variant, property and observation)
and `results/refit.csv`, with a `status` per row and per fit. `results/summary.json` holds the input hashes, the
extension SHA-256, the aggregate metrics and the limits.

## Native conversion comparison

The [before/after record](results/batch4-summary.json) retains strict
prediction-comparison misses at the new optima without changing tolerances.
The historical optima have nonzero exact gradients; the native fits improve
the same-wheel objective. No calibration/holdout ranking, availability,
bound activity or physical conclusion changes. The notebook formats the current
retained operands; [exact metrics](results/batch4-metrics-before-after.csv) and
[fit uncertainties](results/batch4-fit-before-after.csv) own the current values.
Covariance is scaled by the retained variance factor; exact tied derivatives
replace forward differences in uncertainty propagation. Old wheel differences
prevent isolated attribution of every retained prediction shift.
