# Figiel–Yu–Held 2025 catalog review

## Sources

- Figiel, P.; Yu, G.; Held, C. “Predicting Thermodynamic Properties of
  Ions in Single Solvents and in Mixed Solvents Using a Modified Born Term
  within the ePC-SAFT Framework.” *Industrial & Engineering Chemistry
  Research* **2025**, 64, 9406–9418.
- DOI: `10.1021/acs.iecr.5c00475`
- Main article PDF SHA-256:
  `a3c940895c530f72f47f22f8c0f4796b5ad37918a4edd274f85497c6dd81ad1f`
- Supporting-information PDF SHA-256:
  `005b38ed566ec3c09b87e1ca3a9dd6eeafc9ba75e1a30b9322291d770bb93895`
- Gross, J.; Sadowski, G. “Application of the Perturbed-Chain SAFT
  Equation of State to Associating Systems.” *Industrial & Engineering
  Chemistry Research* **2002**, 41, 5510–5515.
- DOI: `10.1021/ie010954d`
- Gross–Sadowski PDF SHA-256:
  `e4991a757cad356c2d0f6ec20972adb74458218b2583cd92b19b825ccf1d89bb`
- Zotero collection: `ePC-SAFT Progression`

## Admitted records

`figiel-yu-held-2025/1` contains all 13 species for which the article reports
model parameters:

- solvents: water, methanol, and ethanol;
- cations: H+, Li+, Na+, K+, VO2+, and V3+; and
- anions: Cl-, Br-, I-, and SO4(2-).

The installed identity uses the explicit component IDs `water`, `methanol`,
`ethanol`, `hydrogen-cation`, `lithium-cation`, `sodium-cation`,
`potassium-cation`, `vanadyl-cation`, `vanadium-iii-cation`,
`chloride-anion`, `bromide-anion`, `iodide-anion`, and `sulfate-anion`.

The runtime inventory is:

- 75 single-component records: five records for each solvent and six for
  each ion;
- 41 `k_ij` records: all 15 nonblank Table 4 values and all 26 nonblank
  Table 5 values;
- six explicit 2B association sites;
- 18 association records: six pure records and 12 unlike-solvent records;
- six correlations for solvent diameter and relative permittivity; and
- three model records for the ionic-region relative permittivity, the
  ion-fraction suppression coefficient, and the Eq. 11 formulation.

Every printed zero in Tables 4 and 5 is retained as an explicit record.
Blank cells are absent. The resolver therefore rejects a selection whenever
the source does not supply every required solvent–ion, opposite-charge ion,
or solvent–solvent dispersion interaction. It does not replace a blank with
zero or infer an interaction from another salt.

## Association transformation

Figiel Table 2 reports each solvent as a 2B associating component and gives
the pure association energy and volume. The EOS association owner admits only
explicit site-pair records; it has no runtime combining rule. The three
unlike-solvent pairs are therefore evaluated once using Gross–Sadowski 2002
Eqs. 2–3 with no binary correction parameter:

| Pair | Cross energy / K | Cross volume |
| --- | ---: | ---: |
| water–methanol | 2662.6 | 0.03952009081880056 |
| water–ethanol | 2539.55 | 0.03797842740938607 |
| methanol–ethanol | 2776.45 | 0.03374752517428578 |

Each value is stored for both donor–acceptor directions. Water’s Table 2
diameter correlation is evaluated at 298.15 K for Eq. 3; methanol and ethanol
use their printed constant diameters. These records make no otherwise-blank
Figiel dispersion interaction available. In particular, methanol–ethanol
alone still fails closed because Table 5 does not report its `k_ij`.

## Supporting-information disposition

The supporting information was inspected for additional model parameters and
for the validation context behind the main article. Its model-input statements
are consistent with the main article. Tables S4 and S5 contain calculated or
experimental validation outputs, including vanadium-electrolyte density and
osmotic data and ion-solvation free energies. They are not EOS input
parameters and are not copied into the runtime bundle. Those tables belong in
installed-artifact validation evidence when a separate Validation campaign is
assigned.

## License, use basis, and tolerance

The inspected article and supporting information are publisher-copyrighted.
The package retains factual numeric parameter values, bibliographic
attribution, hashes, and short source locators; it does not redistribute the
papers, tables, figures, or prose.

Printed Table 1–5 values and explicit zeroes have zero numeric transcription
tolerance. Solvent molar masses are separate NIST identity records and are
converted by the unit system without rounding. The Gross–Sadowski cross
energies use exact decimal arithmetic on the printed pure energies. Cross
volumes are retained at the evaluated double precision shown above; tests use
exact stored decimals for the catalog contract and floating-point tolerances
only for EOS and derivative comparisons.

## Verification contract

Package tests independently fix all 13 component IDs, all 10 ion parameter
tuples, all three solvent parameter tuples, all 41 binary values, all 18
association records, source hashes, source DOIs, applicability bounds, the
bundle fingerprint, and representative selection fingerprints. Public EOS
and reference-state tests cover NaCl, vanadyl sulfate, and V(III) chloride,
including the 1:3 V(III)/chloride formula unit. Native SDK tests compare the
two multivalent phase blocks with the public Helmholtz owner and check
directional gradients and Hessians. A component-order test verifies
permutation invariance, and an unreported H+/I- pair verifies fail-closed
selection behavior.
