# Held–Cameretti–Sadowski 2008 catalog review

## Source

- Held, C.; Cameretti, L. F.; Sadowski, G. “Modeling Aqueous Electrolyte
  Solutions. Part 1. Fully Dissociated Electrolytes.” *Fluid Phase
  Equilibria* **2008**, 270, 87–96.
- DOI: `10.1016/j.fluid.2008.06.010`
- Inspected article PDF SHA-256:
  `20855b3de708fd529b29c7c4d8f6d0a1270802b610408606db1488b32d9c97e8`
- Zotero collection: `ePC-SAFT Progression`, `polar`

## Admitted records

`held-cameretti-sadowski-2008/1` contains water and all 29 ion rows printed in
Table 2:

- univalent cations: H+, Li+, Na+, K+, Rb+, Cs+, NH4+, and choline+;
- univalent anions: F-, Cl-, Br-, I-, OH-, SCN-, ClO4-, ClO3-, H2PO4-,
  BrO3-, and NO3-; and
- multivalent ions: Mg2+, Ca2+, Sr2+, Ba2+, Co2+, Cu2+, Fe2+, Cr3+, SO4(2-),
  and HPO4(2-).

Every ion uses the paper’s one-segment convention and retains its printed
diameter, dispersion energy, and charge. Table 1 water records retain the
segment count, dispersion energy, two-site association topology, association
energy, association volume, and temperature-dependent diameter correlation.
The package also records water molar mass as a separate NIST identity record.

The runtime inventory is:

- 30 components and 120 single-component records;
- 435 explicit pair records: water–ion `k_ij = 0` and ion–ion `k_ij = 1`;
- two water association sites and two association records;
- two formulation records for fully dissociated Debye–Hückel and
  solvent-only relative permittivity; and
- one water diameter correlation and one constant water relative-permittivity
  record.

The concentration-independent solvent dielectric is instantiated at 78.09
from the canonical Figiel–Yu–Held 2025 water record. This is a numerical
source supplement for the Held equation, not a Born contribution. No
ion-pairing, reaction-equilibrium, `K_w`, or Born record is admitted.

## Applicability and validation disposition

Packet version 4 corrects the Debye–Hückel closest-approach length to
`sigma_i`, as stated on p. 88 after Eq. (2) in Held 2008. Versions 2 and 3
stored `0.88 sigma_i` for both this length and hard-chain packing. Packing
remains `0.88 sigma_i`; Held 2014 Eq. (21) explicitly confirms that packing
convention for the earlier Held 2008 calculations. The 29 ion diameters and
energies, water parameters, association topology, constant dielectric, and
dispersion choices retain their values. Versions 1–3 remain historical inputs.

The [divalent chloride diagnosis](../../../analyses/2026-divalent-chloride-activity/README.md)
retains the source comparison and numerical effect. Version 4 is a correction
to the printed equation, not an independently established reproduction of
Held's original calculations or a new fit. Held Table 4 itself reports
MgCl2 activity ARD 11.08% and AAD 0.19, and CaCl2 ARD 26.11% and AAD 0.80.
The original Lobo fitting rows and Goldberg–Nuttall evaluated divalent grid
were unavailable for this diagnosis. No packet version is qualified here to
0.02 decade over the divalent fit range.

The paper reports a water correlation over 0–100 °C and ion-density fitting and
mean-ionic-activity comparisons around 20–30 °C, with the reported MIAC values
at 25 °C. The package’s exact reference domain is 298.15 K and 1 bar; evaluating
the EOS elsewhere is not a claim that the paper validated that state.

Tables 3–5 are comparison and validation material rather than runtime model
inputs. Reproducing those results belongs to installed-wheel numerical
evidence, not to the EOS catalog.

## Numerical verification

Package tests fix all 29 ion identities, charges, one-segment convention,
diameters, dispersion energies including printed zeroes, all water–ion and
ion–ion pair classifications, water association records, source metadata, and
the bundle fingerprint. They construct and evaluate representative 1:1, 2:1,
1:2, and 3:1 electroneutral selections and preserve the existing water/H+/OH-
mechanical behavior. The displaced narrow
`held-2008-water-self-ionization/1` identity is absent and does not resolve.
