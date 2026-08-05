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

The paper reports a water correlation over 0–100 °C and ion-density fitting and
mean-ionic-activity comparisons around 20–30 °C, with the reported MIAC values
at 25 °C. The package’s exact reference domain is 298.15 K and 1 bar; evaluating
the EOS elsewhere is not a claim that the paper validated that state.

Tables 3–5 are comparison and validation material rather than runtime model
inputs. Reproducing those results belongs to installed-artifact Validation
evidence, not to the EOS catalog.

## Verification contract

Package tests fix all 29 ion identities, charges, one-segment convention,
diameters, dispersion energies including printed zeroes, all water–ion and
ion–ion pair classifications, water association records, source metadata, and
the bundle fingerprint. They construct and evaluate representative 1:1, 2:1,
1:2, and 3:1 electroneutral selections and preserve the existing water/H+/OH-
mechanical behavior. The displaced narrow
`held-2008-water-self-ionization/1` identity is absent and does not resolve.
