# Figure 1(a) NaCl + water campaign design

## Claim

Exercise the current native HELD2 workflow over a bounded 298.15 K
NaCl + water vapor-pressure sweep and compare its phase topology with the
SAFT-\(\gamma\) Mie boundary reported by Perdomo et al. The evidence is
black-box installed-artifact validation, not a same-EOS reproduction and not
a global-stability proof.

## Source contract

The official supporting workbook is Zenodo record
`10.5281/zenodo.13646853`, file
`Zenodo_file_reviewed_NoExp (1).xlsx`, SHA-256
`c659847256044fb783c3fe090454103e109cd13cb5696bb715910428f9935834`.
Only compact Figure 1(a), 298.15 K rows are retained.

Seven unique molalities span the curve from trace salt to the paper's
6.0 mol/kg modeling limit. The packet retains both published 5.6 mol/kg
pressures:

- workbook Figure 1(a): 2541 Pa;
- paper Table 3: 2508 Pa.

The 6.0 mol/kg pressure is linearly interpolated between workbook cells
`M28:O29`; the bracketing rows themselves remain recorded in the packet.

## Execution

Run the merged Equilibrium wheel with the frozen Provider NaCl-water bundle.
For every source definition, evaluate pressures at factors 0.99, 1.00, and
1.01. The central row alone carries a source-topology comparison; bracket
rows are solver stress evidence and do not pretend to be published points.

The representative Table 3 central row is run through the native C++
diagnostic first. The retained campaign then calls only the public Python
`tp_flash` API from an isolated environment containing exact wheels.

## Decision axes

The result independently records:

1. artifact and installed `RECORD` integrity;
2. public-route availability;
3. solver and numerical status;
4. material balance, normalization, electroneutrality, and finite-state
   checks;
5. cross-EOS phase-topology comparison at source boundary rows;
6. source coverage and the 5.6 mol/kg discrepancy; and
7. search/root/globality limitations.

An error is retained as a row, never dropped. Any result that lacks
`globality_certificate="not_guaranteed"` fails closed.

## Self-review

The design uses no fitted LiCl-water-butanol inputs, no Provider or
Equilibrium source imports, no duplicated EOS/controller logic, no inferred
experimental points, and no continuous-globality claim. The only derived
source scalar is the named 6.0 mol/kg interpolation.
