# Electrolyte Reference States And MIAC Validation Design

Status: approved for implementation on 2026-07-16

## Outcome

The provider will own one bounded reference-state module built on
`EPCSAFT.evaluate`. It will expose Figiel's x-process infinite-dilution
solvation, transfer between explicit solvent reference states, single-ion
activity on an electroneutral dilution path, and stoichiometric mean ionic
activity coefficients (MIACs). The validation repository will own the broad
comparison for the eight source-complete aqueous salts with the Hamer--Wu data
used by Figiel.

No Helmholtz equation is copied into the reference module or validation
repository. A successful reference-state calculation does not expand the EOS,
solvent, salt, temperature, or pressure domains carried by the selected
parameter records.

## Scientific definitions

For ion `i` in solvent state `S`, Figiel et al. (2025) define

```text
Delta_solv,x G_i^infinity = R T ln(phi_i^infinity,S)              (Eq. 17)
Delta_tr G_i^(infinity,S1->S2)
    = R T ln(phi_i^infinity,S2 / phi_i^infinity,S1)              (Eq. 19)
gamma_i^* = phi_i^S / phi_i^infinity,S                           (Eq. 20)
gamma_pm^*,x = product_i (gamma_i^*)^(nu_i / sum_j nu_j)         (Eq. 21)
gamma_pm^*,m = gamma_pm^*,x /
    [1 + M_solvent m_solute sum_i nu_i]                          (Eq. 22)
```

The single-ion values are tied to the declared electroneutral dilution path.
The MIAC is the stoichiometric combination used for comparison with
literature data. `M_solvent` is the salt-free solvent-mixture molar mass in
`kg/mol`, and `m_solute` is the molality of the undissociated formula unit.

The infinite-dilution result must be numerically checked. The implementation
will evaluate a deterministic sequence of positive formula-unit molalities,
retain electroneutrality exactly, and accept the limiting log fugacity only
when successive values meet the documented numerical tolerance. The accepted
result reports the final reference molality and maximum log-fugacity change.

## Provider interface

`EPCSAFT.reference_state(...)` is the only new entry point:

```python
reference = eos.reference_state(
    temperature=temperature,
    pressure=pressure,
    solvent_mole_fractions=solvent_mole_fractions,
    phase="liquid",
)
activity = reference.at_molality(molality)
transfer = reference.transfer_to(other_reference)
```

The call requires temperature, pressure, phase, and the salt-free solvent
composition. It admits exactly one cation and one anion formula unit and
derives the smallest integer stoichiometric coefficients from their charges.
It rejects a lone ion, multiple cations or anions, a charged solvent reference,
an unsupported phase, or missing fugacity output.

The immutable reference result owns:

- ordered component IDs and ion indices;
- stoichiometric coefficients;
- temperature, pressure, phase, and normalized salt-free solvent composition;
- infinite-dilution log fugacity coefficients;
- x-process infinite-dilution solvation free energies;
- the accepted reference molality and convergence error.

`at_molality(...)` returns the finite-state mole fractions, single-ion activity
coefficients, mole-fraction MIAC, molality MIAC, and, for the single-water
reference path, the dimensionless practical molal osmotic coefficient. The
latter uses formula-unit molality per kilogram of salt-free solvent and the
pure-water activity reference at the same temperature and pressure. A
mixed-solvent reference reports the osmotic scalar as `None` without changing
its MIAC behavior. `transfer_to(...)` returns component-ordered transfer free
energies and requires matching ions, temperature, pressure, and reference
convention.

The implementation lives in `src/epcsaft/reference.py`. `EPCSAFT` remains the
EOS owner and delegates only this reference-state orchestration. There is no
public tolerance knob, backend selector, mutable cache, or compatibility alias.

## Source parameters

One sparse Figiel catalog bundle replaces the narrow water/NaCl bundle. It
contains only the source records needed for water, methanol, and ethanol with
`Li+`, `Na+`, `K+`, `Cl-`, and `Br-`. A distinct aqueous Figiel bundle owns the
source-complete `NaI` and `KI` selections without changing established model
fingerprints. `LiI` remains unsupported because the Figiel Table 4 `Li+/I-`
cell is blank and no zero is inferred. A caller selects one exact solvent/salt
subsystem before constructing `EPCSAFT`. The eight source-complete water salts
are the admitted public reference-property scope. Methanol and ethanol retain
the existing transfer and mixed-solvent behavior; they do not create a broad
nonaqueous validation claim.

## Provider tests

Unit tests own reference-limit convergence, electroneutral composition,
Figiel equation identities, transfer identity/antisymmetry, and rejection
behavior. The provider retains exactly one physical literature anchor:

```text
system: water + Na+ + Cl-
T: 298.15 K
P: 1 bar
formula-unit molality: 1 mol/kg
Hamer--Wu Table 16 gamma_pm^m: 0.657
acceptance: absolute error <= 0.08
```

The same Hamer--Wu table reports `Phi = 0.936` at this state. It is retained as
a comparison diagnostic for the public transformation, not as a fit constant,
model-acceptance threshold, or member of Haslam's undisclosed selected rows.

This is a physical tolerance test, not a frozen implementation golden.

## Validation campaign

`ePC-SAFT-validation` will run only against an installed provider wheel. It
will retain the Hamer--Wu molality and molality-scale MIAC columns for LiCl,
NaCl, KCl, LiBr, NaBr, KBr, NaI, and KI, with DOI, table locators,
transformation record, and hashes. These are the critically evaluated source
tables cited by Figiel Figure 5; they are not digitized points from that figure. The retained
range is `0.001 <= m <= 6 mol/kg`, exactly the overlap covered by the archived
Figiel Figure 5 model curves. Higher-molality Hamer--Wu rows are outside this
campaign rather than silently extrapolated.

The campaign reports every prediction and the per-salt RMSE, MAE, bias, and
maximum absolute error, plus pooled statistics. Acceptance requires:

- no rejected, nonfinite, or nonpositive prediction;
- each salt RMSE `<= 0.35`;
- each salt maximum absolute error `<= 1.25`;
- the first reported low-molality prediction `< 0.98`; and
- pooled RMSE `<= 0.17`.

The retained archive produced pooled RMSE `0.137328` and per-salt RMSE values
from `0.028806` to `0.282501` against these tables. A separately labeled
archive-parity report may compare clean and archived model values, but it is
migration evidence rather than independent physical evidence and cannot
replace the Hamer--Wu acceptance.

The campaign retains one CSV result, one JSON receipt, and one source-data
versus model plot. It does not retain the paper PDF, the old analysis
framework, or an executable EOS copy.
