# Full EOS EOS Scientific Specification

Status: scientific equation authority for the accepted full-EOS scope;
native-SDK authority is determined separately by migration promotion receipt
0015

Design: `../designs/2026-07-16-full-EOS-eos-design.md`

This file owns the scientific equations and conventions consumed by the clean
EOS. The [checked EOS equation mirror](../equations.md) maps those
relations to their current C++ and Python owners and provides synchronized
Markdown and LaTeX views. Runtime source remains executable authority; the
mirror does not create a second equation implementation.

The legacy `equations.tex` is retained only as provenance at review commit
`7920e56dbb0fc1fb418bae44e67ef50f334b75cd`, blob
`2e970476258558de05391d9d77f31d1736d78de8`, and SHA-256
`a0f2f3aff12042a8a50c5254c2af23e04b60af3e886b97b2756e6ecad800bacc`.
It does not override this bounded specification or the current implementation.

## 1. Sources and evidence classes

| ID | Source | Use |
|---|---|---|
| Gross2001 | Gross and Sadowski, *Ind. Eng. Chem. Res.* 2001, 40, 1244-1260, DOI `10.1021/ie0003887` | PC-SAFT hard-chain, dispersion, pressure, chemical potential, fugacity, and residual properties |
| Gross2002 | Gross and Sadowski, *Ind. Eng. Chem. Res.* 2002, 41, 5510-5515, DOI `10.1021/ie010954d` | Association mixing for associating systems |
| Chapman1990 | Chapman et al., *Ind. Eng. Chem. Res.* 1990, 29, 1709-1721, DOI `10.1021/ie00104a021` | Wertheim association energy and site-fraction equations |
| Bulow2019 | Bülow et al., *Fluid Phase Equilibria* 2019, 492, 26-33, DOI `10.1016/j.fluid.2019.03.010` | Debye-Huckel contribution with composition-dependent dielectric treatment |
| Figiel2025 | Figiel, Yu, and Held, *Ind. Eng. Chem. Res.* 2025, 64, 9406-9418, DOI `10.1021/acs.iecr.5c00475` | Fitted Born diameter, SSM+DS Born term, solvation factor, dielectric rule, fugacity relations, and electrolyte reference-state equations |
| Figiel2025SI | Official ACS supporting information, DOI `10.1021/acs.iecr.5c00475.s001`, Figshare article `28888514`, file `54030209`, MD5 `4c3d3964d6e61b908c3ce3b17bb4e6b0` | Comparison of direct, SSM, and SSM+DS formulations and retained validation evidence |
| Figiel2026Correction | Figiel, Yu, and Held, *Ind. Eng. Chem. Res.* 2026, DOI `10.1021/acs.iecr.6c01415` | Authoritative corrected text of Figiel2025 Eq. (7): restores the missing bracket and replaces the second `epsilon_r,bulk` with `epsilon_r,ion` |
| Zuber2014 | Zuber et al., *Fluid Phase Equilibria* 2014, DOI `10.1016/j.fluid.2014.05.037` | Ion-specific dielectric suppression, Eqs. (8)--(10) and Table 2 |
| Yang2024 | Yang, Kontogeorgis, and de Hemptinne, *Fluid Phase Equilibria* 2024, 583, 114103, DOI `10.1016/j.fluid.2024.114103` | Schreckenberg/Uematsu--Frank density-and-temperature dielectric model, Eqs. (1a), (2a), and (2b), and solvent parameters |
| Rueben2024 | Rueben et al., *J. Chem. Eng. Data* 2024, DOI `10.1021/acs.jced.4c00347` | Partial-density perturbation-theory dielectric model, Eqs. (4)--(6), and PT2 ion parameters |
| Perdomo2025 | Perdomo et al., *Comput. Chem. Eng.* 2025, 194, 108977, DOI `10.1016/j.compchemeng.2024.108977` | Read-only NaCl--water vapor/liquid endpoint diagnostic; not an ePC-SAFT parameter source |
| Held2008 | Held, Cameretti, and Sadowski, *Fluid Phase Equilibria* 2008, 270, 87-96, DOI `10.1016/j.fluid.2008.06.010` | Source-complete water plus all 29 Table 2 ion parameters and fully dissociated water-electrolyte formulation; package reference domain 298.15 K and 1 bar |
| Ascani2022 | Ascani, Sadowski, and Held, *J. Chem. Eng. Data* 2022, 67, 1972-1984, DOI `10.1021/acs.jced.1c00866` | Case Study 2 component, association, `k_ij`, `l_ij`, and `k_ij^hb` records at 298.15 K |
| AscaniHeld2021 | Ascani and Held, *Z. Anorg. Allg. Chem.* 2021, 647, 1305-1314, DOI `10.1002/zaac.202100032` | Original Born equation and combined solvent-mass/ion-mole dielectric rule |
| HamerWu1972 | Hamer and Wu, *J. Phys. Chem. Ref. Data* 1972, 1, 1047-1100, DOI `10.1063/1.3253108`; official NIST PDF SHA-256 `8399084b29b4db4ce35ab47ed679446dbedb9bf1686f99f43fce7b2073ee4beb` | Practical molal osmotic-coefficient definition and aqueous NaCl comparison value |

Each equation below has one status:

- **source equation**: a direct mathematical transcription from the cited
  source;
- **exact derivation**: algebra or thermodynamic differentiation of named
  source equations;
- **implementation convention**: an explicit indexing, unit, or API choice;
- **bounded empirical relation**: a source correlation used only in its
  retained domain.

### Archive distillation anchors

The clean implementation was distilled from the read-only lab at
`tannerpolley/ePC-SAFT-lab`; it was not independently reconstructed from this
specification alone. The principal archived owners are:

| Clean family | Archived owner or evidence |
|---|---|
| Hard-chain and dispersion | `packages/epcsaft/src/epcsaft/native/eos/residual/cppad/hard_chain_dispersion.h` |
| Association | `packages/epcsaft/src/epcsaft/native/eos/residual/contributions/association.cpp` and the archived association tests |
| Debye--Huckel and dielectric model | archived residual Helmholtz/electrolyte owners plus the checked local [equation mirror](../equations.md) distilled from the legacy LaTeX provenance |
| Corrected Born SSM+DS | `packages/epcsaft/src/epcsaft/native/eos/residual/contributions/born.cpp`, the correction paper, and the archived Figiel evidence |
| Pressure and density closure | archived state/property owners and branch tests; the clean solver deliberately removes unsafe warm-start relabelling |

For the three retained methane/ethane neutral anchors, the clean and archived
kernels agree within `5e-7` relative at identical `T`, `rho_m`, `x`, and
parameters. The remaining systematic difference is the archive's
`N_A = 6.022140857e23 mol^-1` versus the exact SI value
`6.02214076e23 mol^-1` used by the clean EOS. The clean goldens retain the
exact SI value rather than preserving an obsolete constant for bitwise parity.
Association, Debye--Huckel, Born, property, and density tests additionally use
closed-form, thermodynamic-identity, finite-difference, or branch oracles that
do not call the archived runtime.

## 2. Symbols, bases, and units

The EOS uses absolute temperature `T [K]`, molar density
`rho_m [mol m^-3]`, SI number density `rho_N = N_A rho_m [m^-3]`, and the
Angstrom-basis number density `rho_hat = rho_N / 10^30 [A^-3]`. The identity

```text
rho_m d/d(rho_m) = rho_N d/d(rho_N) = rho_hat d/d(rho_hat)
```

holds because the conversions are constant.

`a^r = A^r/(nRT)` is the dimensionless molar residual Helmholtz energy. Mole
fractions are ordered by `Parameters.component_ids` and satisfy `x_i >= 0`
and `sum_i x_i = 1`. Ordinary ionic states also satisfy
`sum_i x_i z_i = 0` within the declared state tolerance.

Component inputs use:

- segment count `m_i [-]`, segment diameter `sigma_i [A]`, and dispersion
  energy `u_i/k_B [K]`;
- association site topology, `epsilon_AB/k_B [K]`, and `kappa_AB [-]`;
- integer charge `z_i` and source-formulation-specific dielectric and Born
  records; and
- pair values `k_ij` for dispersion energy and, where the selected formulation
  requires it, explicit `l_ij` for cross diameter.

The clean native kernel converts every length to one basis before combining it
with number density. It does not mix SI number density with Angstrom cubed.

## 3. Equation registry and owners

| Family | Governing relation | Status | Source locator | Owner |
|---|---|---|---|---|
| Segment diameter | `d_i(T) = sigma_i [1 - 0.12 exp(-3(u_i/k_B)/T)]` | source equation | Gross2001 Eq. (A.6) | `_native/helmholtz.hpp` |
| Figiel ionic segment diameter | `d_i = 0.88 sigma_i` | source equation | Figiel2025 Eqs. (3), (4) | `_native/helmholtz.hpp` |
| Held water-ionization formulation | water and ion hard-chain/dispersion plus association and Debye--Huckel; `a^Born = 0` | source equation and explicit absence | Held2008 Eqs. (1)--(7), Sections 2--3 | `_native/helmholtz.hpp`, `_native/electrolyte.hpp` |
| Held ion--ion dispersion | `k_(H3O+,OH-) = 1`; water--ion `k_ij = 0` | source equation | Held2008 Eqs. (5), (6) and Section 3 | installed Held bundle and `_native/helmholtz.hpp` |
| Ascani ionic segment diameter | `d_i = sigma_i` in hard-chain, Debye--Huckel, and Born terms | source equation | AscaniHeld2021 model description following Eq. (10) | `_native/helmholtz.hpp` |
| Cross diameter | `sigma_ij = (sigma_i + sigma_j)(1-l_ij)/2` | source equation | Ascani2022SI Eq. (S4) | `_native/helmholtz.hpp`, `_native/association.hpp` |
| Mixture moments | `zeta_n = (pi/6) rho_hat sum_i x_i m_i d_i^n`, `n=0..3` | source equation | Gross2001 Eq. (A.8) | `_native/helmholtz.hpp` |
| Hard chain | hard-sphere mixture Helmholtz minus chain contact terms | source equation | Gross2001 Eqs. (A.9), (A.16), (A.17) | `_native/helmholtz.hpp` |
| Dispersion | first- and second-order perturbation terms using `I_1`, `I_2`, and `C_1` | source equation | Gross2001 Eqs. (A.10)-(A.15) | `_native/helmholtz.hpp` |
| Association | mass-action site fractions and `sum_i x_i sum_A [ln X_iA - X_iA/2 + 1/2]` | source equation | Chapman1990 Eqs. (19), (21) | `_native/association.hpp` |
| Debye-Huckel | finite-size DH residual Helmholtz term | adapted source equation | Bulow2019 Eq. (2), with diameter and screening definitions below | `_native/electrolyte.hpp` |
| Figiel bulk dielectric | salt-free solvent mass mixing followed by ion suppression | bounded empirical relation | Figiel2025 Eqs. (11), (12) | `_native/electrolyte.hpp` |
| Ascani bulk dielectric | combined salt-free-solvent mass and ion-mole mixing | source equation | AscaniHeld2021 Eq. (13) | `_native/electrolyte.hpp` |
| Schreckenberg bulk dielectric | solvent-density and temperature correlation | bounded empirical relation | Yang2024 Eqs. (1a), (2a), (2b) | `_native/electrolyte.hpp` |
| Zuber bulk dielectric | pure-solvent dielectric with ion-specific suppression | bounded empirical relation | Zuber2014 Eq. (9), Table 2 | `_native/electrolyte.hpp` |
| Rueben bulk dielectric | partial-density perturbation-theory relation | source equation | Rueben2024 Eqs. (4)--(6) | `_native/electrolyte.hpp` |
| Born SSM+DS | piecewise-dielectric shell work | corrected source equation | Figiel2026Correction Eq. (7) | `_native/born.hpp` |
| Original Born | direct charging work with ion diameter `sigma_i` | source equation | AscaniHeld2021 Eq. (10) and diameter convention | `_native/born.hpp` |
| Total Helmholtz | `a^r = a^hc + a^disp + a^assoc + a^DH + a^Born` | source equation | Figiel2025 Eq. (1) | `_native/helmholtz.hpp` |
| Pressure | `Z = 1 + rho_m (partial a^r/partial rho_m)_(T,x)` and `P = rho_m R T Z` | source equation | Gross2001 Eq. (A.24) | `_native/derivatives.hpp` |
| Chemical potential | composition derivative of the common `a^r` owner | source equation | Gross2001 Eq. (A.33); Figiel2025 Eq. (13) | `_native/derivatives.hpp` |
| Fugacity | `ln phi_i = mu_i^r/(RT) - ln Z` | source equation | Gross2001 Eq. (A.32) | `_native/derivatives.hpp` |
| Infinite-dilution solvation | `Delta_solv,x G_i^infinity = RT ln(phi_i^infinity)` | source equation | Figiel2025 Eq. (17) | `reference.py` |
| Infinite-dilution transfer | `Delta_tr G_i = RT ln(phi_i^infinity,S2/phi_i^infinity,S1)` | source equation | Figiel2025 Eq. (19) | `reference.py` |
| Ionic activity and MIAC | fugacity ratio, stoichiometric mean, and mole-fraction-to-molality conversion | source equations | Figiel2025 Eqs. (20)--(22) | `reference.py` |
| Practical molal osmotic coefficient | pure-solvent-referenced solvent activity on the formula-unit molality basis | source equation | HamerWu1972 Eqs. (3.4'), (3.5) | `reference.py` |
| Residual caloric properties | temperature and density derivatives of the common `a^r` owner | source equation and exact derivation | Gross2001 Eqs. (A.46)-(A.50) plus identities below | `_native/derivatives.hpp` |
| Density closure | `rho_m R T Z(T,rho_m,x) - P_spec = 0` | exact derivation | pressure identity | `_native/density.hpp` |

The EOS does not implement the lab catalogue's handwritten nonlinear
composition, density, or temperature derivative formulas. CppAD differentiates
the defining equations in this registry. Independent numerical formulas may
exist only in tests.

## 4. Neutral PC-SAFT

The promoted Slice 1 kernel is retained without numerical change. It evaluates
the Gross2001 temperature-corrected diameters, mixture moments, hard-sphere
Helmholtz energy, hard-chain contact correction, and two dispersion
perturbation terms. The coefficient arrays `a_0..a_6` and `b_0..b_6` are the
Gross2001 universal constants already owned by
`src/epcsaft/_native/helmholtz.hpp`.

The pair rules consumed by the admitted bundles are

```text
sigma_ij = [(sigma_i + sigma_j)/2] (1 - l_ij)
(u_ij/k_B) = sqrt[(u_i/k_B)(u_j/k_B)] (1 - k_ij).
```

The `l_ij = 0` Gross rule remains an implementation identity for bundles that
do not admit that family. The Ascani2022 bundle instead requires an explicit
`l_ij` record for every unlike pair and uses its Eq. S4 `(1-l_ij)` sign
convention. A missing applicable unlike-component `k_ij` or Ascani `l_ij`
record is an input error. Same-sign ion dispersion is equation-level
`NOT_APPLICABLE`, not an implicit zero source record. The diagonal identities
`k_ii = l_ii = 0` are not pair records.
An associating neutral subset selected from the Ascani bundle retains its
explicit source `l_ij`; only bundles with no selected `l_ij` family use the
zero-matrix identity.

For charged components, Figiel2025 sets cation--cation and anion--anion
dispersion energies to zero. This gives `d_i = 0.88 sigma_i` from Eqs. (3)--
(4) and sets every same-sign `u_ij` to zero in the dispersion sums. Opposite-
sign and solvent--ion pairs retain the Gross combining equation with their
explicit source-backed `k_ij`. These ionic rules are part of the common scalar
kernel, so Helmholtz values and CppAD composition derivatives consume the same
definition.

## 5. Association

For each explicit site `A` on component `i`, define the site fraction
`X_iA` by the fixed-topology mass-action residual

```text
H_iA(X; T,rho_hat,x) =
    X_iA [1 + rho_hat sum_j x_j sum_B X_jB Delta_iA,jB] - 1 = 0.
```

The association strength is

```text
Delta_iA,jB = g_ij^hs [exp(epsilon_iA,jB/(k_B T)) - 1]
               sigma_ij^3 kappa_iA,jB.
```

These are source equations from Chapman1990 Eqs. (19) and (18), with explicit
mixture/site indexing as an implementation convention. `g_ij^hs` is the
Gross2001 hard-sphere contact function. Site-pair energy and volume come from
explicit `AssociationParameterRecord` values. This first capability admits no
cross-association combining equation. Every allowed site pair therefore has
both an energy and volume record; an omitted unlike-site pair is inactive,
not inferred. String schemes and inferred donor/acceptor pairs are forbidden.

The source-complete bounded 298.15 K Figiel2025 bundle stores both
donor--acceptor cross directions explicitly for water--methanol,
water--ethanol, and methanol--ethanol. Their energy and volume values are
evaluated once from the Gross2002 cross-association equations using the
Figiel2025 Table 2 pure-solvent records and no binary correction parameter.
The runtime still performs no implicit combining or default insertion.

The Ascani2022 water--1-butanol selection likewise stores both cross directions
explicitly. Their energy is the SI Eq. S6 arithmetic mean multiplied by
`1-k_ij^hb` with the source value `k_ij^hb=0.026`; their volume is evaluated
once from the cited Wolbach--Sandler rule. The distinct `k_ij^hb` provenance is
therefore retained without adding a runtime pair matrix or overloading `k_ij`.

The association Helmholtz contribution is

```text
a^assoc = sum_i x_i sum_A [ln(X_iA) - X_iA/2 + 1/2].
```

The value solver runs outside CppAD recording. For independent variables `q`,
CppAD supplies the partial Jacobians of `H(X,q)` and the code solves

```text
dX/dq = -(partial H/partial X)^-1 (partial H/partial q).
```

For `q,r` drawn from `T`, `rho_m`, and independent composition coordinates,
the second sensitivity consumed by `a_TT`, `a_T,rho`, and `a_rho,rho` is

```text
X_qr = -H_X^-1 [H_qr + H_qX X_r + H_Xr X_q + H_XX[X_q,X_r]].
```

CppAD supplies every displayed partial of the defining `H`; the code performs
the linear solves. The EOS rejects `X` outside `(0,1]`, a nonfinite
residual, failed convergence, or an ill-conditioned association Jacobian.

## 6. Bulk dielectric and Debye-Huckel

The relative-permittivity topology is selected by a complete, source-bearing
parameter family. It is not a caller option and there is no phase label or
runtime fallback. Mixed or incomplete topology records fail during model
resolution.

### 6.1 Figiel composition-only and Ascani mole-fraction topologies

For the bounded Figiel2025 solvent/electrolyte model, first form salt-free
solvent mass fractions from component molar masses:

```text
w_s^(salt-free) = x_s M_s / sum_(k in solvents) x_k M_k.
epsilon_r,solvent^(salt-free) = sum_(s in solvents)
                                w_s^(salt-free) epsilon_r,s(T).
```

Figiel2025 Eq. (11) then gives

```text
epsilon_r,bulk = epsilon_r,solvent^(salt-free) /
                 [1 + C_epsilon x_ion],
x_ion = sum_(i in ions) x_i.
```

`C_epsilon = 7.01` is a source-bearing model record. The relation is empirical
for water and methanol correlations, was checked against ethanol, uses data up
to about `x_ion = 0.38`, and was evaluated by Figiel2025 at 298.15 K. Other
temperatures require source-backed pure-solvent dielectric correlations and do
not broaden the empirical ion-suppression evidence.

For the bounded Ascani2022 Case Study 2 model, AscaniHeld2021 Eq. 13 instead
defines

```text
x_sol = sum_(s in solvents) x_s,
epsilon_r,bulk = epsilon_r,solvent^(salt-free) x_sol
                 + sum_(i in ions) epsilon_r,i x_i.
```

The source records set `epsilon_r,i = 8` for the admitted alkali-halide ions.
That mixture quantity is not Figiel's radial `epsilon_r,ion`: the Ascani bundle
contains neither the Figiel radial-region model record nor a dielectric
suppression coefficient. Its water and 1-butanol source correlations are
evaluated exactly at the bundle's sole admitted temperature, 298.15 K. CppAD
differentiates the displayed composition rule in the common owner.

### 6.2 Schreckenberg density topology

Yang2024 Eqs. (1a), (2a), and (2b) give the implemented
Schreckenberg/Uematsu--Frank topology. With `S` the neutral-solvent set,

```text
x_solv = sum_(s in S) x_s,
d_V,mix = sum_(s in S) x_s d_V,s / x_solv,
d_T,mix = sum_(s in S) x_s d_T,s / x_solv,
epsilon_r,bulk = 1 + rho_m x_solv d_V,mix (d_T,mix/T - 1).
```

The admitted water comparison uses `d_V = 0.3777 dm3 mol^-1` and
`d_T = 1403 K` from Yang2024. Because solvent partial density is explicit,
`epsilon_r,bulk -> 1` as `rho_m x_solv -> 0`. Every neutral solvent requires
both parameters; missing records fail closed.

### 6.3 Zuber ion-specific topology

For a single source-supported solvent, Zuber2014 Eq. (9) is

```text
epsilon_r,bulk = epsilon_r,s(T) /
                 [1 + sum_(i in ions) alpha_i,s x_i].
```

The water/Na+/Cl- comparison uses the Table 2 coefficients
`alpha_Na+,water = 3.62` and `alpha_Cl-,water = 6.80`. These parameters are
ion--solvent specific and are not replacements for Figiel's universal
`C_epsilon = 7.01`. Zuber Eq. (8) requires solvent volume fractions for a
mixed-solvent system. The current topology therefore admits exactly one
neutral solvent and rejects a mixed-solvent selection rather than inventing
volume data or a mixing rule.

### 6.4 Rueben perturbation topology

Rueben2024 Eqs. (4)--(6) are implemented in the source's Angstrom number-density
basis. Let `rho_hat_i = rho_m N_A x_i 10^-30` in particles per Angstrom cubed:

```text
y = (4 pi/9) sum_i rho_hat_i
      [mu2_i 10^-19/(k_B T) + 3 alpha_i],
I = 1 + [sum_i a2_i x_i] [exp(-y) - 1],
epsilon_r,bulk = 1 + 3 y {1 + y + [(17/16) I - 1] y^2}.
```

`mu2_i`, `alpha_i`, and `a2_i` are the source-backed dipole-squared,
polarizability, and correlation-integral records. The water/Na+/Cl- comparison
uses water `mu2=5.199 D2`, `alpha=0 A3`, `a2=0.1276`, and Rueben2024 PT2 ion
values Na+ `(5.407 D2, 0 A3, 0.215)` and Cl- `(2.563 D2, 0 A3, 0.138)`.
All three records are required for every selected component. The equation gives
`epsilon_r,bulk -> 1` as all partial densities approach zero.

### 6.5 Common electrostatic owner and bounded comparison

The selected `epsilon_r,bulk` is evaluated once by the scalar-generic native
owner and is consumed unchanged by both Debye--Huckel and the selected Born
equation. Temperature, density, composition, and parameter dependence therefore
remain inside the common CppAD tape. `Mixture.relative_permittivity_formulation` reports the
generic bundle-defined formulation and `State.bulk_relative_permittivity` reports the
evaluated state value.

The bounded comparison uses Perdomo2025 Table 3 at `T=298.15 K` and
`P=2508 Pa`, with explicit-species compositions in EOS order
`(water,Na+,Cl-)`:

```text
x_v = (1-2(1.41063e-116), 1.41063e-116, 1.41063e-116),
x_l = (0.826112, 0.086944, 0.086944).
```

At the corrected retained vapor volume `V=0.9841 m3 mol^-1`, Figiel gives
`epsilon_r=78.09`, while Schreckenberg and Rueben give approximately
`1.00142` and `1.000324`. Thus the density topologies repair the dilute-vapor
limit at the equation level. The retained SAFT-gamma Mie liquid volume
`1.787e-5 m3 mol^-1` is not an ePC-SAFT pressure root: direct ePC-SAFT pressures
are about `-2.08e8` to `-2.16e8 Pa` across the four topologies. Comparisons of
chemical potentials therefore use each topology's own pressure-certified
root, not the retained foreign-EOS liquid volume.

The diagnostic changes only the permittivity records on the installed
Figiel2025 ePC-SAFT/SSM+DS parameter base. Its immutable fingerprints are:

```text
Figiel          sha256:79faf0da75ddb6dfe44d0282f93c5c9497941bee82227d8d76dd1149618838ef
Zuber           sha256:bb00c286a3d36a14485b96ae72577716a3a8813708d7ec5907f1da7cc03bcc1f
Schreckenberg   sha256:07d3034eb210e89c6a2ecb736735c85fafeb6161c03c4febd473428c8883d47a
Rueben          sha256:9d7ae94f2cda937ea541a3d65da3d04c66d62cca0af8fff23896456f3602cbf2
```

The vapor root is `V=0.9862436645 m3 mol^-1` for every topology. The liquid
root and electrostatic results are:

| topology | liquid `rho_m` / mol m^-3 | `epsilon_v` | `epsilon_l` | `a_l^DH` | `a_l^Born` | mean-ion gap / `RT` | water gap / `RT` |
|---|---:|---:|---:|---:|---:|---:|---:|
| Figiel composition-only | 59380.0341 | 78.0900 | 35.1922 | -0.350489 | -24.5552 | -269.7235 | 0.0170 |
| Zuber ion-specific | 59364.9061 | 78.0900 | 40.9716 | -0.292775 | -24.6290 | -269.5921 | 0.0089 |
| Schreckenberg density | 59466.9217 | 1.001419 | 69.7591 | -0.154498 | -24.8146 | -170.8193 | -0.0177 |
| Rueben perturbation | 59456.7965 | 1.000323 | 74.7311 | -0.142058 | -24.8321 | -170.6547 | -0.0231 |

Here the reported phase-equilibrium diagnostic is evaluated from EOS
fugacity coefficients,

```text
Delta_mu_ion/(RT) = (1/2) sum_(i=Na+,Cl-)
                     ln[(x_i^v phi_i^v)/(x_i^l phi_i^l)],
Delta_mu_water/(RT) = ln[(x_water^v phi_water^v)/
                          (x_water^l phi_water^l)].
```

All displayed roots have positive `dP/d(rho_m)` and pressure closure within
`3e-6 Pa`. Zuber's explicit ion coefficients barely change the approximately
`270 RT` mismatch. Schreckenberg and Rueben improve it by about `99 RT`, but a
remaining `171 RT` mean-ion gap is still far from equilibrium. The water gap is
already within `0.024 RT`, which localizes this endpoint disagreement to ionic
solvation rather than solvent phase equality.

As a liquid-property cross-check at `298.15 K`, `1 bar`, and `1 mol kg^-1`
aqueous NaCl, HamerWu1972 Table 16 reports `gamma_+/-=0.657`. With all other
Figiel parameters unchanged, the EOS returns `0.70721` for Figiel,
`0.66422` for Zuber, and `0.55220` for Rueben. Their corresponding practical
osmotic diagnostics are `0.97101`, `0.94640`, and `0.87190`; no acceptance claim
is assigned to those three osmotic numbers here. The Schreckenberg reference
sequence fails closed because its requested liquid branch becomes a lone
stable root during the dilution schedule, and the EOS correctly refuses
to relabel a single root as a requested phase. That is an unresolved
reference-path/topology interaction, not a fitted-property result.

Since the retained ePC-SAFT ion, Born, and interaction parameters were not
fitted with the density topologies, this is equation and topology evidence,
not VLE validation or authority to change a catalog default. The next
scientific step is a Regression-owned sensitivity/refit study with liquid
dielectric-decrement, density, MIAC, and osmotic anchors before cross-package
equilibrium validation. The density-based topology is worth pursuing; changing
only the Zuber suppression coefficients is not supported as the next priority
by this comparison.

The defining Debye screening quantity is

```text
kappa = sqrt[rho_N e^2/(epsilon_0 epsilon_r,bulk k_B T)
             sum_i x_i z_i^2].
```

Only charged components enter the following sum. Their finite-size DH diameter
is the same model-bound diameter used by the PC-SAFT mixture and packing owner:

```text
Figiel: d_i^DH(T) = 0.88 sigma_i,
Ascani: d_i^DH(T) = sigma_i.
```

The Figiel relation follows its Eqs. (3)--(4) and admitted `u_ii = 0`
convention. AscaniHeld2021 explicitly uses `sigma_i` in the hard-chain,
Debye--Huckel, and Born terms. The DH owner consumes the diameter already
resolved by the common segment-diameter owner rather than applying a second
scale. Here `sigma_i` and `d_i^DH` use the same length basis and `d_i^DH` is
converted to metres in the SI DH equation. For each charged component,

```text
chi_i = 3/(kappa d_i)^3
        [3/2 + ln(1+kappa d_i) - 2(1+kappa d_i)
         + (1+kappa d_i)^2/2].
```

The residual contribution is

```text
a^DH = -kappa e^2/(12 pi epsilon_0 epsilon_r,bulk k_B T)
       sum_(i in ions) x_i z_i^2 chi_i.
```

When `sum_i x_i z_i^2 = 0`, the EOS defines `kappa = 0` and `a^DH = 0`
without evaluating `chi_i`; this is the exact neutral limit. The fitted Born
diameter is used only in the Born term. CppAD differentiates the defining
relations above on the positive-ionic-strength branch; the lab
`dkappa/dx_i` catalogue entry is not a production owner.

## 7. Born SSM+DS reconciliation

Figiel2025 defines a fitted `d_i^Born` and the shell thickness

```text
f_mix = sum_k x_k f_k / sum_k x_k,
Delta d_i = [(f_mix - 1)/|z_i|] d_i^Born,
D_i = d_i^Born + Delta d_i.
```

For the admitted systems the sum contains only solvents and ions and therefore
`sum_k x_k = 1`. Ion values `f_i = 1` are source records, not defaults.

### Corrected Eq. (7)

The literal rendered Eq. (7) on article p. 9408 places
`(1 - 1/epsilon_r,bulk)` in both radial terms. That form simplifies to the
direct Born expression at `d_i^Born`; the shell thickness cancels. It therefore
contradicts all of the following evidence in the same source:

- Figure 1 shows `epsilon_r,ion` between the ion and shell and
  `epsilon_r,bulk` outside the shell;
- the preceding prose assigns `epsilon_r,ion = 8` from cited salt data;
- the following prose says SSM+DS depends separately on `d_i^Born` and `D_i`;
  and
- Figiel2025SI compares SSM+DS with direct and SSM models as distinct
  calculations and selects `s-SSM+DS-f`.

Figiel2026Correction confirms that the original Eq. (7) omitted one bracket
and incorrectly named the second dielectric constant `epsilon_r,bulk` instead
of `epsilon_r,ion`. The clean EOS implements the corrected primary-source
equation:

```text
a^(Born,SSM+DS) = -e^2/(4 pi epsilon_0 k_B T) sum_(i in ions) x_i z_i^2
  { (1 - 1/epsilon_r,bulk)/D_i
    + (1 - 1/epsilon_r,ion) [1/d_i^Born - 1/D_i] }.
```

Every `d_i^Born` and `D_i` in this equation is converted from angstrom to metre
before use with the SI factor `e^2/(4 pi epsilon_0 k_B T)`. This equation is a
**corrected source equation**. It matches the radial two-region model, retains
both limiting diameters, and matches the archived candidate kernel.
`epsilon_r,ion = 8` enters through a Figiel/Buelow source-bearing model record.
The native source must not embed 8 or 7.01.

The distinct Ascani bundle uses the original Born expression

```text
a^(Born,original) = -e^2/(4 pi epsilon_0 k_B T)
  sum_(i in ions) x_i z_i^2 (1 - 1/epsilon_r,bulk)/sigma_i.
```

The ion diameter is its source `sigma_i`, converted from angstrom to metre.
SSM and DS are disabled for this bundle: there is no fitted Born diameter,
solvation factor, shell, or radial dielectric input. The internal equation
selection follows the mutually exclusive resolved source schema and immutable
bundle fingerprint. It is not a caller-selectable Born mode or a literature
`solvation_factor`. `State.born` carries whichever Born equation the
immutable bundle selects; no second result owner is created.

## 8. Common derivatives and state properties

CppAD records the common `a^r(T,rho_m,x)` equations for fixed component and
site topology. Association uses the implicit sensitivities in Section 5.

Define

```text
a_T   = (partial a^r/partial T)_(rho_m,x)
a_TT  = (partial^2 a^r/partial T^2)_(rho_m,x)
a_rho = (partial a^r/partial rho_m)_(T,x)
Z     = 1 + rho_m a_rho
P     = rho_m R T Z.
```

The unconstrained composition partials used below follow the Gross2001
constant-`T`, constant-density mixture convention. For component `i`,

```text
mu_i^r/(RT) = a^r + (Z-1) + a_xi - sum_j x_j a_xj
ln(phi_i)   = mu_i^r/(RT) - ln Z.
```

The EOS reports the total chemical potential and fugacity coefficient.
It does not publish the lab's nonunique contribution allocation of `ln Z`.

The exact residual caloric identities at the state density are

```text
u^r/(RT) = -T a_T
h^r/(RT) = -T a_T + (Z-1)
c_v^r/R  = -2 T a_T - T^2 a_TT.
```

With `Z_T = (partial Z/partial T)_(rho_m,x)` and
`Z_rho = (partial Z/partial rho_m)_(T,x)`,

```text
c_p^r/R = c_v^r/R
          + [Z + T Z_T]^2/[Z + rho_m Z_rho] - 1.
```

The `-1` removes the ideal-gas `c_p-c_v = R` contribution. The EOS also
reports the same-density residual entropy

```text
s_V^r/R = -a^r - T a_T.
```

If it reports the departure entropy relative to an ideal gas at the same
`T,P`, the result is named separately and equals `s_V^r/R + ln Z`. No result
named only `residual_entropy` may switch reference convention by state route.

`c_p^r` requires finite, sufficiently conditioned `Z + rho_m Z_rho`. A
specified-pressure state additionally requires a positive value because that
route admits only mechanically stable roots. Explicit-density evaluation does
not reject a finite state solely because the slope is negative; it reports the
mathematical EOS result and does not certify that state as stable.

## 9. Native pure-phase Helmholtz block

The version-1 native consumer contract evaluates one fixed-temperature pure,
nonassociating phase at positive amount `n [mol]` and volume `V [m^3]`. Define

```text
rho_m   = n/V
n_ref   = 1 mol
rho_ref = 1 mol m^-3
n_hat   = n/n_ref

Phi = A/(R T n_ref)
    = n_hat * [ln(rho_m/rho_ref) - 1 + a^r(T,rho_m)].
```

The logarithmic term supplies the complete ideal density dependence needed by
nonreactive phase equilibrium. Its reference contributes only a conserved
species-linear constant; it is not a reactive standard-state or ideal-gas
caloric model. CppAD records this single `Phi(n,V)` owner and returns the value,
gradient, Hessian, and symmetric third tensor in row-major `(n_mol,V_m3)`
order. Pressure and chemical potential diagnostics are derived from that same
gradient:

```text
P [Pa]  = -R T n_ref Phi_V
mu/(RT) = Phi_n.
```

The phase block rejects mixture, association, charge, invalid source domain,
packing failure, and nonfinite states. It does not select liquid or vapor,
solve density, or certify equilibrium. One model-bound capsule is evaluated
independently at every candidate phase state.

For the first pure-methane regression slice, the same defining function accepts
the active coordinate order
`(n_mol,V_m3,m,sigma_angstrom,epsilon_over_k_kelvin)` and returns its gradient
and Hessian. The pressure and chemical-potential parameter Jacobians are exact
rows of that Hessian:

```text
dP/dq_j       = -R T Phi_(V,q_j)
d(mu/RT)/dq_j = Phi_(n,q_j).
```

The regression package may lift each observed state's two phase volumes and
assemble its own pressure, potential, and density residuals without copying an
EOS equation or differentiating a density solver. The EOS does not own
target data, residual scaling, parameter bounds, Ceres, identifiability, or fit
acceptance. Mixture, association, electrolyte, and binary-interaction active
parameters remain outside this version-1 parameterized block.

The first regression formulation fixes each phase amount at `1 mol`, uses
volume as the only row-local phase coordinate, and must reject coalesced or
mechanically unstable candidate phases. Its immutable receipt retains the
actual parameter tuple, units, transforms/scales, bounds, and source
fingerprint; the EOS's returned fingerprint continues to identify the
source model rather than the mutable trial point.

## 10. Electrolyte reference states

`Mixture.reference` holds a salt-free solvent mole composition fixed and
dilutes one cation/anion formula unit along an exactly electroneutral path. It
evaluates the complete EOS at formula-unit molalities `10^-6`, `10^-8`,
`10^-10`, and `10^-12 mol kg^-1`. The limiting ion log fugacities are accepted
only when the maximum change between the final two states is at most `5e-5`;
the accepted molality and observed change are returned as diagnostics.

For ion `i` in solvent state `S`, the x-process quantities are

```text
Delta_solv,x G_i^infinity = RT ln(phi_i^infinity,S)
Delta_tr G_i^(infinity,S1->S2)
    = RT [ln(phi_i^infinity,S2) - ln(phi_i^infinity,S1)]
gamma_i^* = phi_i^S / phi_i^infinity,S.
```

`ElectrolyteReference.infinite_dilution_log_fugacity_coefficient` stores the
two dimensionless `ln(phi_i^infinity,S)` values in `ion_component_ids` order.

For one formula unit with stoichiometric coefficients `nu_i`,

```text
ln(gamma_pm^*,x) = sum_i nu_i ln(gamma_i^*) / sum_i nu_i
gamma_pm^*,m = gamma_pm^*,x /
    [1 + M_solvent m_solute sum_i nu_i].
```

These are Figiel2025 Eqs. (17), (19), and (20)--(22). `M_solvent`
is calculated from the explicitly supplied salt-free solvent mole fractions
and source-backed solvent molar masses. The reference module calls
`Mixture.state`; it owns no Helmholtz or derivative equation. Single-ion
values are path-defined model quantities. The stoichiometric MIAC is the
literature comparison quantity.

For a single active water solvent, the same finite-molality pressure state
also gives

```text
ln(a_water) = ln(x_water) + ln(phi_water) - ln(phi_water^0)
Phi = -ln(a_water) / [m_solute M_water sum_i nu_i].
```

This is Hamer--Wu (1972) Eqs. (3.4') and (3.5) in the SI-equivalent basis:
`m_solute` is formula-unit molality in `mol/kg` of salt-free solvent and
`M_water` is the resolved molar mass in `kg/mol`. The pure-water fugacity
coefficient is evaluated at the same temperature and pressure. Mixed-solvent
references retain their MIAC result and expose no scalar osmotic coefficient.

The first admitted reference scope contains exactly one cation, one anion,
water, methanol, ethanol, or a selected mixture of those solvents at
`T = 298.15 K` within the Figiel parameter domain.
Multi-salt reference paths, non-electroneutral dilution, automatic species
generation, and alternative concentration conventions remain closed.

### 10.1 Installed charge-neutral reference basis

For an admitted aqueous phase model with `C` components, exactly one water
reference carrier, at least one cation, and at least one anion, the native SDK
exposes `C-1` rows spanning the null space of the ordered charge vector. Every
neutral component contributes its unit row. For the charged components, the
first cation is paired with every anion and every remaining cation is paired
with the first anion. Each pair uses the smallest positive integer
stoichiometry obtained from the greatest common divisor of the two charge
magnitudes. If `B` is this row-major basis and `c` its returned contractions,

```text
B q = 0
c_r = sum_i B_ri ln(phi_i^reference).
```

Each neutral contraction uses that component's infinite-dilution fugacity
coefficient in the reported pure-water composition at the same `T,P`; every
ionic formula-unit contraction uses the terminal values from its canonical
four-molality reference sequence. Any caller-owned charge-neutral row `nu` is
represented as `nu = alpha B`, so its gauge-invariant EOS reference is
`alpha c`. No individual ionic reference is exposed. In particular, the
original ordered `(water, H3O+, OH-)` case retains

```text
(-2, 1, 1) = -2 B_water + B_salt
Delta C_fugacity = c_salt - 2 c_water.
```

The result binds the ordered components and charges already advertised by the
model table, parameter fingerprint, exact source-domain `T,P`, pure-water
reference composition, water molar mass, the worst formula-unit convergence
diagnostic, and the Helmholtz basis `A/(RT n_ref)` with `n_ref = 1 mol` and
`rho_ref = 1 mol m^-3`. Existing exact catalog admissions remain intact. A full
phase-admitted `user-provided` model may expose the same reference topology,
but unknown source applicability fails value evaluation and its caller-owned
provenance is not promoted by transport. This is not a simultaneous multi-salt
finite-composition reference path. It returns values only: temperature,
pressure, composition, and parameter derivatives are explicitly unavailable.
Reaction constants, source standard-state conversions, gauges, constraints,
and equilibrium remain consumer-owned.

The Held source-complete aqueous-electrolyte bundle is admitted through this
same generic basis at exactly `298.15 K` and `100000 Pa`; it does not retain a
case-specific standard-reference callback. The ordinary public pure-water state
supplies any liquid molar-volume initialization through the same fixed-pressure
density owner. EOS does not supply `K_w` or transform an external activity
convention into this mechanical basis.

## 11. Density closure

Specified-pressure evaluation solves

```text
R(rho_m) = rho_m R T Z(T,rho_m,x) - P_spec = 0.
```

Candidate roots must have a finite residual and a finite slope

```text
dR/d(rho_m) = R T [Z + rho_m Z_rho].
```

Mechanical stability requires a positive, numerically nonsingular slope. A
deterministic logarithmic scan brackets sign-changing candidates in the
positive-density interval below the packing limit. `phase=None` accepts one
stable root. In a multi-root state, vapor and liquid requests select the
lowest- and highest-density stable roots. An unlabelled lone root is reported
as `single`; a phase-specific request rejects after that branch topology is
lost rather than relabelling the remaining root. Failed closure reports through
`EosError`; this local root calculation is not a global phase-stability or
critical-point certificate.

## 12. Domain and rejection rules

Finite temperature bounds attached to the selected source records are
intersected by the resolver and enforced for every topology. A record whose
source domain is explicitly unknown contributes no invented numerical bound.

The first electrolyte model is bounded to the Figiel source records for water,
methanol, ethanol, `Li+`, `Na+`, `K+`, `Cl-`, and `Br-` at `T = 298.15 K` and
`0 <= x_ion <= 0.38`. A distinct aqueous bundle admits source-complete NaI and
KI at the same conditions without changing established model fingerprints.
LiI remains excluded because Figiel Table 4 publishes no `Li+/I-` interaction.
Each admitted MIAC calculation selects one solvent
composition, one cation, and one anion. A mixed-solvent selection additionally
requires every solvent pair and cross-association record. The ordinary-state
composition is electroneutral. The
retained Figiel records report `p = 1 bar`; evaluation at another pressure is
a mathematical EOS result, not source-backed physical validation at that
pressure. The dielectric suppression relation is not evidence for another
solvent, temperature, component topology, or larger ion fraction. Solver
convergence does not expand that evidence.

The Schreckenberg, Zuber, and Rueben parameter families are source-named
implementation capabilities, not runtime defaults. A custom
`Parameters.from_dictionary(...)` selection must provide one complete
topology. Zuber currently requires one neutral solvent because the EOS has
no source-backed mixed-solvent volume-fraction payload. Rueben requires its
three parameters for every selected component. No new topology broadens the
temperature, pressure, composition, or fitted-parameter evidence of an
existing bundle.

The Ascani Case Study 2 bundle is a separate exact selection of water,
1-butanol, `Na+`, `K+`, and `Cl-` at `T = 298.15 K` and `p = 1 bar`. It requires
all ten unordered `k_ij` and `l_ij` pair classifications explicitly. The
same-sign `Na+/K+` dispersion pair is source-marked not applicable; its
equation-level zero cannot be used as source evidence for another pair. The
bundle also carries the source-derived water/1-butanol cross-association
energy and volume while preserving the distinct `k_ij_hb` provenance. Its
combined dielectric and original-Born equations do not alter the Figiel
SSM+DS bundles or their fingerprints.

The Held2008 bundle is a source-complete selection of water plus all 29 Table 2
ions at `T = 298.15 K` and `p = 1 bar`. It requires explicit water--ion
`k_ij = 0` and ion--ion `k_ij = 1`, so no ion--ion dispersion is silently
introduced. It uses the source's fully dissociated Debye--Huckel formulation
with concentration-independent water permittivity and has no Born term,
ion-pairing, or reaction-equilibrium behavior. Evaluation or the
standard-reference callback outside the recorded source state is not
source-backed, and the callback fails closed outside that state.

Topology is determined by the selected components, not by unrelated records
remaining in their source bundle. A selection containing no nonzero charge is
neutral or associating according to its explicit sites; ionic contributions
activate only when the selection contains both positive and negative charges.

Construction or evaluation fails on:

- missing parameter, model, pair, site, source, unit, or domain records;
- non-electroneutral ordinary ionic composition;
- nonpositive dielectric values or Born/shell diameters;
- association nonconvergence or an ill-conditioned association Jacobian;
- packing-limit violation, density branch ambiguity, or mechanical
  instability of a specified-pressure root; and
- a requested derivative or property whose defining Jacobian is singular.

## 13. Deliberate exclusions

This specification does not admit Bjerrum pairing, ion-pair equilibrium,
reaction equilibrium, caller-selected dielectric or Born modes, unsupported
Zuber mixed-solvent volume mixing, synthetic Rueben parameter completion,
analytic derivative modes, ideal-gas caloric data, phase equilibrium,
regression, multi-salt reference paths, alternate single-ion conventions,
mixture/associating/electrolyte native phase blocks, or reactive standard-state
chemical potentials. ADR 0001 defines only the fixed-temperature pure
nonassociating SDK v1 prefix and its append-only pure-parameter tail. Migration
receipt `promotion-0015-provider-native-sdk-v1` alone determines their
authority; this specification does not admit a consumer or fitted model.
Adding any other exclusion requires a new source-backed capability design
rather than a flag in this model.

The Held source-complete mechanical/reference candidate does not change those
exclusions. Reaction stoichiometry, `K_w`, mixed-standard-state conversion,
conservation constraints, Ipopt, KKT and physical certificates, globality, and
predictive adjudication remain outside EOS and require their own authority.
