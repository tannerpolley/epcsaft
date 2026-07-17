# Full Provider EOS Scientific Specification

Status: implementation authority for the approved full-provider program

Design: `../designs/2026-07-16-full-provider-eos-design.md`

This file owns the scientific equations and conventions consumed by the clean
provider. The lab equation catalogue and native code remain evidence. They do
not override this bounded specification.

## 1. Sources and evidence classes

| ID | Source | Use |
|---|---|---|
| Gross2001 | Gross and Sadowski, *Ind. Eng. Chem. Res.* 2001, 40, 1244-1260, DOI `10.1021/ie0003887` | PC-SAFT hard-chain, dispersion, pressure, chemical potential, fugacity, and residual properties |
| Gross2002 | Gross and Sadowski, *Ind. Eng. Chem. Res.* 2002, 41, 5510-5515, DOI `10.1021/ie010954d` | Association mixing for associating systems |
| Chapman1990 | Chapman et al., *Ind. Eng. Chem. Res.* 1990, 29, 1709-1721, DOI `10.1021/ie00104a021` | Wertheim association energy and site-fraction equations |
| Bulow2019 | Bülow et al., *Fluid Phase Equilibria* 2019, 492, 26-33, DOI `10.1016/j.fluid.2019.03.010` | Debye-Huckel contribution with composition-dependent dielectric treatment |
| Figiel2025 | Figiel, Yu, and Held, *Ind. Eng. Chem. Res.* 2025, 64, 9406-9418, DOI `10.1021/acs.iecr.5c00475` | Fitted Born diameter, SSM+DS Born term, solvation factor, dielectric rule, fugacity and reference-state relations |
| Figiel2025SI | Official ACS supporting information, DOI `10.1021/acs.iecr.5c00475.s001`, Figshare article `28888514`, file `54030209`, MD5 `4c3d3964d6e61b908c3ce3b17bb4e6b0` | Comparison of direct, SSM, and SSM+DS formulations and retained validation evidence |

Each equation below has one status:

- **source equation**: a direct mathematical transcription from the cited
  source;
- **exact derivation**: algebra or thermodynamic differentiation of named
  source equations;
- **implementation convention**: an explicit indexing, unit, or API choice;
- **bounded empirical relation**: a source correlation used only in its
  retained domain.

## 2. Symbols, bases, and units

The provider uses absolute temperature `T [K]`, molar density
`rho_m [mol m^-3]`, SI number density `rho_N = N_A rho_m [m^-3]`, and the
Angstrom-basis number density `rho_hat = rho_N / 10^30 [A^-3]`. The identity

```text
rho_m d/d(rho_m) = rho_N d/d(rho_N) = rho_hat d/d(rho_hat)
```

holds because the conversions are constant.

`a^r = A^r/(nRT)` is the dimensionless molar residual Helmholtz energy. Mole
fractions are ordered by `ParameterSet.component_ids` and satisfy `x_i >= 0`
and `sum_i x_i = 1`. Ordinary ionic states also satisfy
`sum_i x_i z_i = 0` within the declared state tolerance.

Component inputs use:

- segment count `m_i [-]`, segment diameter `sigma_i [A]`, and dispersion
  energy `u_i/k_B [K]`;
- association site topology, `epsilon_AB/k_B [K]`, and `kappa_AB [-]`;
- integer charge `z_i`, relative permittivity `epsilon_r [-]`, fitted Born
  diameter `d_i^Born [A]`, and solvation factor `f_i [-]`; and
- pair values `k_ij`, `l_ij`, and `k_ij^hb` only when their consuming source
  equation requires them.

The clean native kernel converts every length to one basis before combining it
with number density. It does not mix SI number density with Angstrom cubed.

## 3. Equation registry and owners

| Family | Governing relation | Status | Source locator | Planned owner |
|---|---|---|---|---|
| Segment diameter | `d_i(T) = sigma_i [1 - 0.12 exp(-3(u_i/k_B)/T)]` | source equation | Gross2001 Eq. (A.6) | `_native/helmholtz.hpp` |
| Mixture moments | `zeta_n = (pi/6) rho_hat sum_i x_i m_i d_i^n`, `n=0..3` | source equation | Gross2001 Eq. (A.8) | `_native/helmholtz.hpp` |
| Hard chain | hard-sphere mixture Helmholtz minus chain contact terms | source equation | Gross2001 Eqs. (A.9), (A.16), (A.17) | `_native/helmholtz.hpp` |
| Dispersion | first- and second-order perturbation terms using `I_1`, `I_2`, and `C_1` | source equation | Gross2001 Eqs. (A.10)-(A.15) | `_native/helmholtz.hpp` |
| Association | mass-action site fractions and `sum_i x_i sum_A [ln X_iA - X_iA/2 + 1/2]` | source equation | Chapman1990 Eqs. (19), (21) | `_native/association.hpp` |
| Debye-Huckel | finite-size DH residual Helmholtz term | adapted source equation | Bulow2019 Eq. (2), with diameter and screening definitions below | `_native/electrolyte.hpp` |
| Bulk dielectric | salt-free solvent mass mixing followed by ion suppression | bounded empirical relation | Figiel2025 Eqs. (11), (12) | `_native/electrolyte.hpp` |
| Born SSM+DS | piecewise-dielectric shell work | exact source reconciliation | Figiel2025 p. 9408, Eq. (7), Figure 1, surrounding prose; discrepancy recorded below | `_native/born.hpp` |
| Total Helmholtz | `a^r = a^hc + a^disp + a^assoc + a^DH + a^Born` | source equation | Figiel2025 Eq. (1) | `_native/helmholtz.hpp` |
| Pressure | `Z = 1 + rho_m (partial a^r/partial rho_m)_(T,x)` and `P = rho_m R T Z` | source equation | Gross2001 Eq. (A.24) | `_native/derivatives.hpp` |
| Chemical potential | composition derivative of the common `a^r` owner | source equation | Gross2001 Eq. (A.33); Figiel2025 Eq. (13) | `_native/derivatives.hpp` |
| Fugacity | `ln phi_i = mu_i^r/(RT) - ln Z` | source equation | Gross2001 Eq. (A.32) | `_native/derivatives.hpp` |
| Residual caloric properties | temperature and density derivatives of the common `a^r` owner | source equation and exact derivation | Gross2001 Eqs. (A.46)-(A.50) plus identities below | `_native/derivatives.hpp` |
| Density closure | `rho_m R T Z(T,rho_m,x) - P_spec = 0` | exact derivation | pressure identity | `_native/density.hpp` |

The provider does not implement the lab catalogue's handwritten nonlinear
composition, density, or temperature derivative formulas. CppAD differentiates
the defining equations in this registry. Independent numerical formulas may
exist only in tests.

## 4. Neutral PC-SAFT

The promoted Slice 1 kernel is retained without numerical change. It evaluates
the Gross2001 temperature-corrected diameters, mixture moments, hard-sphere
Helmholtz energy, hard-chain contact correction, and two dispersion
perturbation terms. The coefficient arrays `a_0..a_6` and `b_0..b_6` are the
Gross2001 universal constants already owned by
`src/epcsaft/_native/eos.hpp` at Slice 1.

The pair rules consumed by this slice are

```text
sigma_ij = (sigma_i + sigma_j)/2
(u_ij/k_B) = sqrt[(u_i/k_B)(u_j/k_B)] (1 - k_ij).
```

Both are source equations from Gross2001. A missing `k_ij` record remains an
input error under the clean parameter contract; the resolver does not assume
zero.

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
explicit `AssociationParameterRecord` values or from an admitted source
combining equation with all operands present. String schemes and inferred
donor/acceptor pairs are forbidden.

The association Helmholtz contribution is

```text
a^assoc = sum_i x_i sum_A [ln(X_iA) - X_iA/2 + 1/2].
```

The value solver runs outside CppAD recording. For independent variables `q`,
CppAD supplies the partial Jacobians of `H(X,q)` and the code solves

```text
dX/dq = -(partial H/partial X)^-1 (partial H/partial q).
```

Higher derivatives are added only when a listed provider property consumes
them. The provider rejects `X` outside `(0,1]`, a nonfinite residual, failed
convergence, or an ill-conditioned association Jacobian.

## 6. Bulk dielectric and Debye-Huckel

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

The defining Debye screening quantity is

```text
kappa = sqrt[rho_N e^2/(epsilon_0 epsilon_r,bulk k_B T)
             sum_i x_i z_i^2].
```

For each component, with the finite-size ion diameter `d_i` in meters,

```text
chi_i = 3/(kappa d_i)^3
        [3/2 + ln(1+kappa d_i) - 2(1+kappa d_i)
         + (1+kappa d_i)^2/2].
```

The residual contribution is

```text
a^DH = -kappa e^2/(12 pi epsilon_0 epsilon_r,bulk k_B T)
       sum_i x_i z_i^2 chi_i.
```

The clean Figiel formulation uses the source-defined ion diameter relation for
the DH term and the fitted Born diameter only in the Born term. The scientific
implementation plan must bind the exact Figiel diameter relation and source
domain before Task 6 code is admitted. CppAD differentiates the three defining
relations above; the lab `dkappa/dx_i` catalogue entry is not a production
owner.

## 7. Born SSM+DS reconciliation

Figiel2025 defines a fitted `d_i^Born` and the shell thickness

```text
f_mix = sum_k x_k f_k / sum_k x_k,
Delta d_i = [(f_mix - 1)/|z_i|] d_i^Born,
D_i = d_i^Born + Delta d_i.
```

For the admitted systems the sum contains only solvents and ions and therefore
`sum_k x_k = 1`. Ion values `f_i = 1` are source records, not defaults.

### Published Eq. (7) discrepancy

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

The clean provider classifies the printed second `epsilon_r,bulk` as a likely
typesetting error. It implements the exact piecewise-dielectric work implied by
Figure 1 and the prose:

```text
a^(Born,SSM+DS) = -e^2/(4 pi epsilon_0 k_B T) sum_(i in ions) x_i z_i^2
  { (1 - 1/epsilon_r,bulk)/D_i
    + (1 - 1/epsilon_r,ion) [1/d_i^Born - 1/D_i] }.
```

This equation is an **exact source reconciliation**, not a literal
transcription of printed Eq. (7). It matches the radial two-region model,
retains both limiting diameters, and matches the archived candidate kernel.
`epsilon_r,ion = 8` enters through a Figiel/Buelow source-bearing model record.
The native source must not embed 8 or 7.01.

The provider does not expose direct Born, SSM-only, or DS-only runtime modes.
The admitted ionic formulation uses SSM+DS as one fixed equation.

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

The provider reports the total chemical potential and fugacity coefficient.
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

The `-1` removes the ideal-gas `c_p-c_v = R` contribution. The provider also
reports the same-density residual entropy

```text
s_V^r/R = -a^r - T a_T.
```

If it reports the departure entropy relative to an ideal gas at the same
`T,P`, the result is named separately and equals `s_V^r/R + ln Z`. No result
named only `residual_entropy` may switch reference convention by state route.

`c_p^r` requires finite `Z + rho_m Z_rho` on a mechanically stable branch.
The provider rejects singular or ill-conditioned states instead of returning a
finite-looking property.

## 9. Density closure

Specified-pressure evaluation solves

```text
R(rho_m) = rho_m R T Z(T,rho_m,x) - P_spec = 0.
```

Candidate roots must have a finite residual and a finite slope

```text
dR/d(rho_m) = R T [Z + rho_m Z_rho].
```

Mechanical stability requires a positive slope. A deterministic scan and
bracket identifies all candidates in the admitted positive-density interval.
`phase=None` accepts one certified stable root. Vapor and liquid requests
select the lowest- and highest-density certified stable roots. The solver
rejects branch loss, coalescence, poor conditioning, and absent requested
roots.

## 10. Reference-state electrolyte quantities

Figiel2025 Eqs. (17), (19), and (20) define the admitted x-process quantities:

```text
Delta_solv G_i^infinity = R T ln(phi_i^infinity)
Delta_transfer G_i^(S1->S2,infinity)
    = R T ln[phi_i^(S2,infinity)/phi_i^(S1,infinity)]
gamma_i^* = phi_i/phi_i^infinity.
```

For an electrolyte with stoichiometric coefficients `nu_c` and `nu_a`,

```text
ln gamma_pm^* = [nu_c ln gamma_c^* + nu_a ln gamma_a^*]
                /(nu_c + nu_a).
```

The provider approaches infinite dilution along an explicit electroneutral
stoichiometric path. It does not evaluate a lone non-electroneutral bulk ion.
Each API names the x, molality, or other concentration scale; scale conversion
is not implicit.

## 11. Domain and rejection rules

The first electrolyte capability is bounded to the complete source-backed
water/sodium/chloride bundle and its retained conditions. The Figiel dielectric
suppression relation is not evidence for arbitrary solvents, temperatures, or
ion fractions above its reported range. A mathematically converged state
outside the source domain is not physically accepted.

Construction or evaluation fails on:

- missing parameter, model, pair, site, source, unit, or domain records;
- non-electroneutral ordinary ionic composition;
- nonpositive dielectric values or Born/shell diameters;
- association nonconvergence or an ill-conditioned association Jacobian;
- packing-limit violation, density branch ambiguity, or mechanical
  instability; and
- a requested derivative or property whose defining Jacobian is singular.

## 12. Deliberate exclusions

This specification does not admit Bjerrum pairing, ion-pair equilibrium,
reaction equilibrium, alternate dielectric rules, direct/SSM-only Born modes,
analytic derivative modes, ideal-gas caloric data, phase equilibrium,
regression, or a public derivative/SDK contract. Adding one requires a new
source-backed capability design rather than a flag in this model.
