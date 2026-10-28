# Full EOS scientific specification

This document owns symbols, bases, source applicability and EOS conventions.
It states required behavior, not installed availability. The canonical
[equation definitions](../equations.md) are edited once; repeated derivations
below explain their meaning. [Algorithms](algorithms.md) own numerical methods
and the [GitHub issue roadmap](https://github.com/tannerpolley/ePC-SAFT/milestones) owns implementation gaps.
Correct mathematics is not weakened to match an implementation defect.

## 1. Sources and evidence classes

| ID | Source | Use |
|---|---|---|
| Gross2001 | Gross and Sadowski, *Ind. Eng. Chem. Res.* 2001, 40, 1244-1260, DOI `10.1021/ie0003887` | PC-SAFT hard-chain, dispersion, pressure, chemical potential, fugacity, and residual properties |
| Gross2006 | Gross and Vrabec, *AIChE J.* 2006, 52, 1194-1204, DOI `10.1002/aic.10683` | Optional dipole-dipole term, Eqs. (7)-(16) and Table 1 constants; Table 2 PCP-SAFT parameters (ADR 0011) |
| Gross2005 | Gross, *AIChE J.* 2005, 51, 2556-2568, DOI `10.1002/aic.10502` | Optional quadrupole-quadrupole term, Eqs. (8)-(19) and Table 2 constants; Table A1 PCP-SAFT parameters (ADR 0014) |
| Gross2002 | Gross and Sadowski, *Ind. Eng. Chem. Res.* 2002, 41, 5510-5515, DOI `10.1021/ie010954d` | Association mixing for associating systems |
| Chapman1990 | Chapman et al., *Ind. Eng. Chem. Res.* 1990, 29, 1709-1721, DOI `10.1021/ie00104a021` | Wertheim association energy and site-fraction equations |
| Bulow2019 | Bülow et al., *Fluid Phase Equilibria* 2019, 492, 26-33, DOI `10.1016/j.fluid.2019.03.010` | Debye-Huckel contribution with composition-dependent dielectric treatment |
| Figiel2025 | Paul Figiel, Gangqiang Yu, and Christoph Held, *Ind. Eng. Chem. Res.* 2025, 64(18), 9406-9418, DOI `10.1021/acs.iecr.5c00475` | Fitted Born diameter, SSM+DS Born term, solvation factor, dielectric rule, fugacity relations, and electrolyte reference-state equations |
| Ji2012 | Ji, Held, and Sadowski, *Fluid Phase Equilibria* 2012, 335, 64-73, DOI `10.1016/j.fluid.2012.05.029` | Fully dissociated ionic-liquid model with dispersion among all ion pairs |
| Figiel2025SI | Official ACS supporting information, DOI `10.1021/acs.iecr.5c00475.s001`, Figshare article `28888514`, file `54030209`, MD5 `4c3d3964d6e61b908c3ce3b17bb4e6b0` | Comparison of direct, SSM, and SSM+DS formulations and retained validation evidence |
| Figiel2026Correction | Figiel, Yu, and Held, *Ind. Eng. Chem. Res.* 2026, DOI `10.1021/acs.iecr.6c01415` | Authoritative corrected text of Figiel2025 Eq. (7): restores the missing bracket and replaces the second `epsilon_r,bulk` with `epsilon_r,ion` |
| Zuber2014 | Zuber et al., *Fluid Phase Equilibria* 2014, DOI `10.1016/j.fluid.2014.05.037` | Ion-specific dielectric suppression, Eqs. (8)--(10) and Table 2 |
| Yang2024 | Yang, Kontogeorgis, and de Hemptinne, *Fluid Phase Equilibria* 2024, 583, 114103, DOI `10.1016/j.fluid.2024.114103` | Schreckenberg/Uematsu--Frank density-and-temperature dielectric model, Eqs. (1a), (2a), and (2b), and solvent parameters |
| Rueben2024 | Rueben et al., *J. Chem. Eng. Data* 2024, DOI `10.1021/acs.jced.4c00347` | Partial-density perturbation-theory dielectric model, Eqs. (4)--(6), and PT2 ion parameters |
| Perdomo2025 | Perdomo et al., *Comput. Chem. Eng.* 2025, 194, 108977, DOI `10.1016/j.compchemeng.2024.108977` | Read-only NaCl--water vapor/liquid endpoint diagnostic; not an ePC-SAFT parameter source |
| Held2008 | Held, Cameretti, and Sadowski, *Fluid Phase Equilibria* 2008, 270, 87-96, DOI `10.1016/j.fluid.2008.06.010` | Source-complete water plus all 29 Table 2 ion parameters and fully dissociated water-electrolyte formulation; package reference domain 298.15 K and 1 bar |
| Bulow2021 | Bülow, Ascani, and Held, *Fluid Phase Equilibria* 2021, 537, 112989, DOI `10.1016/j.fluid.2021.112989` | Linear ion-mole permittivity rule (§2.3, Table 5) and Born charging work (Eq. 24) for salts in an ionic liquid, with no neutral species |
| Held2014 | Held et al., *Chem. Eng. Res. Des.* 2014, 92, 2884-2897, DOI `10.1016/j.cherd.2014.05.017` | Affine `k_ij(T)` laws of Tables 1-2 (pp. 2888-2889), temperature-law falsifiers; the text after Eqs. (21)-(22), the ion diameter from zero like-ion dispersion (§4) |
| Uyan2015 | Uyan et al., *Fluid Phase Equilibria* 2015, 393, 91-100, DOI `10.1016/j.fluid.2015.02.026` | Quadratic `k_ij(T)` at 298 K, Eqs. (26)-(27), and the Wolbach--Sandler volume correction outside the square root, Eq. (9) and Table 5; Eqs. (10)-(11), the ion diameter from zero like-ion dispersion (§4) |
| Bulow2021MDEA | Bülow et al., *Ind. Eng. Chem. Res.* 2021, 60, 6327-6336, DOI `10.1021/acs.iecr.1c00176` | Cubic water--MDEA `k_ij(T)`, Eq. (26) and Table 3, temperature-law falsifier only |
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

## 3. Equation and implementation ownership

The [equation catalogue](../equations.md) gives each defining relation a stable
ID and source locator. Its YAML index is generated, not independently edited.
Current native owners and numerical falsifiers are linked by the
[GitHub issue roadmap](https://github.com/tannerpolley/ePC-SAFT/milestones). Historical `_native`
paths and pure-version interfaces do not describe current ownership.

ADR 0006 retires residual point-multipole behavior while preserving the
separate dipolar-density-perturbation relative-permittivity formulation; ADR 0011
readmits the optional dipole-dipole term of section 4.1 and ADR 0014 the optional
quadrupole-quadrupole term of section 4.2. The dipole-quadrupole cross term stays retired.
The common Helmholtz definition owns thermodynamic derivatives. Independent
finite-difference formulas belong to numerical verification, not an alternate
production equation owner.

See [Figiel input identity and calibration evidence](#figiel-input-identity-and-calibration-evidence)
for the adopted packet and its source limits.

## 4. Neutral PC-SAFT

The neutral model uses the Gross2001 temperature-corrected diameters, mixture moments, hard-sphere
Helmholtz energy, hard-chain contact correction, and two dispersion
perturbation terms. The coefficient arrays `a_0..a_6` and `b_0..b_6` are the
Gross2001 universal constants retained as `A0`, `A1`, `A2`, `B0`, `B1`, `B2` in
`engine/native/contributions/helmholtz.hpp`, transcribed from Gross2001 Table 1.
The catalogue defines their interpolation; it does not duplicate the arrays.

The pair rules consumed by the admitted bundles are

Defining relation: [EqID: `cross_segment_diameter`](../equations.md#cross_segment_diameter), [EqID: `cross_dispersion_energy`](../equations.md#cross_dispersion_energy), [EqID: `same_sign_ionic_dispersion`](../equations.md#same_sign_ionic_dispersion).

The `l_ij = 0` Gross rule remains an implementation identity for bundles that
do not admit that family. The Ascani2022 bundle instead requires an explicit
`l_ij` record for every unlike pair and uses its Eq. S4 `(1-l_ij)` sign
convention. A missing applicable unlike-component `k_ij` or Ascani `l_ij`
record is an input error. For electrolyte families declaring same-sign ion
exclusion, same-sign dispersion is equation-level `NOT_APPLICABLE`, not an
implicit zero source record. Families declaring all-ion-pair dispersion retain
the Gross combining rule for those pairs. The diagonal identities
`k_ii = l_ii = 0` are not pair records.
An associating neutral subset selected from the Ascani bundle retains its
explicit source `l_ij`; only bundles with no selected `l_ij` family use the
zero-matrix identity.

For charged components, the electrolyte model family declares its source rule.
Figiel2025 sets cation--cation and anion--anion dispersion energies to zero;
Ji2012 retains dispersion among all ionic species. An ion declaring
`fixed.charged_molecule_ion` is a charged molecule: it keeps dispersion with
itself and with every other declared charged molecule, and takes the neutral
Gross2001 diameter with its own `u_i` instead of a packing diameter.
Held2014 and Uyan2015 derive the fixed ion diameter from zero like-ion
dispersion; that a kept like-ion dispersion restores the temperature-dependent
diameter is our inference from that derivation, which the sources state only in
the zero-dispersion direction. Ji2012 is precedent for kept like-ion dispersion,
not for the ion diameter, which its text does not state.
The declaration is opt-in with no default. Opposite-sign and
solvent--ion pairs always retain the Gross combining equation with their
explicit source-backed `k_ij`. These ionic rules are part of the common scalar
kernel, so Helmholtz values and CppAD composition derivatives consume the same
definition.

### 4.1 Optional dipole-dipole term

A model that declares `polar: dipole-dipole` adds the Gross2006 third-order
Padé dipole-dipole contribution to the residual Helmholtz energy. Every component
supplies `dipole_moment` (debye) and `dipole_count` n_μ; zero marks a nonpolar
species. With the reduced moment written
`q_i = n_μ,i μ_i^2/(4π ε_0 k T m_i σ_i^3)`, ε_i appears only in the J2 energy
`ε_ij = (ε_i ε_j)^0.5`, which carries no k_ij. The chain lengths in the
J integrals are capped at 2 (source Eqs. 15-16). The term uses the PC-SAFT
packing fraction and segment diameters of the same state.

Defining relation: [EqID: `dipolar_helmholtz`](../equations.md#dipolar_helmholtz), [EqID: `dipolar_second_order`](../equations.md#dipolar_second_order), [EqID: `dipolar_third_order`](../equations.md#dipolar_third_order), [EqID: `dipolar_integrals`](../equations.md#dipolar_integrals).

When every dipolar species is absent, A2 = A3 = 0 and the Padé factor takes its
limit 1. The value and the first composition derivative then equal the analytic
limit a_dd ≈ A2 ∝ x_μ², so the absent species' dipolar chemical-potential share is
zero. Higher composition derivatives along an absent dipolar direction are
outside the claim; public actions report a support change there. The constants are fitted to two-centre
Lennard-Jones plus point-dipole simulations at 0.55 ≤ T/T_c ≤ 0.99 with
elongation L* ≤ 1. Polar parameters are valid only with their own fitted m, σ,
ε/k and k_ij (ADR 0011).

### 4.2 Optional quadrupole-quadrupole term

A model that declares `polar: quadrupole-quadrupole` adds the Gross2005 third-order
Padé quadrupole-quadrupole contribution instead. Every component supplies
`quadrupole_moment` (debye angstrom) and `quadrupole_count` n_Q; zero marks a
nonquadrupolar species. The reduced moment is
`q_i = n_Q,i Q_i^2/(4π ε_0 k T m_i σ_i^5)`, with the same debye bridge as the dipole,
and the J integrals use the uncapped `m_ij = (m_i m_j)^0.5` and
`m_ijk = (m_i m_j m_k)^(1/3)` of source Eqs. 18-19. The two-body third-order
integral is zero (Eq. 14).

Defining relation: [EqID: `quadrupolar_helmholtz`](../equations.md#quadrupolar_helmholtz), [EqID: `quadrupolar_second_order`](../equations.md#quadrupolar_second_order), [EqID: `quadrupolar_third_order`](../equations.md#quadrupolar_third_order), [EqID: `dipolar_integrals`](../equations.md#dipolar_integrals).

The absent-species limit and the support-change status are those of section 4.1.
The constants are fitted to two-centre Lennard-Jones plus point-quadrupole
simulations at elongations 0 ≤ L* ≤ 0.8 (m up to 1.8447) and 0.5 ≤ T/T_c ≤ 0.99;
larger m, such as benzene's 2.2463, is the source's own extrapolation. The source
finds its mixing rule unsuitable for real mixtures with more than one quadrupolar
component. A model declares either the dipole or the quadrupole term, never both.

## 5. Association

The version-1 parameter document declares exactly one choice for each
base, association, electrolyte, and permittivity family. Unknown or duplicate
family kinds are input errors. Native admission accepts base pc-saft and
association none or general-site. Association none requires no selected sites;
contradictory topology is rejected. General-site with no selected sites has
the exact zero-association limit. Changing a choice does not establish that
an existing fit transfers to the new formulation.

Component IDs must be nonempty strings and unique in the supplied document
before subset selection. Repeated records are ambiguous physical inputs, not
parameter overrides. Every selected component supplies fixed.charge_number,
including zero for a neutral species. Python transports that record; native
admission rejects omission by component identity. No charge is inferred from
a component name, source label, or absent field. A charged component may
supply fixed.charged_molecule_ion with magnitude 1 (§4); admission refuses it
on a neutral component, beside a packing_diameter, or with another magnitude. The version-3 Gross/DEEA
packets state previously implicit neutral charges without changing coefficients.

For each explicit site `A` on component `i`, define the site fraction
`X_iA` by the fixed-topology mass-action residual

Defining relation: [EqID: `association_mass_action`](../equations.md#association_mass_action).

The association strength is

Defining relation: [EqID: `association_strength`](../equations.md#association_strength).

These are source equations from Chapman1990 Eqs. (19) and (18), with explicit
mixture/site indexing as an implementation convention. `g_ij^hs` is the
Gross2001 hard-sphere contact function. Chapman1990 Eq. 18 carries `d^3`; the
strength here carries `sigma_ij^3`, the PC-SAFT convention written explicitly by
Nasrifar and Tafazzol 2010 (DOI `10.1021/ie901181n`, Eq. 34). A parameter set
regressed with `d_ij^3` is a different model: its volume must be converted in the
packet, and no constant conversion is exact because `(d/sigma)^3` depends on
temperature. See [association source conversions](#association-source-conversions)
for packet-specific evidence. Every admitted site pair is an explicit
topology edge whose source is either a direct coefficient record or one of two
named closed-registry rules. For unlike components `i` and `j`, both rules
resolve

```text
epsilon_iA,jB/k_B = [(epsilon_iA,iB/k_B) + (epsilon_jA,jB/k_B)] / 2
kappa_iA,jB = sqrt(kappa_iA,iB kappa_jA,jB).
```

`arithmetic-energy-geometric-volume` uses these values unchanged. `wolbach-sandler`
is Gross2002 Eqs. 2--3, which multiply the volume by the size factor
`[sqrt(sigma_i sigma_j) / ((sigma_i + sigma_j)/2)]^3`. Because Delta already
carries `sigma_ij^3 = [(sigma_i + sigma_j)/2 (1 - l_ij)]^3`, the size factor is
applied exactly by using `[sqrt(sigma_i sigma_j)(1 - l_ij)]^3` in Delta for
that edge. The state-evaluated `sigma(T)` keeps the temperature dependence and
its derivatives, which matter for the water correlation. The size factor is 1
for like pairs.

Both inputs must be donor--acceptor self-edges for the same two components and
must contain the output sites. The output must itself connect one donor and one
acceptor. Resolved outputs are derived, non-selectable coordinates; their
recorded derivative chain is `1/2` for each energy input and
`kappa_iA,jB/(2 kappa_input)` for each volume input. This lets later exact
parameter derivatives propagate through the rule without treating its output
as an independent fit coefficient.

See [association source conversions](#association-source-conversions)
for the immutable packet history. An omitted unlike-site pair remains inactive,
not inferred. String schemes, inferred donor/acceptor pairs, implicit combining,
and default insertion are forbidden.

Every selected pair of components that both carry sites must therefore declare
its cross association: at least one unlike-site edge (explicit or a named
rule), or one `topology.no_cross_association` item naming `component_id_a`,
`component_id_b` and its `provenance`. The declaration makes every unlike-site
edge of that pair inactive; it preserves a model, not evidence that the physical
interaction is absent. `build_parameters` refuses a missing declaration before
any evaluation, naming both components, and refuses a duplicate, a self or
unknown pair, a nonzero native marker and a declaration contradicted by a
resolved active edge. A declaration outside the selection is ignored. The native
entry is `no_cross_association` on the unordered component pair with a zero
dimensionless marker, so Python, direct C++ and ASL records share this owner.

The Ascani2022 water--1-butanol selection likewise stores both cross directions
explicitly. Their energy is the SI Eq. S6 arithmetic mean multiplied by
`1-k_ij^hb` with the source value `k_ij^hb=0.026`; their volume is evaluated
once from the cited Wolbach--Sandler rule. The distinct `k_ij^hb` provenance is
therefore retained without adding a runtime pair matrix or overloading `k_ij`.

**Temperature laws** (#567). A pair `k_ij` or `l_ij` record, an explicit edge's
`energy_over_k` or `volume`, and a combining-rule edge's `volume_correction`
carry exactly one of `value` or `temperature_law`:

```json
{"identity": "uyan2015/water-co2/k_ij", "family": "k_ij",
 "temperature_law": {"form": "polynomial-plus-reciprocal-temperature",
   "reference_temperature": {"magnitude": 298.0, "unit": "kelvin"},
   "terms": [{"power": 0, "coefficient": {"identity": "uyan2015/water-co2/a0", "value": {"magnitude": -0.022, "unit": "dimensionless"}}},
             {"power": 1, "coefficient": {"identity": "uyan2015/water-co2/a1", "value": {"magnitude": 0.00042, "unit": "1 / kelvin"}}}],
   "reciprocal_temperature_slope": {"identity": "...", "value": {"magnitude": 0.0, "unit": "kelvin"}}}}
```

Defining relation: [EqID: `pair_parameter_temperature_law`](../equations.md#pair_parameter_temperature_law).
Powers are unique integers 0-3 and the power-0 term is required; omitted terms
are absent. `T_r > 0` is fixed, not fitted. Coefficient `a_p` has unit
`kelvin ** (e - p)` and `s` unit `kelvin ** (e + 1)`, with `e = 1` for the
association energy and `e = 0` otherwise, spelled `dimensionless`, `kelvin`,
`kelvin ** n`, `1 / kelvin` or `1 / kelvin ** n`. Each coefficient is one
`temperature_law_coefficient` active slot selected by its identity; a law owns
its symmetric pair in both directions, and its stored constant and `k_ij`
reciprocal slope are unconsumed. A record repeated under one identity (for
example one volume correction on both cross edges of a component pair) is one
law and one set of slots. Source validity ranges remain the caller's. An
evaluated `k_ij` must be finite, `l_ij < 1`, association energy and volume
nonnegative, at every state: a law leaving that domain is refused with its
target and temperature, never clamped.

A combining-rule edge that declares a `volume_correction`, or whose self edges
carry a law, is evaluated at each temperature from its explicit self edges:
[EqID: `association_volume_correction`](../equations.md#association_volume_correction),
the Uyan2015 Eq. (9) convention with `(1 - c)` outside the square root and the
size factor applied once by `wolbach-sandler`. Its self-edge constants keep
their derivative chain; its own energy and volume are not coordinates. Rule
edges without these declarations keep construction-time resolution. Held2014's
corrected water/1-butanol cross energy (2420.53610 K) is representable as an
explicit energy; adopting that whole edge, with its Eq. (18) volume, is not.

The association Helmholtz contribution is

Defining relation: [EqID: `association_helmholtz`](../equations.md#association_helmholtz).

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

These are implicit-function identities, not a claim that every displayed
derivative is implemented. An equivalent stationary elimination must reproduce
them on its regular domain. The EOS rejects `X` outside `(0,1]`, a nonfinite
residual, failed convergence, or an ill-conditioned association Jacobian.

## 6. Bulk dielectric and Debye-Huckel

The active record includes charged-slot Debye-Huckel and packing diameters in
Angstrom and dimensionless ion-suppression coefficients. The scalar-generic
permittivity and screening owners consume that record directly. Pressure-closed
permittivity columns include the exact density chain, including the density
topology, and ratios differentiate both states. Model choices, charges, masses
and dielectric source correlations stay fixed. Parameter admission and fitting
execution belong to [ALG-PARAMETERS](algorithms.md#alg-parameters) and
[ALG-FIT](algorithms.md#alg-fit).

The relative-permittivity topology is selected by a complete, source-bearing
parameter family. It is not a caller option and there is no phase label or
runtime fallback. Mixed or incomplete topology records fail during model
resolution.

The `solvent-mole-weighting` declaration selects neutral component identifiers
in `solvents` and requires `normalization` equal to `unnormalized` or
`solvent-normalized`. It consumes only their constant or existing polynomial
relative-permittivity laws. Charged, missing, duplicate or unknown solvent
identifiers and unconsumed dielectric laws are refused. An actually neutral-only
component selection retains the removed-electrolyte convention.
Defining relation: [EqID: `solvent_mole_weighted_bulk_permittivity`](../equations.md#solvent_mole_weighted_bulk_permittivity).
Bülow2021a Eq. (19) supplies the unnormalized relation; normalization is the
explicitly requested alternative. Existing mass and density rules retain their arithmetic.

### 6.1 Figiel composition-only and Ascani mole-fraction topologies

For the bounded Figiel2025 solvent/electrolyte model, first form salt-free
solvent mass fractions from component molar masses:

Defining relation: [EqID: `salt_free_solvent_permittivity`](../equations.md#salt_free_solvent_permittivity).

Figiel2025 Eq. (11) then gives

Defining relation: [EqID: `ion_suppressed_bulk_permittivity`](../equations.md#ion_suppressed_bulk_permittivity), with `alpha_i = C_epsilon` for every ion.

`C_epsilon = 7.01` is a source-bearing model record. The relation is empirical
for water and methanol correlations, was checked against ethanol, uses data up
to about `x_ion = 0.38`, and was evaluated by Figiel2025 at 298.15 K. Other
temperatures require source-backed pure-solvent dielectric correlations and do
not broaden the empirical ion-suppression evidence.

For the bounded Ascani2022 Case Study 2 model, AscaniHeld2021 Eq. 13 instead
defines

Defining relation: [EqID: `component_mixed_bulk_permittivity`](../equations.md#component_mixed_bulk_permittivity).

The source records set `epsilon_r,i = 8` for the admitted alkali-halide ions.
That mixture quantity is not Figiel's radial `epsilon_r,ion`: the Ascani bundle
contains neither the Figiel radial-region model record nor a dielectric
suppression coefficient. Its water and 1-butanol source correlations are
evaluated exactly at the bundle's sole admitted temperature, 298.15 K. CppAD
differentiates the displayed composition rule in the common owner.

With every species charged (an ionic liquid plus a salt), the solvent sum is empty and the rule is `epsilon_r,bulk = sum_i x_i epsilon_r,i` (Bulow2019 Eq. 7; Bulow2021 §2.3, Table 5); the Debye-Huckel term and the Born charging work from vacuum (Bulow2021 Eq. 24) need only this `epsilon_r,bulk > 0`. Only this form admits a parameter set without a neutral species: ion suppression and the density form average neutral-solvent records, and Rueben2024 fit the perturbation form's ion records to water-salt mixtures. Under this form, stability admits ionic trials only on the condensed side of their own composition's loop ([ALG-STABILITY](algorithms.md#supports-seeds-and-trial-solves)).

### 6.2 Schreckenberg density topology

Yang2024 Eqs. (1a), (2a), and (2b) specify the
Schreckenberg/Uematsu--Frank topology, admitted as the permittivity choice
`solvent-density-temperature` with component records `dielectric_volume`
(`liter / mole`) and `dielectric_temperature` (`kelvin`). With `S` the
neutral-solvent set,

Defining relation: [EqID: `density_correlation_mixed_parameters`](../equations.md#density_correlation_mixed_parameters), [EqID: `density_correlated_bulk_permittivity`](../equations.md#density_correlated_bulk_permittivity).

The source water comparison uses `d_V = 0.3777 dm3 mol^-1` and
`d_T = 1403 K` from Yang2024. Because solvent partial density is explicit,
`epsilon_r,bulk -> 1` as `rho_m x_solv -> 0`. Every neutral solvent requires
both parameters; missing records, records on an ion, or a constant solvent
permittivity record under this choice fail closed.

### 6.3 Zuber ion-specific topology

For a single source-supported solvent, Zuber2014 Eq. (9) is

Defining relation: [EqID: `ion_suppressed_bulk_permittivity`](../equations.md#ion_suppressed_bulk_permittivity), with one `alpha_i,s` per ion.

The water/Na+/Cl- comparison uses the Table 2 coefficients
`alpha_Na+,water = 3.62` and `alpha_Cl-,water = 6.80`. These parameters are
ion--solvent specific and are not replacements for Figiel's universal
`C_epsilon = 7.01`. Zuber Eq. (8) requires solvent volume fractions for a
mixed-solvent system. The current topology therefore admits exactly one
neutral solvent and rejects a mixed-solvent selection rather than inventing
volume data or a mixing rule.

### 6.4 Rueben perturbation topology

Rueben2024 Eqs. (4)--(6) are specified in the source's Angstrom number-density
basis, admitted as the permittivity choice `dipolar-perturbation` with
component records `dipole_scaling` (`debye ** 2`), `polarizability_scaling`
(`angstrom ** 3`) and `correlation_integral_parameter` (dimensionless). Let `rho_hat_i = rho_m N_A x_i 10^-30` in particles per Angstrom cubed:

Defining relation: [EqID: `effective_dipolar_density`](../equations.md#effective_dipolar_density), [EqID: `dipolar_correlation_integral`](../equations.md#dipolar_correlation_integral), [EqID: `dipolar_perturbation_bulk_permittivity`](../equations.md#dipolar_perturbation_bulk_permittivity).

`mu2_i`, `alpha_i`, and `a2_i` are the source-backed dipole-squared,
polarizability, and correlation-integral records. The water/Na+/Cl- comparison
uses water `mu2=5.199 D2`, `alpha=0 A3`, `a2=0.1276`, and Rueben2024 PT2 ion
values Na+ `(5.407 D2, 0 A3, 0.215)` and Cl- `(2.563 D2, 0 A3, 0.138)`.
All three records are required for every selected component. The equation gives
`epsilon_r,bulk -> 1` as all partial densities approach zero.

<a id="65-common-electrostatic-owner-and-bounded-comparison"></a>

### 6.5 Common electrostatic owner

The selected `epsilon_r,bulk` is evaluated once by the scalar-generic native
owner and is consumed unchanged by both Debye--Huckel and the selected Born
equation. Temperature, density, composition, and parameter dependence therefore
remain inside the common CppAD tape. The parameter record identifies the generic
bundle-defined formulation, and `State.bulk_relative_permittivity` reports the evaluated
state value.

See [permittivity comparisons on the Figiel parameter base](#permittivity-comparisons-on-the-figiel-parameter-base)
for the retained endpoint and liquid-property evidence.

The defining Debye screening quantity is

Defining relation: [EqID: `debye_screening_parameter`](../equations.md#debye_screening_parameter).

Only charged components enter the following sum. Their finite-size DH diameter
is a separately declared record. Equality to a packing or Born diameter is
allowed only when the selected source prescribes it. Source-specific examples:

```text
Figiel: d_i^DH(T) = 0.88 sigma_i,
Ascani: d_i^DH(T) = sigma_i.
```

The Figiel relation follows its Eqs. (3)--(4) and admitted `u_ii = 0`
convention. AscaniHeld2021 explicitly uses `sigma_i` in the hard-chain,
Debye--Huckel, and Born terms. The DH owner consumes its declared `d_i^DH`; the packing owner consumes
`d_i^pack` and Born consumes `d_i^Born`. A source relation between these
records must be retained explicitly, not inferred universally. Here `sigma_i` and `d_i^DH` use the same length basis and `d_i^DH` is
converted to metres in the SI DH equation. For each charged component,

Defining relation: [EqID: `debye_finite_size_chi`](../equations.md#debye_finite_size_chi).

The residual contribution is

Defining relation: [EqID: `debye_huckel_helmholtz`](../equations.md#debye_huckel_helmholtz).

When `sum_i x_i z_i^2 = 0`, the EOS defines `kappa = 0` and `a^DH = 0`
without evaluating `chi_i`; this is the exact neutral limit. The fitted Born
diameter is used only in the Born term. CppAD differentiates the defining
relations above on the positive-ionic-strength branch; the lab
`dkappa/dx_i` catalogue entry is not a production owner.

## 7. Canonical Born source reconciliation

Charged-slot Born diameters belong to the same scalar active record; the corrected
Born owner consumes them directly. Source correction factors and solvent solvation
factors stay fixed. A declared native hydration hold treats the diameter as a
dependent coordinate through EqID: `dependent_born_diameter_derivative`; the
observed target uses the single-ion x-process J/mol convention. The hold is not
a high-weight residual. Stored-diameter acceptance is separate from the inner
root solve. No cross-family diameter tie or coupled hold is qualified here.

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
of `epsilon_r,ion`. The required canonical primary-source equation is:

Defining relation: [EqID: `born_outer_shell_diameter`](../equations.md#born_outer_shell_diameter), [EqID: `born_residual_helmholtz`](../equations.md#born_residual_helmholtz).

Every `d_i^Born` and `D_i` in this equation is converted from angstrom to metre
before use with the SI factor `e^2/(4 pi epsilon_0 k_B T)`. This equation is a
**corrected source equation**. It matches the radial two-region model, retains
both limiting diameters, and matches the archived candidate kernel.
`epsilon_r,ion = 8` enters through a Figiel/Buelow source-bearing model record.
The native source must not embed 8 or 7.01.

The fixed pairs `(c_shell, c_dielectric) = (0,0), (1,0), (1,1)` recover the
uncorrected, shell-only, and shell-plus-dielectric limits. The redundant
`(0,1)` pair is rejected. Source packets retain the literature-specific
diameter values and provenance, but all resolve through this single equation;
no source or bundle selects a different runtime Born owner.

## 8. Common derivatives and state properties

CppAD records the common `a^r(T,rho_m,x)` equations for fixed component and
site topology. Association uses the implicit sensitivities in Section 5.

Define

Defining relation: [EqID: `compressibility_factor`](../equations.md#compressibility_factor), [EqID: `pressure`](../equations.md#pressure).

The unconstrained composition partials used below follow the Gross2001
constant-`T`, constant-density mixture convention. For component `i`,

Defining relation: [EqID: `residual_chemical_potential`](../equations.md#residual_chemical_potential), [EqID: `log_fugacity_coefficient`](../equations.md#log_fugacity_coefficient).

The chemical-potential convention and fugacity coefficient must be explicit.
It does not publish the lab's nonunique contribution allocation of `ln Z`.

The exact residual caloric identities at the state density are

Defining relation: [EqID: `residual_internal_energy`](../equations.md#residual_internal_energy), [EqID: `residual_enthalpy`](../equations.md#residual_enthalpy), [EqID: `residual_heat_capacity_cv`](../equations.md#residual_heat_capacity_cv).

With `Z_T = (partial Z/partial T)_(rho_m,x)` and
`Z_rho = (partial Z/partial rho_m)_(T,x)`,

Defining relation: [EqID: `residual_heat_capacity_cp`](../equations.md#residual_heat_capacity_cp).

The `-1` removes the ideal-gas `c_p-c_v = R` contribution. The EOS also
reports the same-density residual entropy

Defining relation: [EqID: `residual_entropy`](../equations.md#residual_entropy).

If it reports the departure entropy relative to an ideal gas at the same
`T,P`, the result is named separately and equals `s_V^r/R + ln Z`. No result
named only `residual_entropy` may switch reference convention by state route.

For the current central state result, `c_p^r` is available only when
`D = Z + rho_m Z_rho = P_rho/(RT)` is finite and positive and the resulting
quotient is finite. Positivity is the existing mechanically stable branch rule;
no additional universal slope cutoff is adopted. A numerical accuracy claim at
that state additionally requires a recorded fixed-pressure enthalpy-derivative
ladder on the same branch. A near-zero positive `D` is retained as conditioning
stress outside the qualified claim, not silently rejected or treated as
accurate. P4 owns general per-action conditioning and unavailable behavior.
Explicit-density evaluation does not reject an otherwise finite state solely
because the slope is negative; it reports the other mathematical EOS results
without certifying stability and leaves `c_p^r` unavailable.

On a `derived_properties` request, State also assembles fixed-composition
response and partial molar properties from the same T/V Hessian columns,
with no further derivative evaluation:

Defining relation: [EqID: `isothermal_compressibility`](../equations.md#isothermal_compressibility), [EqID: `isobaric_expansivity`](../equations.md#isobaric_expansivity), [EqID: `heat_capacity_ratio`](../equations.md#heat_capacity_ratio), [EqID: `isentropic_compressibility`](../equations.md#isentropic_compressibility), [EqID: `speed_of_sound`](../equations.md#speed_of_sound), [EqID: `joule_thomson_coefficient`](../equations.md#joule_thomson_coefficient), [EqID: `partial_molar_volume`](../equations.md#partial_molar_volume), [EqID: `partial_molar_enthalpy`](../equations.md#partial_molar_enthalpy), [EqID: `partial_molar_entropy`](../equations.md#partial_molar_entropy).

The same positive-`D` rule admits every one of them. Values that need total
`c_p` or `c_v`, and the partial enthalpy and entropy, follow the
thermochemistry availability of Section 9. Unavailable values stay unavailable;
no floor or residual-only substitute is returned. The six named residual
contributions, dipolar included, sum to `residual_helmholtz`.

## 9. General phase energy and caloric basis

The [phase formulation](../equilibrium/formulation.md) defines amount-volume
coordinates for mixtures, association and electrolytes. Pure nonassociating
states are a restriction of that definition, not a separate thermodynamic
owner. Observable derivative orders and chain rules are specified in
[algorithms](algorithms.md).

Residual enthalpy and internal energy differ by the residual pressure-volume
term; total enthalpy and internal energy differ by the full pressure-volume
term. Missing reference thermochemistry means total calorics are unavailable,
not equal to residual calorics. Standard-pressure and standard-density bases
must be converted consistently; see EqID: `standard_pressure_density_conversion`.
A reaction-consistent reference gauge alone does not determine absolute
enthalpy, entropy or heat capacity. An empty record list declares a species
without ideal-gas data; its Gibbs standard and consumed derivatives are completed
from the reaction chain only when supplied records leave nothing but the ionic
charge gauge (EqID: `reference_species_calorics`). Caloric completion retains
the supplied h/s anchors, independently of the chemical rows' conserved gauge.
Fixed source pressure stays fixed when system pressure moves; system-pressure
references move with it. Their frozen internal energy and heat capacity use the
same Gibbs-conjugate volume (EqID: `gibbs_conjugate_molar_volume`), separately
from EOS mechanical volume. A regular local volume chart does not establish a
globally stable reacting minimum.

## 10. Electrolyte reference states

Figiel2025 Eqs. (17), (19), (20)–(22) define solvation, transfer, single-ion
reference activities and stoichiometric mean ionic activity. Hamer–Wu1972
Eqs. (3.4′), (3.5) define the practical molal osmotic coefficient. The
canonical definitions are EqID: `infinite_dilution_solvation_energy`,
EqID: `transfer_free_energy`, EqID: `single_ion_activity_coefficient`,
EqID: `mean_ionic_activity_coefficient`, EqID: `molality_activity_conversion`
and EqID: `practical_osmotic_coefficient`. The scalar solvent difference uses
the extensive residual Helmholtz Euler identity at prescribed pressure before
the molality division; its exact parameter derivative uses the shared CppAD
fugacity columns. This protects dilute cancellation without changing the
thermodynamic definition or reporting floor.

Declare solvent composition, component order, charges, formula-unit
stoichiometry, molar masses, concentration basis and reference pressure.
Dilute along an exactly electroneutral path with fixed salt-free solvent
composition. Single-ion quantities are path/gauge dependent; charge-neutral
stoichiometric contractions carry the physical comparison. The neutral basis
is EqID: `neutral_reference_stoichiometric_basis`; its existence does not
establish current public availability. Mixed-solvent MIAC does not imply a
unique scalar osmotic coefficient. Multi-salt paths require a separately
specified reference convention.

Properties supplies the raw lower reference partials C_T, C_TT, C_theta,
C_Ttheta and C_TP and the bounded higher-order reference columns in EqID:
`reference_mixed_temperature_pressure_derivatives`, using EqID:
`reference_single_insertion_taylor`. The [F8/F9 notebook](../../../analyses/2026-reactive-caloric-actions/notebook.qmd)
qualifies the stated finite-terminal derivatives and observed dilution changes.
It does not extend reacting caloric/action consumers by implication.

The analytic first insertion is EqID `infinite_dilution_reference`. One ion-free
solvent root and the existing association solution supply its values and exact
T/P/parameter derivatives. Born remains nonzero; DH has zero first insertion.
Neutral solvent mixtures use the same amount insertion. Charged backgrounds
are refused, even if net neutral. The former EqIDs `reference_molality_sequence`
and `reference_limit_terminal_pair` are retained qualification calculations.
No terminal-change error is measured by the analytic method; retain actual
density and site diagnostics without assigning a zero error estimate.

A neutral solute in one pure solvent of an all-neutral model has the
fugacity-based Henry constant of EqID: `henry_constant_px`, read from one unit
reaction row of the same reference owner and analytic solvent insertion. It is
`lim f/x`, not the partial-pressure volatility `p/x`. Normalized
infinite-dilution activity coefficients tend to one and carry no absolute
fugacity scale.

Thermochemical standard states and reaction correlations belong to the common
Engine property/formulation owners; they do not modify the residual EOS.
Their coefficients require complete independent source evidence. Reference
pressure is not the trial pressure unless the declared transformation says so.

## 11. Density closure

Use EqID: `density_residual`, EqID: `pressure_density_derivative` and
EqID: `mechanical_stability`. The root-selection algorithm and failure meaning
are [ALG-DENSITY](algorithms.md#alg-density). A mathematical fixed-density
state can exist on a mechanically unstable branch. A specified-pressure
physical phase must meet the root's mechanical criterion; neither choice
establishes global equilibrium.

## 12. Domain and rejection rules

The optional `cation_mixing` declaration adds the empirical solvent-mass
pair/triple term in EqID `cation_pair_triple_mixing_helmholtz`. It names one
neutral solvent, two distinct cations and one declared anion, with both constant
coefficients supplied explicitly. The term is evaluated in every phase on the
molality basis of that declared neutral solvent, with positive solvent mass and
the declared anion retained. Additional anions are admitted: $\theta_{ij}$ is
independent of the anion, while only the declared $\psi_{ijk}$ is nonzero;
every other anion triple coefficient is zero as a model assumption
(Pitzer and Kim, 1974, Eqs. 8, 9 and 11,
[doi:10.1021/ja00825a004](https://doi.org/10.1021/ja00825a004);
Peiper and Pitzer (1982), J. Chem. Thermodyn. 14, 613–638, Eqs. (7)–(8), Zotero X8GEDSJ2).
The unequal-charge electrostatic $E\theta$ terms are excluded.
Transferring chloride-fitted coefficients to sulfate is unvalidated.
Other neutral components enter this contribution only through
mole fractions and total-mole normalization; their masses are not added to the
declared solvent mass. There is no phase selector or automatic deactivation in
an organic phase. Omitting either cation from
a component selection omits the term; a zero amount in a retained component set
keeps its insertion derivatives when the declared solvent amount is positive.
With no declaration, arithmetic is unchanged.
The [cation-mixing qualification](#cation-mixing-qualification) records
the retained domain and extrapolation limits.

At exact $x_s=0$, if $x_i=0$ or $x_j=0$, define the contribution
as zero on that fixed support. Every derivative within that support, including
thermal, volume and active-coefficient derivatives, is exactly zero. Dense
arrays store zero transverse mixing entries; these are defined storage, not
physical insertion derivatives or limits of the aqueous Hessian. In general
$F_{ij}=2\theta/W+\psi n_k/W^2$ diverges as $W\to0$, so this definition is
not a continuous full-composition extension. Dry states with both cations
present remain refused even with zero coefficients or an underflowing cation
product. The individual exact-zero comparisons introduce no water threshold;
positive-solvent arithmetic and its positive-mass refusal are unchanged.

Existing `IonicInsertion` refusals for absent-ion actions and for ionic values
in an uncharged phase remain. In a dry phase containing one participating
cation, the absent partner's scalar chemical potential and fugacity, including
fugacity actions, are also unavailable with `IonicInsertion`. Separate scalar
insertions when both cations are absent do not authorize mixed insertion
actions. Component support signatures distinguish dry and wet tapes; these
stored zeros cannot authorize continuation across the dry boundary. This
admits neutral pure-organic reference evaluations without qualifying organic
ions or predicting Lithium extraction. Other neutral components retain their
existing effects on the other EOS contributions.

Finite temperature bounds attached to the selected source records are
intersected by the resolver and enforced for every topology. A record whose
source domain is explicitly unknown contributes no invented numerical bound.

See the [Figiel and iodide input domain](#figiel-and-iodide-input-domain)
for the source-specific selection and pressure restriction.

Relative-permittivity behavior is selected by one source-neutral Model Family
choice in the canonical Parameters document. The public builder admits
ion-fraction suppression, component-permittivity mixing, ion-specific suppression,
solvent-only, solvent-density-temperature (§6.2) and dipolar-perturbation (§6.4)
permittivity, plus `none`. A choice consumes only its own records; any other
permittivity record is rejected. Citations to the
papers that support a choice and its coefficients remain in provenance; paper
names do not select runtime implementations. `Parameters.from_mapping(...)`
and `Parameters.from_json(...)` parse the same strict document, which must
provide every coefficient and one complete fixed topology needed by the
selected model. No model choice or topology broadens the temperature,
pressure, composition, or fitted-parameter evidence recorded by its sources.

See [Ascani Case Study 2 input selection](#ascani-case-study-2-input-selection)
for the source-specific pair classifications.

The [Held2008 input selection and reference restriction](#held2008-input-selection-and-reference-restriction)
retains its fail-closed source-state requirement.

The fully dissociated Debye--Huckel model family also accepts explicitly
declared neutral--ion and opposite-charge `k_ij` records when a different
source strategy requires those dispersion interactions. Same-sign ionic
dispersion remains identically suppressed by the EOS and its declared pair
record must be `k_ij = 1`. The Engine does not infer either strategy from a
component identity or application.

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

## 13. Scope and derivative boundaries

EOS supplies residual thermodynamics. Equilibrium, reaction standard states,
reference thermochemistry and fitting have separate responsibilities inside
the same Engine. Their absence in the current implementation is a transfer
gap, not a permanent ban on the accepted target. Solids, critical-point
conditions as an augmented equation system, kinetics and Bjerrum pairing
require separately justified scientific extensions. Residual multipoles other
than the optional dipole-dipole term (ADR 0011) remain retired by ADR 0006.

Residual EOS parameters and reference thermochemistry are separately declared
physical records. Public density/pressure state calls supply temperature and
composition with exactly one of density or pressure; native amount-volume
coordinates retain their defining meaning. Reference calculations additionally require
the declared stoichiometry, basis and standard state. Physical-input resolution
and unit arithmetic belong to the native admission owner.

The model/property owner converts declared reaction correlations to reference
Gibbs offsets and their temperature derivatives; equilibrium consumes those
values consistently. Request value and derivative orders explicitly. Missing
thermochemistry or unsupported actions produce typed unavailable total fields
while preserving valid residual properties. These are required semantics;
concrete class names and installed availability belong to the implementation
plan and public interface, not this mathematical specification.

On an ion-free support, Born and DH are identically zero before reciprocal
permittivity, solvation mixing or screening operations are evaluated. Their
supported temperature, density, neutral-composition and parameter actions,
including mixed actions, are exactly zero. Pure CO₂ has no dielectric solvent
in this model: its zero bulk-permittivity diagnostic supplies no ionic-insertion
limit. The empty electrostatic sums require no replacement permittivity or floor.

With finite positive limiting solvent permittivity, a dilute ionic path has a
linear Born leading term and a DH leading term proportional to ionic strength
to the power 3/2. Born's insertion first derivative is generally nonzero; DH's
insertion second derivative diverges. Neutral-restriction derivatives therefore
do not establish unrestricted composition derivatives. Absent charged ordinary
chemical-potential/fugacity slots and their parameter columns are unavailable
with `IonicInsertion`; present neutral slots and regular scalar calorics remain
available. Pressure closure preserves this admission. Full-species insertion
planes on a zero-inventory face require separate qualification.

Zero amount of an associating species is not permission to add a positive
composition floor. The unweighted mass-action Jacobian may remain regular
where a composition-weighted stationary association Hessian becomes singular.
Retain valid central states and report unsupported derivative directions
separately. See [ALG-ASSOCIATION](algorithms.md#alg-association).

Ordinary physical ionic states require electroneutrality and explicit charges.
Unconstrained composition partials use a declared mathematical extension;
they do not admit a charged physical feed. A neutral subset is selected from
its own components/sites, not unrelated ionic records in its source packet.
Source extrapolation must be identified separately from mathematical domain
validity; convergence does not expand a source's applicability.

Central mixing-energy and entropy terms use the continuous limit of n_i log n_i
at zero component amount. Positive total amount suffices for mole fractions;
logarithmic coordinates and insertion chemical potentials require separate
interior/limit treatment. An incipient phase at zero total amount uses the
independent intensive chart in the equilibrium formulation.

## Appendix: Special cases and numerical qualification

### Figiel input identity and calibration evidence

<a id="adopted-figiel-input-and-model-identity"></a>

The source-complete starting input is the immutable
[figiel-yu-held-2025/3 packet](../../../data/packets/figiel-yu-held-2025/3/packet.toml),
with parameter JSON SHA-256
3d37235cd8217e056f55005bfdf2b65dcf7c4a8621e239bdd40a70c273a9be8c.
Its packet fingerprint is retained in [packet-locks.json](../../../data/packet-locks.json).
Supply the JSON to Parameters.from_json and an ordered component subset;
no paper name selects physics. For the fixed-density equation control, the
subset is (water, sodium-cation, chloride-anion), T = 298.15 K,
rho_m = 56000 mol m^-3, and x = (0.98, 0.01, 0.01) in that order. The
calculated pressure is about 17.1 MPa, outside the packet's reported 1-bar
conditions; it tests admission and equations, not empirical accuracy.

The document declares base: pc-saft with the Gross reference and
dispersion rules of §4; association: general-site with explicit water A/B
sites and no ion sites (§5); electrolyte: born with Debye--Huckel,
same-sign ionic dispersion exclusion, and corrected SSM+DS Born (§§6–7);
and permittivity: ion-fraction-suppression with salt-free solvent mass
averaging and source coefficient 7.01 (§6). Its explicit charge, pair,
association, and correlation records supply units and source locators.
Segment, ionic packing, Debye--Huckel, and fitted Born diameters are distinct
inputs; their source relations are stated in §§4, 6, and 7. Equilibrium and
reported activities use the existing reference conventions in §10.
The packet's reported evidence domain is 298.15 K, 1 bar, and ion mole
fraction up to about 0.38. It does not provide missing ideal caloric inputs.

The [current-wheel admission control](../../../docs/scientific/reviews/modular-baseline-admission-2026-09-22.md)
retains sampled values, derivative outputs, exact input identities, and limits.

Figiel parameter roles and calibration rows remain with the existing
[reproduction source manifest](../../../analyses/reproduction/2025-figiel/inputs/source/figiel-2025-regression-target-ledger.yaml):
the canonical Data packet owns reusable model inputs, while the reproduction
analysis owns paper-specific fitted targets and partitions. Its
[parameter-role ledger](../../../analyses/reproduction/2025-figiel/inputs/source/figiel-2025-parameter-provenance.csv)
(SHA-256 `932e8baa90fcefbaa8c3a8730cdeadd83a4c01f0a3b109f4e4cd0319aee9312b`)
covers 122 source cells from Tables 2–5, their association/segment-count
footnotes, and Eq. 11. It classifies 57 cells as fitted in Figiel, 33 inherited
from cited sources, 13 fixed or assumed, and 19 blank or unpublished. Explicit
published zero interactions stay fitted records; blank cells remain absent.
Table 1 solvent-permittivity values and the fixed ionic-region value are
source-provided inputs in the version-3 Data packet, outside that census and
not regression targets. The parameter document preserves compatible identities:
Table 2 water `f_k` is `component/water/solvation_factor`; Table 3 sodium
`d_i^Born` is `component/sodium-cation/born_diameter`; Table 4 sodium–chloride
`k_ij` is `pair/chloride-anion/sodium-cation/k_ij`; Table 5 sodium–water `k_ij`
is `pair/sodium-cation/water/k_ij`; the Eq. 11 coefficient 7.01 is
`model/relative_permittivity/ion_fraction_suppression_coefficient`; and Table 1
`epsilon_r,ion = 8` is
`model/relative_permittivity/ionic_region_relative_permittivity`. The source
`2B` association choice is materialized as explicit A/B sites, and the Table 3
segment sigma remains distinct from the fitted Born diameter; packing and
Debye–Huckel diameters use the source transformation
`d_i^pack = d_i^DH = 0.88 sigma_i`.

The 407-row [calibration-observation ledger](../../../analyses/reproduction/2025-figiel/inputs/source/figiel-2025-regression-target-ledger.csv)
(SHA-256 `4d6ce439f732d44c6cee08d9d800cef6f844eb297a30ebf935daa1c65403c5e7`)
records each row's role, value and unit, species/solvent, temperature,
pressure, composition basis/value, thermodynamic basis, source locator and
row, extraction or construction, source hash, and uncertainty status. It retains
the published order: SI Table S5 water solvation, NaBr MIAC, aqueous MIAC,
Eq. 19 transfer construction, dielectric correlation, then SI Table S4 density
and osmotic-coefficient data. The current five-ion Born-diameter tracer uses
only five SI Table S5 reported averages as targets and their 27 underlying
literature values as support: `T = 298.15 K`, `P = 1 bar`, `x_ion -> 0`, with
the x-treatment from gas at 1 bar to a hypothetical dilute-ideal solution at
1 bar. The SI reports no pointwise uncertainty.

The ledger classifies 191 rows as direct, 47 table-derived, 90 digitized, 43
constructed, and 36 missing-primary-source. In particular, its 36 NaBr MIAC
rows for methanol and ethanol remain explicitly missing-primary-source: each
is dimensionless at 298.15 K and 1 bar on a formula-unit molality basis, but
the retained lab CSVs have no primary-paper extraction result. Digitized and
constructed values remain marked as such; Eq. 19 constructed solvation rows
are not source observations. The ledger preserves calibration roles and
conditions without asserting a complete replay of the paper's objective,
weights, or uncertainty model. It is source-design evidence, not a new fit or
an extension of the baseline's numerical claim.

### Figiel and iodide input domain

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

### Association source conversions

The MEA nine-species packet from version 5 converts the Fakouri
Baygi 2015 (DOI `10.1016/j.cherd.2014.07.017`, Eq. 4) MEA 2B volume this way;
the [conversion analysis](../../../analyses/2026-mea-association-diameter/README.md)
owns the evidence.

The Figiel2025, Figiel–Yu–Held 2025, Khudaida2026 and MEA nine-species packets
name `wolbach-sandler` from their current versions onward, with its
pure-component self-edge inputs and its provenance for both admitted cross
directions. Their earlier versions name `arithmetic-energy-geometric-volume`
and stay unchanged as immutable records.

### Permittivity comparisons on the Figiel parameter base

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

The diagnostic changes only the permittivity records on the explicit
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
path, which formerly failed on a lone stable root, now completes and gives
`0.56348`. The #27 common-data comparison
([analyses/2026-dielectric-alternatives](../../../analyses/2026-dielectric-alternatives/README.md))
supersedes these single-point diagnostics for water/NaCl at 298.15 K and 1 bar.

Since the retained ePC-SAFT ion, Born, and interaction parameters were not
fitted with the density topologies, this is equation and topology evidence,
not VLE validation or authority to change a caller-owned default. The next
scientific step is a Regression-owned sensitivity/refit study with liquid
dielectric-decrement, density, MIAC, and osmotic anchors before cross-package
equilibrium validation. The density-based topology is worth pursuing; changing
only the Zuber suppression coefficients is not supported as the next priority
by this comparison.

### Cation-mixing qualification

The coefficients are qualified only for retained aqueous water/chloride-salt
cases at 298.15 K and 1 bar; the #658 Yao
calibration covers ionic strength through 8.3 mol/kg. It gives no vapor, dry-salt,
mixed-solvent, nonaqueous, temperature-law, caloric or stability qualification.
In particular, it has no ideal-gas residual limit at fixed composition.
Occupied additional neutrals and organic-phase effects are unqualified
extrapolation.

For an **illustrative**, electroneutral organic-rich composition, take water
$x_s=10^{-3}$, Li and Mg $x_i=x_j=10^{-6}$, chloride $x_k=3\times10^{-6}$,
and the remainder neutral organic. Substituting $M_s=0.01801528$ kg/mol,
$\theta=-0.2205$ kg/mol and $\psi=0.024$ kg²/mol² gives
$\Delta a^r\simeq-2.426\times10^{-8}$, about $-6.0\times10^{-5}$ J/mol
of phase at 298.15 K. The Li chemical-potential contribution is nevertheless
$\Delta\mu_i/(RT)=2\theta m_j+\psi m_jm_k\simeq-0.02426$, about
$-60$ J/mol of inserted Li or $-0.0105$ decade in its activity contribution.
These rounded formula substitutions are neither an Engine run nor a predicted
phase or physical validation. At fixed water molalities the molar energy scales
with $x_s$ while ion chemical potentials stay finite. At fixed nonzero ion
fractions the pair and triple terms instead grow as $x_s^{-1}$ and $x_s^{-2}$
as $x_s\to0$; there is no universal small-organic-phase bound.

### Ascani Case Study 2 input selection

The Ascani Case Study 2 bundle is a separate exact selection of water,
1-butanol, `Na+`, `K+`, and `Cl-` at `T = 298.15 K` and `p = 1 bar`. It requires
all ten unordered `k_ij` and `l_ij` pair classifications explicitly. The
same-sign `Na+/K+` dispersion pair is source-marked not applicable; its
equation-level zero cannot be used as source evidence for another pair. The
bundle also carries the source-derived water/1-butanol cross-association
energy and volume while preserving the distinct `k_ij_hb` provenance. Its
source records preserve the combined-dielectric and uncorrected-Born
provenance, while it and the Figiel records resolve through the same canonical
`born` equation with explicit fixed coefficients. Those source distinctions
remain packet provenance and do not create runtime variants or alter packet
fingerprints.

### Held2008 input selection and reference restriction

The Held2008 bundle is a source-complete selection of water plus all 29 Table 2
ions at `T = 298.15 K` and `p = 1 bar`. It requires explicit water--ion
`k_ij = 0` and ion--ion `k_ij = 1`, so no ion--ion dispersion is silently
introduced. It uses the source's fully dissociated Debye--Huckel formulation
with concentration-independent water permittivity and has no Born term,
ion-pairing, or reaction-equilibrium behavior. Evaluation or the
standard-reference evaluator outside the recorded source state is not
source-backed, and the evaluator fails closed outside that state.
