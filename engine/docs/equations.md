# EOS Equation Registry

Canonical, human-editable equation definitions. Edit each equation here once;
`python3 engine/scripts/docs/sync_equation_registry.py` generates the YAML index.
The index is not an implementation-availability declaration. Symbols, source
applicability and conventions are in [full-eos.md](science/full-eos.md) and
[formulation.md](equilibrium/formulation.md); methods and derivative domains
are in [algorithms.md](science/algorithms.md). The [GitHub issue roadmap](https://github.com/tannerpolley/ePC-SAFT-greenfield/milestones)
links the current implementation and numerical evidence gaps.

Each definition has one stable heading, the five metadata lines, and one
LaTeX display block. Preserve this small format so the existing checker can
reject accidental omissions. Source status describes derivation, not acceptance
of a numerical implementation. Historical method equations do not select a
current numerical method.

## `active_binary_interaction_mapping`

- Section: Native phase coordinates
- Status: Implementation convention
- Source: EOS active-coordinate convention for the symmetric binary pair
- Description: Map one active unlike-pair k_ij coordinate into the symmetric binary matrix.
- Change note: Resolved diagonal values remain fixed.

$$
\boldsymbol{k}(\theta)=
\begin{bmatrix}
k_{11}&\theta\\
\theta&k_{22}
\end{bmatrix}.
$$

## `active_binary_size_interaction_mapping`

- Section: Native phase coordinates
- Status: Implementation convention
- Source: EOS active-coordinate convention for the symmetric binary pair
- Description: Map one active unlike-pair l_ij coordinate into the symmetric binary matrix.
- Change note: Resolved diagonal values remain fixed and k_ij remains unchanged.

$$
\boldsymbol{l}(\theta)=
\begin{bmatrix}
l_{11}&\theta\\
\theta&l_{22}
\end{bmatrix}.
$$

## `active_parameter_coordinate_identity`

- Section: Native phase coordinates
- Status: Implementation convention
- Source: Native active-coordinate convention for explicit storage slots
- Description: Define one linear active coordinate and its derivative identity over the slots it owns.
- Change note: A symmetric unlike pair owns both reciprocal off-diagonal slots; fixed slots have zero coordinate derivative.

$$
p_s(\theta_a)=\theta_a\;(s\in S_a),\qquad
\frac{\partial p_s}{\partial\theta_a}=1\;(s\in S_a),\qquad
\frac{\partial F}{\partial\theta_a}=\sum_{s\in S_a}\frac{\partial F}{\partial p_s}.
$$

## `active_parameter_unit_transformation`

- Section: Native phase coordinates
- Status: Exact derivation
- Source: Multivariate chain rule for linear unit conversion
- Description: Transform active-coordinate gradients and mixed Hessians between source and native units.
- Change note: Supported parameter conversions have a multiplicative scale and no additive offset.

$$
p^{\mathrm{native}}=s_p p^{\mathrm{source}},\qquad
\frac{\partial F}{\partial p^{\mathrm{source}}}=s_p\frac{\partial F}{\partial p^{\mathrm{native}}},\qquad
\frac{\partial^2F}{\partial p^{\mathrm{source}}\partial q^{\mathrm{source}}}
=s_p s_q\frac{\partial^2F}{\partial p^{\mathrm{native}}\partial q^{\mathrm{native}}}.
$$

## `angstrom_number_density`

- Section: Symbols, bases, and units
- Status: Exact derivation
- Source: [Gross2001], unit conversion used by Appendix A
- Description: Convert molar density to molecular number density on the angstrom basis.
- Change note: Uses the current EOS Avogadro constant and angstrom conversion.

$$
\widehat{\rho}=\rho_m\,N_A\,10^{-30}.
$$

## `association_contact`

- Section: Contribution intermediates / Association
- Status: Source equation
- Source: [Gross2001], hard-sphere contact rule used by association
- Description: Unlike-site hard-sphere contact function.
- Change note: Mirrors the reduced-diameter association contact expression.

$$
g_{ij}^{hs}=
\frac{1}{1-\eta}
+\frac{3d_id_j}{d_i+d_j}\frac{\zeta_2}{(1-\eta)^2}
+2\left(\frac{d_id_j}{d_i+d_j}\right)^2
\frac{\zeta_2^2}{(1-\eta)^3}.
$$

## `association_helmholtz`

- Section: Contribution intermediates / Association
- Status: Source equation
- Source: [Chapman1990], Eq.~(21)
- Description: Association residual Helmholtz contribution.
- Change note: Uses explicit site multiplicity from the resolved model.

$$
a^{assoc}=\sum_A x_{i(A)}\nu_A
\left(\ln X_A-\frac{X_A}{2}+\frac{1}{2}\right).
$$

## `association_mass_action`

- Section: Contribution intermediates / Association
- Status: Source equation
- Source: [Chapman1990], Eq.~(19)
- Description: Association site-fraction mass-action equation with explicit multiplicity.
- Change note: Written in the zero-residual form used by the EOS.

$$
F_A=X_A\left[
1+\widehat{\rho}\sum_B x_{i(B)}\,\nu_B\,X_B\Delta_{AB}
\right]-1=0.
$$

## `association_strength`

- Section: Contribution intermediates / Association
- Status: Source equation
- Source: [Gross2002], cross-association rule; [Chapman1990], association strength
- Description: Site-pair association strength using explicit source-resolved cross records.
- Change note: Keeps association energy and volume distinct from k_ij and l_ij.

$$
\Delta_{AB}^{ij}=g_{ij}^{hs}
\left[\exp\left(\frac{\epsilon_{AB}^{ij}/k_B}{T}\right)-1\right]
\sigma_{ij}^3\kappa_{AB}^{ij}.
$$

## `born_outer_shell_diameter`

- Section: Born solvation relations
- Status: Source equation
- Source: [Figiel2025], shell-thickness definition preceding Eq.~(7), parameterized to include its uncorrected limit
- Description: Ion-specific outer Born shell diameter with one explicit fixed correction coefficient.
- Change note: The fitted Born diameter, solvation factor, and model-definition coefficient remain distinct inputs.

$$
f_{\mathrm{mix}}=\frac{\sum_kx_kf_k}{\sum_kx_k},\qquad
D_i=d_i^{Born}\left(1+c_{\mathrm{shell}}\frac{f_{\mathrm{mix}}-1}{|z_i|}\right).
$$

## `born_radial_charging_work`

- Section: Born solvation relations
- Status: Source equation
- Source: Born charging work; [Figiel2026Correction], corrected Eq.~(7)
- Description: Canonical radial charging sum with an explicit dielectric-saturation coefficient.
- Change note: The admitted fixed coefficients express the uncorrected, shell-only, and shell-plus-dielectric limits through one equation.

$$
W_{Born}=
\sum_{i\in\mathcal{I}}x_i z_i^2
\left[
\frac{1-\epsilon_{r,\mathrm{bulk}}^{-1}}{D_i}
+c_{\mathrm{dielectric}}\left(1-\epsilon_{r,i}^{-1}\right)
\left(\frac{1}{d_i^{Born}}-\frac{1}{D_i}\right)
\right].
$$

## `born_residual_helmholtz`

- Section: Born solvation relations
- Status: Source equation
- Source: Born charging work; [Figiel2026Correction], corrected Eq.~(7)
- Description: Dimensionless residual Helmholtz contribution from the sole Born charging equation.
- Change note: Shared by fixed and active-coordinate evaluation paths without a formulation selector.

$$
a^{Born}=-\frac{e^2}{4\pi\epsilon_0k_BT}\,W_{Born}.
$$

## `candidate_phase_criteria`

- Section: Equilibrium
- Status: Source equation
- Source: [Pereira2012] system PP with tolerances of Section 4
- Description: A retained point is a candidate stable phase when its lower value is within the gap tolerance of the upper bound, its reduced gradient agrees with the multipliers on every coordinate away from the lower bound, and it is separated from every other candidate in packing fraction or in composition.
- Change note: Gap 5e-2, relative derivative 0.5, packing and composition distinctness 1e-3, as in the source.

$$
\mathrm{UBD} - L^{m} \le \varepsilon_b,\qquad \frac{|\partial g/\partial \tilde{x}_i^{m} - \alpha_i|}{|\alpha_i|} \le \varepsilon_\alpha\ (\tilde{x}_i^{m} > \tilde{x}^{l}),\qquad |\eta^{m} - \eta^{n}| \ge \varepsilon_\eta\ \text{or}\ \max_i |\tilde{x}_i^{m} - \tilde{x}_i^{n}| \ge \varepsilon_x.
$$

## `chemical_potential_stationarity`

- Section: Equilibrium
- Status: Source equation
- Source: [Pereira2012] Eq. 5; [Gross2001] Eq. (A.33)
- Description: Intensive stationarity in phase composition with balance multipliers λ and per-phase charge multipliers φ_p; the same row for a free and a pinned phase.
- Change note: Written intensively (divided by N_p) so a pinned phase carries the same row; phase equilibrium and reaction affinity follow from it.

$$
\mu_{p,i}-\boldsymbol{\lambda}\cdot\mathbf{B}_i-z_i\varphi_p=0,\quad i\in S_p.
$$

## `close_packing_fraction`

- Section: Density and packing
- Status: Exact derivation
- Source: Geometric close packing of equal spheres
- Description: EOS hard upper packing-domain boundary.
- Change note: Mirrors close_packing_fraction without introducing solver policy.

$$
\eta_{\mathrm{cp}}=\frac{\pi}{3\sqrt{2}}.
$$

## `component_mixed_bulk_permittivity`

- Section: Electrolyte dielectric and Debye--Huckel relations
- Status: Source equation
- Source: [AscaniHeld2021], Eq.~(13)
- Description: Ascani combined solvent-mass and ion-mole bulk relative permittivity.
- Change note: Ion mixture permittivities remain distinct from Figiel radial-region permittivity.

$$
\epsilon_{r,\mathrm{bulk}}=
\epsilon_{r,\mathrm{sf}}\sum_{s\in\mathcal{S}}x_s
+\sum_{i\in\mathcal{I}}x_i\epsilon_{r,i}.
$$

## `component_permittivity_temperature_polynomial`

- Section: Electrolyte dielectric and Debye--Huckel relations
- Status: Implementation convention
- Source: Typed constitutive correlation
- Description: Reference-temperature polynomial for one dielectric solvent's relative static permittivity.
- Change note: Each continuous coefficient is independently selectable while the polynomial powers and reference temperature remain fixed.

$$
\epsilon_{r,s}(T)=c_{0,s}+\sum_{n=1}^{N_s}c_{n,s}
(T-T_{\mathrm{ref},s})^n.
$$

## `compressibility_factor`

- Section: Thermodynamic properties
- Status: Source equation
- Source: [Gross2001], Eq.~(A.24)
- Description: Compressibility factor from the residual Helmholtz density derivative.
- Change note: Uses molar density as the independent density coordinate.

$$
Z=1+\rho_m a_\rho,
\qquad
a_\rho=\left(\frac{\partial a^r}{\partial\rho_m}\right)_{T,\boldsymbol{x}}.
$$

## `cp_ideal_constant`

- Section: Reference thermochemistry / Ideal-gas heat-capacity kinds
- Status: Implementation convention
- Source: Declared simplification; no literature form
- Description: A constant ideal-gas heat capacity for a species with no better data, declared as such; its integrals are linear and logarithmic.
- Change note: New declared kind; the only kind without a literature source, and the packet must say so.

$$
C^{\circ}_p=c,\qquad h^{\circ}(T)=h^{\circ}_{\mathrm{ref}}+c\,(T-T_{\mathrm{ref}}),\qquad s^{\circ}(T,P^{\circ})=s^{\circ}_{\mathrm{ref}}+c\ln\frac{T}{T_{\mathrm{ref}}}
$$

## `cp_ideal_dippr_107`

- Section: Reference thermochemistry / Ideal-gas heat-capacity kinds
- Status: Bounded empirical relation
- Source: Aly and Lee 1981, Fluid Phase Equilibria 6, 169 (DIPPR equation 107)
- Description: Ideal-gas heat capacity of one species, hyperbolic form, with its closed-form enthalpy and entropy integral differences; valid on the declared temperature domain only.
- Change note: The admitted real extension uses $c=|C|/T$ and $e=|E|/T$ in the even terms. Anchored differences use stable equivalents; $C=0$ gives $B\Delta T$ and $B\ln(T/T_{\mathrm{ref}})$, while $E=0$ gives no $D$ contribution.

$$
\begin{aligned} c&=|C|/T,\qquad e=|E|/T,\\ C^{\circ}_p&=A+B\left[\frac{c}{\sinh c}\right]^2+D\left[\frac{e}{\cosh e}\right]^2,\\ \Delta h&=A\Delta T+B|C|\left[\coth c\right]_{T_{\mathrm{ref}}}^{T}-D|E|\left[\tanh e\right]_{T_{\mathrm{ref}}}^{T},\\ \Delta s&=A\ln\frac{T}{T_{\mathrm{ref}}}+B\left[c\coth c-\ln\sinh c\right]_{T_{\mathrm{ref}}}^{T}-D\left[e\tanh e-\ln\cosh e\right]_{T_{\mathrm{ref}}}^{T},\end{aligned}\qquad \Delta T=T-T_{\mathrm{ref}},
$$

## `cp_ideal_nasa7`

- Section: Reference thermochemistry / Ideal-gas heat-capacity kinds
- Status: Bounded empirical relation
- Source: NASA RP-1311 Part I (1994), section 4.2, Eqs. (4.6)-(4.8), https://ntrs.nasa.gov/citations/19950013764
- Description: Ideal-gas seven-coefficient NASA7 (original NASA/Chemkin) standard-state thermochemistry on each source-declared temperature interval; coefficients a0-a4 give heat capacity and a5-a6 are enthalpy and entropy integration constants.
- Change note: Uses dimensionless tau=T/(1 K) while preserving the source coefficients' Kelvin convention; P_src is an explicit source standard pressure and P is the caller-declared standard pressure. Do not infer 298.15 K, 1 atm or interval continuity.

$$
\begin{aligned}
\tau&=T/T_{\mathrm u},\qquad T_{\mathrm u}=1\,\mathrm K,\qquad T>0,\qquad T\in I_r,\qquad P^{\circ}>0,\qquad P^{\circ}_{\mathrm{src},r}>0,\\
\frac{C^{\circ}_{p,r}(T)}{R}&=a_{0r}+a_{1r}\tau+a_{2r}\tau^2+a_{3r}\tau^3+a_{4r}\tau^4,\\
\frac{h^{\circ}_r(T)}{RT}&=a_{0r}+\frac{a_{1r}}{2}\tau+\frac{a_{2r}}{3}\tau^2+\frac{a_{3r}}{4}\tau^3+\frac{a_{4r}}{5}\tau^4+\frac{a_{5r}}{\tau},\\
\frac{s^{\circ}_r(T,P^{\circ})}{R}&=a_{0r}\ln\tau+a_{1r}\tau+\frac{a_{2r}}{2}\tau^2+\frac{a_{3r}}{3}\tau^3+\frac{a_{4r}}{4}\tau^4+a_{6r}-\ln\frac{P^{\circ}}{P^{\circ}_{\mathrm{src},r}},\\
\frac{\mathrm d h^{\circ}_r}{\mathrm dT}&=C^{\circ}_{p,r},\qquad T\left(\frac{\partial s^{\circ}_r}{\partial T}\right)_{P^{\circ},P^{\circ}_{\mathrm{src},r}}=C^{\circ}_{p,r}.
\end{aligned}
$$

## `cp_ideal_nasa9`

- Section: Reference thermochemistry / Ideal-gas heat-capacity kinds
- Status: Bounded empirical relation
- Source: NASA RP-1311 Part I (1994), section 4.2, Eqs. (4.9)-(4.11); McBride, Zehe and Gordon, NASA/TP-2002-211556 (2002), Empirical Equations for Fitting Thermodynamic Functions, Eqs. (1)-(3), https://ntrs.nasa.gov/citations/19950013764 and https://ntrs.nasa.gov/citations/20020085330
- Description: Ideal-gas nine-coefficient NASA9 (new NASA/Glenn) standard-state thermochemistry on each source-declared temperature interval; coefficients a0-a6 give heat capacity and a7-a8 are enthalpy and entropy integration constants.
- Change note: Uses dimensionless tau=T/(1 K) while preserving the source coefficients' Kelvin convention; P_src is an explicit source standard pressure and P is the caller-declared standard pressure. Evaluate only on declared intervals, preserve source boundary behavior and do not infer 298.15 K or 1 atm.

$$
\begin{aligned}
\tau&=T/T_{\mathrm u},\qquad T_{\mathrm u}=1\,\mathrm K,\qquad T>0,\qquad T\in I_r,\qquad P^{\circ}>0,\qquad P^{\circ}_{\mathrm{src},r}>0,\\
\frac{C^{\circ}_{p,r}(T)}{R}&=a_{0r}\tau^{-2}+a_{1r}\tau^{-1}+a_{2r}+a_{3r}\tau+a_{4r}\tau^2+a_{5r}\tau^3+a_{6r}\tau^4,\\
\frac{h^{\circ}_r(T)}{RT}&=-a_{0r}\tau^{-2}+a_{1r}\frac{\ln\tau}{\tau}+a_{2r}+\frac{a_{3r}}{2}\tau+\frac{a_{4r}}{3}\tau^2+\frac{a_{5r}}{4}\tau^3+\frac{a_{6r}}{5}\tau^4+\frac{a_{7r}}{\tau},\\
\frac{s^{\circ}_r(T,P^{\circ})}{R}&=-\frac{a_{0r}}{2}\tau^{-2}-a_{1r}\tau^{-1}+a_{2r}\ln\tau+a_{3r}\tau+\frac{a_{4r}}{2}\tau^2+\frac{a_{5r}}{3}\tau^3+\frac{a_{6r}}{4}\tau^4+a_{8r}-\ln\frac{P^{\circ}}{P^{\circ}_{\mathrm{src},r}},\\
\frac{\mathrm d h^{\circ}_r}{\mathrm dT}&=C^{\circ}_{p,r},\qquad T\left(\frac{\partial s^{\circ}_r}{\partial T}\right)_{P^{\circ},P^{\circ}_{\mathrm{src},r}}=C^{\circ}_{p,r}.
\end{aligned}
$$

## `cp_ideal_polynomial`

- Section: Reference thermochemistry / Ideal-gas heat-capacity kinds
- Status: Bounded empirical relation
- Source: Poling, Prausnitz and O'Connell 2001, The Properties of Gases and Liquids, 5th ed., Appendix A polynomial form
- Description: Ideal-gas heat capacity as a power series in centered temperature about an explicit reference temperature, with closed-form integrals; valid on the declared temperature domain only.
- Change note: The accepted basis is $(T-T_{\mathrm{ref}})^k$, including the constant form; $T_{\mathrm{ref}}$ anchors both integrals.

$$
C^{\circ}_p=\sum_{k=0}^{n}a_k(T-T_{\mathrm{ref}})^{k},\qquad h^{\circ}(T)=h^{\circ}_{\mathrm{ref}}+\sum_{k=0}^{n}\frac{a_k}{k+1}(T-T_{\mathrm{ref}})^{k+1},\qquad s^{\circ}(T)=s^{\circ}_{\mathrm{ref}}+\sum_{k=0}^{n}a_kT_{\mathrm{ref}}^{k}\int_{0}^{(T-T_{\mathrm{ref}})/T_{\mathrm{ref}}}\frac{z^k}{1+z}\,\mathrm dz
$$

## `cp_ideal_shomate`

- Section: Reference thermochemistry / Ideal-gas heat-capacity kinds
- Status: Bounded empirical relation
- Source: NIST Chemistry WebBook SRD 69, Gas Phase Heat Capacity (Shomate Equation), https://webbook.nist.gov/cgi/cbook.cgi?ID=C7732185&Mask=203, equations and unit definitions
- Description: Use NIST tabulated numerical coefficients and their prescribed units: A-E and G in the heat-capacity/entropy convention, F and H_coef in the kJ/mol enthalpy convention. The enthalpy polynomial returns an increment; formation enthalpy/reference must be supplied separately. Valid only on each declared temperature interval.
- Change note: The explicit 1000 factor converts the NIST kJ/mol enthalpy expression to Engine J/mol. H_coef names the NIST coefficient, not the observable; each source interval retains its fitted increment and supplied formation enthalpy.

$$
\begin{aligned}t&=T/(1000\,\mathrm K),\\ C_p^\circ/[\mathrm{J\,mol^{-1}K^{-1}}]&=A+Bt+Ct^2+Dt^3+E/t^2,\\ [h^\circ(T)-h^\circ(298.15\,\mathrm K)]/[\mathrm{J\,mol^{-1}}]&=10^3[At+Bt^2/2+Ct^3/3+Dt^4/4-E/t+F-H_{\mathrm{coef}}],\\ s^\circ/[\mathrm{J\,mol^{-1}K^{-1}}]&=A\ln t+Bt+Ct^2/2+Dt^3/3-E/(2t^2)+G.\end{aligned}
$$

## `cross_association_energy_correction`

- Section: Native phase coordinates
- Status: Implementation convention
- Source: EOS typed cross-association correction applied to the declared arithmetic energy combining rule
- Description: Correct every matching role-aware cross-association energy while preserving association volume.
- Change note: One unordered component-and-role coordinate owns every multiplicity-equivalent matching edge.

$$
\frac{\epsilon_{AB,r}^{ij}}{k_B}
=\left(1-k_{\mathrm{HB},ij,r}\right)
\frac{\epsilon_{AB,r,0}^{ij}}{k_B},
\qquad
\frac{\partial\left(\epsilon_{AB,r}^{ij}/k_B\right)}
{\partial k_{\mathrm{HB},ij,r}}
=-\frac{\epsilon_{AB,r,0}^{ij}}{k_B}.
$$

## `cross_dispersion_energy`

- Section: Parameter setup
- Status: Source equation
- Source: [Gross2001], Eq.~(A.7)
- Description: Unlike-pair dispersion energy for an applicable pair.
- Change note: Keeps k_ij semantically distinct from l_ij and association data.

$$
\frac{\epsilon_{ij}}{k_B}
=\sqrt{\frac{\epsilon_i}{k_B}\frac{\epsilon_j}{k_B}}\,(1-k_{ij}),
\qquad z_iz_j\le0.
$$

## `cross_segment_diameter`

- Section: Parameter setup
- Status: Source equation
- Source: [Gross2001], Eq.~(A.7), with source-resulted l_ij extension
- Description: Unlike-pair segment diameter including the distinct l_ij family.
- Change note: Makes the current l_ij semantics explicit.

$$
\sigma_{ij}=\frac{\sigma_i+\sigma_j}{2}(1-l_{ij}).
$$

## `debye_finite_size_chi`

- Section: Electrolyte dielectric and Debye--Huckel relations
- Status: Source equation
- Source: [Bulow2019], finite-size function in Eq.~(2)
- Description: Exact finite-size DH function and its Taylor series at zero, truncated at degree 24 in the implementation. The polynomial is not the exact function away from the limit.
- Change note: A recording-point reduced diameter y_ref may select a cancellation-safe implementation, but is an algorithm choice, not a thermodynamic variable. With the switch at y_ref = 0.2, chi and its first and second Taylor derivatives match 50-digit references to 3e-14, 7e-13 and 1.1e-11 relative on both branches (2026-09-23).

$$
\chi(y)=\frac{3[\ln(1+y)-y+y^2/2]}{y^3}\quad(y>0),\qquad \chi(0)=1,\qquad y=\kappa d_i^{DH},\qquad \chi(y)=\sum_{k=0}^{24}\frac{3(-y)^k}{k+3}+O(y^{25})\quad(y\to0).
$$

## `debye_huckel_helmholtz`

- Section: Electrolyte dielectric and Debye--Huckel relations
- Status: Source equation
- Source: [Bulow2019], Eq.~(2)
- Description: Finite-size Debye--Huckel residual Helmholtz contribution.
- Change note: d_i^DH is the separately declared finite-size diameter in metres; equality to packing or Born diameters requires source justification.

$$
a^{DH}=-\frac{\kappa e^2}{12\pi\epsilon_0\epsilon_{r,\mathrm{bulk}}k_BT}\sum_{i\in\mathcal I}x_i z_i^2\chi(\kappa d_i^{DH}).
$$

## `debye_screening_parameter`

- Section: Electrolyte dielectric and Debye--Huckel relations
- Status: Source equation
- Source: [Bulow2019], Debye screening definition associated with Eq.~(2)
- Description: SI Debye screening parameter for the resolved bulk dielectric.
- Change note: Uses molar density and the current EOS physical constants.

$$
\kappa=
\sqrt{\frac{\rho_mN_Ae^2S_z}
{\epsilon_0\epsilon_{r,\mathrm{bulk}}k_BT}}.
$$

## `declared_solvent_source_reference_transfer`

- Section: Native phase coordinates
- Status: Documentation-only
- Source: Unqualified source-activity conversion; reference-pressure consistency must be established in Properties reference thermochemistry (M3)
- Description: Trial-pressure fugacity/activity contraction; this expression is not an accepted conversion of a fixed source standard state. P in the displayed formula means P_trial, not P0.
- Change note: Candidate for checking the reference-state discrepancy; not an adopted conversion. Do not infer that a source reference pressure is only provenance, or substitute trial pressure for fixed P0. The accepted pressure-density conversion has a separate ID.

$$
\begin{aligned}
\Delta_r^{\mathrm{src}\rightarrow\mathrm{EOS}}
&=c_r+s_r+\Delta\nu_r\ln\!\left(\frac{P}{\rho^\circ RT}\right),\\
s_r&=\sum_i b_i^{(r)}\ell_i^{\mathrm{src}},\\
\ell_i^{\mathrm{src}}&=
\begin{cases}
\ln x_i^S, & z_i=0,\ x_i^S>0,\\
\ln(\overline M_S m^\circ), & z_i=0,\ x_i^S=0,\\
\ln(\overline M_S m^\circ), & z_i\ne0.
\end{cases}
\end{aligned}
$$

## `density_correlated_bulk_permittivity`

- Section: Electrolyte dielectric and Debye--Huckel relations
- Status: Bounded empirical relation
- Source: Yang, Kontogeorgis, and de Hemptinne (2024), DOI 10.1016/j.fluid.2024.114103, Eqs. (1a), (2a), and (2b)
- Description: Schreckenberg density-and-temperature bulk relative permittivity.
- Change note: The solvent partial density makes the dilute-density limit explicit.

$$
\epsilon_{r,\mathrm{bulk}}=
1+\rho_m x_{\mathrm{solv}}d_{V,\mathrm{mix}}
\left(\frac{d_{T,\mathrm{mix}}}{T}-1\right).
$$

## `density_correlation_mixed_parameters`

- Section: Electrolyte dielectric and Debye--Huckel relations
- Status: Bounded empirical relation
- Source: Yang, Kontogeorgis, and de Hemptinne (2024), DOI 10.1016/j.fluid.2024.114103, Eqs. (2a)--(2b)
- Description: Solvent-fraction-normalized Schreckenberg mixture parameters.
- Change note: Only neutral-solvent components enter these mixtures.

$$
x_{\mathrm{solv}}=\sum_{s\in\mathcal{S}}x_s,\qquad
d_{V,\mathrm{mix}}=\frac{\sum_{s\in\mathcal{S}}x_sd_{V,s}}{x_{\mathrm{solv}}},
\qquad
d_{T,\mathrm{mix}}=\frac{\sum_{s\in\mathcal{S}}x_sd_{T,s}}{x_{\mathrm{solv}}}.
$$

## `density_residual`

- Section: Density closure and mechanical domain
- Status: Exact derivation
- Source: Exact pressure-closure equation using Eq.~`pressure`
- Description: Scalar residual whose roots are specified-pressure EOS states.
- Change note: Root-search mechanics are not thermodynamic equations.

$$
F(\rho_m;T,\boldsymbol{x},P^\star)
=P_{\mathrm{EOS}}(T,\rho_m,\boldsymbol{x})-P^\star=0.
$$

## `density_upper_bound`

- Section: Density closure and mechanical domain
- Status: Exact derivation
- Source: Exact inversion of Eqs.~`packing_fraction_coefficient` and `close_packing_fraction`
- Description: Model-specific density boundary implied by close packing.
- Change note: Reuses the sole EOS packing geometry owner.

$$
\rho_{\max}=\frac{\eta_{\mathrm{cp}}}{c_\eta(T,\boldsymbol{x})}.
$$

## `dipolar_correlation_integral`

- Section: Electrolyte dielectric and Debye--Huckel relations
- Status: Source equation
- Source: Rueben et al. (2024), DOI 10.1021/acs.jced.4c00347, Eq. (5)
- Description: Rueben correlation integral from component correlation parameters.
- Change note: The correlation parameter is an explicit record for every component.

$$
I=1+\left(\sum_i a_{2,i}x_i\right)\left[\exp(-y)-1\right].
$$

## `dipolar_perturbation_bulk_permittivity`

- Section: Electrolyte dielectric and Debye--Huckel relations
- Status: Source equation
- Source: Rueben et al. (2024), DOI 10.1021/acs.jced.4c00347, Eq. (6)
- Description: Rueben perturbation-theory bulk relative permittivity.
- Change note: The relation tends to one as all component partial densities vanish.

$$
\epsilon_{r,\mathrm{bulk}}=
1+3y\left\{1+y+\left[\frac{17}{16}I-1\right]y^2\right\}.
$$

## `dispersion_chain_coefficients`

- Section: Contribution intermediates / Dispersion
- Status: Source equation
- Source: [Gross2001], Eqs.~(A.13)--(A.14)
- Description: Chain-length interpolation of universal dispersion coefficients.
- Change note: Mirrors the current A0/A1/A2 and B0/B1/B2 arrays.

$$
\begin{aligned}
a_k(\overline m)&=a_{0k}
+\frac{\overline m-1}{\overline m}a_{1k}
+\frac{(\overline m-1)(\overline m-2)}{\overline m^2}a_{2k},\\
b_k(\overline m)&=b_{0k}
+\frac{\overline m-1}{\overline m}b_{1k}
+\frac{(\overline m-1)(\overline m-2)}{\overline m^2}b_{2k}.
\end{aligned}
$$

## `dispersion_compressibility_factor`

- Section: Contribution intermediates / Dispersion
- Status: Source equation
- Source: [Gross2001], Eq.~(A.15)
- Description: Chain compressibility correction multiplying the second-order term.
- Change note: Mirrors the current rational expression.

$$
C_1=\left[
1+\overline m\frac{8\eta-2\eta^2}{(1-\eta)^4}
+(1-\overline m)
\frac{20\eta-27\eta^2+12\eta^3-2\eta^4}
{[(1-\eta)(2-\eta)]^2}
\right]^{-1}.
$$

## `dispersion_first_moment`

- Section: Contribution intermediates / Dispersion
- Status: Source equation
- Source: [Gross2001], Eq.~(A.10) mixture sum
- Description: First unlike-pair dispersion mixture moment.
- Change note: Includes current segment counts and cross diameter.

$$
M_1=\sum_i\sum_j x_ix_jm_im_j\sigma_{ij}^3
\frac{\epsilon_{ij}/k_B}{T}.
$$

## `dispersion_helmholtz`

- Section: Contribution intermediates / Dispersion
- Status: Source equation
- Source: [Gross2001], Eq.~(A.10)
- Description: PC-SAFT first- and second-order dispersion Helmholtz contribution.
- Change note: Reconciled to the current angstrom number-density basis.

$$
a^{disp}=-2\pi\widehat{\rho}I_1M_1
-\pi\widehat{\rho}\,\overline m\,C_1I_2M_2.
$$

## `dispersion_integral_one`

- Section: Contribution intermediates / Dispersion
- Status: Source equation
- Source: [Gross2001], Eq.~(A.11)
- Description: First PC-SAFT dispersion integral correlation.
- Change note: Evaluated through the shared polynomial owner.

$$
I_1(\eta,\overline m)=\sum_{k=0}^{6}a_k(\overline m)\eta^k.
$$

## `dispersion_integral_two`

- Section: Contribution intermediates / Dispersion
- Status: Source equation
- Source: [Gross2001], Eq.~(A.12)
- Description: Second PC-SAFT dispersion integral correlation.
- Change note: Evaluated through the shared polynomial owner.

$$
I_2(\eta,\overline m)=\sum_{k=0}^{6}b_k(\overline m)\eta^k.
$$

## `dispersion_polynomial`

- Section: Contribution intermediates / Dispersion
- Status: Exact derivation
- Source: [Gross2001], polynomial form in Eqs.~(A.11)--(A.12)
- Description: Horner-evaluated sixth-order dispersion polynomial.
- Change note: Documents the real arithmetic performed by polynomial.

$$
\mathcal{P}(\boldsymbol{c},\eta)=\sum_{k=0}^{6}c_k\eta^k.
$$

## `dispersion_second_moment`

- Section: Contribution intermediates / Dispersion
- Status: Source equation
- Source: [Gross2001], Eq.~(A.10) mixture sum
- Description: Second unlike-pair dispersion mixture moment.
- Change note: Includes the squared reduced pair energy.

$$
M_2=\sum_i\sum_j x_ix_jm_im_j\sigma_{ij}^3
\left(\frac{\epsilon_{ij}/k_B}{T}\right)^2.
$$

## `dual_inner_problem`

- Section: Equilibrium
- Status: Source equation
- Source: [Pereira2012] Eqs. 9--11
- Description: Global inner infimum on the admitted phase domain; g is reduced molar Gibbs energy, x0 is feed composition. Reactive systems replace composition by conserved-moiety coordinates.
- Change note: The feed multiplier term is required in the lower bound. A local inner minimum supplies no certified global lower bound.

$$
q(\boldsymbol\alpha)=\inf_{\boldsymbol x,v}\left[g(\boldsymbol x,v)-\boldsymbol\alpha\cdot\boldsymbol x\right],\qquad q(\boldsymbol\alpha)+\boldsymbol\alpha\cdot\boldsymbol x^{0}\le g^{*}(\boldsymbol x^{0}).
$$

## `effective_dipolar_density`

- Section: Electrolyte dielectric and Debye--Huckel relations
- Status: Source equation
- Source: Rueben et al. (2024), DOI 10.1021/acs.jced.4c00347, Eq. (4)
- Description: Rueben effective dipolar-density variable in the Angstrom number-density basis.
- Change note: Dipole-squared and polarizability records are explicit for every component.

$$
y=\frac{4\pi}{9}\sum_i \widehat{\rho}_i
\left(\frac{\mu_i^2\,10^{-19}}{k_BT}+3\alpha_i\right),
\qquad
\widehat{\rho}_i=\rho_mN_Ax_i10^{-30}.
$$

## `electrolyte_phase_coordinate_jacobian`

- Section: Native phase coordinates
- Status: Exact derivation
- Source: Exact Jacobian of Eqs.~`phase_mole_fraction`--`phase_molar_density`
- Description: State-coordinate Jacobian used to transform electrolyte Helmholtz gradients.
- Change note: Retained by the Equilibrium EOS workflow historical CAP-07 proposal; current qualification is linked from its GitHub issue; the generic EOS tape differentiates amount-volume coordinates directly.

$$
\frac{\partial\rho_m}{\partial n_j}=\frac{1}{V},
\quad
\frac{\partial\rho_m}{\partial V}=-\frac{n}{V^2},
\quad
\frac{\partial x_i}{\partial n_j}=\frac{\delta_{ij}-x_i}{n},
\quad
\frac{\partial x_i}{\partial V}=0.
$$

## `electrolyte_phase_hessian_chain_rule`

- Section: Native phase coordinates
- Status: Exact derivation
- Source: Exact second-order multivariable chain rule
- Description: Transform state Helmholtz Hessians into explicit amount-volume coordinates.
- Change note: Retained by the Equilibrium EOS workflow historical CAP-07 proposal; current qualification is linked from its GitHub issue; the generic EOS tape differentiates amount-volume coordinates directly.

$$
\frac{\partial^2(na^r)}{\partial q_\alpha\partial q_\beta}
=n_\alpha a^r_{y_p}y_{p,\beta}
+n_\beta a^r_{y_p}y_{p,\alpha}
+n\left(
a^r_{y_py_q}y_{p,\alpha}y_{q,\beta}
+a^r_{y_p}y_{p,\alpha\beta}
\right),
\quad
\boldsymbol{y}=(\rho_m,\boldsymbol{x}).
$$

## `electroneutral_amount_chart`

- Section: Equilibrium
- Status: Implementation convention
- Source: [Pereira2012]; Perdomo et al. 2025 (Perdomo2025) Section 3
- Description: Electroneutral reduced coordinates for the tangent-plane search: one ion eliminated per charged trial phase.
- Change note: Candidate reduced chart for charged search; not a current solver-availability claim. Eliminated ion has nonzero charge and reconstructed amounts must be nonnegative.

$$
n_{p,e}=-\frac{1}{z_e}\sum_{i\neq e} z_i\,n_{p,i}.
$$

## `fixed_pressure_fugacity_correction`

- Section: Thermodynamic properties
- Status: Documentation-only
- Source: EOS Newton correction using Eq.~`pressure_density_derivative`
- Description: Correct a converged-state log fugacity for its bounded pressure residual.
- Change note: First-order residual-error estimate for comparison with the exact extensive-coordinate pressure constraint; not a replacement for solving that constraint.

$$
\begin{aligned}
\Delta\rho_m&=-\frac{P-P^\star}{P_\rho},\\
\ln\phi_i(P^\star)&\approx
\ln\phi_i(P)+(\ln\phi_i)_\rho\,\Delta\rho_m.
\end{aligned}
$$

## `fixed_pressure_total_derivative`

- Section: Thermodynamic properties
- Status: Exact derivation
- Source: Exact implicit derivative of the pressure constraint
- Description: Total derivative of a property with respect to an active coordinate at fixed pressure.
- Change note: Applied by the generic EOS active-coordinate block, including canonical Born, solvation, and permittivity coordinates.

$$
\left.\frac{dQ}{d\theta}\right|_{P}
=Q_\theta+Q_\rho\left.\frac{d\rho_m}{d\theta}\right|_P,
\qquad
\left.\frac{d\rho_m}{d\theta}\right|_P=-\frac{P_\theta}{P_\rho}.
$$

## `formula_unit_composition`

- Section: Electrolyte reference states and reported activities
- Status: Exact derivation
- Source: [Figiel2025], formula-unit molality reference-state basis
- Description: Convert salt formula-unit molality per kilogram of salt-free solvent to mole fractions.
- Change note: Uses one kilogram of salt-free solvent as the construction basis.

$$
\begin{aligned}
\overline M_s&=\sum_{s\in\mathcal{S}}y_sM_s,
&n_s^{(0)}&=\frac{1}{\overline M_s},\\
n_{\mathrm{tot}}&=n_s^{(0)}+m(\nu_++\nu_-),
&x_s&=\frac{y_sn_s^{(0)}}{n_{\mathrm{tot}}},
&x_\pm&=\frac{\nu_\pm m}{n_{\mathrm{tot}}}.
\end{aligned}
$$

## `hard_chain_helmholtz`

- Section: Contribution intermediates / Hard chain
- Status: Source equation
- Source: [Gross2001], Eq.~(A.16)
- Description: PC-SAFT hard-chain residual Helmholtz contribution.
- Change note: Mirrors the site-independent chain term.

$$
a^{hc}=\overline{m}\,a^{hs}
-\sum_i x_i(m_i-1)\ln g_{ii}^{hs}.
$$

## `hard_sphere_contact`

- Section: Contribution intermediates / Hard chain
- Status: Source equation
- Source: [Gross2001], Eq.~(A.17)
- Description: Hard-sphere radial distribution function at segment contact.
- Change note: Uses the current pure-component contact diameter d_i/2.

$$
g_{ii}^{hs}=
\frac{1}{1-\eta}
+\frac{3(d_i/2)\zeta_2}{(1-\eta)^2}
+\frac{2(d_i/2)^2\zeta_2^2}{(1-\eta)^3}.
$$

## `hard_sphere_helmholtz`

- Section: Contribution intermediates / Hard chain
- Status: Source equation
- Source: [Gross2001], Eqs.~(A.8)--(A.9)
- Description: Boublik-Mansoori-Carnahan-Starling-Leland hard-sphere Helmholtz term.
- Change note: Mirrors the current zeta-based expression.

$$
a^{hs}=\frac{1}{\zeta_0}\left[
\frac{3\zeta_1\zeta_2}{1-\zeta_3}
+\frac{\zeta_2^3}{\zeta_3(1-\zeta_3)^2}
+\left(\frac{\zeta_2^3}{\zeta_3^2}-\zeta_0\right)
\ln(1-\zeta_3)\right].
$$

## `helmholtz_lagrangian_dual`

- Section: Equilibrium
- Status: Documentation-only
- Source: [Pereira2012] problem D and Eq. 6
- Description: Lagrangian dual in independent admissible composition-volume coordinates; a and v are molar Helmholtz energy and volume. Equality to the equilibrium convex envelope requires the source assumptions and an exact global inner infimum.
- Change note: A source method, not an installed discovery claim; local search and finite outer cuts cannot certify a closed global gap.

$$
g^D(\tilde{\boldsymbol x}^{0})=\sup_{\boldsymbol\alpha}\inf_{\tilde{\boldsymbol x},v}\left[g(\tilde{\boldsymbol x},v)+\boldsymbol\alpha\cdot(\tilde{\boldsymbol x}^{0}-\tilde{\boldsymbol x})\right],\qquad g=\frac{a+Pv}{RT}.
$$

## `ideal_gas_density`

- Section: Density closure and mechanical domain
- Status: Exact derivation
- Source: Ideal-gas limiting relation
- Description: Ideal-gas density scale used to construct the density search interval.
- Change note: This is a numerical scale, not a replacement EOS.

$$
\rho_{\mathrm{ideal}}=\frac{P^\star}{RT}.
$$

## `ideal_gas_helmholtz`

- Section: Reference thermochemistry / Ideal side
- Status: Exact derivation
- Source: Michelsen and Mollerup 2007, Thermodynamic Models: Fundamentals and Computational Aspects, ch. 2 (ideal-gas contribution); the translational and mixing term is registered under phase_helmholtz_amount_volume
- Description: Molar reference offset added to the dimensionless ideal density/mixing term using fixed rho0. Per-species pressure-standard thermochemistry must be declared; absent reference data makes total calorics unavailable.
- Change note: Includes pressure-to-density standard-state conversion; a_ref denotes a molar offset, not an independently complete ideal Helmholtz function.

$$
\frac{a^{\circ}_{\mathrm{ref}}}{RT}=\sum_i x_i\left[\frac{h^{\circ}_i(T)}{RT}-\frac{s^{\circ}_i(T,P_0)}{R}+\ln\frac{RT\rho_0}{P_0}\right],\qquad h^{\circ}_i(T)=\Delta_f h^{\circ}_i(T_{\mathrm{ref}})+\int_{T_{\mathrm{ref}}}^{T}C^{\circ}_{p,i}(\tau)\,d\tau,\qquad s^{\circ}_i(T,P_0)=s^{\circ}_i(T_{\mathrm{ref}},P_0)+\int_{T_{\mathrm{ref}}}^{T}\frac{C^{\circ}_{p,i}(\tau)}{\tau}\,d\tau.
$$

## `implicit_second_action`

- Section: Equilibrium
- Status: Exact derivation
- Source: Twice differentiating F(y(theta),theta)=0 on a smooth fixed branch with nonsingular F_y
- Description: a, b and c are independent parameter directions; subscripts on F are partial derivatives and y_a is a total directional action. The third relation uses the augmented variable z=(y,theta), tangents t_a=(y_a,a) and second variations u_ab=(y_ab,0). Along a diagonal path theta(t)=theta_0+t w the coefficients of y(t)=sum_k y_k t^k satisfy the same recurrence with the unknown coefficient set to zero; mixed actions follow by half-sign polarisation of the diagonal paths. Solve linear systems without explicitly inverting F_y.
- Change note: Reconciled against transfer audit 2026-09-06; E5 (2026-09-16) adds the third-order relation and the coefficient recurrence that realises A1 to A3; required mathematics, not an availability claim.

$$
F_y y_a=-F_a,\qquad F_y y_{ab}=-\left(F_{ab}+F_{ay}y_b+F_{yb}y_a+F_{yy}[y_a,y_b]\right),\qquad F_y y_{abc}=-\left(D^3F[t_a,t_b,t_c]+D^2F[u_{ab},t_c]+D^2F[u_{ac},t_b]+D^2F[u_{bc},t_a]\right),\qquad F_y y_k=-[t^k]F\bigl(y(t),\theta(t)\bigr)\big|_{y_k=0}.
$$

## `implicit_sensitivity`

- Section: Equilibrium
- Status: Exact derivation
- Source: Implicit-function theorem on the KKT system
- Description: First implicit sensitivity of F(y,theta)=0 on a smooth fixed branch with a square nonsingular row Jacobian; solve the linear system rather than forming its inverse.
- Change note: Rows made of first energy derivatives generally need second energy derivatives here. Higher-order observables/rows require their own order; see ALG-ACTIONS.

$$
\frac{\mathrm d\boldsymbol y}{\mathrm d\boldsymbol\theta}=-\mathbf{F}_{\boldsymbol y}^{-1}\,\mathbf{F}_{\boldsymbol\theta}.
$$

## `infinite_dilution_solvation_energy`

- Section: Electrolyte reference states and reported activities
- Status: Source equation
- Source: [Figiel2025], infinite-dilution fugacity-to-solvation-energy relation
- Description: Single-ion infinite-dilution solvation Gibbs energy on the EOS reference convention.
- Change note: Uses the terminal reference-sequence log fugacity.

$$
\Delta G_{i,\mathrm{solv}}^\infty
=RT\ln\phi_i^\infty.
$$

## `interphase_potential_equality`

- Section: Equilibrium
- Status: Exact derivation
- Source: Consequence of chemical_potential_stationarity
- Description: Equal chemical potential for neutral species and equal electrochemical potential for ions in two supports.
- Change note: Not imposed a second time.

$$
\mu_{p,i}-\mu_{q,i}=z_i\,(\varphi_p-\varphi_q).
$$

## `ion_fraction_suppressed_bulk_permittivity`

- Section: Electrolyte dielectric and Debye--Huckel relations
- Status: Bounded empirical relation
- Source: [Figiel2025], Eq.~(11)
- Description: Figiel-family ion-suppressed bulk relative permittivity.
- Change note: The suppression coefficient is an explicit resolved model record.

$$
\epsilon_{r,\mathrm{bulk}}=
\frac{\epsilon_{r,\mathrm{sf}}}
{1+C_\epsilon x_{\mathrm{ion}}},
\qquad
x_{\mathrm{ion}}=\sum_{i\in\mathcal{I}}x_i.
$$

## `ion_specific_suppressed_bulk_permittivity`

- Section: Electrolyte dielectric and Debye--Huckel relations
- Status: Bounded empirical relation
- Source: Zuber et al. (2014), DOI 10.1016/j.fluid.2014.05.037, Eq. (9)
- Description: Zuber single-solvent ion-specific dielectric suppression.
- Change note: Each ion-solvent suppression coefficient is an explicit resolved record.

$$
\epsilon_{r,\mathrm{bulk}}=
\frac{\epsilon_{r,s}(T)}
{1+\sum_{i\in\mathcal{I}}\alpha_{i,s}x_i}.
$$

## `ion_suppressed_bulk_permittivity`

- Section: Electrolyte dielectric and Debye--Huckel relations
- Status: Source equation
- Source: Zuber et al. 2014 per-ion suppression coefficients (Zuber2014); [Figiel2025] is the special case with one coefficient for every ion
- Description: Bulk relative permittivity of an electrolyte solution as the salt-free solvent-mixture permittivity suppressed by the ion mole fractions with per-ion coefficients.
- Change note: One equation covers ion-specific suppression, ion-fraction suppression (all alpha equal) and solvent-only (alpha = 0); the model layer fills the coefficients.

$$
\varepsilon=\frac{\varepsilon_{\mathrm{sf}}}{1+\sum_i\alpha_i x_i}.
$$

## `ionic_packing_diameter`

- Section: Parameter setup
- Status: Implementation convention
- Source: Declared ion-specific hard-chain diameter; source packets retain the originating convention
- Description: Ionic hard-chain packing diameter, independent of Debye--Huckel and Born diameters.
- Change note: Removes the Born-selector dependency and keeps packing, Debye--Huckel, and Born lengths independently declarable.

$$
d_i=d_i^{\mathrm{pack}},\qquad z_i\ne0,\qquad d_i^{\mathrm{pack}}>0.
$$

## `ionic_strength_charge_sum`

- Section: Electrolyte dielectric and Debye--Huckel relations
- Status: Source equation
- Source: [Bulow2019], charge-weighted composition in the screening definition
- Description: Mole-fraction charge-square sum entering Debye screening.
- Change note: Only charged components contribute.

$$
S_z=\sum_{i\in\mathcal{I}}x_i z_i^2.
$$

## `log_fugacity_coefficient`

- Section: Thermodynamic properties
- Status: Source equation
- Source: [Gross2001], Eq.~(A.32)
- Description: Natural logarithm of the component fugacity coefficient.
- Change note: Available only when the compressibility factor is positive.

$$
\ln\phi_i=\frac{\mu_i^r}{RT}-\ln Z.
$$

## `log_fugacity_density_derivative`

- Section: Thermodynamic properties
- Status: Documentation-only
- Source: Exact derivative of Eqs.~`residual_chemical_potential` and `log_fugacity_coefficient`
- Description: Density derivative of log fugacity used by active fixed-pressure observations.
- Change note: Independent density-coordinate identity for checking the generic extensive-coordinate Hessian.

$$
\begin{aligned}
\left(\ln\phi_i\right)_\rho
&=\left(\frac{\mu_i^r}{RT}\right)_\rho-\frac{Z_\rho}{Z},\\
\left(\frac{\mu_i^r}{RT}\right)_\rho
&=2a_\rho+\rho_ma_{\rho\rho}
+a_{x_i\rho}-\sum_jx_ja_{x_j\rho}.
\end{aligned}
$$

## `material_balance`

- Section: Equilibrium
- Status: Source equation
- Source: [Pereira2012] Eq. 3
- Description: Conserved-moiety balance over the phases, each weighted by its amount; B is the identity when no reaction is declared.
- Change note: Rank of B is the number of balance rows.

$$
\sum_p N_p\,\mathbf{B}\,\mathbf{x}_p=\mathbf{B}\,\mathbf{n}_{\mathrm{feed}}.
$$

## `mean_ionic_activity_coefficient`

- Section: Electrolyte reference states and reported activities
- Status: Source equation
- Source: [Figiel2025], stoichiometric mean ionic activity definition
- Description: Stoichiometric mean ionic activity coefficient on the mole-fraction basis.
- Change note: Supports the one-formula-unit cation/anion topology.

$$
\ln\gamma_\pm^{(x)}
=\frac{\nu_+\ln\gamma_+^{(x)}+\nu_-\ln\gamma_-^{(x)}}
{\nu_++\nu_-}.
$$

## `mean_segment_count`

- Section: Density and packing
- Status: Source equation
- Source: [Gross2001], Appendix A mixture definition
- Description: Composition-averaged segment count.
- Change note: Reconciled to the current residual Helmholtz owner.

$$
\overline{m}=\sum_i x_i m_i.
$$

## `mechanical_stability`

- Section: Density closure and mechanical domain
- Status: Source equation
- Source: Standard isothermal mechanical-stability criterion
- Description: Positive isothermal pressure slope required for a stable density root.
- Change note: Finite positive slope defines local mechanical stability; numerical conditioning is recorded separately.

$$
\left(\frac{\partial P}{\partial\rho_m}\right)_{T,\boldsymbol{x}}
=P_\rho>0.
$$

## `modified_fraction_chart`

- Section: Equilibrium
- Status: Exact derivation
- Source: [Pereira2012] composition chart extended to charged species per Perdomo et al. 2025 (Perdomo2025) Section 3
- Description: Modified mole fractions of a charged phase after eliminating one ion by electroneutrality; the modified fractions sum to one and the eliminated ion is reconstructed from the retained ones.
- Change note: Valid only where each retained factor is nonzero and the reconstructed composition is admissible. Other supports require another regular charge chart; no installed search is implied.

$$
\tilde{x}_i = f_i\,x_i,\qquad f_i = 1 - \frac{q_i}{q_e}\ (i \neq e),\qquad \sum_{i \neq e}\tilde{x}_i = 1,\qquad x_e = -\frac{1}{q_e}\sum_{i \neq e} q_i x_i,\qquad x_i = \frac{\tilde{x}_i}{f_i}.
$$

## `molality_activity_conversion`

- Section: Electrolyte reference states and reported activities
- Status: Source equation
- Source: [Figiel2025], mole-fraction to molality activity conversion
- Description: Convert the mean ionic coefficient to the formula-unit molality basis.
- Change note: Uses salt-free solvent molar mass in kilograms per mole.

$$
\gamma_\pm^{(m)}=
\frac{\gamma_\pm^{(x)}}
{1+\overline M_s\,m(\nu_++\nu_-)}.
$$

## `multiplier_potential_link`

- Section: Equilibrium
- Status: Exact derivation
- Source: [Pereira2012] Eq. 7 in the modified chart of Perdomo et al. 2025 (Perdomo2025)
- Description: Dimensional modified chemical potentials in J/mol are conjugate to modified fractions. Dimensionless dual multiplier alpha is their difference divided by RT; d is the dependent normalized coordinate.
- Change note: The regular chart and admissible composition are required. Dual multipliers are distinct from dimensional balance/charge multipliers and numerical optimizer multipliers.

$$
\tilde\mu_i=\frac{\mu_i-(z_i/z_e)\mu_e}{f_i},\qquad \sum_{i\ne e}\tilde\mu_i\tilde x_i=\sum_i\mu_i x_i,\qquad \alpha_k=\frac{\tilde\mu_k-\tilde\mu_d}{RT}.
$$

## `neutral_reference_stoichiometric_basis`

- Section: Native phase coordinates
- Status: Exact derivation
- Source: Exact charge-neutral basis from explicit stoichiometric declarations
- Description: Native standard-reference rows and log-fugacity contractions.
- Change note: Neutral-component unit rows and a spanning set of cation-anion formula rows are immutable for a resolved model.

$$
\boldsymbol{b}^{(s)}=\boldsymbol{e}_s,
\qquad
\boldsymbol{b}^{(+-)}=\nu_+\boldsymbol{e}_+
+\nu_-\boldsymbol{e}_-,
\qquad
c_r=\sum_i b_i^{(r)}\ln\phi_i^{\mathrm{ref}},
\qquad
\ln\phi_i^{\mathrm{ref}}=
\begin{cases}
\ln\phi_i(\boldsymbol{x}_{\mathrm{neutral}}), & z_i=0,\\
\ln\phi_i^\infty, & z_i\ne 0.
\end{cases}
$$

## `neutral_segment_diameter`

- Section: Parameter setup
- Status: Source equation
- Source: [Gross2001], Eq.~(A.6)
- Description: Temperature-dependent hard-chain diameter for neutral components.
- Change note: Reconciled to the current neutral diameter branch.

$$
d_i(T)=\sigma_i(T)\left[1-0.12\exp\left(-\frac{3(\epsilon_i/k_B)}{T}\right)\right],
\qquad z_i=0.
$$

## `outer_approximation_bound`

- Section: Equilibrium
- Status: Source equation
- Source: [Pereira2012] problem OP; [Pereira2012] Section 2 after Mitsos and Barton 2007 (MitsosBarton2007)
- Description: Finite outer approximation of the semi-infinite dual: maximize a scalar bounded by every retained cut and by the feed free energy; the optimum is an upper bound on the dual value and its multipliers drive the next inner problem.
- Change note: Solved as a linear program over the retained set M; the feed cut has zero slope.

$$
\mathrm{UBD} = \max_{v,\boldsymbol\alpha}\ v\quad \text{s.t.}\quad v \le g(\tilde{\boldsymbol x}^{m}, v^{m}) + \boldsymbol\alpha\cdot(\tilde{\boldsymbol x}^{0} - \tilde{\boldsymbol x}^{m})\ \ \forall m \in \mathcal M,\qquad v \le g(\tilde{\boldsymbol x}^{0}, v^{0}).
$$

## `packing_coefficient_parameter_derivative`

- Section: Native phase coordinates
- Status: Exact derivation
- Source: Chain rule through the packing-fraction coordinate map; the packing fraction definition and the temperature-dependent segment diameter are registered under packing_fraction and the hard-sphere diameter entries
- Description: At fixed T, packing fraction and composition, differentiate c=(pi/6) N sum_i x_i m_i d_i cubed, where N=N_p N_A 10^-30 for phase amount N_p in mol and diameters in angstrom; c has cubic-metre units. The displayed partials specialize the neutral Gross diameter rule.
- Change note: For other diameter correlations or active coefficients differentiate their defining equation as well; these partials do not establish universal zero derivatives for additional families.

$$
\frac{\partial c}{\partial m_i}=\frac{\pi}{6}N x_i d_i^{3},\qquad \frac{\partial c}{\partial \sigma_i}=\frac{\pi}{6}N x_i m_i\,3d_i^{2}\frac{d_i}{\sigma_i},\qquad \frac{\partial c}{\partial(\varepsilon_i/k)}=\frac{\pi}{6}N x_i m_i\,3d_i^{2}\,\sigma_i\frac{0.36}{T}e^{-3\varepsilon_i/(kT)},\qquad \left.\frac{\partial v_p}{\partial\theta_k}\right|_{T,\eta_p,x_p}=\frac{1}{\eta_p}\frac{\partial c}{\partial\theta_k}
$$

## `packing_fraction`

- Section: Density and packing
- Status: Source equation
- Source: [Gross2001], Eq.~(A.8) at n=3
- Description: Packing fraction from molar density and segment geometry.
- Change note: Shared by Helmholtz evaluation and EOS domain owners.

$$
\eta=\zeta_3
=\frac{\pi}{6}\rho_mN_A10^{-30}\sum_i x_i m_i d_i^3.
$$

## `packing_fraction_coefficient`

- Section: Density and packing
- Status: Exact derivation
- Source: Exact specialization of Eq.~`packing_fraction`
- Description: Composition and temperature coefficient satisfying eta=rho_m c_eta.
- Change note: Mirrors the value-only volume-domain helper.

$$
c_\eta(T,\boldsymbol{x})=
\frac{\pi}{6}N_A10^{-30}\sum_i x_i m_i d_i(T)^3,
\qquad \eta=\rho_m c_\eta.
$$

## `packing_moments`

- Section: Density and packing
- Status: Source equation
- Source: [Gross2001], Eq.~(A.8)
- Description: Segment-density moments used by hard-chain and association terms.
- Change note: Uses the current angstrom number-density basis.

$$
\zeta_n=\frac{\pi}{6}\widehat{\rho}\sum_i x_i m_i d_i^n,
\qquad n=0,1,2,3.
$$

## `packing_volume_bounds`

- Section: Native phase coordinates
- Status: Exact derivation
- Source: Exact inversion of eta=c_eta/V_m
- Description: Molar-volume interval corresponding to caller-owned packing-fraction limits.
- Change note: The EOS owns c_eta while the caller owns eta_min and eta_max.

$$
V_{m,\min}=\frac{c_\eta(T,\boldsymbol{x})}{\eta_{\max}},
\qquad
V_{m,\max}=\frac{c_\eta(T,\boldsymbol{x})}{\eta_{\min}}.
$$

## `phase_dimensionless_helmholtz`

- Section: Native phase coordinates
- Status: Exact derivation
- Source: Exact ideal-mixing plus residual Helmholtz construction
- Description: Dimensionless mechanical phase energy. Reference offsets are omitted only for mechanical/nonreactive differences where they cancel; reactive and total caloric calculations require the declared total phase energy.
- Change note: n_ref=1 mol and rho_ref=rho0=1 mol per cubic metre. Phi in mechanical derivative entries denotes this mechanical function.

$$
\Phi_{\mathrm{mech}}\equiv\frac{A_{\mathrm{mech}}}{RTn_{\mathrm{ref}}}=\sum_i\frac{n_i}{n_{\mathrm{ref}}}\left[\ln\frac{n_i/V}{\rho_0}-1\right]+\frac{n}{n_{\mathrm{ref}}}a^r(T,\rho_m,\boldsymbol x).
$$

## `phase_electroneutrality`

- Section: Equilibrium
- Status: Source equation
- Source: [Ascani2022] Eq. 5; Tsanas et al. 2019 Section 2.2.1 and 2022 Section 2.4 (Tsanas2019, Tsanas2022)
- Description: Per-phase electroneutrality for every ion-containing phase except one phase of nonzero amount; the omitted phase's neutrality follows from the balance rows and a neutral feed.
- Change note: One row fewer than the ion-containing phases; if every ion-containing phase is pinned at zero amount, all carry the row.

$$
\mathbf{z}\cdot\mathbf{x}_p=0\quad\text{for all ion-containing }p\text{ except one of nonzero amount}.
$$

## `phase_helmholtz_amount_volume`

- Section: Equilibrium
- Status: Exact derivation
- Source: [Gross2001] Appendix A with the ideal-gas mixing term
- Description: Mechanical phase energy contains ideal density/mixing and residual terms; total phase energy additionally contains the molar reference offset defined by ideal_gas_helmholtz.
- Change note: N_p is mol, V_p is cubic metres, rho0 is the fixed reference molar density. Missing reference thermochemistry does not define total calorics.

$$
\frac{A_{p,\mathrm{mech}}}{RT}=\sum_i n_{p,i}\left[\ln\frac{n_{p,i}}{V_p\rho_0}-1\right]+N_p a^r\left(T,\frac{N_p}{V_p},\frac{\boldsymbol n_p}{N_p}\right),\qquad A_{p,\mathrm{total}}=A_{p,\mathrm{mech}}+N_p a^{\circ}_{\mathrm{ref}}.
$$

## `phase_molar_density`

- Section: Native phase coordinates
- Status: Exact derivation
- Source: Exact extensive-to-intensive coordinate transformation
- Description: Phase molar density from total amount and extensive volume.
- Change note: Uses SI mole and cubic-metre coordinates.

$$
\rho_m=\frac{n}{V}.
$$

## `phase_mole_fraction`

- Section: Native phase coordinates
- Status: Exact derivation
- Source: Exact extensive-to-intensive coordinate transformation
- Description: Phase mole fractions from explicit component amounts.
- Change note: Requires positive total amount n and nonnegative component amounts. Positive individual amounts are required only for their logarithmic interior coordinates; insertion derivatives have separate domains.

$$
x_i=\frac{n_i}{n}.
$$

## `phase_normalization`

- Section: Equilibrium
- Status: Implementation convention
- Source: Composition closure of every declared phase
- Description: Composition normalisation of every phase, free or pinned, so that the coordinates (N_p, x_p, V_p) are uniform and a pin removes exactly one variable.
- Change note: Replaces incipient_normalization, which named a row that is now universal.

$$
\sum_{i\in S_p}x_{p,i}=1\quad\text{for every phase }p.
$$

## `phase_packing_fraction`

- Section: Density and packing
- Status: Exact derivation
- Source: Exact extensive-coordinate form of Eq.~`packing_fraction`
- Description: Packing fraction in explicit component amounts and extensive volume.
- Change note: Mirrors the native phase-domain derivative owner.

$$
\eta(T,\boldsymbol{n},V)=
\frac{\pi}{6}\frac{N_A10^{-30}}{V}\sum_i n_i m_i d_i(T)^3.
$$

## `phase_packing_volume`

- Section: Equilibrium
- Status: Exact derivation
- Source: Packing-fraction definition, [Gross2001] Eq. (A.8)
- Description: Phase volume in cubic metres from amounts in mol, diameters in angstrom and dimensionless packing fraction.
- Change note: Includes molecule-per-mol and cubic-angstrom-to-cubic-metre conversion; matches phase_packing_fraction.

$$
V_p=\frac{\pi N_A10^{-30}}{6\eta_p}\sum_i n_{p,i}m_i d_i^3.
$$

## `phase_pressure`

- Section: Native phase coordinates
- Status: Exact derivation
- Source: Exact extensive Helmholtz derivative
- Description: Phase pressure from the volume gradient of the dimensionless phase Helmholtz function.
- Change note: Returned on the original pascal scale.

$$
P=-RTn_{\mathrm{ref}}
\left(\frac{\partial\Phi}{\partial V}\right)_{T,\boldsymbol{n}}.
$$

## `phase_support`

- Section: Equilibrium
- Status: Implementation convention
- Source: Declared phase support; structural, no source equation
- Description: A species outside a phase's declared support has zero amount in that phase; the coordinate is removed, not bounded.
- Change note: Structural row: it removes a variable rather than adding an equation, so it carries no multiplier and does not enter the count rule.

$$
n_{p,i}=0\quad\text{for }i\notin S_p\ \text{(the variable is absent)}.
$$

## `phase_total_amount`

- Section: Native phase coordinates
- Status: Exact derivation
- Source: Exact extensive-coordinate definition
- Description: Total phase amount from explicit component amounts.
- Change note: Shared by neutral and electrolyte phase blocks.

$$
n=\sum_i n_i.
$$

## `practical_osmotic_coefficient`

- Section: Electrolyte reference states and reported activities
- Status: Source equation
- Source: [HamerWu1972], practical molal osmotic-coefficient definition
- Description: Dimensionless practical molal osmotic coefficient from pure-solvent-referenced solvent activity.
- Change note: Defined only for the admitted single-active-solvent aqueous path.

$$
\Phi=
-\frac{\ln a_s}
{(\nu_++\nu_-)mM_s},
\qquad
\ln a_s=\ln x_s+\ln\phi_s-\ln\phi_s^{\mathrm{pure}}.
$$

## `pressure`

- Section: Thermodynamic properties
- Status: Source equation
- Source: [Gross2001], Eq.~(A.24)
- Description: EOS pressure from compressibility factor.
- Change note: Returns pressure in pascals on the SI molar-density basis.

$$
P=\rho_mRTZ.
$$

## `pressure_density_derivative`

- Section: Thermodynamic properties
- Status: Exact derivation
- Source: Exact derivative of Eq.~`pressure`
- Description: Isothermal pressure derivative used by density closure and stability checks.
- Change note: Transforms CppAD density partials without expanding the EOS derivative.

$$
\left(\frac{\partial P}{\partial\rho_m}\right)_{T,\boldsymbol{x}}
=RT\left(Z+\rho_mZ_\rho\right),
\qquad
Z_\rho=a_\rho+\rho_ma_{\rho\rho}.
$$

## `pressure_state_hessian`

- Section: Thermodynamic properties
- Status: Exact derivation
- Source: Differentiate Eq.~`pressure` with Eq.~`compressibility_factor`, at fixed composition and physical parameters
- Description: Pressure gradient and Hessian on the independent SI temperature/molar-density chart. Here a denotes dimensionless molar residual Helmholtz energy, r denotes molar density, and subscripts denote partial derivatives.
- Change note: Implemented for the bounded single-component neutral slice by `HelmholtzModel::neutral_property` using the existing Helmholtz tape and falsified by `native.properties.issue10` plus the adapter campaign. Association, electrolyte and active-parameter cases remain unavailable; broader use requires a smooth branch and complete implicit chains.

$$
\begin{aligned}
P_T &= R(r+r^2a_r)+RT r^2a_{Tr},\\
P_r &= RT(1+2ra_r+r^2a_{rr}),\\
P_{TT} &= 2Rr^2a_{Tr}+RT r^2a_{TTr},\\
P_{Tr} &= R(1+2ra_r+r^2a_{rr})+RT(2ra_{Tr}+r^2a_{Trr}),\\
P_{rr} &= RT(2a_r+4ra_{rr}+r^2a_{rrr}).
\end{aligned}
$$

## `pressure_equality`

- Section: Equilibrium
- Status: Source equation
- Source: [Pereira2012] Eq. 6; [Gross2001] Eq. (A.24)
- Description: Stationarity in phase volume: every EOS phase is at the system pressure.
- Change note: A stationarity row of the same minimization.

$$
-\frac{\partial A_p}{\partial V_p}=P.
$$

## `pure_active_parameter_mapping`

- Section: Native phase coordinates
- Status: Implementation convention
- Source: EOS parameterized pure-phase coordinate convention
- Description: Construct the pure nonassociating model from active m, sigma, and epsilon-over-k coordinates.
- Change note: No generic mutable parameter overlay is introduced.

$$
\boldsymbol{\theta}_{\mathrm{pure}}
=\left(m,\sigma,\epsilon/k_B\right),
\qquad
\boldsymbol{k}=\boldsymbol{l}=(0).
$$

## `reaction_affinity`

- Section: Equilibrium
- Status: Exact derivation
- Source: Consequence of chemical_potential_stationarity with $\mathbf{B}\boldsymbol\nu=0$
- Description: Vanishing affinity for every reaction in the nullspace of the balance matrix.
- Change note: A problem with B = I has no affinity row.

$$
\sum_i\nu_i\,\mu_{p,i}=\boldsymbol\lambda\cdot(\mathbf{B}\boldsymbol\nu)+\varphi_p\,(\boldsymbol z\cdot\boldsymbol\nu)=0.
$$

## `reaction_activity_quotient`

- Section: Equilibrium / Reaction standard states
- Status: Exact derivation
- Source: Chemical-potential definition on one declared standard-state basis
- Description: Products-positive reaction quotient from dimensionless activities on one identified basis. Activities and the equilibrium constant must use the same basis.
- Change note: Zero activities may have limiting values, but logarithmic insertion derivatives are not implied.

$$
\mu_i=\mu_i^{\circ,b}+RT\ln a_i^{b},\qquad
\ln Q_b=\sum_i\nu_i\ln a_i^{b},\qquad
\ln K_b=-\frac{\sum_i\nu_i\mu_i^{\circ,b}}{RT}.
$$

## `reaction_standard_state_conversion`

- Section: Equilibrium / Reaction standard states
- Status: Exact derivation
- Source: Invariance of chemical potential under a declared standard-state change
- Description: Convert a reaction constant only from explicit per-species standard chemical-potential offsets. Units alone do not define a conversion.
- Change note: Mole-fraction uses explicit dimensionless standard fraction one; molarity and molality require explicit positive standard concentration or molality and the matching coefficient/reference path.

$$
\ln K_c=\ln K_b+\frac{1}{RT}\sum_i\nu_i\left(\mu_i^{\circ,b}-\mu_i^{\circ,c}\right).
$$

## `reaction_standard_state_polynomial_image`

- Section: Equilibrium / Reaction standard states
- Status: Exact derivation
- Source: Integration of $\Delta c_p(T)=\Delta c_{p,0}+\Delta c_p'(T-T_0)$ into Eq.~`reaction_standard_state_conversion`; owner decision 2026-09-21 (#84)
- Description: A standard-state change whose old-minus-new reaction enthalpy, entropy, heat capacity and heat-capacity slope are known at $T_0$ shifts $\ln K$ by a member of the `reaction_constant_ln_polynomial` family, so one correlation form carries the change exactly through every temperature derivative. Zero heat-capacity terms give the constant-enthalpy shift.
- Change note: Its implementation replaces the former constant-enthalpy point-conversion routine; the general identity in `reaction_standard_state_conversion` remains its source. The reference-state pressure path is not part of this identity; a fixed-pressure or system-pressure reference chain contributes separately.

$$
\frac{\Delta g(T)}{RT}=\frac{C_0}{T}+C_1+C_2\ln\frac{T}{T_0}+C_3T,\qquad
C_0=\frac{\Delta h_0-\Delta c_{p,0}T_0+\tfrac12\Delta c_p'T_0^2}{R},\quad
C_1=\frac{\Delta c_{p,0}-\Delta s_0}{R},\quad
C_2=\frac{\Delta c_p'T_0-\Delta c_{p,0}}{R},\quad
C_3=-\frac{\Delta c_p'}{2R},\\
\delta L'=-\frac{\Delta h(T)}{RT^2},\qquad
\delta L''=\frac{2\Delta h(T)}{RT^3}-\frac{\Delta c_p(T)}{RT^2},\qquad
\delta L'''=-\frac{6\Delta h(T)}{RT^4}+\frac{4\Delta c_p(T)}{RT^3}-\frac{\Delta c_p'}{RT^2}.
$$

## `reaction_constant_ln_polynomial`

- Section: Equilibrium / Reaction standard states
- Status: Bounded empirical relation
- Source: Austgen, Rochelle, Peng and Chen 1989, Ind. Eng. Chem. Res. 28, 1060, eq for ln K(T); used for the MEA reactions R1–R3
- Description: Temperature correlation of a reaction equilibrium constant on its declared standard state, with its analytic temperature derivative, which is the van 't Hoff reaction enthalpy on that standard state divided by RT²; valid on the declared temperature domain only.
- Change note: $T_{ref}>0$ is explicit. A legacy $c\ln(T/(1\,\mathrm K))$ coefficient is normalized without changing values by replacing $a_{old}$ with $a=a_{old}+c\ln(T_{ref}/(1\,\mathrm K))$. Standard-state quantities are held fixed in the displayed derivative; temperature-dependent reference conversions contribute separately. Initial Engine regression excludes reaction parameters; application-owned reaction calibration is distinct.

$$
\ln K=a+\frac{b}{T}+c\ln\frac{T}{T_{ref}}+dT,\qquad \frac{d\ln K}{dT}=-\frac{b}{T^2}+\frac{c}{T}+d,\qquad \Delta_r h^{\circ}=RT^2\frac{d\ln K}{dT}
$$

## `reaction_constant_log10_reciprocal`

- Section: Equilibrium / Reaction standard states
- Status: Bounded empirical relation
- Source: The −log10 K = a/T + b + cT form as recorded for MEA reaction R5 in the MEA-Thermodynamics source contract (Aroua 1999, Tong 2012)
- Description: Temperature correlation of a reaction equilibrium constant in decadic-logarithm form, converted to the natural logarithm, with its analytic temperature derivative; valid on the declared temperature domain only.
- Change note: Source notation ln(T) means ln(T/(1 K)); coefficients use the numerical Kelvin convention. Standard-state quantities are held fixed in the displayed derivative; temperature-dependent reference conversions contribute separately. Initial Engine regression excludes reaction parameters; application-owned reaction calibration is distinct.

$$
-\log_{10}K=\frac{a}{T}+b+cT,\qquad \ln K=-\ln 10\left(\frac{a}{T}+b+cT\right),\qquad \frac{d\ln K}{dT}=\ln 10\left(\frac{a}{T^2}-c\right)
$$

## `reduced_hessian_stability`

- Section: Equilibrium
- Status: Exact derivation
- Source: Second-order sufficient condition for a constrained local minimum
- Description: Reduced Lagrangian Hessian over the constraint nullspace is positive definite.
- Change note: Sufficient strict local minimum condition at a regular stationary point; neither necessary in degenerate cases nor a global certificate.

$$
\mathbf{Z}^{\mathsf T}\,\nabla^{2}\mathcal{L}\,\mathbf{Z}\succ 0.
$$

## `neutral_reaction_reference_transfer`

- Section: Equilibrium / Reaction standard states
- Status: Exact derivation
- Source: Chemical-potential equality on the density standard and the dilute molality-to-mole-fraction identity
- Description: Charge-neutral reaction contraction at a declared pure neutral solvent reference. Products have positive stoichiometry; solvent molar mass is in kg/mol and the molality standard in mol/kg. The source basis is explicit.
- Change note: The implemented contraction is a qualified finite terminal estimate at 1e-14 mol/kg, not an exact infinite-dilution limit. The two last ladder differences are observed diagnostics, not rigorous error bounds. Pressure jets and the first temperature derivative differentiate this terminal state; parameter derivatives and caloric reference values are not implied.

$$
\Lambda_i^\rho=\ln(\rho/\rho_0)+\mu_i^{res}/(RT)
=\ln\phi_i+\ln(P/(RT\rho_0)),\quad
C_r=\sum_i\nu_{ri}\Lambda_i^\rho,\quad \sum_i\nu_{ri}z_i=0,\\
A_r=\left(\sum_{i\ne s}\nu_{ri}\right)\ln(M_s m^\circ),\qquad
\ln K_{\rho,r}=\begin{cases}
\ln K_{m,r}+A_r+C_r,&\text{common molality source},\\
\ln K_{x,r}+C_r,&\text{raw mole-fraction source},
\end{cases}\\
n_s=1/M_s,\quad n_i=m r_i\ (i\ne s),\quad \sum_i z_i r_i=0,\quad
\mathcal M=\{10^{-6},10^{-8},10^{-10},10^{-12},10^{-14}\}\,\mathrm{mol\,kg^{-1}},\\
\widehat C_r=C_r(10^{-14}),\qquad
e_r=\max\{|C_r(10^{-14})-C_r(10^{-12})|,|C_r(10^{-12})-C_r(10^{-10})|\}\le 5\times10^{-5}.
$$

## `neutral_reaction_reference_pressure_jet`

- Section: Equilibrium / Reaction standard states
- Status: Exact derivation
- Source: Implicit differentiation of the pressure closure and the neutral reaction reference contraction at fixed temperature and composition
- Description: First and second pressure derivatives of the finite terminal reference. The extensive residual potential is psi=Ares/(RT), the terminal state uses N=1 mol and V=1/rho, and DeltaNu is the stoichiometric mole change.
- Change note: Derivatives follow the same regular positive-slope density branch as the value. A reference at fixed source pressure has zero derivatives with respect to system pressure. E2 and E3 provide the displayed partials; no runtime finite differences are used. Parameter and third pressure derivatives remain unavailable; the first temperature derivative is `standard_state_transfer_temperature_derivative`.

$$
C=-\Delta\nu\ln(V\rho_0)+\sum_i\nu_i\psi_{n_i},\qquad
p_V=RT(-V^{-2}-\psi_{VV}),\quad p_{VV}=RT(2V^{-3}-\psi_{VVV}),\\
C_V=-\Delta\nu/V+\sum_i\nu_i\psi_{n_iV},\qquad
C_{VV}=\Delta\nu/V^2+\sum_i\nu_i\psi_{n_iVV},\\
C_P=\frac{C_V}{p_V},\qquad
C_{PP}=\frac{C_{VV}}{p_V^2}-\frac{C_Vp_{VV}}{p_V^3},\qquad
g_P=-RT C_P,\quad g_{PP}=-RT C_{PP}.
$$

## `standard_state_transfer_temperature_derivative`

- Section: Equilibrium / Reaction standard states
- Status: Exact derivation
- Source: Implicit differentiation of the pressure closure at the fixed terminal composition of the neutral reaction reference; Lambda = ln phi + ln(P/(RT rho0)) and RT^2 d ln phi_i/dT|_P = -hbar_i^res
- Description: First temperature derivative at fixed reference pressure of the finite terminal reference and of the density-basis equilibrium constant. psi=Ares/(RT) is evaluated at N=1 mol and V=1/rho, and the source basis term A is temperature independent. The target is transported as s=-ln K_rho, whose slope maps to species through the same minimum-norm owner as its value.
- Change note: Order one only. The slope is admitted when T times the change over the last three ladder points (1e-10, 1e-12 and 1e-14 mol/kg) is at most 5e-5, the value criterion in RT units. This is an observed convergence criterion, not a limit-error bound: the Debye-Hueckel sqrt(m) form gives a model-dependent multiple M = -3/2(1 + d ln eps/d ln T) - T alpha_P/2 of the value's terminal change (measured 1.54 for the retained aqueous NaCl case), and a steep declared permittivity can exceed the criterion while the value passes. A system-pressure row combines this slope with the pressure jet at first order. Second temperature derivatives, Delta cp, parameter derivatives and caloric reference values are not implied.

$$
p_T=R(V^{-1}-\psi_V)-RT\psi_{VT},\quad V_T=-p_T/p_V,\\
\left.\frac{dC}{dT}\right|_P=\sum_i\nu_i\psi_{n_iT}+C_VV_T,\qquad
RT^2\left.\frac{dC}{dT}\right|_P=-\sum_i\nu_i\left(\bar h_i^{res}+RT\right),\\
\frac{d\ln K_\rho}{dT}=\frac{d\ln K_{src}}{dT}+\left.\frac{dC}{dT}\right|_P,\qquad
\Delta h_\rho=\Delta h_{src}+RT^2\left.\frac{dC}{dT}\right|_P,\qquad
s_r'=-\frac{d\ln K_{\rho,r}}{dT}.
$$

## `repeated_direction_energy_third`

- Section: Helmholtz derivatives / Directional actions
- Status: Exact derivation
- Source: Multivariate Taylor expansion with one direction inserted at order three and the repeated direction at order one
- Description: Cancellation-resistant third directional energy action with two equal directions. Brackets extract a Taylor coefficient, not a finite step evaluation; all derivatives include the same implicit association solution.
- Change note: Pure fifth-order terms cancel and terms with two copies of u start at order six. This evaluates a third derivative, not a public fifth-derivative capability. Normalize directions before coefficient evaluation and restore their three norms afterwards. No association equation or physical parameter changes.

$$
[t^5]F(p+vt\pm ut^3)=\frac{D^5F[v,v,v,v,v]}{120}
\pm\frac{D^3F[u,v,v]}{2},\\
D^3F[u,v,v]=[t^5]\{F(p+vt+ut^3)-F(p+vt-ut^3)\}.
$$

## `reference_limit_terminal_pair`

- Section: Electrolyte reference states and reported activities
- Status: Implementation convention
- Source: Finite terminal-change diagnostic for the declared reference molality sequence.
- Description: Maximum observed terminal change in log fugacity between the last two molality points.
- Change note: Convergence diagnostic; neither an error bound nor an estimate of actual bias without additional convergence assumptions.

$$
e_{\mathrm{ref}}=\max_{i\in\{+,-\}}\left|\ln\varphi_i(m_{k})-\ln\varphi_i(m_{k-1})\right|.
$$

## `reference_molality_sequence`

- Section: Electrolyte reference states and reported activities
- Status: Implementation convention
- Source: EOS convergence sequence for the [Figiel2025] infinite-dilution reference
- Description: Finite molality ladder and terminal estimate of infinite-dilution log fugacity, accepted only under the declared convergence checks.
- Change note: A terminal estimate is not mathematical equality with the infinite-dilution limit; retain observed error diagnostics and reference-path failure.

$$
\mathcal M_{\mathrm{ref}}=\{10^{-6},10^{-8},10^{-10},10^{-12}\}\,\mathrm{mol\,kg^{-1}},\qquad \widehat{\ln\phi_i^\infty}=\ln\phi_i(m=10^{-12}\,\mathrm{mol\,kg^{-1}}).
$$

## `residual_chemical_potential`

- Section: Thermodynamic properties
- Status: Source equation
- Source: [Gross2001], Eq.~(A.33)
- Description: Dimensionless residual chemical potential under unconstrained composition partials.
- Change note: Mirrors the current constant-temperature, constant-density mixture convention.

$$
\frac{\mu_i^r}{RT}
=a^r+(Z-1)+a_{x_i}-\sum_jx_ja_{x_j}.
$$

## `residual_enthalpy`

- Section: Thermodynamic properties
- Status: Source equation
- Source: [Gross2001], Eq.~(A.48)
- Description: Residual molar enthalpy.
- Change note: Includes the residual pressure-volume term Z-1.

$$
h^r=RT\left[-Ta_T+(Z-1)\right].
$$

## `residual_entropy`

- Section: Thermodynamic properties
- Status: Source equation
- Source: [Gross2001], Eq.~(A.46)
- Description: Residual entropy at the same temperature and density.
- Change note: Does not silently change to a same-pressure ideal-gas reference.

$$
s_V^r=-R\left(a^r+Ta_T\right).
$$

## `residual_heat_capacity_cp`

- Section: Thermodynamic properties
- Status: Source equation
- Source: [Gross2001], Eq.~(A.50)
- Description: Residual isobaric molar heat capacity on a conditioned mechanical branch.
- Change note: The final minus one removes the ideal-gas cp-cv contribution.

$$
\frac{c_p^r}{R}
=\frac{c_v^r}{R}
+\frac{(Z+TZ_T)^2}{Z+\rho_mZ_\rho}-1.
$$

## `residual_heat_capacity_cv`

- Section: Thermodynamic properties
- Status: Source equation
- Source: [Gross2001], Eq.~(A.49)
- Description: Residual isochoric molar heat capacity.
- Change note: Transforms CppAD first and second temperature partials.

$$
c_v^r=R\left(-2Ta_T-T^2a_{TT}\right).
$$

## `residual_helmholtz`

- Section: Residual Helmholtz energy
- Status: Implementation convention
- Source: [Gross2001], PC-SAFT residual decomposition extended by admitted electrolyte terms
- Description: Sum of every active residual Helmholtz contribution.
- Change note: Mirrors the single EOS Helmholtz owner and its contribution result.

$$
a^r=a^{hc}+a^{disp}+a^{assoc}+a^{DH}+a^{Born}.
$$

## `residual_internal_energy`

- Section: Thermodynamic properties
- Status: Source equation
- Source: [Gross2001], Eq.~(A.47)
- Description: Residual molar internal energy.
- Change note: Uses the CppAD temperature partial at fixed density and composition.

$$
u^r=-RT^2a_T.
$$

## `salt_free_solvent_permittivity`

- Section: Electrolyte dielectric and Debye--Huckel relations
- Status: Source equation
- Source: [Figiel2025], Eq.~(11) solvent basis; [AscaniHeld2021], Eq.~(13)
- Description: Mass-weighted relative permittivity of the salt-free solvent mixture.
- Change note: Uses component molar masses and includes exactly the neutral components with declared permittivity correlations.

$$
\epsilon_{r,\mathrm{sf}}=
\frac{\sum_{s\in\mathcal{S}}x_sM_s\epsilon_{r,s}}
{\sum_{s\in\mathcal{S}}x_sM_s}.
$$

## `same_sign_ionic_dispersion`

- Section: Parameter setup
- Status: Implementation convention
- Source: [Figiel2025], same-sign ion exclusion; [Ji2012], all-ion-pair dispersion
- Description: Select the source-declared dispersion rule for equally charged ions.
- Change note: The caller selects all ion pairs or same-sign exclusion as part of the electrolyte model family.

$$
\epsilon_{ij}=\begin{cases}\sqrt{\epsilon_i\epsilon_j}(1-k_{ij}),&z_iz_j\leq 0\ \text{or all-ion-pair dispersion},\\0,&z_iz_j>0\ \text{and same-sign exclusion}.\end{cases}
$$

## `segment_sigma_correlation`

- Section: Parameter setup
- Status: Bounded empirical relation
- Source: Source-resulted component correlation represented by the EOS parameter schema
- Description: Evaluate a component segment-size correlation from its constant and exponential terms.
- Change note: Mirrors the general sparse correlation implemented by segment_sigma.

$$
\sigma_i(T)=\sigma_{i,0}+\sum_{q\in\mathcal{Q}_i} A_q\exp(b_qT).
$$

## `single_ion_activity_coefficient`

- Section: Electrolyte reference states and reported activities
- Status: Source equation
- Source: [Figiel2025], electrolyte reference-state fugacity relation
- Description: Single-ion mole-fraction activity coefficient relative to the infinite-dilution fugacity reference.
- Change note: Exposes the exponential of the finite-minus-reference log fugacity.

$$
\gamma_i^{(x)}=
\exp\left(\ln\phi_i-\ln\phi_i^\infty\right).
$$

## `standard_pressure_density_conversion`

- Section: Reference thermochemistry
- Status: Exact derivation
- Source: Exact derivation from the ideal chemical potential and h = -T squared times the temperature derivative of g/T at fixed standard basis
- Description: Fixed P0 and rho0; Delta nu is the sum of reaction stoichiometric coefficients, products positive. The density-basis enthalpy coefficient is an offset; total physical calorics remain basis invariant.
- Change note: Reconciled against transfer audit 2026-09-06; required mathematics, not an availability claim.

$$
g^{\circ}_{i,\rho}=g^{\circ}_{i,P}+RT\ln\frac{RT\rho_0}{P_0},\qquad h^{\circ}_{i,\rho}=h^{\circ}_{i,P}-RT,\qquad s^{\circ}_{i,\rho}=s^{\circ}_{i,P}-R\left(1+\ln\frac{RT\rho_0}{P_0}\right),\\
\ln K_\rho=\ln K_P-\Delta\nu\ln\frac{RT\rho_0}{P_0},\qquad \frac{d\ln K_\rho}{dT}=\frac{d\ln K_P}{dT}-\frac{\Delta\nu}{T}.
$$

## `system_gibbs_helmholtz_coordinates`

- Section: Equilibrium
- Status: Exact derivation
- Source: [Pereira2012] Eqs. 2--4; exact derivation from $G=\min_V[A+PV]$
- Description: System Gibbs objective in Helmholtz coordinates; density-standard offsets appear exactly once and are consistent with ideal_gas_helmholtz.
- Change note: Use either mechanical energy plus offsets or total energy, never both reference terms.

$$
\frac{G}{RT}=\sum_p\left[\frac{A_{p,\mathrm{mech}}}{RT}+\frac{PV_p}{RT}+\sum_i n_{p,i}\frac{g^{\circ}_{i,\rho}}{RT}\right].
$$

## `tangent_plane_distance`

- Section: Equilibrium
- Status: Source equation
- Source: [Pereira2012] Eq. 12; Perdomo et al. 2025 (Perdomo2025) Eq. 40
- Description: Tangent-plane distance in volume-composition space; nonnegativity for every trial is the stability criterion.
- Change note: A negative value is a one-sided instability witness.

$$
\mathrm{TPD}(\boldsymbol x,v)=g(\boldsymbol x,v)-g(\boldsymbol x^{0},v^{0})-\nabla_{\boldsymbol x}g(\boldsymbol x^{0},v^{0})\cdot(\boldsymbol x-\boldsymbol x^{0}).
$$

## `thermochemistry_from_reactions`

- Section: Reference thermochemistry / Declared derivations
- Status: Implementation convention
- Source: Minimum-norm conserved gauge as practised in MEA-Thermodynamics (evaluate_direct_absorption_heat.py, build_thermochemistry), made a declared native derivation
- Description: For an independent reaction-row matrix nu, the minimum-norm solution defines one declared gauge of reaction enthalpies. Adding a conserved-moiety vector preserves reaction enthalpies but does not determine absolute calorics.
- Change note: Gauge is a convention, not replacement evidence for missing enthalpies, entropies or heat capacities. Equilibrium invariance additionally requires consistent Gibbs reference functions.

$$
\boldsymbol{\nu}\,\mathbf{h}^{\circ}=\Delta_r\mathbf{h}^{\circ},\qquad \mathbf{h}^{\circ}=\boldsymbol{\nu}^{\mathsf T}\left(\boldsymbol{\nu}\boldsymbol{\nu}^{\mathsf T}\right)^{-1}\Delta_r\mathbf{h}^{\circ},\qquad \boldsymbol{\nu}\left(\mathbf{h}^{\circ}+\mathbf{B}^{\mathsf T}\mathbf{c}\right)=\Delta_r\mathbf{h}^{\circ}\ \text{for every }\mathbf{c}
$$

## `ideal_enthalpy`

- Section: Reference thermochemistry / Ideal calorics
- Status: Exact derivation
- Source: Ideal-gas mixture calorics on declared component pressure-standard records
- Description: Molar ideal enthalpy is the linear mixture of component standard enthalpies.
- Change note: Every component record must be available on its selected source interval.

$$
h^{id}=\sum_i x_i h_i^{\circ}(T,P_0).
$$

## `ideal_internal_energy`

- Section: Reference thermochemistry / Ideal calorics
- Status: Exact derivation
- Source: Ideal-gas identity $h-u=RT$
- Description: Molar ideal internal energy on the same declared thermochemical basis as ideal enthalpy.
- Change note: No trial-pressure or density-standard offset is inserted into the physical value.

$$
u^{id}=h^{id}-RT.
$$

## `ideal_entropy`

- Section: Reference thermochemistry / Ideal calorics
- Status: Exact derivation
- Source: Ideal-gas mixture entropy at the phase's actual molar density
- Description: Molar ideal entropy from component pressure-standard entropies plus ideal mixing and density terms.
- Change note: Use $x_i\ln x_i\to0$ at zero amount; no logarithmic insertion derivative is promised there.

$$
s^{id}=\sum_i x_i\left[s_i^{\circ}(T,P_0)-R\ln\frac{x_i\rho RT}{P_0}\right].
$$

## `ideal_heat_capacity_cp`

- Section: Reference thermochemistry / Ideal calorics
- Status: Exact derivation
- Source: Ideal-gas mixture calorics
- Description: Molar ideal isobaric heat capacity is the linear mixture of component source correlations.
- Change note: Each source is evaluated only on its selected interval.

$$
c_p^{id}=\sum_i x_i c_{p,i}^{\circ}(T).
$$

## `ideal_heat_capacity_cv`

- Section: Reference thermochemistry / Ideal calorics
- Status: Exact derivation
- Source: Ideal-gas identity $c_p-c_v=R$
- Description: Molar ideal isochoric heat capacity on the same selected source intervals.
- Change note: It inherits ideal-$c_p$ availability.

$$
c_v^{id}=c_p^{id}-R.
$$

## `total_enthalpy`

- Section: Reference thermochemistry / Totals
- Status: Exact derivation
- Source: Michelsen and Mollerup 2007, ch. 2, residual-plus-ideal decomposition at fixed (T, V, n)
- Description: Enthalpy of a phase as the declared ideal reference part plus the residual part the residual EOS definition; reported on the P° basis of the declared thermochemistry. Unavailable without declared reference thermochemistry.
- Change note: Required total-property definition; registration does not establish implementation.

$$
H=\sum_i n_i\,h^{\circ}_i(T)+H^{\mathrm{res}}(T,V,\mathbf{n}),\qquad H^{\mathrm{res}}=A^{\mathrm{res}}-T\left(\frac{\partial A^{\mathrm{res}}}{\partial T}\right)_{V,\mathbf{n}}+PV-nRT
$$

## `total_entropy`

- Section: Reference thermochemistry / Totals
- Status: Exact derivation
- Source: Michelsen and Mollerup 2007, ch. 2, ideal-gas mixture entropy at the phase volume plus the residual at fixed (T, V, n)
- Description: Entropy of a phase as the ideal-gas mixture entropy at the phase's own (T, V, n) plus the residual EOS entropy; reported on the P° basis of the declared thermochemistry, which enters only through the per-species s°. Unavailable without declared reference thermochemistry.
- Change note: Required total-property definition; registration does not establish implementation.

$$
S=\sum_i n_i\left[s^{\circ}_i(T,P^{\circ})-R\ln\frac{n_iRT}{VP^{\circ}}\right]+S^{\mathrm{res}}(T,V,\mathbf{n}),\qquad S^{\mathrm{res}}=-\left(\frac{\partial A^{\mathrm{res}}}{\partial T}\right)_{V,\mathbf{n}}
$$

## `total_internal_energy`

- Section: Reference thermochemistry / Totals
- Status: Exact derivation
- Source: Residual-plus-ideal decomposition at fixed $(T,V,\mathbf n)$
- Description: Molar total internal energy is the ideal value plus the residual EOS value and obeys the full pressure-volume identity.
- Change note: Unavailable without declared reference thermochemistry.

$$
u=u^{id}+u^{res}=h-\frac{P}{\rho}.
$$

## `total_heat_capacity_cp`

- Section: Reference thermochemistry / Totals
- Status: Exact derivation
- Source: Residual-plus-ideal decomposition on a regular fixed-pressure branch
- Description: Molar total isobaric heat capacity combines ideal and residual values at the same state.
- Change note: It is unavailable when either input is unavailable, including a nonpositive mechanical denominator for the residual term.

$$
c_p=c_p^{id}+c_p^{res}.
$$

## `total_heat_capacity_cv`

- Section: Reference thermochemistry / Totals
- Status: Exact derivation
- Source: Residual-plus-ideal decomposition at fixed density
- Description: Molar total isochoric heat capacity combines ideal and residual values at the same state.
- Change note: It does not inherit the fixed-pressure mechanical-stability condition.

$$
c_v=c_v^{id}+c_v^{res}.
$$

## `transfer_free_energy`

- Section: Electrolyte reference states and reported activities
- Status: Source equation
- Source: [Figiel2025], transfer Gibbs-energy relation
- Description: Ion transfer free energy between compatible solvent reference states.
- Change note: Source and target must share temperature, pressure, phase, ions, and stoichiometry.

$$
\Delta G_{i,\mathrm{tr}}^{A\rightarrow B}
=RT\left(\ln\phi_{i,B}^{\infty}-\ln\phi_{i,A}^{\infty}\right).
$$

## `tunneling_transform`

- Section: Equilibrium
- Status: Documentation-only
- Source: [Pereira2012] Eq. 10 and 11
- Description: Exponential tunneling transform of the tangent-plane distance used by the source for the stage-I search from 10 nc random starts. Source search candidate; not an adopted discovery algorithm.
- Change note: Evaluate only within a separately reviewed search design using the shared thermodynamic owner.

$$
f_T(\tilde{\boldsymbol x}, v) = \big(\mathrm{TPD}(\tilde{\boldsymbol x}, v) - d^{*}\big)\exp\!\Big(\frac{p^{*}}{\|\tilde{\boldsymbol x} - \tilde{\boldsymbol x}^{*}\|}\Big),\qquad p^{*} = 10^{-3}.
$$

## `weighted_lagrangian_hessian`

- Section: Numerical methods
- Status: Exact derivation
- Source: Exact differentiation of the constrained objective used by the Ipopt eval_h callback
- Description: All primal coordinates and arbitrary callback multipliers are included. Lambda_opt denotes optimizer constraint multipliers, separate from embedded chemical balance/charge coordinates. Objective scaling never scales constraint curvature.
- Change note: Reconciled against transfer audit 2026-09-06; required mathematics, not an availability claim.

$$
\nabla_y^2\mathcal L_{\mathrm{opt}}=\mathrm{obj\_factor}\,\nabla_y^2 f(y)+\sum_j\Lambda_j^{\mathrm{opt}}\nabla_y^2 g_j(y).
$$
