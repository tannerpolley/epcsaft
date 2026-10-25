# EOS Equation Registry

Canonical, human-editable equation definitions. Edit each equation here once;
`python3 engine/scripts/docs/sync_equation_registry.py` generates the YAML index.
The index is not an implementation-availability declaration. Symbols, source
applicability and conventions are in [full-eos.md](science/full-eos.md) and
[formulation.md](equilibrium/formulation.md); methods and derivative domains
are in [algorithms.md](science/algorithms.md). The [GitHub issue roadmap](https://github.com/tannerpolley/ePC-SAFT/milestones)
links the current implementation and numerical evidence gaps.

Each definition has one stable heading, the five metadata lines, and one
LaTeX display block. Preserve this small format so the existing checker can
reject accidental omissions. Source status describes derivation, not acceptance
of a numerical implementation. Historical method equations do not select a
current numerical method.

## `errors_in_variables_objective`

- Section: Parameter estimation
- Status: Exact derivation
- Source: Zavala (2008), thesis chapter 4, §4.2, Eq. (4.2), printed p. 53 (PDF p. 65); independent physical-unit Gaussian errors selected in #510; the known joint Gaussian covariance (generalized least-squares negative log-likelihood) selected in #572.
- Description: Joint least squares for model parameters v and one standardized correction u per measured input shared by outputs of experiment k, with an optional known correlation R over the positive-weight output errors followed by the active input corrections.
- Change note: sigma_q is finite and nonnegative; zero eliminates its coordinate. Objective sigma_y is finite and positive. Difference residuals use sqrt(w)=a/sigma_y exactly once. Report-only outputs and zero-sigma inputs occupy no row of e or R. A declared covariance C converts through the sigmas; R is fixed, symmetric positive definite with unit diagonal within 64 m epsilon, and never repaired. R = I recovers the independent sum.

$$
q_{kj}=q^M_{kj}+\sigma_{q,kj}u_{kj},\quad
r^x_{kj}=u_{kj},\quad r^y_{ki}=\frac{f_{ki}(\theta(v),X_k(q_k))-y^M_{ki}}{\sigma_{y,ki}},\qquad
\mathbf e=\begin{bmatrix}\mathbf r^y\\\mathbf r^x\end{bmatrix},\quad
R=\Sigma^{-1}C\,\Sigma^{-1}=LL^T,\quad
\min_{v,u}\frac12\lVert L^{-1}\mathbf e\rVert^2
\;\overset{R=I}{=}\;\frac12\sum_{k,i}(r^y_{ki})^2+\frac12\sum_{k,j}(r^x_{kj})^2.
$$

## `heat_capacity_ratio`

- Section: Thermodynamic properties
- Status: Exact derivation
- Source: Exact derivation from `pressure`, `residual_chemical_potential` and the T/V Hessian columns of State; #566 design
- Description: Ratio of the total isobaric and isochoric molar heat capacities.
- Change note: Computed only on request; unavailable unless the total c_p (`total_heat_capacity_cp`) and c_v (`total_heat_capacity_cv`) are available, finite and positive, so it needs ideal thermochemistry.

$$
\gamma=\frac{c_p}{c_v}.
$$

## `henry_constant_px`

- Section: Equilibrium / Reaction standard states
- Status: Exact derivation
- Source: Fugacity definition with the raw density-standard reference of `neutral_reaction_reference_transfer`; #566 design
- Description: Fugacity-based Henry constant in Pa of a neutral solute i in one pure solvent of an all-neutral model: H = lim f_i/x_i, not the partial-pressure volatility p/x.
- Change note: Admitted per neutral solute on one unit reaction row of `neutral_reaction_reference_transfer`, using the analytic solvent insertion of `infinite_dilution_reference`. Mixed solvents and charged models are refused (NonNeutralDeclaration). The molality form M_s H^px is not exposed.

$$
H_i^{px}=\lim_{x_i\to0}\frac{f_i}{x_i}=P\phi_i^\infty=RT\rho_0\exp L_i^\infty,\qquad L_i=\ln\frac{\rho_m}{\rho_0}+\Phi_{n_i},\qquad\rho_0=1\ \mathrm{mol\,m^{-3}}.
$$

## `independent_feed_measurement_coordinates`

- Section: Parameter estimation
- Status: Exact derivation
- Source: Normalization and electroneutrality constraints in the equilibrium formulation; #510 independent species-fraction measurement chart.
- Description: Caller-named independent species fractions and dependent species preserve a fixed feed total N in mol, normalization and charge; B has identity rows for independent species.
- Change note: Fractions are mol/mol; charges z are component charge numbers. One neutral dependent species or two dependent species with distinct charges solve the constraint block. The entire caller-declared box must map into [0,1], with fixed support; no trial normalization or input clipping is permitted.

$$
\mathbf x_k=\mathbf b_k+B_k\mathbf q_k,\quad \mathbf n_k=N_k\mathbf x_k,\quad
\mathbf1^TB_k=\mathbf z^TB_k=0,\quad \mathbf1^T\mathbf b_k=1,\quad\mathbf z^T\mathbf b_k=0,\qquad
\frac{\partial\mathbf n_k}{\partial u_{kj}}=N_kB_{k,:,j}\sigma_{q,kj}.
$$

## `errors_in_variables_residual_derivative`

- Section: Parameter estimation
- Status: Exact derivation
- Source: Differentiation of errors_in_variables_objective and independent_feed_measurement_coordinates; existing exact pressure-closed and solved-equilibrium actions.
- Description: Native first derivatives of standardized measurement residuals; D maps optimizer model parameters to physical parameters; a declared correlation applies the same triangular solve to every complete column.
- Change note: Input-prior rows are exact identities before whitening. T/P actions include finite and reference states; molality actions hold the reference ladder fixed. Ionic T/P reference columns are the terminal-pair limit of `reference_limit_terminal_pair`. Raw prediction columns dy/dv and sigma_q dy/dq come from the physical actions, never from whitened rows.

$$
\frac{\partial r^y_{ki}}{\partial v_a}=\frac1{\sigma_{y,ki}}\sum_b\frac{\partial f_{ki}}{\partial\theta_b}D_{ba},\quad
\frac{\partial r^y_{ki}}{\partial u_{kj}}=\frac{\sigma_{q,kj}}{\sigma_{y,ki}}\frac{\partial f_{ki}}{\partial X_k}\frac{\partial X_k}{\partial q_{kj}},\qquad
\frac{\partial r^x_{kj}}{\partial v_a}=0,\quad\frac{\partial r^x_{kj}}{\partial u_{lh}}=\delta_{kl}\delta_{jh},\qquad
\frac{\partial(L^{-1}\mathbf e)}{\partial(v,u)}=L^{-1}\frac{\partial\mathbf e}{\partial(v,u)}.
$$

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

## `affine_parameter_coordinate_map`

- Section: Parameter regression / Affine coordinates
- Status: Exact derivation
- Source: Multivariate chain rule applied to `active_parameter_coordinate_identity`; issue #209 revision 1.
- Description: A normalized optimizer value v maps to a base value q and proportional physical targets across parameter families. Bounds apply to q. Coefficients carry target-unit/base-unit; reciprocal pairs count once.
- Change note: Nonzero signed coefficients and positive numerical scales are explicit inputs. Row order is independent EOS targets, dependent Born diameters, then reaction coefficients. Reaction coordinates have no EOS targets. Fixed slots remain outside the target list.

$$
q_j=o_j+s_jv_j,\qquad \theta_k=c_kq_j,\qquad
D_{kj}=\frac{\partial\theta_k}{\partial v_j}=c_ks_j,
\qquad s_j>0,\quad c_k\ne0,
$$
$$
J_{\mathrm{opt}}=\operatorname{diag}\!\left(\sqrt{w}\,\frac{dr}{dy}\right)J_{\mathrm{physical}}D,
\qquad C_{\mathrm{physical}}=D C_{\mathrm{optimizer}}D^T.
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
- Source: [Chapman1990], Eq.~(15)
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
- Description: Site-pair association strength using explicit source-resolved cross records. It carries `sigma_ij^3` where Chapman1990 Eq. 18 carries `d^3` (the PC-SAFT parameter convention); volumes fitted using `d^3` are converted in their packets. For a `wolbach-sandler` rule edge, the Gross2002 Eq. 3 size factor makes `sigma_ij` the geometric mean `sqrt(sigma_i sigma_j)(1-l_ij)`.
- Change note: Keeps association energy and volume distinct from k_ij and l_ij.

$$
\Delta_{AB}^{ij}=g_{ij}^{hs}
\left[\exp\left(\frac{\epsilon_{AB}^{ij}/k_B}{T}\right)-1\right]
\sigma_{ij}^3\kappa_{AB}^{ij}.
$$

## `association_volume_correction`

- Section: Contribution intermediates / Association
- Status: Source equation
- Source: [Uyan2015], Eqs.~(8)--(9) and Table 5
- Description: Volume of a combining-rule edge evaluated at the state temperature from its explicit self edges and a declared volume correction c; the Wolbach--Sandler size factor stays in `association_strength`.
- Change note: (1 - c) multiplies the volume once, outside the square root (#567); c = 0 recovers the existing rule.

$$
\frac{\epsilon_{AB}^{ij}}{k_B}=\frac12\left[\frac{\epsilon_{AB}^{ii}(T)}{k_B}+\frac{\epsilon_{AB}^{jj}(T)}{k_B}\right],
\qquad
\kappa_{AB}^{ij}=\sqrt{\kappa_{AB}^{ii}(T)\,\kappa_{AB}^{jj}(T)}\,\bigl[1-c_{ij}(T)\bigr].
$$

## `born_outer_shell_diameter`

- Section: Born solvation relations
- Status: Source equation
- Source: [Figiel2025], Eqs.~(8)--(10), parameterized to include its uncorrected limit
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
+c_{\mathrm{dielectric}}\left(1-\epsilon_{r,\mathrm{ion}}^{-1}\right)
\left(\frac{1}{d_i^{Born}}-\frac{1}{D_i}\right)
\right].
$$

## `born_residual_helmholtz`

- Section: Born solvation relations
- Status: Source equation
- Source: Born charging work; [Figiel2026Correction], corrected Eq.~(7)
- Description: Dimensionless residual Helmholtz contribution from the sole Born charging equation.
- Change note: Shared by fixed and active-coordinate evaluation paths. On an ion-free support the empty sum is identically zero before permittivity or solvation operations; this restriction does not supply ionic insertion derivatives.

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
- Source: [Ascani2022] Eqs. 9--10; [Pereira2012] Eqs. 5--7; [Gross2001] Eq. (A.33)
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
- Source: [AscaniHeld2021], Eq.~(13); all-ion limit [Bulow2019], Eq.~(7), and [Bulow2021], §2.3 and Table 5
- Description: Ascani combined solvent-mass and ion-mole bulk relative permittivity; with no neutral solvent, the solvent sum is empty and the rule is linear in ion mole fractions.
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
- Source: [Gross2001], Eq.~(A.22)
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
- Source: [Gross2001], Eq.~(A.15)
- Description: Unlike-pair dispersion energy for an applicable pair.
- Change note: Keeps k_ij semantically distinct from l_ij and association data. A pair with a declared reference temperature $T_{\mathrm{ref},ij}>0$ takes the record's linear-in-$1/T$ law; its slope $b_{ij}$ (K) is an active-parameter family like $k_{ij}$ (#152). The applicable pairs are the first case of `same_sign_ionic_dispersion`, with C the declared charged molecules (#641).

$$
\frac{\epsilon_{ij}}{k_B}
=\sqrt{\frac{\epsilon_i}{k_B}\frac{\epsilon_j}{k_B}}\,\bigl(1-k_{ij}(T)\bigr),
\qquad k_{ij}(T)=k_{ij}+b_{ij}\Bigl(\frac{1}{T}-\frac{1}{T_{\mathrm{ref},ij}}\Bigr),
\qquad z_iz_j\le0,\ \text{all-ion-pair dispersion, or } i,j\in C.
$$

## `cross_segment_diameter`

- Section: Parameter setup
- Status: Source equation
- Source: [Gross2001], Eq.~(A.14), with source-resulted l_ij extension
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
- Change note: The current reduced diameter y may select a cancellation-safe implementation, but is an algorithm choice, not a thermodynamic variable. With the switch at y = 0.2, chi and its first and second Taylor derivatives match 50-digit references to 3e-14, 7e-13 and 1.1e-11 relative on both branches (2026-09-23).

$$
\chi(y)=\frac{3[\ln(1+y)-y+y^2/2]}{y^3}\quad(y>0),\qquad \chi(0)=1,\qquad y=\kappa d_i^{DH},\qquad \chi(y)=\sum_{k=0}^{24}\frac{3(-y)^k}{k+3}+O(y^{25})\quad(y\to0).
$$

## `debye_huckel_helmholtz`

- Section: Electrolyte dielectric and Debye--Huckel relations
- Status: Source equation
- Source: [Bulow2019], Eq.~(2)
- Description: Finite-size Debye--Huckel residual Helmholtz contribution.
- Change note: d_i^DH is the separately declared finite-size diameter in metres; equality to packing or Born diameters requires source justification. The ion-free support restriction is identically zero before screening operations; the dilute three-halves power does not admit a finite ionic insertion Hessian.

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

## `dipolar_helmholtz`

- Section: Contribution intermediates / Dipole-dipole
- Status: Source equation
- Source: [Gross2006], Eq.~(7)
- Description: Optional third-order Padé dipole-dipole Helmholtz contribution; zero when every dipolar species is absent.
- Change note: Admitted as an opt-in family by ADR 0011.

$$
a^{dd}=\frac{A_2}{1-A_3/A_2},\qquad a^{dd}=0\ \text{if}\ A_2=0.
$$

## `dipolar_integrals`

- Section: Contribution intermediates / Dipole-dipole
- Status: Source equation
- Source: [Gross2006], Eqs.~(10)--(16) and Table 1
- Description: Pair and triplet correlation integrals with chain-length-interpolated coefficients, shared by the dipole (Table 1) and the quadrupole ([Gross2005], Eqs.~(12)--(19) and Table 2).
- Change note: The unlike energy in J2 is the geometric mean without k_ij. Only the dipole caps its chain lengths at 2 (ADR 0011); the quadrupole uses the uncapped $m_{ij}$, $m_{ijk}$ (ADR 0014).

$$
J_{2,ij}=\sum_{n=0}^{4}\left(a_{n,ij}+b_{n,ij}\frac{\epsilon_{ij}}{kT}\right)\eta^n,\quad
J_{3,ijk}=\sum_{n=0}^{4}c_{n,ijk}\eta^n,\quad
a_{n}=a_{0n}+\frac{m-1}{m}a_{1n}+\frac{m-1}{m}\frac{m-2}{m}a_{2n},\quad
\epsilon_{ij}=\sqrt{\epsilon_i\epsilon_j},\\
\text{dipole: } m_{ij}=\min\!\left(\sqrt{m_im_j},2\right),\ m_{ijk}=\min\!\left((m_im_jm_k)^{1/3},2\right);\qquad
\text{quadrupole: } m_{ij}=\sqrt{m_im_j},\ m_{ijk}=(m_im_jm_k)^{1/3}.
$$

## `dipolar_second_order`

- Section: Contribution intermediates / Dipole-dipole
- Status: Exact derivation
- Source: [Gross2006], Eq.~(8) with $q_i=(\epsilon_i/kT)n_{\mu,i}\mu_i^{*2}$
- Description: Second-order dipole-dipole perturbation term.
- Change note: The reduced moment $q_i=n_{\mu,i}\mu_i^2/(4\pi\varepsilon_0kTm_i\sigma_i^3)$ absorbs $\epsilon_i$.

$$
A_2=-\pi\widehat{\rho}\sum_i\sum_j x_ix_j\frac{\sigma_i^3\sigma_j^3}{\sigma_{ij}^3}q_iq_jJ_{2,ij}.
$$

## `dipolar_third_order`

- Section: Contribution intermediates / Dipole-dipole
- Status: Exact derivation
- Source: [Gross2006], Eq.~(9) with the same $q_i$
- Description: Third-order dipole-dipole perturbation term.
- Change note: Uses $\sigma_{ij}=(\sigma_i+\sigma_j)/2$ for every pair.

$$
A_3=-\frac{4\pi^2}{3}\widehat{\rho}^{\,2}\sum_i\sum_j\sum_k x_ix_jx_k
\frac{\sigma_i^3\sigma_j^3\sigma_k^3}{\sigma_{ij}\sigma_{ik}\sigma_{jk}}q_iq_jq_kJ_{3,ijk}.
$$

## `dispersion_chain_coefficients`

- Section: Contribution intermediates / Dispersion
- Status: Source equation
- Source: [Gross2001], Eqs.~(A.18)--(A.19)
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
- Source: [Gross2001], Eq.~(A.11)
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
- Source: [Gross2001], Eq.~(A.12)
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
- Source: [Gross2001], Eq.~(A.16)
- Description: First PC-SAFT dispersion integral correlation.
- Change note: Evaluated through the shared polynomial owner.

$$
I_1(\eta,\overline m)=\sum_{k=0}^{6}a_k(\overline m)\eta^k.
$$

## `dispersion_integral_two`

- Section: Contribution intermediates / Dispersion
- Status: Source equation
- Source: [Gross2001], Eq.~(A.17)
- Description: Second PC-SAFT dispersion integral correlation.
- Change note: Evaluated through the shared polynomial owner.

$$
I_2(\eta,\overline m)=\sum_{k=0}^{6}b_k(\overline m)\eta^k.
$$

## `dispersion_polynomial`

- Section: Contribution intermediates / Dispersion
- Status: Exact derivation
- Source: [Gross2001], polynomial form in Eqs.~(A.16)--(A.17)
- Description: Horner-evaluated sixth-order dispersion polynomial.
- Change note: Documents the real arithmetic performed by polynomial.

$$
\mathcal{P}(\boldsymbol{c},\eta)=\sum_{k=0}^{6}c_k\eta^k.
$$

## `dispersion_second_moment`

- Section: Contribution intermediates / Dispersion
- Status: Source equation
- Source: [Gross2001], Eq.~(A.13)
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
- Source: [Perdomo2025] Section 3.1 (elimination of one charged species before Eq. 18; reduced gradient Eq. 22)
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

## `fixed_pressure_heat_capacity_parameter_derivative`

- Section: Thermodynamic properties
- Status: Exact derivation
- Source: Composition of `fixed_pressure_total_derivative`, `total_heat_capacity_cp` and `residual_heat_capacity_cp`
- Description: First parameter derivative of total molar isobaric heat capacity at fixed T, P and normalized composition on a regular declared density branch. The ideal-gas thermochemistry coefficients are fixed inputs independent of the residual-EOS parameters. Parameter regression calls the Properties pressure-chart action once per parameter direction; the common energy order is three, including the implicit association-site chain (ALG-ACTIONS).
- Change note: Qualified by the total-value, derivative-ladder and density-lift checks in `engine/tests/native/regression_j2_inventory_checks.cpp` and the retained `analyses/2026-heat-capacity-fitting/` evidence. No parameter Hessian or reacting-equilibrium heat-capacity derivative is claimed.

$$
\left.\frac{dc_p}{d\theta}\right|_{T,P,\mathbf{x}}
=\left.\frac{\partial c_p^{res}}{\partial\theta}\right|_{T,\rho,\mathbf{x}}
-\left.\frac{\partial c_p^{res}}{\partial\rho}\right|_{T,\theta,\mathbf{x}}
 \frac{P_\theta}{P_\rho},
\qquad \frac{dc_p^{id}}{d\theta}=0,\qquad P_\rho>0.
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
- Source: [Gross2001], Eq.~(A.4)
- Description: PC-SAFT hard-chain residual Helmholtz contribution.
- Change note: Mirrors the site-independent chain term.

$$
a^{hc}=\overline{m}\,a^{hs}
-\sum_i x_i(m_i-1)\ln g_{ii}^{hs}.
$$

## `hard_sphere_contact`

- Section: Contribution intermediates / Hard chain
- Status: Source equation
- Source: [Gross2001], Eq.~(A.7) at i=j
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
- Source: [Gross2001], Eq.~(A.6)
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
- Description: Exact density closure for the declared neutral ideal-gas base; also the density scale for residual-model searches.
- Change note: The ideal-gas base has one stable vapor root with pressure slope RT and no packing bound.

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
- Change note: Removes the Born-selector dependency and keeps packing, Debye--Huckel, and Born lengths independently declarable. A declared charged molecule (i in C) takes `neutral_segment_diameter` and declares no packing diameter (#641).

$$
d_i=d_i^{\mathrm{pack}},\qquad z_i\ne0,\ i\notin C,\qquad d_i^{\mathrm{pack}}>0.
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

## `isentropic_compressibility`

- Section: Thermodynamic properties
- Status: Exact derivation
- Source: Exact derivation from `pressure`, `residual_chemical_potential` and the T/V Hessian columns of State; #566 design
- Description: Isentropic compressibility at fixed composition.
- Change note: Computed only on request; availability of `isothermal_compressibility` and `heat_capacity_ratio`.

$$
\kappa_S=-\frac{1}{v}\left(\frac{\partial v}{\partial P}\right)_{s,\boldsymbol{x}}=\kappa_T\frac{c_v}{c_p}.
$$

## `isobaric_expansivity`

- Section: Thermodynamic properties
- Status: Exact derivation
- Source: Exact derivation from `pressure`, `residual_chemical_potential` and the T/V Hessian columns of State; #566 design
- Description: Isobaric expansivity at fixed composition.
- Change note: Computed only on a derived_properties request; unavailable unless (dP/drho)_T,x is finite and positive. A positive slope does not certify phase stability.

$$
\alpha_P=\frac{1}{v}\left(\frac{\partial v}{\partial T}\right)_{P,\boldsymbol{x}}=\frac{P_T}{\rho_m P_\rho},\qquad P_T=\frac{P}{T}-RT\Phi_{VT}.
$$

## `isothermal_compressibility`

- Section: Thermodynamic properties
- Status: Exact derivation
- Source: Exact derivation from `pressure`, `residual_chemical_potential` and the T/V Hessian columns of State; #566 design
- Description: Isothermal compressibility at fixed composition.
- Change note: Computed only on a derived_properties request; unavailable unless (dP/drho)_T,x is finite and positive. A positive slope does not certify phase stability.

$$
\kappa_T=-\frac{1}{v}\left(\frac{\partial v}{\partial P}\right)_{T,\boldsymbol{x}}=\frac{1}{\rho_m P_\rho},\qquad P_\rho=\left(\frac{\partial P}{\partial\rho_m}\right)_{T,\boldsymbol{x}}.
$$

## `joule_thomson_coefficient`

- Section: Thermodynamic properties
- Status: Exact derivation
- Source: Exact derivation from `pressure`, `residual_chemical_potential` and the T/V Hessian columns of State; #566 design
- Description: Joule-Thomson coefficient (dT/dP at constant molar enthalpy) in K/Pa; negative values are physical.
- Change note: Computed only on request; unavailable unless total c_p is available and positive. Numerator and denominator are both molar.

$$
\mu_{JT}=\left(\frac{\partial T}{\partial P}\right)_{h,\boldsymbol{x}}=\frac{Tv_T-v}{c_p},\qquad v_T=\left(\frac{\partial v}{\partial T}\right)_{P,\boldsymbol{x}}=-\frac{P_T}{P_v}.
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
- Source: [Pereira2012] Eq. 1 (without reactions) and Eq. 3
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
- Description: Convert the mean ionic coefficient to formula-unit molality; define total solute-species molality for a single solvent.
- Change note: Uses solvent molar mass in kilograms per mole. Each occupied non-solvent molecular or ionic species contributes once to q.

$$
\gamma_\pm^{(m)}=
\frac{\gamma_\pm^{(x)}}
{1+M_s q},\qquad q=m(\nu_++\nu_-)\ \text{for one binary salt}.
\qquad m_i=\frac{x_i}{x_s M_s},\quad q=\sum_{i\ne s}m_i,\quad D=M_s q.
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
- Source: [Gross2001], Eq.~(A.9); for i in C, [Held2014] Eqs.~(21)--(22) and the text after them, [Uyan2015] Eqs.~(10)--(11)
- Description: Temperature-dependent hard-chain diameter for neutral components and declared charged molecules.
- Change note: Reconciled to the current neutral diameter branch. #641: C is the set of ions declaring `fixed.charged_molecule_ion`. Held2014 and Uyan2015 obtain the fixed ion diameter 0.88 sigma from this equation with the like-ion dispersion energy set to zero; an ion that keeps its own like-pair well depth takes this equation with its own u_i (an inference from that derivation; the sources state only the zero-dispersion direction). Its u_i derivative is `packing_coefficient_parameter_derivative`'s neutral rule.

$$
d_i(T)=\sigma_i(T)\left[1-0.12\exp\left(-\frac{3(\epsilon_i/k_B)}{T}\right)\right],
\qquad z_i=0\ \text{or}\ i\in C.
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

## `pair_parameter_temperature_law`

- Section: Parameter setup
- Status: Implementation convention
- Source: [Held2014], Tables 1--2; [Uyan2015], Eqs.~(26)--(27); [Bulow2021MDEA], Eq.~(26)
- Description: One declared temperature law for k_ij, l_ij, explicit association energy over k_B or volume, or a combining-rule volume correction; the family covers the cited affine, quadratic and cubic laws and the existing reciprocal-temperature term.
- Change note: Coefficients are active slots; T_r and the term powers are fixed declarations (#567). Undeclared k_ij keeps `cross_dispersion_energy`'s constant and reciprocal-slope arithmetic.

$$
q(T)=\sum_{p=0}^{3}a_p\,(T-T_r)^p+s\left(\frac1T-\frac1{T_r}\right),\qquad T,\ T_r>0.
$$

## `partial_molar_enthalpy`

- Section: Thermodynamic properties
- Status: Exact derivation
- Source: Exact derivation from `pressure`, `residual_chemical_potential` and the T/V Hessian columns of State; #566 design
- Description: Partial molar enthalpy of component i at fixed T and P on the supplied ideal thermochemistry basis.
- Change note: Computed only on request; availability of `partial_molar_volume` plus Available thermochemistry for the whole record. Component h_i^o comes from the existing interval selection and correlation owner at the record's standard pressure P^o. d_i is the fixed-P, fixed-x temperature derivative of ln phi_i. The x-weighted sum equals the total molar enthalpy.

$$
\bar h_i=h_i^\circ-RT^2d_i,\qquad d_i=\Phi_{n_iT}+\Phi_{n_iV}v_T+\frac{1}{T}-\frac{v_T}{v}.
$$

## `partial_molar_entropy`

- Section: Thermodynamic properties
- Status: Exact derivation
- Source: Exact derivation from `pressure`, `residual_chemical_potential` and the T/V Hessian columns of State; #566 design
- Description: Partial molar entropy of component i at fixed T and P on the supplied ideal thermochemistry basis.
- Change note: Availability of `partial_molar_enthalpy`, and unavailable at x_i = 0, where ln x_i diverges. It keeps the same-density residual entropy of State; the x-weighted sum equals the total molar entropy.

$$
\bar s_i=s_i^\circ-RL_i-RTd_i,\qquad L_i=\ln\frac{x_i\rho_mRT}{P^\circ}+\Phi_{n_i}.
$$

## `partial_molar_volume`

- Section: Thermodynamic properties
- Status: Exact derivation
- Source: Exact derivation from `pressure`, `residual_chemical_potential` and the T/V Hessian columns of State; #566 design
- Description: Partial molar volume of component i at fixed T and P, the amount derivative of extensive volume, evaluated at N = 1 mol.
- Change note: Computed only on request, in component order; unavailable for an IonicInsertion component or a nonpositive slope. A zero-fraction neutral species gives its finite insertion limit. For present ions only electroneutral combinations carry a physical claim.

$$
\bar v_i=-\frac{P_{n_i}}{P_V},\qquad P_{n_i}=RT\left(\frac{1}{V}-\Phi_{Vn_i}\right),\qquad P_V=-RT\left(\frac{1}{V^2}+\Phi_{VV}\right).
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
- Source: [HamerWu1972], practical molal osmotic-coefficient definition; common-pressure identity from Eqs.~`log_fugacity_coefficient` and `pressure`
- Description: Dimensionless practical molal osmotic coefficient from pure-solvent-referenced solvent activity.
- Change note: One neutral solvent, evaluated at common prescribed T/P with the same parameters and permittivity model. The scalar residual chemical-potential difference uses the extensive Helmholtz Euler identity in extended precision, with Z=P/(RT rho). Subscripts f and p denote the finite solution and pure solvent, a^r is per mole, and the sum traverses occupied non-solvent species. For one salt, q=(nu_++nu_-)m. At q=0 activity is exactly one and Phi is undefined; below the total-species molality floor Phi is unavailable. The one-salt helper also retains its formula-molality floor. Liquid-side regularity is required; this does not establish global phase stability.

$$
\Phi=
-\frac{\ln a_s}
{M_s q},
\qquad
\ln a_s=\ln x_s+\ln\phi_s-\ln\phi_s^{\mathrm{pure}}
=\ln x_s+\frac{\mu_s^r-\mu_s^{r,\mathrm{pure}}}{RT}
+\ln\frac{\rho_m}{\rho_m^{\mathrm{pure}}}.
\qquad
\Delta\mu_s^r/(RT)=
\frac{a_f^r-a_p^r+\frac{P}{RT}(\rho_f^{-1}-\rho_p^{-1})
-\sum_{i\ne s}x_i[\mu_i^{r,f}/(RT)-\mu_s^{r,p}/(RT)]}{x_s}.
$$

## `pressure`

- Section: Thermodynamic properties
- Status: Source equation
- Source: [Gross2001], Eq.~(A.23)
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
- Source: [Pereira2012] Eq. 8; [Gross2001] Eq. (A.23)
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

## `quadrupolar_helmholtz`

- Section: Contribution intermediates / Quadrupole-quadrupole
- Status: Source equation
- Source: [Gross2005], Eq.~(8)
- Description: Optional third-order Padé quadrupole-quadrupole Helmholtz contribution; zero when every quadrupolar species is absent.
- Change note: Admitted as an opt-in family by ADR 0014; it shares the dipole's Padé form and J integrals.

$$
a^{qq}=\frac{A_2}{1-A_3/A_2},\qquad a^{qq}=0\ \text{if}\ A_2=0.
$$

## `quadrupolar_second_order`

- Section: Contribution intermediates / Quadrupole-quadrupole
- Status: Exact derivation
- Source: [Gross2005], Eq.~(9) with $q_i=(\epsilon_i/kT)n_{Q,i}Q_i^{*2}$
- Description: Second-order quadrupole-quadrupole perturbation term.
- Change note: The reduced moment $q_i=n_{Q,i}Q_i^2/(4\pi\varepsilon_0kTm_i\sigma_i^5)$ absorbs $\epsilon_i$.

$$
A_2=-\frac{9\pi}{16}\widehat{\rho}\sum_i\sum_j x_ix_j\frac{\sigma_i^5\sigma_j^5}{\sigma_{ij}^7}q_iq_jJ_{2,ij}.
$$

## `quadrupolar_third_order`

- Section: Contribution intermediates / Quadrupole-quadrupole
- Status: Exact derivation
- Source: [Gross2005], Eq.~(10) with the same $q_i$ and $J_{3,ij}=0$ (Eq.~(14))
- Description: Third-order quadrupole-quadrupole perturbation term.
- Change note: $(4\pi^2/3)(3/4)^3=9\pi^2/16$; the sign is positive, as in Eq.~(A11).

$$
A_3=\frac{9\pi^2}{16}\widehat{\rho}^{\,2}\sum_i\sum_j\sum_k x_ix_jx_k
\frac{\sigma_i^5\sigma_j^5\sigma_k^5}{\sigma_{ij}^3\sigma_{ik}^3\sigma_{jk}^3}q_iq_jq_kJ_{3,ijk}.
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

## `ideal_mole_fraction_reaction`

- Section: Equilibrium / Reaction standard states
- Status: Source equation
- Source: Nasrifar and Tafazzol (2010), DOI 10.1021/ie901181n, Eq. (11), Table 1 and section 6, pp.7621,7627; ideal activity assumption
- Description: Dimensionless mole-fraction activities with standard fraction one; positive products stoichiometry and caller-declared constants on the ideal-mole-fraction basis. Pressure does not enter these liquid chemistry rows.
- Change note: Specializes reaction_activity_quotient. ReactionLogPolynomial owns the temperature law; its finite positive logarithm reference is independent of its evaluation interval. No EOS density or infinite-dilution correction applies.

$$
a_i=x_i=\frac{n_i}{N},\quad N=\sum_i n_i,\qquad
F_r=\sum_i\nu_{ir}\ln x_i-\ln K_r^x(T)=0.
$$

## `ideal_molarity_reaction`

- Section: Equilibrium / Reaction standard states
- Status: Exact derivation
- Source: IAPWS TGD9-18, Table 2, p.15, 25 °C concentration constants; reaction_activity_quotient and ideal_reaction_amount_coordinates
- Description: Fixed declared solution volume V in m³ and absolute amounts in mol, with standard concentration c° = 1000 mol/m³. A declared neutral solvent has unit activity; no identifier selects this exception.
- Change note: Products-positive stoichiometry. The formal proton convention H3O+ = H+ + H2O with unit water activity preserves the numerical source constants. This is ideal concentration chemistry, without EOS density closure or real-solution validation.

$$
a_i^c=\frac{n_i}{V c^\circ},\quad a_s^c=1\text{ for a declared unit-activity solvent},\qquad
F_r=\sum_i\nu_{ir}\ln a_i^c-\ln K_r^c(T)=0,\quad c^\circ=1000\ \mathrm{mol\,m^{-3}}.
\qquad m_s=0,\quad m_i=1\ (i\ne s),\qquad
\ln a_i^c=m_i[y_i+\ln(N_{\rm ref}/\mathrm{mol})-\ln(Vc^\circ/\mathrm{mol})],\qquad
\partial_j F_r=\nu_{jr}m_j,\quad\partial_{jk}F_r=0.
$$

## `ideal_reaction_amount_coordinates`

- Section: Equilibrium / Amount coordinates
- Status: Exact derivation
- Source: ideal_mole_fraction_reaction, ideal_molarity_reaction and material_balance; logarithm and exponential chain rules
- Description: Log amounts on either ideal activity basis, scaled by positive feed mol, with independent conserved null-space rows. There are s coordinates and r+(s-r) equations; charge conservation follows from the neutral feed and charge-conserving reactions.
- Change note: The feasibility objective is zero. No density, phase-amount, normalization or balance-dual coordinate is added. CppAD differentiates the same row expression for the existing Ipopt driver and exact-Jacobian polish.

$$
y_i=\ln(n_i/N_{\rm ref}),\quad N_{\rm ref}=N_{\rm feed},\quad
\ln x_i=y_i-\operatorname{logsumexp}_j y_j,\qquad
C_a=\sum_i B_{ai}e^{y_i}-b_a/N_{\rm ref},\quad B\nu=0,\quad b=Bn^0,
\qquad
\partial_j F_r=\nu_{jr}-\Delta\nu_r x_j,\quad
\partial_{jk}F_r=-\Delta\nu_r(\delta_{jk}x_j-x_jx_k),\quad
\partial_j C_a=B_{aj}e^{y_j},\quad\partial_{jk}C_a=\delta_{jk}B_{aj}e^{y_j}.
$$

## `mdea_reaction_conserved_amounts`

- Section: Equilibrium / Reaction inventories
- Status: Exact derivation
- Source: Null space of Nasrifar and Tafazzol (2010), Eqs. (1)-(6), p.7621; hydronium stoichiometry
- Description: Physical inventory checks for the caller-declared aqueous MDEA/CO2/H2S example. The Engine computes a general null-space basis, independent of these example species names.
- Change note: Water equivalent includes hydronium. The paper's printed Eq. (12) omits it and does not conserve R1; these conserved equations do not reproduce that printed balance.

$$
W=n_w+n_{OH^-}+n_{H_3O^+}+n_{HCO_3^-}+n_{CO_3^{2-}},\quad
M=n_{MDEA}+n_{MDEAH^+},\quad C=n_{CO_2}+n_{HCO_3^-}+n_{CO_3^{2-}},\quad
S=n_{H_2S}+n_{HS^-}+n_{S^{2-}},\quad
Z=n_{H_3O^+}+n_{MDEAH^+}-n_{OH^-}-n_{HCO_3^-}-2n_{CO_3^{2-}}-n_{HS^-}-2n_{S^{2-}}.
$$

## `effective_water_amounts`

- Section: Equilibrium / Amount conversion
- Status: Source equation
- Source: Nasrifar and Tafazzol (2010), section 6, p.7627, effective-water passage
- Description: Lump every charged species amount into declared neutral water, preserving every remaining free molecular species and the current true-species mole total. This molecular representation does not preserve every elemental inventory.
- Change note: Native conversion requires this ideal model's criterion-meeting liquid result and a complete disjoint water/ion/molecular partition. It does not calculate molecular phase equilibrium or reproduce Figure 15.

$$
\widetilde n_w=n_w+\sum_{i\in I}n_i,\qquad
\widetilde n_j=n_j\ (j\in J),\qquad
\widetilde N=\widetilde n_w+\sum_{j\in J}\widetilde n_j=N.
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
- Source: Austgen, Rochelle, Peng and Chen 1989, Ind. Eng. Chem. Res. 28, 1060, ln K(T); Parkhurst and Appelo (2013), USGS Techniques and Methods 6-A43, SOLUTION_SPECIES, printed p.200, analytical expression (https://doi.org/10.3133/tm6A43).
- Description: Temperature correlation of a reaction equilibrium constant on its declared standard state, with its analytic temperature derivative, which is the van 't Hoff reaction enthalpy on that standard state divided by RT²; valid on the declared temperature domain only.
- Change note: $T_{ref}>0$ is explicit. A legacy $c\ln(T/(1\,\mathrm K))$ coefficient is normalized without changing values by replacing $a_{old}$ with $a=a_{old}+c\ln(T_{ref}/(1\,\mathrm K))$. Standard-state quantities are held fixed in the displayed derivative; temperature-dependent reference conversions contribute separately. Engine fitting admits stored natural-log coefficients through `reaction_constant_coordinate_derivative`, at first order for state actions and for referenced finite-addition enthalpy actions at fixed temperature. Heat-capacity coefficient actions, coefficient Hessians and referenced free-temperature fitting are not implied; parameter adoption remains application-owned. Coefficients $(a,b,c,d,e,f)$ have units $(1,\mathrm K,1,\mathrm K^{-1},\mathrm K^2,\mathrm K^{-2})$; temperatures are in kelvin and $K$ is dimensionless on its declared source activity standard. All coefficients are finite and $0<T_{min}\le T\le T_{max}$; the interval is inclusive. The two added coefficients default to zero and remain fixed during fitting. PHREEQC's $\log_{10}K=A_1+A_2T+A_3/T+A_4\log_{10}[T/(1\,\mathrm K)]+A_5/T^2+A_6T^2$ maps, for the same reaction direction and source standard, by the displayed coefficient conversion.

$$
\begin{gathered}
L(T)\equiv\ln K=a+\frac bT+c\ln\frac{T}{T_{ref}}+dT+\frac e{T^2}+fT^2,\qquad \Delta_r h^{\circ}=RT^2 L'(T),\\
L'=-\frac b{T^2}+\frac cT+d-\frac{2e}{T^3}+2fT,\qquad L''=\frac{2b}{T^3}-\frac c{T^2}+\frac{6e}{T^4}+2f,\qquad L'''=-\frac{6b}{T^4}+\frac{2c}{T^3}-\frac{24e}{T^5},\\
a=\ln(10)A_1+A_4\ln[T_{ref}/(1\,\mathrm K)],\quad b=\ln(10)A_3,\quad c=A_4,\quad d=\ln(10)A_2,\quad e=\ln(10)A_5,\quad f=\ln(10)A_6.
\end{gathered}
$$

## `reaction_constant_coordinate_derivative`

- Section: Equilibrium / Reaction standard states
- Status: Exact derivation
- Source: Chain rule for `reaction_constant_ln_polynomial`, the minimum-norm standard state in `chemical_potential_stationarity`, and `implicit_sensitivity`; issue #187 accepted design
- Description: A fitted source-basis coefficient changes its density-basis reaction target by the same basis function. The source molality scale and reference contraction have no reaction-coefficient dependence. The minimum-norm image supplies chemical-row derivatives, and the existing defining Jacobian supplies the solved-state derivative.
- Change note: Coefficients a, b, c, d have units 1, K, 1, K⁻¹, with explicit positive Tref. Scalar ln K admits only a with basis 1. Evaluate the coefficient law at the path temperature T(t) for state actions, including its free-temperature chain. For referenced enthalpy actions at fixed temperature, coefficient directions also change the source temperature slope by $\phi_k'(T)$ and the reaction enthalpy by $RT^2\phi_k'(T)$. The bounded caloric actions consume the same analytic thermal coefficient terms; outer coefficient Hessians remain outside this consumer set. Source/reference declarations and non-fitted coefficients remain fixed. The stoichiometric matrix N has reactions as rows and species as columns. The constant source reaction-enthalpy shift $h$ (J/mol) at a fixed declared pivot $T_p$ is one derived coefficient direction (#484), preserving the supplied correlation. The existing coefficient chain consumes the combined direction $\mathbf w$. For the `reaction_constant_log10_reciprocal` view, a native law with $c=e=f=0$ replays as $(A,B,C)=-(b,a,d)/\ln 10$; source metadata remains caller-owned.

$$
\phi_k(T)=\left(1,\frac1T,\ln\frac{T}{T_{ref}},T\right)_k,
\qquad
\frac{\partial\ln K_{\rho,r}}{\partial\theta_{r,k}}=\phi_k(T),
\qquad
\frac{d}{dT}\frac{\partial\ln K_{src,r}}{\partial\theta_{r,k}}=\phi'_k(T)=\left(0,-\frac1{T^2},\frac1T,1\right)_k,
\qquad
\frac{\partial\Delta_r h_\rho}{\partial\theta_{r,k}}=RT^2\phi'_k(T)=(0,-R,RT,RT^2)_k,
\qquad
\frac{\partial s_i}{\partial\theta_{r,k}}=-\big(N^+\big)_{ir}\phi_k(T),
\qquad
\frac{d\mathbf y}{d\theta_{r,k}}=-\mathbf F_{\mathbf y}^{-1}\mathbf F_{\theta_{r,k}},
\quad \mathbf F_{\theta_{r,k}}=\text{chemical rows only},
\qquad
\frac{\partial(a,b,c,d)}{\partial h}=\mathbf w=
\left(\frac{1}{RT_p},-\frac1R,0,0\right),\qquad
\delta\ln K=\frac hR\left(\frac1{T_p}-\frac1T\right),\qquad
\frac{\partial\Delta_r h^\circ}{\partial h}=1,\qquad
\delta\frac{d\ln K}{d(1/T)}=-\frac hR.
$$

## `source_log_equilibrium_constant_observation`

- Section: Parameter regression / Reaction observations
- Status: Exact derivation
- Source: `reaction_constant_ln_polynomial` (Austgen et al. 1989, ln K(T) correlation); coefficient differentiation and existing Difference residual, issue #211
- Description: A direct observation is the named reaction's declared source-basis natural-log equilibrium constant at fixed positive temperature. It does not contain the EOS reference contraction or an equilibrium-state derivative.
- Change note: Fixed source declarations and Tref; a/b/c/d units are 1, K, 1, K⁻¹. Scalar K admits ln K with the a basis only. Reuse named-coordinate tying and rejection of fixed thermochemistry, affinity targets and dependent sets. Only Difference residuals are admitted. An absolute ln K uncertainty σ uses scale=1 and weight=σ⁻²; a chosen nondimensionalization scale requires the corresponding weight transformation. Covariance retains ALG-FIT's local rank/bound/error assumptions.

$$
y_r(T)=\ln K_{src,r}(T)=a_r+\frac{b_r}{T}+c_r\ln\frac{T}{T_{ref,r}}+d_rT,
\qquad
\frac{\partial y_r}{\partial\theta_{s,k}}=\delta_{rs}\left(1,\frac1T,\ln\frac{T}{T_{ref,r}},T\right)_k,
\qquad
\frac{\partial y_r}{\partial\theta_{EOS}}=0,
\qquad r_i=\frac{y_i-y_{obs,i}}{scale_i}.
$$

## `reaction_constant_log10_reciprocal`

- Section: Equilibrium / Reaction standard states
- Status: Bounded empirical relation
- Source: The −log10 K = a/T + b + cT form as recorded for MEA reaction R5 in the MEA-Thermodynamics source specification (Aroua 1999, Tong 2012)
- Description: Temperature correlation of a reaction equilibrium constant in decadic-logarithm form, converted to the natural logarithm, with its analytic temperature derivative; fitting uses the explicitly converted native coefficients; valid on the declared temperature domain only.
- Change note: The three-coefficient view requires native $c=e=f=0$; otherwise it is refused. Source notation ln(T) means ln(T/(1 K)); coefficients use the numerical Kelvin convention. Standard-state quantities are held fixed in the displayed derivative; temperature-dependent reference conversions contribute separately. Engine fitting admits stored natural-log coefficients through `reaction_constant_coordinate_derivative`, at first order only; parameter adoption remains application-owned.

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
- Description: Charge-neutral reaction contraction with a declared neutral solvent and optional per-species pure-liquid, pure-vapor or infinite-dilution references. Products have positive stoichiometry; solvent molar mass is in kg/mol and the molality standard in mol/kg. The source basis is explicit.
- Change note: The solvent standard uses the exact ion-free liquid state. Pure overrides require neutral species and a positive-slope root on the declared side of the pure isotherm's mechanical loop; a single root with no sampled loop serves either branch. Infinite-dilution species use `infinite_dilution_reference` at the same source pressure. Every returned T/P/parameter column differentiates that reference, with independent numerical availability. Nondefault mixed common-molality declarations are refused. Uyan (2015), Eqs. 22–23, specify pure water and solute infinite dilution; the historical CO2 vapor selection at e12a5ff8 is not an unambiguous paper attribution (#274). Per-species activities use raw mole fractions; uniform common molality retains A_r. Declarations are fixed under system composition actions. Fixed-source system-pressure derivatives are exact zeros. Ionic references qualify charge-neutral contractions, without an absolute single-ion standard.

$$
\Lambda_i^\rho=\ln(\rho/\rho_0)+\mu_i^{res}/(RT)
=\ln\phi_i+\ln(P/(RT\rho_0)),\quad
L_i^{b_i}=\begin{cases}
\Lambda_i^\rho(T,P_{src},\mathbf e_i;\mathrm{liquid}),&b_i=\mathrm{pure\ liquid},\\
\Lambda_i^\rho(T,P_{src},\mathbf e_i;\mathrm{vapor}),&b_i=\mathrm{pure\ gas},\\
\lim_{m\to0}\Lambda_i^\rho(T,P_{src},\mathbf x^s(m);\mathrm{liquid}),&b_i=\mathrm{infinite\ dilution},
\end{cases}\\
\mu_i^{\circ,b}=g_i^\rho+RT L_i^{b_i},\qquad
\ln a_i^b=\ln x_i+\Lambda_i^\rho-L_i^{b_i},\qquad
C_r^b=\sum_i\nu_{ri}L_i^{b_i},\quad \sum_i\nu_{ri}z_i=0,\\
\delta_i=L_i^{c_i}-L_i^{b_i},\qquad
\ln K_{c,r}=\ln K_{b,r}-\sum_i\nu_{ri}\delta_i,\qquad
\ln K_{c,r}+C_r^c=\ln K_{b,r}+C_r^b,\\
A_r=\left(\sum_{i\ne s}\nu_{ri}\right)\ln(M_s m^\circ),\qquad
\ln K_{\rho,r}=\begin{cases}
\ln K_{m,r}+A_r+C_r^b,&\text{common molality source},\\
\ln K_{x,r}+C_r^b,&\text{raw mole-fraction source},
\end{cases}\\
n_s=1/M_s,\quad n_i=m r_i\ (i\ne s),\quad \sum_i z_i r_i=0,\quad
\mathbf x^s(0)=\mathbf e_s,\\
\left.\partial_T L_i\right|_{P_{src}}=-\frac{h_i^{res}+RT}{RT^2},\qquad \partial_{P_{src}}L_i=\frac{v_i}{RT}\quad\text{(pure neutral species)},\\
L_i-\ln[P_{src}/(RT\rho_0)]\to0,\quad T L_{i,T}+1\to0,\quad P_{src}L_{i,P}-1\to0\quad(P_{src}\to0,\mathrm{pure\ gas}),\\
C_r^{(k)}=C_r(m_k),\quad m_k\in\{10^{-10},10^{-12},10^{-14}\}\ \text{(rung calculations)},\\
\widehat C_r=\frac{10\,C_r(10^{-14})-C_r(10^{-12})}{9},\qquad
\widehat C_r^{-}=\frac{10\,C_r(10^{-12})-C_r(10^{-10})}{9},\qquad
|\widehat C_r-\widehat C_r^{-}|\le 5\times10^{-5},\\
D=\partial_T^a\partial_{P_{src}}^b\partial_\theta^p C_r(T,P_{src};\mathbf x^s(0)).
$$

## `neutral_reaction_reference_pressure_jet`

- Section: Equilibrium / Reaction standard states
- Status: Exact derivation
- Source: Implicit differentiation of the pressure closure and the neutral reaction reference contraction at fixed temperature and composition
- Description: First and second source-pressure derivatives of the analytic solvent-limit contraction. The extensive regular residual potential is psi=Ares/(RT), at N=1 mol and V=1/rho_s; its single ionic insertion includes Born and the zero Debye–Hückel limit. DeltaNu is the stoichiometric mole change.
- Change note: Derivatives follow the same regular positive-slope density branch as the value. A reference at fixed source pressure has zero derivatives with respect to system pressure. E2 and E3 provide the displayed single-state partials; for mixed references, sum the masked dilution contribution and each selected pure-state contribution through `reference_single_insertion_taylor`. No runtime finite differences are used. Higher source-pressure derivatives are supplied by `reference_mixed_temperature_pressure_derivatives` in Properties; downstream use is separately qualified. The first temperature and parameter derivatives are `standard_state_transfer_temperature_derivative` and `standard_state_transfer_parameter_derivative`, and the mixed temperature-pressure derivative is `standard_state_transfer_temperature_pressure_derivative`.

$$
C=-\Delta\nu\ln(V\rho_0)+\sum_i\nu_i\psi_{n_i},\qquad
p_V=RT(-V^{-2}-\psi_{VV}),\quad p_{VV}=RT(2V^{-3}-\psi_{VVV}),\\
C_V=-\Delta\nu/V+\sum_i\nu_i\psi_{n_iV},\qquad
C_{VV}=\Delta\nu/V^2+\sum_i\nu_i\psi_{n_iVV},\\
C_P=\frac{C_V}{p_V},\qquad
C_{PP}=\frac{C_{VV}}{p_V^2}-\frac{C_Vp_{VV}}{p_V^3},\qquad
g_P=-RT C_P,\quad g_{PP}=-RT C_{PP}.
$$

## `speed_of_sound`

- Section: Thermodynamic properties
- Status: Exact derivation
- Source: Exact derivation from `pressure`, `residual_chemical_potential` and the T/V Hessian columns of State; #566 design
- Description: Thermodynamic speed of sound in m/s: the composition-frozen acoustic limit without transport, chemical relaxation or phase redistribution.
- Change note: Computed only on request; availability of `isentropic_compressibility`, plus a positive molar mass for every occupied species. M is in kg/mol; zero-fraction species do not contribute.

$$
w=\frac{1}{\sqrt{\rho_m M\kappa_S}}=\sqrt{\frac{c_p}{c_v}\frac{P_\rho}{M}},\qquad M=\sum_i x_iM_i.
$$

## `standard_state_transfer_temperature_derivative`

- Section: Equilibrium / Reaction standard states
- Status: Exact derivation
- Source: Implicit differentiation of the pressure closure at the fixed ion-free solvent composition of the neutral reaction reference; Lambda = ln phi + ln(P/(RT rho0)) and RT^2 d ln phi_i/dT|_P = -hbar_i^res
- Description: First temperature derivative at fixed reference pressure of the returned reference limit and of the density-basis equilibrium constant. psi=Ares/(RT) is evaluated at N=1 mol and V=1/rho, and the source basis term A is temperature independent. The target is transported as s=-ln K_rho, whose slope maps to species through the same minimum-norm owner as its value.
- Change note: Differentiate the returned solvent insertion with `reference_single_insertion_taylor`, retaining the implicit density and site chains. No finite dilution or terminal-change criterion is evaluated. A system-pressure row combines this slope with the pressure jet. The second T derivative and active-parameter derivatives have their own requested columns.

$$
p_T=R(V^{-1}-\psi_V)-RT\psi_{VT},\quad V_T=-p_T/p_V,\\
\left.\frac{dC}{dT}\right|_P=\sum_i\nu_i\psi_{n_iT}+C_VV_T,\qquad
RT^2\left.\frac{dC}{dT}\right|_P=-\sum_i\nu_i\left(\bar h_i^{res}+RT\right),\\
\frac{d\ln K_\rho}{dT}=\frac{d\ln K_{src}}{dT}+\left.\frac{dC}{dT}\right|_P,\qquad
\Delta h_\rho=\Delta h_{src}+RT^2\left.\frac{dC}{dT}\right|_P,\qquad
s_r'=-\frac{d\ln K_{\rho,r}}{dT}.
$$

## `standard_state_transfer_parameter_derivative`

- Section: Equilibrium / Reaction standard states
- Status: Exact derivation
- Source: Implicit differentiation of the pressure closure in an active EOS coordinate at fixed temperature, reference pressure and ion-free solvent composition of the neutral reaction reference
- Description: First derivative in an active EOS parameter theta, at fixed temperature and reference pressure, of the returned reference limit. psi=Ares/(RT) is evaluated at N=1 mol and V=1/rho. ln K_src and the molality scale A are parameter independent, so the target slope is -dC/dtheta. It maps to species through the same minimum-norm owner as the value and enters each chemical row linearly in theta at the solved state.
- Change note: One active coordinate at fixed T/P and ion-free solvent fractions. The parameter scale used in polarization is max(abs(theta),1) in native units, with max(abs(d_V),1e-3 m3/mol) for dielectric volume. Scaling is restored after extraction; it introduces no approximation or column admission tolerance. Referenced enthalpy additionally consumes `standard_state_transfer_temperature_parameter_derivative`; heat capacity consumes `reference_mixed_temperature_pressure_derivatives`.

$$
p_\theta=-RT\psi_{V\theta},\quad V_\theta=-p_\theta/p_V,\qquad
\left.\frac{dC}{d\theta}\right|_{T,P}=\sum_i\nu_i\psi_{n_i\theta}+C_VV_\theta,\\
\frac{\partial\ln K_\rho}{\partial\theta}=\left.\frac{dC}{d\theta}\right|_{T,P},\qquad
\frac{\partial s_r}{\partial\theta}=-\left.\frac{dC_r}{d\theta}\right|_{T,P}.
$$

## `standard_state_transfer_temperature_parameter_derivative`

- Section: Equilibrium / Reaction standard states
- Status: Exact derivation
- Source: Mixed implicit differentiation of the isobaric reference contraction and pressure closure along one active EOS coordinate, at fixed temperature, reference pressure and ion-free solvent composition
- Description: Mixed derivative $\partial^2C/\partial T\partial\theta$ of the returned reference limit, with $\psi=A^{res}/(RT)$ at $N=1$ mol. It supplies the EOS-coordinate derivative of referenced reaction enthalpy; $C$ is independent of the solved equilibrium state.
- Change note: Evaluate only for a consumed referenced-enthalpy parameter column. The shared single-insertion Taylor calculation supplies the mixed partial at the ion-free solvent root, with exact Born, association and density chains. A failed column retains its action cause and does not erase values or sibling directions. Higher consumed derivatives follow `reference_mixed_temperature_pressure_derivatives` and `reference_species_calorics`.

$$
\mathbf w=\mathbf e_T+V_T\mathbf e_V,\quad
\mathbf u=\mathbf e_\theta+V_\theta\mathbf e_V,\quad
V_\theta=-\frac{p_\theta}{p_V},\\
D^2C[\mathbf w,\mathbf u]=\sum_i\nu_iD^3\psi[\mathbf e_{n_i},\mathbf w,\mathbf u]
+\Delta\nu\frac{V_TV_\theta}{V^2},\\
D^2p[\mathbf w,\mathbf u]=-R\left(\psi_{V\theta}+V_\theta\psi_{VV}+\frac{V_\theta}{V^2}\right)
+\frac{2RTV_TV_\theta}{V^3}-RTD^3\psi[\mathbf e_V,\mathbf w,\mathbf u],\\
V_{T\theta}=-\frac{D^2p[\mathbf w,\mathbf u]}{p_V},\quad
\frac{\partial^2C}{\partial T\,\partial\theta}=D^2C[\mathbf w,\mathbf u]+C_VV_{T\theta},\quad
\frac{\partial\Delta_rh_\rho}{\partial\theta}=RT^2\frac{\partial^2C}{\partial T\,\partial\theta},\\
D^3\psi[\mathbf q,\alpha\mathbf w+\mathbf u,\alpha\mathbf w+\mathbf u]
-D^3\psi[\mathbf q,\alpha\mathbf w-\mathbf u,\alpha\mathbf w-\mathbf u]
=4\alpha D^3\psi[\mathbf q,\mathbf w,\mathbf u],\quad
\mathbf q\in\{\boldsymbol\nu,\mathbf e_V\},\quad
\alpha=\frac{\|\mathbf u\|_s}{\|\mathbf w\|_s}.
$$

## `standard_state_transfer_temperature_second_derivative`

- Section: Equilibrium / Reaction standard states
- Status: Exact derivation
- Source: Second implicit differentiation of the pressure closure at the fixed ion-free solvent composition; for identical species RT^2 dLambda/dT = -(h^res + RT) gives d2Lambda/dT2 = -(c_p^res + R)/(RT^2) + 2(h^res + RT)/(RT^3)
- Description: Second temperature derivative at fixed reference pressure of the returned reference limit along w = e_T + V_T e_V, with psi = Ares/(RT) at N = 1 mol. The shared single-insertion Taylor owner evaluates this lower partial, including implicit association and pressure closure. It supplies the EOS-basis reaction heat capacity consumed by reference species calorics.
- Change note: Evaluate only for a consumed heat-capacity value, enthalpy temperature action or second temperature action. The analytic solvent-limit column supplies s_TT=-(d2 ln K_src/dT2 + C_TT). Numerical failures retain their causes independently. Qualification retains the existing 0.01 J/(mol K) scaled comparison with the pinned finite calculation; this comparison is outside runtime.

$$
\mathbf{w}=\mathbf{e}_T+V_T\mathbf{e}_V,\qquad
D^2C[\mathbf{w},\mathbf{w}]=\sum_i\nu_iD^3\psi[\mathbf{e}_{n_i},\mathbf{w},\mathbf{w}]+\Delta\nu\,V_T^2/V^2,\\
D^2p[\mathbf{w},\mathbf{w}]=-\frac{2RV_T}{V^2}+\frac{2RTV_T^2}{V^3}
-R\left(2(\psi_{VT}+V_T\psi_{VV})+T\,D^3\psi[\mathbf{e}_V,\mathbf{w},\mathbf{w}]\right),\\
V_{TT}=-\frac{D^2p[\mathbf{w},\mathbf{w}]}{p_V},\qquad
\left.\frac{d^2C}{dT^2}\right|_P=D^2C[\mathbf{w},\mathbf{w}]+C_VV_{TT},\qquad
\Delta c_{p,\rho}=2RT\frac{d\ln K_\rho}{dT}+RT^2\frac{d^2\ln K_\rho}{dT^2}.
$$

## `standard_state_transfer_temperature_pressure_derivative`

- Section: Equilibrium / Reaction standard states
- Status: Exact derivation
- Source: Implicit differentiation of C_P = C_V/p_V (`neutral_reaction_reference_pressure_jet`) along the isobaric temperature direction w of `standard_state_transfer_temperature_second_derivative`, at the fixed ion-free solvent composition
- Description: Mixed temperature-pressure derivative of the returned reference limit of a system-pressure reference, with psi = Ares/(RT) at N = 1 mol. The shared single-insertion Taylor owner evaluates the raw source-pressure partial at the resolved terminal state. It supplies s_TP = -d2C/dTdP to the chemical rows of second-order actions that move temperature and pressure, and the pressure slope of the EOS-basis reaction enthalpy, d(Delta_r h_rho)/dP = RT^2 d2C/dTdP, to reference species calorics.
- Change note: Evaluate when pressure and temperature columns are both requested. A fixed source pressure gives exact zero system-pressure derivatives. Differentiate the same returned solvent insertion; a failed mixed column leaves values and sibling columns available. Higher orders are supplied only by `reference_mixed_temperature_pressure_derivatives`.

$$
\mathbf{w}=\mathbf{e}_T+V_T\mathbf{e}_V,\qquad
DC_V[\mathbf{w}]=\sum_i\nu_iD^3\psi[\mathbf{e}_{n_i},\mathbf{e}_V,\mathbf{w}]+\Delta\nu\,V_T/V^2,\\
Dp_V[\mathbf{w}]=R\left(-V^{-2}-\psi_{VV}\right)+\frac{2RTV_T}{V^3}-RT\,D^3\psi[\mathbf{e}_V,\mathbf{e}_V,\mathbf{w}],\\
\frac{\partial^2C}{\partial T\,\partial P}=\frac{1}{p_V}\left(DC_V[\mathbf{w}]-\frac{C_V\,Dp_V[\mathbf{w}]}{p_V}\right),\qquad
\frac{\partial^2s_r}{\partial T\,\partial P}=-\frac{\partial^2C_r}{\partial T\,\partial P},\qquad
\frac{\partial\Delta_rh_\rho}{\partial P}=RT^2\frac{\partial^2C}{\partial T\,\partial P}.
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

## `reference_single_insertion_taylor`

- Section: Helmholtz derivatives / Reference derivatives
- Status: Exact derivation
- Source: Taylor expansion of the neutral reaction contraction and implicit pressure closure; issue #225, section 3
- Description: One fixed-volume amount insertion in the extensive reduced residual Helmholtz energy supplies gradient coefficients along a pressure-closed, fixed-composition base path. Temperature and EOS coordinates can move. Orders one and two also supply the lower reference slopes. The degree-seven/nine evaluations are bounded reference calculations, not unrestricted energy derivative surfaces.
- Change note: The same central unweighted site Jacobian eliminates association through degree nine. At an ion-free reference, a separate tape reuses the regular scalar contributions and Born radial work while excluding the analytically zero DH insertion; ordinary state admission is unchanged. Two amount insertions start above the extracted degree. No zero-concentration ionic Hessian is evaluated. Coefficients retain factorials until mixed extraction, and physical normalization is restored once. Mask pure-selected species from the solvent insertion and add their pure-state insertions. Fixed-source system-pressure derivatives are zero.

$$
z(t)=z_0+\sum_{j=1}^{r}z_jt^j,\quad r\in\{0,1,2,3,4\},\qquad
[t^k]D\psi(z(t))[u]=\frac12[t^{r+1+k}]
\{\psi(z(t)+ut^{r+1})-\psi(z(t)-ut^{r+1})\},\quad 0\le k\le r,\\
D^4\psi[u,v,v,v]=\frac{3!}{2}[t^7]\{\psi(z_0+vt+ut^4)-\psi(z_0+vt-ut^4)\},\\
D^5\psi[u,v,v,v,v]=\frac{4!}{2}[t^9]\{\psi(z_0+vt+ut^5)-\psi(z_0+vt-ut^5)\},\\
p(t)=RT(t)\{V(t)^{-1}-D\psi(z(t))[e_V]\},\qquad
V_k=\frac{[t^k]P_{src}(t)-[t^k]p(z(t))|_{V_k=0}}{p_V(z_0)},\\
C(t)=-\Delta\nu\ln(V(t)\rho_0)+D\psi(z(t))[\nu],\qquad
z_{j,n}=0\quad(j>0).
$$

## `reference_mixed_temperature_pressure_derivatives`

- Section: Properties / Reaction references
- Status: Exact derivation
- Source: Differentiation of `neutral_reaction_reference_transfer` with the single-insertion construction; issue #225, sections 2–3 and F8/F9
- Description: Raw SI derivatives through fourth temperature/pressure order, and through third temperature/pressure order with one active EOS-coordinate partial. No repeated parameter derivatives or dilution-composition derivatives are consumed.
- Change note: Each column differentiates the returned analytic reference at the ion-free solvent state using `reference_single_insertion_taylor`. Columns retain independent action diagnostics and causes. The existing scaled comparison with `reference_limit_terminal_pair` is qualification evidence only; no terminal sequence is evaluated in runtime.

$$
C_{ab}=\partial_T^a\partial_{P_{src}}^b C,\quad a+b\le4,\qquad
C_{ab,\theta}=\partial_T^a\partial_{P_{src}}^b\partial_\theta C,\quad a+b\le3,\\
C_T=\frac{h^{res}+RT}{RT^2},\qquad C_P=-\frac{v}{RT}\quad\text{(identical neutral species, }C=-\Lambda\text{)}.
$$

## `reference_limit_terminal_pair`

- Section: Electrolyte reference states and reported activities
- Status: Implementation convention
- Source: Finite terminal-change diagnostic for the declared reference molality sequence.
- Description: Pinned qualification comparison for the former finite-molality reference. It extrapolates the 1e-12 and 1e-14 mol/kg points and reports its change from the 1e-10/1e-12 pair. It is not the runtime reference after #823.
- Change note: Convergence diagnostic; neither an error bound nor an estimate of actual bias without additional convergence assumptions. The neutral limit removes the Debye-Hueckel c sqrt(m) tail exactly; an O(m) remainder d m leaves -d 1e-13. For Held 2014 NaCl in water at 298.15 K, e_f is 1.3e-11 against raw rung changes of 1.06e-5.

$$
\Delta_k=f(m_k)-f(m_{k-1}),\qquad
\hat f_k=\frac{10f(m_k)-f(m_{k-1})}{9},\quad m_k=m_{k-1}/100,\quad e_f=s\,|\hat f_{10^{-14}}-\hat f_{10^{-12}}|.
$$

## `reference_molality_sequence`

- Section: Electrolyte reference states and reported activities
- Status: Implementation convention
- Source: EOS convergence sequence for the [Figiel2025] infinite-dilution reference
- Description: Historical finite-molality inputs retained only for qualification of `infinite_dilution_reference`. No runtime ladder remains.
- Change note: The paired baseline declaration and numerical record of #823 identify the old inputs and extrapolation. Terminal changes are convergence diagnostics, not certified error bounds. The runtime molality is now exactly zero and terminal-change fields are unavailable; solvent-root and association diagnostics remain defined.

$$
\mathcal M_{\mathrm{ref}}=\{10^{-6},10^{-8},10^{-10},10^{-12},10^{-14}\}\,\mathrm{mol\,kg^{-1}},\qquad \widehat{\ln\phi_i^\infty}=\frac{10\ln\phi_i(10^{-14}\,\mathrm{mol\,kg^{-1}})-\ln\phi_i(10^{-12}\,\mathrm{mol\,kg^{-1}})}{9}.
$$

## `infinite_dilution_reference`

- Section: Electrolyte reference states and reported activities
- Status: Exact derivation
- Source: [Gross2001], A.32–A.33; [Chapman1990], Eqs. 21–24; `born_radial_charging_work`; #823 accepted design D1–D12
- Description: Fixed-volume first insertion into the extensive residual potential at one ion-free solvent pressure root, for a pure solvent or declared neutral solvent mixture. The reference is dimensionless and uses one total log Z correction.
- Change note: Hard chain, dispersion, association and admitted polar energy reuse their scalar owners. Association retains solvent-contact response even for a site-free ion. Born has a nonzero insertion; DH and cation pair/triple mixing have zero first insertions. Bulk permittivity and Born shell diameters must be positive, with a regular density root and association Jacobian. A net-neutral charged background is refused as IonicReferenceBackground. T/P and active-coordinate derivatives retain the implicit density/site chains. This defines no ionic-concentration Hessian. No factor of solvent molar mass enters b_i or ln(phi_i); the existing molality conversion remains separate. The retained NaCl comparisons have a case-specific 1e-9 log-fugacity criterion, not a universal error bound or physical validation.

$$
b_i=\left.\partial_{n_i}\left[\frac{A^r}{RT}\right]\right|_{T,V,n_{j\ne i},\mathbf n=\mathbf y_s},\qquad
b_i=a_s^r+\rho_s a^r_{\rho}+a^r_{x_i}-\sum_s y_s a^r_{x_s},\\
P(T,\rho_s,\mathbf y_s)=P_{src},\qquad
\ell_i^\infty=b_i-\ln\frac{P_{src}}{\rho_s RT},\qquad
L_i^\infty=b_i+\ln\frac{\rho_s}{\rho_0}.
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
- Source: [Gross2001], Eq.~(A.46)
- Description: Residual molar enthalpy.
- Change note: Includes the residual pressure-volume term Z-1.

$$
h^r=RT\left[-Ta_T+(Z-1)\right].
$$

## `residual_entropy`

- Section: Thermodynamic properties
- Status: Source equation
- Source: [Gross2001], Eq.~(A.48) without the ln Z term (same temperature and density)
- Description: Residual entropy at the same temperature and density.
- Change note: Does not silently change to a same-pressure ideal-gas reference.

$$
s_V^r=-R\left(a^r+Ta_T\right).
$$

## `residual_heat_capacity_cp`

- Section: Thermodynamic properties
- Status: Exact derivation
- Source: Temperature and density derivatives of [Gross2001] Eqs.~(A.22) and (A.46)
- Description: Residual isobaric molar heat capacity on a conditioned mechanical branch.
- Change note: The final minus one removes the ideal-gas cp-cv contribution.

$$
\frac{c_p^r}{R}
=\frac{c_v^r}{R}
+\frac{(Z+TZ_T)^2}{Z+\rho_mZ_\rho}-1.
$$

## `residual_heat_capacity_cv`

- Section: Thermodynamic properties
- Status: Exact derivation
- Source: Temperature derivative of [Gross2001] Eq.~(A.46)
- Description: Residual isochoric molar heat capacity.
- Change note: Transforms CppAD first and second temperature partials.

$$
c_v^r=R\left(-2Ta_T-T^2a_{TT}\right).
$$

## `cation_pair_triple_mixing_helmholtz`

- Section: Residual Helmholtz energy
- Status: Exact derivation
- Source: Solvent-mass integration of the constant short-range terms in Song et al. (2010), §2, p. 1287, Eqs. (1), (2), (10-1), (10-2), Zotero W3AX2GHK; #658 design D1–D7. Multi-anion indexing: Pitzer and Kim (1974), Eqs. (8), (9), (11), [doi:10.1021/ja00825a004](https://doi.org/10.1021/ja00825a004); like-sign θ independent of the counter-ion and ψ per ion triple as written in Peiper and Pitzer (1982), J. Chem. Thermodyn. 14, 613–638, Eqs. (7)–(8), Zotero X8GEDSJ2.
- Description: One declared neutral solvent s, distinct cations i and j and declared anion k; M_s is kg/mol, W is kg solvent, theta is kg/mol and psi is kg²/mol². The reduced molar contribution is Delta a^r = F/N, with N = sum n.
- Change note: Evaluated in every phase on the declared solvent's molality basis, with positive solvent mass and the declared anion retained. Additional anions are admitted: theta_ij is independent of the anion; only the declared psi_ijk is nonzero, with psi_ija = 0 for every other anion as a model assumption. The unequal-charge electrostatic Etheta terms are excluded. Transfer of chloride-fitted coefficients to sulfate is unvalidated. Other neutral components enter this contribution only through mole fractions and total-mole normalization; their masses are not added to W. The constant coefficients are qualified only for retained aqueous water/chloride-salt cases at 298.15 K and 1 bar. Occupied additional neutrals and organic-phase effects are unqualified extrapolation, with no phase selector or automatic deactivation; at exact zero declared solvent with either participating cation exactly absent, the contribution is defined as zero on that fixed support. All derivatives within that support are zero; transverse entries in dense arrays have defined zero storage and claim no physical insertion derivative. Dry states containing both cations remain refused, including zero coefficients. Existing absent-ion action refusals remain; in a dry phase containing one participating cation, the absent partner’s scalar chemical potential and fugacity are also unavailable with IonicInsertion. This is not a continuous full-composition extension or organic-ion qualification. No density dependence, pressure correction, unequal-charge electrostatic mixing, vapor, dry-salt, mixed-solvent or caloric qualification. For W > 0, faces retain nonzero absent-ion insertion derivatives. Coefficients correct this EOS and are not imported Pitzer coefficients. Full EOS §12 gives an illustrative organic-phase substitution and its limits.

$$
F=\frac{\Delta A^r}{RT}
=Wm_im_j\left(2\theta_{ij}+\sum_{a\in\mathcal A}\psi_{ija}m_a\right)
=\frac{2\theta_{ij}n_i n_j}{W}+\frac{\psi_{ijk}n_i n_j n_k}{W^2},\quad
\psi_{ija}=0\ (a\ne k),\quad
W=M_s n_s>0,\quad m_l=\frac{n_l}{W},\qquad
\Delta a^r=\frac{2\theta x_i x_j}{M_s x_s}
+\frac{\psi x_i x_j x_k}{(M_s x_s)^2},\quad
\Delta P=0,\qquad
\Delta\Phi=\frac{2m_im_j}{\sum_{l\ne s}m_l}(\theta+\psi m_k).
$$

## `residual_helmholtz`

- Section: Residual Helmholtz energy
- Status: Implementation convention
- Source: [Gross2001], PC-SAFT residual decomposition extended by admitted electrolyte terms
- Description: Sum of every active residual Helmholtz contribution; $a^{dd}$ only when the optional polar family is declared, and $a^{mix}$ only with `cation_mixing`.
- Change note: Mirrors the single EOS Helmholtz owner and its contribution result.

$$
a^r=a^{hc}+a^{disp}+a^{assoc}+a^{DH}+a^{Born}+a^{dd}+a^{mix}.
$$

## `residual_internal_energy`

- Section: Thermodynamic properties
- Status: Exact derivation
- Source: Temperature derivative of [Gross2001] Eq.~(A.46) with u=h-Pv
- Description: Residual molar internal energy.
- Change note: Uses the CppAD temperature partial at fixed density and composition.

$$
u^r=-RT^2a_T.
$$

## `salt_free_solvent_permittivity`

- Section: Electrolyte dielectric and Debye--Huckel relations
- Status: Source equation
- Source: [Figiel2025], Eq.~(12); [AscaniHeld2021], Eq.~(13)
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
- Source: [Figiel2025], same-sign ion exclusion; [Ji2012] §2.1 and Eq.~(4), all-ion-pair dispersion and ion-based species with dispersion among all ion pairs; [Held2014] text after Eqs.~(21)--(22)
- Description: Select the source-declared dispersion rule for equally charged ions.
- Change note: The caller selects all ion pairs or same-sign exclusion as part of the electrolyte model family. #641: C is the opt-in set of charged molecules (`fixed.charged_molecule_ion`, no default). The pair i, j in C includes i = j, so a declared ion keeps its self pair with k_ii = 0; two different declared molecules of the same sign keep their pair under the combining rule and need a declared k_ij like any applicable pair.

$$
\epsilon_{ij}=\begin{cases}\sqrt{\epsilon_i\epsilon_j}(1-k_{ij}),&z_iz_j\leq 0,\ \text{all-ion-pair dispersion, or } i,j\in C,\\0,&z_iz_j>0,\ \text{same-sign exclusion, and not } i,j\in C.\end{cases}
$$

## `segment_sigma_correlation`

- Section: Parameter setup
- Status: Bounded empirical relation
- Source: Source-resulted component correlation represented by the EOS parameter fields
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
- Description: System Gibbs objective in Helmholtz coordinates; density-standard offsets appear exactly once and are consistent with ideal_gas_helmholtz. The ideal-gas base uses q = ln[(1/v)/rho0], v = exp(-q)/rho0 for unit-amount volume v; at an unnormalized trial S = sum_i x_i, physical density is S/v.
- Change note: Use either mechanical energy plus offsets or total energy, never both reference terms. At fixed q, v_q = -v, v_qq = v and composition/temperature partials vanish. Do not normalize trial compositions.

$$
\frac{G}{RT}=\sum_p\left[\frac{A_{p,\mathrm{mech}}}{RT}+\frac{PV_p}{RT}+\sum_i n_{p,i}\frac{g^{\circ}_{i,\rho}}{RT}\right].
$$

## `tangent_plane_distance`

- Section: Equilibrium / Tangent-plane testing
- Status: Exact derivation
- Source: [Michelsen1982] Part I; [Pereira2012] Eq. 12; [Perdomo2025] §3.3 Eqs. 37–40; [Ascani2022] §3.3 Eqs. 33–35; conserved-plane extension adopted in #223
- Description: Reduced molar tangent-plane distance at fixed positive T/P, normalized charge-neutral composition w on support S and molar volume v in m³/mol; e is the all-ones vector and z contains charge numbers. The mechanical reduced energy includes residual and ideal-density Helmholtz energy plus Pv/RT; s is the dimensionless density-standard offset and π the chemical balance plane.
- Change note: Freeze standards including reference transfer and gauge at the solved snapshot. Chemical λ is dimensionless here, unlike dimensional chemical_potential_stationarity. A negative plane value needs feasible donor withdrawal and error qualification; finite search supplies no global certificate.

$$
\begin{aligned}
\widetilde g(\mathbf w,v)&=\frac{a_{\mathrm{mech}}(\mathbf w,v)}{RT}+\frac{Pv}{RT},\qquad
s_i=\frac{g^\circ_{i,\rho}(T,P)}{RT},\qquad \boldsymbol\pi=\mathbf B^{\mathsf T}\boldsymbol\lambda,\\
t(\mathbf w,v)&=\widetilde g(\mathbf w,v)+\mathbf w^{\mathsf T}(\mathbf s-\boldsymbol\pi),\qquad
\mathbf e^{\mathsf T}\mathbf w=1,\quad \mathbf z^{\mathsf T}\mathbf w=0,\quad
w_i>0\ (i\in S),\ w_i=0\ (i\notin S).
\end{aligned}
$$

## `tangent_plane_gibbs_dual_identity`

- Section: Equilibrium / Tangent-plane testing
- Status: Exact derivation
- Source: Homogeneity of system_gibbs_helmholtz_coordinates and chemical_potential_stationarity; #223 adopted derivation
- Description: Gibbs-dual identity for an exact feasible all-free fixed-T/P stationary candidate and any charge-neutral split with the same conserved totals b. N_q is mol; extensive reduced energies are on the request feed basis.
- Change note: Pinned feasibility cells and approximate numerical states do not inherit exact stationary equality. The split identity alone does not construct a feasible withdrawal.

$$
\frac{G_0}{RT}=\boldsymbol\lambda^{\mathsf T}\mathbf b,\qquad
\frac{G_{\mathrm{split}}}{RT}-\boldsymbol\lambda^{\mathsf T}\mathbf b
=\sum_q N_q\,t(\mathbf x_q,v_q),\qquad
\sum_q N_q\mathbf B\mathbf x_q=\mathbf b.
$$

## `tangent_plane_trial_balance`

- Section: Equilibrium / Tangent-plane testing
- Status: Exact derivation
- Source: Basis transformation of normalization/electroneutrality and Euler stationarity; #223 adopted derivation
- Description: Internal unit-amount trial balance basis on S; R contains total-amount and charge rows, omitting charge for neutral support. H is the retained full-rank row transformation; λ_t is conjugate to B_t, not an optimizer multiplier.
- Change note: No feed-moiety-proportion constraint or reaction is imposed. Charged density-temperature Vapor trials use H=I and the numerical log-charge row h below, evaluated by maximum-shifted log sums. Only that trial charge balance or a declared Vapor Electro row is replaced in Newton correction and classification. Original physical criteria, Ipopt callbacks and solved-state derivatives remain unchanged; full transformed rank does not establish branch or recovery accuracy. Qualification remains incomplete for #771. Other trials retain their existing basis; use λ_t dot b_t in every basis.

$$
\begin{aligned}
\mathbf R_S&=\begin{pmatrix}\mathbf e_S^{\mathsf T}\\\mathbf z_S^{\mathsf T}\end{pmatrix},\quad
\mathbf B_t=\mathbf H\mathbf R_S,\quad
\mathbf b_t=\mathbf B_t\mathbf w^0=\mathbf H\mathbf R_S\mathbf w^0,\qquad
\mathbf R_S\mathbf w^0=\begin{pmatrix}1\\0\end{pmatrix}
\quad\text{(omit charge row and second entry for neutral support)},\\
N\mathbf B_t\mathbf w&=\mathbf b_t,\quad \mathbf e_S^{\mathsf T}\mathbf w=1
\ \Longrightarrow\ N=1,\ \mathbf z_S^{\mathsf T}\mathbf w=0,\\
\frac{\boldsymbol\mu_{\mathrm{mech}}}{RT}+\mathbf s-\boldsymbol\pi
&=\mathbf B_t^{\mathsf T}\boldsymbol\lambda_t
=k\mathbf e_S+\varphi_t\mathbf z_S,\qquad
 t=\boldsymbol\lambda_t^{\mathsf T}\mathbf b_t=k
\quad\text{at exact feasible stationary density and composition}.
\end{aligned}
\qquad
\begin{aligned}
Q_+&=\sum_{z_i>0}z_i e^{y_i},\quad Q_-=\sum_{z_i<0}|z_i|e^{y_i},\\
h&=\ln Q_+-\ln Q_-,\quad R=|\tanh(h/2)|,\\
\frac{\partial h}{\partial y_i}&=
\begin{cases}z_i e^{y_i}/Q_+,&z_i>0,\\
z_i e^{y_i}/Q_-,&z_i<0,\\0,&z_i=0.\end{cases}
\end{aligned}

$$

## `michelsen_modified_tangent_plane`

- Section: Equilibrium / Tangent-plane testing
- Status: Source equation
- Source: [Michelsen1982] Part I Eqs. 13–19 and 25–27; T,P Jacobian derived for #647 from the IsothermalPhaseHessian columns
- Description: Michelsen's modified tangent-plane objective g* over positive trial mole numbers W on support S, at fixed T/P and the frozen shift s = standard − plane of `tangent_plane_distance`; w = W/ΣW and v is the trial label's density-owner molar volume (m³/mol). ψᵢ = Fᵢ − ln v = ln φᵢ + ln(P/RT), with F the residual reduced Helmholtz energy at unit amount (Fᵢ dimensionless, F_iV in m⁻³, F_VV in mol m⁻⁶, and Fᵢⱼ in mol⁻¹, using SI derivative columns). Jᵢⱼ = n ∂ln φᵢ/∂nⱼ at fixed T,P is dimensionless and symmetric with Σᵢ wᵢJᵢⱼ = 0.
- Change note: All-neutral, non-reacting systems only (#647). A stationary point has tpd = k = −ln ΣW; the shared tail reevaluates tpd and keeps every acceptance check.

$$
\begin{aligned}
g^*(\mathbf W)&=1+\sum_{i\in S}W_i\left(\ln W_i+\psi_i(\mathbf w)+s_i-1\right),\qquad
\frac{\partial g^*}{\partial W_i}=g_i=\ln W_i+\psi_i+s_i,\\
\ln W_i^{(k+1)}&=-s_i-\psi_i^{(k)}\quad\text{(Eq. 19)},\qquad
\alpha_i=2\sqrt{W_i},\qquad
H_{ij}=\delta_{ij}+\frac{\sqrt{W_iW_j}}{\sum_l W_l}J_{ij}+\tfrac12\delta_{ij}g_i\quad\text{(Eqs. 25–27)},\\
J_{ij}&=F_{ij}+1+\frac{p_ip_j}{p_V},\qquad p_i=\frac1v-F_{iV},\qquad p_V=-F_{VV}-\frac1{v^2}.
\end{aligned}
$$

## `tangent_plane_shift_log_curvature`

- Section: Equilibrium / Tangent-plane testing
- Status: Exact derivation
- Source: Exact log-composition chain rule for the frozen standard shift in system_gibbs_helmholtz_coordinates; #223
- Description: Shift coefficients c=s−π are fixed, but supported fractions w_i=exp(y_i) and free amount N are coordinates. Curvature of the extensive reduced shift must remain in shared exact callbacks at arbitrary points.
- Change note: Linear amount-space energy does not imply zero log-chart curvature; direct amount curvature is zero.

$$
L_c=N\sum_{i\in S}c_i e^{y_i},\qquad
\partial_{y_i}L_c=N c_iw_i,\quad
\partial_{y_i y_j}L_c=\delta_{ij}N c_iw_i,\quad
\partial_{N y_i}L_c=c_iw_i,\quad \partial_{NN}L_c=0.
$$

## `electroneutral_trial_reconstruction`

- Section: Equilibrium / Tangent-plane testing
- Status: Implementation convention
- Source: Normalization and two linear charge/amount equations; #223 adopted numerical admission
- Description: Reconstruct finite positive supported u=exp(y). Ionic supports require both charge signs. Rescale charge groups and normalize; structural zeros remain zero. Liquid trials then correct arithmetic residuals on the largest positive and negative entries (ties by index). Charged Vapor trials use multiplicative log rescaling and omit the additive ionic correction to preserve trace coordinates.
- Change note: Reject loss of positivity/support or domain. Neutral supports normalize directly; a singleton is exactly one. Reevaluate energy at v=c(T,w)/η. This is numerical admission, not exact-real feasibility.

$$
\begin{aligned}
Q_+&=\sum_{z_i>0}z_i u_i,\quad Q_-=\sum_{z_i<0}|z_i|u_i,\\
\widehat u_i&=\begin{cases}u_i\sqrt{Q_-/Q_+},&z_i>0,\\u_i\sqrt{Q_+/Q_-},&z_i<0,\\u_i,&z_i=0,\end{cases}
\qquad w_i=\frac{\widehat u_i}{\sum_{j\in S}\widehat u_j},\\
\begin{pmatrix}1&1\\z_+&z_-\end{pmatrix}
\begin{pmatrix}\delta w_+\\\delta w_-\end{pmatrix}
&=\begin{pmatrix}1-\sum_i w_i\\-\sum_i z_iw_i\end{pmatrix},\qquad
(w_+,w_-)\leftarrow(w_++\delta w_+,w_-+\delta w_-).
\end{aligned}
$$

## `normalized_trial_reaction_affinity`

- Section: Equilibrium / Tangent-plane testing
- Status: Exact derivation
- Source: tangent_plane_trial_balance contracted with charge-conserving Bν=0; #223 adopted derivation
- Description: Normalized stationary trial affinity for a mole-changing reaction, products positive; μ_mech/RT+s is the total reduced chemical potential. The trial's normalized reaction direction differs from the amount-space reaction vector.
- Change note: π annihilates ν in amount space, not generally its normalized direction. The added full declaration restores original zero-affinity rows; the search trial need not have zero affinity.

$$
\mathbf B\boldsymbol\nu=0,\quad \mathbf z^{\mathsf T}\boldsymbol\nu=0
\ \Longrightarrow\
\boldsymbol\nu^{\mathsf T}\left(\frac{\boldsymbol\mu_{\mathrm{mech}}}{RT}+\mathbf s\right)
=k\,\mathbf e^{\mathsf T}\boldsymbol\nu,\qquad
\mathbf d_w=\boldsymbol\nu-\mathbf w\,\mathbf e^{\mathsf T}\boldsymbol\nu.
$$

## `present_phase_tangent_plane_identity`

- Section: Equilibrium / Tangent-plane testing
- Status: Exact derivation
- Source: Euler homogeneity of phase_helmholtz_amount_volume with chemical and pressure rows; #223 adopted derivation
- Description: Full present-phase identity on normalized compositions. a_V,p is the SI molar-volume column of residual reduced Helmholtz energy at unit amount; pressure-row error cannot be omitted. μ_mech contains ideal-density and residual chemical potentials without s.
- Change note: The charge term is retained until charge-neutral admission; approximate chemical equality alone does not force zero TPD.

$$
\begin{aligned}
r_{p,i}&=\frac{\mu_{\mathrm{mech},p,i}}{RT}+s_i-\pi_i-z_i\varphi_p,\qquad
r_{Pp}=v_p\left(\frac{1}{v_p}-a_{V,p}-\frac{P}{RT}\right),\\
t(\mathbf x_p,v_p)&=\mathbf x_p^{\mathsf T}\mathbf r_p
+\varphi_p\mathbf z^{\mathsf T}\mathbf x_p-r_{Pp}.
\end{aligned}
$$

## `tangent_plane_donor_descent`

- Section: Equilibrium / Tangent-plane testing
- Status: Exact derivation
- Source: First variation of system_gibbs_helmholtz_coordinates along conserved donor withdrawal; #223 adopted derivation
- Description: Withdraw βw from strictly positive donor p containing S and create βw in the trial, holding donor and trial molar volumes fixed. β and its capacity are mol; d_p is the per-mole slope of extensive G/RT.
- Change note: Compute the donor slope directly from physical gradients; chemical plane uncertainty cancels. A negative TPD without finite positive donor capacity is an unresolved plane trial.

$$
\begin{aligned}
\beta_{\max,p}&=\min_{i\in S}\frac{n_{p,i}}{w_i}>0,\\
d_p=\left.\frac{d(G_{\mathrm{split}}/RT)}{d\beta}\right|_{0}
&=\widetilde g(\mathbf w,v_t)+\mathbf s^{\mathsf T}\mathbf w
-\mathbf w^{\mathsf T}\left(\frac{\boldsymbol\mu_{\mathrm{mech},p}}{RT}+\mathbf s\right)+r_{Pp}\\
&=t(\mathbf w,v_t)-\mathbf w^{\mathsf T}\mathbf r_p
-\varphi_p\mathbf z^{\mathsf T}\mathbf w+r_{Pp}.
\end{aligned}
$$

## `tangent_plane_descent_uncertainty`

- Section: Equilibrium / Tangent-plane testing
- Status: Implementation convention
- Source: First-order propagation, Euler homogeneity and final-summation rounding; #223 adopted estimated-uncertainty criterion
- Description: ε_tg estimates trial reduced molar energy error; ε_pμi donor reduced chemical error; ε_pP pressure-row error; ε_f feasibility/representation correction; ε_fp arithmetic. Include supplied standard-term uncertainty separately in the applicable energy/chemical terms. F_n and F_V denote native association derivative-floor estimates in unit-mol amount and SI volume columns. For the single full candidate amount correction, freeze ν_j^0=μ_j^0/(RT) at the parent; d=m¹−m⁰ is the stored update, u=δ−d its update-rounding discrepancy, m²=N′x′ the physical representation and m³=N′exp(ℓ′) its reported log representation. C_n is the applied and representation displacement in mol; ε_f is extensive G/(RT) on the request amount basis. Update rounding contributes only to ε_f, retaining nonzero uncertainty when a stored update rounds away. Candidate and donor consumers share the final totals.
- Change note: Report association residual/conditioning and Euler defect at reevaluated points. A maximum floor may replace each column floor with its scale retained. γ_m bounds final summation rounding only, not the EOS evaluator. Local first-order estimates are not certified enclosures.

$$
\begin{aligned}
q_j&=|m_j^2-m_j^1|+|m_j^3-m_j^2|,\qquad
C_n=\sum_j(|d_j|+q_j),\qquad
\epsilon_f=\sum_j(|d_j|+|u_j|+q_j)|\nu_j^0|,\\
E_d&=\epsilon_{tg}+\sum_i|w_i|\epsilon_{p\mu i}
+\epsilon_{pP}+\epsilon_f+\epsilon_{fp},\\
\epsilon_{tg,\mathrm{assoc}}&\approx\sum_i|w_i|F_{t,n_i}+|v_t|F_{t,V},\qquad
\epsilon_{pP,\mathrm{assoc}}\approx|v_p|F_{p,V},\\
\epsilon_{\mathrm{Euler}}&=\left|a^{\mathrm{res}}-\sum_i w_i a^{\mathrm{res}}_{n_i}-v a^{\mathrm{res}}_V\right|,\qquad
\epsilon_{\mathrm{sum}}\le\gamma_m\sum_j|\mathrm{term}_j|,\quad
\gamma_m=\frac{m\varepsilon_{\mathrm{mach}}}{1-m\varepsilon_{\mathrm{mach}}},\\
t&<0,\quad d_p<-E_d\quad\text{admit numerical slope descent only with all terms available;}
\quad d_p+E_{\mathrm{bound}}<0\quad\text{requires justified total bounds for proof}.
\end{aligned}
$$

## `conserved_finite_phase_split`

- Section: Equilibrium / Tangent-plane testing
- Status: Exact derivation
- Source: Component-conserving donor split of system_gibbs_helmholtz_coordinates; #223 adopted finite-descent check
- Description: Finite split preserving every component and thus moieties and charge. The withdrawn donor is reclosed at the request pressure on its declared density branch; the trial keeps its own molar volume. Other phases are unchanged. ΔG/RT and E_Δ are extensive on the request's mol feed basis.
- Change note: Try β_max/2 and at most eight halvings, reevaluating physical points. Numerical descent requires actual conservation and a Gibbs decrease exceeding its propagated estimate; proof additionally requires justified feasibility and error bounds. Issue #390: a fixed donor volume is admissible for the minimized A+PV, but its mechanical penalty (g_vv>0 on a stable root) blocked useful finite descent within the bounded search; the reclosed donor volume also seeds the added declaration. With exact pressure closure and a smoothly continuing root, g_v=0 makes N_p g_v dv_p/dβ vanish, so the first-order donor slope is unchanged; no smooth path is claimed across a root-selection change.

$$
\begin{aligned}
\mathbf n'_p&=\mathbf n_p-\beta\mathbf w,\quad
\mathbf n_t=\beta\mathbf w,\quad N'_p=N_p-\beta,\quad
\mathbf x'_p=\frac{\mathbf n'_p}{N_p-\beta},\quad 0<\beta<\beta_{\max,p},\\
P_{\mathrm{EOS}}(T,\mathbf x'_p,v'_p)&=P,\quad V'_p=(N_p-\beta)v'_p,\quad V_t=\beta v_t,\quad
\mathbf n'_p+\mathbf n_t=\mathbf n_p,\\
\frac{\Delta G}{RT}&=(N_p-\beta)g(\mathbf x'_p,v'_p)+\beta g(\mathbf w,v_t)-N_p g(\mathbf x_p,v_p),
\quad g=\widetilde g+\mathbf s^{\mathsf T}\mathbf x,\\
\frac{\Delta G}{RT}&<-E_\Delta\quad\text{for numerically resolved finite descent with estimated uncertainty}.
\end{aligned}
$$

## `thermochemistry_from_reactions`

- Section: Reference thermochemistry / Declared derivations
- Status: Implementation convention
- Source: Minimum-norm conserved gauge as practised in MEA-Thermodynamics (evaluate_direct_absorption_heat.py, build_thermochemistry), made a declared native derivation
- Description: For an independent reaction-row matrix nu, the minimum-norm solution defines one declared gauge of reaction enthalpies. Adding a conserved-moiety vector preserves reaction enthalpies but does not determine absolute calorics.
- Change note: Gauge is a convention, not replacement evidence for missing enthalpies, entropies or heat capacities. Equilibrium invariance additionally requires consistent Gibbs reference functions. Native calorics on referenced reactions use `reference_species_calorics`, which anchors this gauge on supplied records and leaves only the ionic charge gauge.

$$
\boldsymbol{\nu}\,\mathbf{h}^{\circ}=\Delta_r\mathbf{h}^{\circ},\qquad \mathbf{h}^{\circ}=\boldsymbol{\nu}^{\mathsf T}\left(\boldsymbol{\nu}\boldsymbol{\nu}^{\mathsf T}\right)^{-1}\Delta_r\mathbf{h}^{\circ},\qquad \boldsymbol{\nu}\left(\mathbf{h}^{\circ}+\mathbf{B}^{\mathsf T}\mathbf{c}\right)=\Delta_r\mathbf{h}^{\circ}\ \text{for every }\mathbf{c}
$$

## `equilibrium_isobaric_heat_capacity`

- Section: Equilibrium / Caloric observables
- Status: Exact derivation
- Source: Total extensive enthalpy differentiated along the fixed-pressure equilibrium branch; [issue #225](https://github.com/tannerpolley/ePC-SAFT/issues/225), sections 1 and 5, owner-selected D3
- Description: Apparent isobaric heat capacity includes phase and reaction relaxation, normalized by the declared positive feed sum. Native actions use physical increments on a dimensionless path parameter; the auxiliary thermal increment is 1 K.
- Change note: Equilibrium owns the phase-amount/enthalpy Taylor products. Parameter fitting consumes the value and first mixed thermal actions, dividing out physical increments before the affine optimizer map. Driven T/P and a smooth regular branch are required. No second apparent heat-capacity action or physical calorimetry is established.

$$
H=\sum_p n_p h_p,\qquad N_f=\sum_i n_{f,i},\qquad
c_{P,eq}=\frac{DH[w_T]}{N_f\delta_T},\qquad \delta_T=1\ \mathrm{K},
$$
$$
Dc_{P,eq}[w]=\frac{D^2H[w_T,w]}{N_f\delta_T}
-\frac{DH[w_T]DN_f[w]}{N_f^2\delta_T},\qquad
\partial_k c_{P,eq}=\frac{D^2H[w_T,w_k]}{N_f\delta_T\delta_k}\quad(DN_f[w_k]=0).
$$

## `reference_species_calorics`

- Section: Reference thermochemistry / Declared derivations
- Status: Exact derivation
- Source: Gibbs-Helmholtz and Gibbs derivatives, using standard_pressure_density_conversion and the owner-selected anchored Gibbs construction in [issue #225](https://github.com/tannerpolley/ePC-SAFT/issues/225), sections 1–3, D1(a)/D2/broadened D5
- Description: Supplied species records anchor density-standard Gibbs energy, enthalpy and entropy. Record-less species and every consumed derivative are completed through the same stoichiometric system. Its remaining null space must be the charge gauge, which cancels in electroneutral phases and neutral directions. The chemical rows' minimum-norm representative preserves affinities but cannot replace these caloric anchors.
- Change note: The completion checks consistency with supplied records on the selected smooth interval. Each raw derivative retains its own admission status. The bounded reference inventory is T/P order at most four and T/P order at most three with one EOS-coordinate derivative. Source coefficients are evaluated at the path temperature. Missing anchors, an undetermined neutral combination, or a singular chart cannot be supplied by a gauge. Numerical availability is qualified by the retained notebook, not by this equation definition.

$$
g_{\rho,i}^{\circ}=h_{P,i}^{\circ}-T s_{P,i}^{\circ}
 +RT\ln\frac{RT\rho_0}{P_0}\quad(i\in S),\qquad \rho_0=1\ \mathrm{mol/m^3},\\
\boldsymbol\nu_D\mathbf f_D^{\circ}=\Delta_r\mathbf f_\rho-\boldsymbol\nu_S\mathbf f_S^{\circ},
\qquad\mathcal N(\boldsymbol\nu_D)\subseteq\operatorname{span}(\mathbf z_D),\\
L=\ln K_\rho=\ln K_{src}+A+C,\qquad
\Delta_r g_\rho=-RTL,\quad\Delta_r h_\rho=RT^2L_T,\quad
\Delta_r s_\rho=R(L+TL_T),\quad\Delta_r c_{p,\rho}=R(2TL_T+T^2L_{TT}),\\
h=\sum_i x_i(g_{\rho,i}^{\circ}-Tg_{\rho,i,T}^{\circ}+RT)+h^{res},\\
s=\sum_i x_i\left[-g_{\rho,i,T}^{\circ}+R-R\ln\frac{x_i\rho}{\rho_0}\right]+s_{TV}^{res},
\qquad c_p=\sum_i x_i(-Tg_{\rho,i,TT}^{\circ}+R)+c_p^{res}.
$$

For a raw T/P derivative of a reaction reference, all lower product terms remain:

$$
\partial_T^a\partial_P^b(-RT C)
=-R\left(TC_{ab}+aC_{a-1,b}\right).
$$

## `gibbs_conjugate_molar_volume`

- Section: Reference thermochemistry / Declared derivations
- Status: Exact derivation
- Source: Frozen-composition Gibbs derivatives and the anchored completion in [issue #225](https://github.com/tannerpolley/ePC-SAFT/issues/225), section 1, D1(a)
- Description: The volume conjugate to system pressure includes the record-anchored standard pressure derivative. It is distinct from the EOS mechanical volume. Internal energy and frozen heat capacity use that same Gibbs volume.
- Change note: Frozen raw T/P partials eliminate density and association at fixed composition before following the solved composition path. The local volume chart requires finite v_G and finite nonzero v_G,P, with a regular EOS root and site Jacobian. Record v_G>0 and v_G,P<0 as stable-chart diagnostics alongside quotient conditioning; these signs do not restrict mathematical admission or establish global reacting stability. Fixed-source-pressure standards have zero system-pressure derivatives and reduce to the mechanical definitions. A zero or unresolved denominator is diagnosed before a negative disposition; it is never replaced by the mechanical volume.

$$
v_G=g_P=\rho^{-1}+\sum_i x_i g_{\rho,i,P}^{\circ},\qquad
V_G=\sum_p n_p v_{G,p},\qquad u_G=h-Pv_G,\\
c_{v,G}=\left(\frac{\partial u_G}{\partial T}\right)_{v_G,\mathbf x}
=c_p+T\frac{v_{G,T}^2}{v_{G,P}},\qquad
v_{G,T}=\left(\partial_Tv_G\right)_{P,\mathbf x},\quad
v_{G,P}=\left(\partial_Pv_G\right)_{T,\mathbf x}.
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

## `finite_addition_enthalpy`

- Section: Parameter regression / Enthalpy observations
- Status: Exact derivation
- Source: Closed-system enthalpy balance for a declared finite feed addition, using EqIDs `total_enthalpy`, `ideal_enthalpy` and `residual_enthalpy`; accepted observation definition in issue #210
- Description: Signed total-system enthalpy increment per mole of total feed added between two declared equilibrium states, less the ideal-plus-residual feed enthalpy at its declared temperature, pressure, composition and branch. Each phase enthalpy uses EqID `total_enthalpy` on the declared thermochemistry's $P^{\circ}$ basis. Exothermic absorption is negative. A zero-pinned phase contributes zero and is not evaluated.
- Change note: Both endpoint problems share parameters, thermochemistry, ordered reactions and source reference. Their feed amounts differ by $n_{add}\mathbf y$ within the admission tolerance. Each endpoint may have its own specification and temperature, with referenced free-temperature endpoints consuming the shared native first caloric actions at their solved temperatures. Every feed species needs an ideal-gas record. The exact first action holds feed state and amount fixed, differentiates both endpoint phase amounts and enthalpies, and subtracts the feed residual action; ideal feed enthalpy is independent of fitted parameters. The observation is an enthalpy increment for free-pressure endpoints; interpreting it as measured heat is an energy-balance choice unless endpoints and feed are isobaric at one pressure.

$$
\Delta h_{add}=\frac{H^{(2)}-H^{(1)}}{n_{add}}-h_{feed},\qquad
H^{(k)}=\sum_{p\notin Z_k}n_p^{(k)}h_p^{(k)},\\
h_{feed}=h^{id}(T_f,\mathbf y)+h^{res}(T_f,P_f,\mathbf y;\theta),\\
\frac{d\Delta h_{add}}{d\theta}=\frac{1}{n_{add}}\sum_{k=1}^{2}\sigma_k\sum_{p\notin Z_k}
\left(\dot n_p^{(k)}h_p^{(k)}+n_p^{(k)}\dot h_p^{(k)}\right)
-\left.\frac{\partial h^{res}}{\partial\theta}\right|_{T_f,P_f,\mathbf y},\qquad
\sigma_1=-1,\quad\sigma_2=1.
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
- Change note: Unavailable without declared reference thermochemistry. These mechanical definitions apply to fixed-source-pressure standards. System-pressure references use the same record-anchored Gibbs volume in `gibbs_conjugate_molar_volume`.

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
- Change note: For fixed-source-pressure standards this mechanical quantity does not inherit the fixed-pressure mechanical-stability condition. System-pressure references use frozen c_v,G on the regular local Gibbs-volume chart in `gibbs_conjugate_molar_volume`; it is not a reacting heat capacity at fixed total volume.

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

## `log_fugacity_coefficient_parameter_derivative`

- Section: Properties / Pressure-closed parameter actions
- Status: Exact derivation
- Source: Implicit differentiation of `pressure`, `compressibility_factor` and `log_fugacity_coefficient`; extensive residual Helmholtz Hessian including implicit association.
- Description: First physical-parameter columns at fixed temperature, pressure and component amounts.
- Change note: Requires positive finite pressure-density slope. Parameters use their native units; the volume is m³ for the declared amounts. Numerical qualification belongs to the electrolyte parameter-gradient analysis.

$$
V_a=-\frac{\phi_{Va}}{n/V^2+\phi_{VV}},\qquad
 g_{i,a}=\phi_{ia}+\left(\phi_{iV}-\frac1V\right)V_a,
\quad \phi=A^r/(RT),\quad g_i=\ln\varphi_i.
$$

## `electrolyte_activity_parameter_derivative`

- Section: Properties / Electrolyte reference actions
- Status: Exact derivation
- Source: `mean_ionic_activity_coefficient`, `molality_activity_conversion`, `practical_osmotic_coefficient` and `log_fugacity_coefficient_parameter_derivative`.
- Description: Molality-scale mean ionic activity, solvent activity and practical osmotic first parameter columns.
- Change note: Ionic reference columns are the terminal-pair limit of `reference_molality_sequence`, the exact derivative of the returned reference. Solvent columns subtract the same-T/P pure-solvent column, including the implicit density, association and dielectric/Born derivatives. Physical coordinates, composition, T/P and molar masses are held fixed. At pure solvent activity columns are zero; unavailable Phi has unavailable columns. Energies multiply ionic reference columns by RT in J/mol.

$$
\frac{\partial\gamma_\pm^{(m)}}{\partial\theta_a}
=\gamma_\pm^{(m)}\frac{\sum_k\nu_k(g_{k,a}-g^\infty_{k,a})}{\sum_k\nu_k},
\qquad \frac{\partial\Phi}{\partial\theta_a}
=-\frac{g_{s,a}-g^\mathrm{pure}_{s,a}}{M_s q},\qquad
\frac{\partial a_s}{\partial\theta_a}=a_s(g_{s,a}-g^\mathrm{pure}_{s,a}).
$$

## `dependent_born_diameter_derivative`

- Section: Parameter regression / Dependent hydration coordinate
- Status: Exact derivation
- Source: Implicit-function theorem applied to `infinite_dilution_solvation_energy` and `log_fugacity_coefficient_parameter_derivative`.
- Description: A Born diameter holds a declared single-ion solvation energy, with exact first columns contracted through the independent cross-family affine chart.
- Change note: Independent scalar holds only; nonzero scaled H_d, a declared bracket and both energy and diameter-correction criteria are required. Final acceptance uses the stored diameter without solving the equality again. Target uncertainty is fixed, not estimated.

$$
H=RT\ln\varphi_i^\infty-\Delta G_{i,\mathrm{target}}^\infty=0,
\qquad\frac{d d_i}{d\theta_a}=-\frac{H_a}{H_d},
\qquad J^\mathrm{opt}_{jc}=\sqrt{w_j}\,\frac{dr_j}{dy_j}\sum_r\frac{\partial y_j}{\partial\theta_r}D_{rc},
\quad D_{rc}=\frac{d\theta_r}{dv_c}.
$$

For the explicit targets $(\sigma_i,d_{\mathrm{packing},i},d_{\mathrm{DH},i})=(q,0.88q,0.88q)$,
the implicit-function contraction includes every independent raw column:

$$
\frac{d d_{\mathrm{Born},i}}{d v_j}
=-\frac{s_j\left(H_{\sigma_i}+0.88H_{d_{\mathrm{packing},i}}+0.88H_{d_{\mathrm{DH},i}}\right)}{H_{d_{\mathrm{Born},i}}}.
$$

## `extensive_enthalpy_specification`

- Section: Equilibrium caloric specifications
- Status: Exact derivation
- Source: `total_enthalpy`, extensive phase summation and resolved-feed normalization; issue #224 slice A
- Description: Fixed pressure and one extensive enthalpy target H_spec [J] determine free temperature and phase amounts on declared neutral supports.
- Change note: Same-density total calorics use the supplied component anchors and standard pressure; the denominator uses the current feed and temperature.

$$
r_H=\frac{\sum_p N_p h_p^{tot}-H_{spec}}{N_{feed}RT}.
$$

## `extensive_entropy_specification`

- Section: Equilibrium caloric specifications
- Status: Exact derivation
- Source: `total_entropy`, extensive phase summation and resolved-feed normalization; issue #224 slice A
- Description: Fixed pressure and one extensive entropy target S_spec [J/K] determine free temperature and phase amounts on declared neutral supports.
- Change note: Same-density ideal and residual entropy share the supplied component anchors and standard pressure; zero pinned amounts contribute zero.

$$
r_S=\frac{\sum_p N_p s_p^{tot}-S_{spec}}{N_{feed}R}.
$$

## `critical_equations`

- Section: Equilibrium critical and envelope conditions
- Status: Source equation
- Source: Heidemann and Khalil 1980, AIChE J. 26, 769–779 (doi:10.1002/aic.690260510); Michelsen 1980, Fluid Phase Equilib. 4, 1–10, pp. 7–8, Eqs. 24–26; normalized amount perturbation as in the FeOs "Stability and critical points" documentation; issue #224 slice B
- Description: Necessary critical conditions of one neutral nonreactive phase on a positive feed z at fixed (T, V), with F = A/(RT), n the phase amounts, unknown mode u, D = diag(sqrt z) and w = D u; T and P are free and the pressure row closes P.
- Change note: Ideal parts Q_ij = δ_ij/n_i and F_ijk = −δ_ijk/n_i²; residual parts are the extensive residual contractions E2–E5. ∂c/∂u_j = 3 sqrt(z_j) F_nnn[w,w,e_j] and ∂²c/∂u_i∂u_j = 6 sqrt(z_i z_j) F_nnn[w,e_i,e_j]. A solution is a critical-equation candidate: the conditions certify no stability.

$$
M=DQD,\quad Q=F_{nn},\quad Mu=0,\quad u^{T}u-1=0,\quad c=F_{nnn}[w,w,w]=0.
$$

## `envelope_local_condition`

- Section: Equilibrium critical and envelope conditions
- Status: Source equation
- Source: Michelsen 1980, Fluid Phase Equilib. 4, 1–10, pp. 2–4, Eqs. 9 and 14–18 (selectable specification and local parametrization), adapted to the raw T/P chart; issue #224 slice B
- Description: One scalar row of an envelope cell (a finite feed phase and an incipient phase pinned at zero amount, T and P free): the selected comparison coordinate minus its target s, with T_unit = 1 K, P_unit = 1 Pa and K_k = x_k^inc/x_k^feed.
- Change note: ∂row/∂s = −1. The ln T and ln P rows have second derivatives −1/T² and −1/P²; the ln K row is a difference of log-composition coordinates and is affine.

$$
\ln(T/T_{unit})-s=0,\quad \ln(P/P_{unit})-s=0,\quad \ln x_k^{inc}-\ln x_k^{feed}-s=0.
$$

## `azeotrope_condition`

- Section: Equilibrium critical and envelope conditions
- Status: Source equation
- Source: equal phase compositions at a homogeneous azeotrope (Gibbs–Konovalov); FeOs `PhaseEquilibrium.binary_azeotrope` as the independent implementation; issue #571
- Description: The C − 1 rows of a binary azeotrope cell (a free feed phase and an incipient phase pinned at zero amount, exactly one of T and P free) that replace the first C − 1 species balances; the remaining balance is the total amount, and the declared feed composition is only the seed.
- Change note: The rows are affine in the log-composition chart and add no Hessian. The C-th equality follows from the two normalization rows.

$$
\ln x_k^{inc}-\ln x_k^{feed}=0,\quad k=\text{support}[0..C-2],\qquad \sum_p N_p\sum_i x_{p,i}-N_{feed}=0.
$$

## `solid_standard_gibbs`

- Section: Equilibrium pure stoichiometric solids
- Status: Exact derivation
- Source: dissolution S(s) = Σ_i ν_{s,i} A_i with Σ_i ν_{s,i} μ_i° − μ_s° = −RT ln K_s; the density-basis constant from `neutral_reaction_reference_transfer` at the source standard pressure P°; issue #573, ADR 0013
- Description: Standard Gibbs energy of a pure stoichiometric solid at unit activity on its pure-solid standard state, with g_i the compiled fluid standard vector, ln K_s^ρ the sourced dissolution constant converted on the charge-neutral fluid-side stoichiometry at the declared reference pressure P°, and V_s the declared solid molar volume.
- Change note: P° is required, so ∂g_s/∂P = V_s exactly at fixed T; a transfer at system pressure would add the dissolved ions' infinite-dilution partial molar volume. Only ln K_s is physical: a shift of g_i along the balance row space moves the objective by a constant. ln K_s^ρ is linear in T about the solved temperature (first-order T actions).

$$
g_s(T,P)=\sum_i\nu_{s,i}\,g_i(T)+RT\ln K_s^{\rho}(T,P^\circ)+V_s\,(P-P^\circ).
$$

## `solid_saturation_index`

- Section: Equilibrium pure stoichiometric solids
- Status: Exact derivation
- Source: affinity of dissolution at the solved fluid state; with `chemical_potential_stationarity`, Σ_i ν_{s,i} μ_i = RT λ·(Bν_s) because Σ_i ν_{s,i} z_i = 0; issue #573
- Description: Saturation index of solid s in decades; at the optimizer endpoint its scaled affinity is the bound multiplier, z_s = −A_s/RT, and the row g_s/RT − λ·Bν_s = 0 holds for a present solid.
- Change note: A Pinned(0) solid reports SI_s at the solved fluid state without a coordinate or row.

$$
\mathrm{SI}_s=\frac{A_s}{RT\ln 10},\qquad \frac{A_s}{RT}=\frac{\sum_i\nu_{s,i}\mu_i-g_s}{RT}=\boldsymbol\lambda\cdot(\mathbf B\boldsymbol\nu_s)-\frac{g_s}{RT}.
$$

## `solid_amount_complementarity`

- Section: Equilibrium pure stoichiometric solids
- Status: Exact derivation
- Source: KKT conditions of the Gibbs minimum with the bound n_s ≥ 0 and no stationarity row for n_s; issue #573
- Description: Complementarity between the solid amount n_s (mol) and its saturation index; Ipopt's returned bound multiplier of n_s on the extensive objective G/RT is z_s = −A_s/RT, so (n_s/N_feed) z_s is the barrier parameter at an interior-point endpoint.
- Change note: One endpoint classification makes the declaration square: present (n_s/N_feed > z_s, saturation row) or absent (n_s = 0). The rule replaced A_s/RT ≥ −ε after the barrier left a stable solid at −A_s/RT above ε for small excess feeds.

$$
n_s\ge 0,\qquad \mathrm{SI}_s\le 0,\qquad n_s\,\mathrm{SI}_s=0.
$$

## `solid_material_balance`

- Section: Equilibrium pure stoichiometric solids
- Status: Exact derivation
- Source: `material_balance` with the solids' moiety columns; issue #573
- Description: Conserved-moiety balance over the fluid phases and the free solid amounts; the solid term adds n_s g_s/RT to the objective and is linear in n_s, so it adds no curvature.
- Change note: No upper bound on n_s: the balance bounds it.

$$
\sum_p N_p\,\mathbf{B}\,\mathbf{x}_p+\sum_s n_s\,\mathbf{B}\boldsymbol\nu_s=\mathbf{B}\,\mathbf{n}_{\mathrm{feed}}.
$$
