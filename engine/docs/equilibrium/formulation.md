# Equilibrium formulation

This document specifies required phase, reaction and reactive-phase mathematics.
[Equation definitions](../equations.md) own the displayed relations;
[algorithms](../science/algorithms.md) own their numerical solution.
[GitHub issue roadmap](https://github.com/tannerpolley/ePC-SAFT-greenfield/milestones) records implementation gaps.
A proposed or registered equation is not a capability claim. Current completion
covers declared-phase specifications. PH/PS/UV and a phase-stability/discovery
algorithm are deferred in the plan; their mathematics below remains reference
material. Mechanical/local stability checks and truthful global-claim limits
continue to apply to current results.

## Energy, reference states and coordinates

At fixed temperature and pressure, fluid equilibrium minimizes the sum of
phase Helmholtz energies plus pressure-volume work subject to conserved
amounts and phase electroneutrality. Sources: Gross–Sadowski2001 Appendix A,
Pereira2012 Eqs. 2–8, Ascani2022 Eq. 5, and the retained Tsanas2019/2022
studies under `research/`. Use EqID: `phase_helmholtz_amount_volume` and
EqID: `system_gibbs_helmholtz_coordinates`. All phase values and derivatives
come from the same defining energy.

A phase carries amount (mol), normalized composition, and volume (m³), or an
equivalent packing-fraction coordinate (EqID: `phase_packing_volume`). Evaluate
intensive properties at unit amount before scaling extensive contributions.
An incipient phase has zero amount but a defined intensive composition and
volume per mole; do not evaluate a logarithm of its zero extensive amount.
The positive-amount logarithmic coordinate transformation is not defined at zero.

Eliminating volume by a global admissible minimum of Helmholtz energy plus
pressure-volume work gives the Gibbs formulation. A mechanically stable
stationary density root identifies only a local branch; it does not establish
the global minimum. Joint volume/composition coordinates avoid an inner root
choice but still require stability assessment of the final candidate.

The ideal mixing/density term supplies nonreactive equilibrium at fixed T.
Reactive or caloric work additionally requires declared standard-state offsets
and their temperature dependence. Use EqID: `standard_pressure_density_conversion`,
EqID: `ideal_gas_helmholtz`, EqID: `total_enthalpy`, and EqID: `total_entropy`.
Reaction correlations constrain a reaction subspace. A chosen gauge may define
a reference convention but cannot supply missing absolute caloric evidence.

## Multiplier units and signs

Chemical potentials are dimensional (J/mol). Balance multipliers lambda are
conjugate to the declared conserved-moiety row basis; lambda dot B_i has units
J/mol. Charge numbers z_i are dimensionless, so the charge multiplier varphi_p
in chemical stationarity has units J/mol. With stationarity written as
mu_i minus lambda dot B_i minus z_i varphi_p = 0, a physical electric potential
psi_p in volts corresponds to varphi_p = -F psi_p, where F is the Faraday
constant. Only potential differences are determined without a gauge.

The dual of reduced molar energy uses dimensionless multipliers alpha, with
potential differences divided by RT. Optimizer callback multipliers Lambda_opt
multiply the implemented constraint functions (including their declared scales);
they are distinct from chemical multipliers embedded in the primal coordinates.
The weighted-Hessian falsifier must vary both independently.

## Constraints and specification

The balance matrix spans conserved moieties; its nullspace is the admitted
reaction span. Reactions must conserve charge. Independent balance rows replace
reaction extents as primal coordinates; extents may be reconstructed afterward.
At each phase use EqID: `chemical_potential_stationarity`, EqID: `pressure_equality`,
EqID: `phase_normalization`, EqID: `phase_support`, EqID: `phase_electroneutrality`
and EqID: `material_balance`.

One charge row can be omitted only for an ion-containing phase of nonzero
amount when neutral feed, charge-conserving reactions and balances imply it.
A zero-amount phase does not carry that implication. All returned charged
phases must satisfy electroneutrality. Support removes coordinates; explicitly
empty support must not silently mean all species. A support unable to form a
neutral ionic phase is inadmissible.

For the retained fixed-phase specification with k pinned phase amounts and f
free intensives, the structural count is rows minus unknowns = k − f. A square
system requires k = f. The admitted phase-count condition is Π − k ≤ rank B.
Reject dependent reactions, charged feed, impossible supports, inconsistent
pins, invalid domains and deficient structural rank with a named reason.
Counts and a generic-point rank check do not prove existence, uniqueness or
regularity at the computed solution. Energy specifications add rows/unknowns
and require their own count derivation; they cannot inherit this rule blindly.

The resolved feed n_feed is the amount basis of the whole problem, and returned
phase amounts are on that same basis. An amounts feed is absolute mol and is
never normalized. A mole-fraction feed resolves to n_i = x_i, mol per mol of
feed; a mass-fraction feed resolves to n_i = w_i / M_i, mol per kg of feed. A
fraction feed must sum to one within the composition-sum policy value; it is
rejected otherwise and never renormalized. A molality feed is absolute mol:
solute entries are mol per kg of solvent times the stated solvent mass, and
solvent entries are relative molar amounts rescaled so the solvent mass equals
the stated value. One declared solvent needs no amount. Several declared solvents
with no positive solvent amount are rejected, because their ratio is a physical
input; a declared solvent left at zero beside a supplied one stays zero.

A pinned amount is feed-relative, N_p = v_p Σ n_feed with 0 ≤ v_p ≤ 1. When the
declared reaction space conserves total moles, rank([B; e]) = rank(B) with
e ∝ (1, …, 1), the pins additionally satisfy Σ v_p ≤ 1, and Σ v_p = 1 when every
phase is pinned; a nonconserving reaction space admits no universal total-mole
sum rule. Reaction participation is sign independent: a species with nonzero
stoichiometry is retained in support regardless of orientation. An all-pinned
nonreactive declaration compiles square by count but its retained balance and
normalization rows are dependent, Σ_i balance_i = Σ_p v_p normalization_p, so it
reports generic rank loss and no nonlinear solve is run.

Ordinary reaction targets are a constant affinity, a scalar K, or the P3
correlation ln K(T) on the EOS density basis, activity (ρ x_i / 1 mol m⁻³)
exp(μ_i^res / RT). A correlation on another basis requires the explicit neutral
reference declaration below. An ordinary scalar K is a
temperature-independent target: its series is −RT ln K, −R ln K, 0, 0, so it
states zero density-basis reaction enthalpy and differs from a correlated target
by Δ_r h / T at first temperature order. A temperature-dependent source constant
is therefore stated as a correlation on the EOS basis, not as its value at one
temperature; the engine does not detect that misuse. The two native conversions,
an ideal-gas pressure basis (a − Δν ln(R T_ref ρ₀ / P₀), c − Δν) and a
standard-state shift with temperature-independent Δh and Δs (a − Δs / R,
b + Δh / R, old minus new), have an exact polynomial image of
a + b/T + c ln(T/T_ref) + dT. A condensed-phase reference (molality or infinite
dilution, solution molarity, a pure-liquid mole-fraction reference) transfers
through the residual potential of the reference state and has none. A supplied
constant standard Gibbs vector must satisfy every reaction target at the
evaluation temperature by the same criterion as the inferred vector, and a
free temperature admits it only when the targets are constant. Supplied species
thermochemistry is compared with each target through Gibbs, enthalpy and heat
capacity on the density basis; a mismatch is a reference inconsistency, not a
preference for either input.

The implemented neutral-reference slice uses one pure neutral solvent and a
positive charge-neutral all-solute dilution path. An explicit source-basis enum
selects common-molality or raw mole-fraction correlations; a nonempty source ID
records provenance. EqID: `neutral_reaction_reference_transfer` converts each
source at its declared fixed or system reference pressure, then the existing
reaction/SVD owner compiles the density-basis target. At fixed temperature,
system-pressure references follow each current nonlinear pressure; a fixed
source pressure has zero derivatives with respect to system pressure.
The liquid selector chooses the dense branch when several stable roots exist;
a unique stable root follows the density owner's unspecified-label semantics.
Reference qualification must therefore establish the physical solvent state
for the application; the selector supplies no universal phase classification.
This qualified finite terminal reference supports fixed-T/P central solves and
feed actions. Its exact first and second pressure derivatives also support
square pinned free-pressure cells and pressure/feed actions through order two
(EqID: `neutral_reaction_reference_pressure_jet`). Its first temperature
derivative at the reference pressure supports order-one temperature and mixed
T/P actions (EqID: `standard_state_transfer_temperature_derivative`). Feed-only
actions at fixed T/P retain order three. Declared-branch continuation along
fixed T, fixed P and feed uses the same first-order chains for its tangent.
Third-order actions involving pressure, temperature actions above order one,
free temperature, raw lever derivatives and parameter actions remain unavailable. Ideal/Total
enthalpy and heat capacity follow EqID: `reference_species_calorics` at fixed
pressure; other caloric quantities remain unavailable. Ordinary unreferenced reactions
retain their existing consistency checks and derivative behavior.

Bubble, dew, pure saturation, fixed-fraction and ordinary flash can be
constructed as data specifications of these rows. Names may help users
construct a problem; they must not select application-specific physics.
A single declared phase and one-component coexistence are ordinary cells of
the same specification. Admission and numerical qualification must cover
phase/support, reaction, charge, pin and free-intensive combinations explicitly;
structural acceptance of a cell is not evidence that it solves.

## Stationarity, feasibility and stability

For free phase amounts at fixed T/P, the rows are first-order conditions of the
constrained Gibbs minimum. Phase potential equality and zero reaction affinity
follow from stationarity: EqID: `interphase_potential_equality` and
EqID: `reaction_affinity`. The positive-definite reduced Hessian
(EqID: `reduced_hessian_stability`) is a sufficient strict local condition on a
regular constraint manifold, not a necessary condition for every local minimum
or a global certificate.

Pinned amounts and free intensives define a feasibility specification using
intensive chemical-potential and pressure rows. They are not all the KKT rows
of the original fixed-T/P minimization. Their algorithm must satisfy the actual
constraint Jacobian and, when advertised, weighted Hessian requirements in
[ALG-EQUILIBRIUM](../science/algorithms.md#alg-equilibrium).

For free phase amounts the classifying stationarity residuals are the
first-order conditions in the solved logarithmic chart, N_p x_i μ̃_i / N_feed for
composition and N_p / N_feed times the pressure row for packing, which are the
quantities the constrained minimizer controls; the unscaled chemical and
pressure rows are retained as diagnostics and are certified only to the
criterion times N_feed / (N_p x_i), which is unbounded where an amount or a
fraction vanishes. The amount-stationarity component (phase Gibbs energy equal
to the balance dual) is not a classified row; it is bounded by the optimizer's
dual tolerance and checked by the retained Gibbs-dual identities. A pinned specification classifies its unscaled rows, which
are its constraints. Report achieved balance, charge, normalization, pressure,
stationarity/affinity and phase-distinctness diagnostics on declared scales. A stationary candidate
is not automatically a stable equilibrium. A requested distinct phase set
that coalesces must report the degeneration explicitly, including a pair whose
identities swap between consecutive states of a followed branch (passage
through a critical point).

## Discovery and global claims

EqID: `tangent_plane_distance` defines stability relative to the candidate
supporting plane, over every admissible trial composition and volume with
appropriate charge/conservation coordinates. A feasible negative trial is an
instability witness. Failure to find one in finite search is inconclusive.
Sources: Michelsen1982; Pereira2012; Perdomo2025. Reactive systems use conserved
moiety multipliers in the plane, consistent with the declared reference gauge.

The exact global dual infimum (EqID: `helmholtz_lagrangian_dual`, EqID:
`dual_inner_problem`) can supply a lower bound. A local inner minimizer is an
upper bound on that infimum; substituting it as a certified lower bound is
invalid. A finite outer approximation alone does not close the global gap.
EqID: `outer_approximation_bound` and the historical tunneling/candidate
relations describe source methods, not a selected current discovery algorithm.
Search/controller selection and a justified bound are deferred future work in
GitHub future-work issues, not prerequisites for declared-phase acceptance. Budget exhaustion never means proven stable.

## Sensitivity and validity

EqID: `implicit_sensitivity` applies on a smooth fixed branch and support with
a square full-rank row Jacobian. Solve its linear system rather than forming
an inverse. EqID: `implicit_second_action` defines second directional actions,
the third-order relation and the diagonal coefficient recurrence that realises
them; admission at an accepted result uses the classified residual of the
solved chart. All direct and implicit
terms of the requested observable must be included. Energy derivative order
is determined by the observable, as tabulated in the algorithms specification.

At rank loss, phase coalescence, support changes, branch switches or a
nonsmooth ionic boundary, report the unavailable action and its reason while
preserving valid central values. Conditioning and residuals qualify each
linear solve. Do not return an unexplained empty array as a derivative.

Solid phases, kinetics and augmented critical-point equations remain separately
justified extensions; this fluid formulation does not imply their support.
