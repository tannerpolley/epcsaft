# Equilibrium formulation

This document specifies required phase, reaction and reactive-phase mathematics.
[Equation definitions](../equations.md) own the displayed relations;
[algorithms](../science/algorithms.md) own their numerical solution.
[TRANSFER.md](../../../docs/greenfield/TRANSFER.md) records implementation gaps.
A proposed or registered equation is not a capability claim.

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

Bubble, dew, pure saturation, fixed-fraction and ordinary flash can be
constructed as data specifications of these rows. Names may help users
construct a problem; they must not select application-specific physics.

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

Report achieved balance, charge, normalization, pressure, stationarity/affinity
and phase-distinctness diagnostics on declared scales. A stationary candidate
is not automatically a stable equilibrium. A requested distinct phase set
that coalesces must report the degeneration explicitly.

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
Search/controller selection and a justified bound remain open work in the
transfer plan. Budget exhaustion never means proven stable.

## Sensitivity and validity

EqID: `implicit_sensitivity` applies on a smooth fixed branch and support with
a square full-rank row Jacobian. Solve its linear system rather than forming
an inverse. EqID: `implicit_second_action` defines second directional actions;
third actions differentiate the same equation again. All direct and implicit
terms of the requested observable must be included. Energy derivative order
is determined by the observable, as tabulated in the algorithms specification.

At rank loss, phase coalescence, support changes, branch switches or a
nonsmooth ionic boundary, report the unavailable action and its reason while
preserving valid central values. Conditioning and residuals qualify each
linear solve. Do not return an unexplained empty array as a derivative.

Solid phases, kinetics and augmented critical-point equations remain separately
justified extensions; this fluid formulation does not imply their support.
