# Equilibrium formulation

This document specifies required phase, reaction and reactive-phase mathematics.
[Equation definitions](../equations.md) own the displayed relations;
[algorithms](../science/algorithms.md) own their numerical solution.
[GitHub issue roadmap](https://github.com/tannerpolley/ePC-SAFT/milestones) records implementation gaps.
A proposed or registered equation is not a capability claim. Declared-phase
specifications and the explicit native tangent-plane test/fixed-T/P discovery
path adopted in #223 have retained numerical qualification for the stated
neutral and admitted ionic/reacting cases. The independent delivered review
supports these retained cases and stated limits; PH/PS/UV remain deferred.
Mechanical/local stability checks and truthful global-claim limits apply to
both declared and discovered results.

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
sum rule. Reaction participation is sign independent. A participant is removed
only when a globally nonnegative conserved functional has an exactly zero feed
inventory and a positive coefficient for that species. Exact native dot products
verify the functional against every original reaction and feed component; an LP
objective rounded to zero alone does not establish absence. Positive trace
inventories retain their species. An all-pinned
nonreactive declaration compiles square by count but its retained balance and
normalization rows are dependent, Σ_i balance_i = Σ_p v_p normalization_p, so it
reports generic rank loss and no nonlinear solve is run.

On a certified face, restrict state coordinates to the surviving species A and
select independent echelon row operations L from B[:, A]. Keep the full-column
map L B for balances and feed controls; only state rank and total-mole conservation
use its restriction to A. Preserve discarded inventory restrictions for action
admission. The feasible reaction space is ker(B[:, A]), embedded with zero
amounts outside A, rather than reactions with stoichiometric columns deleted.
Full reaction declarations and species reference/caloric completion stay fixed
before restriction to A. Removed species return literal zero amounts and
fixed-face actions. Feed directions must preserve vanished conserved inventories;
individual removed feed slots need not independently vanish.
Supplied phase composition guesses follow the same support restriction before
the existing supported-composition normalisation.

For a declared liquid/vapor boundary, an incipient phase means the limit of an
emerging positive phase on the same exact inventory face. For each nonnegative
certificate w, conservation gives 0 = Σ_p N_p wᵀx_p. Every positive phase therefore
has wᵀx_p = 0, and the closed composition face contains its incipient limit.
A zero amount alone does not imply this restriction for an arbitrary intensive
phase. Restrict both finite and incipient supports before differentiation;
the continuous energy extension x ln x → 0 is differentiated only along surviving
coordinates. No chemical row for a removed species survives. A discarded balance
row vanishes on A, so its potential contributes nothing to surviving stationarity.
The original reaction/reference completion and ker(B[:, A]) therefore define
the same constrained energy and stationary variations on this face.

Automatic incipient restriction is admitted only with exactly two phases of
opposite liquid/vapor kinds, fixed positive T, free P, exactly one `Pinned(0)`
amount and one free amount, and a nonempty effective incipient support whose
every species is uncharged. Outside this domain, an intersection with certified
absence retains `zero_inventory_incipient_support`, naming the species, phase
and excluded condition. Other empty supports remain refused. The reduced square
Jacobian must be regular at the solved endpoint; generic rank alone does not
qualify it. Singular or coalesced endpoints and charged incipient extensions
are unqualified. This face equivalence neither proves uniqueness nor establishes
global stability. Incipient supports without certified absence retain their
existing admission, and finite additions may connect separately solved faces.

The conserved-inventory bound call admits no additional charge, energy or chemical
conditions. On the union of compiled supports it maximizes species amount under
$Bn=b$, $n\ge0$, with structural zero outside that union. Amount results are in
mol; their separate feed-normalized value is $u_i/N_{\mathrm{feed}}$.
When reactions change total moles this value differs from a mole-fraction maximum.
For an explicit supported fraction floor $f\ge0$, $S f\le1$, fraction bounds use
$y=n/N$, $\tau=N_{\mathrm{feed}}/N$: $By=(b/N_{\mathrm{feed}})\tau$,
$\sum y_i=1$, $y_i\ge f$ on support, $y_i=0$ outside, $\tau\ge0$.
These are aggregate inventory fractions. They bound an individual liquid fraction
only for a sole finite liquid with a zero-amount vapor, as in the #237 Born cases.
They establish neither thermodynamic equilibrium nor stability, and no derivative
of the potentially nonsmooth maximum is supplied.

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
reaction/SVD owner compiles the density-basis target. References refresh at
each current nonlinear temperature, and system-pressure references also
follow each current pressure. A fixed source pressure retains thermal motion
and has zero derivatives with respect to system pressure.
The liquid selector chooses the dense branch when several stable roots exist;
a unique stable root follows the density owner's unspecified-label semantics.
Reference qualification must therefore establish the physical solvent state
for the application; the selector supplies no universal phase classification.
This qualified finite terminal reference supports fixed-T/P central solves and
feed actions. Its exact first and second pressure derivatives also support
square pinned free-pressure cells and pressure/feed actions through order two
(EqID: `neutral_reaction_reference_pressure_jet`). Its first temperature
derivative at the reference pressure supports order-one temperature and mixed
T/P actions (EqID: `standard_state_transfer_temperature_derivative`); its second
and mixed temperature-pressure derivatives extend temperature actions, alone or
mixed with feed or pressure, to order two (EqIDs:
`standard_state_transfer_temperature_second_derivative`,
`standard_state_transfer_temperature_pressure_derivative`). Feed-only
actions at fixed T/P retain order three. Declared-branch continuation along
fixed T, fixed P and feed, or a straight physical EOS/Born parameter path, uses
the same first-order chains for its tangent. All other declarations and
source conventions stay fixed on the parameter path.
Third-order actions involving pressure or temperature remain unavailable.
Raw referenced lever requests retain their existing refusal. Referenced
one-free-intensive pinned cells use current-point T/P row jets and
A1/A2 with at most one parameter-bearing direction, including T/P–parameter
terms. No external direction may drive its solved intensive; simultaneous free
T/P still lacks the additional constraint for a regular square chart.
Record-anchored Gibbs h/s/u, frozen c_p/c_v and Gibbs-conjugate volume use the
shared caloric composition and its consumed reference columns. The retained
prescribed-T/P evidence and subsequent free-intensive F6 qualification have
separate source identities and domains.
Ordinary unreferenced reactions
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

EqID: `tangent_plane_distance` defines the reduced molar distance from the
candidate's conserved supporting plane. Sources: Michelsen1982 Part I;
Pereira2012 Eq. 12; Perdomo2025 §3.3 Eqs. 37–40; Ascani2022 §3.3 Eqs. 33–35.
The reactive extension follows the homogeneous energy and conserved-moiety
identities, EqID: `tangent_plane_gibbs_dual_identity`. Temperature is in K,
pressure in Pa, amounts in mol and molar volume in m³/mol; reduced energies and
chemical multipliers are in RT units. The plane uses the result's chemical
balance multipliers, recovered for all-free cells or primal for pinned cells,
not optimizer constraint multipliers. Freeze actual standard terms, including
reference transfer and gauge, at the solved T/P snapshot.

Only explicit `Discover()` requests run discovery. Only explicit
`stability_test(...)` calls test an accepted declared result. Declared solves,
continuation, parameter fitting and ASL never invoke either implicitly.
Discover once at initialization, then declare the returned phases for subsequent
calculations. Explicit testing of pinned/free-intensive results releases pins
only for diagnosing Gibbs descent at their solved fixed T/P; it says nothing
about feasibility along the original pinned path. Zero-amount phases cannot
supply a donor. Discovery with free T/P remains excluded.

The native test compiles a one-free-phase unit-amount trial through the shared
exact-curvature driver, without reactions or feed-moiety proportions. Its
internal reduced balance space spans normalization and charge, retaining the
row transformation (EqID: `tangent_plane_trial_balance`). The stationary trial
identity uses the multiplier contraction with its transformed right-hand side;
species ordering need not leave a total-amount row first. The frozen standard
shift has constant coefficients but retains exact log-chart curvature
(EqID: `tangent_plane_shift_log_curvature`). A normalized reacting trial is a
search composition, not necessarily reaction equilibrium: EqID:
`normalized_trial_reaction_affinity`. The full added declaration restores every
original reaction, reference input and affinity row.

Raw optimizer points are not physically feasible by construction. Native
normalization and ionic group rescaling followed by a two-pivot arithmetic
correction admit a nearby trial (EqID: `electroneutral_trial_reconstruction`).
Recompute its molar volume from composition and packing and reevaluate shared
energy. Retain raw residuals, displacement and reconstructed point separately.
Repair remaining candidate balance/charge defects in supported amount variables,
then check positivity, support and domain and reevaluate; keep this corrected
candidate separate from the public result. No failed initialization supplies a
trial point. Arithmetic reconstruction does not assert exact-real feasibility.

A negative plane value requires a strictly positive donor with the entire trial
support and finite withdrawal capacity. The actual fixed-molar-volume donor
slope includes chemical and pressure corrections (EqID:
`tangent_plane_donor_descent`); present phases obey the full chemical-plus-pressure
identity (EqID: `present_phase_tangent_plane_identity`). Association derivative
floors, Euler defect, feasibility correction, supplied standard-term uncertainty
and arithmetic contribute separately to an estimated slope uncertainty
(EqID: `tangent_plane_descent_uncertainty`). A resolved negative slope and plane
value also require a retained finite inventory-conserving split whose computed
Gibbs decrease exceeds its propagated estimate (EqID: `conserved_finite_phase_split`).
Optimizer convergence or positive local trial inertia is not necessary for this
descent evidence.

The adopted claim is **numerically resolved descent with estimated uncertainty**.
A **proved instability witness** additionally requires justified feasibility
and total error bounds; current association floors are not such enclosures.
Missing error terms, near-zero margins and failed reconstruction remain
indeterminate. A negative trial without an eligible donor remains an unresolved
plane trial. Completed finite testing reports only **no witness in the recorded
trial set**, with unresolved trials named; it never certifies global stability.
Restricted supports and density-label starts limit that set, and positive local
curvature does not exclude a remote split. [ALG-STABILITY](../science/algorithms.md#alg-stability)
owns enumeration, conserved addition/removal and separate work budgets.

Charged trial supports are admitted only when the declared bulk-permittivity
form reaches the vacuum limit, $\epsilon_r\to1$ as $\rho_m\to0$.
The existing `DensityTemperature` (Schreckenberg) and `DipolarPerturbation`
(Rueben) forms satisfy this condition by their definitions in
[full EOS §6.2 and §6.4](../science/full-eos.md#6-bulk-dielectric-and-debye-huckel),
EqID: `density_correlated_bulk_permittivity` and EqID:
`dipolar_perturbation_bulk_permittivity`. `IonSuppression` (including
solvent-only, ion-fraction and ion-specific choices) and `ComponentMixing`
are density independent and are refused for charged trials. Their electrolyte
terms are defined for condensed phases; a free-volume trial search cannot infer
that domain from a Liquid label or mechanical pressure-root admission.

One native support-admission rule reads the model's declared form before any
trial root or optimizer invocation. A requested family containing a charged
support under a refused form raises `ParameterError` with
"permittivity model has no low-density limit; declare phases". The same rule
runs before discovery declaration solves. Caller-supplied neutral-only supports
remain admitted on these models. No density cutoff, optimizer-success filter,
parameter substitution or model change is introduced. Declared solves,
continuation, fitting and ASL retain their existing domain and calculations.
The vacuum limit is a necessary admission condition for this search, not
physical validation or a global stability certificate. The previously evaluated
solvent-only gas-density descent remains retained numerical evidence and is
classified as refused under this owner-adopted rule.

The exact global dual infimum (EqID: `helmholtz_lagrangian_dual`, EqID:
`dual_inner_problem`) can supply a lower bound. A local inner minimizer is an
upper bound on that infimum; substituting it as a certified lower bound is
invalid. A finite outer approximation alone does not close the global gap.
EqID: `outer_approximation_bound` and the historical tunneling/candidate
relations describe source methods, not a selected current discovery algorithm.
The adopted finite native controller does not implement those global-bound
methods. Justified global bounds, replacement at the balance-rank phase limit
and generalized multi-donor/reactive withdrawal remain separately scoped work.
Budget exhaustion never means proven stable. The [qualification notebook](../../../analyses/2026-tpd-phase-discovery/notebook.qmd)
retains neutral comparisons, admitted density-dependent ionic/reacting checks,
density-independent charged-search refusals and unchanged declared work.
These are numerical verification; physical validation and a global certificate
do not follow from this doctrine.

Reuse of discovered declarations is numerically qualified for the retained
owner-enumerated cases. It is not a guarantee for every amount scale or initial
point. A 100 mol neutral discovery reproduces the one-mol intensive split and
extensive Gibbs scaling, while both cold declared solving and supplied-result
public reuse coalesce and are rejected by the distinct-phase check. This is a
known declared-driver conditioning limit, not phase nonexistence or failure of
the discovered-state scaling. [#243](https://github.com/tannerpolley/ePC-SAFT/issues/243)
owns investigation; no shared-driver change is selected in #223.

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
