# Equilibrium formulation

This document specifies required phase, reaction and reactive-phase mathematics.
[Equation definitions](../equations.md) own the displayed relations;
[algorithms](../science/algorithms.md) own their numerical solution.
[GitHub issue roadmap](https://github.com/tannerpolley/ePC-SAFT/milestones) records implementation gaps.
A proposed or registered equation is not a capability claim. Declared-phase
specifications and the explicit native tangent-plane test/fixed-T/P discovery
path adopted in #223 have retained numerical qualification for the stated
neutral and admitted ionic/reacting cases. The independent delivered review
supports these retained cases and stated limits; PH/PS use the bounded neutral formulation below; UV remains deferred.
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
equivalent packing-fraction coordinate (EqID: `phase_packing_volume`). The declared
neutral ideal-gas base uses the log inverse unit-amount volume coordinate from
EqID: `system_gibbs_helmholtz_coordinates`; its occupied fraction is zero. It
admits one vapor phase at fixed positive T/P. Caloric conditions, solved-state
caloric actions, source-reference transfer and phase searches are deferred and
refused before calculation. Evaluate
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

## Ideal liquid chemistry

`IdealMoleFraction` declares species identifiers and charges without EOS parameters.
At fixed positive T/P, one free liquid uses EqID: `ideal_mole_fraction_reaction`
and EqID: `ideal_reaction_amount_coordinates`. Scalar constants are dimensionless
mole-fraction constants; laws identify `ideal-mole-fraction` explicitly. EOS-basis
laws, reaction/source references, supplied species Gibbs energies and affinity
targets are refused. The positive logarithm reference in `ReactionLogPolynomial`
may lie outside its evaluation interval; evaluation outside the interval remains refused.

The default `activity_basis="mole_fraction"` retains the equations above.
With `activity_basis="molarity"`, EqID: `ideal_molarity_reaction` requires a finite
positive `solution_volume_m3` and an absolute `Amounts` feed. Scalar K is the
dimensionless concentration constant on the 1 mol/L standard; correlations declare
`ideal-molarity-1-mol-per-litre`. No basis conversion is inferred.
`unit_activity_solvent` optionally names one neutral species on the positive liquid
support. Its activity is one; its amount and stoichiometry remain in the balances.
Volume and unit-solvent declarations are refused on the default basis.

A neutral nonnegative mol or mole-fraction feed defines the default conserved inventories.
Every reacting species requires a positive attainable inventory; an inventory face
forcing one absent is refused. The native inventory bounds limit log amounts,
including amounts greater than the feed reference where attainable. A supplied
positive normalized composition guess is a native cold seed. Previous-result
starts are refused before compilation or reaction-law evaluation.

One native row expression supplies values and exact derivatives to the existing
feasibility/Ipopt lifecycle and post-solve polish. The returned phase amount is
the current species total in mol, with normalized true species fractions. The
objective is zero; it is not physical Gibbs energy. Density and packing arrays
are empty. Endpoint Jacobian rank is checked. EOS-only solved-state actions,
stability searches and caloric conditions are unavailable for this model.

The explicit `reaction_activities(model, result)` call reads dimensionless species
activities and their logarithms at each returned phase composition and density,
including a finite candidate that misses its requested criterion. It performs
no calculation or storage during an unrequested solve, continuation or fit.
The model must match the result, which must retain its original declaration
and physical arrays. Records follow phase/reaction order; `species` arrays use
component order. A zero phase amount does not remove its defined composition.

EqID: `reaction_activity_quotient` governs the products-positive `log_quotient`.
Unreferenced Helmholtz reactions, including ideal gas and species-derived
targets, return density-standard activities. Ideal mole-fraction reactions
return exactly the mole fractions. Molarity reactions instead use returned phase
amounts times fractions and the shared fixed-volume expression, with unit solvent
activity where declared. Concentration pH is an analysis quantity under a declared
individual-ion convention, not an Engine result. Source references retain their declared
pressure and pure/dilute species standards; common molality subtracts the
solvent-mass/standard-molality factor for every solute. `basis`,
`reference_pressure_pa` and the `species` solvent/molality fields identify
these conventions. Ionic values use the EOS single-ion convention; only
charge-neutral combinations have convention-independent meaning. Values come
from the state and reference owner, independently of targets and residuals.

`species.causes` names unavailable component quantities; `cause` names an
unavailable quotient. A finite log value survives activity exponentiation
underflow or overflow, with the activity absent. A structural zero has activity
zero and log activity negative infinity when its reference exists. Zero
stoichiometric coefficients are skipped, and conflicting infinite terms make
the quotient unavailable. This call supplies no pressure-basis quotient, pH,
derivative, solid activity, persistent cache or fitting payload.

EqID: `effective_water_amounts` converts a criterion-meeting result to named
molecular amounts after checking model identity and the complete species partition.
For the MDEA example, EqID: `mdea_reaction_conserved_amounts` supplies independent
inventory checks. These are numerical verification of conserved hydronium chemistry,
not physical validation of ideal ionic activities, reproduction of Figure 15 or
agreement with the paper's nonconserved printed water balance. Hot-condition
extrapolation remains deferred; only declared evaluation intervals are admitted.

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
An explicit charged Vapor declaration that would select Vapor as the
omitted-charge phase is refused with `charged vapor cannot be the omitted-charge phase`.
Internal charged Vapor trials instead carry their directly constructed pure-charge
balance. Their numerical log-charge representation is confined to Newton
correction and classification (EqID `tangent_plane_trial_balance`); physical
constraints and Ipopt derivatives are unchanged.

A zero-amount phase does not carry that implication. All returned charged
phases must satisfy electroneutrality. Support removes coordinates; explicitly
empty support must not silently mean all species. A support unable to form a
neutral ionic phase is inadmissible.

A neutral feed is one whose mole fractions x_feed = n_feed / N meet the state
rule of full-eos.md §12, |Σ z_i x_i| ≤ 10⁻¹² (absolute); any other feed is
refused as charged. The row-free phase p₀ carries the feed charge Q amplified by
N/n₀, which is unbounded, so an admitted feed is solved with its ions of the
excess sign scaled by min(C, A)/max(C, A), where C = Σ_{z_i>0} z_i n_i and
A = Σ_{z_i<0} |z_i| n_i. A feed with C = A, including one with no ions, is left
unchanged. The scaling preserves signs and zero amounts, changes the charge by
exactly −Q and the amounts by Σ|Δn_i| ≤ |Q| ≤ 10⁻¹² N (every ion has
|z_i| ≥ 1); each balance row against the declared feed moves by at most
max|B|·10⁻¹². Ions of a single sign at that charge are scaled to zero. The
charge left in a returned phase is the solver's own: the charge component of
the balance-row residual in p₀ and each charge row's residual elsewhere. A
charge row bounds its phase's charge by the row criterion. In p₀ the balance
residual is normalized by N, so its charge is that criterion amplified by
N/n₀, like the feed charge. Neither is bounded by the state rule; the largest
measured is 3.2×10⁻¹⁴, 31 times under the state rule.

For k pinned phase amounts, f free intensives and e independent caloric rows,
rows minus unknowns = k + e − f; a square system requires f = k + e.
Ordinary specifications retain Π − k ≤ rank B. PH/PS adds one independent
balance, permitting Π − k ≤ rank B + 1. Generic-point rank does not prove
existence, uniqueness or regularity at the solution. Reject dependent reactions,
charged feed, impossible supports, inconsistent pins and invalid domains.

PH/PS declares exactly one extensive H_spec [J] or S_spec [J/K], free T and
fixed P, free phase amounts and neutral declared supports. EqIDs:
`extensive_enthalpy_specification`, `extensive_entropy_specification`.
Every consuming state requires complete component records and interval admission.
Ordinary density-basis reactions retain Δh = t − Tt′ and Δcp = −Tt″:
constant affinity t = a has Δh = a, while constant K has Δh = 0.
Charged calorics and source-referenced free-T reactions are excluded.


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

When no target is supplied, participating species thermochemistry supplies the
density-standard reaction Gibbs series through order three, preserving missing
records and excluded-interval causes. This source admits ordinary fixed-T/P
problems without solids, caloric conditions or source-reference transfer.
Attaching records alone does not select absolute species standards: with neither
reactions nor a supplied vector the compiler uses zero species standards.

Ordinary supplied reaction targets are a constant affinity, a scalar K, or the P3
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

### Pure stoichiometric solids

A pure solid s (ADR 0013) is not a fluid phase: it has no composition, density
root, pressure, normalization or charge row. It declares a charge-neutral
dissolution stoichiometry ν_s over model components, exactly one sourced ln K_s
or P3 law with its reaction reference, the reference's source standard pressure
P° (required) and its molar volume V_s. EqIDs: `solid_standard_gibbs`,
`solid_saturation_index`, `solid_amount_complementarity`,
`solid_material_balance`. A free solid adds one amount coordinate n_s ≥ 0 with
no upper bound and no stationarity row in the optimizer; it is counted as a
complementarity pair, and one endpoint classification adds either its
saturation row or n_s = 0, so the classified declaration is square. A
`Pinned(0)` solid adds no coordinate and reports SI_s only. Solids require fixed
T and P and a neutral reference; PH/PS, critical, envelope, azeotrope, pinned
free-intensive cells, discovery, stability trials, continuation, fit
observations and ASL records refuse them. Solid solutions and kinetics are not
expressible.

## Stationarity, feasibility and stability

For free phase amounts at fixed T/P, the rows are first-order conditions of the
constrained Gibbs minimum. Phase potential equality and zero reaction affinity
follow from stationarity: EqID: `interphase_potential_equality` and
EqID: `reaction_affinity`. The positive-definite reduced Hessian
(EqID: `reduced_hessian_stability`) is a sufficient strict local condition on a
regular constraint manifold, not a necessary condition for every local minimum
or a global certificate.

Pinned amounts or an independent caloric row, with free intensives, define a feasibility specification using
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
dual tolerance and checked by the retained Gibbs-dual identities. A feasibility specification classifies its unscaled rows, which
are its constraints. Report achieved balance, charge, normalization, pressure,
stationarity/affinity and phase-distinctness diagnostics on declared scales. A stationary candidate
is not automatically a stable equilibrium. A requested distinct phase set
that coalesces must report the degeneration explicitly, including a pair whose
identities swap between consecutive states of a followed branch (passage
through a critical point) other than an envelope's requested ln K straddle.

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

An explicit `stability_test(...)` with
`StabilityOptions.evaluate_supplied_trials = true` instead evaluates each
supplied `composition_guess`, without minimizing composition or creating a
finite split. Entries, including repetitions and unavailable points, retain
caller order. It uses the same unit-amount trial, frozen standards, chemical
plane and EqID: `tangent_plane_distance`. Its point is the supplied
full-species vector after native reconstruction, log/exp chart preparation and
final reconstruction; supported entries must be positive and structural zeros
remain zero. The returned composition is the actual energy-evaluation point.
The total input displacement and final reconstruction displacement are separate.

A supplied packing coordinate is replayed through the chart and retains its
pressure defect; without one, the existing Liquid/Vapor density root selects
packing at the initially reconstructed composition. Reconstruction preserves
that packing when recomputing molar volume, without a second root solve. The
returned packing, density and pressure row describe the final point. Packing
outside the domain or nonpositive supported composition yields an unavailable
row with its existing cause, never an invented TPD. Charged-family admission
below applies even when packing is supplied. An explicit charged Vapor support
requires `DensityTemperature` and must match its own lowest stable pressure root
after final reconstruction, within the existing root window and residual criterion,
with positive finite pressure slope. A sole root on the Liquid side is unavailable;
a sole root with no loop evidence retains the density owner's `Unspecified` meaning.
This admission supplies no charged-Vapor discovery or physical vapor-validation claim. Donor slope and estimated
uncertainty remain available when a positive full-support donor exists, but
optimizer-only values are absent. Exhausting the evaluation sequence, a zero
or negative TPD, or a donor estimate establishes no stability conclusion.
[ALG-STABILITY](../science/algorithms.md#alg-stability) owns these operations.

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

Charged trial supports are admitted under every declared bulk-permittivity
form. The existing `DensityTemperature` (Schreckenberg) and
`DipolarPerturbation` (Rueben) forms reach the vacuum limit,
$\epsilon_r\to1$ as $\rho_m\to0$, by their definitions in
[full EOS §6.2 and §6.4](../science/full-eos.md#6-bulk-dielectric-and-debye-huckel),
EqID: `density_correlated_bulk_permittivity` and EqID:
`dipolar_perturbation_bulk_permittivity`; their charged trials count at every
density. `IonSuppression` (including solvent-only, ion-fraction and
ion-specific choices) and `ComponentMixing` are density independent. Their
electrolyte terms are defined for condensed phases, and a free-volume trial
search cannot infer that domain from a Liquid label or mechanical pressure-root
admission alone. Under these forms a charged trial point counts only when it is
a stable density root of its own composition at the test T and P that the
density owner's side classification (#423) places on the liquid side of that
composition's own mechanical loop, the classification declared Liquid phases
already satisfy (#568).

[ALG-STABILITY](../science/algorithms.md#supports-seeds-and-trial-solves) owns
the seed and end-point checks. A vapor-side root, a sole root without a loop
sample and a point away from every stable root are unavailable with the cause
"charged trial outside the condensed domain of density-independent
permittivity"; they are never witnesses or stability evidence. Neutral trials,
vacuum-limit forms, declared solves, continuation, fitting and ASL retain their
existing domain and calculations. No density cutoff, new tolerance,
optimizer-success filter, parameter substitution or model change is introduced.
The condition is necessary for this search, not physical validation or a global
stability certificate, and it establishes condensed density only, not the
correlation's fitted range. The previously evaluated solvent-only gas-density
descent at 2.32 mol/m³ lies on a sole root without a loop and stays refused.

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
extensive Gibbs scaling. Since [#243](https://github.com/tannerpolley/ePC-SAFT/issues/243)
(PR #400) cold declared solving and supplied-result public reuse at that scale
also keep the split.

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

Neutral nonreactive critical equations (EqID: `critical_equations`) are
necessary conditions: a solution is a critical-equation candidate with T, P,
density and composition, not a stability certificate, and its pressure-closed
quantities that divide by P_rho = 0 are unavailable. Envelope cells (EqID:
`envelope_local_condition`) trace a feed with a zero-pinned incipient phase
(ALG-EQUILIBRIUM). Solid phases and kinetics remain separately
justified extensions; this fluid formulation does not imply their support.
