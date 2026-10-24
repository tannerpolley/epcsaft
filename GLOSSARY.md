# ePC-SAFT Engine

## Language

**Solvent activity**:
Dimensionless solvent fugacity divided by its pure-liquid value at the same temperature and pressure, with the same physical parameters and permittivity model. The practical osmotic coefficient is its negative logarithm divided by solvent molar mass times total solute-species molality.
Source: [practical osmotic coefficient](engine/docs/equations.md#practical_osmotic_coefficient), extensive Helmholtz Euler identity and Hamer and Wu (1972).

**Errors-in-variables estimation**:
Parameter estimation that treats measured inputs as unknown true inputs and penalizes their deviations, together with output deviations, using declared measurement uncertainties.
Source: Zavala (2008), thesis chapter 4, §4.2, Eq. (4.2), printed p. 53; [measurement objective](engine/docs/equations.md#errors_in_variables_objective).

**Reconciled input**:
An estimated physical input consistent with the model and the declared measurement-error objective; its correction is measured from the supplied input value.
Source: Zavala (2008), thesis chapter 4, §4.2, Eq. (4.2), printed p. 53; [measurement objective](engine/docs/equations.md#errors_in_variables_objective).

**Species reference**:
The declared pure-liquid, pure-vapor or solvent infinite-dilution chemical-potential standard for one species. `SpeciesReference.PureGas` and `PureLiquid` select the native vapor or liquid branch at finite pressure and are refused where only the other side's root exists; a single root on an isotherm whose ALG-DENSITY scan finds no mechanical loop (a finite search, in practice above the EOS critical temperature) serves either branch; the solvent's default liquid reference uses the exact ion-free solvent state.
Source: chemical-potential standard-state invariance in `neutral_reaction_reference_transfer`; Uyan et al. (2015), Eqs. 22–23 for pure/dilute conventions; historical generator e12a5ff8 for pure CO2 vapor; issue #274.

**Support**:
Species with coordinates in a phase; omitted species have structural zero amount.
Source: equilibrium formulation

**Inventory face**:
Nonnegative amount states with species forced absent by an exactly zero globally nonnegative conserved inventory; state balances restrict to surviving coordinates while the full feed map and original reaction/reference completion remain.
Source: equilibrium formulation, #231

**Admission**:
Checking declared physical inputs and the requested mathematical domain before evaluation.
Source: full EOS §12

**Derivative action**:
A requested directional derivative of a named observable on a stated branch and coordinate chart.
Source: [ALG-ACTIONS](engine/docs/science/algorithms.md#alg-actions)

**Convergence record**:
Achieved residual, comparison criterion, termination reason and available conditioning/error diagnostics.
Source: [ALG-TERMINATION](engine/docs/science/algorithms.md#alg-termination)

**Numerical floor**:
An estimate of attainable accuracy propagated from input/inner-solve error; not automatically a rigorous bound.
Source: [ALG-TERMINATION](engine/docs/science/algorithms.md#alg-termination)

**Tolerance gated**:
A verification check whose error exceeds its criterion but not a numerical floor derived from reported quantities; a recorded diagnosis pending an attempted root-cause fix, never an end state or a hand-picked looser tolerance.
Source: [ALG-TERMINATION](engine/docs/science/algorithms.md#alg-termination)

**Fingerprint**:
A recorded identity of parameter, topology or installed wheel contents; not scientific evidence of accuracy.
Source: [packet requirements](data/docs/packet-contract.md) and root `AGENTS.md`

**Typed failure**:
An unavailable result carrying a specific cause after diagnosis. A finite density search that cannot establish the requested branch after a resolved sampled contradiction retains its cause without asserting physical nonexistence. No physical value is invented.
Source: [ALG-TERMINATION](engine/docs/science/algorithms.md#alg-termination), fix before a negative disposition; [ALG-DENSITY](engine/docs/science/algorithms.md#alg-density).

**Claim boundary**:
Observable, domain, evidence and uncertainty to which a scientific conclusion applies.
Source: CSE scientific method

**Tangent-plane distance**:
Reduced molar energy above the candidate conserved supporting plane on normalized electroneutral composition/volume trials.
Source: Michelsen1982 Part I; Perdomo2025 §3.3; equilibrium formulation

**Numerically resolved descending split**:
Positive-donor inventory-conserving split with negative computed slope and finite Gibbs decrease exceeding propagated uncertainty estimates; numerical verification, not proof.
Source: equilibrium formulation; #223 adopted D9

**Proved instability witness**:
Feasible Gibbs descent with justified total numerical error and feasibility bounds; association floors alone are insufficient.
Source: equilibrium formulation

**Stability certificate**:
A justified bound or theorem for the stated domain; finite search without a witness is insufficient.
Source: equilibrium formulation

**Critical-equation candidate**:
A solved state of the necessary Heidemann–Khalil conditions M u = 0, u'u = 1 and c = 0 with one lowest zero eigenvalue of M; it certifies neither stability nor a unique critical point.
Source: EqID `critical_equations`; Heidemann and Khalil (1980); Michelsen (1980)

**Envelope cell**:
A free feed phase with an incipient phase pinned at zero amount, T and P free and one selected ln T, ln P or ln K condition; its continuation passes a mixture critical point by one ln K straddle.
Source: EqID `envelope_local_condition`; Michelsen (1980), pp. 2–4

**Azeotrope cell**:
A binary bubble/dew cell whose first C − 1 species balances are replaced by equal phase compositions and whose remaining balance is the total amount; the declared feed composition is only the seed.
Source: EqID `azeotrope_condition`; FeOs `PhaseEquilibrium.binary_azeotrope`

**Component**:
A major responsibility with an implementation or records owner: Engine, Data, Analyses, Applications, Tooling or Documentation.
Source: `ARCHITECTURE.yaml` project map

**Layer**:
An architectural responsibility inside a component; Engine layers have enforceable native include rules.
Source: `ARCHITECTURE.yaml`

**Capability**:
A named behavior such as density evaluation or equilibrium solution; it may use several layers without changing their ownership.
Source: `ARCHITECTURE.yaml` header (Python is organised by capability) and `python.capabilities`

**Milestone**:
A numbered GitHub milestone (for example `04 - Equilibrium`) grouping leaf issues with explicit acceptance; its number expresses delivery priority, not a scientific dependency.
Source: `docs/scientific/CONTEXT.md` "Implementation completion order"; ADR 0010

**File role**:
Why a file exists: runtime, specification, test, data, evidence, tooling, configuration, generated or external.
Source: `ARCHITECTURE.yaml` project map

**Qualification**:
Reviewed evidence supports a scientific capability within stated inputs, domains and limitations; never implied by inventory or a passing check.
Source: `docs/scientific/CONTEXT.md` "The project map" (leaf issues record unqualified capability) and "Review rules"

**Ideal mole-fraction activity**:
The dimensionless activity $a_i=x_i$ with standard fraction one and unit activity coefficient.
Source: Nasrifar and Tafazzol (2010), Eq. (11) and section 6; `engine/docs/equations.md` ideal_mole_fraction_reaction.

**Ideal molarity activity**:
The dimensionless activity $a_i=n_i/(V c^\circ)$ at declared fixed solution volume,
with $c^\circ=1$ mol/L; an explicitly declared neutral solvent may instead have unit activity.
Source: IAPWS TGD9-18, Table 2, p.15 (ideal concentration constants); `engine/docs/equations.md` ideal_molarity_reaction.

**Effective-water amount**:
The free water amount plus all declared ionic species amounts, used with preserved free molecular amounts.
Source: Nasrifar and Tafazzol (2010), section 6, p.7627; `engine/docs/equations.md` effective_water_amounts.
