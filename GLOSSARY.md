# ePC-SAFT Engine

## Language

**Species reference**:
The declared pure-liquid, pure-vapor or solvent infinite-dilution chemical-potential standard for one species. `SpeciesReference.PureGas` and `PureLiquid` select the native vapor or liquid branch at finite pressure and are refused where only the other side's root exists; a single root on an isotherm whose ALG-DENSITY scan finds no mechanical loop (a finite search, in practice above the EOS critical temperature) serves either branch; the solvent's default liquid reference retains the near-pure terminal realization.
Source: chemical-potential standard-state invariance in `neutral_reaction_reference_transfer`; Uyan et al. (2015), Eqs. 22–23 for pure/dilute conventions; historical generator e12a5ff8 for pure CO2 vapor; issue #274.

**Independent reference case**:
Analytic identity, measurement, evaluated reference, or source model calculation; source class and uncertainty remain distinct.
Source: CSE terminology and transfer audit §7

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
Source: algorithms

**Convergence record**:
Achieved residual, comparison criterion, termination reason and available conditioning/error diagnostics.
Source: algorithms

**Numerical floor**:
An estimate of attainable accuracy propagated from input/inner-solve error; not automatically a rigorous bound.
Source: algorithms

**Tolerance gated**:
A verification check whose error exceeds its criterion but not a numerical floor derived from reported quantities; a recorded diagnosis pending an attempted root-cause fix, never an end state or a hand-picked looser tolerance.
Source: algorithms

**Fingerprint**:
A recorded identity of parameter, topology or installed wheel contents; not scientific evidence of accuracy.
Source: [packet requirements](data/docs/packet-contract.md) and root `AGENTS.md`

**Typed failure**:
An unavailable result with a specific reason, final only for genuinely singular or undefined mathematics or a state outside the stated physical domain, after an attempted root-cause fix; not an invented physical value.
Source: algorithms, fix before a negative disposition

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

**Component**:
A major responsibility with an implementation or records owner: Engine, Data, Analyses, Applications, Tooling or Documentation.
Source: `ARCHITECTURE.yaml` project map

**Layer**:
An architectural responsibility inside a component; Engine layers have enforceable native include rules.
Source: `ARCHITECTURE.yaml`

**Capability**:
A named behavior such as density evaluation or equilibrium solution; it may use several layers without changing their ownership.
Source: sourceless

**Milestone**:
Bounded work with explicit acceptance requirements, owned by a component and, where applicable, a layer; GitHub uses the established CORE1, C1, P1, … identifiers.
Source: sourceless

**File role**:
Why a file exists: runtime, specification, test, data, evidence, tooling, configuration, generated or external.
Source: `ARCHITECTURE.yaml` project map

**Qualification**:
Reviewed evidence supports a scientific capability within stated inputs, domains and limitations; never implied by inventory or a passing check.
Source: sourceless
