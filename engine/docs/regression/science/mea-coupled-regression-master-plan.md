# MEA Coupled Regression Master Plan

Status: Regression transport ready; reduced tracer blocked on a certified
stable homogeneous state and frozen rank-sufficient coordinates.

This document is the single Regression-owned plan for fitting ePC-SAFT
parameters to aqueous monoethanolamine observations. It coordinates the
application contract owned by MEA-Thermodynamics, thermodynamic state and
sensitivity results owned by Equilibrium, thermodynamic values and exact
parameter partials owned by EOS, and the Ceres workflow owned by Regression.

It does not admit a parameter set, authorize a fit, persist values to the
EOS catalog, or make a predictive claim.

## 1. Governing decision

The first coupled tracer shall use:

1. one admitted CO2 partial-pressure observation; and
2. one eligible liquid-speciation equality or linear aggregate.

CO2 loading is a fixed state input, not the measured response. Heat of
absorption is excluded from the first tracer because no source-complete heat
row or observation definition is admitted.

Under the user-approved low-pressure first-tracer convention, the pressure
prediction is

```text
p_CO2 = f_CO2_liquid
```

for neutral CO2 in the certified homogeneous reacting liquid. This is the
liquid-fugacity equivalent of the measured low-pressure gas partial pressure
with vapor fugacity coefficient fixed to one. It is an explicit application
convention for this reduced engineering tracer, not an experimental
uncertainty, a general vapor-phase approximation, or a predictive claim.
There is no vapor composition, bubble root, or phase-equilibrium solve in the
first tracer. Reactive-bubble capability remains a separate later path for
applications that require nonideal vapor or phase-boundary results.

The old proposed tracer—CO2-loading equality plus heat equality—is rejected as
unexecutable from the admitted source packet. The old historical 12-parameter
block and 435 reserved-observation count are also stale.

## 2. Ownership

| Owner | Owns | Must not own |
|---|---|---|
| MEA-Thermodynamics | species and reactions; equilibrium-constant sources; standard-state convention; feed/state construction; source rows; target eligibility; measurement mapping; parameter selection; bounds, scales, starts, sharing, regularization, partitions, and scientific gates | EOS equations, Equilibrium algorithms, Ceres, or fitted-value persistence |
| EOS | typed parameters and topology; Helmholtz/reference calculations; phase and caloric primitives; exact explicit parameter partials; domain/applicability and artifact identity | reaction equilibrium, Ceres, application datasets, or promotion |
| Equilibrium | source-to-EOS reference transformation; reacting-liquid solves; exact implicit state sensitivities; solver/numerical/physical certificates; separately admitted later reactive-bubble solves | parameter selection, residual weights, regression objectives, or application promotion |
| Regression | source-bound observation validation; parameter transforms and sharing execution; residual/Jacobian assembly; Ceres; rank, conditioning, active-bound and confirmation diagnostics; immutable authority-neutral results | chemistry defaults, a second equilibrium formulation, copied EOS/reaction equations, or EOS catalog mutation |
| Validation campaign | installed-artifact black-box replay and durable cross-package evidence authored by the accountable package task | production algorithms or private source imports |

Organization decision
[`ePC-SAFT/.github#1`](https://github.com/ePC-SAFT/.github/issues/1) admits the
transport. The downstream application constructs the model-bound evaluator
from exact installed EOS and Equilibrium artifacts and supplies the
versioned process-local handle described by Regression issue 15. Regression
does not import or link Equilibrium.

## 3. Current source contract

The current exact Gate-0 application subject is
`MEA-Thermodynamics@269c954230b73bffe19d157137143a52d9c685f6`, tree
`7d58f7b50b2d3e0682a862b92e8f1c998501cf80`, merged through PR 51. Later
exploratory SciPy fits are rejected evidence and do not define this contract.

### 3.1 Chemistry

Canonical species order:

```text
CO2, MEA, H2O, MEAH+, MEACOO-, HCO3-, CO3^2-, H3O+, OH-
```

Canonical charges:

```text
0, 0, 0, +1, -1, -1, -2, +1, -1
```

Canonical reactions:

```text
2 H2O                <=> H3O+ + OH-
CO2 + 2 H2O          <=> HCO3- + H3O+
HCO3- + H2O          <=> CO3^2- + H3O+
MEACOO- + H2O        <=> HCO3- + MEA
MEAH+ + H2O          <=> MEA + H3O+
```

The application source contract is
`data/reference/MEA/manifests/chemical_reaction_source_contract.json`,
SHA-256
`77c6c0a705aa2f5bcc5f9ff90722ffadfe9b2ef0ae3cf41b8a1ac8fecbd286d8`.
Its EOS-basis transform is not yet complete. The source spelling
`CO3--` and runtime spelling `CO3^2-` must be reconciled by one explicit
identity mapping in the installed application contract, not by a Regression
alias.

### 3.2 Pressure and speciation partition

The frozen grouped split is
`data/reference/MEA/manifests/grouped_split_manifest.csv`, SHA-256
`af205ad5968667cf25dc9205d780738035769664a94cc9a421cd3c67148ff804`.

| Observation family | Training | Reserved |
|---|---:|---:|
| CO2 partial-pressure rows | 89 | 167 |
| Speciation states | 58 | 53 |
| Total states | 147 | 220 |
| Eligible speciation scalar records | 208 | 100 on 16 states |
| Maximum eligible pressure plus speciation residuals | 297 | 267 |

The reserved total is 267, not 435. The stale value counted 168 speciation
records whose manifest explicitly marks them ineligible.

Training memberships:

- pressure: Aronu `vle_obs_0034` through `vle_obs_0069`; Hilliard
  `vle_obs_0119` through `vle_obs_0142`; Idris `vle_obs_0162` through
  `vle_obs_0171`; Ma'mun `vle_obs_0246` through `vle_obs_0264`;
- speciation: Böttinger states 030 through 068; Matin states 001 through 019.

Reserved memberships:

- pressure: Aronu 0001–0033 and 0070–0106; Hilliard 0107–0118 and
  0143–0161; Jou 0186–0233; Xu 0265–0272 and 0279–0288;
- speciation: Böttinger 001–029 and Jakobsen 001–024. Only Jakobsen
  009–024 currently carries eligible quantitative validation observations.

The speciation membership file is
`data/reference/MEA/manifests/speciation_target_membership.csv`, SHA-256
`1e267adc62a6572d029caf92cfa39529cd38c9c0475fdf2f66a9546dc92c8f9b`.

### 3.3 Volumetric partition

The separately frozen volumetric campaign contains:

- 153 future training densities: 63 Amundsen, 56 Augusto, and 34 Dhage;
- 78 reserved densities: 40 Amundsen, 28 Augusto, and 10 Dhage.

The 44 derived excess-molar-volume values are diagnostics derived from parent
density rows. They must not be added as simultaneous independent residuals.
The volumetric campaign is a parameter-identification stage, not evidence that
the coupled pressure/speciation path is ready.

## 4. Observation contracts

### 4.1 State inputs

Every state must bind:

```text
row_id
source_id and immutable source hash
partition and group identity
temperature_K
pressure specification and pressure role
unloaded_MEA_mass_fraction
CO2_loading_mol_per_mol_MEA
apparent carbon, MEA, and water totals
ordered species/reaction/reference identities
installed artifact fingerprints
```

Loading, temperature, and unloaded MEA fraction are state specifications. They
are not residuals unless a future source explicitly reports them as uncertain
observations and admits a data-reconciliation contract.

### 4.2 CO2 partial pressure

Each pressure row must add:

```text
observed_p_CO2
unit
measurement_origin:
  direct_partial_pressure |
  calibration_derived_partial_pressure |
  total_pressure_derived |
  model_derived
measured_primitive
derivation and solvent-vapor correction
gas-analysis method and reported vapor composition
uncertainty and covariance when reported
source table/figure/row locator
residual transform, scale, and rationale
```

`model_derived` rows are never fit targets. The merged Gate-0 packet classifies
Hilliard `vle_obs_0137` as
`calibration_derived_partial_pressure` and freezes it as the first tracer row;
its `574 Pa` value comes from the source's calibrated gas composition and
measured total pressure. That classification does not turn the derived value
into a direct partial-pressure measurement or assign experimental uncertainty.
Only `vle_obs_0171` presently carries a numeric pressure uncertainty,
`0.077 kPa`; it is not substituted for the preregistered tracer row.

### 4.3 Liquid speciation

Admitted roles are:

- positive species equality;
- source-defined linear aggregate equality, principally `MEA + MEAH+`; and
- one-sided upper censoring for a reported zero or below-detection result.

Training contains 136 positive equalities, 58 aggregates, and 14 censored
records. Reserved evidence contains 83 positive equalities, 16 aggregates,
and one censored record.

Aggregate coefficients, measurement identity, and covariance must be explicit
to avoid double counting. A reported zero without a numerical detection limit
cannot become `log(0)`, an invented zero equality, or an invented censor bound.
The current 14 training and one reserved zero records are therefore
non-executable until MEA acquires source-backed bounds or demotes them and
recomputes the executable observation counts and hashes. The 297/267 counts
remain maxima until that censor and pressure-origin adjudication finishes; the
grouped split hash changes only if state membership changes.
Balance-inferred, ambiguous, contextual, calibration-derived, and
`target_eligible=no` records do not enter Ceres.

### 4.4 Heat of absorption

Heat is `NOT_ADMITTED_NO_SOURCE_CONTRACT`. There are no qualified rows.
Before heat can join the objective, MEA must freeze:

- row identities and primary-source locators;
- differential, integral, or semi-differential definition;
- calorimetric versus derived measurement role;
- sign convention;
- initial and final loading/reference states;
- constant-pressure or other process path;
- energy unit and per-mole-absorbed-CO2 basis;
- uncertainty, covariance, partition, residual scale, and promotion gate.

EOS must then expose the required caloric primitives and exact temperature
and active-parameter partials, and Equilibrium must expose the corresponding
total sensitivity through the reacting state. This capability is outside
Equilibrium issues 36–38.

## 5. Thermodynamic and derivative contract

For a fixed application state and physical fitted parameters `theta`,
Equilibrium owns the reacting-liquid optimization and its complete active-set
KKT system:

```text
H(w, theta) =
  [gradient_z(L),
   c(z, theta),
   g_active(z, theta)] = 0
```

where `w` contains the reduced amount chart, `log(V)`, multipliers, and the
active-set variables required by the Equilibrium formulation; `c` is the
implemented equality-constraint block and `g_active` is the fixed active-bound
block. The equality block contains the compiled independent material balances;
the ionic amount chart enforces electroneutrality identically and charge is
independently certified after the solve. Reaction and pressure equilibrium
appear once through the implemented stationarity formulation—the fixed-pressure
objective contains the `P*V/(R*T)` term—and are not appended again as duplicate
constraints. On a fixed smooth active set:

```text
H_w * dw/dtheta = -H_theta
dw/dtheta = -H_w^{-1} H_theta
```

Equilibrium then extracts the state-coordinate rows `dz/dtheta`. For any
primitive observation `y(z, theta)`, the exact total derivative consumed by
Regression is:

```text
dy/dtheta = partial(y)/partial(theta)
          + partial(y)/partial(z) * dz/dtheta
```

EOS-parameter sensitivities require exact EOS mixed partials of the
state gradient with respect to every active physical parameter, exact
parameter derivatives of packing/domain functions, and exact derivatives of
the transformed reference vector. Missing terms make the Jacobian unavailable;
neither Equilibrium nor Regression substitutes finite differences or frozen
speciation.

The first-tracer CO2 pressure primitive is the liquid fugacity of neutral CO2.
For the EOS Helmholtz basis `Phi = A/(R T n_ref)` with
`n_ref = 1 mol` and `rho_ref = 1 mol/m^3`,

```text
ln(f_CO2_liquid / (rho_ref R T)) = partial(Phi)/partial(n_CO2)
```

at fixed temperature. Its exact total derivative is

```text
d ln(f_CO2_liquid)/dtheta_j
  = sum_l Phi_(n_CO2,n_l) * dn_l/dtheta_j
  + Phi_(n_CO2,V) * dV/dtheta_j
  + Phi_(n_CO2,theta_j)
```

EOS supplies the gradient, Hessian, and explicit active-parameter
chemical-potential partial. Equilibrium supplies exact amount and volume
sensitivities. This composition uses first-order implicit sensitivities and
EOS second derivatives; it needs no third derivative, density root,
vapor-incidence residual, or bubble solve. A missing, reordered, nonfinite, or
unsupported term makes that observation Jacobian unavailable.

## 6. Equilibrium value and certificate contract

An installed composed evaluator receipt must retain:

- `T` in K and `P` in Pa;
- ordered species IDs, charges, amounts in mol, mole fractions, and volume in
  m3;
- EOS component order, parameter fingerprint, reference identity, and
  artifact identity;
- ordered reaction and balance rows;
- source standard-state, activity-scale, conversion, and reaction-constant
  provenance;
- amount-chart topology and all active lower, upper, trace, and domain bounds;
- independent solver, numerical, physical, local-minimum, search, root,
  predictive, and globality statuses.

The current public homogeneous-liquid result already exposes the core state and
certificate values but does not expose all ordered balance/reaction matrices,
artifact hashes, sensitivity active sets/chart topology, or
search/root-completeness fields above. For the first homogeneous-liquid
tracer, search and root fields are `not_applicable`. Later reactive-bubble
results require their own complete search and root evidence.

Current Equilibrium value gates independently recompute:

- balance and charge closure at `1e-9` absolute;
- pressure closure at `1e-8` relative;
- reaction affinity at `1e-7`;
- KKT stationarity and complementarity at `1e-7`;
- positivity, packing, EOS domain, and reduced-Hessian local-minimum
  status.

Sensitivity availability is separate from value validity. A Jacobian fails
closed for an uncertified primal state, singular system, infinity-norm
condition estimate above `1e6`, bound activity within the `1e-7` margin,
trace/structural-face activity, topology change, missing EOS mixed
partial, or missing transformed-reference derivative.

Equilibrium issues 36 and 38 are complete, and Regression issue 15 is complete
through PR 27. The public sensitivity result includes exact active
EOS-parameter amount and volume columns with conditioning and failure evidence,
and Regression owns the typed process-local positive-observation transport.
Issue 37 is closed as not planned; reactive-bubble work is not a tracer gate.

The retained issue-38 Hilliard/Böttinger state is negative evidence, not an
admissible value sentinel. It is `FEASIBLE_ONLY`: solver, numerical, physical,
and EOS-domain gates pass, while the local-minimum gate fails. An independent
reduced-Hessian calculation confirms negative curvature. Exact balance
retraction left no predeclared seed that passed the descent/admission screen.
The local-minimum gate remains unchanged.

Before issue 16 executes, the downstream application must bind a fresh
immutable installed evaluator handle that composes the public EOS and
Equilibrium contracts into ordered primitive values and complete exact
parameter columns for a certified stable state. The handle binds artifact
identity, units, state and parameter ordering, reference identity, chart
topology, and EOS fingerprints. Regression must not consume a private
underscored API.

## 7. Regression residual and Jacobian assembly

Regression reuses the existing ordered `[N fitted | Q lifted]` Ceres core,
result family, native module, and native target. Equilibrium eliminates the
reacting state, so the coupled MEA path has `Q = 0`.

The application chooses one admitted transform per observation:

```text
positive equality: identity or log transform
linear aggregate:  a^T y
upper/lower censor: frozen differentiable one-sided rule
```

The log base, scale, family normalization, and any uncertainty interpretation
remain MEA inputs. Equal weights of 1.0 and the historical regularization
weight 0.003 are not source-backed defaults and shall not be inferred.

For the evaluator's ordered primitive values `y`, exact physical-parameter
Jacobian `Y_theta`, sharing map `S`, and affine solver transform `D`,
Regression assembles:

```text
J = W * h_y * Y_theta * S * D
```

Every requested fitted column must be exact, finite, ordered, and complete.
Censored rows that are locally inactive contribute no local rank.

## 8. Parameter-identification stages

The old 12-parameter pressure/speciation block is not the production plan.

| Stage | Active coordinates | Current role |
|---|---|---|
| Shared-cation analog | `MEAH+::sigma`, `MEAH+::epsilon_over_k` | identify a common cation contribution from admitted direct analog densities after immutable counterion parameters/covariance or declared joint regularization exist |
| Carbamate correction | `MEACOO-::sigma` | smallest correction from electroneutral reactive-liquid evidence after an immutable transfer prior and scale exist; isolated-ion states and directional Maiti context are not residuals |
| Reactive refinement | the three coordinates above | frozen objective is Amundsen density plus eligible speciation; pressure is a rejection guard unless a pre-fit amendment freezes metrology, scaling, sensitivities, and the revised objective |

Initial safeguard contract:

| Parameter | Historical start | Bounds |
|---|---:|---:|
| `MEAH+::sigma` | 3.48508556586 angstrom | 2.0–5.8 angstrom |
| `MEAH+::epsilon_over_k` | 232.687201645 K | 50–950 K |
| `MEACOO-::sigma` | 3.53543525721 angstrom | 2.0–5.8 angstrom |

These starts are provisional seeds, not scientific authority. Affine coordinate
scales remain unfrozen. Segment counts, Born diameters,
`MEACOO-::epsilon_over_k`, reaction constants, trace-ion parameters, and all
binary interactions remain fixed initially.

`k_ij(MEAH+,MEACOO-)` may replace, not augment, the carbamate correction only
after a preregistered sensitivity amendment. A quadratic transfer prior for
`MEACOO-::sigma` remains blocked on an immutable uncertainty scale and generic
regularization support.

Frozen density scaling is application input:

- Augusto uses the greater of the row standard deviation and reported relative
  uncertainty times density;
- Dhage uses 0.0005 g/cm3;
- Amundsen uses 0.0005 g/cm3 for unloaded rows and 0.002 g/cm3 for loaded rows;
- each source/salt family is normalized by its admitted row count.

## 9. Executable sequence

### Gate 0 — application contract

MEA must:

- complete the EOS-basis reaction transformation and carbonate identity;
- enforce `target_eligible`;
- propagate eligible aggregates into the executable payload;
- acquire source-backed detection bounds for censored rows or demote them, then
  recompute executable counts and observation hashes;
- classify pressure metrology;
- freeze source-backed row pressure or an explicitly bounded pressure
  convention and a 2-oxazolidone inclusion/exclusion rule before selecting the
  sentinel or tracer state;
- freeze affine and residual scales, starts, bounds, sharing, and exact tracer
  rows;
- remove heat from the active tracer or admit a complete calorimetric packet.

### Gate 1 — certified installed reacting-liquid value sentinel

Equilibrium issue 38 retained one source-complete fixed-`T,P` aqueous MEA
liquid state using the nine-species/five-reaction bundle and an exact installed
EOS artifact. The retained state is `FEASIBLE_ONLY`, not a certified local
minimum, so it is falsification evidence and cannot be consumed by Ceres. A
distinct Equilibrium continuation/multistart design and immutable artifact must
return a certified stable state before this gate passes.

### Gate 2 — exact active-parameter sensitivities

EOS supplies the exact explicit partials for only the preregistered coordinate
subset. Equilibrium issue 36 supplies the implicit columns and source-reference
transform. Value-valid but Jacobian-unavailable states remain unusable for
Ceres, and derivative completeness does not override a failed local-minimum
gate.

### Gate 3 — homogeneous-liquid observable composition

The installed downstream evaluator composes the certified reacting-liquid
state, neutral CO2 liquid fugacity, eligible speciation primitive, and every
exact selected parameter column. The two first-tracer rows are independent
source observations; they are not relabeled a paired experiment. A retained
installed-artifact receipt must prove the values, Jacobian columns,
fingerprints, status propagation, and failure behavior.

### Gate 4 — cross-package evaluator transport

Organization decision
[`ePC-SAFT/.github#1`](https://github.com/ePC-SAFT/.github/issues/1) admits the
inverted process-local transport. Regression issue 15 and PR 27 bind the
downstream-composed, model-bound, versioned evaluator through that transport.
The first executable slice supports positive identity and log equalities,
value-only and value-plus-Jacobian responses, deterministic ordering, complete
fingerprints, per-row failures, and no Regression-to-Equilibrium dependency.
Linear aggregates and censoring remain frozen designs until a real
source-bound row requires their runtime surface.

### Gate 5 — reduced coupled tracer

Regression issue 16 fits one or two predeclared coordinates against exactly:

- Hilliard `vle_obs_0137`, observed `p_CO2 = 574 Pa`, modeled as neutral
  `f_CO2_liquid` at `T = 313.15 K`, total pressure `7326.7 Pa`, unloaded MEA
  mass fraction `0.30`, and loading `0.466 mol CO2/mol MEA`; and
- Böttinger `cheq_canon_00194`, observed `x_MEACOO- = 0.0502` at the same
  temperature, unloaded MEA mass fraction, and loading.

The Böttinger source reports no pressure. Its evaluation pressure is the
application-declared fixed state anchor, not a Böttinger measurement. Both
rows are positive equalities; no loading residual, heat residual, aggregate,
censor, bubble calculation, or vapor composition is present.

The exact merged MEA subject binds the row payloads through:

- `data/reference/MEA/observations/vapor_liquid_equilibrium/Canonical_VLE_Observations.csv`,
  SHA-256
  `9e7d9ba5fead8bfa83a311dad341e3e2e8df1806d5249642a23562e99a72cb73`,
  and `data/reference/MEA/manifests/pco2_metrology_manifest.csv`, SHA-256
  `0d14803873a60534ec5d7df382cfbd0ae03e4aaeba68bb5d54be7e4def8397cc`;
- `data/reference/MEA/observations/liquid_speciation/Canonical_Combined_ChEq.csv`,
  SHA-256
  `8c07df9efd1c1ecbd775ccdd42791e0cef1880b3837e5749a60d2142aa85809e`,
  and `data/reference/MEA/manifests/speciation_target_membership.csv`,
  SHA-256
  `a89a3f0373a86813482158f180939cf57f74be038cd59f244dfadcb689923190`.

Regression consumes the installed application contract or immutable campaign
copy of these rows; it does not read a sibling MEA source checkout at runtime.

The two-row sensitivity matrix must have rank `N` at every declared start
before optimization. If no source/physics-selected subset with `N <= 2` has
that rank, the tracer remains blocked; rows, weights, priors, or parameters are
not changed after observing the result.

### Gate 6 — staged identification and reactive refinement

Run the shared-cation analog and carbamate-correction stages only through exact
installed derivatives and their frozen source contracts. Then fit the
three-coordinate reactive block against the preregistered training
observations. Pressure may enter the objective only after metrology and scaling
are frozen; otherwise it remains a coupled rejection guard.

The full staged campaign uses 32 deterministic scrambled-Sobol starts with seed
`390035`, plus the provisional historical seed. Materially distinct best
basins, incomplete rank without profile-supported identifiability, material
start/prior/bound dependence, or a safeguard-bound solution rejects promotion.

### Gate 7 — reserved evaluation

Evaluate without refitting:

- the 220 reserved pressure/speciation states, reporting 37 context-only
  speciation states separately from the 16 quantitatively eligible states; and
- the 78 reserved volumetric rows.

Do not combine derived excess volumes with parent density residuals. Preserve
every failed, skipped, inapplicable, and evaluated row.

Any failed eligible training state rejects the fit. Every failed reserved state
counts as a failed prediction and cannot be relabeled skipped or inapplicable
after execution.

### Gate 8 — later heat expansion

Only after the source, caloric EOS, Equilibrium sensitivity, and
observation-scale gates in section 4.4 are complete may heat join the same
parameter vector and Ceres objective. Its addition requires a new rank and
conditioning preflight; it does not reinterpret the earlier tracer.

## 10. Performance and execution shape

One evaluator batch groups observations only when their complete
specifications and equilibrium operations are identical:

1. perform one top-level evaluation per logical state and operation;
2. reuse its certified values for every attached pressure, speciation,
   aggregate, censor, density, or later heat observation;
3. request a residual-only response, without returned fitted-parameter
   columns, for Ceres residual-only calls;
4. request all exact fitted columns once for Jacobian calls;
5. return deterministic row order and one status for every input.

A fixed-`T,P` liquid, later reactive-bubble boundary, and later calorimetric
path are different operations even when temperature, loading, and composition
match. A
residual-only response may omit returned fitted-parameter columns, but
Equilibrium still consumes the exact state derivatives and EOS tensors
required internally by Ipopt and certification; it cannot replace them with a
value-only phase endpoint.

No nested solver parallelism is allowed. Bounded concurrency requires an
explicit installed reentrancy capability. Per-fit budgets and cancellation
must be task-owned. The historical IDAES/eNRTL workflow is architectural
corroboration only; it is not copied, imported, or used as a runtime
dependency.

## 11. Result and acceptance semantics

Every result separately reports:

- upstream Equilibrium solver, numerical, physical, local-minimum, search,
  root-completeness, and globality status;
- Ceres termination and Regression numerical validity;
- full Jacobian rank, conditioning, active bounds, KKT/bound diagnostics, and
  primary/confirmation-start agreement;
- pressure, each speciation species/aggregate/censor, density, and later heat
  residual summaries separately;
- exact input, evaluated, skipped, failed, training, and reserved accounting;
- workflow validity, scientific comparison, and predictive status.

Optimization success never converts a rejected or uncertified state into a
valid row. A density improvement accompanied by degraded reserved pressure or
speciation fails. Current legacy median-error gates are application policies,
not experimental-uncertainty limits, and require explicit preregistration
before reuse.

A promotion candidate additionally retains one complete successor receipt:
exact EOS, Equilibrium, and Regression versions, commits/trees, installed
wheel or ABI hashes, and capability fingerprints; ordered parameter identities,
units, values, bounds, sources, transforms, and fitted domain; observation and
split hashes; reaction/reference identities; all multistart diagnostics;
Equilibrium certificates; complete training/reserved residual records; and
final bundle fingerprint.

No result establishes global parameter uniqueness, uncertainty, prediction,
EOS-catalog authority, or promotion by itself.

## 12. Current readiness

| Dependency | Status | Smallest missing work |
|---|---|---|
| Regression multi-parameter Ceres core | ready for admitted exact evaluators | none |
| Generic positive-observation equations | implemented in PR 27 | none |
| MEA chemistry/source transform | Gate-0 contract frozen | retain exact merged application subject and source hashes |
| MEA pressure metrology | first row selected | Hilliard `vle_obs_0137`; no model-derived row |
| MEA speciation payload | first row selected | Böttinger `cheq_canon_00194`; direct positive equality |
| MEA heat | not admitted | source-complete packet plus caloric value/sensitivity contract |
| Equilibrium reacting-liquid values and sensitivities | exact transport complete; retained state `FEASIBLE_ONLY` | certified stable homogeneous state from a distinct continuation/multistart design |
| Equilibrium reactive bubble | issue 37 closed as not planned | not a tracer gate |
| Installed MEA campaign | frozen state retained as negative evidence | re-freeze exact artifacts after a stable state is available |
| Regression downstream transport | implemented | PR 27 / issue 15 complete |
| Reduced coupled fit | blocked | certified stable state plus frozen rank-sufficient `N <= 2` coordinates |
| Full mixed campaign | not ready | successful tracer, frozen three-coordinate contract, and rank-sufficient admitted rows |

The next Regression-side work that does not depend on Equilibrium is the
generic strict dataset/control and preparation/preflight work in issues 21 and
20. It must remain chemistry-free and must not fabricate missing MEA
observations or derivatives.

## 13. Current issue reconciliation

- Regression issue 16 names the Hilliard pressure and Böttinger speciation
  rows, uses the liquid-fugacity-equivalent convention, excludes heat/loading
  and bubble work, retains the corrected reserved maximum of 267, and limits
  the tracer to `N <= 2`.
- Regression issue 15 and PR 27 completed the exact positive-observation
  transport; it must still not be bypassed by a Python callback or source
  checkout.
- MEA PR 51 owns the exact first-tracer source/state preparation.
- Equilibrium issue 36 is complete for exact active-parameter state
  sensitivities. Issue 38 is closed with retained `FEASIBLE_ONLY`
  negative-curvature evidence. Issue 37 is closed as not planned and remains
  outside the tracer.
- Resume only after an installed Equilibrium artifact certifies a stable
  homogeneous state and the application freezes a rank-sufficient `N <= 2`
  coordinate block. Re-freeze all artifact identities on resume.

This document controls the coordinated sequence. Application manifests,
installed artifact receipts, and owner-specific implementation contracts
remain authoritative for their own content.
