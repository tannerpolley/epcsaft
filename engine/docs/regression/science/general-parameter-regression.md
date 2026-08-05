# General ePC-SAFT Parameter Regression

Status: general nonpolar/electrolyte families fit-ready for their typed
observation contracts; descriptor-driven fixed-topology association block
executable; polar regression excluded

Date: 2026-07-30

## Active installed EOS binding

The active frontend dependency is EOS `0.2.0.dev0`, merged commit
`d0c2575908323db6c487d4e541cf497c93c89d9e`, tree
`e29386139e7ba81115b0f454c79e0c869a262d37`, wheel SHA-256
`bc7e637de084330ebded4ddfd52e02bc1ce5451221128692972ebba8856d098e`, and
installed-header SHA-256
`2adb85c3dc0502384b5ae5fce69cc46f8874415ca7e6fad85b1a5d7de3a8ea53`.
The installed static-library SHA-256 is
`e4e6c58aa14bf285a9cc347eb6b8b34fa8a7ee563576522223a7aa5e67930198`.
The previous frontend subject's final EOS JUnit receipt artifact is
`provider-python-frontend-v0.2/14fa374/provider-tests.xml`,
SHA-256
`56d1d3e9fcf47bb700df4223fa2d9a20444dc97aa22b5c39525d446a296ba3cf`.
That receipt remains historical evidence for its exact prior subject and is
not rebound to the active ordinary-sigma artifact. Artifact names are relative
to an executor-supplied immutable artifact root. Campaign-specific EOS
subjects below remain historical evidence for the named reference runs.

## Objective

Regression shall support source-bound fitting of any admitted continuous
ePC-SAFT parameter for any EOS-supported component set when the caller
supplies an informative dataset and a complete numerical contract.

Published methane, ethane, Figiel, and later literature workflows are
reference campaigns. They provide positive or negative evidence for an exact
parameter-family domain but do not restrict the runtime to those components,
mixtures, or papers. One campaign does not prove scientific accuracy for every
mixture.

Every literature-related result follows
`docs/science/literature-reproduction-contract.md`. Regression reproduces a
complete declared numerical contract; it does not infer missing author rows,
objectives, weights, starts, bounds, topology, or solver choices from a
published tuple. Exact author-method replay, source-faithful reconstruction,
published-tuple property replay, and modern refitting are distinct evidence
classes.

This design supersedes the earlier assumption that every paper-specific
campaign must remain a separate regression surface. It does not claim that
every EOS parameter is fit-ready today.

## Ownership boundary

EOS remains the only owner of ePC-SAFT equations, resolved parameter
records, thermodynamic states, density closure, and nonlinear derivatives.
Regression owns:

- source-bound observation contracts;
- active parameter selection, bounds, scales, transforms, starts, and sharing;
- residual construction and Ceres execution;
- exact residual-Jacobian assembly from EOS derivatives;
- rank, conditioning, active-bound, confirmation, and row-accounting evidence;
- immutable fit results and fitted-versus-reference comparisons.

Regression must not accept arbitrary executable residual plugins, copied EOS
equations, a numerical production derivative backend, or a mutable parameter
registry. The supported vocabulary for any fit is a closed, versioned typed
schema. EOS may advertise a versioned superset. Regression accepts only
descriptor members its installed version understands, reports the remainder
as unsupported, and rejects any request that references an unknown member
rather than registering it dynamically.

## Capability meanings

The package shall distinguish these states for every parameter family:

1. `REPRESENTED`: EOS can resolve and evaluate a model containing the
   parameter.
2. `DERIVATIVE_READY`: an installed EOS artifact exposes the exact active
   derivative required by at least one observation contract.
3. `FIT_READY`: Regression can construct and solve the corresponding
   source-bound problem for caller-supplied components and rows.
4. `REFERENCE_VALIDATED`: at least one source-backed installed-artifact
   campaign has passed its declared gates and retained producer, reviewer,
   artifact, command, and evidence identities.
5. `ADMITTED`: an accepted promotion receipt binds the exact capability,
   producer and distinct reviewer, source, installed artifacts, evidence
   subject, and user approval. Documentation, a merge, or a converged fit
   cannot set this state.

An example dataset may advance validation evidence without narrowing the
component domain of a `FIT_READY` family.

## Typed parameter identity

An active coordinate is identified by one closed identity variant, not by a
paper or array position:

```text
ParameterCoordinate
  family: ParameterFamily
  identity:
    ComponentParameter(component_id)
    | PairParameter(component_id_a, component_id_b)
    | ModelParameter(model_parameter_id)
    | AssociationParameter(quantity, site_pair_endpoints...)
  capability_id | None                   # None for descriptor-bound slots
  provider_parameter_fingerprint
  provider_topology_fingerprint
  provider_artifact_fingerprint
  transform: Affine(origin, scale)       # v1; p = origin + scale*z
  lower_bound
  upper_bound
```

Each family permits exactly one identity variant. Pair, component, and
association-site endpoints use canonical lexical order; duplicate canonical
identities are rejected. One association identity may contain several site
pairs to declare intentional many-to-one parameter sharing. EOS's immutable
descriptor supplies the finite slot list, canonical physical units, component
and site order, multiplicities, allowed bonds, physical domains, tensor layout,
and exact mapping to derivative columns. Model, topology, parameter, artifact,
identity, order, unit, or bound mismatches fail before optimization.

Bounds and affine scales are mandatory source or workflow inputs; the runtime
supplies no chemistry defaults. `RegressionProblem.start_vectors` is an
ordered collection of complete physical-coordinate vectors. Every vector has
length `N` and uses the exact `parameters` order. The first is primary and all
remaining vectors are confirmations. Per-coordinate start lists and Cartesian
products are forbidden. Multiple observations may share one coordinate, and
multiple coordinates may be fitted together. A future non-affine transform
requires a new closed transform kind and its exact `dp/dz`; v1 must not infer
one from parameter sign or bounds.

Discrete model choices, charge number, component identity, site topology, and
electrolyte formulation are specifications, not continuous fit coordinates.

### Complete EOS-record inventory

The current EOS record vocabulary is the starting inventory, not proof of
derivative readiness:

| EOS record family | Regression treatment |
|---|---|
| `molar_mass` | fixed measured input; not an adjustable ePC-SAFT coordinate |
| `charge_number` | fixed discrete identity; not a continuous coordinate |
| `electrolyte_formulation` | discrete model selection; not fitted |
| `segment_count`, `segment_diameter`, `dispersion_energy_over_k` | continuous component coordinates |
| `relative_permittivity`, `born_diameter`, `solvation_factor` | continuous component coordinates when an exact observation derivative is advertised |
| `schreckenberg_dielectric_volume`, `schreckenberg_dielectric_temperature` | continuous component/correlation coordinates when exact derivatives are advertised |
| `zuber_ion_suppression_coefficient` | continuous component coordinate when exact derivatives are advertised |
| `rueben_dipole_scaling`, `rueben_polarizability_scaling`, `rueben_correlation_integral_parameter` | EOS-represented polar-model inputs explicitly excluded from the Regression roadmap by user decision; no Regression capability is planned |
| `k_ij`, `l_ij` | continuous unordered component-pair coordinates |
| `association_energy_over_k`, `association_volume` | component/site-pair coordinates selected from the installed EOS fixed-topology ordinary-`sigma_ij^3*kappa` descriptor; topology is immutable EOS input, not a fitted coordinate |
| `ion_fraction_suppression_coefficient`, `ionic_region_relative_permittivity` | continuous model coordinates when exact derivatives are advertised |
| correlation terms for `segment_diameter`, `relative_permittivity`, and `solvation_factor` | individually named coefficient coordinates; no whole-correlation opaque fit |

New EOS record families require an explicit schema revision in both
packages. Regression must never infer fit support merely because a record can
be resolved.

## Definitive parameter-family and data-requirement matrix

`FIT_READY` below means that caller-supplied rows can execute today when the
installed EOS advertises the named exact capability. It does not mean
that every dataset is informative, that a reference value will be reproduced,
or that the result is predictive. Every problem still needs declared units,
source identity, partitions, residual scales, parameter bounds/scales/starts,
complete row accounting, full projected parameter rank, acceptable
conditioning, non-bound diagnostics, and confirmation-start agreement.

| Parameter family | Identity | Minimum source-bound data contract | Exact executable contract | Current state and MEA role |
|---|---|---|---|---|
| `segment_count` (`m`) | component | Multiple pure rows with `T`, observed saturation pressure and liquid density, molar mass, phase-volume bounds/starts, partitions, and scales | EOS `(n,V,m)` value/gradient/Hessian; Regression lifted liquid/vapor volumes with pressure, chemical-potential, and density residuals | `FIT_READY`, one family at a time. Standalone recovery is optional for MEA unless `m` is selected in MEA's application-owned parameter block. |
| `segment_diameter` (`sigma`) | component | Same pure saturation pressure/liquid-density contract, over a range that gives nonzero independent sensitivity | EOS `(n,V,sigma)` value/gradient/Hessian; same lifted Ceres owner | `FIT_READY` for a constant coordinate. Named temperature-correlation coefficients are separately `NOT_READY`. Optional for MEA unless selected. |
| `dispersion_energy_over_k` (`epsilon/k`) | component | Same pure saturation pressure/liquid-density contract with enough vapor-pressure sensitivity | EOS `(n,V,epsilon/k)` value/gradient/Hessian; same lifted Ceres owner | `FIT_READY`, one family at a time. Optional for MEA unless selected. |
| `k_ij` | unordered component pair | One supported domain: fixed measured `T,P,x,y` VLE rows; fixed-composition liquid-mixture-density rows with ordered `T,P,x` and a caller-supplied mixture molar mass; source-bound MIAC rows with formula molality; or single-ion solvation-Gibbs targets. Rows must vary enough to identify the selected pair. | Exact EOS Hessian for lifted neutral VLE or liquid-mixture-density closure, or exact first derivative for the admitted aqueous/solvation direct observable; one typed pair coordinate | `FIT_READY` in the advertised neutral-VLE, fixed-composition liquid-density, aqueous-MIAC, and organic-ion-solvation domains. MEA needs exact pair derivatives only if a preregistered amendment selects a pair coordinate; the current staged three-coordinate block selects none. |
| `l_ij` | unordered component pair | Fixed measured `T,P,x,y` VLE rows or fixed-composition liquid-mixture-density rows sensitive to cross diameter, with explicit source sign convention and caller-supplied mixture molar mass | Exact EOS `(n_1,...,n_C,V,l_ij)` Hessian for the installed neutral phase capability and its lifted pressure-density closure | `FIT_READY` for an installed neutral phase capability that advertises the row's ordered component set; no association or electrolyte domain is implied. |
| `born_diameter` | ion component | One or more source-defined single-ion solvation-Gibbs targets with exact x-process convention, state, component order, and numerical scale | Exact EOS solvation-Gibbs value/first derivative for the active ion; direct-observable Ceres row | `FIT_READY`; five Figiel ions are reference evidence. It is parallel parameter groundwork, not coupled-MEA readiness. |
| `solvation_factor` | component | Source-bound MIAC rows at declared `T,P` and formula-unit molality, with solvent/ion identities and scale | Exact EOS MIAC value/first derivative for one active factor; direct-observable Ceres rows | `FIT_READY` for the advertised constant factor. A temperature correlation is separately `NOT_READY`; not on the MEA critical path unless selected. |
| `relative_permittivity` | solvent component | Source-bound single-ion solvation-Gibbs targets with fixed other model inputs, exact solvent identity, x-process convention, state, and scale | Exact EOS solvation-Gibbs value/first derivative for the active solvent permittivity; direct-observable Ceres rows | `FIT_READY`. Corrected Validation subject `e4cb7af` gives five independent rank-1 water fits returning `78.0899937514462` through `78.08999375166104` versus fixed `78.09`. This is implementation evidence, not a paper-fitted target or MEA prerequisite. |
| `ion_fraction_suppression_coefficient` | model | Salt-free-normalized relative-permittivity observations spanning enough total-ion mole fraction to identify one coefficient | Exact EOS relative-permittivity ratio and first derivative; direct-observable Ceres rows | `FIT_READY`; 36 digitized Figiel water/methanol rows are reference evidence. Optional standalone recovery. |
| `ionic_region_relative_permittivity` | model | Source-defined single-ion solvation-Gibbs targets with fixed Born diameters and other inputs | Exact EOS SSM+DS solvation-Gibbs value/first derivative; direct-observable Ceres rows | `FIT_READY`; five independent Figiel fits recover fixed `8` within `2.2e-9`. Parallel groundwork, not coupled-MEA readiness. |
| `association_energy_over_k` | one or more canonical component/site pairs in an EOS-advertised fixed topology | Association-sensitive pure rows or fixed-composition liquid-mixture-density rows sufficient for the selected subset and sharing contract, with exact component order, identities, units, topology, bounds, scales, and starts | One variable-width EOS request returns `(n_1,...,n_C,V,selected slots...)` value/gradient/full Hessian; Regression maps `S` slots to `N` fitted coordinates and sums shared derivative columns | `FIT_READY_FIXED_TOPOLOGY`; energy-only, joint, subset, sharing, pure 2B/non-2B, and anonymous associating-mixture mechanics are routine-tested. Case-specific identifiability and acceptance remain Validation evidence. |
| `association_volume` | one or more canonical component/site pairs in the same immutable descriptor | The same source-complete requirement applies; volume slots may be selected independently or jointly with ordinary-component and energy slots | The same descriptor/request/result path; no scalar association callback, topology-name parser, fixed-component-count branch, or paper-specific runtime | `FIT_READY_FIXED_TOPOLOGY`; no literature-specific tuple or topology is accepted or packaged. |
| `k_hb_ij` | source-defined cross-association combining-rule coordinate | A source-defined combining rule and sign convention plus composition/temperature-varying association-sensitive observations that separate cross from pure association | New EOS record/transform identity, versioned fingerprint, exact chain-rule derivative, and an association-endpoint Regression identity are required | `NOT_READY`. Ascani's fixed `0.026` is provenance, not a recovery dataset. Do not alias it to resolved association energy or invent zero defaults. Not the next MEA investment. |
| `schreckenberg_dielectric_volume`, `schreckenberg_dielectric_temperature` | component/correlation coordinates | Multi-temperature electrolyte relative-permittivity or other source observables that independently identify the selected coefficient | New exact EOS active-parameter callback and Regression correlation identity/surface | `REPRESENTED_NOT_DERIVATIVE_READY`; no retained rank-sufficient recovery series. Optional, not on the MEA critical path. |
| `zuber_ion_suppression_coefficient` | ion component | Relative-permittivity, MIAC, osmotic/activity, or solvation observations spanning ion fraction/molality and preferably temperature | New exact EOS active-parameter callback and matching typed Regression observation contract | `REPRESENTED_NOT_DERIVATIVE_READY`; the retained one-row osmotic oracle is insufficient. Optional. |
| Named correlation coefficients for `segment_diameter`, `relative_permittivity`, or `solvation_factor` | component, correlation family, coefficient kind, and term index | Multi-temperature rows covering a range sufficient to distinguish constant, amplitude, and exponent coefficients | New closed correlation-coordinate descriptor, exact EOS derivative including the correlation chain rule, and Regression identity | `REPRESENTED_NOT_DERIVATIVE_READY`; current constant-parameter evidence does not validate coefficient recovery. |
| Temperature-dependent `k_ij(T)` coefficients | pair plus named correlation coefficient | Multi-temperature mixture observations and the exact source correlation form, reference temperature, units, and sign convention | New EOS pair-correlation record/capability and exact derivative; Regression correlation identity | `NOT_REPRESENTED_NOT_READY`. Ascani provides a source law but the retained bounded case is only 298.15 K. |
| Polar parameter families | source-specific | none in this roadmap | none in Regression | `EXCLUDED_FROM_REGRESSION_ROADMAP_BY_USER_DECISION`. EOS may still represent and evaluate polar physics for other owners. |

`molar_mass`, `charge_number`, component identity, association-site topology,
association scheme, and electrolyte formulation remain fixed measured or
discrete model inputs. They are not continuous Regression coordinates.
Experimental uncertainty may define a residual weight only when the source
reports it and the submitted contract declares that interpretation.

## Downstream decision: MEA before optional family recovery

Standalone parameter recovery and coupled reactive regression are different
deliverables. The Figiel/Born, pure-association, dielectric, and correlation
campaigns can demonstrate individual EOS/Regression contracts, but they
do not establish the coupled state sensitivities needed by
MEA-Thermodynamics.

The complete cross-owner sequence is controlled by
`docs/science/mea-coupled-regression-master-plan.md`. Its corrections are:

- loading is a fixed state input, not the pressure residual;
- the first tracer uses one admitted CO2 partial-pressure observation and one
  eligible speciation equality or aggregate with `N <= 2`;
- heat remains excluded because no source-complete heat contract exists;
- the application-selected production candidate is a staged three-coordinate
  ion block, not the stale historical 12-parameter block;
- the frozen pressure/speciation partitions contain 297 maximum eligible
  training residuals and 267, not 435, reserved residuals; and
- the first tracer uses the explicitly application-approved low-pressure
  `p_CO2 = f_CO2_liquid` convention; a reactive bubble is a later independent
  capability, not a prerequisite.

Regression still reuses its one native Ceres engine/result/target and consumes
only complete exact downstream-composed value and Jacobian blocks. It does not
copy reaction/EOS equations, run a second equilibrium formulation, or
finite-difference a black-box solve. Issue 15 and PR 27 completed the
model-bound transport. Runtime is now blocked on a certified stable installed
homogeneous state and a frozen rank-sufficient `N <= 2` coordinate block. EOS
exact parameter partials and Equilibrium active-parameter state sensitivities
are implemented inputs to that composition; reactive-bubble issue 37 is closed
as not planned and is not a tracer gate.

### `k_hb_ij` rule

`k_hb_ij` must not be inferred from a resolved cross-association energy and
silently fitted under another name. A fit may use it only when EOS exposes
the source-defined combining rule, its component/site identity, a versioned
transform identifier and fingerprint, and its exact total derivative including
the transform chain rule. Otherwise the caller may fit the resolved
`association_energy_over_k` coordinate, and the result must retain that
different meaning.

## Source-bound observations

Dataset-level provenance is stored once:

```text
SourceDescriptor
  source_id
  citation and durable locator
  source artifact SHA-256
  canonical transformed-dataset SHA-256
  transformation record
  units and composition/property bases
  license or use basis
  tolerance and residual-scale rationale
```

Every observation row references that descriptor and contains:

```text
row_id
observation_type
component_ids and exact order
T and P specifications
composition and basis
phase or reference-state convention
observed value and unit
residual scale and its rationale
partition: TRAINING | HELD_OUT | STRESS
source_id and row locator
```

For aqueous mean-ionic-activity fits, one scalar `k_ij` is active per problem.
Each row carries the ordered EOS model identity and the explicit fixed
`(k_water,cation, k_water,anion, k_cation,anion)` context. The active pair
replaces exactly one entry before the installed EOS batch callback is
evaluated. The residual and exact Jacobian are

```text
r_i = (1 - gamma_model,i / gamma_observed,i) / s_i
dr_i/dz = -(gamma_model,i / gamma_observed,i)
           * d ln(gamma_model,i)/d k_ij * (d k_ij/dz) / s_i.
```

The fixed context is a required workflow input, not a catalog default inferred
by Regression. Water--cation, water--anion, and cation--anion are separate
closed EOS capabilities; fitting one does not silently activate the other
two.

For an advertised organic ion-solvation capability, every row carries the
ordered `(solvent, cation, anion)` model, active ion, active unordered pair,
and explicit fixed
`(k_solvent,cation, k_solvent,anion, k_cation,anion)` context. The active pair
replaces exactly one entry. With solver coordinate `z`, physical pair
parameter `k = k_origin + k_scale z`, and solvation-energy scale `s_G`,

```text
r_G = (G_solv(k) - G_solv,observed) / s_G
dr_G/dz = (d G_solv/d k_ij) * k_scale / s_G.
```

EOS owns the pure-solvent infinite-dilution reference sequence and only
returns success after its value and derivative limits converge. Regression
does not copy the reference sequence, Born term, or EOS. Solvent--cation,
solvent--anion, and cation--anion coordinates are separate closed capability
IDs. The active component must be the cation or anion, never the solvent.

For the model-level Figiel dielectric ion-suppression fit, every row supplies
the total ion mole fraction and the observed dimensionless ratio to the
salt-free solvent permittivity. The residual and exact Jacobian are

```text
r_i = (epsilon_model,i / epsilon_saltfree,i - y_i) / s_i
dr_i/dz = (d epsilon_model,i / d alpha)
           * (d alpha/dz) / (epsilon_saltfree,i * s_i).
```

EOS supplies both permittivities and the exact first derivative;
Regression does not copy the dielectric equation. Solvent identity remains
source provenance even where normalization makes the scalar model response
independent of the salt-free solvent value.

The canonical transformed-dataset hash binds every row to exactly one
partition. Row IDs are globally unique, partitions are disjoint, and only
`TRAINING` rows may create Ceres residual blocks. `HELD_OUT` and `STRESS` rows
are evaluated after the primary fit without refitting. Any overlap,
relabeling, or attempt to pass an untouched row into Ceres fails closed.
Without a predeclared acceptance cutoff, their errors remain descriptive and
predictive status remains not adjudicated.

Initial typed observation contracts are:

- pure saturation, vapor-pressure-only, and fixed-pressure-density rows;
- fixed-composition VLE rows with observed phase compositions and lifted phase
  volumes;
- fixed-composition liquid-mixture-density rows with caller-supplied ordered
  component IDs, mole fractions, mixture molar mass, pressure/density scales,
  lifted-volume origins/starts/bounds, provenance, and partitions;
- mean-ionic-activity and staged aqueous-`k_ij` mean-ionic-activity rows;
- solvation-Gibbs and ion-solvation-`k_ij` rows; and
- salt-free-normalized relative-permittivity-ratio rows.

Generic activity, fugacity, osmotic, aggregate, and censored observations are
not current executable contracts.

The fixed-composition liquid-density row is generic at the Regression boundary:
it does not name a chemical, infer a missing parameter, select association
sites, or provide a catalog fallback. The caller supplies the ordered model
component IDs and composition, while the installed EOS capability supplies the
matching model identity, active coordinate, topology, and exact native
derivatives. The active artifact currently advertises the neutral binary
`k_ij`/`l_ij` phase callbacks; a wider component set becomes usable only when a
matching installed EOS descriptor and callback advertise it.

For one row with mixture molar mass `M`, lifted molar volume `V`, observed
pressure `P_obs`, observed mass density `rho_obs`, and affine active parameter
`theta(z)`, Regression sends the caller's composition to EOS and uses

```text
r_P   = (P_EOS(T, x, V, theta(z)) - P_obs) / s_P
r_rho = (M / V - rho_obs) / s_rho
```

The native EOS Hessian supplies the pressure derivatives with respect to `V`
and `theta`; Ceres owns `z` and the log-volume nuisance coordinate. No density
root optimizer or finite-difference production derivative is introduced.

The minimum future downstream observation vocabulary is frozen as follows for
a source-declared positive scale `s`:

```text
identity equality:       r = (y - y_observed) / s
positive log equality:   r = log(y / y_observed) / s
linear aggregate:        r = (a^T y - b) / s
lower censor:            r = max(0, (L - y) / s)^2
upper censor:            r = max(0, (y - U) / s)^2
```

The log equality fails closed unless both values are positive. Aggregate
coefficients and `b` carry source-defined units compatible with `s`. The
squared hinge is continuously differentiable: its derivative is zero on the
inactive side and at the boundary, and is the exact analytic derivative of
the active quadratic branch. It is a censor penalty, never an experimental
uncertainty model. These semantics may enter runtime only with a real
source-bound observation and an authorized evaluator transport.

The package dependency remains EOS-only. Organization decision
[`ePC-SAFT/.github#1`](https://github.com/ePC-SAFT/.github/issues/1) additionally
permits a caller to supply one typed process-local evaluator composed from
installed owner artifacts without creating a Regression-to-Equilibrium
dependency. Only observations explicitly advertised by that exact handle are
authorized. Regression still must not import Equilibrium, run a second
equilibrium implementation, infer phase behavior, or finite-difference a
black-box solve. Predicted bubble/dew/flash and other unadvertised operations
remain unsupported.

The loader rejects missing units, ambiguous composition bases, duplicate row
identities, unsupported phases, nonfinite values, missing source identities,
nonfinite or nonpositive residual scales, scales dimensionally incompatible
with their observed values, source/hash mismatches, and
parameter/observation combinations that the installed artifact does not
advertise. In canonical value units, every scale must satisfy
`0 < s_q < infinity`.

## Exact derivative contract

EOS `1e571ab0a84603a51ed6994b14286f683fb12b88` supplies the first two
general model-bound capability descriptors for neutral binary `k_ij` and
`l_ij`. EOS `86983ff` adds three component-identity descriptors for a
neutral, nonassociating, constant-diameter pure model, backed by the exact
`(n,V,active_parameter)` value/gradient/Hessian callback. EOS
`2d1816cf376294156684fee85611a93fc41d0970` adds complete typed descriptors
and topology fingerprints for the existing active-Born and aqueous-
solvation-factor callbacks. Regression selects the requested known descriptor
from that finite EOS capability set, reports unknown
descriptors as unsupported, and rejects an unknown capability request. Later
families must supply the same closed metadata:

- schema and capability identifiers;
- parameter family and identity shape;
- accepted observation or phase-block contract;
- full state and active coordinate order and units;
- row-major tensor layout and exact buffer dimensions;
- identity-to-coordinate mapping;
- returned value units;
- exact derivative order;
- fixed topology and domain limits;
- local/global scope and explicit non-global claims;
- phase/branch stability and discovery status where applicable;
- certificate and tolerance-set identity;
- explicit unsupported-model, domain, topology, branch, derivative, and
  evaluation failure boundaries;
- EOS source, parameter, topology, transform, and artifact fingerprints;
- accepted-artifact status and receipt identity, or an explicit
  authority-neutral status.

Missing or unknown descriptor metadata prevents `DERIVATIVE_READY` and
`FIT_READY`; Regression does not fill it from local defaults.

The direct-observable descriptors have derivative order one because EOS
returns the observable and its exact total derivative with respect to the sole
active physical parameter. Regression's typed direct-observation contracts
make those exact domains `FIT_READY`; this does not admit fitted values to the
EOS catalog or generalize their fixed reference paths.

Direct observable residuals consume EOS values and first total
derivatives. For a solvation-Gibbs target, solver coordinate `z`, physical
parameter `p = p_origin + p_scale z`, and residual scale `s_G`:

```text
r_G = (G_solv(p) - G_solv,observed) / s_G
J_G = (dG_solv/dp) * p_scale / s_G
```

The mean-ionic-activity contract preserves the source-frozen relative
residual rather than substituting a log residual:

```text
gamma_model = exp(log_gamma_model)
r_gamma = (1 - gamma_model/gamma_observed) / s_relative
J_gamma = -(gamma_model/gamma_observed)
          * d(log_gamma_model)/dp * p_scale / s_relative
```

The implemented direct-observable EOS paths each accept exactly one
active parameter per problem. Parameter sharing is therefore explicit across
rows, while the five Figiel Born targets are five independent one-parameter
problems. Organic ion-solvent pair fits are likewise independent; other pair
values remain explicit row inputs. The separate installed joint pure callback
is the first multi-active seam and supplies all three ordered pure-parameter
columns in one phase evaluation.

The reference organic-ion campaign uses four constructed nearest-pure
endpoints retained by Validation: K+/methanol, Br-/methanol, Na+/ethanol, and
Cl-/ethanol. Three digitized compositions are only near unity. Therefore this
campaign proves exact-derivative fit execution and provides a descriptive
comparison with rounded Table 5 parameters; it is not exact reconstruction of
the paper's pure-organic fitting data, and the remaining mixed-solvent rows
are not duplicated as training residuals.

Lifted phase-equilibrium residuals consume EOS phase-potential values,
gradients, and Hessians. The Hessian supplies exact derivatives of pressure
and chemical-potential residuals with respect to lifted state variables and
active parameters; the capability's coordinate order and row-major layout are
authoritative. No third derivative is required for a Ceres Jacobian whose
residuals are values or first gradients of that phase potential.

For a scalar pure-saturation row with lifted
`V_L = V_L,origin exp(u_L)` and `V_V = V_V,origin exp(u_V)`, the four
residuals are

```text
r_P,L   = (P_L - P_observed) / s_P
r_P,V   = (P_V - P_observed) / s_P
r_mu    = (Phi_n,L - Phi_n,V) / s_mu
r_rho,L = (M/V_L - rho_L,observed) / s_rho
```

The exact parameter column uses
`dP/dtheta = -RT Phi_(V,theta)` and
`d(Phi_n)/dtheta = Phi_(n,theta)`. The lifted-volume columns use the
corresponding `Phi_(V,V)` and `Phi_(n,V)` entries times `V`; the density
column is `-(M/V_L)/s_rho`. Thus the EOS Hessian is sufficient and no
third derivative or density-root sensitivity is introduced. Regression rejects
out-of-bound or inverted volumes and nonpositive mechanical-stability
curvature before accepting an evaluation.

The fixed-topology association surface re-closes the pressure roots for
vapor-pressure, pure fixed-pressure-density, combined saturation
pressure/liquid-density, and fixed-composition liquid-mixture-density problems
at each independent parameter start. A
safeguarded Newton correction begins from the row's retained volume; a bounded
log-volume scan and bisection is used only when that correction cannot remain
inside the declared phase interval. Only mechanically stable roots are
admitted. This prevents a confirmation start from inheriting
primary-parameter phase volumes while preserving the same residual,
derivative, and fail-closed contracts. The immutable EOS build and public
state replay verify this mechanics contract; they do not accept a fitted
parameter tuple.

EOS may use diagnosed density closure for a direct density/property
contract. Branch identity, closure residual, and conditioning then remain
EOS diagnostics returned with the observation. Regression does not
reimplement the closure.

Unsupported family/observation pairs fail before Ceres starts.

## General fitted and lifted coordinate contract

For one problem, let:

```text
z in R^N          ordered fitted solver coordinates
q in R^Q          ordered Regression-owned lifted coordinates
p = p0 + D z      physical fitted coordinates; D is diagonal and nonsingular
theta = S p       ordered evaluator parameter slots; S has shape K_theta x N
y in R^K          ordered evaluator primitive outputs
r = W h(y, q)     ordered scaled residuals in R^R
```

`S` is an immutable sparse sharing map. It can route one fitted coordinate to
multiple declared evaluator slots, but it cannot merge duplicate parameter
identities implicitly. EOS or an authorized downstream evaluator returns
exact primitive derivative blocks `Y_theta = dy/dtheta` and `Y_q = dy/dq`.
With `H_y = dh/dy` and the explicit lifted-coordinate derivative
`H_q = partial h/partial q`, Regression assembles

```text
J_parameter = W H_y Y_theta S D             # R x N
J_lifted    = W (H_y Y_q + H_q)             # R x Q
J_full      = [J_parameter | J_lifted]       # R x (N + Q)
```

The `H_q` term is mandatory when a residual depends directly on a lifted
coordinate, including liquid density and common reporting pressure. All
buffers are row-major. Parameter columns follow `RegressionProblem.parameters`;
lifted columns follow observation order and then the observation contract's
declared local order. Residuals follow observation order and then the
contract's declared residual order. A missing, reordered, nonfinite, or
unsupported derivative column fails closed.

The structural preflight requires `R >= N + Q`; this is necessary, not
sufficient. Let `J_q` be the lifted block and let `U_q` contain the left
singular vectors spanning its numerical column space. The nuisance-projected
parameter matrix is

```text
J_projected = (I - U_q U_q^T) J_parameter   # R x N
```

For `Q = 0`, `J_projected = J_parameter`. Full numerical validity requires
rank `N + Q` for `J_full` and rank `N` for `J_projected`. Both use

```text
tau = 100 * epsilon_binary64 * max(R, N + Q) * sigma_max
rank = count(sigma_i > tau)
```

and report condition numbers in scaled Ceres coordinates. The nonsingular
physical transform `D` preserves exact rank but can materially change
conditioning. Raw physical parameter columns, scaled solver columns, and
nuisance-projected columns are distinct diagnostics. A coordinate on an
active bound is reported separately and is not described as an identified
interior parameter.

Training rows contribute residual blocks and lifted coordinates. Held-out and
stress rows never enter Ceres and are evaluated without refitting. Pure-
saturation reporting keeps fixed fitted parameters and solves exactly three
reporting coordinates `(log V_L, log V_V, log P_common)` against
`(P_L-P_common, P_V-P_common, mu_L/RT-mu_V/RT)`; liquid density remains a
descriptive prediction.

## One engine and result

The intended public surface is:

```text
fit_parameters(problem: RegressionProblem) -> RegressionResult
parameter_capabilities(model) -> tuple[ParameterCapability, ...]
```

`RegressionProblem` owns one ordered parameter block, typed observations,
source identities, numerical controls, and declared primary and confirmation
starts. It contains no backend selector.

`RegressionResult` owns:

- resolved problem, EOS, source, and artifact fingerprints;
- ordered fitted parameters with start, final value, unit, movement, bounds,
  scale, and active-bound distance;
- per-row observed/model values, raw/scaled residuals, derivative status, and
  evaluated/skipped/failed accounting;
- Ceres termination and cost diagnostics;
- full and parameter-projected singular values, rank, and conditioning;
- primary/confirmation agreement;
- separate solver, numerical, physical-validity, workflow, scientific, and
  predictive statuses.

The current installed positive-observation adapter is deliberately narrower
than this general result vocabulary: it accepts only solved-state,
KKT-certified training rows (`N >= 1`, one ordered slot per fitted
parameter), rejects held-out and stress rows at construction, and either
evaluates every required row or fails closed. Its row diagnostics preserve the
source locator and expose the exact scaled Ceres-coordinate Jacobian; they do
not claim reserved-partition evaluation or partial-row accounting. Future
reserved evaluation, lifted coordinates, sharing, aggregate, and censored
observations are separate admitted slices.

The canonical general result exposes only the ordered `parameters`
collection; it has no scalar `parameter` compatibility alias. Accepted
paper-specific entry points may retain their established presentation wrappers
while their internals migrate to this engine, but they must not retain separate
solvers or scientific logic. Each completed migration deletes its superseded
paper-specific native cost path. New paper-named fit functions are forbidden
when the paper can be represented as data plus a `RegressionProblem`.

## Identifiability and numerical preflight

Every observation capability declares its residual and nuisance-variable
dimensions. Before Ceres, Regression verifies that every active parameter has
at least one exact finite derivative path to an observation. A missing or
structurally zero column is a fail-closed input error.

Rank uses the scaled matrices and tolerance defined by the general-coordinate
contract above. For direct observations, the parameter Jacobian is tested
directly. Lifted problems report both the full matrix and the
nuisance-projected parameter matrix. Active-bound coordinates are reported
separately and are not hidden as identified interior parameters.

Primary and confirmation starts receive the same deterministic derivative
preflight. Initial numerical rank deficiency is reported but need not prevent a
nonlinear solve unless all declared starts contain a missing/zero parameter
path. A completed fit cannot be numerically or scientifically valid unless the
final projected parameter rank equals the number of independent active
coordinates. The numerical contract must declare a maximum accepted condition
number and its rationale; the engine supplies no universal scientific cutoff.

This bounded preflight establishes structural and local numerical readiness,
not practical or global identifiability. The retained fixed-2B literature
campaign additionally requires issue 28's nuisance-reoptimized profiles and
declared accepted region before it may report practical-identifiability
evidence; other topology/parameter selections require their own Validation
campaigns. A dense set of correlation-generated rows does not become
independent experimental information merely because the local Jacobian has
full rank.

## Routine Regression tests and Validation campaigns

Regression retains compact deterministic numerical-tolerance tests for the
main generic public workflows and supported active-coordinate variations. The
routine matrix covers accepted methane/ethane joint pure fits, representative
scalar and pair adapters, joint `(m,sigma,epsilon/k)`, fixed-topology
association with 2B and non-2B descriptors, subsets and sharing,
fixed-composition `k_ij`/`l_ij`,
the composed positive-observation transport, and deterministic result records.
Related cases are parameterized when they share one invariant.

These tests use manufactured analytic cases or the smallest source-backed
anchors, exact production derivatives, independent derivative oracles where
needed, and declared
`abs(value-reference) <= atol + rtol*abs(reference)` comparisons. Tolerances
are justified by scale, oracle accuracy, source precision, and deterministic
numerical behavior; expected values are not generated by the implementation
path under test.

The sibling Validation repository owns full literature rows and objectives,
broad temperature/property/composition matrices, profile surfaces,
bootstrap/data-perturbation or parameter-distribution studies, competing
topology and induced/cross-association strategies, transfer evidence, durable
receipts, and independent review. Validation consumes immutable installed EOS
and Regression artifacts through public contracts and does not replace the
compact Regression mechanics sentinels.

The routine package command is:

```bash
.venv/bin/pytest -q
```

`pyproject.toml` excludes the explicit `campaign` marker from that command.
Campaign-marked package replays are opt-in retained evidence, not routine CI or
a merge gate. New full-grid or paper-complete reproductions belong in
Validation as installed-artifact campaigns; they must not be added to
Regression's default suite.

## Acceptance and failure rules

A fit is not scientifically successful merely because Ceres terminates.
Results separately report:

- solver convergence;
- numerical convergence and confirmation-start agreement;
- complete exact Jacobian columns;
- local parameter rank and conditioning;
- active bounds and parameter movement;
- observation-domain, phase, topology, and source validity;
- predictive status based only on a predeclared untouched partition and
  approved criteria.

Rank deficiency, unsupported derivatives, EOS fingerprint mismatch,
incomplete row evaluation, invalid physical state, or confirmation failure
returns a diagnostic result or a fail-closed input error according to whether
Ceres began. No runtime automatically persists fitted parameters into a
EOS catalog.

## Performance contract

Broad input support must not move the Ceres evaluation loop into Python.
`RegressionProblem` is validated and serialized once; C++ owns residual
evaluation, Jacobian assembly, Ceres, and diagnostics. EOS calls are
native, model-bound, and batched by compatible topology/observation contract.

The implementation shall:

- tape or prepare a fixed EOS topology once per active problem and reuse
  it when the EOS capability permits;
- reuse row buffers and avoid per-iteration Python objects;
- preserve the exact sparse parameter-to-observation dependency structure;
- choose a deterministic dense or sparse Ceres policy internally from the
  problem shape, with no public backend selector;
- propagate EOS evaluation budgets, deadlines, cancellation, and failed
  row identities;
- retain a benchmark comparing a migrated reference campaign with its legacy
  native path before deleting that path.

Residual-only trial evaluations request values only. Jacobian evaluations
request the complete exact derivative blocks. An evaluator that cannot honor
the requested mode or buffer contract fails closed; Regression does not
silently calculate a missing derivative or pretend that a derivative-bearing
call was value-only.

### Downstream evaluator transport

Organization decision
[`ePC-SAFT/.github#1`](https://github.com/ePC-SAFT/.github/issues/1) authorizes
this exact inverted process-local direction. Each concrete installed evaluator
ABI still requires package evidence before Regression may consume it. The
transport is model-bound and immutable for the lifetime of a solve. It declares
schema version; evaluator capability ID;
ordered parameter-slot, primitive-output, and lifted-coordinate identities
with units; fixed `S` dimensions and entries; row and batch shapes; row-major
buffer sizes; supported `VALUES_ONLY` and `VALUES_AND_JACOBIAN` requests;
threading/reentrancy rules; applicability, topology/branch, and sensitivity
conditioning statuses; evaluated/skipped/failed row accounting; and EOS,
Equilibrium, model, parameter, topology, and installed-artifact fingerprints.

The caller owns input and output buffers for the duration of a synchronous
call. An evaluator may batch only rows with one compatible immutable
topology/observation contract, and it must preserve canonical row order in its
output. Unsupported requests, identity or unit mismatch, short buffers,
nonfinite values, incomplete columns, topology/branch changes, failed
sensitivities, and cancellation return explicit failure with the affected row
identities. No partial batch is accepted as a complete residual evaluation.
Production transfer remains blocked only until the exact installed artifact
contract and its replay evidence exist.

Performance evidence cannot relax derivative, rank, or physical-validity
gates.

The native generalized core is one contiguous `[N fitted | Q lifted]`
parameter block. Starts supplied by the public general problem are complete
ordered fitted-parameter vectors of length `N`; lifted starts remain owned by
the source-bound observation rows and their declared confirmation policy.
Confirmation agreement gates the fitted vector and relative cost, not equality
of nuisance coordinates that may legitimately converge to different
representations of the same closure. Its build-testing-only closed-form
fixture has `N = 2`,
`Q = 1`, `R = 4`, full rank 3, and nuisance-projected parameter rank 2; it
recovers the same solution from two complete starts, observes value-only and
Jacobian requests separately, exercises the slot-sharing map `(0,1,0)`,
rejects `R < N + Q`, rejects an incompletely written exact Jacobian, and
reports a zero parameter column as rank deficient. It is neither installed
nor exposed as a EOS capability.

Accepted methane and ethane training are the first installed multi-active
cutover. Each four-row campaign uses `N = 3`, `Q = 8`, `R = 16`, the EOS
coordinate order `(n,V,m,sigma,epsilon/k)`, a `16 x 11` exact Jacobian, full
rank 11, and nuisance-projected parameter rank 3. The shared core preserves
the accepted parameter, residual, cost, confirmation, reporting, and status
contracts. Their established public result wrapper remains a presentation
projection; the duplicate pure training Ceres loop has been deleted. The
three-coordinate common-pressure reporting solve remains intentionally
separate from training because it owns reporting closure rather than fitted
parameters.

For this cutover, the primary fitted-parameter vector is the immutable
component specification start. The deterministic confirmation vector adds
`0.1` solver unit to each fitted coordinate (equivalently `0.1` times its
declared physical transform scale, clamped to the existing bounds). Both runs
use the observation-declared lifted-volume starts. This is a start
perturbation only; it changes no source row, weight, bound, scale, or
acceptance threshold.

The public `fit_parameters` transport admits two exact multi-active shapes.
The first is the three ordered pure parameters `(m,sigma,epsilon/k)` on
pure-saturation observations. This closed EOS-backed adapter maps the identity
sharing order `(m,sigma,epsilon/k)` to the EOS
`(n,V,m,sigma,epsilon/k)` callback and delegates training to the generalized
core.

The second is the EOS descriptor-driven fixed-topology association shape. It
accepts any finite ordered active subset advertised by that descriptor,
including ordinary-component and association slots, multiple site pairs, and
explicit many-to-one sharing. The native bridge keeps `N` fitted coordinates
and `S` EOS slots distinct and accumulates the full Hessian columns according
to the declared sharing map. Vapor-pressure, fixed-pressure-density, and
combined saturation rows use its pure model domain. Fixed-composition liquid-
mixture-density rows use the same slot and sharing contracts through the
generic associating-mixture model domain, preserving caller component order and
composition. Exact Jacobian replay, pure 2B/non-2B references, and an anonymous
associating-mixture sentinel establish mechanics only; no association tuple,
topology, practical identifiability, or prediction is accepted. Reactive
association remains outside this neutral descriptor and fail-closed.

The installed callbacks return value, gradient, and Hessian in one inseparable
ABI call. The core propagates Ceres's residual-only request and
does not copy the assembled Regression Jacobian for that call, but the adapter
still receives EOS derivatives. This is an artifact limitation, not a
claim that trial evaluations are sensitivity-free.

The installed pair-family campaigns evaluate all 17 audited May 2015
methane/ethane rows as training data. Each `68 x 35` solve converges with full
rank 35 and projected parameter rank 1. The independent fitted values are
non-bound `k_ij = -0.00843032298906253` and
`l_ij = -0.002774426668544412`, each confirmed from two perturbed starts.
This is in-sample reproduction evidence only; the retained pressure-closure
result remains negative under its historical gate and predictive status remains
`NOT_ADJUDICATED_NO_APPROVED_HELD_OUT_CUTOFF`.

### Constant-`k_ij` maturity hardening

Constant `k_ij` is the priority pair family. Maturity work stays inside the
existing one-active-pair contract and does not add a second fitter, a generic
target registry, simultaneous multi-pair recovery, or a temperature-dependent
`k_ij(T)` coordinate.

One compact campaign matrix shall exercise each currently admitted observation
domain:

1. neutral fixed-composition VLE through the lifted-volume Hessian path;
2. aqueous mean ionic activity through the exact first-derivative path; and
3. single-ion solvation Gibbs energy through the exact first-derivative path.

For every applicable case the campaign must verify the exact EOS
derivative against an independent analytic or automatic-differentiation
oracle where one exists, projected parameter rank one, the declared
conditioning gate, active-bound reporting, agreement of the primary and
perturbed starts, complete row accounting, and deterministic results. Neutral
VLE also verifies invariance to row order and to the caller's
ordering of the unordered pair identity. Direct aqueous cases verify that the
two inactive pair values supplied by each row remain fixed context rather than
additional fitted coordinates. Negative controls must fail closed for an
unsupported EOS fingerprint or observation domain, incomplete EOS
row evaluation, and a terminal projected-rank-zero or conditioning failure.

The 17-row methane/ethane solve is the performance sentinel. On the current
development host it completes in `0.55 s` wall time under the explicit campaign
command; this measurement is diagnostic, not a portable acceptance cutoff.
The campaign shall retain Ceres residual- and Jacobian-evaluation accounting
so a future slowdown can be diagnosed without adding a user-selectable backend.
Ceres remains single-threaded `DENSE_QR` for this `68 x 35` problem because the
measured solve is small and fast. A bounded native solver-time budget is
required so a malformed or unexpectedly expensive case returns a diagnostic
result instead of running indefinitely; the budget does not interrupt an
individual EOS callback and cannot replace EOS-side evaluation
deadlines.

### May methane/propane source-backed reference

May et al. (2015), Table 6, supplies 22 direct methane/propane VLE rows in
`evidence/may-2015-methane-propane-vle.csv`. The source PDF is retained by
SHA-256 `53fd1bdd55dc6807ec76cf88626438d8dfceb3ec09149d4405ea36cfbe6b842a`;
the exact CSV artifact is
`97a07b274dc4da6a281614f3fd39c520ebd6678776413746b13bc8665113c529`.
The source values are kept as `T/K`, `p/kPa`, methane liquid fraction, and
propane vapor fraction with both standard and combined uncertainty columns.
Regression converts pressure by `P_Pa = 1000 p_kPa` and derives the methane
vapor fraction as `1 - y_propane`. Uncertainties are descriptive provenance;
they are not residual weights or acceptance cutoffs.

The test-only bundle is `purpose="user-provided"`: it copies the installed
Gross--Sadowski methane and propane pure records, adds one explicit zero
methane/propane active-pair initialization, and selects the canonical order
`("methane", "propane")`. EOS remains the EOS and exact Hessian owner;
no EOS catalog is changed. Every row uses the observed pressure scale,
unit `mu/RT` scales, liquid start/origin `4.0e-5 m3/mol`, ideal-gas vapor
start/origin `RT/P`, and the declared liquid `[3.0e-5, 2.0e-4]` and vapor
`[5.0e-5, 1.0e-2]` bounds.

The source-backed in-sample campaign fits exactly one unordered constant pair
coordinate. Its lifted residual/Jacobian shape is `88 x 45`, with full rank
`45/45`, projected parameter rank `1`, full condition number
`3894.3041379063716`, no active bound, and two confirmation runs after the
primary start. The three declared starts converge to the same fitted
`k_ij = 0.0038919335722629794`; final cost is `0.03734758119771876`, all 22
rows are evaluated, and every row reports `EXACT_PROVIDER_HESSIAN`. A complete
second campaign rerun was bitwise identical for fitted value and cost (the
observed repeatability deltas are zero; the test retains a `1e-14` diagnostic
tolerance). Solver, numerical, workflow, and physical-validity statuses are
reported separately. May supplies the observations, not this fitted
parameter: the result is in-sample Regression evidence only, with no
prediction, EOS-catalog persistence, production, or scientific authority
claim.

An additional literature reference may be admitted only when its raw
observations, component order, units, `k_ij` convention, fitted comparison
value, and installed EOS capability are all independently available.
Published parameters without their fitting rows remain comparison provenance,
not a regression benchmark. A further neutral VLE case may add component and
state transfer evidence; an already source-bound aqueous case remains valid
cross-domain evidence but is not a substitute for missing neutral rows.

The scalar pure implementation replays the four accepted methane training
rows independently for each family. Every `16 x 9` solve has full rank 9 and
projected parameter rank 1, with a non-bound result and confirmation-start
agreement. The local anchors are `m = 1.0001569260577763`,
`sigma = 3.7063548743836034 angstrom`, and
`epsilon/k = 150.00325287725062 K`. They are deterministic in-sample
implementation evidence, not a replacement for the accepted joint fit,
prediction evidence, or EOS-catalog authority.

## Implementation sequence

1. **Complete.** Introduce the typed capability, parameter, observation,
   problem, and result contracts without adding a second native module or
   target.
2. **Complete for fixed-composition VLE.** Generalize neutral,
   non-associating binary `k_ij` from paper-named input ownership to
   caller-supplied component pairs and rows, matching the advertised
   phase-block domain and retaining reference campaigns as evidence.
3. **Complete for fixed-composition VLE.** Add EOS active `l_ij` support
   and admit it through the same phase observation contract and Ceres owner.
4. **Complete for independent scalar pure saturation.** Admit `m`, `sigma`,
   and `epsilon/k` component coordinates through the same contract, result,
   Ceres owner, and native target. The accepted joint methane/ethane public
   result and reporting presentation remain unchanged while their training
   now uses the same ordered three-parameter Ceres owner.
5. **Complete for the advertised direct-observable domains.** Admit one Born
   diameter from one or more solvation-Gibbs targets, or one solvation factor
   from one or more mean-ionic-activity targets. Each problem remains one
   parameter at a time, consumes the EOS value/first derivative, and
   uses no lifted phase variables.
6. **Complete for the advertised model-level dielectric domain.** Admit one
   ion-suppression coefficient from normalized relative-permittivity rows.
   The 36-row Figiel water/methanol reference fit is a `36 x 1`, rank-1,
   non-bound training solve and recovers `7.067350349980952` versus the
   paper's descriptive `7.01`; the digitized rows define no scientific or
   predictive cutoff.
7. **Complete for model-level ionic-region relative permittivity.** Admit one
   positive SSM+DS coordinate from source-bound solvation-Gibbs targets on an
   installed model-bound callback. Validation subject `81e9a3f` fits all five
   admitted Table S5 ion models independently; every `1 x 1` problem is rank
   1, non-bound, and confirmed from starts 4 and 12, recovering the fixed
   input `8` within `2.2e-9`. This is in-sample mechanics evidence, not
   recovery of a parameter that Figiel fitted or a shared multi-model value.
8. **Complete for component solvent relative permittivity.** Admit one
   positive solvent coordinate from source-bound single-ion solvation-Gibbs
   targets. EOS subject `7cadaad`, Regression subject `177890a`, and
   corrected Validation subject `e4cb7af` bind the exact installed-artifact
   path. Five
   independent `1 x 1` fits are rank 1, non-bound, and confirmed from starts
   50 and 110. Recovery of fixed `78.09` is mechanics evidence, not a
   paper-fitted or predictive result.
9. **Mechanics complete for fixed-topology association; campaign evidence pending.**
   One ordinary-sigma descriptor/request/result path supports arbitrary
   advertised subsets, multiple site pairs, explicit sharing, and any neutral
   component count advertised by the installed EOS. The retained pure 2B and
   anonymous associating-mixture cases are references through that path, not
   special capabilities.
   Practical-identifiability profiles, competing topology or induced-
   association strategies, and source-complete installed-artifact campaigns
   remain Validation work. Polar families remain excluded.
10. **Selected coupled MEA path blocked upstream.** The exact positive-
    observation transport is implemented. The frozen Hilliard/Böttinger
    homogeneous state remains `FEASIBLE_ONLY`. Its local-minimum gate fails,
    and an independent reduced-Hessian calculation confirms negative
    curvature. The reduced pressure-plus-speciation Ceres tracer waits for a
    certified stable installed state and a frozen rank-sufficient `N <= 2`
    coordinate block. The exact sequence is in
    `docs/science/mea-coupled-regression-master-plan.md`.

Every standalone scalar family admitted in steps 1--8 must be independently
fit-ready and reference-validated. Step 9 records executable fixed-topology
mechanics without parameter or topology acceptance. Step 10 is a separate coupled capability: it
does not promote optional standalone families and cannot use their evidence as
a substitute for exact Equilibrium sensitivities.

## Explicit exclusions

- arbitrary Python or C++ residual plugins;
- mutable runtime registration;
- copied EOS or Equilibrium equations;
- numerical production derivatives;
- hidden bounds, scales, starts, units, chemistry, or source defaults;
- automatic experimental-uncertainty interpretation;
- parameter-catalog persistence;
- polar parameter regression;
- claims of global uniqueness, uncertainty, predictive validity, or
  downstream readiness without their separate evidence.

## Public usability and deterministic records

The executable caller path is documented in the
[general Regression quickstart](../general-regression-quickstart.md).
`ObservationDataset.from_records` delegates row validation and transformed-row
hashing to the existing observation and hash owners. `prepare_fit` resolves
installed capability identity and constructs the existing `RegressionProblem`.
Preflight evaluates exact Jacobians at each declared start without invoking
Ceres. `PreparedFit.evaluate` exposes one immutable, read-only assembled
residual/Jacobian evaluation at a declared fitted and lifted solver point.
An omitted lifted point uses EOS root-resolved starts for multi-parameter
pure-phase contracts and row-declared starts otherwise; the returned record
contains the exact selected coordinates. It
preserves ordered fitted, lifted, row, and residual-component identities,
declares row-major layout, reports full and nuisance-projected diagnostics,
and binds parameter, topology, installed-artifact, preparation, source,
dataset, and observation-order fingerprints. It calls the same exact EOS
derivative and Regression assembly owners as preflight, never enters Ceres,
and grants no new authority. `PreparedFit.fit` calls the existing
`fit_parameters` owner.

`RegressionResult.to_record` and `to_json_bytes` expose
`epcsaft-regression-result` schema version 1 for both direct-EOS and
composed-positive producers. Export is deterministic and non-mutating.
Literature context keeps fitted full-precision numbers separate from
source-printed strings and represents absent profile, bootstrap, uncertainty,
and Validation-campaign artifacts explicitly. Literature export requires the
exact prepared-fit provenance. Failed diagnostics identify unavailable
numerical fields explicitly; evaluated data reject nonfinite values. Solver,
numerical, workflow/physical, scientific, predictive, and authority statuses
remain separate. A result record is not a Validation receipt or authority
change.

Fixed-topology capabilities may advertise a mathematically open upper
descriptor bound. Canonical result export omits that nonnumeric field, lists
its path under `unavailable_fields`, and identifies it under `open_bounds`
with the upper direction. Lower descriptor bounds remain finite. Export does
not encode infinity as JSON, silently drop the field, or permit nonfinite
calculated residuals, Jacobians, fitted values, or diagnostics.
