# Full Provider EOS Design

Status: approved design; local implementation complete through density closure

Approved: 2026-07-16

Doctrine: `../ePC-SAFT-organization/GOVERNANCE.md`, revision 2

## 1. Objective

Extend the clean provider from its promoted explicit-density neutral PC-SAFT
baseline to the smallest complete, source-backed provider implementation that
supports:

- hard-chain, dispersion, association, Debye--Huckel, and Born residual
  Helmholtz contributions;
- the Figiel--Yu--Held solvation-shell plus dielectric-saturation (SSM+DS)
  Born formulation;
- explicit-density and specified-pressure state evaluation;
- residual chemical potentials, fugacity coefficients, and residual caloric
  properties; and
- exact nonlinear derivatives required by those provider calculations.

The optimization objective remains minimum maintained surface subject to those
capabilities producing source-faithful, reproducible results. This is a
source-faithful distillation of the lab's proven equations, parameters, and
numerical evidence. It is neither a bulk port nor an independent rewrite: each
clean owner is traced to an archived owner or independent retained oracle,
while legacy configuration and compatibility surface is deliberately omitted.

## 2. Chosen architecture

Use one deep, scalar-generic residual Helmholtz kernel and expand it through
independently testable capability leaves. Every state result, derivative, and
solver consumes that kernel.

Rejected alternatives:

- A bulk port followed by pruning would import legacy configuration modes,
  duplicate derivative routes, compatibility seams, and tests that protect the
  old implementation rather than the required behavior.
- A generic contribution-plugin framework would add speculative abstractions,
  registration, and configuration without an admitted external contribution
  owner.

The provider retains one Python package, one native module, and one native
target unless an observed build constraint proves that another is necessary.

## 3. Capability boundary

### Included

- Neutral, associating, and ionic mixtures represented by source-bearing
  `ParameterSet` records.
- Association with explicit site topology and site-pair parameters.
- Debye--Huckel electrostatics.
- Born SSM+DS with fitted Born diameters, component solvation factors, bulk
  relative permittivity, and an explicit ionic-region relative permittivity.
- The Figiel 2025 salt-free solvent mass-fraction dielectric rule and ion-
  fraction suppression rule for the bounded systems to which their evidence
  applies.
- CppAD derivatives in temperature, molar density, and composition through the
  orders consumed by the admitted provider outputs.
- Branch-safe density closure for requested vapor and liquid roots.
- Contribution-resolved and total residual Helmholtz energy, compressibility,
  pressure, residual chemical potentials, fugacity coefficients, residual
  entropy, internal energy, enthalpy, isochoric heat capacity, and isobaric
  heat capacity.

### Excluded

- Ideal-gas caloric properties and total caloric properties; the current
  provider does not own ideal-gas heat-capacity data.
- Chemical reactions, dissociation or ion-pair equilibrium, electroneutrality
  repair, and automatic species generation.
- Phase equilibrium, phase discovery, global stability, Ipopt, regression,
  Ceres, and fitted-parameter workflows.
- A public high-order derivative tensor API or native consumer SDK before a
  real extension consumer contract is approved.
- Single-ion activity, infinite-dilution, solvation, and transfer workflows
  before their reference paths and literature campaign are separately
  admitted.
- Alternative derivative backends, model/backend selectors, compatibility
  aliases, mutable scientific dictionaries, and undocumented defaults.

## 4. Public contract

The lifecycle remains:

```text
ParameterBundle -> ordered ParameterSet -> EPCSAFT -> EosResult
```

`ParameterBundle` and `ParameterSet` remain the only parameter lifecycle.
`EPCSAFT` remains the only model owner. `EPCSAFT.evaluate` remains the state
entry point and accepts exactly one state specification:

```text
(temperature, molar_density, mole_fractions)
or
(temperature, pressure, mole_fractions, phase=None)
```

`phase` is rejected for explicit-density evaluation. For specified-pressure
evaluation, `None` accepts exactly one resolved mechanically stable root and
rejects an ambiguous multi-root state. The admitted explicit values are
`"vapor"` and `"liquid"`; they select the lowest- and highest-density stable
roots when multiple roots exist. If only one stable root exists, its diagnostic
branch is `"single"`; a hint never relabels it as liquid or vapor.

`EosResult` remains the single state result. It gains only fields produced by
the admitted full-provider capability. A small immutable density diagnostic
record is returned for pressure-specified states and is absent for explicit-
density states.

No public argument selects a derivative backend or contribution
implementation. Model activation follows the resolved record topology:

- no sites and no charged components: hard-chain plus dispersion;
- sites with complete association pair records: add association;
- nonzero charges with the complete admitted electrolyte records: add
  Debye--Huckel and Born SSM+DS.

Partial topology is an input error. An ordinary charged composition must be
electroneutral within a declared numerical tolerance; the provider does not
invent counterions or alter composition.

## 5. Resolved scientific input

The existing record model remains the sole source of values, units, source
locators, domains, and fingerprints. One bundle-level model-parameter record
is added for a scientific value that belongs to the formulation rather than
to one component or pair. It follows the same provenance, unit, domain,
immutability, and fingerprint rules as the existing records. The consuming
resolver performs model-specific completeness checks and creates one immutable
native model.

The existing families cover segment count, segment diameter, dispersion
energy, charge, relative permittivity, Born diameter, solvation factor,
association energy and volume, and pair interactions. The clean electrolyte
model adds only bundle-level source-bearing model records demonstrated by its
equations, notably:

- the ionic-region relative permittivity used by dielectric saturation; and
- any fixed coefficient in the admitted dielectric suppression relation that
  is not a universal physical constant.

The archived `eps_ion = 8.0` implementation is not copied as a hidden
constant. Figiel--Yu--Held 2025 states that 8 is the assumed relative
permittivity of the ionic region and attributes its origin to Buelow et al.; it
must enter the model through a cited record and participate in fingerprints.
Likewise, the coefficient 7.01 in the Figiel dielectric suppression rule is
source-bearing model data, not a universal constant.

Sparse semantics are preserved: omission means absent, while a sourced zero is
a present value. The resolver never substitutes zero, one, a combining rule,
or another component's value for a missing unlike-component record unless that
exact relation is the admitted source equation. The diagonal `k_ii = 0` is the
exact self-interaction identity and is not an empirical pair record.

## 6. Governing kernel

The scalar-generic kernel evaluates

```text
a_res = a_hc + a_disp + a_assoc + a_DH + a_Born,SSM+DS
```

in reduced molar residual Helmholtz form. It retains the promoted Gross--
Sadowski hard-chain and dispersion implementation and distills only the
required association and electrolyte equations from lab evidence and primary
sources.

The slice-local scientific specification must identify, for every consumed
equation, whether it is a primary-source equation, exact derivation,
implementation convention, or bounded approximation. The broad lab equation
catalogue is evidence, not an editable governing specification for this
provider.

### Association

Association sites and allowed site pairs come from `SiteRecord` and
`AssociationParameterRecord`; there is no string association scheme. Site
fractions solve the Wertheim mass-action residual for the exact fixed
topology. The value solver runs outside CppAD recording. CppAD supplies the
partials of the defining residual and Helmholtz term, and implicit-function
linear solves supply the temperature, density, and composition sensitivities
consumed by provider properties.

The implementation rejects nonconvergence, nonphysical site fractions,
singular or ill-conditioned association Jacobians, and topology changes.

### Debye--Huckel and dielectric model

The Debye--Huckel term uses SI number density, explicit integer charge, and the
admitted bulk relative-permittivity function. Solvent dielectric values are
combined on the salt-free solvent mass-fraction basis and then suppressed by
the admitted total ion mole-fraction relation. Its source domain is retained;
the provider does not silently generalize the relation beyond the reported
solvents, temperature treatment, or ion-fraction range.

### Born SSM+DS

The clean Born owner implements Figiel--Yu--Held 2025 equations 7--10 as
resolved against the article text and supporting evidence:

- each ion has a source-bearing fitted Born diameter;
- each component has the required solvation factor, with ion values of one
  represented as records rather than defaults;
- the shell thickness is charge-normalized;
- the bulk and ionic-region dielectric regions are distinct; and
- SSM and DS are both part of this admitted formulation, not independent
  runtime switches.

Before coding this term, the slice scientific specification must reconcile the
typeset equation-7 extraction with the article's explicit ionic-region
permittivity prose. The retained primary article and any supporting material
are the authority; the archived C++ implementation is only a comparison
candidate.

## 7. Derivatives and properties

CppAD is the sole nonlinear production derivative authority. It records the
defining scalar equations for fixed component and association topology; it
does not record density or association iterations.

The internal derivative payload is limited to what the public provider uses:

- first density derivative for compressibility and pressure;
- first and second density derivatives and temperature cross derivatives for
  density diagnostics and residual heat capacities;
- first composition derivatives for residual chemical potentials and
  fugacity coefficients; and
- the association implicit sensitivities required by those derivatives.

The implementation follows the thermodynamic definitions from the common
Helmholtz owner. Contribution-level property views are diagnostic decompositions
of that same derivative result, not separate equation kernels.

`C_P` is admitted only where the required pressure derivatives are finite and
the density branch is mechanically stable. Near a singular density or
association Jacobian, the provider reports failure instead of a finite-looking
property.

A successful fugacity calculation does not imply a validated single-ion
activity convention or literature-scale conversion.

## 8. Density closure

Density closure solves

```text
R(rho) = rho * R * T * Z(T, rho, x) - P_target = 0
```

over a physically bounded positive-density domain derived from the packing
fraction limit. A deterministic logarithmic scan identifies sign-changing
brackets. Each candidate is refined without recording the solver history.

Every accepted root records its bracket, residual, `dR/d(rho)`, iteration
count, and mechanical-stability status. Failed closure reports its reason
through `EosError`. In a multi-root state, vapor selects the lowest-density
stable root and liquid selects the highest. There is no warm-start or cache
path. Nonpositive or numerically singular slopes are rejected.

Cross-seeded tests must prove that a vapor request cannot return a liquid root
and vice versa. The provider makes no critical-point or global phase-stability
claim from this local density closure.

## 9. Errors and diagnostics

Parameter/model construction errors continue to use `ParameterError` and its
stable context. State and solver failures use one public `EosError` subclass of
`ValueError` with a stable code and immutable context. Failed states raise
`EosError`; successful pressure states carry density-solver diagnostics.
Source and fit domains remain explicit on `ParameterSet` and bound admitted
state inputs; a returned EOS value is not a separate physical-validation
certificate. Message parsing is not an API.

Diagnostics distinguish:

1. solver convergence;
2. numerical convergence under tolerance, scaling, and initialization checks;
3. physical validity under source domain, electroneutrality, topology, and
   branch checks.

No layer implies the next.

## 10. Evidence and validation

Production tests retain only compact evidence reachable from this capability.
The lab retains papers, broad tables, old implementation variants, and failed
experiments. The validation repository tests immutable installed provider
artifacts and never imports the lab or sibling source tree.

Required package evidence includes:

- primary-source equation and parameter locators;
- independent high-precision or numerical checks of each contribution;
- CppAD directional checks in temperature, density, and composition;
- association mass-action residuals and implicit-sensitivity checks;
- neutral and zero-charge reduction to the promoted Slice 1 result;
- low-density, pure-component, composition-boundary, and component-permutation
  behavior;
- dielectric, Debye-length, Born-diameter, shell-thickness, and charge-scaling
  checks;
- explicit rejection of non-electroneutral, incomplete, out-of-domain, and
  singular states;
- branch-adversarial density tests; and
- wheel-build and isolated installed-artifact parity.

Required validation evidence includes at least one bounded associating system
and one bounded electrolyte system with retained literature data and a plotted
model comparison. Tolerances must be justified from source precision,
independent numerical convergence, and the quantity being tested. Agreement
with archived implementation output can characterize a migration but cannot
by itself establish physical validity.

## 11. Continuous implementation leaves

Implementation proceeds without routine approval stops through these local,
independently verifiable leaves:

1. Extend the typed model resolver and write the slice-local scientific
   specification.
2. Add association values and implicit derivatives.
3. Add the bulk dielectric model and Debye--Huckel contribution.
4. Add the reconciled Born SSM+DS contribution with source-bearing ionic-region
   permittivity and dielectric coefficient.
5. Add composition and temperature derivatives plus residual chemical,
   fugacity, and caloric properties.
6. Add branch-safe density closure and specified-pressure evaluation.
7. Reconcile documentation and update the architecture cost ledger.
8. Validate the
   immutable installed artifact, and prepare promotion receipts.

Each leaf starts from tests, removes any displaced owner in the same commit,
and keeps one kernel and one public lifecycle. A leaf may be locally committed
without changing runtime authority. Publication, release, and promotion remain
separate explicit gates.

The program stops only for an unresolved primary-source equation,
missing parameter/license evidence, a numerical result that invalidates the
design, or an external-write/authority decision. It does not stop merely
because one internal leaf completed.

## 12. Completion criteria

The design is complete when:

- every included capability is reachable through the one provider lifecycle;
- every retained production file, dependency, and public symbol serves that
  capability or its direct proof;
- all five Helmholtz contributions and consumed derivatives have unique
  owners;
- neutral, associating, and bounded electrolyte fixtures pass source-backed
  package and installed-artifact validation;
- pressure-specified states preserve requested density branch identity;
- no lab, migration, sibling source, alternate derivative backend, hidden
  scientific constant, compatibility seam, or duplicate input/result owner is
  present; and
- the architecture ratchet records and justifies the final maintained surface.
