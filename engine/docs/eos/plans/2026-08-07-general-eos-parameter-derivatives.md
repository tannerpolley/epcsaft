# General EOS Parameter-Derivative Plan

- Status: completed owner contract retained after CAP-01 through CAP-05
- Date: 2026-08-07
- Owner: EOS and `Parameters`
- Parent: [General scientific capability hard cutover](../../plans/2026-08-07-general-scientific-capability-hard-cutover.md)
- Decision: [ADR 0004](../../scientific/adr/0004-one-active-parameter-and-observable-system.md)

This document is the retained scientific owner contract and CAP provenance. It
is not a live execution queue; current work is selected from GitHub Issues.

## Outcome

Create one EOS derivative path that accepts any ordered set of supported,
continuous, fixed-topology parameter identities for any admitted component
topology. Remove the fixed binary, fixed ternary aqueous, named-solvent, and
single-parameter derivative routes. EOS owns values, partial derivatives,
reference conventions, and pure-saturation sensitivities; consumers never
reimplement these equations.

## Scientific model decisions

### Pair coordinates

`k_ij`, `l_ij`, and the admitted `k_HB,ij` correction are unordered component
pairs. The first two correct cross dispersion energy and cross segment size.
`k_HB,ij` corrects the typed cross-association-energy combining rule; it is not
an association-volume parameter and is not a directed coefficient. Association
site roles select eligible edges, while the pair correction remains reciprocal.
Its identity therefore contains an unordered component pair and an unordered
association-role pair. One coordinate applies to every native edge matching
that declared role pair, including multiplicity-equivalent edges. A request
without a role pair is accepted only when topology resolves exactly one
eligible role pair; otherwise it fails as ambiguous. Independent role pairs
require independent coordinates.
See [the primary-source review](../../research/2026-08-07-pair-interaction-parameter-symmetry.md).

Direct constant `k_ij` values are the first implementation slice. A typed
temperature correlation later exposes each continuous coefficient through the
same parameter identity and tape, not through a new evaluator.

### Canonical Born formulation

The runtime formulation is named `born`. It includes explicit, physically
named correction coefficients for the solvent-shared and dielectric-saturation
terms currently described in the literature as SSM+DS. Generated parameter
templates emit every coefficient. The fixed model topology admits exactly
`(c_shell, c_dielectric) = (0, 0)`, `(1, 0)`, and `(1, 1)` through one
equation; no separate `original-born`, SSM, or SSM+DS runtime branch exists.

The source-neutral equation-registry identities are
`born_outer_shell_diameter`, `born_radial_charging_work`, and
`born_residual_helmholtz`. For ion `i`, the sole equation is parameterized as

```text
D_i = d_i^Born [1 + c_shell (f_mix - 1) / |z_i|]
W_Born = sum_i x_i z_i^2 [
  (1 - 1/epsilon_bulk) / D_i
  + c_dielectric (1 - 1/epsilon_i) (1/d_i^Born - 1/D_i)
]
```

with the existing physical prefactor converting `W_Born` to residual
Helmholtz energy. `d_i^Born` is a required positive length (angstrom at the
document boundary); `f_mix`, `epsilon_bulk`, and `epsilon_i` are required
positive dimensionless quantities. `c_shell` and `c_dielectric` are required
dimensionless model-definition fields, not continuous fit coordinates.
`(0, 0)` is the canonical uncorrected Born limit, `(1, 0)` is the shell-only
SSM model, and `(1, 1)` admits shell and dielectric-saturation corrections.
Parameters resolution rejects `(0, 1)` because `c_shell = 0` makes
`D_i = d_i^Born` and nulls the dielectric increment at every state; admitting
it would create a globally redundant, observationally indistinguishable model
identity. It is not silently canonicalized or aliased. The old
`original-born`, `ssm-ds`, source-named equation IDs, and omitted-field forms
are rejected rather than migrated or defaulted.

Fully dissociated Debye--Huckel remains a distinct physical contribution. It is
not folded into Born and is not selected by a literature source.

### Reference-dependent quantities

Solvent composition, selected ion pair, ion stoichiometry, phase, and reference
convention are explicit typed inputs. Component names and component positions
carry no behavior. A reference-dependent quantity reports its convention with
its value and derivative block.

## Current defects to retire

### Priority 0

- Dormant fixed-aqueous route:
  `active_aqueous_kij_observation`,
  `thermodynamic_derivatives_active_aqueous_kij`,
  `evaluate_ion_solvation_kij`, `evaluate_aqueous_miac_kij`,
  `evaluate_aqueous_miac_kij_batch`, `NativeIonSolvationKij`, and
  `NativeAqueousMiacKij`.
- `active_aqueous_k_ij` optional plumbing in native Helmholtz, derivative, and
  polar structures.
- `evaluate_mixture_phase_active_parameter`, which assumes exactly two
  components and one pair coordinate.
- Python capability/derivative dispatch that admits only a binary single-pair
  coordinate.

### Priority 1

- MIAC solvation-factor differentiation tied to `{0,+1,-1}`, solvent index 0,
  ion indices 1/2, and one active scalar.
- Ion-fraction permittivity differentiation tied to exactly three components
  and one coefficient.
- Literal `"water"` selection in reference construction.
- Separate `original-born` formulation and source-named correction symbols.
- Fragmented optional active slots in the native derivative input.

### Priority 2

- Exactly-one-neutral-solvent restrictions where the physical equations admit
  an explicit solvent composition.
- Author-named polar constants or runtime families; citations remain in the
  equation registry and research notes.

## Target interface

### Resolved active coordinate

Each coordinate contains:

- stable physical family and source-neutral identity;
- scope: component, unordered pair, association edge/rule, model coefficient,
  or typed correlation coefficient;
- canonical component/site identities and resolved native indices;
- physical value, unit, open/closed domain, and optional affine optimizer scale;
- transformation/application rule, including reciprocal pair writes;
- topology and model applicability; and
- correlation family and coefficient identity when applicable.

The descriptor contains no observed target, data path, author name, optimizer
policy, or fit-specific alias.

### Native tape

Record one scalar-generic EOS defining function with coordinate order:

```text
[state coordinates..., active parameter coordinates...]
```

The state prefix is chosen by the owning EOS observable. The active suffix is
the exact order of `ActiveParameterSet`. Before evaluation, all active values
are applied to a local resolved model through one generic dispatcher. Pair
updates write both matrix orientations. Association corrections resolve onto
all eligible role-aware edges and propagate the exact combining-rule chain.

The tape core is established first for coordinates whose physical application
rules are already resolved. Born/solvation and reference-dependent coordinates
are admitted only by the canonical Born/reference cutover; they then use this
same dispatcher and tape rather than a specialized callback.

The result includes value, gradient/Jacobian, required higher partial tensors,
coordinate identities, units, model fingerprint, and structured status. It
also declares each tensor's state and active-coordinate axes, derivative order,
shape, symmetry convention, and layout. An owner must request every mixed or
higher-order CppAD partial consumed by its residual, certificate, implicit
linearization, or exact-Hessian route; unsupported orders fail before solving.
It does not expose a public C++ ABI or installed header.

### EOS observable blocks

The initial owners are:

1. fixed-state EOS properties with explicit state and basis;
2. reference-dependent electrolyte quantities with complete reference input;
3. bulk auxiliary-model quantities such as relative permittivity; and
4. pure-saturation state values and exact total selected-parameter
   sensitivities.

Each descriptor is versioned and declares its applicability, physical domain,
rows, state/reference metadata, and required derivative order. The returned
batch block preserves those declarations and distinguishes row failure from a
shared-solve failure. Regression cannot select a native target index or
contribution-specific hidden value.

A pure-saturation block also records the liquid and vapor branch/bracket
identities, root residuals, `dR/d(rho)` conditioning, and rejection status for
root swaps, coalescence, or branch loss. A finite converged solve without this
branch evidence is not an evaluable saturation row.

## Implementation sequence

1. **Identity consolidation.** Replace family-specific native slots with the
   resolved descriptor and ordered value vector. Reject duplicates, unsupported
   topology changes, invalid domains, and unresolved identities before taping.
2. **Generic application.** Implement arbitrary-width component, pair,
   association, model, and correlation coefficient application where the
   physical rule is already resolved. Start with direct `k_ij`; keep one typed
   extension seam without dormant callbacks. Born/reference coordinates enter
   only in steps 5-6 after their sole physical semantics are fixed.
3. **Pair correctness.** Canonicalize unordered pairs, write reciprocal matrix
   entries, and add typed `k_HB,ij` resolution to selected cross-association
   edges with exact chain rule.
4. **Tape unification.** Route existing fixed-state and pure-parameter
   derivatives through one active suffix. Delete scalar/binary/aqueous callback
   variants once the numerical comparisons pass.
5. **Born cutover.** Name the sole formulation `born`, express correction
   coefficients physically, require explicit fields, and delete the original
   selector and source-named branches.
6. **Reference cutover.** Replace component-name and positional assumptions
   with explicit solvent composition, ion-pair identity, stoichiometry, phase,
   and reference convention.
7. **Observable blocks.** Expose fixed-state, reference-dependent,
   permittivity, and pure-saturation values/Jacobians through the common block.
8. **Dead-surface deletion.** Remove all symbols listed above, generated docs
   entries for deleted callbacks, narrow tests, and stale guidance in the same
   atomic cutover.

## Numerical verification

Persistent EOS checks must falsify numerical claims, including:

- arbitrary-N mixtures with one and several active unordered pairs;
- reordered components producing permuted but identical values/Jacobians;
- a ternary electrolyte and a multisolvent/multi-ion composition without fixed
  positions;
- a mixed active set spanning component, pair, association, Born/solvation,
  and permittivity coefficients where derivatives are implemented;
- `k_HB,ij` reciprocal storage and role-aware cross-edge application;
- each analytic/CppAD column against an independent central difference at
  interior points using tolerance justified by conditioning;
- one parameterized Born matrix at a nondegenerate state with `f_mix != 1`
  compares `(0, 0)`, `(1, 0)`, and `(1, 1)` against independent original-Born,
  shell-only SSM, and SSM+DS equations and confirms that `(0, 1)` cannot
  resolve to an evaluable thermodynamic model;
- explicit-reference invariance under component reordering; and
- structured failures at domain, topology, density, and association-solve
  boundaries without fabricated values.

Named literature datasets and broad matrices move to installed-wheel analysis
sets. Tests must not assert callback names, exception class or prose,
loader mechanics, schema wording, compatibility, or provenance labels.

## Failure semantics

Resolution errors are deterministic preflight failures. EOS evaluation errors
identify the observable row, state, active coordinate set, physical domain, and
native cause. No missing field is inferred, no invalid point is clipped, no
finite-difference derivative is substituted, and no large residual is returned.

## Issue-ready decomposition

Do not publish these units until the four canonical plans are approved.

| ID | Work unit | Acceptance |
| --- | --- | --- |
| EOS-01 | Resolved active-coordinate descriptor | Mixed families resolve in one stable order with complete domains |
| EOS-02 | Arbitrary-width native application and tape | Supported pre-Born N-component/P-parameter cases return declared tensor orders without specialized callbacks |
| EOS-03 | Pair and cross-association correction handling | Reciprocal `k_ij`/`l_ij` and role-aware `k_HB,ij` checks pass |
| EOS-04 | Canonical Born formulation | The four-input nested-model matrix passes, `(0, 1)` fails before evaluation, and no alternate Born selector remains |
| EOS-05 | Explicit electrolyte reference model | No component name/position controls behavior |
| EOS-06 | Generic reference/permittivity derivatives | Multi-ion/multisolvent Jacobians match independent differences |
| EOS-07 | EOS observable derivative blocks | Versioned fixed-state and pure-saturation blocks preserve active-set order, applicability, tensor order/layout, and typed row/solve status |
| EOS-08 | Specialized-route deletion and audit | Listed symbols, tests, and generated references are absent |

## Done definition

EOS has one active-parameter derivative machinery, one Born equation, explicit
reference inputs, arbitrary supported component and parameter counts, exact
owner-provided blocks, and no source-specific or fixed-aqueous runtime path.
