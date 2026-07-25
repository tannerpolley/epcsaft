# Provider Equation Documentation Mirror Design

Status: approved design; implementation not started

Approved: 2026-07-24

Doctrine: `../ePC-SAFT-organization/GOVERNANCE.md`, revision 4

Legacy evidence:

- repository:
  `~/Workspaces/Engineering/ePC-SAFT-legacy/repositories/ePC-SAFT-legacy-monorepo`;
- pinned review commit:
  `7920e56dbb0fc1fb418bae44e67ef50f334b75cd`;
- source path: `docs/latex/equations.tex`;
- Git blob: `2e970476258558de05391d9d77f31d1736d78de8`; and
- local SHA-256:
  `a0f2f3aff12042a8a50c5254c2af23e04b60af3e886b97b2756e6ecad800bacc`.

Migration classifies that file as retained lab evidence. It does not carry
current runtime or scientific authority into Provider.

## 1. Objective

Add a complete, readable equation mirror for the mathematics implemented by
the current Provider C++ owners. Preserve the useful organization and
source-tracing metadata of the legacy `equations.tex`, but reconstruct the
content from the current Provider implementation and scientific contract.

The mirror must make it inexpensive to answer:

- which defining relation Provider implements;
- which source, derivation, or convention supports it;
- which C++ location owns it;
- which thermodynamic basis and units it uses; and
- whether the documentation and implementation have drifted.

This work changes no equation, parameter, numerical algorithm, public API,
native ABI, runtime dependency, capability, or authority.

## 2. Authority and ownership

The three surfaces have distinct roles:

1. `src/epcsaft/_native/**` remains the sole production implementation of the
   EOS and native defining equations; the small exact reference/property
   transformations currently owned by `src/epcsaft/reference.py` remain
   Python-owned.
2. `docs/science/full-provider-eos.md` remains the normative contract for
   admitted capability, sources, units, conventions, derivative policy, and
   failure boundaries.
3. `docs/latex/equations.tex` becomes the exhaustive human-readable mirror of
   the mathematics implemented by the current Provider.

The LaTeX document is not a runtime input and cannot override C++ or broaden
the admitted capability. CppAD remains the sole nonlinear production
derivative authority.

## 3. Chosen architecture

Use stable equation identifiers on both sides of the documentation boundary.

Each defining relation in LaTeX has metadata of the form:

```tex
% EqID: residual_helmholtz
% Source: \cite{Figiel2025}, Eq.~(1)
% Status: Source equation
% Description: Sums the active residual Helmholtz contributions.
% Change note: Reconciled to the current scalar-generic Provider owner.
```

The corresponding C++ owner carries a nearby comment:

```cpp
// EqID: residual_helmholtz
```

A Python-owned transformation uses the same identifier with its native comment
syntax:

```python
# EqID: practical_osmotic_coefficient
```

A Provider-owned synchronization script parses the LaTeX and both implementation
surfaces and generates:

- `docs/equations_registry.yaml`, the machine-readable traceability registry;
- `docs/equations.md`, the searchable GitHub-rendered equation index.

`docs/latex/references.bib` contains only sources cited by the current
Provider mirror. A PDF may be built locally with a documented `latexmk`
command, but no generated PDF is tracked.

Rejected alternatives:

- A registry-first design would make YAML the scientific editing surface and
  weaken direct review of notation and equation structure.
- Generating LaTeX from implementation comments would place presentation
  markup and bibliography concerns inside production source.
- Human-maintained file and line references without stable identifiers would
  permit silent drift and break during ordinary refactoring.

## 4. Equation coverage

The document retains the legacy catalogue's broad progression while including
only current Provider mathematics:

1. symbols, ordering, thermodynamic bases, constants, and unit conversions;
2. segment parameters, temperature-dependent diameters, and formulation-bound
   ionic diameter rules;
3. pair-diameter, dispersion, and association mixing relations, including
   distinct `k_ij`, `l_ij`, and hydrogen-bond semantics;
4. packing moments and packing fraction;
5. hard-sphere, hard-chain, and dispersion intermediates and residual
   Helmholtz contributions;
6. association strengths, mass-action residuals, and site-fraction state
   definition;
7. admitted bulk-relative-permittivity rules;
8. finite-size Debye--Huckel contribution;
9. original Born and corrected SSM+DS Born contributions;
10. total residual Helmholtz energy;
11. exact thermodynamic identities for compressibility, pressure, chemical
    potential, fugacity coefficient, and residual caloric properties;
12. density residual and admitted root-domain definition;
13. electrolyte reference sequence, activity, mean ionic activity,
    solvation/transfer energy, and practical osmotic coefficient; and
14. native phase-coordinate and exact algebraic transformations that are
    Provider thermodynamics rather than transport mechanics.

The inventory includes defining equations and exact thermodynamic identities.
It excludes ordinary validation checks, buffer handling, serialization,
iteration histories, and other control flow.

## 5. Derivative boundary

The new document does not retain the legacy catalogue's expanded analytical
nonlinear gradients, Hessians, dielectric derivatives, association
derivatives, or higher-order density/composition/temperature derivatives.

Its derivative section records only:

- independent variables and fixed discrete topology;
- dependent thermodynamic quantities;
- the exact property identity or implicit-function relation consumed by an
  admitted result; and
- CppAD as the owner of the nonlinear partial tensors.

Exact linear transformations and exact implicit-function linear algebra may
be documented when the implementation consumes them. They are not alternate
derivative kernels. The documentation must not suggest that an expanded
handwritten derivative formula is a production owner.

## 6. Synchronization contract

The checker fails closed when:

- a LaTeX `EqID` is duplicated;
- an implementation `EqID` has no LaTeX entry;
- an implemented LaTeX equation has no C++ or Python owner comment;
- an equation lacks required source, status, description, or change metadata;
- an identifier uses an unknown status or malformed name;
- a generated Markdown or YAML file differs from current source; or
- a prohibited handwritten nonlinear derivative entry is introduced.

An equation may have multiple implementation-owner references when one
defining relation is reused. Generated records include current file, line, and
nearby symbol context. Authored documentation does not depend on line numbers.

Documentation-only thermodynamic identities must use an explicit
documentation or exact-derivation status. They cannot be counted as direct
C++ coverage.

## 7. Reconciliation method

Reconstruct one equation family at a time:

1. identify the current scalar-generic C++ owner or exact Python transformation
   owner;
2. identify the current scientific-contract entry and source locator;
3. compare the legacy equation only as provenance and structural guidance;
4. write the current equation using Provider units, bases, and conventions;
5. attach the stable `EqID` beside the defining implementation expression;
6. generate and inspect the registry and Markdown view; and
7. verify that the family contains no obsolete derivative or model variant.

The legacy file is not copied wholesale. Historical equations that are not
implemented by current Provider remain in the legacy repository.

## 8. Proposed file delta

Add:

- `docs/latex/equations.tex`;
- `docs/latex/references.bib`;
- `docs/latex/README.md`;
- `docs/equations_registry.yaml`;
- `docs/equations.md`;
- `scripts/docs/sync_equation_registry.py`; and
- one focused equation-registry test module.

Modify only:

- current native equation-owner files and `reference.py` to add `EqID`
  comments;
- `docs/science/full-provider-eos.md` to point to the local checked mirror and
  state its non-runtime role;
- the smallest existing contributor or test entry point needed to document
  the synchronization command.

No native or Python function body, build target, packaging input, or public
symbol is changed by the documentation slice.

## 9. Verification

The implementation gate requires:

1. bidirectional `EqID` validation;
2. generated Markdown and YAML freshness;
3. unique and complete metadata;
4. declared equation-family coverage;
5. negative space for handwritten nonlinear derivative entries;
6. an optional successful local LaTeX build when `latexmk` is available;
7. the existing complete Provider test suite, proving runtime results are
   unchanged;
8. lint and `git diff --check`; and
9. the repository cleanup audit.

The LaTeX toolchain is documentation tooling only. It does not become a
runtime or wheel dependency.

## 10. Explicit exclusions

This design does not authorize:

- importing all 172 legacy equations;
- retaining obsolete analytical derivative catalogues;
- copying legacy runtime code or generated registries;
- changing current equations to match legacy notation;
- adding a second EOS, derivative, density, or association owner;
- documenting Equilibrium-, Regression-, Ceres-, or Ipopt-owned mathematics
  as Provider equations;
- committing generated PDFs or broad paper archives;
- changing capability, receipt, promotion, or release status; or
- adding a dependency on the legacy or lab repositories.
