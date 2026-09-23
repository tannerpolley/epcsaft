# Complete model specification and parameter compatibility

[Overview](README.md) · [Physics](physics-and-families.md) · [Derivative contract](derivatives-and-evidence.md)

## What the user selects

A complete specification applies to the modeled mixture, including its
cross-interactions. It is more precise than a list of enabled energy terms or
a paper name. Individual contribution choices remain inspectable and editable
inside that specification.

| Part of the specification | Information to retain |
|---|---|
| Equation choices | Exact reference, dispersion, association, multipolar, electrostatic and Born formulation identifiers and revisions; correction/erratum adoption. |
| Species and chemistry | Species identities, charge numbers, declared supports, full-dissociation or explicit-speciation convention, and chemical standard-state inputs when needed. |
| Association topology | Site types, component mapping, multiplicities, allowed edges, solvation versus pair labels, and construction of unlike strengths. |
| Shared constitutive rules | Packing, contact, electrostatic and Born diameters; their relations; permittivity; solvent averaging; cross parameters and temperature dependence. |
| Parameter roles | Fixed physical inputs, inherited effective parameters, newly fitted coefficients and active fitting coordinates. |
| Reference conventions | Ideal thermal inputs, solvation/transfer conventions, activity scale, infinite-dilution definition and electrolyte concentration basis. |
| Calibration attachment | Exact fitting formulation, fixed/fitted partition, observation identities, weights, constraints, objective and transformations. |
| Evidence attachment | Independently tested species, states, observables, derivatives, failures and uncertainty. |

Parameter and evidence records refer to the specification; they need not be
stored inside a single runtime object. Numerical policy is separately recorded:
solver tolerances, stopping criteria, association method and derivative method
affect numerical reproducibility but do not define a different molecular
association topology. The existing ownership boundaries should be retained.

## Example selections, expressed conceptually

These are examples of information, not a proposed callable API or accepted
input schema.

| Choice | Figiel baseline | A proposed association research variant |
|---|---|---|
| Reference/dispersion | Declared PC-SAFT reference and source mixing/ionic pair rules | Initially preserved explicitly; reassessed during calibration. |
| Association | Source solvent sites, ions without sites | Solvent sites plus separately labeled ion-solvation and ion-pair sites. |
| Electrostatic/solvation | DH and corrected Born SSM+DS | A deliberately declared combination; feedback from pairing must be resolved. |
| Dielectric | Figiel composition/decrement relation and source thermal inputs | Unchanged or deliberately replaced, with the changed dependencies recorded. |
| Parameters | Identified source parameter bundle | New variant identity; inherited values are starting guesses until revalidated. |
| Claims | Only the baseline's retained evidence | New evidence required for every affected claim. |

The second column is a baseline description, not automatic certification of
every state. The third is a new hybrid hypothesis, not a reproduction of Yang's
full ePPC-SAFT model. Reproducing Yang also requires its polar, MSA, solvent,
mixing and applicable NAHS choices.

## Double counting and missing interactions

| Failure mode | Example | Required response |
|---|---|---|
| Direct duplication | Adding full DH and full MSA descriptions of the same ion correlations | Ordinary selection is incompatible; a special partition requires its own derivation and evidence. |
| Inconsistent decomposition | Adding pairing without deciding which ions screen the remaining charges | Specify the combined model, internal equations and complete derivative chain. |
| Hidden fitted compensation | Adding explicit dipolar attraction to a nonpolar fit | Revalidate and, as necessary, refit affected parameters, including shared base parameters. |
| Missing cross-interaction | Enabling dipole–dipole and quadrupole–quadrupole only | Declare whether dipole–quadrupole interaction is included, neglected or unavailable, and justify the approximation. |

Concrete compatibility rules follow:

1. Use one coherent reference-fluid/dispersion/contact package.
2. Use one long-range ionic description for a given interaction unless an
   explicit combined derivation partitions it.
3. Select one Born formulation. Direct, shell-only and shell-plus-saturation
   are alternatives/limits of a solvation model, not three energies to sum.
4. Specify association topology and every active unlike interaction.
5. Do not describe the same pair as both a chemical species and a site bond
   without reconciling balances, standard states and electrostatic populations.
6. Distinguish a residual polar term from dielectric response. Both may use
   dipole inputs while representing different quantities.
7. Born plus ion–solvent association is not inherently inconsistent. Published
   models use both; their physical partition and joint fit must be assessed.

Biazi explicitly acknowledges an approximation when separately adding solvent
hydrogen bonding and ion–solvent association. That is evidence of a declared
modeling assumption, not proof of exact independence.
[Biazi2025](references.md#biazi2025)

## Four levels of compatibility

**Mathematical compatibility** requires consistent units, reduced-energy basis,
reference conventions, required inputs, topology and regularity of requested
derivatives. The Engine can enforce declared mathematical requirements and
report missing inputs or numerical failures.

**Physical compatibility** concerns overlap, omitted interactions and the
meaning of coupled approximations. Known incompatibilities can be explained
by selection rules. Metadata cannot prove that physical double counting is
absent. New combinations require source-based reasoning and evidence.

**Calibration compatibility** asks whether the bundle was fitted for these
equations, sites, mixing rules and reference states. The data layer can compare
the recorded fitting specification with the selected one. A mismatch identifies
a research variant, not a valid transfer by default.

**Empirical applicability** asks where the formulation/parameterization pair
has useful predictive evidence. Fitted ranges and failures remain visible in
data and analyses. In accordance with current policy, publication ranges do
not become arbitrary Engine prohibitions on exploratory calculations.

## Parameters belong to a formulation

Segment count $m$, diameter $\sigma$, and attraction $\epsilon/k_B$ are
usually effective model parameters. A measured molecular dipole moment, a
fitted Born diameter and an ion-specific dielectric scaling coefficient have
different meanings. The last need not be a literal dipole moment of an ion.
An unusually small effective size is not by itself a transferable molecular
dimension. [Esper2023](references.md#esper2023), [YangPair2024](references.md#yangpair2024)

Record the fitted/fixed partition rather than assuming every value was fitted
or every inherited value is immutable. Changing a formulation does not prove
that every coefficient needs refitting. It does require checking the affected
parameters and claims; even shared base parameters may have compensated for
the old physics.

For fair comparisons, retain both a fixed-parameter diagnostic and an
appropriately constrained refit. The former explains the immediate effect of
an equation change; the latter tests an alternative model's achievable
predictive quality. Keep fitting and validation observations distinct, and
compare on the same definitions, data and uncertainty basis.

## Smallest integration with the existing repository

| Existing owner | Future responsibility |
|---|---|
| Native input builder and Python mapping | Resolve the explicit specification and parameters; validate mathematical requirements once. |
| Native Contributions assembly | Consume those choices and reuse common state intermediates. |
| Native Properties, equilibrium and regression | Reuse one differentiated thermodynamic model; no paper-specific property or equilibrium implementations. |
| Data parameter/provenance records | Attach fitting formulation identity, fixed/fitted roles and source transformations. |
| Analyses | Retain comparison partitions, numerical qualification and independent physical validation. |
| Applications or a future selection interface | Offer complete starting selections that expand into explicit physics and compatible parameter bundles. |

Current source anchors are `engine/native/contributions/helmholtz.hpp`,
`engine/native/core/parameters.hpp`, `engine/native/surface/parameter_builder.cpp`
and [the packet contract](../../../../data/docs/packet-contract.md). These are
an inspected starting point, not promises that every proposed field is admitted.

Prefer extending the existing parameter records and admission owner. Do not
prebuild a plugin registry, dynamic loading system, parallel catalogs or
general compatibility theorem engine. A physical-value fingerprint can identify
normalized inputs while citations remain separate provenance, as required by
[ADR 0003](../../adr/0003-explicit-parameter-inputs.md). A hash is identity, not
evidence of scientific accuracy.
