---
status: accepted
date: 2026-08-07
---

# Use one active-parameter and owner-observable system

## Context

EOS, Equilibrium, and Regression currently contain overlapping parameter
capability maps and observation dispatch. Several routes are constrained by a
paper, fixed component count/order, one pair, or one target. That duplication
makes a scientifically valid parameter derivative depend on which consumer
asks for it and encourages Regression to own thermodynamic calculations.

## Decision

`Parameters` owns one ordered `ActiveParameterSet` of continuous,
fixed-topology physical parameter identities. EOS accepts it unchanged for
partial derivatives. Equilibrium accepts it unchanged and propagates total
solved-state sensitivities. Regression consumes owner-provided observable
values and Jacobians in the same order and owns only residual and optimization
policy.

Owners advertise source-neutral observable descriptors with explicit inputs,
units/basis, domain, row identity, derivative support, and structured failure.
Runtime branches may represent a distinct physical formulation or numerical
algorithm, but not an author, source, campaign, application, fingerprint,
literal component, fixed component position, or special species tuple.

The sole solvation formulation is `born`; explicit neutral correction
coefficients reduce it algebraically to the original Born expression. Pair
parameters are unordered component-pair identities. A typed `k_HB,ij`
correction applies to eligible role-aware cross-association-energy edges and is
not an association-volume or directed parameter.

This is a hard cutover. Superseded parameter slots, callbacks, registries,
schemas, aliases, and fallbacks are deleted with their replacement.

## Consequences

- One fitted coordinate has one identity from input through `FitResult`.
- Mixed component, pair, association, Born/solvation, auxiliary-model, and typed
  correlation coefficients can share one active set when EOS supports them.
- Regression cannot implement EOS or Equilibrium observation families.
- Governance D-034 and D-035 admit the Equilibrium-owned fixed-T,P homogeneous
  reactive and fixed-topology reactive bubble/VLE families through the same
  fit/result route. D-031 continues to defer nonreactive phase-equilibrium
  observations.
- Generated templates provide explicit required neutral identity values;
  loaders do not infer omissions.
- Engine tests falsify numerical claims; literature acceptance remains in
  Validation.

The execution authority is
[`../plans/2026-08-07-general-scientific-capability-hard-cutover.md`](../plans/2026-08-07-general-scientific-capability-hard-cutover.md).
