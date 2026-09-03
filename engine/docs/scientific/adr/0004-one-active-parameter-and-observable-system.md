---
status: superseded
superseded_by: greenfield rebuild, engine/docs/greenfield/TRANSFER.md
date: 2026-08-07
migration_source_path: engine/docs/adr/0004-one-active-parameter-and-observable-system.md
migration_source_commit: 9f8728b3bffbd548aebb0356a32e31798b36c385
migration_source_blob: dc2546945b7c6dbf15ea291ad5260f1a1afce766
---

> **Superseded by the greenfield rebuild, 2026-09-02.** The decision recorded here was taken
> and is not being rewritten; what it decided no longer describes the tree. Owner compilation leaves the fit path entirely; Stage 10 is `fit(model, observations)` on Ceres with no observable-owner system.
> See `engine/docs/greenfield/TRANSFER.md` for the replacement and its stages. The plan this
> record links to was deleted with the architecture it specified; `git log` carries it.


# Use one active-parameter and owner-observable system

Kind: computational-architecture

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

- One fitted coordinate has one identity from input through `LeastSquaresResult`.
- Mixed component, pair, association, Born/solvation, auxiliary-model, and typed
  correlation coefficients can share one active set when EOS supports them.
- Regression cannot implement EOS or Equilibrium observation families.
- D-034 and D-035, retained by ADR 0007, admit the Equilibrium-owned fixed-T,P homogeneous
  reactive and fixed-topology reactive bubble/VLE families through the same
  owner-neutral Ceres route. D-031 continues to defer unqualified
  phase-equilibrium observations.
- Generated templates provide explicit required neutral identity values;
  loaders do not infer omissions.
- Engine tests falsify numerical claims; literature acceptance remains in
  analysis sets.

The execution authority is
[`../../plans/2026-08-07-general-scientific-capability-hard-cutover.md`](../../plans/2026-08-07-general-scientific-capability-hard-cutover.md).
