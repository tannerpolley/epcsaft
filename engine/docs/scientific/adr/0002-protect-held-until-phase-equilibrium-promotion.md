---
status: fulfilled
date: 2026-08-10
amends: 0001-one-general-phase-equilibrium-owner.md
fulfilled_by: governance D-030 and CAP-09 promotion receipt
migration_source_path: engine/docs/equilibrium/adr/0002-protect-held-until-phase-equilibrium-promotion.md
migration_source_commit: 9f8728b3bffbd548aebb0356a32e31798b36c385
migration_source_blob: 077021af554ca17daaaa3c9170b63a8e5a8982b5
---

# Protect HELD until the phase-equilibrium promotion gate passes

Kind: numerical-method

## Context

ADR 0001 set the correct destination—one source-neutral nonreactive
phase-equilibrium owner—but prescribed an atomic replacement. CAP-09 then
removed the mature `tp_flash` HELD/HELD2 workflow and introduced a separate
Gibbs-minimization implementation. That was not the intended lineage. The
requested change was to evolve and eventually rename the validated HELD/HELD2
workflow while adding explicit phase selection.

HELD2 already carries expensive numerical and literature evidence, including
the accepted eight-feed Khudaida Figure 2 campaign. Deleting that incumbent
before a candidate proves the same workflows makes the migration itself the
largest uncontrolled scientific change.

## Decision

The exact pre-CAP-09 `tp_flash` HELD/HELD2 route is restored as the protected
incumbent. `phase_equilibrium` is a HELD-derived candidate interface, not a
separate physical solver. During the transition:

- both public operations remain explicit and neither automatically falls back
  to the other;
- both use the same HELD/HELD2 numerical core and canonical EOS owner;
- `tp_flash` retains its result, diagnostics, active-parameter sensitivity,
  source-domain admission, and HELD2 certificate behavior unchanged;
- `phase_equilibrium` may add typed feed amounts, phase roles, and constrained
  semantics without changing an incumbent discovery calculation;
- a declared role mismatch fails closed and never becomes an initialization
  hint or a different topology claim; and
- historical HELD/HELD2 names remain live only for the protected incumbent and
  its exact evidence, not as selection keys for new numerical algorithms.

The candidate implementation now carries exact protected incumbent discovery,
declared one-to-three-phase all-component local solves with at most one vapor,
fixed-T coexistence with all-component or neutral-component partner support,
non-authoritative start ordering, and exact fixed-topology sensitivities. This implementation status does
not itself promote the candidate. Until the installed-wheel campaign, two
independent reviews, remaining-tree audit, and explicit user approval are bound
to one accepted receipt, `phase_equilibrium` is not the sole owner and must not
claim the completed ADR 0001 boundary.

## Promotion gate

Retiring `tp_flash` requires a reviewed, user-approved receipt over one exact
installed wheel that proves:

1. every admitted HELD/HELD2 incumbent workflow is accepted by the candidate
   with state, objective, residual, certificate, and sensitivity differences
   inside the incumbent's declared tolerances;
2. the complete accepted Khudaida Figure 2 eight-feed sweep is no worse than
   `tp_flash` against the retained descriptive paper comparison;
3. constrained mode actually solves and certifies each supported declared
   liquid/vapor topology rather than post-filtering discovery;
4. one-, two-, and claimed three-phase boundaries, neutral and charged cases,
   hints, permutation behavior, and failure semantics have proportional
   numerical falsifiers; and
5. independent phase-science and solver-architecture reviewers accept the
   implementation and remaining-tree audit.

A failed comparison is retained evidence and blocks promotion. It never
selects `tp_flash` automatically. ADR 0001 remains the destination after this
gate; its atomic-retirement sequencing is superseded by this decision.

## Consequences

- The unrelated CAP-09 scratch phase solver is removed rather than retained as
  a parallel implementation.
- Existing HELD2 certificates and validation evidence remain executable.
- The temporary extra public interface is deliberate migration evidence, not
  a compatibility promise.
- The one-owner target remains unchanged, but current capability claims are
  limited to what the installed-wheel evidence proves.

## Fulfillment

On 2026-08-12 the complete installed-wheel Khudaida promotion campaign passed
under the user-approved `1e-4` lineage-only tolerance while every stricter
physical and derivative gate remained unchanged. Governance D-030 therefore
made `phase_equilibrium` the sole public route, mapped omitted fixed-T,P
topology to the retained internal discovery controller, prohibited fallback
for declared topology, and retired the transitional Python `tp_flash` surface.
The HELD/HELD2 numerical core remains live internally as required by this ADR.
