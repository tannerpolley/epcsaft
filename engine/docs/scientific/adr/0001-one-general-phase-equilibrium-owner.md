---
status: accepted
date: 2026-08-07
amended_by: 0002-protect-held-until-phase-equilibrium-promotion.md, 0005-native-engine-ownership-and-generated-evidence.md
migration_source_path: engine/docs/equilibrium/adr/0001-one-general-phase-equilibrium-owner.md
migration_source_commit: 9f8728b3bffbd548aebb0356a32e31798b36c385
migration_source_blob: c285ce8150732c2c706bbe1503248e335e8a5930
---

# Use one general nonreactive phase-equilibrium owner

Kind: computational-architecture

The destination is implemented. [ADR 0002](0002-protect-held-until-phase-equilibrium-promotion.md)
protected the HELD/HELD2 transition until Governance D-030 accepted the
promotion receipt and made `phase_equilibrium` the sole public route.

## Context

The former decision proposed one facade that delegated to separate neutral,
electrolyte, and reactive source-derived owners while retaining old public
operations. That preserves the route split and its fixed binary, fingerprint,
and source-case assumptions.

## Decision

`phase_equilibrium` is the sole nonreactive phase-equilibrium owner. Omitted
fixed-T,P topology runs the source-neutral Phase-Stability Discovery sequence
selected by D-041. Declared fixed-T,P topology enters GREPE directly for one to
three ordered liquid/vapor phases and never falls back to discovery. Hints may
reorder declared-topology work but cannot change the scientific claim. Neutral
and charged systems share the same discovery stages and GREPE owner; charge
topology selects only the required coordinate chart and constraints.

Physical formulations and specialized numerical kernels may remain distinct
internally when selected by typed topology or conditioning. Author, paper,
campaign, fingerprint, component name, and fixed component position never
select runtime behavior. Success reports the exact topology and stability or
globality claim certified.

Homogeneous chemical equilibrium remains a separate owner. Reactive phase
equilibrium is deferred; only a shared state, phase, parameter, and sensitivity
seam is designed now. The cutover retains no facade aliases or old neutral and
electrolyte routes.

## Consequences

- HELD, HELD2, Mitsos--Barton, and Perdomo remain citations and historical
  lineage, not runtime keys or public types.
- Callers may supply expected phase count/types when known or request discovery.
- Equilibrium consumes the canonical `ActiveParameterSet` and obtains EOS
  partials for exact implicit solved-state sensitivities.
- Named literature cases remain Validation evidence, not runtime identities.

The execution authority is
[`../../equilibrium/plans/2026-08-07-general-phase-equilibrium.md`](../../equilibrium/plans/2026-08-07-general-phase-equilibrium.md).
