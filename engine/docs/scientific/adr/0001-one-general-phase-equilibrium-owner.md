---
status: accepted
date: 2026-08-07
amended_by: 0002-protect-phase-stability-until-promotion.md
migration_source_path: engine/docs/equilibrium/adr/0001-one-general-phase-equilibrium-owner.md
migration_source_commit: 9f8728b3bffbd548aebb0356a32e31798b36c385
migration_source_blob: c285ce8150732c2c706bbe1503248e335e8a5930
---

# Use one general nonreactive phase-equilibrium owner

Kind: computational-architecture

The destination is implemented. [ADR 0002](0002-protect-phase-stability-until-promotion.md)
protected the HELD/HELD2 transition until the D-030 decision retained by
[ADR 0007](0007-retire-governance-and-validation-owners.md) accepted the
promotion result and made `phase_equilibrium` the sole public route.

## Context

The former decision proposed one facade that delegated to separate neutral,
electrolyte, and reactive source-derived owners while retaining old public
operations. That preserves the route split and its fixed binary, fingerprint,
and source-case assumptions.

## Decision

`phase_equilibrium` is the sole nonreactive phase-equilibrium owner. Omitted
fixed-T,P topology retains the admitted HELD/HELD2 discovery boundary: binary
neutral mixtures and one- or two-phase charged results. Declared fixed-T,P
topology supports one to three ordered liquid/vapor phases. Hints may reorder
work but cannot change the scientific claim. Neutral and charged systems use
the same physical owner; neutral topology simply omits charge/electrolyte
operations.

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

- HELD/HELD2 remain internal topology-appropriate numerical kernels; their
  fingerprint-qualified Python route and duplicate public schema are removed.
- Callers may supply expected phase count/types when known or request discovery.
- Equilibrium consumes the canonical `ActiveParameterSet` and obtains EOS
  partials for exact implicit solved-state sensitivities.
- Named literature cases remain analysis evidence, not runtime identities.

The execution authority is
[`../../equilibrium/plans/2026-08-07-general-phase-equilibrium.md`](../../equilibrium/plans/2026-08-07-general-phase-equilibrium.md).
