---
status: accepted
date: 2026-08-07
---

# Use one general nonreactive phase-equilibrium owner

## Context

The former decision proposed one facade that delegated to separate neutral,
electrolyte, and reactive source-derived owners while retaining old public
operations. That preserves the route split and its fixed binary, fingerprint,
and source-case assumptions.

## Decision

`phase_equilibrium` is the sole nonreactive phase-equilibrium owner. It supports
one to three liquid/vapor phases in discovery and constrained-specification
modes. Hints may reorder work but cannot change the scientific claim. Neutral
and charged systems use the same physical owner; neutral topology simply omits
charge/electrolyte operations.

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

- The neutral HELD1/electrolyte HELD2 split and fingerprint-qualified flash
  routes are removed during implementation.
- Callers may supply expected phase count/types when known or request discovery.
- Equilibrium consumes the canonical `ActiveParameterSet` and obtains EOS
  partials for exact implicit solved-state sensitivities.
- Named literature cases remain Validation evidence, not runtime identities.

The execution authority is
[`../plans/2026-08-07-general-phase-equilibrium.md`](../plans/2026-08-07-general-phase-equilibrium.md).
