# Retire residual Polar while preserving dipolar relative permittivity

Kind: computational-architecture

- Status: accepted
- Date: 2026-08-27

## Context

The residual point-multipole formulation added fixed dipole and quadrupole
inputs, a model-family selector, an inlined native contribution, and a public
`ares_polar` result. The project retired that formulation from Engine scope.
The separately sourced `dipolar-density-perturbation` relative-permittivity
equation also uses dipolar terminology but is not a residual Polar term.

## Decision

Delete the residual Polar family, formulation, physical-moment inputs, native
equations, aggregate contribution, public result, tests, and current packet and
documentation selections. Do not retain an alias, parser, dormant zero, or
compatibility route.

Keep `dipolar-density-perturbation`, `dipole_scaling`,
`polarizability_scaling`, `correlation_integral_parameter`,
`bulk_relative_permittivity`, and their temperature and density derivative
behavior in the electrolyte relative-permittivity owner. Historical results,
source transcriptions, accepted plans, and ADR 0003 remain provenance.

This record supersedes only the final consequence of ADR 0003 that required
open polar leaves to report one aggregate `ares_polar` contribution. Every
other ADR 0003 decision and consequence remains accepted.

## Verification and validity boundary

Public installed-wheel EOS tests falsify the retained contribution closure,
relative-permittivity equation and derivatives, state identities, density
closure, saturation, reference behavior, and active parameter derivatives.
Source, include, configured Meson/link, and symbol scans prove that the deleted
header and residual surface are absent. These checks do not validate broad
fluid or literature domains, which remain analysis evidence.
