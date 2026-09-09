# Equation and convention corrections — 2026-09-06

Scope: the uncommitted documentation reconciliation on attached branch
`greenfield`, baseline `ff614901be22727095947e71fb4626ef2dae7a7c`.
This review changes documentation and local readiness tooling only. It does
not qualify numerical Engine behavior or execute the next evidence campaign.

## Corrections from the deeper review

The earlier reconciliation left several inherited contradictions. They were
corrected in the canonical definitions and regenerated YAML index:

- Separate packing, DH finite-size and Born diameter records; equality requires
  a source-specific relation. The current native owners already read separate
  packing and DH records.
- Identify mechanical phase energy separately from reference-complete energy.
  Define dimensional chemical/charge multipliers, dimensionless dual multipliers
  and optimizer callback multipliers independently, including the RT scaling.
- Mark the historical trial-pressure source-activity transfer as unresolved
  provenance. It cannot override the fixed-P0 standard-state conversion. Its
  derivation belongs to reference-thermochemistry qualification; changing every pressure symbol would not resolve it.
- Make Shomate enthalpy's kJ/mol-to-J/mol factor explicit and distinguish the
  NIST H coefficient from an enthalpy observable. Verified against the
  [NIST WebBook unit definitions](https://webbook.nist.gov/cgi/cbook.cgi?ID=C7732185&Mask=203).
- Permit zero component amounts in the composition definition at positive total
  amount; distinguish interior logarithms, central limits and insertion actions.
- Identify the DH polynomial as a finite Taylor approximation, and terminal
  reference molality as a limit estimator. Neither is an exact finite-point
  identity. Point to the retained dispersion coefficient arrays/source table.
- Remove remaining live stage-status/evidence-class/formulation claims from
  `engine/ARCHITECTURE.yaml`; it points to the permanent doctrine and current
  transfer plan.

These are documentation/derivation findings. They are not newly measured code
failures. The dated transfer audit remains unchanged, including its execution
limitations and original evidence tables.
