# Provider Repository Context

This clean repository is the production home for ePC-SAFT provider code. Its
optimization objective is minimum maintained surface subject to source-backed,
reproducible thermodynamic results.

`governance_doctrine_revision: 2`

Canonical local doctrine: `../ePC-SAFT-organization/GOVERNANCE.md`.

The local candidate contains one typed parameter lifecycle and one scalar-
generic residual Helmholtz kernel for hard-chain, dispersion, association,
Debye--Huckel, and corrected Born SSM+DS contributions. It returns residual
chemical potentials, fugacity coefficients when defined, residual caloric
properties, pressure, and branch diagnostics for pressure-specified states.
The installed catalogs are Gross--Sadowski methane/ethane and the bounded
Figiel--Yu--Held water/methanol/ethanol alkali-halide source bundle. One reference-state
module owns x-process infinite-dilution solvation and transfer plus single-ion
and stoichiometric mean ionic activity on one explicit electroneutral formula-
unit path.

Pressure closure accepts exactly one of molar density or pressure. An
unlabelled pressure state must have one mechanically stable root. A genuine
multi-root state requires `phase="vapor"` or `phase="liquid"`; a phase hint
never renames a lone root.

The currently accepted authority receipt is still
`promotion-0009-provider-slice-1`, limited to the neutral explicit-density
baseline. The full local candidate has not been pushed, independently
validated, promoted, or released.

Phase equilibrium, phase discovery/global stability, reaction and species
equilibrium, multi-salt reference paths, regression, native consumer transport,
and ideal-gas caloric data remain outside this repository's current capability.

The architecture and cost ratchet are recorded in `ARCHITECTURE.yaml`; the
scientific equation owner is `docs/science/full-provider-eos.md`.
