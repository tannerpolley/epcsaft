# Generic Relative-Permittivity Formulation Design

Status: implemented; parameter-catalog portions superseded by ADR 0003

Approved: 2026-07-28

The generic formulation naming and routing described here are implemented.
Its installed-catalog migration and provenance-scanner steps are historical;
EOS now accepts explicit caller parameters under ADR 0003.

Doctrine: `../ePC-SAFT-organization/GOVERNANCE.md`, revision 4

## 1. Objective

Replace every author-named relative-permittivity runtime selector, parameter
family, native symbol, serialized key, and equation identifier with a generic
physical name.

Every parameter bundle will carry one immutable
`relative_permittivity_formulation` model record. Parameter resolution and
native evaluation will route exclusively on that record and on generic
parameter-family names. Author names, publication years, citations, source
identifiers, and locators remain provenance; they do not select equations or
alter runtime control flow.

This is a hard schema and API cut. It preserves the admitted equations,
parameter values, native integer formulation codes, numerical behavior,
derivative authority, and Provider capability.

## 2. Design constraints

The implementation must preserve these boundaries:

- native C++ remains the sole production owner of nonlinear
  relative-permittivity equations;
- CppAD remains the sole nonlinear production derivative authority;
- `ParameterBundle -> ParameterSet -> EPCSAFT.evaluate -> EosResult` remains
  the state/property lifecycle;
- bundle records remain immutable;
- no mutable scientific dictionary, backend selector, hidden physical
  default, compatibility alias, or duplicate equation owner is introduced;
- runtime selection must not inspect a bundle name, source identifier,
  locator, citation, author, or publication year; and
- source-complete caller-provided bundles remain possible through the same
  typed record boundary.

## 3. Formulation record

`ModelParameterFamily` gains the family
`relative_permittivity_formulation`. Its value is exactly one of:

| Value | Physical rule |
|---|---|
| `none` | No electrolyte relative-permittivity model is active. |
| `ion-fraction-suppression` | Suppress solvent permittivity by total ionic mole fraction. |
| `component-permittivity-mixing` | Mole-fraction mixing of component relative permittivities. |
| `density-volume-temperature` | Solvent-density rule using mixed dielectric volume and temperature parameters. |
| `ion-specific-suppression` | Suppress solvent permittivity by ion-specific coefficients. |
| `dipolar-density-perturbation` | Dipolar-density perturbation formulation. |
| `solvent-only` | Use the solvent relative permittivity without ionic suppression. |

Every installed bundle must contain exactly one such record in `model.csv`.
Neutral and associating bundles use the explicit `none` value. Complete
caller-provided bundles have the same requirement.

Absence, duplication, or an unknown value is an invalid bundle. Resolution
must fail before native-model construction. There is no inferred fallback
from the presence of parameter families.

The formulation record uses the existing immutable `ModelParameterRecord`
shape, including `record_id`, `source_id`, `locator`, and `domain_id`.
Those fields document provenance only.

## 4. Generic parameter families

The schema replaces author-bearing family names with physical names:

| Previous family | Generic family |
|---|---|
| `dielectric_ion_suppression_coefficient` | `ion_fraction_suppression_coefficient` |
| `schreckenberg_dielectric_volume` | `dielectric_volume` |
| `schreckenberg_dielectric_temperature` | `dielectric_temperature` |
| `zuber_ion_suppression_coefficient` | `ion_specific_suppression_coefficient` |
| `rueben_dipole_scaling` | `dipole_scaling` |
| `rueben_polarizability_scaling` | `polarizability_scaling` |
| `rueben_correlation_integral_parameter` | `correlation_integral_parameter` |

Already-generic families retain their names, including
`relative_permittivity`, `solvation_factor`, and
`ionic_region_relative_permittivity`.

For each formulation, resolution validates the exact required topology and
parameter-family coverage. A family required by the selected formulation
must be complete for its prescribed component role. A formulation-specific
family is rejected when a different formulation is selected. This makes the
formulation record the selector and the family records its typed inputs,
rather than allowing family presence to become an implicit selector.

Scientific units, finite-value checks, component-role rules, mixing rules,
and domain checks remain unchanged.

## 5. Generic electrolyte topology

Author-bearing electrolyte topology labels are also runtime selectors and
therefore become generic:

| Previous value | Generic value |
|---|---|
| `ionic-ascani-original-born` | `ionic-original-born` |
| `ionic-held-2008-debye-huckel` | `ionic-fully-dissociated-debye-huckel` |
| `held-2008-debye-huckel` | `fully-dissociated-debye-huckel` |

`ionic-ssm-ds` remains unchanged. These renames do not change the ionic
diameter, Debye--Huckel, Born, dissociation, or reference-state equations.

## 6. Resolution and runtime flow

Resolution follows one explicit path:

1. Load the single `relative_permittivity_formulation` model record.
2. Validate the selected electrolyte topology and component roles.
3. Validate the exact generic parameter-family coverage for that
   formulation.
4. Resolve typed values into `_ResolvedModel`.
5. Map the generic formulation to the existing native integer code.
6. Evaluate the unchanged scalar-generic native equation.

The native integer codes remain stable:

| Code | Formulation |
|---:|---|
| 0 | `none` |
| 1 | `ion-fraction-suppression` |
| 2 | `component-permittivity-mixing` |
| 3 | `density-volume-temperature` |
| 4 | `ion-specific-suppression` |
| 5 | `dipolar-density-perturbation` |
| 6 | `solvent-only` |

Native constants, fields, function-local names, diagnostics, comments, and
Python/native serialization keys use only generic terminology. Integer
stability is retained for the internal native boundary; the old strings are
not accepted.

The public Python property becomes
`EPCSAFT.relative_permittivity_formulation`. The existing
`EPCSAFT.permittivity_model` property is removed without an alias or
deprecation shim.

The separately meaningful quantities `relative_permittivity`,
`bulk_relative_permittivity`, and
`ionic_region_relative_permittivity` retain their existing names and
semantics.

## 7. Bundle and admission migration

Every installed bundle and every test fixture is migrated atomically:

- add exactly one explicit formulation record;
- rename formulation-specific family records;
- rename generic electrolyte topology values where applicable; and
- retain the original source identifiers, citations, publication years, and
  locators in bundle provenance.

Bundle directory names may remain source identifiers because
`src/epcsaft/data/bundles/` is the provenance boundary.

The narrow staged aqueous admissions must stop loading or comparing a
particular catalog bundle, source identifier, or domain identifier. They
instead admit a resolved parameter set by generic structural facts:

- exact component identities and roles required by the callback;
- admitted electrolyte topology;
- selected generic relative-permittivity formulation;
- complete generic parameter-family coverage;
- exact active-parameter ordering; and
- the existing finite-value, unit, and domain constraints.

This preserves the callback's bounded scientific surface without making
provenance metadata executable.

## 8. Equation identifiers and generated documentation

Equation identifiers beside runtime owners become generic:

| Previous EqID | Generic EqID |
|---|---|
| `figiel_bulk_permittivity` | `ion_fraction_suppressed_bulk_permittivity` |
| `ascani_bulk_permittivity` | `component_mixed_bulk_permittivity` |
| `schreckenberg_mixed_parameters` | `density_correlation_mixed_parameters` |
| `schreckenberg_bulk_permittivity` | `density_correlated_bulk_permittivity` |
| `zuber_bulk_permittivity` | `ion_specific_suppressed_bulk_permittivity` |
| `rueben_effective_dipolar_density` | `effective_dipolar_density` |
| `rueben_correlation_integral` | `dipolar_correlation_integral` |
| `rueben_bulk_permittivity` | `dipolar_perturbation_bulk_permittivity` |
| `figiel_ionic_segment_diameter` | `corrected_born_ionic_segment_diameter` |
| `ascani_ionic_segment_diameter` | `original_born_ionic_segment_diameter` |

`docs/latex/equations.tex` remains the equation-registry source. Source
citations stay in documentation metadata and bibliography, while equation
identity and runtime ownership use generic physical language.
`scripts/docs/sync_equation_registry.py` regenerates
`docs/equations.md` and `docs/equations_registry.yaml` after the rename.

## 9. Runtime source-name policy

Add `scripts/check_runtime_provenance_names.py`. It scans `src/epcsaft/`
while excluding `src/epcsaft/data/bundles/` and rejects configured
provenance-name tokens in Python, C++, headers, binding declarations,
comments, diagnostics, and string literals.

The token policy data lives at
`src/epcsaft/data/bundles/provenance_names.toml`, inside the excluded bundle
provenance boundary, so the checker itself does not embed author names.
Matching is token-aware and case-insensitive. The initial policy covers every
author-bearing runtime name removed by this migration.

The check is deterministic, has a nonzero exit status on a violation, and
prints the file, line, and prohibited token. It is a static repository gate,
not a scientific unit test. `tests/test_runtime_provenance_names.py` covers
its path exclusion, token matching, and diagnostic behavior.

## 10. Error behavior

The migration keeps structured fail-closed errors. It adds or updates errors
for:

- missing or duplicate `relative_permittivity_formulation`;
- unknown formulation values;
- missing, incomplete, or extraneous formulation-specific families;
- incompatible formulation and electrolyte topology;
- old family, topology, public-property, or serialized-key names; and
- source-name policy violations.

No old value is silently translated. Git history is the recovery path.

## 11. Test-first implementation strategy

Implementation proceeds in observable slices:

1. Add failing record and bundle tests for the required formulation record,
   generic values, and rejection of legacy values.
2. Add failing resolver tests for explicit selection, exact family coverage,
   neutral `none`, and provenance-independent routing.
3. Add failing public/native-boundary tests for the new property, generic
   export keys, topology labels, and stable integer codes.
4. Add failing independent numerical-oracle tests for all seven formulation
   values, then rename the native implementation without changing equations.
5. Add failing staged-admission tests proving source identifiers do not
   control acceptance.
6. Add failing equation-registry and runtime source-policy checks.
7. Migrate installed bundles, fixtures, documentation, and generated
   registries until all gates pass.

At each slice, the focused test must be observed failing for the intended
reason before production code changes.

## 12. Verification and acceptance

The implementation is complete only when:

1. every installed and fixture bundle has exactly one immutable explicit
   formulation record;
2. runtime resolution uses only the generic formulation and generic
   parameter-family names;
3. the seven native integer codes and all scientific equations are unchanged;
4. independent relative-permittivity numerical oracles pass at their existing
   tolerances;
5. native derivative and fixed-pressure callback tests pass;
6. staged admissions are invariant to provenance-name changes;
7. the legacy public property and every legacy runtime selector, family,
   native symbol, export key, diagnostic, and EqID are absent;
8. the source-name policy check passes outside the bundle provenance
   boundary;
9. equation Markdown and YAML are freshly generated and registry tests pass;
10. the complete Provider test suite, lint, type checks, and
    `git diff --check` pass; and
11. the repository cleanup audit reports no task-owned residue.

## 13. Explicit exclusions

This design does not authorize:

- changing a relative-permittivity, ionic-diameter, Debye--Huckel, Born, or
  reference-state equation;
- changing parameter values or numerical tolerances;
- adding runtime-configurable equations or a generic parameter overlay;
- adding compatibility aliases or accepting legacy strings;
- moving citations out of provenance and scientific documentation;
- renaming bundle directories solely to hide their provenance;
- broadening the staged aqueous callbacks or catalog persistence;
- changing CppAD derivative authority; or
- changing Provider capability, receipt, promotion, or release status.
