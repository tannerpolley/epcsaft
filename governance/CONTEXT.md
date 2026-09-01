# Governance Context

`governance/` is the doctrine and control-record directory inside the private
`tannerpolley/ePC-SAFT` monorepo. `GOVERNANCE.md` revision 26 is canonical.

The repository root coordinates work. `engine/` owns the sole `epcsaft`
runtime; `data/` owns reusable scientific packets; `validation/` owns
installed-wheel campaigns.

Source repositories and their former names are migration provenance after
archive. The organization-profile repository remains separate and thin.
Artifacts and rollback bundles live in external local state, not in Git.

## Current scientific program

The current program matures General Reactive Equilibrium and Phase Equilibrium
(GREPE) as one generic, configured-physics equilibrium owner and applies it in
two independent downstream programs:

1. aqueous MEA--H2O--CO2 model selection, regression, validation, and
   manuscript revision in `tannerpolley/MEA-Thermodynamics`, followed by a
   separately scoped predictive lane in `tannerpolley/MEA-Absorption-Column`;
2. configured reactive-electrolyte LLE, application parameterization, and an
   equation-oriented IDAES surrogate in `tannerpolley/Lithium_Extraction`.

Neither application selects an Engine route by component, paper, or operation
name. Declared state roles, phases, species support, reactions, observations,
and model capabilities compile the problem.

The canonical generic design is
`engine/docs/equilibrium/plans/2026-08-14-general-reactive-equilibrium-phase-equilibrium.md`.
Current work is tracked only in native GitHub issues. GREPE owns declared
temperature-solved Bubble-T and Dew-T observations through the generic
fixed-pressure coexistence formulation. The simultaneous sparse backend from
closed issue #51 is historical; issue #54 owns the active Bülow 2021
reproduction and integration. Historical plans, completed issues, campaign
reports, and issue mirrors retain evidence but are not parallel work queues.

## Regression architecture

Regression consumes one or more opaque EOS or anchored D-034/D-035 Equilibrium
owner handles. One call owns one physical start, bounds, affine coordinates,
ordered observations, optional diagonal discrepancy scales, and native Ceres
controls. The result retains fresh final physical predictions and Jacobians,
the weighted optimizer Jacobian, typed failures, counters, and local covariance
status. Regression owns no thermodynamic model construction, topology,
multistart, campaign stopping rule, parameter adoption, or Python iteration
callback. Equilibrium retains its inner Ipopt solves and exact Hessians.

## Cross-repository authority

MEA-Thermodynamics, MEA-Absorption-Column, and Lithium_Extraction are official
downstream applications under D-038. They remain separate sibling Git
repositories. An accountable Engine task may work across the applicable seam,
but every application consumes a hash-identified non-editable Engine wheel and
keeps its application-specific inputs, analyses, and claims locally.

- Engine owns equations, generic equilibrium compilation, exact derivatives,
  generic regression mechanics, and compact numerical evidence.
- MEA-Thermodynamics owns MEA chemistry hypotheses, source data, model
  hierarchy, fitting campaigns, validation, parameter adoption, and its
  thermodynamics manuscript.
- MEA-Absorption-Column owns absorber integration, column validation, and its
  manuscript. Its fixed-chemistry benchmark remains distinct from any future
  predictive MEA lane.
- Lithium_Extraction owns HBTA/TOPO chemistry hypotheses, source translation,
  application fitting, surrogate-domain decisions, flowsheet integration, and
  its manuscript.

The configured Better BibTeX export
`~/Zotero/exports/references.bib` is the single bibliography
authority. Repository bibliographies are one-way projections from that export;
they do not become independent citation databases.
