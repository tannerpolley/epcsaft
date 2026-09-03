---
status: accepted
date: 2026-09-01
---

# Retire Governance and Validation as directory owners

Kind: computational-architecture

## Context

The former `governance/` tree duplicated current execution rules, architecture,
scientific decisions, and GitHub state. The former `validation/` tree mixed
scientific inputs, scripts, results, figures, checks, and administrative files
under one miscellaneous owner. Git history already preserves their obsolete
plans and administrative records.

## Decision

Current repository execution rules live in `AGENTS.md`; executable topology in
`engine/ARCHITECTURE.yaml`; scientific vocabulary and boundaries in
`docs/scientific/CONTEXT.md`; and accepted decisions in this sole ADR
sequence. Data owns reusable source packets and `data/packet-locks.json`.
Scientific work is retained as one-question analysis sets under `analyses/`.
Reproductions use the documented local extension
`analyses/reproduction/<short-id>/`; no alias, redirect, compatibility tree, or
second decision ledger survives.

The former decision identifiers map as follows:

- D-001 through D-028 are superseded topology, administration, or migration
  history except the still-live CppAD, derivative, evidence, and one-wheel
  principles now stated by ADRs 0003--0005, `AGENTS.md`, and
  `engine/ARCHITECTURE.yaml`.
- D-029 and D-030 remain represented by ADRs 0001 and 0002 and the retained
  CAP-09 Khudaida result.
- D-031, D-034, D-035, D-036, D-037, and D-042 retain their exact observation,
  topology, derivative, failure, and claim limits through ADRs 0004 and 0005,
  Engine source, and the surviving owners below. Their old identifiers remain
  historical locators, not a second authority sequence.
- D-032 and D-033 remain represented by ADR 0005 and the build-first fast root
  checks named for CAP-11 above.
- D-038 remains the named downstream-application rule in root `AGENTS.md`.
  D-039 is superseded by D-042; D-040 is complete; D-041 remains represented by
  ADR 0005 and the sole `NativeModel` EOS owner.

The current evidence owners for the former capability claims are:

- CAP-11: `tools/check-analysis-layout.py`, `tools/check-all.sh`, and the
  current architecture/documentation files; no standalone administrative
  result is an active owner.
- CAP-12: the homogeneous row in
  `analyses/2026-grepe-hard-cut/results/validation-gates.csv`.
- CAP-13: the reactive-bubble row in
  `analyses/2026-grepe-hard-cut/results/validation-gates.csv`.
- CAP-14: `engine/ARCHITECTURE.yaml` and
  `engine/tests/equilibrium/numerical/test_grepe.py`.

The 1998 Chase validation ladder is aggregated by
`analyses/reproduction/1998-chase/results/chemical-equilibrium-validation.json`;
its detailed N2O4/NO2 and water-gas-shift results are owned by
`analyses/reproduction/1998-chase/results/n2o4-no2-equilibrium.json` and
`analyses/reproduction/1998-chase/results/water-gas-shift-equilibrium.json`.
The corresponding IAPWS and Belov result owners are
`analyses/reproduction/2019-iapws/results/water-self-ionization.json` and
`analyses/reproduction/2023-belov/results/test-01-equilibrium.json`.
Former per-case `result.json`, receipt, and migration paths are historical
provenance only and are not promised as current evidence owners.

## Consequences and validity boundary

The cutover changes file ownership and active paths only. It changes no
equation, parameter value, unit, basis, standard state, solver control,
tolerance, public callable, numerical result, failure, or scientific claim.
Historical repository names and old paths inside retained results remain
provenance. Compact Engine checks and explicit analysis commands keep their
distinct failure meanings; repository readiness does not establish scientific
correctness.
