# Unified ePC-SAFT Engine direct cutover

Status: canonical plan; G2-G4 complete; G5 blocked

Execution note (2026-08-04): explicit G0 start approval was received. The
non-authoritative candidate and local gate evidence are recorded in
`receipts/unified-engine-g2-source-build.json` through
`receipts/unified-engine-g5-downstream.json`. G6 remains closed.

Approved architecture: 2026-08-04

Doctrine: revision 8, decision D-023

Plan owner: Governance documentation and Management coordination

Runtime owner after cutover: `ePC-SAFT/ePC-SAFT-engine`

## Purpose

Consolidate EOS, Equilibrium, and Regression into one Git repository, one
`epcsaft` distribution, one build, and one release identity. Assign one
vertical-slice task to each production capability from scientific input through
Engine implementation, installed-wheel Validation, and downstream adoption.

The cutover removes cross-wheel SDK transport, independent version selection,
Git-source application pins, and live compatibility paths. Data and Validation
remain separate because scientific input ownership and black-box evidence form
useful repository seams.

This plan does not authorize code moves, GitHub mutations, artifact promotion,
repository archival, release publication, or downstream writes. Execution
starts only after gate G0 receives a separate user instruction.

## Approved decisions

1. Rename `ePC-SAFT/ePC-SAFT-eos` to `ePC-SAFT/ePC-SAFT-engine` at cutover.
2. Transfer accepted Equilibrium and Regression implementations into Engine.
3. Build one `epcsaft` wheel for Python 3.13 or newer.
4. Publish these Python imports:
   - `epcsaft`
   - `epcsaft.equilibrium`
   - `epcsaft.regression`
5. Retire these distribution and import identities without aliases:
   - `epcsaft-equilibrium`
   - `epcsaft-regression`
   - `epcsaft_equilibrium`
   - `epcsaft_regression`
6. Keep EOS, Equilibrium, and Regression as internal modules with separate
   tests and private native targets where needed.
7. Keep Data, Validation, Lab, Governance, and downstream applications in
   separate Git repositories.
8. Give each consumer its own environment and lock it to one immutable Engine
   wheel. Do not create a shared mutable project environment.
9. Preserve historical repositories, plans, and receipts as non-executable
   provenance. Remove their names and paths from live setup, tests, builds, and
   current documentation.
10. Change runtime authority once, after every acceptance gate passes. Do not
    run a mixed active release with both the Engine wheel and the separate
    Equilibrium or Regression wheels.

## Current baseline

The baseline below records the planning inspection. Gate G1 must capture a new
machine-readable baseline before implementation.

| Subject | Planning state |
|---|---|
| EOS | `ePC-SAFT/ePC-SAFT-eos`, `main`, `d0c2575908323db6c487d4e541cf497c93c89d9e`, clean |
| Equilibrium | `ePC-SAFT/ePC-SAFT-equilibrium`, `main`, `ec73669cbb832d2309281153344e5f0dadc7b101`, clean |
| Regression | `ePC-SAFT/ePC-SAFT-regression`, `main`, `d62ea7fffe249604952b716454c0853897ba3e37`, clean |
| Validation | `codex/10-perdomo-figure3-regressed-bundle`, `a6896353c20ac97946ccd6b03e8bff1462cc3c3a`, clean |
| Data | `codex/repo-setup`, `c096285415d4d3198b9d00fc75af48b837dd1305`, dirty user work |
| Governance | `codex/governance-identity-follow-up`, `d99bdfe4e02b683d50e6d5c6c251b03fdb1f8c92` before this plan |

The current package metadata declares `0.2.0.dev0` for all three
distributions, but the setup scripts select different EOS wheel objects.
Equilibrium resolves the installed EOS CMake target. Regression transports an
include directory through `EPCSAFT_INCLUDE_DIR`. The version constraint cannot
prove that the selected native interface matches the consumer.

The current downstream application locks contain Git dependencies on the
retired `tannerpolley/ePC-SAFT` repository:

| Consumer | Current source |
|---|---|
| Lithium Extraction | `b4144c72e580eaa2d2b7ebf19c2040063f1eddd7` |
| MEA Thermodynamics | `9f51afd0f9c11a6497ddca05c8b2dd0ea0ffa785` |
| MEA Absorption Column | `9f51afd0f9c11a6497ddca05c8b2dd0ea0ffa785` |

Data contains user work that this cutover must preserve. No cutover task may
rewrite, stash, reset, or absorb that work without a Data-owned handoff.

Planning validation found these current-state failures:

- the Data checkout commit differs from the Governance packet lock;
- the current preflight mixes generic scientific packets with parameter
  ownership and reports NIST and Perdomo packets as unregistered parameters;
- declared EOS catalog packets are absent from the current EOS checkout;
- Equilibrium uses a temporary EOS wheel instead of the recorded active object;
- Regression uses a different EOS object from the recorded active object;
- the external and tracked artifact indexes are stale;
- candidate wheels remain outside the active channel; and
- mutable EOS `dist/` wheels remain in repository and worktree paths.

The planning pass must not normalize those subjects into a new three-wheel
baseline. G1 records their exact dispositions. G3 separates generic scientific
packet inventory from parameter ownership. G4 creates the replacement Engine
artifact state after one candidate cohort passes.

## Target repository and package

```text
ePC-SAFT-engine/
├── .github/workflows/
├── src/epcsaft/
│   ├── eos/
│   ├── equilibrium/
│   ├── regression/
│   └── _native/
├── tests/eos/
├── tests/equilibrium/
├── tests/regression/
├── tests/integration/
├── docs/
├── CMakeLists.txt
├── pyproject.toml
├── README.md
├── AGENTS.md
└── CONTEXT.md
```

The directory names may follow the accepted source after implementation
inventory. The public imports, one-distribution rule, and test ownership do not
vary.

One root build owns the dependency graph. Equilibrium and Regression link EOS
native targets inside that build. They do not discover an installed EOS header,
read `EPCSAFT_INCLUDE_DIR`, install an EOS provider wheel, or load a sibling
source tree.

The wheel may contain more than one private native extension when Ipopt, Ceres,
or build-time isolation makes that useful. Those extensions share one wheel
version and expose no independently installable SDK.

## Runtime interface

Downstream Python callers use one namespace:

```python
from epcsaft import Mixture, Parameters
from epcsaft.equilibrium import equilibrate
from epcsaft.regression import fit
```

Engine owns a small public interface. Backend selection flags, provider-wheel
paths, sibling repository paths, temporary directories, and artifact-store
layout do not cross that interface.

The first cutover wheel must install from a hash-addressed artifact. Each
application records the wheel version, SHA-256 digest, source commit, Python
version, and platform tag in its lock or receipt. A version string by itself
does not identify an accepted Engine build.

## Scientific data and parameter transfer

Data owns canonical packet bytes. Engine and Validation consume explicit
materializations. No installed package imports or fetches Data at runtime.

Every materialized packet record contains:

- Data repository and immutable commit;
- packet identifier, version, and path;
- full packet fingerprint;
- materialization mode and destination;
- selected files and their hashes when Engine embeds a runtime subset; and
- the task or receipt that accepted the transfer.

Engine loads parameter bytes into one typed parameter representation. A
parameter transfer receipt adds:

- component and species order;
- units and normalization conventions;
- model, association, charge, and reference-state topology;
- source and parameter fingerprints;
- fitted parameter values and covariance or diagnostic payload when present;
- Engine wheel hash and source commit; and
- the consuming EOS, Equilibrium, or Regression route.

Parameter fingerprints identify normalized scientific values. They do not
replace packet provenance. Acceptance requires both the Data packet record and
the Engine parameter fingerprint.

The source-reference transfer interface must supply every value and derivative
that an admitted Equilibrium or Regression route consumes. The cutover must do
one of these for each declared field:

- implement and validate the derivative payload; or
- remove the field and reject the capability before execution.

Empty derivative placeholders do not satisfy the contract.

## Vertical-slice ownership

One Engine task owns a production capability through these stages:

1. identify or materialize the Data packet;
2. implement the EOS behavior and derivative payload;
3. implement the Equilibrium or Regression use;
4. test the module and Engine integration interface;
5. build and verify the Engine wheel;
6. author the bounded installed-wheel Validation campaign; and
7. adapt each affected downstream application.

The Management task in the non-Git parent coordinates repository writers and
records exact subjects. Governance stores doctrine and receipts. Validation
stores evidence. An independent reviewer accepts promotion when the receipt
protocol requires one.

The same task may edit more than one repository for its slice, but it uses one
branch and commit series per repository. It verifies each repository root,
branch, status, and local instructions before writing. It does not treat the
parent coordination folder as a Git repository.

## Execution gates

### G0: plan review and start approval

Entry condition: doctrine revision 8, D-023, this plan, and
`manifests/runtime-consolidation.json` agree.

Actions:

- review the complete Governance diff;
- decide whether to commit and publish the planning change; and
- obtain a separate user instruction to start implementation.

Exit evidence: planning commit and explicit start approval bound to that
commit.

### G1: source and consumer baseline

Actions:

- fetch and record the current remote default branches without integrating
  unreviewed work;
- inventory every branch, worktree, dirty path, open migration subject, active
  wheel, SDK header, native target, public import, setup path, and live test;
- classify retained historical evidence and executable current material;
- inventory all downstream imports, dependency declarations, lock entries,
  and scientific routes;
- classify Data packets by scientific role so generic property and observation
  packets do not enter the parameter-ownership ledger;
- reconcile the Governance Data lock only after the Data owner resolves or
  hands off the dirty checkout;
- disposition each candidate, mutable `dist/`, and mismatched active wheel
  without promoting another independent package cohort;
- update active repository `AGENTS.md` and `CONTEXT.md` pointers to doctrine
  revision 8 through repository-owned changes; and
- freeze independent EOS, Equilibrium, and Regression feature work.

Stop when any source commit, dirty path, or active consumer cannot be assigned
to an owner.

Exit evidence: a machine-readable baseline receipt with exact commits, trees,
artifact hashes, worktrees, dirty-state dispositions, and consumer list.

### G2: unified source and build

Actions:

- create one Engine cutover branch from the accepted EOS source commit;
- transfer accepted Equilibrium and Regression implementation, tests, and
  current documentation with source-commit provenance;
- create one root `pyproject.toml`, CMake build, version, and Python 3.13
  environment;
- move public imports to `epcsaft.equilibrium` and
  `epcsaft.regression`;
- link EOS, Equilibrium, and Regression native targets in the same build;
- remove provider-wheel installation, installed-header discovery,
  `EPCSAFT_INCLUDE_DIR`, and temporary wheel paths; and
- delete superseded tests, setup scripts, and duplicate build files after the
  replacement tests pass.

No compatibility modules or aliases may enter the Engine tree.

Exit evidence: one source tree builds one wheel in a clean disposable staging
directory, and all Engine module tests run from one environment.

### G3: module and transfer contracts

Actions:

- replace sibling Data discovery with explicit packet materialization;
- complete the packet and parameter transfer records defined above;
- complete or narrow source-reference derivative transfer;
- add an integration slice that materializes one exact Data packet, fits a
  parameter through Regression, consumes the result through EOS and
  Equilibrium, and records the resulting fingerprints;
- remove `from_catalog`, sibling path defaults, stale fingerprints, provider
  terminology, and tests that assert retired setup behavior; and
- run static searches that fail on forbidden live names and paths.

Exit evidence: unit and integration tests prove the same typed scientific
objects cross internal module seams without filesystem or package discovery.

### G4: single-wheel validation

Actions:

- build the candidate Engine wheel from an immutable commit;
- inspect wheel metadata, RECORD entries, native libraries, public headers,
  imports, version, Python tag, and platform tag;
- install it into a new empty environment with no editable packages;
- prove all three modules originate from that wheel;
- run Engine tests against the installed artifact;
- update current Validation runners to accept one Engine wheel and one hash;
- execute EOS, Equilibrium, Regression, transfer, and cross-module campaigns;
  and
- retain failed outcomes with their exact subject hashes.

Exit evidence: one candidate wheel and a passing Validation receipt. The wheel
remains non-authoritative until G7.

### G5: downstream consumer validation

Actions for Lithium Extraction, MEA Thermodynamics, and MEA Absorption Column:

- replace the Git dependency with the immutable candidate Engine wheel source
  and SHA-256 hash;
- regenerate each lock without unrelated dependency upgrades;
- update imports to the final Engine namespace;
- remove repository-path, editable-install, and fallback setup logic;
- run the application's focused thermodynamic and integration tests; and
- record the Engine wheel and Data packet identities used by each result.

Search the Engineering workspace again before exit. Add any newly discovered
consumer to the manifest and execute the same gate.

Exit evidence: every active consumer passes in its own environment with the
same Engine wheel bytes.

### G6: external-operation preflight and approval

Prepare one operation record covering:

- rename `ePC-SAFT/ePC-SAFT-eos` to `ePC-SAFT/ePC-SAFT-engine`;
- preserve repository identity, branch protection, issues, pull requests,
  Actions settings, and redirect behavior;
- publish the accepted wheel as a versioned GitHub Release asset and mirror it
  in the Governance artifact object store;
- archive `ePC-SAFT/ePC-SAFT-equilibrium` and
  `ePC-SAFT/ePC-SAFT-regression` after downstream acceptance;
- update local checkout names and remotes;
- update organization links, workflow references, and project fields; and
- state rollback limits and recovery ownership.

Exit evidence: read-only preflight exports and explicit user approval for the
exact external operations. Doctrine approval alone does not satisfy this gate.

### G7: authority and artifact cutover

Execute the approved operation in this order:

1. merge the accepted Engine commit into the EOS repository default branch;
2. rename the repository to `ePC-SAFT-engine`;
3. build or verify the release wheel from that exact commit;
4. publish the versioned release asset and verify its downloaded SHA-256 hash;
5. replace the artifact store with one active `epcsaft` Engine alias;
6. apply reviewed downstream lock changes;
7. run the installed-wheel and downstream smoke suite from fresh environments;
8. accept the runtime authority receipt;
9. archive Equilibrium and Regression repositories; and
10. update Governance manifests from transition to completed state.

The three current wheels remain authoritative until step 8. The Engine
candidate remains non-authoritative before step 8. If a pre-receipt step fails,
stop and retain the old authority. Do not accept partial authority.

Exit evidence: one accepted Engine commit and wheel, passing consumers, updated
Governance state, archived old repositories, and a complete external-operation
receipt.

### G8: retirement and cleanup

Actions:

- remove obsolete local environments, active aliases, SDK wheel candidates,
  setup paths, ignored build output, and clean retired checkouts through
  recorded cleanup operations;
- preserve historical receipts and repository history as non-executable
  provenance;
- remove old names from live documentation, current tests, build files,
  dependency locks, and automation;
- update all active `AGENTS.md` and `CONTEXT.md` files to the completed Engine
  topology; and
- run repository cleanup audits in every changed repository.

Exit evidence: no live reference can install, import, build, or discover a
retired distribution or repository.

## Acceptance checks

The cutover cannot pass unless all checks below succeed.

### Repository and packaging

- one active runtime repository: `ePC-SAFT/ePC-SAFT-engine`;
- one active runtime distribution: `epcsaft`;
- one accepted wheel hash and source commit;
- no active `epcsaft-equilibrium` or `epcsaft-regression` distribution;
- no editable or Git-source Engine dependency in Validation or applications;
- no cross-repository Engine build input; and
- no mutable repository `dist/` wheel used as an installation source.

### Imports and native build

- `epcsaft`, `epcsaft.equilibrium`, and `epcsaft.regression` import from the
  same installed wheel;
- old import names fail;
- Equilibrium and Regression link Engine-owned native targets;
- no `EPCSAFT_INCLUDE_DIR` or installed EOS-provider discovery remains; and
- wheel inspection finds the expected native libraries and no duplicate SDK
  identity.

### Scientific transfer

- packet materialization verifies repository, commit, path, version, and full
  fingerprint;
- parameter loading verifies units, order, topology, source, and parameter
  fingerprints;
- Regression output crosses the typed Engine interface into an EOS and
  Equilibrium consumer;
- required source-reference derivatives exist and pass independent checks, or
  the unsupported capability is absent from the interface; and
- no runtime or package test discovers a sibling Data checkout.

### Validation and applications

- Engine module and integration suites pass;
- Validation installs only the immutable Engine wheel;
- current EOS, Equilibrium, Regression, and cross-module campaigns pass;
- Lithium Extraction, MEA Thermodynamics, and MEA Absorption Column pass their
  focused suites with the same wheel hash; and
- each receipt records exact commands, environment, Data inputs, Engine commit,
  and wheel hash.

### Legacy-path audit

Live files must contain none of these unless a test asserts rejection:

```text
epcsaft-equilibrium
epcsaft-regression
epcsaft_equilibrium
epcsaft_regression
EPCSAFT_INCLUDE_DIR
tannerpolley/ePC-SAFT
ePC-SAFT/ePC-SAFT-eos
ePC-SAFT/ePC-SAFT-equilibrium
ePC-SAFT/ePC-SAFT-regression
from_catalog
../ePC-SAFT-data
```

Historical plans, receipts, archived repositories, and compressed artifact
manifests are excluded from this live-path audit.

## Risks and controls

| Risk | Control |
|---|---|
| Source transfer drops accepted behavior | Capture exact source commits and run old and replacement contract tests before deleting superseded files |
| One wheel increases build time or dependency size | Keep private native targets modular, use build caching, and measure the accepted wheel |
| Python 3.13 excludes an active consumer | Inventory consumer interpreters in G1 and stop before build changes if a required consumer cannot move |
| Validation imports the Engine source checkout | Install into a fresh environment and verify module origins and wheel RECORD hashes |
| Dirty Data work is overwritten | Treat the Data checkout as externally owned and require a Data handoff before any write |
| Downstream results change during dependency migration | Record baseline outputs and compare with source-backed tolerances |
| GitHub rename or archival loses settings | Run G6 preflight, export metadata, and require exact operation approval |
| Agents create mixed compatibility paths | Assign one vertical-slice owner, ban aliases, and run the live-path audit at each gate |

## Out of scope

- merging Data, Validation, Lab, Governance, or downstream applications into
  Engine;
- preserving separate package installation identities;
- maintaining old import aliases or setup fallbacks;
- changing scientific equations without a capability-specific decision and
  evidence;
- deleting historical Git history, accepted receipts, or archived scientific
  evidence; and
- executing any external operation during planning.

## Plan completion

Planning completes when doctrine revision 8, D-023, this document,
`manifests/runtime-consolidation.json`, and the live Governance pointers agree;
the Governance diff passes local validation; and G0 remains pending.

Implementation begins only after the user reviews this plan and issues a
separate start instruction.
