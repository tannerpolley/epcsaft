# EOS Naming, SDK, Artifact, and Worktree Normalization Plan

Status: planned; execution blocked on ecosystem quiescence

Plan owner: Organization policy, with Migration controlling execution

Created: 2026-07-29

Doctrine baseline: revision 4

Target doctrine: revision 5

## 1. Goal

Coordinate one ecosystem migration that:

- canonizes **EOS** as the repository role and removes the legacy `Provider`
  role label from every live authoritative surface;
- renames the GitHub repository to `ePC-SAFT/ePC-SAFT-eos`;
- preserves the Python distribution and import identity `epcsaft`;
- removes legacy role terminology from live Python, C, C++, serialized,
  command-line, status, documentation, tracker, and worktree contracts;
- limits the installed EOS parameter catalog to an explicit allowlist of
  model-defining PC-SAFT and ePC-SAFT reference parameterizations;
- relocates case-study, fitted, alternative, and paper-reproduction parameter
  bundles to their accountable Validation campaigns;
- normalizes native SDK ownership and consumer checks;
- makes wheel and installed-environment identity reproducible;
- classifies and preserves scientific evidence artifacts without treating
  accumulated development wheels as current releases; and
- installs the repository-owned worktree access hierarchy selected by the
  user.

This is a living execution plan. It does not authorize implementation,
cleanup, publication, or a GitHub rename by itself.

## 2. Success contract

The migration is complete only when all of the following are true:

1. The repository is named `ePC-SAFT/ePC-SAFT-eos`, while its distribution and
   import remain `epcsaft`.
2. Live code, public APIs, schemas, filenames, command-line flags, statuses,
   worktree labels, canonical documentation, and active tracker surfaces use
   EOS terminology.
3. No compatibility aliases retain the retired terminology.
4. Every remaining case-insensitive match for the retired term belongs to an
   explicit historical-evidence allowlist.
5. Accepted receipts, immutable result identities, and Git history remain
   truthful and replayable.
6. EOS remains the sole owner of its thermodynamic kernel, native capsule,
   public headers, capability vocabulary, and any admitted standalone native
   model transport.
7. Equilibrium and Regression compile and run against the exact EOS wheel and
   header recorded by their build subject.
8. No consumer installs EOS from a mutable sibling `dist/` path.
9. Every retained wheel has an owner, exact identity, evidence reference, and
   retention class.
10. Development environments are disposable, reproducible, and free of
    unrelated sibling packages.
11. Cross-repository worktrees are owned by the repository being changed and
    implement the approved access matrix.
12. The EOS distribution contains only allowlisted model-defining reference
    parameterizations; every other retained parameter bundle is evidence data
    in an accountable Validation campaign and enters EOS only through a public
    parameter-input contract.
13. Validation parameter bundles are never imported by production packages,
    exposed as a production catalog, or treated as runtime authority.
14. Package, installed-wheel, cross-package, scientific, negative-space, and
    historical-evidence checks pass.
15. GitHub repository identity, redirects, branch protection, issues, pull
    requests, milestones, and organization Project references are verified
    after the rename.

## 3. Decisions fixed by this plan

### 3.1 Names

- GitHub repository: `ePC-SAFT/ePC-SAFT-eos`
- Python distribution: `epcsaft`
- Python import: `epcsaft`
- Canonical repository role: `EOS`
- Equilibrium distribution/import:
  `epcsaft-equilibrium` / `epcsaft_equilibrium`
- Regression distribution/import:
  `epcsaft-regression` / `epcsaft_regression`

The repository name and Python package identity are independent. Renaming the
repository does not require renaming the distribution or import.

### 3.2 Purge boundary

The target is a zero-match live surface, not rewritten history.

The live surface includes:

- public and private source identifiers;
- exported types, fields, functions, methods, and constants;
- current serialized keys and status vocabularies;
- CLI flags and help;
- active tests and fixtures, except fixtures whose bytes are historical
  evidence;
- current filenames and directory names;
- `AGENTS.md`, `CONTEXT.md`, READMEs, current ADRs, current plans, active
  manuscript sources, and generated current documentation;
- active worktree names and routing tables; and
- open tracker titles, descriptions, labels, milestones, and Project fields.

The historical allowlist may contain only:

- accepted receipt IDs and receipt bodies;
- receipt-bound artifact, capability, result, and path identities;
- immutable historical result files whose hashes or schemas require the
  original spelling;
- exact historical repository identities;
- historical Git commits, pull requests, closed issues, and comments; and
- the minimum validator fixtures required to replay those records.

Historical evidence must not be presented as current naming authority.

### 3.3 Compatibility

This is a breaking cleanup. Do not add aliases, duplicate fields, dual-read
serialization, deprecated spellings, forwarding wrappers, or alternate CLI
flags.

### 3.4 Package and artifact identity

Package version identifies a compatible contract line. Commit, tree, wheel
SHA-256, and embedded-header SHA-256 identify the exact executable subject.

The provisional coordinated development line after the breaking rename is:

- `epcsaft 0.3.0.dev0`;
- `epcsaft-equilibrium 0.3.0.dev0`; and
- `epcsaft-regression 0.3.0.dev0`.

The refresh gate must confirm or amend this line before edits begin.
Experimental builds do not require a new package version for every commit, but
they must be stored at immutable commit- and hash-addressed paths. A path
already consumed by another environment or campaign must never be overwritten.

### 3.5 Parameter-bundle ownership

EOS owns generic parameter schemas, validation, resolution, and evaluation. It
may ship only parameterizations admitted by an explicit, reviewed allowlist of
model-defining PC-SAFT or ePC-SAFT references. The allowlist is based on a
source's architectural role, not merely on publication in a paper. Initial
classification examples selected by the user are Gross and Sadowski (2001)
and Figiel et al. (2025); the migration inventory must verify the exact records
and editions before admission.

Validation owns parameter values used for:

- a substance-specific or application-specific case study;
- a literature reproduction or comparison campaign;
- a regression candidate, fitted result, sensitivity study, or alternative
  parameterization;
- a rejected, superseded, or non-admission result; or
- any published parameter set not admitted to the model-defining allowlist.

Such values are test and evidence data. They must be source-traceable,
hash-bound, scoped to a campaign, and supplied through EOS's public
parameter-input API. They must not be installed as an EOS catalog entry,
imported by Equilibrium or Regression as package data, or promoted into a
Validation runtime catalog.

Regression may own a transient or committed fit result while it is being
produced and diagnosed. Durable cross-package or literature acceptance
evidence, including the exact accepted candidate values, belongs to the
accountable Validation campaign. Promotion to the EOS catalog requires a
separate Organization decision that adds the source and exact records to the
model-defining allowlist; successful validation alone is insufficient.

The just-landed Crespo et al. DEEA subject is preserved as an exact migration
input. Its generic ordinary-`sigma^3 kappa` EOS and SDK capability remains an
EOS responsibility, while its substance-specific DEEA catalog bundle is
scheduled for relocation into its Validation campaign. Existing commits,
wheels, and receipts remain historical evidence and are not rewritten.

Implementation is tracked by EOS issue
[`ePC-SAFT/ePC-SAFT#48`](https://github.com/ePC-SAFT/ePC-SAFT/issues/48)
and Validation issue
[`ePC-SAFT/ePC-SAFT-validation#9`](https://github.com/ePC-SAFT/ePC-SAFT-validation/issues/9).

## 4. Current audit snapshot

This snapshot explains the plan but is not an execution baseline. Refresh it
after current work lands.

### 4.1 Wheels

The 2026-07-29 snapshot contains 105 wheels totaling about 75.9 MB, including:

- 103 wheels in the non-Git sibling artifact store;
- one mutable EOS `dist/` wheel; and
- one Lab experiment wheel.

The artifact store contains about 70 MB and has no complete manifest,
current pointer, retention class, or owner/reference index.

The active Validation 0.2 subjects exist and match their documented hashes:

| Distribution | Commit | Wheel SHA-256 |
| --- | --- | --- |
| `epcsaft` | `14fa3745264db66b8e59c12268737d694c706f2f` | `48f3a75c9fc16ba71616aa703b526f41c2dcf89a7e00eebe23f75fcb8fa24594` |
| `epcsaft-equilibrium` | `1213ddcde31a55a0cd36e667dea72933fd132781` | `040bc759e52b4be2e707f1da5da4271c517418f1f0ca62401f0e1ef3292b8a99` |
| `epcsaft-regression` | `f6982a8f0661555dfe91c5be38aec118e37ab83c` | `c03dc63304c228df0f9f0e288e5fae3b744ca3e1486089b94b3c8a52e8dd8b38` |

One accepted historical wheel from commit
`a10cdb27566dbe5aa622f684379f577a34fc1cdb` is absent locally. Its receipt
states that it is rebuildable. Rebuilding it must create a new state record;
the accepted receipt must not be changed.

### 4.2 Version and header collision

Different EOS wheel and SDK-header subjects share the same `0.2.0.dev0`
distribution version:

| Subject | Identity |
| --- | --- |
| Current EOS source header | `4de9acf02639f67e7a83267c5adac96104d53e76e8caa8c5ba678ff235d18688` |
| Current-commit candidate wheel | `6b02e36c8b6c090d44e6933909214a68ecfb5913d04e0cdc59547ffa794dde25` |
| Mutable EOS `dist/` wheel | `4cee10a9158576307cda93f611b6ade3a7cf8819df44f83efe8cbc61ab038789` |
| Header embedded in mutable `dist/` wheel | `bb46f103b7116efcc31f000cd07082f2d79c36a3236265d3e1d9ab4d27cc733f` |
| Header installed in consumer environments | `2cd2b73b83c65936dff21155fd800a87b56e81cce977df7b8491ccfb2bf4c50b` |

The Equilibrium and Regression environments point to the mutable `dist/`
filename, but their installed headers differ from the bytes now at that path.
Dependency checks still pass because all subjects report
`0.2.0.dev0`.

### 4.3 Native contracts

The current ecosystem has:

1. the public EOS capsule `epcsaft.native_sdk.v1`, ABI 1, with an append-only,
   size-gated function table and size-tagged results;
2. a public standalone EOS native-model header, static library, and CMake
   transport not yet reconciled with ADR 0001; and
3. a private Regression evaluator capsule
   `epcsaft.regression.evaluator.v1`, ABI 1.

The EOS capsule's prefix and size checks are a sound compatibility design.
The principal risks are indistinguishable wheel subjects, an undocumented
second transport, static-library linkage where a header-only consumer target
would suffice, and repeated consumer validation logic.

### 4.4 Environments

- EOS is installed editable in its repository environment.
- The EOS environment also contains Regression from a cross-repository
  worktree, which is not part of the EOS environment's ownership.
- Equilibrium is editable and consumes an EOS wheel from the mutable EOS
  `dist/` path.
- Regression is editable and consumes an EOS wheel from the same mutable
  path.
- System Python contains none of the three distributions.
- Validation uses isolated installed-wheel environments rather than a
  permanent resident environment.

## 5. Update protocol for incoming work

### 5.1 Landing register

Before execution, create a register with one row for every open pull request,
active branch, dirty checkout, linked worktree, candidate wheel, and pending
receipt:

| Field | Meaning |
| --- | --- |
| Repository | Owning Git repository |
| Subject | PR, branch, worktree, wheel, receipt, or local change |
| Exact identity | Full commit/tree, path, or SHA-256 |
| Owner | Agent or human responsible for disposition |
| Status | Active, merged, retained, superseded, or abandoned |
| Touched surface | API, SDK, build, docs, evidence, tracker, or worktree |
| Plan effect | Rename mapping, artifact ledger, allowlist, rebuild, or none |

### 5.2 Refresh rule

Run refresh gate R0 after all current agents report completion and again
before each major phase. Compare the new inventory with the preceding
baseline.

- A new live API or SDK field extends the rename and compatibility map.
- A new wheel extends the artifact ledger.
- A new or changed parameter bundle extends the parameter-ownership ledger.
- A new accepted receipt extends the immutable historical allowlist.
- A new current document joins the live purge.
- A superseded experiment receives an evidence disposition before removal.
- An architecture conflict stops execution until the plan and decision ledger
  are amended.

Do not rebase, rewrite, clean, or absorb active work merely to fit an older
version of this plan.

## 6. Execution ownership and concurrency

- Migration owns sequencing, exact-subject binding, the landing register,
  artifact ledger, external-write gates, and state receipts.
- Organization owns doctrine and shared policy.
- EOS owns its runtime, package, parameters, native SDK, headers, build
  metadata, and assigned Validation campaigns.
- Equilibrium and Regression remain top operational agents and may work across
  repositories through target-owned access worktrees.
- Validation is an evidence home with no resident agent and no outbound
  access.
- EOS may work in its own repository and Validation only.
- One primary Equilibrium or Regression executor owns migration writes one
  repository at a time.
- A distinct reviewer acts only on stable commits.
- No two writers mutate the same repository or artifact subject concurrently.

## 7. Phased execution

### Phase 0: Wait and preserve

1. Allow active package, Validation, Lab, manuscript, and pull-request work to
   finish.
2. Record every remaining dirty or untracked path.
3. Preserve unfinished work without automatic stashing, rebasing, cleaning,
   or branch deletion.
4. Stop package builds and Validation writers.

Exit gate:

- every active subject is merged, retained at an exact commit, or excluded by
  an owner decision;
- no relevant writer or build remains active; and
- no target checkout contains unexplained changes.

### Phase 1: Refresh gate R0

Capture:

- repository HEADs, trees, branches, remotes, protections, and worktrees;
- open issues, pull requests, milestones, and Project items;
- public APIs and serialized contracts;
- all native SDK declarations, sizes, fields, callbacks, and consumers;
- package metadata, CMake exports, build floors, and dependency constraints;
- every installed and campaign-local parameter bundle, its source, exact
  records, consumers, and proposed ownership class;
- all wheel paths, hashes, embedded headers, and producing commits;
- all environment distributions, origins, native modules, and headers; and
- a fresh case-insensitive terminology inventory.

Update this plan, the landing register, and the proposed version line before
implementation.

### Phase 2: Governance and decisions

1. Adopt Doctrine revision 5.
2. Add Organization and Migration decisions for:
   - EOS vocabulary and repository identity;
   - the historical evidence boundary;
   - breaking API treatment with no aliases;
   - cross-repository access;
   - artifact identity and retention;
   - the model-defining parameter allowlist and Validation evidence boundary;
   - SDK transport ownership; and
   - the exact GitHub cutover protocol.
3. Update every repository's doctrine pointer.

Exit gate:

- user approval binds the exact doctrine amendment;
- Migration records the implementation and external-write gates; and
- no runtime or remote change has occurred.

### Phase 3: Worktree topology

After all affected worktrees are clean and inactive:

1. inspect every worktree's `.git` pointer and common Git directory;
2. retain or retire each old worktree according to reachability and approved
   disposition;
3. create target-owned access worktrees:

```text
ePC-SAFT-eos/.worktrees/equilibrium-access
ePC-SAFT-eos/.worktrees/regression-access
ePC-SAFT-equilibrium/.worktrees/regression-access
ePC-SAFT-regression/.worktrees/equilibrium-access
ePC-SAFT-validation/.worktrees/equilibrium-access
ePC-SAFT-validation/.worktrees/regression-access
ePC-SAFT-validation/.worktrees/eos-access
```

4. verify that each `.git` file points to the target repository's common Git
   directory and remote;
5. verify the access matrix rejects prohibited origins.

Worktree creation, removal, and branch cleanup remain separate Git mutations
with exact-target checks.

### Phase 4: Artifact ledger

Generate one machine-readable artifact inventory keyed by SHA-256. Each record
contains:

- distribution, version, Python ABI, and platform;
- full commit and tree;
- wheel SHA-256;
- embedded public-header SHA-256;
- producing repository and command when known;
- every current path;
- receipt, result, and plan references;
- current versus historical status;
- rebuildability; and
- one retention class.

Retention classes are:

- current development subject;
- active Validation subject;
- accepted receipt-bound evidence;
- retained failed or non-admission evidence;
- stable non-production checkpoint;
- superseded replay;
- rebuildable missing artifact; and
- unreferenced candidate pending review.

Do not move or delete artifacts in this phase.

### Phase 5: Live terminology and API migration

Prepare a reviewed symbol and path map. Expected changes include:

- public phase types;
- public result and diagnostic fields;
- native context and metadata identifiers;
- current serialization keys;
- CLI flags;
- status strings;
- source and test filenames;
- worktree and routing labels;
- current package documentation;
- active Validation documentation;
- active Lab material;
- Migration's living controls; and
- Organization policy and manuscript sources.

Prepare a parallel parameter-bundle disposition map with one row per bundle:

- exact source and record hash;
- current repository and consumers;
- model-defining allowlist decision;
- target EOS catalog entry or Validation campaign;
- Regression provenance when fitted;
- receipt and artifact references; and
- migration action with replay test.

Apply the migration in this order:

1. Organization doctrine and Migration controls;
2. EOS;
3. Equilibrium;
4. Regression;
5. Validation;
6. active Lab material; and
7. remaining Organization documentation and manuscript files.

Do not add aliases or dual-read compatibility.

### Phase 6: SDK ownership and packaging

1. Keep `epcsaft.native_sdk.v1` when append-only ABI compatibility remains
   proven.
2. Retain generic parameter parsing, validation, resolution, and evaluation in
   EOS.
3. Generate the EOS installed catalog from the reviewed model-defining
   allowlist only.
4. Move non-allowlisted bundles, beginning with the Crespo DEEA case, into
   accountable Validation campaigns without changing their scientific bytes
   or provenance.
5. Make every migrated campaign construct or load public parameter inputs
   explicitly; no production package may import Validation data.
6. Rename live legacy-role source identifiers and metadata without needless
   binary-layout changes.
7. Make EOS the sole owner of its public headers, capsule contract,
   capability vocabulary, and transport metadata.
8. Adopt or supersede an ADR for the standalone native-model transport.
9. Export:
   - `epcsaft::headers` as an interface target for capsule consumers; and
   - `epcsaft::native_model` for the standalone native-model implementation.
10. Stop Equilibrium linking the implementation archive merely to obtain
   declarations.
11. Centralize Regression's EOS capsule validation in one C++ owner.
12. Keep the Regression evaluator capsule private unless another package must
   produce it. If it crosses a package boundary, move the declaration to its
   producer and assign a contract version.
13. Align or justify Python, CMake, scikit-build-core, pybind11, and C++ floors.

### Phase 7: Coordinated package line and immutable wheels

1. Confirm the coordinated development version.
2. Build EOS from an exact clean commit.
3. Store the EOS wheel at an immutable commit- and hash-addressed path.
4. Record the embedded header hash.
5. Build Equilibrium and Regression against that exact installed EOS wheel.
6. Store their wheels immutably.
7. Generate one exact coordinated artifact subject containing all three
   wheels, commits, trees, hashes, header identities, and dependency
   constraints.

Mutable `dist/` directories remain disposable package-local output and are
never cross-repository dependency sources.

### Phase 8: Environment recreation

After explicit cleanup approval:

1. replace each repository environment atomically from its exact subject;
2. install only the local editable owner and exact immutable dependencies;
3. remove unrelated sibling distributions;
4. verify:
   - installed distribution version;
   - direct artifact origin;
   - wheel SHA-256;
   - embedded and installed header equality;
   - native-module origin;
   - capsule name and ABI;
   - required table terminal field and result sizes;
   - dependency compatibility; and
   - absence of sibling-source imports and include paths.

Validation continues to create temporary isolated environments for each exact
campaign.

### Phase 9: Verification

Run the narrowest checks first, then broaden:

1. package build and discrete API tests;
2. SDK layout, prefix, offset, size, malformed-input, and capability tests;
3. installed-wheel tests;
4. EOS-to-Equilibrium consumer tests;
5. EOS-to-Regression consumer tests;
6. serialization and error-contract tests;
7. active Validation campaigns;
8. solver, numerical, physical, source, and claim-specific scientific gates;
9. clean-clone and offline exact-wheel replay;
10. import-origin and forbidden-dependency checks;
11. artifact-ledger completeness;
12. worktree ownership and access enforcement;
13. historical receipt replay; and
14. parameter-catalog allowlist and ownership negative-space checks; and
15. terminology negative-space checks.

Terminology verification requires:

- zero matches in the declared live surface;
- every residual match to be present in the reviewed historical allowlist;
- no old-term public alias;
- no old-term current filename or worktree label; and
- no old-term open tracker title, label, milestone, or Project field.

Parameter-ownership verification requires:

- every EOS catalog record to map to an approved model-defining allowlist
  entry;
- every non-allowlisted retained bundle to map to one accountable Validation
  campaign;
- no Validation parameter bundle to be imported by a production package;
- installed-wheel tests to prove EOS functions with explicit campaign
  parameters; and
- exact pre-migration bundle bytes and provenance to remain reachable.

### Phase 10: Independent review

A reviewer distinct from the implementation writer checks:

- completeness of the API and path migration;
- binary and source compatibility claims;
- artifact and environment identity;
- evidence preservation;
- SDK ownership and dependency direction;
- worktree access;
- absence of duplicate EOS ownership; and
- GitHub cutover readiness.

Any material correction returns to the owning phase and invalidates downstream
subjects that depend on it.

### Phase 11: GitHub cutover

At a separate external-write approval gate:

1. capture immutable machine-readable exports for repository metadata at risk;
2. capture issue, pull request, milestone, Project, branch, protection,
   Actions, release, package, deployment, environment, Pages, ruleset, and
   webhook state;
3. verify the destination name remains available;
4. bind the exact repository node ID, numeric ID, default-branch commit, and
   protection state;
5. submit the exact repository rename;
6. verify the repository object retains its identity;
7. verify issues, pull requests, milestones, Project references, branches, and
   protection;
8. verify old web and Git URLs redirect;
9. update all local remotes;
10. update active tracker and Project terminology;
11. reserve the old repository name permanently; and
12. write a post-state receipt.

Do not rely on redirects for Actions or reusable workflows. Do not reuse the
old repository name.

### Phase 12: Publication and cleanup

Publish in dependency and authority order:

1. Organization and Migration authority;
2. EOS;
3. Equilibrium and Regression;
4. Validation;
5. active Lab cleanup;
6. remaining Organization documentation; and
7. Migration state records.

Only after every commit is reachable and verified:

- remove approved unreferenced candidate artifacts;
- remove stale mutable package build output;
- retain receipt-bound and non-admission evidence;
- remove only clean merged worktrees or worktrees approved as abandoned;
- verify canonical checkouts are clean and synchronized; and
- run each repository's cleanup audit.

This plan must then be rewritten or distilled so that it no longer presents
the retired role term as live vocabulary. Its retained pre-cutover revision is
historical planning evidence.

## 8. Stop conditions

Stop the affected phase when:

- a repository, worktree, artifact, public API, SDK layout, or receipt changes
  after its baseline;
- an active agent or writer reappears;
- a dirty path lacks an owner or disposition;
- a source/header/wheel/installed identity differs unexpectedly;
- an accepted receipt would require rewriting;
- a parameter source cannot be classified reproducibly as model-defining or
  campaign evidence;
- a Validation bundle is required as a production runtime dependency;
- the historical allowlist expands without review;
- a consumer needs a second EOS implementation or sibling source;
- a build requires a mutable cross-repository `dist/` path;
- a scientific or numerical gate fails;
- the GitHub target name is occupied;
- remote protection or metadata differs from the cutover record; or
- recovery would require force push, history rewrite, deletion, or an
  unapproved external mutation.

## 9. Recovery principles

- Preserve source commits and exact artifact hashes before mutation.
- Prefer new corrective commits and state records over rewriting accepted
  history.
- A repository rename rollback does not guarantee restoration of redirects or
  metadata behavior.
- Do not reuse the old GitHub repository name.
- Do not delete an artifact until its retention class and evidence references
  are proven.
- Do not remove a worktree until it is clean and its branch disposition is
  approved.
- Do not treat version equality as proof of SDK or scientific-subject
  equality.

## 10. Execution entry point

When the user reports that the other agents are finished:

1. do not begin the refactor;
2. run Phase 0 quiescence verification;
3. run refresh gate R0;
4. update this document and the task plan with newly landed work;
5. present material changes to assumptions and the exact first implementation
   subject; and
6. begin only after the required authority and external-write gates are
   satisfied.
