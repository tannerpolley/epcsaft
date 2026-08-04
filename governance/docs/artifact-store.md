# Artifact Store

Governance owns the project-root `artifacts/` store policy. Package
repositories own builds, but no consumer may install another repository's
mutable `dist/` output.

## Layout

```text
artifacts/
  active/<distribution>/<wheel-name>.whl
  index.json
  objects/sha256/<full-sha256>/<wheel-name>.whl
  archive/legacy-sdk-wheels-20260730.tar.zst
  archive/candidate-wheel-dispositions-20260731.tar.zst
  archive/candidate-wheel-dispositions-20260804.tar.zst
```

`active/` is the only live SDK channel. It contains exactly one immutable
wheel for each admitted distribution: `epcsaft`, `epcsaft-equilibrium`, and
`epcsaft-regression`. `objects/sha256/` contains the hard-linked,
full-SHA-256 object for every indexed active wheel and may retain immutable
objects that have no live alias. Package setup scripts may install only the
active subjects' object paths. The external `artifacts/index.json` and
tracked `manifests/artifacts.json` contain the same three-entry index.

Consumers resolve an admitted distribution through Governance rather than
embedding an alias, filename, or hash in repository-local setup code:

```bash
python tools/artifact_store.py resolve --distribution epcsaft
```

The command checks the complete store policy and prints the absolute,
full-SHA-256 object path. Unknown distributions or an inconsistent store fail
closed without returning an install path.

## Approved Engine cutover

D-023 replaces the three-distribution channel with one `epcsaft` Engine wheel.
The current three entries remain the live channel until the direct-cutover
receipt passes every gate in
`docs/plans/2026-08-04-unified-engine-direct-cutover.md`. Governance MUST NOT
promote an Engine wheel beside the three current distributions or update one
current distribution as a partial cutover.

Correcting a stale current SDK binding to the exact source-bound wheel already
required by a current package is maintenance of the three-distribution channel,
not an Engine cutover. Such a repair still updates the policy, active alias,
tracked index, package setup scripts, and validation evidence as one operation.

At cutover, Governance will:

- admit one Engine wheel built from the accepted Engine commit;
- verify that the wheel contains the EOS, Equilibrium, and Regression modules;
- replace the current `epcsaft` alias with that wheel;
- remove the `epcsaft-equilibrium` and `epcsaft-regression` active aliases
  through recorded dispositions; and
- regenerate the external and tracked artifact indexes as one-entry indexes.

Downstream applications install the full-SHA-256 Engine object into their own
environments. They do not share a mutable project environment or install from
a Git repository.

Superseded wheels are not alternate SDKs. They are retained, when needed for
scientific provenance, only inside the compressed archive recorded by
`manifests/artifact-archive.json`. The archive includes the pre-cleanup
indexes and a path/hash manifest, but is not an installation source.
Later candidate dispositions use a separate compressed archive and are
recorded in `manifests/artifact-dispositions.json`; those archives are also
never installation sources.

Normalize after intentionally adding a wheel:

```bash
python tools/artifact_store.py normalize
python tools/artifact_store.py check
```

The check fails when:

- a wheel exists in a repository `dist/`;
- an alias is missing from the index;
- a wheel exists outside `active/` (or an uncompressed wheel appears in the
  historical archive);
- an object is absent or has the wrong hash;
- more than one active wheel exists for a distribution;
- an active distribution is absent; or
- the retention policy or index is stale.

Normalization intentionally updates the tracked index. Commit that diff
through Governance review; an artifact-store mutation without a corresponding
Governance change is invalid.

New wheels must be built in disposable staging, verified, copied into the
single active alias for their distribution, normalized, and then consumed
only through the full SHA-256 object path. A version string alone is never
sufficient SDK identity. Replacing a wheel updates the one active subject;
the replaced bytes may be added to a new compressed historical archive only
when a receipt requires replay provenance.
