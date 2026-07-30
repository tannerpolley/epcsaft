# Artifact Store

Management owns the project-root `artifacts/` store policy. Package
repositories own builds, but no consumer may install another repository's
mutable `dist/` output.

## Layout

```text
artifacts/
  active/<distribution>/<wheel-name>.whl
  index.json
  objects/sha256/<full-sha256>/<wheel-name>.whl
  archive/legacy-sdk-wheels-20260730.tar.zst
```

`active/` is the only live SDK channel. It contains exactly one immutable
wheel for each admitted distribution: `epcsaft`, `epcsaft-equilibrium`, and
`epcsaft-regression`. `objects/sha256/` contains exactly one hard-linked,
full-SHA-256 object for each of those three wheels; package setup scripts may
install only those object paths. The external `artifacts/index.json` and
tracked `manifests/artifacts.json` contain the same three-entry index.

Superseded wheels are not alternate SDKs. They are retained, when needed for
scientific provenance, only inside the compressed archive recorded by
`manifests/artifact-archive.json`. The archive includes the pre-cleanup
indexes and a path/hash manifest, but is not an installation source.

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
through Management review; an artifact-store mutation without a corresponding
Management change is invalid.

New wheels must be built in disposable staging, verified, copied into the
single active alias for their distribution, normalized, and then consumed
only through the full SHA-256 object path. A version string alone is never
sufficient SDK identity. Replacing a wheel updates the one active subject;
the replaced bytes may be added to a new compressed historical archive only
when a receipt requires replay provenance.
