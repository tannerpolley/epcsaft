# Artifact Store

Management owns the project-root `artifacts/` store policy. Package
repositories own builds, but no consumer may install another repository's
mutable `dist/` output.

## Layout

```text
artifacts/
  index.json
  objects/sha256/<full-sha256>/<wheel-name>.whl
  <retained evidence aliases>
```

`objects/sha256/` is the immutable installation surface. Existing
receipt-bound paths remain as evidence aliases until Management records a
specific disposition. The external `artifacts/index.json` and tracked
`manifests/artifacts.json` contain the same index, map every alias to one
object, and assign a retention class.

Normalize after intentionally adding a wheel:

```bash
python tools/artifact_store.py normalize
python tools/artifact_store.py check
```

The check fails when:

- a wheel exists in a repository `dist/`;
- an alias is missing from the index;
- an object is absent or has the wrong hash;
- a current admitted subject is absent; or
- the retention policy or index is stale.

Normalization intentionally updates the tracked index. Commit that diff
through Management review; an artifact-store mutation without a corresponding
Management change is invalid.

New wheels must be built in disposable staging, verified, copied into an
evidence alias when required, normalized, and then consumed only through the
full SHA-256 object path. A version string alone is never sufficient SDK
identity.
