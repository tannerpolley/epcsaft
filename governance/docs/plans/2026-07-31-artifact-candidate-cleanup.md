# Candidate wheel and environment cleanup

Status: executed locally; publish through the Management change review.

## Actions

- Reinstalled Regression's `epcsaft` distribution from the immutable active
  object at
  `artifacts/objects/sha256/66b7ea8fb29e0a268b555cbdf401c3502517c088669a4157e8f64ab985b59ce9/`.
- Verified the installed `direct_url.json` points to that object and contains
  archive metadata rather than a mutable or temporary source path.
- Removed the two uncompressed candidate wheel aliases from the live artifact
  tree after verifying their bytes:
  - the Khudaida Figure 2 Equilibrium candidate, SHA-256
    `4f2e0bb534072f628239202dc0ae854c1d50a07989dd8b825eb3e2bd86339970`;
  - the superseded EOS native-SDK candidate, SHA-256
    `9f61f6835b2435a3d9649f6f27ebae107a46e6bafc107b263e50681033342370`.
- Preserved both exact wheel bytes in the compressed, non-installable archive
  recorded by `manifests/artifact-dispositions.json`.

## Verification

- `python3 tools/artifact_store.py check` — clean.
- `python3 tools/management_preflight.py check` — clean.
- Management tests — 22 passed.
- The active channel still contains exactly one wheel per admitted
  distribution.

The existing 2026-07-30 historical archive was not rewritten; this cleanup
uses a supplemental archive so its prior hash and receipt remain immutable.
