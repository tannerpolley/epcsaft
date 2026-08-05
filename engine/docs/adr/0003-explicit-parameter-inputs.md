# ADR 0003: Explicit parameter inputs replace the installed catalog

- Status: accepted
- Date: 2026-07-31

## Decision

EOS accepts parameters through a strict source-neutral dictionary or an
explicit external bundle path. It ships no paper-specific parameter catalog,
does not expose catalog discovery or typed-record construction, and does not
depend on ePC-SAFT-data at runtime.

`Mixture.from_dictionary(...)` is the convenience composition entry point. It
delegates to the same strict `Parameters.from_dictionary(...)` loader and does
not create a second parameter representation or validation path.

Parameter fingerprints are computed from normalized scientific values,
component order, options, and validity. Source IDs, citations, locators, record
IDs, bundle IDs, and packet membership do not affect equation identity. Internal
native capabilities are admitted from resolved topology, required fields,
component order, and applicability—not allowlisted fingerprints.

Reusable parameter/property packets and provenance are owned by Data.
Validation may retain paper-, table-, or figure-specific extracts when they are
local evidence for a campaign. The immutable transfer inventory is recorded in
`docs/receipts/parameter-catalog-cutover.json`.

## Consequences

- A cloned EOS checkout runs its complete test suite without Data, Validation,
  network access, or a sibling source tree.
- Applications choose whether to use in-memory raw values or materialize a
  source-traceable external packet.
- Old catalog IDs fail because no catalog API exists.
- This decision does not change equations or implement the open polar leaves.
  Those leaves must report one aggregate `ares_polar` contribution.
