# Scientific packet contract

Status: canonical

Version: 3

## Purpose

A packet is the smallest reusable, reviewed scientific input with enough
metadata to interpret, audit, hash, and reuse it without relying on a package
source checkout. Packet storage does not admit runtime capability or validation
evidence.

## Layout

Packets live at `packets/<packet-id>/<version>/packet.toml`. Optional content
appears only when needed under `parameters/`, `observations/`, `sources/`, or
`provenance/`. Empty placeholder directories and unlisted files are invalid.
There is no parallel top-level `species/`, `properties/` or mutable catalog
directory: packets are the unit of review, hashing, versioning and consumer
locking. A packet may hold several species that share a source product and
license boundary; split packets when source identity, versioning or use basis
differs.

Full papers, broad extraction corpora, rejected experiments, and exploratory
notebooks are outside Data. Generated results and plots remain in
their owning analysis set.

## Property-data admission boundary

An observation file belongs in Data only when it is a reusable, source-backed
property record that can be consumed independently of one named validation
campaign. It should retain the component identity, independent conditions,
explicit units, source locator or source ID, and source uncertainty when the
source provides it. Source-native wide tables are acceptable when their column
names and packet metadata make the units and meaning unambiguous. Observation
files preserve source values as decimal strings and retain source row roles.
Missing source uncertainty stays missing; it is never replaced with zero or
turned into a model-acceptance threshold. Model predictions and fitted values
are not observations.

Paper- or figure-specific matrices, digitized curves, direct experimental rows
kept only to reproduce a named publication, fitted targets, campaign-specific
partitions, and generated results remain under `analyses/`. Normalizing a table or
copying it into CSV does not change that ownership. Promotion of a campaign
record requires a separate review and a new packet version or family; existing
immutable packets and consumer locks are never rewritten in place.

## Manifest

`packet.toml` requires:

- `schema_version = 2` for current packets. Schema version 1 is retained for
  accepted historical packets and does not gain fields retroactively;
- `packet_id`, matching the packet directory;
- `packet_version`, matching the version directory;
- a short `title` and `description`;
- `license_or_use_basis`, including any redistribution restriction;
- `tolerance_rationale`, explaining the retained source uncertainty or stating
  explicitly why no numerical acceptance tolerance is claimed. This field is
  required for schema version 2 packets;
- one or more `[[sources]]` entries with `source_id`, `citation`, and
  `use_basis`;
- one or more `[[domains]]` entries with `domain_id` and `description`; and
- one `[[files]]` entry for every retained file.

Each file entry requires:

- a safe packet-relative `path` under an allowed content directory;
- `role`: `parameters`, `observations`, `source`, `provenance`, or `schema`;
- `media_type`;
- lowercase SHA-256 `sha256` of the exact bytes;
- `source_refs`, naming manifest source IDs;
- `domain_refs`, naming manifest domain IDs; and
- `units_basis`, explaining units or explicitly stating why the content is
  dimensionless or non-numeric.

The manifest lists every retained file exactly once. Paths may not be absolute,
escape the packet, or traverse symlinks. Hash changes require a new commit; a
semantic or source change should normally create a new packet version.

Provenance manifests under `provenance/` are TOML and are validated with the
packet. `repository-migration-v1` is for a packet copied from one immutable
repository source and records one source bundle manifest. Mixed-source
property packets use `source-record-provenance-v1`; each retained non-
provenance file then records its own origin repository, commit, tree, path, and
SHA-256. This distinction prevents a single synthetic origin record from
concealing a cross-repository migration.

Project-derived parameter documents use `derived-record-provenance-v1`. The
source record hashes every derived packet file, names the transformations, and hashes
each retained repository input so mixed-source fits remain reproducible without
claiming a false single-repository origin.

Origin identities are kept exactly as recorded. On 2026-10-22 the repositories cited as `ePC-SAFT/<name>` moved to
`tannerpolley/<name>`, keeping their commits, issues and releases. The archived `ePC-SAFT-eos` is
`tannerpolley/ePC-SAFT-engine` (its renamed successor). A recorded commit therefore resolves at
`github.com/tannerpolley/<name>/commit/<sha>`, for example
`ePC-SAFT/ePC-SAFT-validation@dc1841fc` at `tannerpolley/ePC-SAFT-validation`.

## Consumer lock

The live `data/packet-locks.json` manifest records the monorepo identity, repository-relative
packet path, packet ID and version, and a fingerprint over the complete packet
directory. The smallest lock record is:

```toml
repository = "tannerpolley/ePC-SAFT"
packet_path = "data/packets/<packet-id>/<version>"
packet_id = "<packet-id>"
packet_version = "<version>"
packet_fingerprint = "<64-hex-fingerprint>"
materialization = "byte-for-byte"
```

The fingerprint is the raw 64-hex `packet_fingerprint` value emitted after
`packet` by `python3 tools/validate_packets.py` (the CLI also prints the
display-only `sha256:` prefix). A consumer may add its campaign or runtime
purpose and the exact materialization destination, but must not omit the
repository-relative path or fingerprint. `materialization` is
`byte-for-byte` for analysis or Regression packet materialization,
`subset-copy` for a self-contained admitted EOS subset copied from the packet,
or `subset-generated` when EOS generates an admitted self-contained subset.
EOS runtime copies must be self-contained. Analysis sets and Regression
materialize the packet before execution. Runtime fetching is forbidden.

The live lock manifest deliberately does not contain a commit field. A tracked
manifest cannot name the commit that first contains itself, and inventing or
retaining an earlier commit would falsely claim that new packet bytes exist in
that revision. The containing `tannerpolley/ePC-SAFT` commit binds the lock,
packet and in-repository materialization record atomically. An exported or
post-commit record may additionally record that containing commit after it
exists; it never replaces the packet fingerprint.
