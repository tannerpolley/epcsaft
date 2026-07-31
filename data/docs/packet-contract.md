# Scientific packet contract

Status: canonical

Version: 2

## Purpose

A packet is the smallest reusable, reviewed scientific input with enough
metadata to interpret, audit, hash, and reuse it without relying on a package
source checkout. Packet storage does not admit runtime capability or validation
evidence.

## Layout

Packets live at `packets/<packet-id>/<version>/packet.toml`. Optional content
appears only when needed under `parameters/`, `observations/`, `sources/`, or
`provenance/`. Empty placeholder directories and unlisted files are invalid.

Full papers, broad extraction corpora, rejected experiments, and exploratory
notebooks remain in Lab. Generated campaign results, plots, and receipts remain
in Validation.

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

## Consumer lock

A consumer records the Data repository, immutable commit, packet path, packet
ID and version, and a fingerprint over the complete packet directory. The
smallest lock record is:

```toml
data_repository = "ePC-SAFT/ePC-SAFT-data"
data_commit = "<40-hex-commit>"
packet_path = "packets/<packet-id>/<version>"
packet_id = "<packet-id>"
packet_version = "<version>"
packet_fingerprint = "<64-hex-fingerprint>"
materialization = "byte-for-byte"
```

The fingerprint is the raw 64-hex `packet_fingerprint` value emitted after
`packet` by `python3 tools/validate_packets.py` (the CLI also prints the
display-only `sha256:` prefix) for the pinned commit. A consumer may add
its campaign or runtime purpose and the exact materialization destination, but
must not omit the immutable commit or fingerprint. `materialization` is
`byte-for-byte` for Validation or Regression packet materialization,
`subset-copy` for a self-contained admitted EOS subset copied from the packet,
or `subset-generated` when EOS generates an admitted self-contained subset.
EOS runtime copies must be self-contained. Validation and Regression
materialize the packet before execution. Runtime fetching is forbidden.
