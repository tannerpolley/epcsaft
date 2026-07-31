# Scientific packet contract

Status: canonical

Version: 1

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

- `schema_version = 1`;
- `packet_id`, matching the packet directory;
- `packet_version`, matching the version directory;
- a short `title` and `description`;
- `license_or_use_basis`, including any redistribution restriction;
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
ID and version, and a fingerprint over the complete packet directory. EOS
runtime copies must be self-contained. Validation and Regression materialize
the packet before execution. Runtime fetching is forbidden.
