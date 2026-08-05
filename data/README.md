# ePC-SAFT Data

Canonical reusable scientific packets for the ePC-SAFT ecosystem.

This repository stores reviewed parameter bundles, source-derived datasets,
normalized observations, provenance, units, domains, hashes, and license or use
bases. It does not store runtime algorithms, generated Validation results,
package artifacts, full papers, or exploratory and rejected material.

Each packet is self-contained at:

```text
packets/<packet-id>/<version>/
  packet.toml
  parameters/     # only when present
  observations/   # only when present
  sources/        # only when distributable and present
  provenance/     # only when present
```

Read [the packet contract](docs/packet-contract.md) before adding data.
Read [the data-organization guide](docs/data-organization.md) before deciding
whether a source record belongs in Data or remains campaign evidence in
Validation.

Validate packet schema, provenance, units, domains, and hashes with:

```bash
python3 tools/validate_packets.py
```

Consumers must pin an immutable Git commit and packet fingerprint, materialize
the packet before execution, and remain self-contained. Installed packages must
not fetch or import this repository at runtime.
