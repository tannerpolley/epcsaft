# ePC-SAFT Data Context

`GOVERNANCE.md` revision 6 in `ePC-SAFT/ePC-SAFT-management` is the canonical
ecosystem doctrine.

Repository role: `scientific_asset`

GitHub home: `ePC-SAFT/ePC-SAFT-data`

Data owns canonical reusable scientific packet bytes and their source,
transformation, units, domain, hash, and license or use-basis metadata.

Storage authority is independent of runtime and evidence authority:

- `storage_owner = data` for reusable reviewed packets;
- `runtime_owner = eos` only for a separately admitted self-contained EOS copy;
- `evidence_owner = validation` for generated campaign evidence.

Lab remains the home for exploratory, rejected, broad, or incomplete material.
Management remains the home for policy, registries, operation receipts, wheels,
and SDK artifact objects.

Consumers pin a Data commit, packet path, and packet fingerprint. No installed
package has a runtime dependency on this repository.
