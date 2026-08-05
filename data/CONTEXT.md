# ePC-SAFT Data Context

`GOVERNANCE.md` revision 10 in `ePC-SAFT/ePC-SAFT-governance` is the canonical
ecosystem doctrine.

Repository role: `scientific_asset`

GitHub home: `ePC-SAFT/ePC-SAFT-data`

Data owns canonical reusable scientific packet bytes and their source,
transformation, units, domain, hash, and license or use-basis metadata.

For the ePC-SAFT ecosystem, this repository is the canonical home for
reusable, source-traceable ePC-SAFT parameter packets or bundles. A parameter
packet is a versioned, internally compatible model formulation; its values
must not be treated as universal species properties.

The reusable-data boundary has two admitted forms: source-traceable ePC-SAFT
parameter bundles and reusable property-observation packets. Property packets
retain organized source values, conditions, units, source locators, and source
uncertainty where available. Paper- or figure-specific matrices, direct
experimental records kept only for a named campaign, fitted alternatives, and
generated evidence remain in Validation even when they are stored as CSV.
See `docs/data-organization.md` for the current inventory and examples.
The Khudaida parameter-only successor is packet version 4; versions 2 and 3
remain immutable historical packets until downstream locks move.

Storage authority is independent of runtime and evidence authority:

- `storage_owner = data` for reusable reviewed packets;
- `runtime_owner = engine` only for a separately admitted self-contained copy;
- `evidence_owner = validation` for generated campaign evidence.

Lab remains the home for exploratory, rejected, broad, or incomplete material.
Governance remains the durable repository home for policy, registries,
operation receipts, wheels, and Engine artifact objects. Management is the
coordination project, not a Git repository.

Consumers pin a Data commit, packet path, and packet fingerprint. No installed
package has a runtime dependency on this repository.
