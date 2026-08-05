# Data Directory Agent Contract

Repository Profile: scientific-computing

Read `../governance/GOVERNANCE.md`, the root contract, `CONTEXT.md`, and
`docs/packet-contract.md`.

This directory owns canonical reusable scientific packets: parameters,
source-derived observations, provenance, units, domains, hashes, and use bases.
It owns no runtime algorithm, generated Validation result, package artifact,
full paper, or rejected experiment. Runtime code never imports or fetches Data.

Packets live at `packets/<packet-id>/<version>/packet.toml`. The fail-closed
`tools/validate_packets.py` command is a readiness check, not a scientific
test. Data owns no pytest or unittest suite. An Engine-accountable vertical
slice may revise a packet and its exact Engine/Validation consumers together,
but every materialization remains hash-bound and self-contained.

Use the sole root Codex environment and root checks. Preserve historical packet
provenance and source repository identities. No CI/CD or compatibility paths.
