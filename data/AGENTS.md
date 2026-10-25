# Data Directory Agent Contract

Read `docs/packet-contract.md` and the sections of `../docs/scientific/CONTEXT.md`
the task needs.

This directory owns canonical reusable scientific packets: parameters,
source-derived observations, provenance, units, domains, hashes, and use bases.
It owns no runtime algorithm, generated analysis result, package artifact,
full paper, or rejected experiment. Runtime code never imports or fetches Data.

Packets live at `packets/<packet-id>/<version>/packet.toml`. The fail-closed
`tools/validate_packets.py` command is a readiness check, not a scientific
test. Data owns no pytest or unittest suite. An Engine-accountable vertical
slice may revise a packet and its exact Engine/analysis consumers together,
but every materialization remains hash-bound and self-contained.

Preserve historical packet provenance and source repository identities.
