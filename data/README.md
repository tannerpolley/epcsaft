# ePC-SAFT Data

Canonical reusable scientific packets for the ePC-SAFT project.

Packets live at `packets/<packet-id>/<version>/packet.toml` and include their
source, provenance, units, domains, hashes, and use basis. Validate all packets
from the repository root:

```bash
.venv/bin/python data/tools/validate_packets.py
```

Data contains no runtime package or scientific test suite. Consumers
materialize exact packets and never import or fetch this directory at runtime.
