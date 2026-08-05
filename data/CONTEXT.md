# Data Context

Governance doctrine revision 11 applies. `data/` is the canonical home for
reusable reviewed scientific packets inside `tannerpolley/ePC-SAFT`.

A packet is self-contained, versioned, hash-bound, source-traceable, and
licensed or use-basis declared. It is scientific input, not executable runtime
authority. Engine and Validation materialize exact packet versions and remain
self-contained; no installed package fetches or imports this directory.

Packet readiness is enforced by `tools/validate_packets.py`. Numerical model
agreement belongs in Engine or Validation.
