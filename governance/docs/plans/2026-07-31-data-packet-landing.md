# Canonical Data packet landing and consumer locks

Status: Data packet landing is merged; consumer lock and replay work remains
with the owning package tasks.

## Accepted landing

`ePC-SAFT/ePC-SAFT-data` PR #7 was merged to `main` at commit
`cc197d686e03f79403e50db4200f227543912e82`. It contains source-traceable
packets for every reusable bundle currently listed by Management. The
Management lock manifest at
`manifests/data-packet-locks.json` binds each bundle to an immutable Data
commit, packet path, packet version, and packet fingerprint.

The Data validator accepts schema version 1 only for the two accepted
historical Khudaida packet versions. All newly landed packets use schema
version 2 and carry an explicit tolerance rationale. The landing preserved
the accepted Khudaida v2 and v3 bytes and fingerprints unchanged.

## Consumer handoff

EOS, Validation, and Regression must create consumer-owned lock records from
the Management/Data lock after their current work is reviewed. The records
must pin the Data commit and packet fingerprint, declare materialization, and
bind any self-contained EOS subset or campaign materialization. Existing EOS
and Validation copies remain until those locks and replay evidence are
accepted; no cleanup step may delete them merely because Data now has the
canonical packet.

Validation issues #9, #15, and #17 are the current evidence/replay handoff
subjects. They should reference the immutable Data lock and retain generated
campaign results in Validation. Management does not become a runtime or
campaign owner.

## Verification

- `python3 tools/validate_packets.py` is clean on the merged Data commit.
- The Data packet/unit test suite passes (22 tests).
- Management tests pass (21 tests using the repository's `.venv`).
- Management's full preflight remains sensitive to unrelated local wheel and
  environment drift; those findings stay with the existing artifact cleanup
  work and do not alter the Data lock.
