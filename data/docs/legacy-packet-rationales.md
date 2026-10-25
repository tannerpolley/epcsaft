# Legacy packet tolerance rationales

Schema-version-1 packet directories are accepted historical records. Their
bytes and lock fingerprints are intentionally immutable; the schema-version-2
`tolerance_rationale` field is not backfilled into them.

The retained Khudaida versions use the source-derived uncertainty record as
their tolerance rationale:

| Packet | Accepted fingerprint | Rationale |
|---|---|---|
| `khudaida-2026-figure-2-electrolyte-lle/2` | `ed173d78fcaeca527c89a49af1f34c6880fdc642fc1c08c2d395f131991e53af` | `observations/table-3-standard-uncertainties.csv` retains the reported composition standard uncertainties. No additional acceptance tolerance is asserted; each consuming campaign declares its own comparison tolerance. |
| `khudaida-2026-figure-2-electrolyte-lle/3` | `73c48ac9fb97f06185bd78d27ceaa7ffe6308abe9feee34b15eadff31cc985eb` | Same source-derived uncertainty record and no additional packet-level acceptance tolerance; consuming campaigns declare their own comparison tolerance. |

New packets and new versions use schema version 2 and must carry the rationale
directly in `packet.toml`.
