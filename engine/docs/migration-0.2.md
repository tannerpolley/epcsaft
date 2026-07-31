# Python frontend migration: explicit parameters

The installed EOS parameter catalog and public typed-record constructor were
removed in `0.2.0.dev0`. This is a hard cutover with no compatibility alias.

| Removed call | Replacement |
| --- | --- |
| `Parameters.from_catalog(...)` | Materialize the desired Data packet, then call `Parameters.from_bundle(path)`; or pass compact values to `Parameters.from_dictionary(mapping)`. |
| `Parameters.from_records(...)` | `Parameters.from_dictionary(mapping)` |
| implicit installed bundle lookup | explicit caller-owned path |

The dictionary route requires schema version 1, ordered components, parallel
vectors, square pair matrices, options, and explicit validity. It accepts the
documented abbreviation or descriptive spelling for each field but rejects
both together. See [the Python API guide](python-api.md).

Reusable source-traceable packets now live in ePC-SAFT-data. EOS does not import
that package or locate a sibling checkout, so deployment remains a single
installed EOS artifact plus whatever explicit parameter material the caller
chooses to provide.

`Parameters -> Mixture -> State` and the property API are otherwise unchanged.
Parameter fingerprints now identify normalized scientific inputs rather than
bundle provenance. SDK admission likewise uses resolved topology, fields,
order, and validity rather than catalog allowlists.
