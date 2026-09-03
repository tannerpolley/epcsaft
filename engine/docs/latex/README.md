# Bibliography

This directory holds the project bibliography only. The equation registry lives
in `docs/equations_registry.yaml` (hand-edited authority) and its rendered view
`docs/equations.md`; there is no LaTeX equation source any more. Check the pair
from `engine/`:

```bash
python3 scripts/docs/sync_equation_registry.py --check
```

Zotero owns bibliographic metadata. `references.bib` is a generated snapshot of the
keys in `cited_keys.txt`, taken from `~/Zotero/exports/references.bib`. **Never
hand-edit it**; add or re-key the record in Zotero, then run
`engine/scripts/docs/sync_bibliography.py` (`--check` fails when the snapshot is stale).
The script refuses to emit local attachment paths and names any cited key the export
does not carry. The equation registry no longer reads this file: its `source` field is a
plain citation string and its check is Markdown-only (owner ruling, 2026-09-03). This
snapshot serves manuscripts.
