# Read-only status report

The status report answers three questions from existing owners:

1. **What exists?** `ARCHITECTURE.yaml` owns components, Engine layers, path
   ownership, and file roles. `tools/project-map.py` resolves only files
   that exist in the checkout and reports unmapped or ambiguous paths.
2. **What remains?** Open GitHub issues ([ADR 0010](decisions/0010-github-issue-roadmap-authority.md));
   the report does not copy their mutable state.
3. **What happened?** Ordinary JUnit XML, key-solve comparison records and
   retained analysis results own execution outcomes, dates, inputs, values,
   failures, and build provenance.

These owners do not depend on the report. Tests and checks never import
report code or change their verdicts for presentation. Python and native
numerical checks may emit `EPCSAFT_COMPARE_JSON` lines to their ordinary
captured output; those lines are useful on their own and retain
model/reference operands, units, inputs, criteria, and limitations.

## One-way flow

```text
test / check / key-solve comparison / analysis
        ↓
JUnit XML, summary.json, analyses/README.md index and retained figures
        ↓
tools/report.py (read-only projection, standard library only)
        ↓
build/report.html (one static page, opened from disk)
```

`tools/check-all.sh` writes the page at its end, pass or fail, without changing
its exit status; `tools/report` rewrites it on demand. Writing it never builds
or imports the Engine, runs a test or solve, or renders a notebook. Each
`analyses/**/notebook.qmd` stays its analysis's record and renders on its own.

## Freshness

Each section states the commit and time of its data:

| Section | Source | Commit and time |
|---|---|---|
| Tests | `build/test-results/native.xml`, `pytest.xml` | The `.source` line written beside each XML (commit, `engine/` clean or dirty); file modification time |
| Key solves | main checkout `build/wheels/key-solves-*/summary.json`, written by `tools/refresh-main-wheel` | `source_commit`; the directory's timestamp |
| Analyses | `analyses/README.md` index | Last commit touching the leaf; uncommitted changes in the leaf |

Data from a commit other than HEAD, or from a dirty tree, is flagged. A flag
says only that the record may not describe the current code. An old pass never
becomes a current-code pass; a failed run stays failed. A missing or unreadable
record shows "no data", never an empty green row.
