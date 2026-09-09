# Companion cleanup and project restructuring record

Date: 2026-09-09. Baseline commit `8b6835c` (the uncommitted test-system and
companion work as found); cleanup commits follow it on `greenfield`. This record
replaces the migration brief formerly at `docs/scientific/project-map-design.md`;
its durable vocabulary now lives in [CONTEXT.md](../CONTEXT.md), its signal table
in the [implementation plan](../implementation-plan.md), and its reduced viewer
boundary in the [test-system design](../test-system-design.md). No scientific
milestone is qualified by this work.

## Owner rulings applied

- Automatic execution on source edits is retired. The viewer is read-only; the
  systemd unit starts only `tools/serve-roadmap.py`.
- The `/coverage`, `/render`, `/source` and `/artifact` endpoints had no caller
  and are removed.
- `engine/native/physics/` → `engine/native/contributions/`,
  `engine/native/model/` → `engine/native/properties/`, with
  `contributions_checks`, `properties_checks`, `properties_floor` and
  `properties_support`; `engine/ARCHITECTURE.yaml` → `ARCHITECTURE.yaml`.

## Declaration-ID changes

Native declaration IDs carried the old layer names. `native.physics.*` became
`native.contributions.*` and `native.model.*` became `native.properties.*`
(`F1` alias now resolves to `native.properties.f1`). Both suites were rerun
through the shared recorder under the new IDs: 8 Contributions and 9 Properties
cases pass. Records under the old IDs are retained as history; they are orphaned
from the current inventory and cannot appear as current entries.

## Preservation observations

- Full installed gate (`tools/check-all.sh`): passed, 86.4 s wall. Stage times:
  recorded-evidence 53.1 s, mypy 7.6 s, analysis-layout 2.9 s, architecture
  2.2 s, packets 2.2 s, ruff 2.1 s.
- Integrity harness (`tools/check-test-system.py`): all 15 retained checks pass.
- F1 (`native.properties.f1`): Passed, 11 complete comparison rows including
  the 4-row temperature ladder, 3 distinct criteria.
- F5 (`native.equilibrium.f5`): still Failed, 1,159 complete rows, 198 false
  verdicts, 4 distinct criteria. The recorded failure is unchanged.
- Architecture check: 30 native sources, 0 rules failing; the wheel rebuilt
  after the path migration and installed non-editable.
- Companion: warm snapshot assembly 0.72 s (input generation 0.27 s twice,
  run reading 0.19 s, project map 0.02 s); payload 3.2 MB, of which the map is
  138 KB and recorded runs 3.0 MB. Map, Plan and Tests render; `#owner=` carries
  a selection across the three pages; `/coverage` returns 404.
- Project map: 15 owners, every project file resolved, no unmapped or
  ambiguous paths (`tools/project-map.py`, now a check-all gate row).

## Line accounting (physical lines, baseline → after)

| Group | Files | Before | After |
|---|---|---:|---:|
| Runtime/tooling (audited) | `test_records.py`, `run-numerical-evidence.py`, `serve-roadmap.py`, `companion-data.py`, `numerical-inventory.py`, `test-worker.py`, `check-all.sh`, `companion.js` | 5,076 | 4,616 |
| Integrity harness (audited) | `check-test-system.py` | 1,858 | 1,618 |
| Unchanged (audited) | `pytest_recording.py`, `record-command.py` | 744 | 744 |
| Audited eleven files | | 7,678 | 6,978 |
| Also touched | `check-roadmap.py`, `check-affected.sh` | 263 | 275 |
| New | `project-map.py` | 0 | 133 |

Native test sources changed by rename only (107 lines replaced in 9 files);
no tolerance, formula or criterion changed. The Map view added about 130 lines
of browser code and a 60-line declaration in `ARCHITECTURE.yaml`.

## What was removed and why

- `tools/test-worker.py` and its supervision in the server, worker fields in the
  projection, the UI badge, and the worker restart/coalescing harness check.
- Refresh-time environment verification: display freshness now means the input
  generation matches and the recorded wheel/binaries are the files on disk.
  Build/run-time verification keeps the full runtime, compiler, dependency and
  file check; the verification cache and runtime-context probe that only served
  repeated refreshes are gone.
- Dead `_inputs`/`_number` text-input helpers in the runner.
- The background PID/wait scheduler in `check-all.sh` (foreground loop) and the
  duplicated roadmap/budget/registry runs when `check-affected.sh` hands over to
  the full gate.
- Duplicate milestone-table parsing: `check-roadmap.rows` is the single parser;
  validation, summary and the companion derive from it, and rows now carry
  prerequisites and evidence links.

## Not done

- Finding 9 (finalization and interruption branches across the runner, pytest
  integration and command recorder) was investigated only to the extent of
  confirming the harness still exercises each entry point; no branch was deleted
  because no redundant path was established.
- Application repositories under `downstream/` are summarized by presence, not
  inventoried.
