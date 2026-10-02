# Historical companion preservation observations

Date: 2026-09-09. Baseline commit `8b6835c`. These observations concern
the retired recorder and viewer, not current execution requirements. The
[test-system design](../test-system-design.md) owns current viewer behavior;
the [notebook-cutover record](notebook-cutover-2026-09-09.md) retains its checks.
No scientific milestone is qualified by this record.

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
