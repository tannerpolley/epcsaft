# Read-only project notebook

The project notebook answers three questions from existing owners:

1. **What exists?** `ARCHITECTURE.yaml` owns components, Engine layers, path
   ownership, and artifact roles. `tools/project-map.py` resolves only files
   that exist in the checkout and reports unmapped or ambiguous paths.
2. **What remains?** `docs/scientific/implementation-plan.md` owns milestone
   status, prerequisites, requirements, and qualification evidence in
   Component → Layer → Milestone order.
3. **What happened?** Ordinary pytest XML, Meson logs, retained analysis
   tables, and the pre-cutover numerical manifest own execution outcomes,
   dates, inputs, values, failures, and build provenance.

These owners do not depend on the notebook. Tests and checks never import
notebook code, collect metadata for it, invoke Quarto, or change their verdicts
for presentation. Python and native numerical checks may emit
`EPCSAFT_COMPARE_JSON` lines to their ordinary captured output; those lines are
useful on their own and retain model/reference operands, units, inputs,
criteria, and limitations.

## One-way flow

```text
test / check / analysis
        ↓
ordinary XML, Meson log, JSON, CSV, or retained legacy record
        ↓
tools/render-notebook.py (read-only projection)
        ↓
Quarto Markdown and HTML
```

`tools/render-notebook.py` uses only the Python standard library and the two
existing owner parsers. It validates JSON, confines legacy paths to the legacy
result root, verifies retained comparison-table SHA-256 values, and reports
missing, corrupt, partial, or unavailable inputs explicitly. It neither runs
nor imports a test suite.

The root Quarto project owns navigation, page rendering, static serving, and
preview reload. Its declared resources include the project map, implementation
plan, ordinary test reports, and retained legacy manifest, so Quarto 1.10.18
re-renders when those files change outside the page-source directory. There is
no custom HTTP endpoint, browser application, poller, reload server, worker,
scheduler, cache selector, retry loop, or recovery state.

## Freshness and qualification

The notebook keeps these independent:

| Signal | Meaning |
|---|---|
| Execution | Passed, failed, skipped, errored, interrupted, or unavailable in the recorded run. |
| Produced | Timestamp written by the result owner, or the result-file modification time when that is all the ordinary format provides. |
| Source / inputs | Results newer than the relevant Engine source/input modification boundary are **Unverified**: modification time does not prove matching inputs. Older results are **Outdated by modification time**. Recorded fingerprints/provenance remain directly inspectable but are not continuously recertified. |
| Scientific qualification | Comes only from the implementation-plan milestone row and its reviewed evidence. |

Historical values remain rendered and downloadable when outdated. An old pass
never becomes a current-code pass. A failed or interrupted run remains failed
or incomplete. Missing and corrupt records are warnings, not empty green rows.
The retained F1 values and established F5 disagreement stay visible with their
execution dates, physical inputs, wheel provenance, and applicability limits.

Rendering does not build or install the Engine, execute tests, collect pytest,
launch native executables, rerun calculations, or continuously certify the
compiler and dependency environment. Build provenance remains in
`build/environment-wheel/preparation.json`; scientific physical inputs remain
in the numerical rows and retained analysis owners.

## Operation

Run the viewer directly:

```bash
quarto preview . --host 127.0.0.1 --port 8766 --no-browser --render all --no-execute
```

The persistent user service runs exactly that command and retains its existing
enabled, restart-on-failure startup behavior. `quarto render --no-execute`
produces the static site under `build/companion/`.

Normal verification remains independent:

```bash
bash tools/check-affected.sh HEAD
bash tools/check-all.sh
```

Those commands produce ordinary test output and never start or render the
notebook. Removing `_quarto.yml`, the three root `.qmd` files,
`tools/render-notebook.py`, and `build/companion/` leaves every test, check,
result table, and scientific owner usable.
