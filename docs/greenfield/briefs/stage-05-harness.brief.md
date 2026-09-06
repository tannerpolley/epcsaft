# Stage 05 harness — every check retains its numbers

Branch `stage-05-harness`, based on `greenfield@a518039`, which carries the Stage 5 repair merge `271fb39`; engine sources are pinned there.
Machine-readable claims are in `briefs/stage-05-harness.brief.json`; that file governs where
the two disagree. This is the owner's step 3: a check that passes and retains nothing is
blind, and the fourth audit found every Stage 5 pass blind.

## What this stage is

Every named check in the three native binaries (`specification_checks`, `equilibrium_checks`,
`boundary_checks`) writes what it solved to one results file per binary: its inputs, its
solved state, its convergence record, its deviations, and its status. Pass and fail are
unchanged. The results file is what acceptance reads and what the measured-data campaign
consumes, so its location and schema are stated once, here, and nowhere else.

The ledger is **26 registrations for 24 unique physical claims**: nine `C` registrations,
ten in `equilibrium_checks`, and seven in `boundary_checks`. The repeated physical claims
are intentional: `S10` is registered in both binaries and `S15(c)` is the boundary
subcase of `S15`. They remain separate `(binary, name)` entries in the results; the claims
checker must not collapse them or only check that a name appears.

**Location.** `${CHECK_RESULTS_DIR}/<binary>.json`, with `CHECK_RESULTS_DIR` set by
`meson.build` per test to `<build>/check-results`; a binary run outside meson with the
variable unset writes `./check-results/<binary>.json` and says so on stderr. The file is
rewritten from scratch on every run; nothing appends.

**Schema.** `engine/tests/native/check_results.schema.json` (JSON Schema, draft 2020-12),
validated by `tools/summarize-check-results.py --validate`. One file is one object:

    { "binary": "...", "commit": "<git sha or unknown>", "wheel_sha256": "<or unknown>",
      "started": "<ISO-8601 UTC>", "checks": [ Check, ... ] }
    Check = { "name": "check_s9", "status": "pass" | "fail" | "skipped",
              "skipped_reason": "<name of the failed certificate or missing declaration>" | null,
              "cells": [ Cell, ... ], "failures": [ "<file:line: condition>", ... ] }
    Cell  = { "label": "<the check's own name for the cell>",
              "inputs":  { "fixture": "...", "phases": [ {"name","kind","support","amount"} ],
                           "T": {"fixed"|"free", value}, "P": {...}, "feed": [...],
                           "reactions": [...], "active_parameters": [...] },
              "solved":  { "T", "P", "amounts": [...], "mole_fractions": [[...]], "molar_densities": [...],
                           "packing_fractions": [...], "multipliers": [...], "residuals": [...] } | null,
              "record":  { "residual", "criterion": {"name","value"}, "floor", "amplification",
                           "stop", "iterations", "evaluations" },
              "failure": { "name", "phase", "quantity", "value" } | null,
              "deviations": [ {"quantity": "...", "value", "reference", "deviation", "bar", "within": true|false} ],
              "events": [ {"span","name","iteration","values":[...]} ]   // the trace, always, not only on failure
            }

Every number is `double` written with 17 significant digits; nothing is rounded. A check
with no solve (the compile tier) writes `solved: null` and its `inputs` and `deviations`.
`deviations` is where a check puts what it compared: the centred difference and the analytic
column, the two pressures of S3, the recovered temperature of S5, each with its `bar` (the
tolerance it was held to) and `within`. The campaign's oracle cells use the same field with
`reference` the measured value and `bar` the stated bar.

**The trace is retained on every cell**, pass or fail, so an accepted run can be read later
for what it solved; the printer stays for failures.

**One more event.** When a free intensive has no guess (`FreeValue` absent), the start rule
seeds it, and the span records `intensive_seeded` with the lever index and the value used;
today `seed_or` fills 300 K or 1e5 Pa silently, which is the sentinel the repair removed on
the surface surviving in the trace.

## The order it is written in

1. `check_results.schema.json`, then the writer in `equilibrium_support.hpp`: `record_cell`
   called by a check with a label, the problem, the result and its deviations; `run_check`
   writes the check's entry with its status and failures; `main` writes the file. Hand-written
   JSON, no dependency: `std::ofstream`, `%.17g`, one function that escapes strings. The
   serializer recomputes `within` from the signed deviation and bar. Keep these mutation
   anchors as written in the serializer: `const double emitted_deviation =
   entry.deviation;` and `const bool emitted_within = std::abs(emitted_deviation) <=
   entry.bar;`, then emit both locals. At least one compile-tier structural row is zero on
   the green baseline with `bar: 0`; this is an equality already asserted by the check, not
   a scientific expected value.
2. The three check files: all 26 registrations write one status-bearing result entry. Every
   check has at least one cell; a check that solves records one cell for **each actual**
   `solve_equilibrium` call (including a failed solve's inputs, failure, convergence record,
   deviations, and trace), while a compile-only check records its inputs and structural
   deviations. Use fully qualified `equilibrium_support::record_cell` calls. In `check_s3`,
   the bubble cell is the first recording call and the dew cell is the second. A static
   lower bound from the literal solve call sites plus the runtime cell count is the harness
   coverage check; loop iterations may add cells. No assertion changes.
3. `meson.build`: `CHECK_RESULTS_DIR` per native test.
4. `phase_boundary.cpp`: the `intensive_seeded` event.
5. `tools/summarize-check-results.py`: `--validate` against the schema, a table of every
   cell's record and deviations for the reading, and `--assertions` that rejects a finite
   `within: true` row with `abs(deviation) > bar` (and a `within: false` row from a passing
   check). The claims checker reads all three results files, matches every `(binary, name)`
   registration, requires nonempty cells and statuses, and verifies that each check's cell
   count covers its literal solve-call lower bound. It recognizes `S15(c)` as the S15
   boundary registration without merging its result with S15.

## Installed-wheel readiness

`tools/ensure-engine-wheel.sh` remains the sole build/install caller and already invokes
`tools/verify-installed-wheel.py` after installing the non-editable wheel. This stage makes
that verifier describe the surface that exists now: public modules `epcsaft` and
`epcsaft.equilibrium`, and the bound native entry point `epcsaft._core.solve_equilibrium`.
The verifier continues to require exactly one extension, no transient wheel runpath, one
`epcsaft` distribution, a non-editable `direct_url.json`, and an installed origin equal to
the wheel path. Do not change `ensure-engine-wheel.sh`, add a regression module, add a fake
entry point, loosen provenance checks, or introduce another package route. The later
regression owner must extend and repin this verifier when that surface is real; this
stage's readiness gate is only `bash tools/ensure-engine-wheel.sh`.

## Traps

- Retaining a number is not asserting it. No tolerance, no assertion and no check moves; a
  check that records a deviation records the bar it already held it to.
- The results file is written by the binary from the values it computed, never parsed
  back from its own stderr (`HANDOFF.md`: compute a measured quantity, never parse the log
  that printed it).
- A `skipped` status carries the name of what was missing: the failed certificate, the
  missing declaration. Silent skips are what the fourth audit found.
- Read the tree through Serena before writing: `get_symbols_overview` on
  `tests/native/equilibrium_support.hpp` and the three check files; `find_symbol` for
  `run_check`, `print_trace`, `seed_temperature`, `seed_pressure`; `find_referencing_symbols`
  for `run_check`. Record what it found.

## Canon edits, the coordinator's

- `CANON.md` "What is tested": one paragraph, "Every check retains its numbers", pointing at
  the schema file and the results location above.
- `EQUILIBRIUM.md` §3.4 and S15's row: the failure of a declared topology is
  `declared_phase_absent`, with the mechanism the driver met as its detail, `phase_vanished`
  and the deficit or `phases_coincide` and the distance; on the subcooled side with
  unrestricted supports coincidence is the correct outcome, since `G` is extensive and the
  trivial manifold and the vanished-phase boundary minimum share one `G`; the deficit is
  asserted where a support restriction forces the vanishing (the electrolyte cell) and
  reported elsewhere.

## Checks, named

| id | asserts | mutation (must go red) |
|---|---|---|
| installed-wheel-verifies | `ensure-engine-wheel.sh` builds/installs the current non-editable wheel and its verifier accepts the current public surface and provenance | rename the required `solve_equilibrium` entry point |
| results-valid | every binary's results file validates against the schema | the writer drops the `record` field |
| every-check-present | all 26 `(binary, registration)` entries have a status and nonempty cells; every actual solve has a retained cell with its inputs, output/failure, record and trace; compile-only entries retain inputs and deviations | silence the first `equilibrium_support::record_cell` call in `boundary_checks.cpp`; S3's second solve remains, but its two-call lower bound and cell ledger go red |
| deviations-are-the-assertions | for every finite cell with `within: true`, `|deviation| ≤ bar`; for every `within: false`, the check's status is `fail`, and the green baseline contains a zero-bar structural equality | change the finite emitted deviation by `+1.0` and force the emitted `within` to `true`; JSON remains valid but the zero-bar row is a real contradiction |
| intensive-seeded | a boundary cell with `P` free and no guess carries one `intensive_seeded` event with lever index 1 and the value the start used | the event is not written |

## Stopping

`budget_ceiling`, `gate_impossible`, `source_contradicted`, `environment_broken`,
`mutation_invalid`. Ceilings only, measured as first written.
