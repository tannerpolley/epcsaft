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
   JSON, no dependency: `std::ofstream`, `%.17g`, one function that escapes strings.
2. The three check files: every check records its cells. Twenty-four checks; a check that
   solves records each solve, a check that compiles records its inputs. No assertion changes.
3. `meson.build`: `CHECK_RESULTS_DIR` per native test.
4. `phase_boundary.cpp`: the `intensive_seeded` event.
5. `tools/summarize-check-results.py`: `--validate` against the schema, and a table of every
   cell's record and deviations for the reading; the claims checker gains one clause, every
   named check present in its binary's results file.

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
| results-valid | every binary's results file validates against the schema | the writer drops the `record` field |
| every-check-present | every named check of `EQUILIBRIUM.md` §4 appears in its binary's results file with a status | one check stops recording |
| deviations-are-the-assertions | for every cell with `within: true`, `|deviation| ≤ bar`; for every `within: false`, the check's status is `fail` | a check writes `within` as always true |
| intensive-seeded | a boundary cell with `P` free and no guess carries one `intensive_seeded` event with lever index 1 and the value the start used | the event is not written |

## Stopping

`budget_ceiling`, `gate_impossible`, `source_contradicted`, `environment_broken`,
`mutation_invalid`. Ceilings only, measured as first written.
