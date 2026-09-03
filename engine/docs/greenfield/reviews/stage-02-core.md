# Review of Stage 2 (the core) at `greenfield@31914d9`, against `CORE.md`

Planner's pass, 2026-09-02. Read the tree, not the worker's report. Verdict first: **the
three written headers are what §4, §5 and §6 specify**, and the population of the two
existing solvers is faithful. One behaviour change slipped in that no gate can see, one
lifetime bug in the bindings, two incompleteness items on the trace's root, and one dead pair
of enum members the new guard exposed. None blocks 04.repair; the first should be fixed
before it, since the repair rewrites the loop it lives in.

## 1. The headers (question 1)

| header | against | verdict |
|---|---|---|
| `core/declared.hpp` | §4.1 | exact: `declared()`, `is()` never throws, `require(name)` throws `MissingDeclaration` naming the field |
| `core/convergence.hpp` | §5.1–5.2 | exact record and `classify` order; `ResidualHistory` on a `std::deque` rather than a ring, same statistic, absent below the window; `describe` prints floor and amplification too, which is better than specified |
| `core/trace.hpp` | §6.1 | exact: one thread-local read in the `Scope` constructor, probes test a member, `steady_clock` read only when on, `optional` parent, `Capture` shadows and restores |

Two hardening nits, neither a deviation: `Scope` should delete copy and move (a copied
`Scope` closes its span twice; `Capture` already deletes them); `ConvergenceFailure` could be
`final` like everything else in the layer.

`numerical_policy.hpp`: 29 entries, all `inline constexpr`, all self-named (measured by
regex); `derived()` as specified; `DENSITY_ROOT_RESIDUAL` carries the accuracy sentence but
not the hypothesis sentence (*"5 is a declared hypothesis about how association noise reaches
the reduced residual, to be read against the measured floor / ASSOCIATION_RESIDUAL"*) — the
reader of the entry cannot tell why 5 without it. Two lines.

## 2. One behaviour change no gate sees — fix before 04.repair

`density.cpp:56-61`, `association_failure`, now accepts only `left_domain` and
`budget_exhausted`. The string it replaced matched *every* association non-convergence. The
association solver's budget exit classifies (`make_record(..., classify = true)`), so at the
budget with contraction ≥ 1 it throws `stalled`, and the density scan then **propagates the
exception instead of skipping the scan point** — the whole solve fails where it used to
continue. That is the failure mode the old string existed for, rare (needs a non-contracting
site solve at 50 steps) and invisible to `trace_changes_no_answer`, which compares captured
against uncaptured, not new against old.

Fix: the predicate should ask *who threw*, not *how it stopped*:

```cpp
return failure != nullptr && failure->record.criterion == &policy::ASSOCIATION_RESIDUAL;
```

One line, and it is exactly what `criterion` on the record is for. The repair rewrites this
file's anchored path; land this first so the repair's base is behaviour-identical to
`21b8695`.

Related, cosmetic: the anchored loop's own budget exit sets `stalled` by hand
(`density.cpp`, `if (!finished)`) rather than calling `classify(..., budget_hit = true)`, and
keeps no `ResidualHistory`, so its record's `contraction` is always absent. Spec §5.2 says a
solver never names a residual stop itself. Six lines; 04.repair touches that block anyway.

## 3. `Declared<T>` and the builder (question 2)

The guard sits where §4.2 put it: the type in `core/`, `.is()` at the fourteen read sites,
the two `require(...)` lines in the constructor, and the surface builder as the one place a
document's choice becomes a `Declared` value. The worker's deviation — a packet choice named
`none` leaves the selector undeclared — is **forced by the constructor rule** *declared iff a
charged component exists*, and it is the right reading: for a neutral set there is nothing to
declare, and "no electrostatic term" is a consequence, not a choice.

What it exposes: `ElectrolyteForm::None` and `PermittivityForm::None` are now **unreachable
from the surface** (measured: the four packets declaring `none` all have zero ions; no ion
packet declares it), and reachable only from a hand-built struct, which the constructor then
rejects. Dead members inside a selector are the `EOS-AUDIT.md` C.4 shape. Recommend deleting
both members: "undeclared" then carries the only meaning they had, the two reads that mention
them (`helmholtz_model.cpp:209` becomes `permittivity_form.declared()`; the builder's `none`
branches already do nothing) shrink, and rule 11's job gets simpler. If a future source ever
wants ions with no electrostatic term, that is a new declared choice with a paper behind it,
not a resurrected `None`. Not a Stage 2 defect; a follow-up of one commit.

## 4. The `state_at_pressure` span (question 3)

**A faithful reading, and the general rule it implies is only half applied.** Without it a
capture around `Mixture.state(P=...)` has two roots — `solve_density`, then the trailing
`state()` → `evaluate` — and §8.2's tree property fails. §6.1 says one span per instrumented
function call and §6.1's activation is "a lexical scope around one method call", so *every
public entry point opens a span* is the rule, and `state_at_pressure` is one instance of it.
Two instances are missing, measured:

- `HelmholtzModel::state` (`helmholtz_model.cpp:519`) has no span and calls `evaluate` and
  then `solve_site_fractions` directly (`:531`), so a capture around `Mixture.state(T, rho=...)`
  has **two roots**: `evaluate` and `association_newton`. The tree check runs only on
  `state_at_pressure`, which is why it passed.
- `electrolyte_activity` (`infinite_dilution_reference.cpp:165`) has no span; a capture
  around `Reference.activity(m)` has as many roots as fixed-pressure states it solves.

Fix: `core::Scope scope("state");` and `core::Scope scope("electrolyte_activity");` at the
top of each, and `trace_invariance` extended to `state` and `electrolyte_activity` — the same
property on two more entries, not a new check. Rule to write into `CORE.md` §6.1 so no later
stage misses it: *every function bound in `bindings/` opens a `Scope`.* That is a grep the
checker can do (a bound function name must appear beside a `core::Scope` line in its
definition), and is the sixth one-directional check of the week if it is not done.

## 5. Bindings

`diagnostics_bindings.cpp:97-107`: the `Capture` wrapper's `trace` property casts a pointer
into `wrapper.captured` with `return_value_policy::reference`, and `spans`/`events`
dereference `captured` without checking it. Consequences: a `Trace` object outliving its
`Capture` dangles (should be `reference_internal`), and `captured.spans` accessed inside the
`with` block, before `__exit__`, is undefined behaviour (should return an empty list or raise
by name). The Python test reads `spans` after the block, so neither path is exercised. Two
lines.

`state.density_convergence` and `ConvergenceError`/`DeclarationError` are bound as specified;
the three ratios are native properties; the exception order is right.

## 6. Checker

Rules 8–12 are as specified and their mutations were shown red. Two notes: rule 12 writes
`__import__("ast")` where `ast` is already imported (cosmetic); rule 9's `Convergence`
test is a substring and would accept a comment-free file that merely names the type — it is
the *record on the result* that matters, and the record-consistency check is what carries
that, so no change, but the rule is presence, not correctness, and `CORE.md` §7 should say so.

## 7. Reported numbers, read

From the floor suite on this tree: case G, nine-species cold `Liquid`, uncaptured 8,880 µs,
captured 9,440 µs — **6 % on-cost** for 221 spans, 1,443 events, 35.8 KB of values; case C
association converged in 6 iterations, contraction 6e-4, pivot ratio 3.0; case D's anchored
solve converged in one iteration with `contraction absent`, as the window rule says.

One correction to my own brief: §4 said the associating binary at 313.15 K, 1e5 Pa `Liquid`
would print `budget_exhausted` from the floor harness "on this tree". The harness has no such
density case (its density cases are D and E), so nothing prints it. The claim was about a
case the harness does not run; the record is still right where it is exercised.

## 8. The one-directional-check family

The coordinator's `types_consistent` finding is the fifth of a family: a rule or a record
with a reader on one side only. `CORE.md` §11 named the preflight as a missing fact; the
family deserves its own line there and in `CATEGORIES.md`'s canon section, because the
sweep of the rules that nobody has done yet is a sweep for exactly this shape: *for each
enforcer, name its consumer; for each record, name its enforcer.* Rule 10's "nothing reads
the trace" is the same idea inverted, and it is why the trace can be a debug mode at all.

## Disposition

| item | severity | when |
|---|---|---|
| `association_failure` propagates `stalled` | behaviour change, no gate sees it | before 04.repair, one line |
| `state` and `electrolyte_activity` have no span | tree property holds on one entry of three | with the above, or in 04.repair's delta |
| `Capture` wrapper lifetime and unchecked dereference | bindings bug, unexercised | same commit |
| anchored loop names `stalled` by hand, no history | spec §5.2, cosmetic | 04.repair touches the block |
| dead `None` members | exposed by the guard, not caused by the stage | follow-up commit |
| hypothesis sentence on `DENSITY_ROOT_RESIDUAL` | documentation | any time |
| `Scope` copy/move, `final`, `__import__` | nits | any time |
