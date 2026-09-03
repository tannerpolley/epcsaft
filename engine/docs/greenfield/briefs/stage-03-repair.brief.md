# Stage 03.repair — the surface's three findings

Branch `stage-03-repair`, based on `greenfield@1dd5db9`. Machine-readable claims are in
`briefs/stage-03-repair.brief.json`; that file governs where the two disagree.

Stage 3 closed at `b14851f` with nine gates green and nine mutations red, and reopened on
three findings that no gate could have caught. Each is the same shape as the defect the EOS
repair fixed: **something asserted where nothing looks.** You are not redesigning the
surface. You are making three rules true, and putting a check behind two of them so they
cannot silently come untrue again.

## 1. Two selector members the surface cannot reach

`core/parameters.hpp` gives `PermittivityForm` and `ElectrolyteForm` a `None = 0` member.
The parameter builder never constructs either: a packet whose `electrolyte_choice` or
`permittivity_choice` is `"none"` falls through the `else if (... != "none")` guard and
leaves the field **undeclared**. Four packets say `"none"` and all four have zero ions.

So `None` is dead inside a live selector, and `declared()` already carries its only meaning.
Delete both members.

One consumer must change with them. `helmholtz_model.cpp` reads:

```cpp
if (parameters_.permittivity_form.is(PermittivityForm::None) == false) {
```

`is()` on an undeclared value returns false, so that condition is **true for every parameter
set that exists** — undeclared or not. It is a guard that never guards. It becomes:

```cpp
if (parameters_.permittivity_form.declared()) {
```

which is what it was always trying to say. Note this changes behaviour for neutral sets: they
stop being required to carry a dielectric solvent. That is correct and is the point.

## 2. The Scope rule is half applied

`CORE.md` §6.1 intends every function reachable from `bindings/` to open a `core::Scope`, so
a capture has exactly one root. Two do not: `HelmholtzModel::state` and
`electrolyte_activity`. A capture around `Mixture.state(T, rho=…)` therefore has two roots.

It survived nineteen gates and nineteen mutations because `trace_case` in
`diagnostics_checks.cpp` drives **only** `state_at_pressure`. The test and the defect were
written against the same blind spot.

Add the two `Scope` lines, following `core::Scope scope("state_at_pressure");` exactly. Add
two `trace_invariance` cases that exercise `state` and `electrolyte_activity`. Then write the
rule into `CORE.md` §6.1 **and put a checker rule behind it** —
`rule_bound_functions_open_a_scope` in `tools/check-engine-architecture.py`: for each symbol
bound in `native/bindings/`, its definition opens a `core::Scope`. A rule in prose with
nothing reading it is the defect, not the fix.

## 3. Three residual stops named by hand

`CORE.md` §5.2 is unambiguous:

> `left_domain`, `singular`, `infeasible` are set by the solver that detects them. A solver
> never computes a residual stop itself; it calls `classify`.

The Stage 3 review recorded one violation, in the anchored loop. Writing the rule as a grep
finds **three**:

| site | stop | verdict |
|---|---|---|
| `density.cpp` far-side-clear exit | `converged` | violation |
| `density.cpp` budget exit | `stalled` | violation, the one recorded |
| `density.cpp` domain exit | `left_domain` | **legal** — a solver stop, leave it |
| `association.hpp` after its own residual test | `converged` | violation, found by the rule |

All three call `classify` with the arguments they already hold. Where the solver had it right,
`classify` returns the same stop; where the residual sits under the floor margin, `classify`
returns `converged_to_floor` and the record becomes more honest than it was. Do **not** change
what `classify` returns, and do not touch a tolerance — that is a `scope_expansion` stop.

The positional `core::Stop::stalled` inside the two `core::Convergence` aggregate initialisers
is a placeholder that `classify` overwrites on the next statement. It is not an assignment and
not a violation. `rule_no_hand_named_residual_stop` must match `.stop = core::Stop::<residual>`
and must not match those initialisers.

## Budgets

Every write is small; the two checker rules are the largest at ~40 lines. Ceilings are in the
JSON. Nothing here should approach them — if a fix is growing, you have found a design
question, and a design question is a stop, not a judgement call.

## Stopping

Stop and report rather than adapting if: a `must_contain` marker is absent; a write would pass
its ceiling; a gate that was green at base goes red for a reason other than the write it names;
or a fix appears to require changing `classify`'s contract or a tolerance value. Report the
marker, the file, and what you observed. Do not repair your environment, and do not accept
your own work — the coordinator re-runs every gate and every mutation independently.
