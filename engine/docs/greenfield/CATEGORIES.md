# What is engine, what is core, what is peripheral

Status: inventory and proposal, 2026-09-02. Prepared as input for the agent that will own
the core, the input surface and everything that is not the engine. **The classification
below is measured; the taxonomy is a proposal and that agent decides it.**

## The measurement

The equation-registry annotation turns out to separate the layers almost perfectly. Every
file in `native/physics/` either carries many `// EqID:` annotations or none at all:

| file | lines | EqIDs |
|---|---|---|
| `helmholtz.hpp` | 280 | **23** |
| `electrolyte.hpp` | 346 | **16** |
| `association.hpp` | 316 | **3** |
| `born.hpp` | 87 | **3** |
| `numerical_policy.hpp` | 248 | 0 |
| `parameters.hpp` | 86 | 0 |
| `physical_constants.hpp` | 28 | 0 |
| `scalar.hpp` | 41 | 0 |

Four files carry all 45 equations. Four carry none, and they are exactly: the tolerance
policy, the parameter container, the measured constants, and the transcendental dispatch.
**None of those four is physics.** `target-layout.md:86` already suspected one of them —
*"`numerical_policy.hpp` in `physics/`: every layer reads it, and `physics/` is the lowest,
so it is reachable; but it is not physics"* — and rejected a `native/policy/` directory
because it would hold one file, which the no-single-file-directory rule forbids. It would
hold four, and with the diagnostics work more.

`native/model/` splits the same way, though less sharply: `helmholtz_model.cpp` (7) and
`infinite_dilution_reference.cpp` (6) carry equations; `helmholtz_tape.hpp` (275 lines, 0)
is automatic-differentiation machinery, `density.cpp` (218, 0) is a solver, and
`parameter_builder.{hpp,cpp}` (425, 0) is the input contract — it resolves a declared
document into arrays and is where the `ionic_dispersion` selector should have lived.

## Proposed taxonomy

**Engine — the science.** Files that carry registered equations, and the solvers that
evaluate them. `helmholtz`, `electrolyte`, `association`, `born`; the model's Helmholtz
assembly and reference; equilibrium and regression when they exist. A file here answers a
thermodynamic question.

**Core — what every layer depends on and no layer owns.** Numerical policy and its
categories; the diagnostics and convergence record; measured constants; the transcendental
dispatch; the automatic-differentiation tape. A file here answers *how the engine computes
and reports*, never *what is true of the fluid*. It sits below physics in the include
direction, and it may not know a phase, a route or a problem kind.

**Surface — how a caller states a problem and reads an answer.** `native/bindings/`,
`src/epcsaft/`, and — this is the reclassification that matters — `parameter_builder`. The
input surface is where a model choice is declared or lost: `stoichiometry=(1, 1)` at
`eos.py:175` is a hidden physical default, and the same-sign dispersion selector vanished
because nothing in the surface carried it. Treating the builder as part of `model/` is what
made it look like a mapping detail rather than a contract.

**Peripheral — real work, outside the package.** The architecture checker (227 lines), the
packet validator (572) and promoter (455), the registry and bibliography sync (233), the
data packets, the analyses, the docs. These are not lesser; two of this week's findings came
from the validator and the registry. They are simply not what a caller imports.

## What this proposal implies, stated so it can be argued against

1. **Four files leave `physics/`.** It stops being "the lowest layer" and becomes what its
   name says: the equations. The checker's include rule gains a layer beneath it.
2. **`parameter_builder` leaves `model/`.** With it goes the argument that the input surface
   is a mapping detail — it is the contract that says which model is being asked for.
3. **`helmholtz_tape.hpp` is arguable and I have not decided it.** It is AD machinery, which
   is core by the test above; but it is shaped entirely around one Helmholtz model, which is
   engine. Named as an open question rather than filled in.
4. **`physical_constants.hpp` is arguable.** Measured constants are facts about the world,
   which sounds like engine; but they are reference data rather than equations, and every
   layer reads them. Also named rather than decided.

## The question this exists to serve

"No hidden defaults" and "every solver reports its convergence" are the same discipline
pointed at inputs and at outputs. If the core owns the convergence contract, the argument
that it should also own the declared-choice contract is that both are promises the engine
makes to a caller, enforced in one place, and both were broken this week in the same way —
by a decision that lived in a layer nobody owned.
