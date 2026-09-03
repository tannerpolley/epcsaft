# Canon — the tree, the rules, and what may be tested

Stage 1's design document, and the only stage whose product is the canon itself.

The layer tree, one static archive and one extension, and the checker's rules live in
`tools/check-engine-architecture.py` — 227 lines, each rule mechanical, each shown to fail
on a deliberate violation. `ORCHESTRATION.md` owns how a stage is run and accepted.

What may be tested follows, unchanged: it is a restriction, not a coverage target.

---


The tree this replaces was hard to change because of its tests, not despite them.
`test_internal_certificates.py` asserted on certificate struct fields; the gate
carried route-specific flags; four preservation tests sat outside the default
run. Each was written in good faith and each one welded the code to one
implementation of itself.

So this is not a coverage target. It is a restriction.

## The diagnostic

Before writing any check, answer one question:

> **Would this test survive a legitimate reimplementation of the thing it tests?**

If someone rewrote the module from its equations, with different internals, a
different call sequence, different intermediate types — would the test still be
correct and still pass?

- **Yes** — it tests the claim. Write it.
- **No** — it tests this code. Do not write it. There is no version of it that is
  worth having, and "but it caught a regression once" is how the old tree got
  its 4,404 lines.

## The two categories that pass

**Oracle rows.** Real mixtures, real measured data, in `engine/tests/oracle/`.
NIST light-alkane saturation, Hamer-Wu mean ionic activity coefficients, the
Khudaida LLE tie-lines, the reactive VLE baseline. These are facts about the
world; no implementation can be right and fail them.

**Mathematical properties.** Facts about the equations that hold for every
parameter set: closed forms (Carnahan-Starling where the hard-sphere logarithm
cancels), exact symmetries (relabelling components changes nothing), limits
(every residual term vanishes at zero density), and derivative agreement (a
recorded tape against a centred difference). No stored numbers, no golden files.

## What is banned outright

- Asserting on an intermediate value, a struct field, a call sequence, or an
  iteration count as a **correctness** claim.
- Golden files of computed output. A number that came out of this code is not
  evidence that this code is right.
- Route-specific or downstream-specific flags and gates. A test that knows which
  solver path it is on is testing the path, not the physics.
- Any test that exists to raise coverage.

## Counts and timings are not correctness

`TRANSFER.md` gives each stage evaluation counts, recordings, Hessian calls and
wall times. **These are performance acceptance, not tests.** A correct but slower
implementation fails them, which is exactly what a correctness test must never
do. Report them, compare them to the floor, and treat a miss as a conversation.
Never make one red-first, and never let one block a change that is right.

## Red-first, and why it is not TDD

Before a stage is accepted, its qualifying checks must be shown to **fail on a
deliberately broken kernel**, then pass when it is restored. This is not
test-first development and it does not shape the design: the code is already
written when the mutation runs. It answers one question that a green suite
cannot -- *is this suite checking anything at all?*

It is not hypothetical. Stage 4's first suite passed while checking nothing: the
build sets `b_ndebug=if-release` with `buildtype=release`, `NDEBUG` erased every
`assert`, and `meson test` reported the physics suite green. The mutation caught
it; nothing else would have. See `STAGE-LOG.md`.

Because red-first applies only to checks that pass the diagnostic above, it
cannot produce the bloat it might in a coverage-driven suite -- there is nothing
to be thorough about. The oracle has 71 rows across seven cases and the properties are
finitely many.

## Budget

`TRANSFER.md` predicted 1,650 lines of tests for the whole engine, against 4,404 in
the tree it replaces. The measured figure is 1,926 after Stage 4 and about 2,400 once
Stage 5 lands (`EQUILIBRIUM.md` §8). The number is recorded, never trimmed to fit. If a
stage wants more than its share, the question is whether the checks pass the diagnostic,
not whether the budget should grow. A line count is never completion evidence in either
direction: a ceiling bounds a stage, a minimum proves nothing, and a check is accepted by
reading it against its named property and mutation.
