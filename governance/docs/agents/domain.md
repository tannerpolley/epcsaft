<!-- project-truss:setup:start -->
## Project Truss domain contract

Layout: `single-context`.

Use one root `CONTEXT.md` and root `docs/adr/`.

Read the nearest applicable context before shaping or implementation. Keep durable terminology, invariants, ownership, and accepted decisions there; do not mirror lifecycle state.

<!-- project-truss:setup:end -->

## Runtime terms

Doctrine revision 8 uses these runtime terms:

- **Engine**: the single `epcsaft` distribution and repository that owns EOS,
  Equilibrium, and Regression production modules.
- **Vertical slice**: one task-owned capability change from scientific input
  through Engine implementation, Engine tests, installed-wheel Validation, and
  affected downstream applications.
- **Transition source**: one of the current EOS, Equilibrium, or Regression
  repositories before the D-023 cutover receipt. A transition source owns no
  post-cutover release line.
- **Consumer environment**: an application-owned environment that installs the
  immutable Engine wheel. It is not an SDK source or shared project runtime.
