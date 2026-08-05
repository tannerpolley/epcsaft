# Test suite index

- `numerical/`: equation contributions, thermodynamic invariants, density closure, and
  independent finite-difference checks.

Tests neither locate nor import a sibling repository. Native transport is an
Engine-internal implementation seam; cross-module scientific evidence belongs
in the Equilibrium, Regression, and integration suites rather than ABI/layout
tests.
