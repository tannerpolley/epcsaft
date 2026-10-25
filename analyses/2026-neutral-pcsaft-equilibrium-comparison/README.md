# Neutral propane PC-SAFT equilibrium comparison

This folder preserves the original cross-package run. The
[current Engine refresh](../2026-open-source-equilibrium-comparison/README.md)
rechecks the same case on wheel `4268f839…`: pressure is 998660.8957399124 Pa;
the original external rows below were not rerun. Historical build commands are
not required for the bounded Engine-only refresh.

This analysis compares five independent implementations on one deliberately
small equilibrium case:

- pure propane, no association, no polarity, no electrolyte terms;
- Gross–Sadowski 2001 PC-SAFT parameters (`m = 2.0020`,
  `sigma = 3.6184 Å`, `epsilon/k = 208.11 K`);
- liquid–vapor saturation at `T = 300 K`;
- reported quantities: saturation pressure and liquid/vapor molar density.

The result is an implementation comparison, not an experimental validation or
an external oracle. All five codes use the same neutral PC-SAFT physics and
the same parameter values; small differences therefore expose conventions,
rounding, and solver implementation details.

## Reproduction

Set `FORKS_ROOT` to the directory containing the package forks.

Before collecting values, build each fork using its own native path:

```bash
# Engine: use the repository's normal installed-wheel setup.
bash tools/check-all.sh

# FeOs
cd "${FORKS_ROOT:?Set FORKS_ROOT to the package forks directory}/feos"
cargo test --release --features pcsaft

# teqp (the checkout's declared submodules must be initialized once)
git submodule update --init --recursive
cmake -S . -B /tmp/teqp-local-build -DCMAKE_BUILD_TYPE=Release
cmake --build /tmp/teqp-local-build -j2
ctest --test-dir /tmp/teqp-local-build --output-on-failure

# ThermoPack
git submodule update --init external/pFUnit
cmake -S . -B /tmp/thermopack-local-build -DCMAKE_BUILD_TYPE=Release -Dtest=ON
cmake --build /tmp/thermopack-local-build -j2
cmake --install /tmp/thermopack-local-build

# Clapeyron: use the local project and a Julia installation.
julia --project="$FORKS_ROOT/Clapeyron.jl" -e 'using Pkg; Pkg.test()'
```

The package wrappers used by `scripts/collect_results.py` must be built from
the local FeOs, teqp, and ThermoPack checkouts into the same temporary Python
environment. The exact environment paths are supplied through
`COMPARISON_PYTHON` and `CLAPEYRON_JULIA`; the collection script does not
silently substitute PyPI packages.

```bash
COMPARISON_PYTHON=/path/to/python \
CLAPEYRON_JULIA=/path/to/julia \
.venv/bin/python analyses/2026-neutral-pcsaft-equilibrium-comparison/scripts/collect_results.py

.venv/bin/python analyses/2026-neutral-pcsaft-equilibrium-comparison/scripts/render_figure.py
```

`results/values.csv` is the exact plotted table. `results/provenance.json`
records the source commits, dirty-state flags, commands, and numerical
outputs. `figures/neutral-propane-pcsaft-comparison.png` and `.svg` are
rendered only from that retained table.

## Collected values

| implementation | pressure [kPa] | liquid density [kmol m⁻³] | vapor density [mol m⁻³] |
| --- | ---: | ---: | ---: |
| ePC-SAFT | 998.660895740 | 11.100251233210 | 482.512126815 |
| FeOs | 998.660895548 | 11.100251234177 | 482.512126702 |
| ThermoPack | 998.660895733 | 11.100251233209 | 482.512126811 |
| teqp | 998.660895740 | 11.100251233210 | 482.512126815 |
| Clapeyron | 998.660894643 | 11.100251233072 | 482.512126157 |

The largest spread across implementations is `1.10 × 10⁻⁹` relative in
pressure and `1.36 × 10⁻⁹` relative in vapor density. The Engine solve reports
`max(abs(residuals)) = 1.834 × 10⁻¹³`; FeOs and teqp report liquid/vapor
pressure mismatches of `1.77 × 10⁻¹³` and `1.03 × 10⁻¹⁴`, respectively.

## Verification evidence

- FeOs local source: `cargo test --release --features pcsaft`; 76 tests passed
  (47 unit, 28 integration, 1 doctest). Its Python wheel was built from the
  local checkout.
- ThermoPack local source: optimized Fortran library and vendored pFUnit
  built through CMake; the direct pFUnit runner passed 126 tests. The legacy
  `make unittests_optim` path was also attempted but requires `funitproc` on
  `PATH`.
- teqp local source: CMake build completed and the focused compiled Catch2
  filter `[PCSAFT]` passed 42 assertions in 8 test cases. The unfiltered CTest
  invocation needs the checkout's expected relative `fluiddata` working
  directory and was therefore not used as the acceptance result.
- Clapeyron local source: the direct `PCSAFT("propane")` saturation call
  passed. The long full `Pkg.test()` attempt reached and passed the package's
  PCSAFT sections (`218/218` and `76/76`) before the session terminated during
  later tests; no whole-suite pass claim is made.

## Scope limits

This is a pure-fluid saturation point, so it does not test mixture phase
splits, electrolyte chemistry, association, reaction equilibrium, or the
Engine's full derivative-driven multiphase application layer. It is useful as
a neutral cross-implementation consistency case and as a first regression
oracle, not as evidence that all five packages have equivalent capability.

## Named-release rerun (#680)

The [research notebook](notebook.qmd) reads the new results under
[results/main-8298c602/](results/main-8298c602/). The earlier retained files remain
under their original identities. The new provenance records the downloaded
wheel, its verified hash, the actual installed Engine identity and the declared
producer. Unavailable portions retain their stopping states and causes.
Independent delivered review and owner promotion remain open; the manuscript
is unchanged.
