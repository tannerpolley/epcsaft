## Target layout (authoritative for #186–#192)

**Principle.** C++ is organized by dependency layer, because the only structural fact the compiler and the checker can enforce is include direction. Python is organized by capability (`eos`, `equilibrium`, `regression`), because that is what a user imports. The two do not mirror each other, deliberately: a layer answers "what may this include", a capability answers "what does a caller reach for". Tests and docs follow the Python side (capability), because a test or a document is written against a user-facing claim, not against a compilation unit.

Consequences, stated once:

- `_native/` dies as a name and as a location. All C++ moves under one root, `engine/native/`, outside `src/epcsaft/`. `meson.build` loses `exclude_directories: ['_native']`; `install_subdir('src/epcsaft')` ships only Python; the wheel contents are unchanged (`epcsaft/*.py`, `epcsaft/equilibrium/*.py`, `epcsaft/_core.so`).
- One static library, `epcsaft_native`, from every `native/**/*.cpp` except `bindings/`; one extension, `_core`, from `native/bindings/*.cpp`, linked against it. The three-archive split and its link-order rule go.
- Bindings live in one directory, `native/bindings/`, not beside the code they bind: they are the only files permitted to include pybind11, and one directory makes that a path rule the checker can test. The earlier proposal of `equilibrium_bindings.cpp` beside `module.cpp` in two places was two answers; this is one.
- Saturation is not a file. `solve_saturation` is the one-component, two-phase call of `solve_phase_boundary` (formulation.md §7), so it lives in `native/equilibrium/phase_boundary.cpp`; `Pure.saturation` in Python calls that binding. `helmholtz_tape.hpp` and `density.cpp` live in the `model` layer (EOS state without any equilibrium concept); `phase_boundary.cpp` in the `equilibrium` layer. The four files that were in two trees are now in two *layers* of one tree with a one-way include.
- Permitted include direction: `physics` ← `model` ← `equilibrium` ← `regression` ← `bindings`, downward only (`regression` may include `model` and `equilibrium`; `bindings` may include anything; nothing includes `bindings`). Python includes nothing; it imports `_core`.
- Regression rule, by symbol rather than by module: `regression/` may include `equilibrium/problem.hpp` and call `solve_equilibrium` and `solve_phase_boundary` on a problem with declared phases; it may not include `phase_stability.hpp`, call `find_phases`, or construct a problem with phases omitted. This replaces today's module-level ban (`_regression_discovery_family` in the checker; `regression_core: []` in `ARCHITECTURE.yaml`), which the layout contradicted. The reconciliation, the move, the single archive, and the gate routing are owned by #195, sequenced before #187.

**What the checker validates after the change** (`tools/check-engine-architecture.py`, replacing the link-order and three-archive rules): exactly one static archive and one extension in the Meson products; `nm` on the archive shows no `PyInit`/pybind11 symbols; every `#include "..."` in `native/<layer>/` resolves to the same or a lower layer (a 40-line scan over the include lines); no `constexpr` name defined in two translation units and no category-1 constant outside `physics/physical_constants.hpp` (#186); `ARCHITECTURE.yaml` `native_core.static_libraries` becomes one entry and `permitted_internal_links` becomes the layer order above.

```text
engine/
  meson.build                                   120   one library, one extension
  native/
    physics/                                    950   no CppAD types, no solver
      physical_constants.hpp                     30   category 1, defined once
      numerical_policy.hpp                      250   categories 2–4, Tolerance entries (#186)
      helmholtz.hpp                             200   hard chain, dispersion, packing volume closed form
      association.hpp                           150   strength, site equations, Newton
      electrolyte.hpp                           260   Debye–Hückel; four permittivity forms
      born.hpp                                   60
    model/                                     1810   EOS state and exact derivatives; no phase concept
      helmholtz_tape.hpp                        350   one cached ADFun per topology; value, gradient, Hessian, parameter Jacobian
      helmholtz_model.hpp                       150   HelmholtzModel: parameter arrays, domains, parameter_hash
      helmholtz_model.cpp                       600   construction, validation once, state evaluation
      density.cpp                               250   bracketed scan + Newton on the pressure row for one phase
      infinite_dilution_reference.cpp           450   dilution sequence, MIAC, osmotic coefficient
    equilibrium/                               5270   one minimization, its rows, its stability
      problem.hpp                               300   EquilibriumProblem, EquilibriumResult, specification patterns
      stoichiometry.hpp                         120
      stoichiometry.cpp                         800   balance/reaction rank (SVD), standard-state transform, nullspace
      request_compile.cpp                       600   specification → active rows; deterministic starts; count check
      gibbs_minimizer.hpp                       150
      gibbs_minimizer.cpp                      1200   GibbsMinimizer (Ipopt), residuals, reduced-Hessian check
      phase_boundary.cpp                        300   solve_phase_boundary: T or P unknown; one-component case is saturation
      implicit_sensitivity.cpp                  250   KKT system, rank, condition, output rows
      phase_stability.hpp                       150
      phase_stability.cpp                      1400   find_phases: electroneutral chart, TPD search, dual bound loop
    regression/                                 890
      least_squares.hpp                          80
      least_squares.cpp                         410   Ceres kernel: ResidualBlock, PredictionCache, covariance
      eos_observable.cpp                        150   state → value and Jacobian rows
      equilibrium_observable.cpp                250   problem → value and Jacobian rows
    bindings/                                  1100   the only pybind11 includes
      module.cpp                                150   PYBIND11_MODULE(_core), submodule registration
      eos_bindings.cpp                          400   HelmholtzModel, State, Saturation, Reference as py::class_
      equilibrium_bindings.cpp                  400   Problem, Result, solve_equilibrium, find_phases, solve_phase_boundary
      regression_bindings.cpp                   150   fit, Observation, Fit
  src/epcsaft/                                 2430   Python: maps inputs to _core and _core results to types
    __init__.py                                  60   11 exports
    _core.pyi                                   250
    units.py                                     80
    parameters.py                               400   Parameters.from_packet / from_mapping; ActiveParameters
    eos.py                                      600   Pure, Mixture, State, derivatives, saturation
    reference.py                                100   InfiniteDilutionReference mapping only
    regression.py                               200   fit, Observation, Fit (a module: one capability, one file)
    equilibrium/
      __init__.py                                40   9 exports
      problem.py                                450   Problem, Phase, Reaction, Fixed, Solved, Result
      loader.py                                 250   problem_from_mapping; mass-fraction and molality feeds
  tests/                                       1650   by capability; no numerical/behavior split
    packets.py                                  150   read data/packets tables for fixtures
    fixtures/                                   250   two synthetic tables (immiscible ternary, associating binary)
    eos/                                        520   saturation vs NIST; γ± vs Hamer–Wu (Held, Figiel); closed-form, identity, limit, permutation, closure, anchor
    equilibrium/                                530   LLE vs Figure 2 tie-lines; analytic uptake and reaction; conservation; boundary Jacobians; discovery; values survive failure; validation
    regression/                                 200   NIST parameter recovery; Jacobian vs centered differences; permutation
  docs/                                         —    by capability; unchanged files not listed
    science/full-eos.md                          —
    equilibrium/formulation.md                   —    reverted 2026-09-02; skeleton lattice being designed
    equilibrium/phase-stability.md               —    with #188
    latex/, equations_registry.yaml, equations.md, python-api.md, scientific/adr/   —
  scripts/docs/, tools/                          —    unchanged
```

Totals: native 10,020; Python 2,430; tests 1,650; build 120 → **14,220**, inside the ~14.5k envelope.

**Deleted directories:** `engine/src/epcsaft/_native/`, `engine/cpp/`, `engine/cpp/interfaces/`, `engine/tests/*/numerical/`, `engine/tests/*/behavior/`. No new directory holds one file: `physics` 6, `model` 5, `equilibrium` 10, `regression` 4, `bindings` 4, `equilibrium/` (Python) 3, `tests/fixtures` 2.

**Homes I could not decide (`investigate`):**

- `numerical_policy.hpp` in `physics/`: every layer reads it, and `physics/` is the lowest, so it is reachable; but it is not physics. If #186 prefers a `native/policy/` directory it would hold one file, which the constraint forbids; the alternative is `native/numerical_policy.hpp` at the root of the tree, outside any layer. Decide in #186.
- The closed-form packing volume: 60 lines in `helmholtz.hpp` versus its own header. Folded into `helmholtz.hpp` here; #189 may split it if the tape needs it without the rest of the Helmholtz terms.
- `infinite_dilution_reference.cpp` in `model/`: it runs a dilution sequence of single-phase states, so it needs nothing from `equilibrium/`; if the target later expresses the reference as a boundary problem it moves up one layer.
- `stoichiometry.cpp` standard-state transform: today `transform_source_standard_state` reads EOS reference quantities; if that pulls `infinite_dilution_reference` into the equilibrium layer the include is still downward, so no move, but the line budget for `stoichiometry.cpp` is the least certain (600–900).
- `tests/packets.py` versus a `Parameters.from_packet` call in each test: if the packet loader is thin enough, this helper disappears and its 150 lines go with it.
