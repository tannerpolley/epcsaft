# Native Engine ownership and generated evidence
Kind: computational-architecture

- Status: superseded in part by the greenfield rebuild, 2026-09-03
- Date: 2026-08-27
- Amended: 2026-09-01 for issues #156 and #153

> **Superseded in part, 2026-09-03.** Three decisions below no longer hold on `greenfield`:
> the three Python-free libraries are one archive, `epcsaft_native` (`docs/greenfield/CANON.md`,
> rule 1); HELD and HELD2 are not retained, and `find_phases` is written fresh to Pereira 2012
> (`docs/greenfield/PHASE-STABILITY.md`, Stage 7); `cpp/interfaces/owner_evaluator.hpp` no
> longer exists, and regression reaches the engine through `equilibrium/problem.hpp` only.
> What survives: one extension boundary, Python maps and never computes, no Python inside the
> archive, no iterative Python callbacks, and the lifetime and thread decision with its evidence.

## Context

The Engine formerly separated EOS and Equilibrium across private extension
modules and transported EOS operations through a C table and capsule. The
accepted direct-C++ cutover removed that surface. Regression then carried a
second thermodynamic/IPOPT lifecycle beside its owners. Retaining separate
extensions, a renamed function table, a public C++ ABI, iterative Python
callbacks, or generated documentation as authority would recreate the
displaced ownership.

## Decision

`epcsaft._core` is the sole Python/native extension boundary. Python validates
and maps one owner request, crosses that boundary once, and maps one final
result. The Python-free native libraries are `eos_core`, `equilibrium_core`,
and the Ceres-only `regression_core`. Equilibrium may consume EOS;
Regression consumes only `cpp/interfaces/owner_evaluator.hpp`; `_core` consumes
all three. Reverse and sibling edges are forbidden. Python C-API, pybind11,
GIL, capsule, private-transport ABI, and iterative Python callback edges are
forbidden inside the three libraries.

Within that graph, immutable `epcsaft::NativeModel` is the sole EOS evaluation
owner. `eos_core` validates EOS domains and metadata, selects pure, neutral,
associating, or electrolyte kernels, returns named physical EOS quantities and
requested exact partials, performs source-reference transfer, and reports
primitive evaluation counts. Equilibrium combines those named results into
HELD, HELD2, Problem (67), GREPE, conservation, KKT, continuation,
certification, and solved-state sensitivity calculations; it does not own a
second EOS façade, route selector, physical projection, or mutable EOS counter.

The compiled GREPE request structure is source-linked and lives for one public
solve or one distinct concrete Regression evaluator. Owner source values and
partials are evaluated once per exact optimizer coordinate. Direct physical
parameter coordinates, bounds, and diagnostics remain the optimizer
representation. Solver state lives within its solve, and a fresh final owner
value/Jacobian validation bypasses trial caches. Source files remain the
concrete lifetime authority.

HELD and HELD2 remain distinct retained phase-discovery algorithms inside the
one Equilibrium owner. No diagram, API page, or other generated view may merge,
rename, or supersede them.

Generated diagrams and API pages are disposable derived evidence. Meson owns
the configured native target graph, source owns implementation and lifetimes,
accepted ADRs own decisions and rationale, and executable tests and campaigns
own numerical claims.

The lifetime and thread decision remains no-change. The accepted evidence is
the measured result at
https://github.com/tannerpolley/ePC-SAFT/issues/86#issuecomment-5442040189:
the unchanged installed wheel passed three reproductions with a 1.93 s median
wall time; mutex wait was effectively zero; tape construction was smaller than
the observed outlier excess; all temporary instrumentation was removed.

## Verification and validity boundary

`tools/check-engine-architecture.py` checks the configured Meson products,
backend archive link order, static-library symbol fields, production source
identifiers, install plan, and identified installed wheel. Existing numerical
falsifiers preserve values, derivatives, certificates, tolerances, units, and
HELD/HELD2 behavior. These checks expose declared structural drift; they do not
prove numerical correctness, physical validity, dynamic call count, or complete
runtime reachability.
