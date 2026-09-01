# Native Engine ownership and generated evidence
Kind: computational-architecture

- Status: accepted
- Date: 2026-08-27
- Amended: 2026-08-31 by Governance D-041

## Context

The Engine formerly separated EOS and Equilibrium across private extension
modules and transported EOS operations through a C table and capsule. The
accepted direct-C++ cutover removed that surface, then made the native
Equilibrium and Regression libraries own their complete calculations. Retaining
separate extensions, a renamed function table, a public C++ ABI, iterative
Python callbacks, or generated documentation as authority would recreate the
displaced ownership.

## Decision

`epcsaft._core` is the sole Python/native extension boundary. Python validates
and maps one public request, crosses that boundary once, and maps one final
result. The Python-free native libraries are `eos_core`, `equilibrium_core`,
and `regression_core`. Equilibrium may consume EOS; Regression may consume EOS
and Equilibrium; `_core` consumes all three. Reverse and sibling edges are
forbidden. Python C-API, pybind11, GIL, capsule, private-transport ABI, and
iterative Python callback edges are forbidden inside the three libraries.

The compiled GREPE request structure is source-linked and lives for one public
solve or one distinct concrete Regression experiment structure. Retained
source-block structure lives with that compiled Regression experiment;
source-standard-state values and partials are evaluated once per exact optimizer
coordinate. Direct physical parameter coordinates, bounds, and diagnostics
remain the optimizer representation. Solver state lives within its solve;
simultaneous full-space Regression retains each accepted experiment state and
finalizes its values, sensitivities, certificates, and receipt once after
convergence. Source files remain the concrete lifetime authority.

One source-neutral Phase-Stability Discovery sequence owns omitted-topology
nonreactive discovery. Its retained Stages 1--7 prepare and qualify neutral or
electroneutral-electrolyte candidate phase cards, then GREPE is the sole
configured finite-phase Gibbs/KKT minimizer for each admitted candidate
topology. Two-phase discovery remains the maximum unless the independent
three-basin oracle and all three-phase gates pass. GREPE's separately admitted
declared three-phase capability is unchanged. Declared topology bypasses discovery and
enters GREPE directly without fallback. HELD and HELD2 remain historical
lineage and citations, not live runtime identities.

The public request and resolved Provider inputs own every native object needed
for one solve. The compiled zero-reaction request, Provider context, candidate
cards, and explicit GREPE starts remain value-owned and alive through all
subset solves, physical audits, and fixed-topology sensitivity evaluation.
References never escape that scope. Only accepted result values and structured
failure evidence cross the native/Python boundary through one result mapper.
The later EOS-access cutover may replace the Provider access vocabulary,
concrete object, and internal lifetime arrangement. It must still keep every
object alive for the complete request-owned solve and cannot change discovery
equations, budgets, candidate decisions, certificates, result values, or allow
a reference to escape the solve scope.

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
phase-stability discovery behavior. These checks expose declared structural drift; they do not
prove numerical correctness, physical validity, dynamic call count, or complete
runtime reachability.
