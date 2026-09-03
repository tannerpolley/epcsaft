# Engine scientific context
## Question and intended use

The Engine evaluates source-traceable ePC-SAFT states, configured equilibrium,
and owner-provided regression observables through one installed `epcsaft`
wheel. The intended use is generic scientific computation for the admitted
domains and formulations recorded by the existing owners below. Application
chemistry, parameter adoption, and downstream claims remain outside Engine.

## Scientific map

- Engine terminology and capability boundary: `engine/CONTEXT.md`
- Stable semantic facts: `engine/ARCHITECTURE.yaml`
- Formulation authority: `engine/docs/science/full-eos.md`
- Equation registry and synchronized views: `engine/docs/equations_registry.yaml`
  and `engine/docs/latex/`
- Public API guide: `engine/docs/greenfield/TRANSFER.md` (replaced the deleted plan)
- Phase-equilibrium method and certificates:
  `engine/docs/greenfield/TRANSFER.md` (replaced the deleted plan)
- Accepted scientific decisions: `engine/docs/scientific/adr/`
- Compact numerical falsifiers: `engine/tests/`
- Installed-wheel analyses: `analyses/`
- Accepted decisions: `engine/docs/scientific/adr/`

## Source authority

Primary-source equations, conventions, units, and validity bounds are owned by
the formulation specification and equation registry. Data packets own reusable
source values and provenance. Accepted ADRs own consequential choices and
rationale. Runtime source owns concrete state lifetimes and implementation.
Tests and installed-wheel campaigns own numerical claims.

## Claim boundaries

Generated documentation and architecture views are derived evidence. They may
expose structural drift but do not prove numerical correctness, physical
validity, complete runtime reachability, or a downstream scientific claim.
