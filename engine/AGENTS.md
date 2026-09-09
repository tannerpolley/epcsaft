# Engine-specific instructions

Follow root `AGENTS.md`. Before changing Engine behavior, read the relevant
canonical mathematics and algorithms mapped in `../docs/scientific/CONTEXT.md`;
`../docs/scientific/implementation-plan.md` distinguishes required capability
from installed availability.

`engine/` owns the sole `epcsaft` runtime. Preserve one Python build, one Meson
graph, one version, one Python-free native archive and one extension. Follow
`ARCHITECTURE.yaml` for dependency direction. Native iterative work has no Python
callbacks; Python maps inputs and results. Do not add compatibility imports,
alternate package identities, hidden physical defaults or a public C++ ABI
without an approved scope change and a real external consumer.

Runtime and tests never discover sibling source dynamically. The monorepo commit
binds copied Data fixtures to their packet lock; each materialization also retains
the packet fingerprint and exact file hashes. Analysis reproductions exercise
the installed wheel through public interfaces.

Numerical checks may cover manufactured solutions, symmetries and justified
solver/equilibrium certificates. Internal layouts, solver labels, call sequences
and iteration counts are not numerical correctness criteria; evaluation counts
and timings are performance evidence for equivalent accepted work.
