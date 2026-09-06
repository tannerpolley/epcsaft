# Engine context

The [scientific context](../docs/scientific/CONTEXT.md) owns the vocabulary and
scientific authority map. Read it and the relevant accepted ADR before changing
scientific behavior. Engine owns one runtime and one non-editable wheel;
applications own chemistry selections, fit adoption and application evidence.

Read [equations](docs/equations.md), [EOS conventions](docs/science/full-eos.md),
[equilibrium formulation](docs/equilibrium/formulation.md) and
[algorithms](docs/science/algorithms.md) for required mathematics. Read the
[current transfer plan](../docs/greenfield/TRANSFER.md) for implemented,
incorrect, missing and unqualified behavior. A definition in the equation
index does not establish that the installed Engine implements it.

Existing research documents and dated evidence retain their original source
locators and historical software names. They inform scientific choices but do
not replace the current specifications. Numerical policy values stay with
`native/core/numerical_policy.hpp`; code and numerical evidence must be updated
or revalidated together with the affected specification, as described in the
[scientific records guide](../docs/scientific/README.md).
