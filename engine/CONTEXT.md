# Engine Context

`governance_doctrine_revision: 11`

`engine/` is the sole runtime owner inside `tannerpolley/ePC-SAFT`. It
builds one Python 3.13 distribution, `epcsaft`, containing EOS, equilibrium,
and regression. All native modules share one CMake graph and wheel.

Reusable packets remain owned by `data/` and are materialized with exact
commit, packet, fingerprint, and file hashes. Regression transfers accepted
values through typed `RegressionResult.to_parameters()`; EOS and equilibrium
consume the resulting `Parameters` directly. No package, SDK, filesystem, or
legacy compatibility handoff exists between the three runtime modules.

Engine retains compact numerical evidence. Installed-wheel physical and
literature campaigns belong in `validation/`; independent or manufactured
oracles belong in `lab/`. Historical plans and receipts preserve former
repository and package identities as provenance.
