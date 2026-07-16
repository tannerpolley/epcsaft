# Provider Repository Agent Contract

Repository Profile: scientific-computing

Before work, read doctrine revision 2 at
`../ePC-SAFT-organization/GOVERNANCE.md`. Its future published home is
`https://github.com/ePC-SAFT/.github/blob/main/GOVERNANCE.md`. Local policy may
narrow that doctrine but must not contradict it.

The current functional boundary is the source-backed parameter contract and
neutral, nonassociating, nonionic PC-SAFT evaluation at explicit molar density.
The package root has six public symbols: `ParameterBundle`, `ParameterSet`,
`ParameterError`, `unit_registry`, `EPCSAFT`, and `EosResult`. Native C++
owns the hard-chain and dispersion equations; CppAD owns the density derivative
consumed by pressure; Pint remains the only Python runtime dependency.

This repository is limited to provider parameters, EOS/state/property
evaluation, density closure, CppAD, and a separately admitted resolved-input
SDK. It does not own Ipopt, Ceres, equilibrium dispatch, regression workflows,
or cross-package validation. Specified-pressure density closure is the next
functional slice and requires its own approved implementation plan.
