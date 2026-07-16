# ePC-SAFT Provider

This package currently provides a source-traceable parameter-input contract.
Its root surface is exactly `ParameterBundle`, `ParameterSet`, `ParameterError`,
and `unit_registry`. Load an exact bundle folder or construct immutable typed
records as a `ParameterBundle`, then select exact component IDs to obtain an
ordered `ParameterSet`.

The installed catalog contains only the Gross and Sadowski 2001 methane/ethane
parameter bundle. The compact Figiel/Yu/Held 2025 and Ascani/Held 2021 data in
`tests/data` is a disconnected parser fixture: it is not installed, cataloged,
or a physically complete mixture.

There is no EOS, model, state, derivative, solver, or thermodynamic capability
in this slice, and no thermodynamic runtime authority has transferred from the
lab archive.
