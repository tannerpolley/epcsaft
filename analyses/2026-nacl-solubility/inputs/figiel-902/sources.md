# Binary halite input basis

The frozen observations and source hashes are exact copies of the lead's
`build/issue-902/frozen-observations.json` and `source-hashes.json`.
Their unrelated neutral and hydrate entries remain unexecuted provenance.
The independent readiness review inspected all source hashes and the halite
Tables 17A and 34 in Clarke–Glew (1985), doi:10.1063/1.555730, and
Pinho–Macedo (2005), doi:10.1021/je049922y, ThermoML dataset 1 row 1.

The liquid selection is water, sodium-cation, chloride-anion from
`data/packets/figiel-2025-chloride-mixture-refit/2/parameters/parameters.json`.
The packet supplies the water molar mass, 0.0180153 kg/mol. Explicit salt
molar mass is 0.05844277 kg/mol, the existing Pinho conversion in
`../sources.md`. The selected packet has no ionic molar masses; they are
unnecessary for the salt/water mass-fraction calculation.
Solid molar volume is 27.1e-6 m³/mol (the frozen USGS pitzer.dat input).
The reference pressure is 101325 Pa, independent of the system pressure.
Standard molality is 1 mol/kg water; the solid has unit activity.
No coefficient or standard is fitted.

`results/figiel-902/run-identity.json` identifies the exact input, executable
source, generated declaration, loaded extension and native shared-library bytes.
The generated C++ header transports the Engine's own ParameterDocument entries;
it contains declarations only. Python performs no thermodynamic calculation.
All predicted quantities, source concentration conversions, errors, signed
perturbations and numerical falsifiers are evaluated in the native study.

Starts have liquid composition corresponding to 1 and 12 mol/kg and use the
existing driver's free solid amount initialization. Each central pressure has
both starts at both feeds. Each separate ΔG or temperature perturbation has
both starts at the 12 mol feed. Perturbations are finite separate realizations,
not a covariance estimate or a joint confidence interval. No absent-solid
endpoint occurs in this excess-feed set; the absent-solid inequality is
implemented conditionally but is not independently exercised here.
