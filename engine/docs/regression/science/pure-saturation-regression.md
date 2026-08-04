# Pure-Saturation Regression

The package fits segment count, segment diameter, and dispersion energy for
one closed pure-component record per call. Methane and ethane are accepted
reproducible workflows. Propane Checkpoint A is retained authority-neutral,
non-production evidence. Migration D-023 closes its 120 K admission effort as
an active gate; the pressure-resolution protocol below remains provenance only.
The installed `epcsaft` provider owns Helmholtz energy and nonlinear
derivatives.

## Source records

The methane table contains NIST Chemistry WebBook SRD 69 saturation properties
from 100 through 180 K in 10 K increments. The fit uses 110, 130, 150, and
170 K. The other five rows remain descriptive held-out evidence. The retained
lab file has SHA-256
`a5e16df3bf8ec78483fc340782cddc89ab8b284a9f6dfaecd6cda3ffde579227`.
The package normalizes CRLF to LF and has SHA-256
`dec64d5a6cac414a4a92393a0d728fa27c02135c6a159d0d1881d7b6dde6d26c`.

The ethane table contains NIST SRD 69 saturation properties from 100 through
280 K in 20 K increments. The fit uses 140, 180, 220, and 260 K. The held-out
set uses 120, 160, 200, and 240 K. The package reports 100 and 280 K as domain
stress. The retained CRLF file has SHA-256
`ed09b8781acfb7025ca505878b884f6353ddd9f3f4bd7aae2e6df88bbe847a67`.
The LF-normalized package file has SHA-256
`b01333e827933c0a7148672c8ae3eef78393320c0d18f2c4d5a0fc40d9bef6b2`.

The propane table is the exact 24-row target from the approved Glos,
Kleinrahm, and Wagner 2004 direct-experimental packet at Validation commit
`7e51590757f1cb85f51df98e9fe1f88cd4255a88`, tree
`05af9e948c786ddfcf43dba701970f1cbb6435a2`. The packaged target SHA-256 is
`ccd1cfa15ec44432b06cbf22316d168c61b282631c9b1e1591e497b8d48b5676`;
the governing packet YAML is
`ba31448989f565d05d63908076e836977780aa87199f208310e9b80b03f64697`;
and the 63-row source receipt is
`ed5eb703ccd3e6bb4c4cfa82ecd58c58f9da0c93ab07a204dee94d8b0ae8d081`.
Training uses 150, 210, 270, and 330 K. The 18 remaining interior rows are
held out, while 110 and 340 K are stress rows. Pointwise pressure, liquid-
density, and available vapor-density uncertainties are retained on the source
observations. They do not weight residuals or define model-acceptance cutoffs.

Both tables keep temperature in K, saturation pressure in Pa, and saturated
liquid mass density in kg/m3. The records retain the exact NIST query, decimal
strings, retrieval date, and use basis.

## Provider contract

New installed-artifact replay binds Provider frontend `0.2.0.dev0` to commit
`14fa3745264db66b8e59c12268737d694c706f2f`, tree
`eb04f10f445957cc768bad1ef4f330038c69a293`, wheel SHA-256
`48f3a75c9fc16ba71616aa703b526f41c2dcf89a7e00eebe23f75fcb8fa24594`, and
installed-header SHA-256
`2cd2b73b83c65936dff21155fd800a87b56e81cce977df7b8491ccfb2bf4c50b`.
Its Provider test receipt artifact is
`provider-python-frontend-v0.2/14fa374/provider-tests.xml`,
SHA-256
`56d1d3e9fcf47bb700df4223fa2d9a20444dc97aa22b5c39525d446a296ba3cf`.
The artifact name is relative to an executor-supplied immutable artifact root.
The older Provider artifacts referenced by the retained propane diagnostic are
historical evidence and are not active replay inputs.

Each phase evaluation fixes `n = 1 mol` and sends

```text
(n_mol, V_m3, m, sigma_angstrom, epsilon_over_k_kelvin)
```

to `evaluate_pure_phase_parameters`. The provider returns
`Phi = A/(RT n_ref)`, its five-coordinate gradient and Hessian, pressure,
chemical potential divided by `RT`, and its immutable parameter fingerprint.
Regression evaluates

```text
dP/dq_j       = -R T Phi_(V,q_j)
d(mu/RT)/dq_j = Phi_(n,q_j)
```

and owns no EOS equation.

## Component specifications

Methane keeps the accepted start and scale:

```text
start = (1.08, 3.555744 angstrom, 157.5315 K)
scale = (0.1, 0.1 angstrom, 10 K)
```

Ethane starts from the Gross 2001 provider parameters:

```text
start = (1.6069, 3.5206 angstrom, 191.42 K)
scale = (0.1, 0.1 angstrom, 10 K)
```

Propane starts from the provider's source-backed Gross 2001 propane bundle:

```text
start = (2.0020, 3.6184 angstrom, 208.11 K)
scale = (0.1, 0.1 angstrom, 10 K)
provider fingerprint = sha256:1194db349d0608c89419e70c56ccec9ada2ae0884dd8e64e519e9560e7e8ae42
```

Both specifications use bounds `[0.5, 3.5]`, `[2, 5] angstrom`, and
`[50, 400] K`. Methane fixes molar mass at `0.016043 kg/mol`; ethane fixes it
at `0.030070 kg/mol`; propane fixes it at `0.044096 kg/mol`. Each
specification also fixes its source, dataset, training partition, provider
fingerprint, volume bounds, pressure bounds, solver controls, and confirmation
gate.

The ethane vapor-volume upper bound is `100 m3` and its reporting-pressure
lower bound is `1 Pa`. Those values contain the audited low-temperature rows.
They do not relax residual closure or topology thresholds.

Propane uses liquid-volume bounds `[2e-5, 1.2e-4] m3`, vapor-volume bounds
`[1.5e-4, 2000] m3`, and reporting-pressure bounds `[0.1, 1e7] Pa`. These
domains contain every source-derived liquid-volume and ideal-gas vapor-volume
start, including the 0.6 Pa stress row. They are solver domains, not accuracy
tolerances. The propane iteration budget is 5000 because the unchanged Ceres
tolerances converge at iterations 1090 and 1165 for the primary and perturbed
starts; the accepted methane and ethane budgets remain 500 and their numerical
results are unchanged.

## Fit formulation

The physical parameters use affine dimensionless coordinates:

```text
p = start + scale*z
```

Each training row carries liquid-like and vapor-like volumes. Log transforms
make both positive. Disjoint volume bounds enforce liquid volume below vapor
volume. Observed liquid density and ideal-gas volume provide starts.

Each row contributes liquid pressure, vapor pressure, liquid-vapor chemical
potential equality, and liquid mass-density residuals. Pressure uses observed
pressure as scale. Chemical potential uses scale 1, and density uses observed
density. Each residual receives weight 0.25 for numerical row normalization.
The weights make no uncertainty claim.

The native cost assembles the `16 x 11` Jacobian from the provider Hessian and
the affine, logarithmic, density, and scaling chain rules. An installed-capsule
test independently reconstructs every methane and ethane Jacobian entry from
the Provider coordinate order `(n,V,m,sigma,epsilon/k)`, the returned `5 x 5`
Hessian, the analytic pressure and chemical-potential identities, the
log-volume chain rule, density derivative, and structural zeros. Production
exposes no numerical derivative option.

## Acceptance layers

Solver convergence requires Ceres `CONVERGENCE`, a usable finite solution,
complete Jacobian columns, full rank 11, projected parameter rank 3, parameter movement
inside the declared bounds, and no provider error. A confirmation solve
perturbs all eight nuisance-volume starts and requires agreement of the three
fitted parameter coordinates below `1e-5` plus relative cost agreement below
`1e-8`. Nuisance volumes are not required to agree coordinate by coordinate.
Cost agreement is the
symmetric relative difference `abs(primary - confirmation) /
max(abs(primary), abs(confirmation), numeric_limits<double>::min())`, so the
gate remains relative when both fitted costs are below one.

Physical acceptance requires positive stability slopes, ordered distinct
phases, relative phase-volume separation above `1e-3`, and reporting closure
on training and held-out rows. Pressure closure uses `1e-8` after scaling by
observed pressure. Chemical-potential closure uses `1e-8`.

Stress rows retain their termination, closure, topology, predictions, and
failure reasons. They do not enter solver or physical acceptance and do not
receive aggregate errors unless the reporting state passes the physical
gates. In the accepted ethane workflow, 280 K passes. The 100 K solve terminates
at the `1 Pa` lower pressure bound and fails pressure and chemical-potential
closure, so the receipt records no valid 100 K prediction.

For propane, the primary and perturbed training solves converge with full rank
and agree inside the unchanged confirmation gates. The 110 K stress solve
reaches the vapor-volume ceiling and has no acceptance effect. The 120 K
held-out reporting solve converges and is usable, but its liquid-pressure
closure is approximately `1.0540e-7 Pa`, or `3.29e-8` after division by the
observed `3.2 Pa`; this exceeds the frozen `1e-8` numerical closure gate.
The exact installed-callback diagnostic reproduces the retained pressure and
residual bit-for-bit. The continuous correction is `0.42624401466815054`
liquid-volume ULP, but the two adjacent representable volumes both change the
callback pressure and differ from the first-order prediction. The result is
therefore a binary64 representability finding, not a single-ULP plateau,
provider defect, or model/data decision. The propane result remains
solver-converged and numerically confirmed but not physically valid under the
current criterion. Predictive status remains
`NOT_ADJUDICATED_NO_APPROVED_HELD_OUT_CUTOFF`; neither the source uncertainty
nor this numerical failure is converted into a model-accuracy decision.

The diagnostic consumed provider wheel SHA-256
`f92f79c8d6f614660e5c201b7061c9b02b5cd1a25a4ed8c8fee0b59adaabf2bf`
and installed header SHA-256
`414c257d28322d6be41809a8dc6b98023859dd156202a5205079d674b28b4070`.
At liquid volume `6.296633932725173e-5 m3`, the callback returns
`2.9409473148611442 Pa` and
`dP/dV = -1.824582692583465e13 Pa/m3`. The immediately lower and higher
volumes return `2.940948085635517 Pa` and `2.940946595244363 Pa`; their
observed changes are `+7.707743727713989e-7 Pa` and
`-7.196167812750787e-7 Pa`, versus first-order changes
`+/-2.472770648973055e-7 Pa`. No tracked diagnostic program or replacement
receipt is retained.

The result reports full and parameter-column singular values, ranks, and
condition numbers. Those local diagnostics provide no global identifiability
or uncertainty claim. No source-row error has an approved pass threshold.

## Pressure-resolution correction protocol

This section is a design/evidence protocol. It changes no runtime, result,
receipt, artifact, gate, or authority. The current accepted criterion remains

```text
abs(P_phase - P_report) <= 1e-8 * abs(P_observed).
```

Historical review verdict `CORRECTION_DESIGN_JUSTIFIED` permits one candidate
correction to be evaluated:

```text
abs(P_phase - P_report)
    <= atol_resolution + 1e-8 * abs(P_observed).
```

The relative term is exactly `1e-8`. `P_observed` is only a magnitude scale;
neither it nor any reported experimental uncertainty determines
`atol_resolution`. No numeric absolute term is selected or approved here.

### Candidate derivations

Three derivations are considered:

1. A fixed number of ULPs of the returned pressure is rejected. At propane
   120 K the adjacent-volume callback changes span more than one billion ULPs
   of the approximately `3 Pa` returned pressure, so output spacing alone does
   not measure callback or input representability.
2. The largest observed `nextafter` pressure jump is retained as a diagnostic
   but rejected as the tolerance derivation. It combines real `dP/dV`
   sensitivity with binary64 evaluation error and would make the tolerance
   depend on the chosen state neighborhood.
3. The recommended minimum design is a certified high-precision discrepancy plus
   input-lattice bound. It uses an independently transcribed, non-installed
   oracle and exact installed callbacks. It does not use the propane 120 K
   closure miss or any experimental uncertainty.

For phase `a` at calibration state `s`, let `P64` be the installed callback,
`P*` and `(dP/dV)*` the independently evaluated high-precision values at the
exact binary64 inputs, and `V-`, `V`, and `V+` the adjacent representable
volume triplet. Define

```text
E_callback(s,a) = max_q abs(P64(q) - P*(q)), q in {V-, V, V+}
half_step_V     = 0.5 * max(V - V-, V+ - V)
E_lattice(s,a)  = half_step_V
                  * max_q abs((dP/dV)*(q)), q in {V-, V, V+}
E_transform(s,a)= certified propagated pressure error from binary64
                  V_reference * exp(u_phase)
E_report(s)     = certified absolute error plus half-ULP lattice bound for
                  binary64 P_observed * exp(u_pressure)
E_subtract(s,a) = 0.5 * ulp(P64(V) - P_report64)
B(s,a)          = E_callback + E_lattice + E_transform
                  + E_report + E_subtract.
```

The oracle must use directed rounding or a containing interval so each term is
an upper bound. The single proposed `atol_resolution` is the smallest
binary64 value not less than `max B(s,a)` over the predeclared calibration
cohort. The campaign must verify round-to-nearest-even and record the CPU,
libc/libm, compiler, and installed-wheel identities. No multiplier or
fit-to-pass adjustment is allowed. The candidate applies only to those exact
artifacts; a later provider or Regression artifact requires a fresh campaign
and approval.

The calibration cohort contains every immutable methane and ethane reporting
state and every immutable propane state except
`glos2004-propane-sat-120-k`, using primary reporting states and the existing
confirmation outputs. The excluded 120 K row is the locked pure challenge and
cannot enlarge the candidate. The binary extension in
`docs/science/neutral-hydrocarbon-next-slice.md` may test the same frozen value
but cannot participate in its derivation or widen it. If the 120 K certified
bound exceeds the candidate, or if the candidate does not separate the
declared negative controls, the protocol returns blocked and no replacement
value is chosen in the same campaign.

The current Provider suite has an independent 70-digit active-`kij` tensor
point but no pure-parameter pressure-resolution oracle over this finite state
matrix. Provider must add the smallest test-only independent pure
`(n,V,m,sigma,epsilon/k)` pressure and `dP/dV` oracle without adding a runtime
entry. Validation then owns the installed-artifact replay and the proposed
numeric receipt. Regression cannot select the number from its blocked propane
receipt.

### Replay and status invariants

For each methane, ethane, and propane training, held-out, and stress row, the
campaign records the primary reporting state, plus the existing confirmation
termination, usability, scaled-parameter delta, and cost delta. Each primary
reporting state has seven pressure probes: center; `nextafter` minus/plus one
ULP independently for liquid volume, vapor volume, and common reporting
pressure. The campaign records raw and scaled phase pressure and
chemical-potential closures, topology, stability, volume and pressure bounds,
callback fingerprints, oracle intervals, and every term in `B`. Axis probes
are sufficient; no Cartesian sweep or runtime ULP-probing surface is admitted.

The low-pressure/boundary sentinels are ethane 100 K, propane 110 K, propane
120 K, and one-ULP inside/outside probes for each declared liquid-volume,
vapor-volume, and reporting-pressure bound. Outside probes must reject before
scientific interpretation. Negative pressure-closure controls perturb the
common pressure in both directions by

```text
max(100 * atol_resolution, 1e-6 * abs(P_observed))
```

and must fail the mixed criterion. Wrong source, row order, component,
provider fingerprint, wheel, header, nonfinite input, topology, and stability
controls retain their existing rejection meanings.

Methane and ethane parameters, costs, ranks, reporting values, accepted
receipts, and authority decisions must remain numerically unchanged. Ethane
100 K and propane 110 K remain failed for their existing boundary and closure
reasons. Propane 120 K may change status only after a candidate value, its full
receipt, and the exact corrected artifact receive independent review.

Solver status is unchanged. A reporting row is marked
`pressure_resolution_limited = true` only when the original relative
criterion fails, the mixed criterion passes, and the independent oracle
certifies its bound. This one immutable row boolean is the minimum future
diagnostic addition because the existing failure reasons cannot truthfully
encode a successful but representability-limited result. No new result type or
generic tolerance object is allowed. The existing specification may gain only
one approved `reporting_pressure_resolution_atol_pa` scalar.

Physical validity can be true only when both liquid and vapor pressure
equalities pass the mixed criterion and chemical-potential closure, topology,
phase ordering, stability, bounds, source identity, and provider identity all
pass their existing gates. A representability-limited numerical pass does not
waive another gate. Predictive status remains
`NOT_ADJUDICATED_NO_APPROVED_HELD_OUT_CUTOFF`.

## Failure boundaries

The package rejects source, component, unit, row-order, partition,
specification, and provider-fingerprint mismatches before Ceres starts. Native
callbacks retain compact provider, topology, derivative, solve, and reporting
failures.

The compiled problem validates the full source and specification identity and
returns it to Python. Runtime code claims no provider commit or wheel digest.
The candidate runner binds and byte-checks the installed wheels before import.

The private native test surface exposes the training residual and Jacobian
block. It does not expose the reporting Jacobian, and this slice adds no test-
only runtime seam.

The package owns no binary interaction, association, electrolyte, reactive,
generic-family, persistence, or provider-catalog capability.
