# Downstream consumer preparation review — 2026-09-08

Status: **preparation inventory complete; Case 3C and the implemented Lithium
producer route are selected bounded routes; no Engine or application numerical
acceptance claimed**. This memo is the root review of the three independent
draft consumer specifications. It was prepared from the current draft clones
and retained provenance, not from guessed labels in historical reports.

## Evidence boundary

The consumer calls, input assumptions, output quantities and derivative needs
below were inspected in the independent clones. Their current migration
specifications own execution details; the parent implementation plan owns the
execution queue. This inventory supports design, not numerical qualification.

| application | actual current consumer route | inputs, units, species order and phases | outputs and derivative actions | failure behavior | Engine owner and concrete gap |
|---|---|---|---|---|---|
| MEA-Thermodynamics | analyses/mea_parameter_bundle/scripts/shared_evaluation.py drives 123 state rows. VLE rows solve pressure; speciation rows fix T/P; thermal rows use EquilibriumEnthalpy. | T (K), P (Pa), finite amounts (mol), molar masses (kg/mol), balance totals, reaction standard-state and parameter fingerprints. Nine liquid species in fixed order: CO2, MEA, H2O, MEAH+, MEACOO-, HCO3-, CO3--, H3O+, OH-. VLE adds an incipient neutral vapor with CO2/MEA/H2O only. | 79 system pressures and vapor outputs (Pa/dimensionless), 106 true-species liquid targets and 25 linear aggregate targets (dimensionless), phase density (mol/m3), amounts (mol), packing/composition and thermochemistry (J, J/mol). Ordinary rows request values. Anchor/recovery and fit calls consume state-input blocks and exact active-parameter Jacobian rows in caller order. | Bounded child timeout, exception, non-evaluable root and failed numerical/physical status are retained with code, diagnostic and recovery attempts. No zero/default output or silent row drop. | Engine owns model/state, phase/root/stability, fugacity/activity, enthalpy and native implicit actions. App owns source conversion, aggregate basis, reaction-coordinate values, thermochemical references and objective weighting. Gap: fresh replay and independent qualification on one installed wheel; heat Jacobian and broader pH/ionic/dielectric consumers remain unqualified. |
| MEA-Absorption-Column | BVP/ABS_Column.py consumes ReactiveLiquid.solve and neutral_vapor_state at column nodes; analyses/nccc_validation/scripts/run_reactive_column.py provides Case 3C; BVP/reactive_jacobian.py assembles A2. | Physical BVP state is [Fl_CO2, Fl_H2O, Fv_CO2, Fv_H2O, Hlf, Hvf, P] with flows mol/s, temperatures K or enthalpy flows J/s by mode, pressure Pa. Apparent liquid is [CO2, MEA, H2O], apparent vapor [CO2, H2O, N2, O2]. Reactive liquid uses the same nine-species order as MEA; neutral vapor support is CO2/MEA/H2O; two phase roles are finite liquid and finite neutral vapor. | RHS flow/energy/pressure derivatives (mol/(s m), J/(s m), Pa/m), fugacity and density (Pa, mol/m3), transport/flux/enthalpy outputs, plus 7x7 scaled A2 state Jacobian. A2 native vector is 9 amounts (mol), liquid density (mol/m3), liquid fCO2 (Pa), neutral-vapor fCO2 (Pa). Native state-input derivatives are mapped to T/P/apparent-input coordinates; current A2 requires raw temperature, bidirectional CO2 and no gas-velocity exponent. | Invalid trial state is raised or captured by guard diagnostics; failed Engine action/derivative, density mismatch, unsupported mode, ideal-Henry initialization timeout or BVP nonconvergence remains a failed run with message and partial diagnostics. Reactive branch has no silent ideal fallback. | Engine owns liquid/vapor state, fugacity, phase certificate, enthalpy and exact actions. Absorber owns geometry, hydraulic/transport/heat-transfer correlations, flux/enhancement closure, BVP and scaling. Gap: fresh wheel-qualified A1/A2 and full Case 3C replay; empirical-vs-true density basis and derivative coverage outside the supported mode need explicit disposition. 6c99 cannot be accepted while absent. |
| Lithium_Extraction | Current MS-2 IDAES stage model consumes the retained loaded-organic and reactive-LLE application surrogates. The direct Engine producer is the calibrated generator's `independent_network.solve_condition`; the configured oracle is a separate manufactured capability check. | Stage inputs: Li/Na (mol/L), pH, carrier (mol/L), O/A, T (K), aqueous flow (L/h), stage count. Loaded-state contract bounds include O/A 0.5--1.5, T 293.15--313.15 K and labeled Li/Na/H/carrier mol/L inputs. Direct Engine conditions use T (K), fixed P 101325 Pa, aqueous/organic volumes (L), Li/Na/H/Cl and carrier amounts (mol/L), saponification and incoming loadings. Independent-network order is water, diluent, Li+, Na+, H+, Cl-, free-carrier, sodium-carrier, lithium-carrier with aqueous and organic declared supports; the seven-species manufactured oracle is a separate order water, diluent, Li+, Na+, Cl-, sodium-carrier, lithium-carrier. Both phases are finite liquids. | Stage model outputs Li/Na extraction/loading/distribution/separation fractions and max balance error (dimensionless or mol/h). Direct Engine outputs phase amounts (mol), molar density (mol/m3), extraction, D/selectivity, phase volume and balance/affinity/KKT/Hessian/absence diagnostics. It consumes an exact active pair/diluent/lithium-carrier/k_ij Jacobian and derives Li-extraction sensitivity. The manufactured oracle derivative centered check has retained scaled error 1.85e-10. | Surrogate schema/domain mismatch or out-of-domain input raises; IDAES DOF/nonoptimal solve and fixed-point nonconvergence fail. Engine solver/domain/certificate failures retain each initialization attempt; only explicitly bounded application-accepted candidates can be reported as such. Rezaee 2025/2026 paths are retained negative endpoints, not forced fits. | Engine owns generic two-liquid finite phases, declared supports, amounts/density and exact parameter actions. Lithium owns HBTA/TOPO chemistry, reaction constants, diluent mapping, source interpretation, stage connectivity, surrogates and process limits. Current runs exercise the producer and its surrogate using the selected inputs and installed wheel; old run identities do not block implementation. Broader pure/pair/tie-line, single-ion and complex inputs remain claim limits, not automatic blockers for this route. |

### Required future consumers

The inventory includes consumers that are not all active in the current
acceptance path. They must stay explicit rather than being inferred from a
successful value solve:

- MEA needs native heat-quantity derivatives (or an owner-approved deferral),
  pH/ionic activity, dielectric/permittivity and ionic-volumetric outputs with
  source standard states and domains.
- The absorber needs A2 actions for enthalpy-state, unidirectional/absorption-only,
  gas-velocity and alternate vapor-convention modes. Each needs its own output
  order, units, chain rule and falsifier; selecting Case 3C does not silently
  defer these required gaps.
- Lithium's direct producer needs the existing fixed topology, dodecane
  assumption, feed/condition bounds and application fit inputs to produce a
  bounded result. Its current source feed has blank pH, chloride and TOC, so
  pH/H and charge-closed chloride are explicit modeling assumptions. Source
  qualified pure/pair/tie-line, single-ion and complex data are required only
  for broader physical identification or predictive transfer. The retained
  MS-2 response artifacts remain bounded to their stated feed, fresh dodecane,
  298.15 K, one-stage and O/A domain.

## Use and limits

Current application choices and execution gates are owned by the
[implementation plan](../implementation-plan.md) and each clone's migration
specification. This inventory supplies consumer evidence, not a second decision
queue. Case 3C is selected; Lithium uses its implemented producer and adopted
assumptions. The seven source-readiness rows constrain broader identification
claims, not the current bounded calculation.

No numerical replay was run for this inventory. Qualification requires actual
values and failures from the selected inputs and installed wheel, with justified
numerical limits and independent review. A build or readiness check alone does
not establish agreement.
