# General Regression quickstart

The public direct-EOS path is `prepare_fit -> preflight -> fit -> export`.
Preparation resolves capability IDs, EOS units, parameter/topology
fingerprints, and transformed-row hashes from the installed model. The caller
still supplies every scientific choice: source, rows, objective, bounds,
affine scales, complete starts, partitions, and lifted-volume contracts.

```python
from epcsaft import Mixture, Parameters
from epcsaft.regression import (
    AcquisitionClass, AffineParameterTransform, ConfirmationControls,
    FixedCompositionVleObservation, ObjectiveContract, ObservationDataset,
    PairParameterIdentity, ParameterFamily, ParameterRequest, RankControls,
    ResultContext, RowProvenance, SolverControls, SourceInput,
    parameter_capabilities, prepare_fit,
)

model = Mixture(Parameters.from_bundle(
    "/materialized/gross-2001-methane-ethane/1/parameters",
    components=("methane", "ethane"),
))
assert any(
    capability.installed_ready
    and getattr(capability, "family", None) is ParameterFamily.K_IJ
    for capability in parameter_capabilities(model)
)
records = ({
    "row_id": "may2015-ch4-c2h6-002",
    "source_id": "may-2015",
    "source_locator": "May et al. (2015), retained row 002",
    "component_ids": ["methane", "ethane"],
    "temperature_k": 203.22,
    "pressure_pa": 2_124_000.0,
    "liquid_mole_fraction_first": 0.3653,
    "vapor_mole_fraction_first": 0.8667,
    "pressure_scale_pa": 2_124_000.0,
    "chemical_potential_scales": [1.0, 1.0],
    "liquid_volume_origin_m3_per_mol": 6.0e-5,
    "liquid_volume_start_m3_per_mol": 6.5e-5,
    "liquid_volume_bounds_m3_per_mol": [2.0e-5, 1.0e-4],
    "vapor_volume_origin_m3_per_mol": 9.0e-4,
    "vapor_volume_start_m3_per_mol": 1.0e-3,
    "vapor_volume_bounds_m3_per_mol": [1.0e-4, 1.0e-2],
    "partition": "training",
},)
dataset = ObservationDataset.from_records(
    FixedCompositionVleObservation, records,
    source=SourceInput(
        "may-2015", "May et al. (2015), methane/ethane VLE",
        "retained source artifact",
        "5cd1e74925a3c6504f5106dcf911f2cae2d6e99a5133fccc20454d8991bdbc7f",
        "No transformation for this retained row.",
        "T/K, P/Pa, mole fractions",
        "Regression derivative-contract evidence",
        "pressure by observed P; mu/RT dimensionless",
    ),
    objective=ObjectiveContract(
        "fixed_composition_vle", "native_scaled_least_squares",
        "observation_residual_scales",
        "independent_no_covariance", "squared", (), "fail_fit",
    ),
    row_provenance={"may2015-ch4-c2h6-002": RowProvenance(
        AcquisitionClass.DIRECT_MEASUREMENT, "unique source row", "included",
        "outside critical-region exclusion", "not censored",
        "retained; no outlier rule applied",
    )},
)
prepared = prepare_fit(
    model, datasets=(dataset,),
    parameters=(ParameterRequest(
        ParameterFamily.K_IJ, PairParameterIdentity("methane", "ethane"),
        AffineParameterTransform(0.0, 0.01), -0.15, 0.10,
    ),),
    parameter_slot_indices=(0,),
    start_vectors=((0.0,), (-0.05,), (0.05,)),
    solver=SolverControls(50, 30.0, 1e-12, 1e-12, 1e-12),
    rank=RankControls(1e10),
    confirmation=ConfirmationControls(1e-5, 1e-8),
)
report = prepared.preflight()  # exact evaluations; no Ceres solve
if not report.ready:
    raise RuntimeError(report.reasons)
evaluation = prepared.evaluate((0.0,))  # exact read-only assembly; no Ceres
assert evaluation.jacobian_layout == "row_major"
assert len(evaluation.jacobian) == (
    evaluation.jacobian_diagnostics.residual_count
    * evaluation.jacobian_diagnostics.variable_count
)
result = prepared.fit()
canonical_json = result.to_json_bytes(prepared=prepared, context=ResultContext())
assert canonical_json == result.to_json_bytes(
    prepared=prepared, context=ResultContext()
)
record = result.to_record(prepared=prepared)
assert record["status"]["solver_converged"] == result.solver_converged
assert record["row_accounting"]["failed"] == result.failed_row_count
```

Preflight reports `N`, `Q`, `R`, partitions, derivative completeness, full and
nuisance-projected local rank/conditioning, start-bound status, and failure
reasons. It does not report practical/global identifiability, uncertainty,
prediction, or acceptance. Repeated or insensitive row designs therefore fail
rank preflight before Ceres.

`PreparedFit.evaluate(physical_parameter_point, lifted_solver_point=None)` is
the public audit seam for one declared training-residual point. When the
lifted point is omitted, multi-parameter pure-phase contracts use the EOS
root-resolved phase start; other contracts use their row-declared starts. The
immutable result records the exact selected point. It contains the physical
fitted point, transformed fitted and lifted solver points, lifted physical
values, ordered fitted coordinates and lifted-variable IDs,
ordered row/component residual IDs, residual vector, complete exact row-major
Jacobian, and full and nuisance-projected diagnostics. It also binds parameter,
topology, capability-artifact, installed-EOS-artifact, preparation, dataset,
source, and observation-order fingerprints. Passing the returned
`lifted_solver_point` back explicitly permits independent perturbation of every
fitted and nuisance column. Evaluation never enters Ceres or changes authority;
domain failures, incomplete derivatives, contract mismatches, and nonfinite
values fail closed with stable classified prefixes.

`parameter_capabilities(model)` is the machine-readable support view. It
returns the installed EOS capability records directly,
including units, identity and observation contracts, derivative order,
`installed_ready`, and `unsupported_reason`. The human-readable family matrix
is in
[`general-parameter-regression.md`](science/general-parameter-regression.md#definitive-parameter-family-and-data-requirement-matrix).

## Joint pure and rank-deficient examples

This installed-artifact example uses six explicitly copied NIST WebBook rows.
It does not call a paper-specific Regression helper. The repeated-state variant
then demonstrates a projected-rank failure without entering Ceres.

```python
from epcsaft import Mixture, Parameters
from epcsaft.regression import (
    AcquisitionClass, AffineParameterTransform, ComponentParameterIdentity,
    ConfirmationControls, ObjectiveContract, ObservationDataset,
    ParameterFamily, ParameterRequest, PureSaturationObservation,
    RankControls, RowProvenance, SolverControls, SourceInput, prepare_fit,
)

SATURATION = (
    (100.0, 34375.892, 438.88524),
    (110.0, 88130.038, 424.77725),
    (120.0, 191430.08, 409.90234),
    (130.0, 367319.94, 394.03734),
    (140.0, 641181.43, 376.86505),
    (150.0, 1039961.3, 357.89846),
)

def saturation_record(index, temperature, pressure, density):
    liquid_volume = 0.016043 / density
    vapor_volume = 8.31446261815324 * temperature / pressure
    return {
        "row_id": f"nist-methane-{index}",
        "source_id": "nist-webbook-methane-saturation",
        "source_locator": f"installed copied grid row {index}",
        "component_id": "methane",
        "temperature_k": temperature,
        "pressure_pa": pressure,
        "liquid_density_kg_per_m3": density,
        "molar_mass_kg_per_mol": 0.016043,
        "pressure_scale_pa": 2.0 * pressure,
        "chemical_potential_scale": 2.0,
        "liquid_density_scale_kg_per_m3": 2.0 * density,
        "liquid_volume_origin_m3_per_mol": liquid_volume,
        "liquid_volume_start_m3_per_mol": liquid_volume,
        "liquid_volume_bounds_m3_per_mol": (2.0e-5, 1.0e-4),
        "vapor_volume_origin_m3_per_mol": vapor_volume,
        "vapor_volume_start_m3_per_mol": vapor_volume,
        "vapor_volume_bounds_m3_per_mol": (1.5e-4, 0.1),
        "partition": "training",
    }

rows = tuple(
    saturation_record(index, *values)
    for index, values in enumerate(SATURATION, start=1)
)
source = SourceInput(
    "nist-webbook-methane-saturation",
    "NIST Chemistry WebBook SRD 69, methane saturation table",
    "https://webbook.nist.gov/cgi/fluid.cgi?ID=C74828",
    "dec64d5a6cac414a4a92393a0d728fa27c02135c6a159d0d1881d7b6dde6d26c",
    "Copied the declared 10 K grid; SI values retained.",
    "T/K, P/Pa, saturated-liquid density/(kg/m3)",
    "Public database rows; routine implementation example only.",
    "Equal four-residual mechanics: 2P, 2 for mu/RT, and 2rho.",
)
objective = ObjectiveContract(
    "pure_saturation", "native_scaled_least_squares",
    "observation_residual_scales", "independent_no_covariance",
    "squared", (), "fail_fit",
)

def provenance(records):
    return {
        row["row_id"]: RowProvenance(
            AcquisitionClass.DATABASE_RECORD, "unique copied row", "included",
            "below the declared critical-region cutoff", "not censored",
            "retained; no outlier rule applied",
        )
        for row in records
    }

dataset = ObservationDataset.from_records(
    PureSaturationObservation, rows, source=source, objective=objective,
    row_provenance=provenance(rows),
)
model = Mixture(Parameters.from_bundle(
    "/materialized/gross-2001-methane-ethane/1/parameters",
    components=("methane",),
))
requests = (
    ParameterRequest(
        ParameterFamily.SEGMENT_COUNT, ComponentParameterIdentity("methane"),
        AffineParameterTransform(1.0, 0.1), 0.5, 3.5,
    ),
    ParameterRequest(
        ParameterFamily.SEGMENT_DIAMETER,
        ComponentParameterIdentity("methane"),
        AffineParameterTransform(3.7039, 0.1), 2.0, 5.0,
    ),
    ParameterRequest(
        ParameterFamily.DISPERSION_ENERGY_OVER_K,
        ComponentParameterIdentity("methane"),
        AffineParameterTransform(150.03, 10.0), 50.0, 400.0,
    ),
)
controls = dict(
    solver=SolverControls(200, 30.0, 1e-10, 1e-10, 1e-10),
    rank=RankControls(1e12),
    confirmation=ConfirmationControls(1e-5, 1e-8),
)
joint = prepare_fit(
    model, datasets=(dataset,), parameters=requests,
    parameter_slot_indices=(0, 1, 2),
    start_vectors=((1.08, 3.555744, 157.5315), (1.0, 3.7, 150.0)),
    **controls,
)
assert joint.preflight().ready
assert len(joint.fit().parameters) == 3

repeated_rows = tuple(
    {**rows[0], "row_id": f"repeated-state-{index}"}
    for index in range(3)
)
rank_dataset = ObservationDataset.from_records(
    PureSaturationObservation, repeated_rows, source=source,
    objective=objective, row_provenance=provenance(repeated_rows),
)
rank_deficient = prepare_fit(
    model, datasets=(rank_dataset,), parameters=requests,
    parameter_slot_indices=(0, 1, 2),
    start_vectors=((1.08, 3.555744, 157.5315), (1.0, 3.7, 150.0)),
    **controls,
).preflight()
assert not rank_deficient.ready
assert any("rank_deficient" in reason for reason in rank_deficient.reasons)
```

## Fixed-topology association and constrained-subset example

The installed ethanol model happens to fix a 2B topology (one A and one B
site), but Regression consumes it through the same descriptor used for any
supported finite topology. This example does not infer a topology. The rows are a declared synthetic
installed-model mechanics fixture, so they are not literature or validation
evidence. The full block reports issue #28 as an unresolved practical-
identifiability gate. The second preparation deliberately fits only
`epsilon_AB/k`; `kappa_AB` remains fixed at the source value for this example,
which is a case-specific constraint rather than a universal strategy.

```python
from epcsaft import Mixture, Parameters, unit_registry
from epcsaft.regression import (
    AcquisitionClass, AffineParameterTransform,
    ConfirmationControls, CorrelationProvenance, FixedTopologyAssociationCapability,
    ObjectiveContract, ObservationDataset, ParameterFamily, ParameterRequest,
    PureDensityObservation, RankControls, RowProvenance, SolverControls,
    SourceInput, parameter_capabilities, prepare_fit,
)

model = Mixture(Parameters.from_bundle(
    "/materialized/figiel-2025-reference-electrolytes/1/parameters",
    components=("ethanol",),
))
volumes = (4.0e-5, 4.5e-5, 5.0e-5, 5.2e-5, 5.5e-5)
records = []
for index, volume in enumerate(volumes, start=1):
    state = model.state(
        T=298.15 * unit_registry.kelvin,
        rho=(1.0 / volume) * unit_registry.mole / unit_registry.meter**3,
        x=(1.0,),
    )
    pressure = float(state.pressure.to("pascal").magnitude)
    records.append({
        "row_id": f"ethanol-mechanics-{index}",
        "source_id": "installed-epcsaft-wheel-ethanol-mechanics-fixture",
        "source_locator": f"generated volume {volume}",
        "component_id": "ethanol",
        "temperature_k": 298.15,
        "pressure_pa": pressure,
        "density_kg_per_m3": 0.046069 / volume,
        "molar_mass_kg_per_mol": 0.046069,
        "pressure_scale_pa": pressure,
        "density_scale_kg_per_m3": 0.046069 / volume,
        "volume_origin_m3_per_mol": volume,
        "volume_start_m3_per_mol": volume,
        "volume_bounds_m3_per_mol": (3.5e-5, 6.0e-5),
        "partition": "training",
    })
records = tuple(records)
source = SourceInput(
    "installed-epcsaft-wheel-ethanol-mechanics-fixture",
    "Synthetic values generated from the installed Figiel ethanol model",
    "epcsaft==0.2.0.dev0 wheel artifact used by this repository",
    "bc7e637de084330ebded4ddfd52e02bc1ce5451221128692972ebba8856d098e",
    "Evaluated fixed-volume public EOS states without row selection.",
    "T/K, P/Pa, density/(kg/m3)",
    "Implementation mechanics only; no scientific or predictive authority.",
    "Pressure and density divided by their finite generated magnitudes.",
)
objective = ObjectiveContract(
    "pure_density", "native_scaled_least_squares",
    "observation_residual_scales", "independent_no_covariance",
    "squared", (), "fail_fit",
)
fixture_correlation = CorrelationProvenance(
    "installed EOS fixed-volume evaluation",
    (("declared_fixture", 1.0),),
    "Pa and kg/m3",
    "exactly the five declared states",
    ("temperature_k",),
    tuple((298.15,) for _ in records),
    "No transformation beyond the public EOS state evaluation.",
)
dataset = ObservationDataset.from_records(
    PureDensityObservation, records, source=source, objective=objective,
    row_provenance={
        row["row_id"]: RowProvenance(
            AcquisitionClass.RECONSTRUCTED_CORRELATION,
            "unique generated state", "included", "not a critical-region claim",
            "not censored", "retained without outlier selection",
            correlation=fixture_correlation,
        )
        for row in records
    },
)
descriptor = next(
    capability for capability in parameter_capabilities(model)
    if isinstance(capability, FixedTopologyAssociationCapability)
)
origins = (2.3827, 3.1771, 198.24, 2653.4, 0.03238)
scales = (0.5, 0.2, 50.0, 500.0, 0.02)
bounds = (
    (1.0, 8.0), (2.0, 6.0), (50.0, 800.0),
    (100.0, 8000.0), (1.0e-5, 0.2),
)
requests = tuple(
    ParameterRequest(
        slot.family, slot.identity, AffineParameterTransform(origin, scale), *bound,
    )
    for slot, origin, scale, bound in zip(
        descriptor.slots, origins, scales, bounds, strict=True,
    )
)
controls = dict(
    solver=SolverControls(100, 30.0, 1e-10, 1e-10, 1e-10),
    rank=RankControls(1e12),
    confirmation=ConfirmationControls(1e-5, 1e-8),
)
full_block = prepare_fit(
    model, datasets=(dataset,), parameters=requests,
    parameter_slot_indices=(0, 1, 2, 3, 4),
    start_vectors=(origins, (2.4, 3.18, 200.0, 2600.0, 0.035)),
    **controls,
)
full_report = full_block.preflight()
assert full_report.practical_identifiability_gate is not None

energy_only = prepare_fit(
    model, datasets=(dataset,), parameters=(requests[3],),
    parameter_slot_indices=(0,),
    start_vectors=((2653.4,), (2600.0,)),
    **controls,
)
assert energy_only.preflight().ready
```

## Direct-observable example

This compact model-level dielectric fit uses three explicitly digitized rows.
It has no lifted phase variables.

```python
from epcsaft import Mixture, Parameters
from epcsaft.regression import (
    AcquisitionClass, AffineParameterTransform, ConfirmationControls,
    ModelParameterIdentity, ObjectiveContract, ObservationDataset,
    ParameterFamily, ParameterRequest, RankControls,
    RelativePermittivityRatioObservation, RowProvenance, SolverControls,
    SourceInput, prepare_fit,
)

source_rows = (
    ("figiel-dielectric-008", 0.010880829, 71.61417323),
    ("figiel-dielectric-009", 0.022927461, 66.65354331),
    ("figiel-dielectric-010", 0.033808290, 62.87401575),
)
records = tuple({
    "row_id": row_id,
    "source_id": "figiel-figure-2-dielectric",
    "source_locator": f"Figure 2 digitized row {row_id}",
    "solvent_id": "water",
    "component_ids": ("water", "sodium-cation", "bromide-anion"),
    "temperature_k": 298.15,
    "pressure_pa": 100_000.0,
    "total_ion_mole_fraction": ion_fraction,
    "observed_relative_permittivity_ratio": permittivity / 78.09,
    "residual_scale": 1.0,
    "partition": "training",
} for row_id, ion_fraction, permittivity in source_rows)
source = SourceInput(
    "figiel-figure-2-dielectric",
    "Figiel, Yu, and Held (2025), Figure 2",
    "declared digitization ledger",
    "09e1e820a55861b835fdab27df5134451ecc9329c6d512dcf26565b267b387a6",
    "Digitized three NaBr rows and divided epsilon_r by 78.09.",
    "T/K, P/Pa, mole fraction, dimensionless ratio",
    "Digitized in-sample implementation evidence only.",
    "Equal raw dimensionless residuals; no pointwise uncertainty.",
)
dataset = ObservationDataset.from_records(
    RelativePermittivityRatioObservation, records, source=source,
    objective=ObjectiveContract(
        "relative_permittivity_ratio", "native_scaled_least_squares",
        "observation_residual_scales", "independent_no_covariance",
        "squared", (), "fail_fit",
    ),
    row_provenance={
        row["row_id"]: RowProvenance(
            AcquisitionClass.DIGITIZED_VALUE, "unique digitized row",
            "included", "not a critical-region observable", "not censored",
            "retained without outlier selection",
        )
        for row in records
    },
)
model = Mixture(Parameters.from_bundle(
    "/materialized/figiel-2025-reference-electrolytes/1/parameters",
    components=("water", "sodium-cation", "bromide-anion"),
))
prepared = prepare_fit(
    model, datasets=(dataset,),
    parameters=(ParameterRequest(
        ParameterFamily.ION_FRACTION_SUPPRESSION_COEFFICIENT,
        ModelParameterIdentity(), AffineParameterTransform(7.0, 1.0),
        0.01, 30.0,
    ),),
    parameter_slot_indices=(0,), start_vectors=((2.0,), (12.0,)),
    solver=SolverControls(50, 30.0, 1e-12, 1e-12, 1e-12),
    rank=RankControls(1e10),
    confirmation=ConfirmationControls(1e-5, 1e-8),
)
assert prepared.preflight().lifted_variable_count == 0
assert prepared.fit().workflow_valid
```

## Parameter blocks and association

- Scalar fits use one request and one slot.
- Joint pure `(m, sigma, epsilon/k)` uses those families and slots in that
  order, complete three-value starts, and pure-saturation rows.
- Fixed-topology association uses any ordered subset of the installed EOS
  descriptor's ordinary-component and component/site-pair slots. Multiple EOS
  slots may explicitly share one fitted coordinate. Pure vapor-pressure,
  fixed-pressure-density, and combined saturation rows use the pure callback;
  associating-mixture density rows use the matching mixture callback.
  Topology is immutable EOS input. Issue #28 owns nuisance-reoptimized profile
  evidence for the retained 2B literature case. Any constrained subset must
  declare its source and validity domain; it is not universal.
- Fixed-composition liquid-mixture-density rows use the same lifted Ceres
  owner for an ordered caller-supplied component set, mole-fraction vector,
  mixture molar-mass basis, pressure/density scales, and liquid-volume
  origins/starts/bounds. No chemistry, missing-parameter defaults, or
  association topology is inferred; the installed EOS must advertise the
  matching phase capability and component order.
- Direct-observable families use their MIAC, solvation-Gibbs, or relative-
  permittivity-ratio rows and no lifted variables.

Induced and resolved cross association are not current pure-fit coordinates.
Issue #34 owns future fixed alternatives after EOS advertises their exact
identities, transforms, topology fingerprints, and derivatives.

## Record templates

`ObservationDataset.from_records` requires every listed row field and rejects
unknown fields; there are no optional row fields in these contracts. All
direct-EOS rows also require `row_id`, `source_id`, `source_locator`, and
`partition`.

| Observation type | Additional required row fields |
|---|---|
| `FixedCompositionVleObservation` | `component_ids`, `temperature_k`, `pressure_pa`, `liquid_mole_fraction_first`, `vapor_mole_fraction_first`, `pressure_scale_pa`, `chemical_potential_scales`, liquid/vapor `volume_origin_m3_per_mol`, `volume_start_m3_per_mol`, and `volume_bounds_m3_per_mol` |
| `PureSaturationObservation` | `component_id`, `temperature_k`, `pressure_pa`, `liquid_density_kg_per_m3`, `molar_mass_kg_per_mol`, `pressure_scale_pa`, `chemical_potential_scale`, `liquid_density_scale_kg_per_m3`, and both phase volume origins/starts/bounds |
| `PureVaporPressureObservation` | `component_id`, `temperature_k`, `pressure_pa`, `pressure_scale_pa`, `chemical_potential_scale`, and separated liquid/vapor volume origins/starts/bounds |
| `PureDensityObservation` | `component_id`, `temperature_k`, `pressure_pa`, `density_kg_per_m3`, `molar_mass_kg_per_mol`, `pressure_scale_pa`, `density_scale_kg_per_m3`, `volume_origin_m3_per_mol`, `volume_start_m3_per_mol`, `volume_bounds_m3_per_mol` |
| `FixedCompositionMixtureDensityObservation` | ordered `component_ids`, `mole_fractions`, `temperature_k`, `pressure_pa`, `density_kg_per_m3`, `molar_mass_kg_per_mol`, pressure/density scales, and liquid-volume origin/start/bounds |
| `MeanIonicActivityObservation` | ordered `component_ids`, `active_component_id`, `temperature_k`, `pressure_pa`, `formula_unit_molality_mol_per_kg`, `observed_mean_ionic_activity_coefficient`, `relative_residual_scale` |
| `AqueousKijMeanIonicActivityObservation` | ordered `component_ids`, `active_pair_component_ids`, complete `fixed_k_ij`, `temperature_k`, `pressure_pa`, `formula_unit_molality_mol_per_kg`, `observed_mean_ionic_activity_coefficient`, `relative_residual_scale` |
| `SolvationGibbsObservation` | ordered `component_ids`, `active_component_id`, `temperature_k`, `pressure_pa`, `observed_solvation_gibbs_j_per_mol`, `residual_scale_j_per_mol` |
| `IonSolvationKijObservation` | ordered `component_ids`, `active_component_id`, `active_pair_component_ids`, complete `fixed_k_ij`, `temperature_k`, `pressure_pa`, `observed_solvation_gibbs_j_per_mol`, `residual_scale_j_per_mol` |
| `RelativePermittivityRatioObservation` | `solvent_id`, ordered `component_ids`, `temperature_k`, `pressure_pa`, `total_ion_mole_fraction`, `observed_relative_permittivity_ratio`, `residual_scale` |

The specialized composed-positive transport constructs
`PositiveScalarObservation` directly. It requires `row_id`, `state_id`,
`state_schema_id`, `source_id`, `source_locator`, `primitive_id`,
`primitive_unit`, `transform`, `reference_id`, `reference_fingerprint`,
`observed_value`, `residual_scale`, `residual_scale_unit`, and `partition`.
It does not use `ObservationDataset` or the direct-EOS preparation surface.
Consequently it is not eligible for literature-result export: that export
requires the complete `PreparedFit` acquisition, row-decision, and objective
provenance contract. Composed-positive rows remain source-bound evaluator
transport, not a substitute literature-data schema.

Every dataset requires all `SourceInput` fields: `source_id`, `citation`,
`durable_locator`, `source_artifact_sha256`, `transformation_record`,
`units_and_bases`, `use_basis`, and `residual_scale_rationale`. Every row then
requires a `RowProvenance` with `acquisition`, `duplicate_decision`,
`exclusion_decision`, `critical_region_decision`, `censoring_decision`, and
`outlier_decision`. Author- and reconstructed-correlation rows additionally
require `CorrelationProvenance(equation, coefficients, units,
validity_interval, sampling_fields, sampling_grid, transformation_record)`;
its ordered, potentially multidimensional grid must exactly match the named
numeric fields on its transformed rows.

The current executable `ObjectiveContract` is deliberately finite:

- `residual_family` must exactly match the observation contract;
- `interpretation="native_scaled_least_squares"`;
- `row_weighting="observation_residual_scales"`, binding the explicit scales
  carried by each validated row;
- `covariance_interpretation="independent_no_covariance"`;
- `loss="squared"` with `loss_parameters=()`; and
- `failed_row_treatment="fail_fit"`.

An objective requiring covariance, a robust loss, censoring, aggregation, or
another unsupported meaning fails closed; it is not retained as ignored prose.
Dense correlation samples are not independent experiments.

Export returns data only. It writes no receipt, mutates no EOS catalog, and
grants no production, predictive, scientific, or authority status. The record
keeps those status axes separate. A literature export requires the exact
`PreparedFit` so row provenance and the caller's objective cannot be detached
from the result. `ResultContext` has explicit immutable-reference slots for
profile, bootstrap, uncertainty, and Validation-campaign artifacts; each is
either an `ArtifactReference(artifact_id, sha256)` or `None`. Failed diagnostic
records omit unavailable numerical fields and list their exact paths under
`unavailable_fields`; evaluated rows and all other fields reject nonfinite
values. Fixed-topology descriptor limits are metadata. An intentional open
upper limit omits the numeric `upper_bound_exclusive` field, lists that path
under `unavailable_fields`, and records
`{"field":"upper_bound_exclusive","direction":"upper"}` under `open_bounds`.
Lower limits remain finite. NaN or an infinity in any other position remains
invalid, and this representation does not relax finite-number checks for
evaluated scientific values.
