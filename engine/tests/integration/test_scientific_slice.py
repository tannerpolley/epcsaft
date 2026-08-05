from __future__ import annotations

import hashlib
import json
import math
from dataclasses import asdict
from pathlib import Path

import pytest

from epcsaft import Mixture, Parameters, SourceRecord, ValidityDomain, unit_registry
from epcsaft import equilibrium
from epcsaft.regression import (
    AcquisitionClass,
    AffineParameterTransform,
    ConfirmationControls,
    FixedCompositionMixtureDensityObservation,
    ObjectiveContract,
    ObservationDataset,
    ObservationPartition,
    PairParameterIdentity,
    ParameterFamily,
    ParameterRequest,
    RankControls,
    RowProvenance,
    SolverControls,
    SourceInput,
    prepare_fit,
)


_FIXTURE = Path(__file__).parents[1] / "fixtures" / "data" / "gross-2001-methane-ethane"
_PACKET = _FIXTURE / "1"


def _packet_fingerprint(packet: Path) -> str:
    digest = hashlib.sha256()
    for path in sorted(item for item in packet.rglob("*") if item.is_file()):
        relative = path.relative_to(packet).as_posix().encode()
        payload = path.read_bytes()
        digest.update(len(relative).to_bytes(8, "big"))
        digest.update(relative)
        digest.update(len(payload).to_bytes(8, "big"))
        digest.update(payload)
    return digest.hexdigest()


def _observations(model: Mixture) -> tuple[FixedCompositionMixtureDensityObservation, ...]:
    rows = []
    fractions = (0.35, 0.65)
    molar_mass = 0.35 * 0.0160428 + 0.65 * 0.030069
    for index, (temperature, volume) in enumerate(
        zip((220.0, 240.0, 260.0, 280.0), (4.0e-4, 5.0e-4, 6.0e-4, 7.0e-4), strict=True)
    ):
        pressure = float(
            model.state(
                T=temperature * unit_registry.kelvin,
                rho=(1.0 / volume) * unit_registry.mole / unit_registry.meter**3,
                x=fractions,
            )
            .pressure.to("pascal")
            .magnitude
        )
        rows.append(
            FixedCompositionMixtureDensityObservation(
                row_id=f"materialized-density-{index}",
                source_id="materialized-density",
                source_locator=f"manufactured-from-materialized-packet:{index}",
                component_ids=("methane", "ethane"),
                mole_fractions=fractions,
                temperature_k=temperature,
                pressure_pa=pressure,
                density_kg_per_m3=molar_mass / volume,
                molar_mass_kg_per_mol=molar_mass,
                pressure_scale_pa=pressure,
                density_scale_kg_per_m3=molar_mass / volume,
                volume_origin_m3_per_mol=volume,
                volume_start_m3_per_mol=volume,
                volume_bounds_m3_per_mol=(0.8 * volume, 1.2 * volume),
                partition=ObservationPartition.TRAINING,
            )
        )
    return tuple(rows)


def test_materialized_packet_fit_transfers_to_eos_and_equilibrium() -> None:
    materialization = json.loads((_FIXTURE / "materialization.json").read_text())
    assert materialization["data_commit"] == "9185bce8452191e6425e164c51dd6f8db983fc9b"
    assert _packet_fingerprint(_PACKET) == materialization["packet_fingerprint"]
    for expected in materialization["files"]:
        assert (
            hashlib.sha256((_PACKET / expected["path"]).read_bytes()).hexdigest()
            == (expected["sha256"])
        )

    parameters = Parameters.from_bundle(
        _PACKET / "parameters",
        components=("methane", "ethane"),
    )
    model = Mixture(parameters)
    rows = _observations(model)
    records = []
    provenance = {}
    for row in rows:
        record = asdict(row)
        record["partition"] = row.partition.value
        records.append(record)
        provenance[row.row_id] = RowProvenance(
            AcquisitionClass.DIRECT_MEASUREMENT,
            "deterministic materialized-packet state",
            "included",
            "outside critical-region screening",
            "not censored",
            "no outlier rule",
        )
    dataset = ObservationDataset.from_records(
        FixedCompositionMixtureDensityObservation,
        tuple(records),
        source=SourceInput(
            "materialized-density",
            "Deterministic density states generated from the pinned Data packet",
            materialization["packet_path"],
            materialization["packet_fingerprint"],
            "all four declared states",
            "SI; methane then ethane",
            "cross-module transfer slice",
            "pressure and mass density",
        ),
        objective=ObjectiveContract(
            "fixed_composition_mixture_density",
            "native_scaled_least_squares",
            "observation_residual_scales",
            "independent_no_covariance",
            "squared",
            (),
            "fail_fit",
        ),
        row_provenance=provenance,
    )
    prepared = prepare_fit(
        model,
        datasets=(dataset,),
        parameters=(
            ParameterRequest(
                ParameterFamily.K_IJ,
                PairParameterIdentity("methane", "ethane"),
                AffineParameterTransform(0.0, 0.05),
                -0.1,
                0.1,
            ),
        ),
        parameter_slot_indices=(0,),
        start_vectors=((-0.03,), (0.02,)),
        solver=SolverControls(150, 8.0, 1.0e-12, 1.0e-12, 1.0e-12),
        rank=RankControls(1.0e12),
        confirmation=ConfirmationControls(1.0e-6, 1.0e-8),
    )
    assert prepared.preflight().ready
    fit = prepared.fit()
    assert fit.workflow_valid
    assert fit.parameters[0].final == pytest.approx(0.0, abs=2.0e-8)

    transferred = fit.to_parameters(
        parameters,
        bundle_id="gross-2001-methane-ethane-slice-fit",
        bundle_version=1,
        source=SourceRecord(
            "slice-regression",
            "Unified Engine deterministic Data-to-Regression transfer slice",
            "regression output accepted only by this test slice",
        ),
        domain=ValidityDomain(
            "slice-domain",
            "reported-conditions",
            temperature_min=220.0 * unit_registry.kelvin,
            temperature_max=280.0 * unit_registry.kelvin,
        ),
        locator="tests/integration/test_scientific_slice.py",
    )
    fitted_model = Mixture(transferred)
    fitted_pair = next(record for record in transferred.records if record.family == "k_ij")
    assert float(fitted_pair.value) == pytest.approx(fit.parameters[0].final)
    assert fitted_model.parameter_fingerprint != model.parameter_fingerprint

    state = fitted_model.state(
        T=243.61 * unit_registry.kelvin,
        rho=10_000.0 * unit_registry.mole / unit_registry.meter**3,
        x=(0.5627, 0.4373),
    )
    assert math.isfinite(float(state.pressure.to("pascal").magnitude))
    flash = equilibrium.tp_flash(
        fitted_model,
        243.61 * unit_registry.kelvin,
        6_691_000.0 * unit_registry.pascal,
        (0.5627, 0.4373),
    )
    assert flash.parameter_fingerprint == fitted_model.parameter_fingerprint
