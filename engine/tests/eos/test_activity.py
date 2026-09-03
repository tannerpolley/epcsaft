import json
import math

from epcsaft import Mixture, Parameters

from tests.packets import ORACLE, parameters_json

# Figiel, Yu & Held 2025 SI S2: retained NaBr candidates require AARD < 10%.
SOURCE_ADEQUACY_AARD = 0.10


def test_aqueous_activity_reproduces_source_rows() -> None:
    rows = json.loads(ORACLE.read_text())["cases"]["aqueous_alkali_halide_activity"]
    grouped: dict[str, list[dict[str, object]]] = {}
    for row in rows:
        grouped.setdefault(row["salt"], []).append(row)
    for salt, salt_rows in grouped.items():
        assert len(salt_rows) == 6
        first = salt_rows[0]
        assert all(row["temperature_k"] == first["temperature_k"] for row in salt_rows)
        assert all(row["pressure_pa"] == first["pressure_pa"] for row in salt_rows)
        cation = str(first["cation"])
        anion = str(first["anion"])
        parameters = Parameters.from_json(
            parameters_json("figiel-2025-reference-electrolytes", 3),
            components=("water", cation, anion),
        )
        temperature = first["temperature_k"]
        pressure = first["pressure_pa"]
        assert isinstance(temperature, (int, float))
        assert isinstance(pressure, (int, float))
        reference = Mixture(parameters).reference(
            float(temperature),
            float(pressure),
            solvent_fractions=(1.0, 0.0, 0.0),
            ions=(cation, anion),
            stoichiometry=(1, 1),
        )
        deviations: list[float] = []
        for row in salt_rows:
            molality = row["molality_mol_kg"]
            measured_value = row["measured_mean_ionic_activity_coefficient"]
            assert isinstance(molality, (int, float))
            assert isinstance(measured_value, (int, float))
            calculated = reference.activity(float(molality)).mean_ionic_molality_scale
            measured = float(measured_value)
            assert math.isfinite(calculated) and calculated > 0.0
            deviations.append(abs(calculated / measured - 1.0))
            print(salt, row["molality_mol_kg"], measured, calculated, deviations[-1])
        aard = sum(deviations) / len(deviations)
        assert aard <= SOURCE_ADEQUACY_AARD, f"{salt}: {deviations}"


def test_published_model_curve_reproduces_source() -> None:
    rows = json.loads(ORACLE.read_text())["cases"]["aqueous_lithium_bromide_published_model"]
    assert len(rows) == 27
    first = rows[0]
    parameters = Parameters.from_json(
        parameters_json("figiel-2025-reference-electrolytes", 3),
        components=("water", "lithium-cation", "bromide-anion"),
    )
    reference = Mixture(parameters).reference(
        first["temperature_k"], first["pressure_pa"],
        solvent_fractions=(1.0, 0.0, 0.0),
        ions=("lithium-cation", "bromide-anion"), stoichiometry=(1, 1),
    )
    deviations: list[float] = []
    uncertainty_ratios: list[float] = []
    for row in rows:
        calculated = reference.activity(row["molality_mol_kg"]).mean_ionic_molality_scale
        published = row["published_mean_ionic_activity_coefficient"]
        uncertainty = row["digitization_uncertainty"]
        deviation = abs(calculated / published - 1.0)
        deviations.append(deviation)
        uncertainty_ratios.append(abs(calculated - published) / uncertainty)
        print(row["molality_mol_kg"], published, calculated, uncertainty_ratios[-1])
    assert sum(deviations) / len(deviations) <= sum(
        row["digitization_uncertainty"] / row["published_mean_ionic_activity_coefficient"]
        for row in rows
    ) / len(rows)
