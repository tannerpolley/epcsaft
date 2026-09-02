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
            parameters_json("figiel-2025-reference-electrolytes", 2),
            components=("water", cation, anion),
        )
        reference = Mixture(parameters).reference(
            float(first["temperature_k"]),
            float(first["pressure_pa"]),
            solvent_fractions=(1.0, 0.0, 0.0),
            ions=(cation, anion),
        )
        deviations: list[float] = []
        for row in salt_rows:
            calculated = reference.activity(float(row["molality_mol_kg"])).mean_ionic_molality_scale
            measured = float(row["measured_mean_ionic_activity_coefficient"])
            assert math.isfinite(calculated) and calculated > 0.0
            deviations.append(abs(calculated / measured - 1.0))
            print(salt, row["molality_mol_kg"], measured, calculated, deviations[-1])
        aard = sum(deviations) / len(deviations)
        assert aard <= SOURCE_ADEQUACY_AARD, f"{salt}: {deviations}"
