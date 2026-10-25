#!/usr/bin/env python3
"""Verify the Issue #118 reaction decision table with Python's standard library."""

from __future__ import annotations

import json
import math
from pathlib import Path

AUDIT = Path(__file__).resolve().parents[1] / "results/analysis/reaction-equilibria/reaction-audit.json"
ELEMENTS = {
    "CO2": (1, 0, 2, 0, 0),
    "MEA": (2, 7, 1, 1, 0),
    "H2O": (0, 2, 1, 0, 0),
    "MEAH+": (2, 8, 1, 1, 1),
    "MEACOO-": (3, 6, 3, 1, -1),
    "HCO3-": (1, 1, 3, 0, -1),
    "CO3^2-": (1, 0, 3, 0, -2),
    "H3O+": (0, 3, 1, 0, 1),
    "OH-": (0, 1, 1, 0, -1),
}
EXPECTED_STOICHIOMETRY = {
    "R1": [0, 0, -2, 0, 0, 0, 0, 1, 1],
    "R2": [-1, 0, -2, 0, 0, 1, 0, 1, 0],
    "R3": [0, 0, -1, 0, 0, -1, 1, 1, 0],
    "R4": [0, 1, -1, 0, -1, 1, 0, 0, 0],
    "R5": [0, 1, -1, -1, 0, 0, 0, 1, 0],
}
EXPECTED_LN_K_298_15 = {
    "R1": -32.23229024933365,
    "R2": -14.642290490877887,
    "R3": -23.79214802087026,
    "R4": -3.031961596511823,
    "R5": -21.86574617291778,
}
EXPECTED_DECISIONS = {
    "R1": (
        {"A": 132.899, "B_K": -13445.9, "C": -22.4773, "D_per_K": 0.0},
        [273.15, 498.15],
        [("source-qualified", 293.15, True, 393.15, True)],
    ),
    "R2": (
        {"A": 231.465, "B_K": -12092.1, "C": -36.7816, "D_per_K": 0.0},
        [273.15, 498.15],
        [("source-qualified", 293.15, True, 393.15, True)],
    ),
    "R3": (
        {"A": 216.049, "B_K": -12431.7, "C": -35.4819, "D_per_K": 0.0},
        [273.15, 498.15],
        [("source-qualified", 293.15, True, 393.15, True)],
    ),
    "R4": (
        {"A": 2.151, "B_K": -1545.3},
        [293.0, 323.0],
        [
            ("source-qualified", 293.15, True, 323.0, True),
            ("extrapolated", 323.0, False, 393.15, True),
        ],
    ),
    "R5": (
        {"A_K": 2677.91, "B": 0.3869, "C_per_K": 0.0004277},
        [273.15, 323.15],
        [
            ("source-qualified", 293.15, True, 323.15, True),
            ("extrapolated", 323.15, False, 393.15, True),
        ],
    ),
}
EXPECTED_CLASSIFICATIONS = {
    "R1": ["source-qualified"] * 8,
    "R2": ["source-qualified"] * 8,
    "R3": ["source-qualified"] * 8,
    "R4": ["source-qualified"] * 3 + ["extrapolated"] * 5,
    "R5": ["source-qualified"] * 4 + ["extrapolated"] * 4,
}


def ln_k(reaction: dict, temperature_k: float) -> float:
    parameters = reaction["correlation"]["parameters"]
    form = reaction["correlation"]["form"]
    if form == "A + B_K/T + C*ln(T/K) + D_per_K*T":
        value = (
            parameters["A"]
            + parameters["B_K"] / temperature_k
            + parameters["C"] * math.log(temperature_k)
            + parameters["D_per_K"] * temperature_k
        )
    elif form == "A + B_K/T":
        value = parameters["A"] + parameters["B_K"] / temperature_k
    elif form == "-ln(10)*(A_K/T + B + C_per_K*T)":
        value = -math.log(10.0) * (
            parameters["A_K"] / temperature_k
            + parameters["B"]
            + parameters["C_per_K"] * temperature_k
        )
    else:
        raise AssertionError(f"unsupported correlation form: {form}")
    return value + reaction["source_to_common_ln_k_offset"]


def classify(reaction: dict, temperature_k: float) -> str:
    for segment in reaction["target_segments"]:
        above_lower = temperature_k > segment["lower_k"] or (
            segment["lower_closed"] and temperature_k == segment["lower_k"]
        )
        below_upper = temperature_k < segment["upper_k"] or (
            segment["upper_closed"] and temperature_k == segment["upper_k"]
        )
        if above_lower and below_upper:
            return segment["classification"]
    raise AssertionError(f"{reaction['id']} leaves {temperature_k} K unclassified")


def main() -> None:
    audit = json.loads(AUDIT.read_text())
    assert audit["target_temperature_k"] == [293.15, 393.15]
    assert audit["species_order"] == list(ELEMENTS)
    assert [reaction["id"] for reaction in audit["reactions"]] == list(
        EXPECTED_STOICHIOMETRY
    )

    water = audit["common_standard_state"]
    offset = -math.log(
        water["water_molar_mass_kg_per_mol"] * water["standard_molality_mol_per_kg"]
    )
    expected_offsets = {
        "R1": 2 * offset,
        "R2": offset,
        "R3": offset,
        "R4": 0.0,
        "R5": 0.0,
    }

    for reaction in audit["reactions"]:
        reaction_id = reaction["id"]
        stoichiometry = reaction["stoichiometry"]
        expected_parameters, expected_domain, expected_segments = EXPECTED_DECISIONS[
            reaction_id
        ]
        assert stoichiometry == EXPECTED_STOICHIOMETRY[reaction_id]
        assert reaction["correlation"]["parameters"] == expected_parameters
        assert reaction["source_domain_k"] == expected_domain
        actual_segments = [
            (
                segment["classification"],
                segment["lower_k"],
                segment["lower_closed"],
                segment["upper_k"],
                segment["upper_closed"],
            )
            for segment in reaction["target_segments"]
        ]
        assert actual_segments == expected_segments
        for property_index in range(5):
            assert (
                sum(
                    coefficient * ELEMENTS[species][property_index]
                    for species, coefficient in zip(
                        audit["species_order"], stoichiometry
                    )
                )
                == 0
            )
        assert math.isclose(
            reaction["source_to_common_ln_k_offset"],
            expected_offsets[reaction_id],
            abs_tol=1e-14,
        )
        assert math.isclose(
            ln_k(reaction, 298.15), EXPECTED_LN_K_298_15[reaction_id], abs_tol=1e-12
        )
        classifications = []
        for temperature_k in audit["audit_temperatures_k"]:
            assert math.isfinite(ln_k(reaction, temperature_k))
            classifications.append(classify(reaction, temperature_k))
        assert classifications == EXPECTED_CLASSIFICATIONS[reaction_id]
        assert not reaction["fit"]["eligible"]
        assert reaction["fit"]["coordinates"] == []

    r1 = audit["reactions"][0]
    assert r1["correlation"]["parameters"]["B_K"] == -13445.9
    assert r1["rejected_transcription"]["value"] == 13445.9
    assert not audit["fit_decision"]["pressure_only_refit_sufficient"]
    assert audit["fit_decision"]["eligible_coordinates"] == []
    assert audit["fit_decision"]["source_backed_priors_or_bounds"] == []
    print(
        "PASS: R1-R5 balances, standard states, correlations, domains, and fit decisions verified"
    )


if __name__ == "__main__":
    main()
