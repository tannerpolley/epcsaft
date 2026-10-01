import hashlib
import json
import math
from pathlib import Path
from typing import Any

from epcsaft import Mixture, Parameters, Stop
from tests.packets import ORACLE, parameters_json

# Figiel, Yu & Held 2025 SI S2: retained NaBr candidates require AARD < 10%.
SOURCE_ADEQUACY_AARD = 0.10

REPOSITORY_ROOT = Path(__file__).resolve().parents[3]
HAMER_WU_PACKET = REPOSITORY_ROOT / "data/packets/hamer-wu-1972-aqueous-alkali-halides/1"
HAMER_WU_OBSERVATIONS = HAMER_WU_PACKET / "observations/hamer-wu-1972-aqueous-alkali-halides.csv"
HAMER_WU_SOURCE = HAMER_WU_PACKET / "sources/hamer-wu-1972-aqueous-alkali-halides.yaml"
FIGIEL_DIGITIZED = (
    REPOSITORY_ROOT
    / "analyses/reproduction/2025-figiel/figures/figure-05/source/libr-model-digitized.csv"
)
HAMER_WU_TABLES = {"LiCl": 9, "LiBr": 10, "NaCl": 16, "NaBr": 17, "KCl": 28, "KBr": 29}


def _file_identity(path: Path) -> dict[str, Any]:
    return {
        "path": str(path.relative_to(REPOSITORY_ROOT)),
        "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "size_bytes": path.stat().st_size,
    }


def _parameter_inputs(
    packet_id: str, version: int, components: tuple[str, ...]
) -> tuple[Parameters, dict[str, Any]]:
    path = parameters_json(packet_id, version)
    document = json.loads(path.read_text())
    parameters = Parameters.from_mapping(document, components=components)
    model_choices = {
        str(item["kind"]): str(item["choice"])
        for item in document.get("model_families", ())
        if isinstance(item, dict) and "kind" in item and "choice" in item
    }
    return parameters, {
        "packet_id": packet_id,
        "packet_version": version,
        "parameter_file": _file_identity(path),
        "effective_parameters": {
            "supplied_mapping": document,
            "selected_components": list(components),
            "model_choices": model_choices,
        },
    }


def _hamer_wu_reference(salt: str) -> dict[str, Any]:
    return {
        "source_id": "hamer-wu-1972",
        "citation": (
            "Hamer and Wu, Osmotic Coefficients and Mean Activity Coefficients "
            "of Uni-univalent Electrolytes in Water at 25 C (1972)"
        ),
        "doi": "10.1063/1.3253108",
        "domain": "aqueous-alkali-halides-298.15-k",
        "classification": "critically evaluated reference data",
        "conditions": "298.15 K, 1 bar, water",
        "quantity": "stoichiometric mean ionic activity coefficient on the molality scale",
        "units": "formula-unit molality in mol/kg; gamma_pm_m dimensionless",
        "uncertainty": "not reported",
        "source_table": HAMER_WU_TABLES[salt],
        "owner_files": {
            "packet": _file_identity(HAMER_WU_PACKET / "packet.toml"),
            "observations": _file_identity(HAMER_WU_OBSERVATIONS),
            "source_metadata": _file_identity(HAMER_WU_SOURCE),
            "retained_case_owner": _file_identity(ORACLE),
        },
    }


def _figiel_model_reference() -> dict[str, Any]:
    return {
        "source_id": "figiel-yu-held-2025",
        "citation": (
            "Figiel, Yu, and Held, Predicting Thermodynamic Properties of Ions "
            "in Single Solvents and in Mixed Solvents Using a Modified Born Term "
            "within the ePC-SAFT Framework (2025)"
        ),
        "doi": "10.1021/acs.iecr.5c00475",
        "domain": "figiel-298-k-1-bar",
        "conditions": "298.15 K, 1 bar, LiBr in water",
        "quantity": "published LiBr model curve, molality-scale mean ionic activity coefficient",
        "units": "formula-unit molality in mol/kg; gamma_pm_m dimensionless",
        "classification": "published-model digitization",
        "uncertainty": "per-point digitization uncertainty retained",
        "figure": "Figure 5, right panel, orange LiBr model line",
        "comparison_scope": "curve reproduction; equation matching is not established",
        "equation_matching_established": False,
        "owner_files": {
            "digitized_curve": _file_identity(FIGIEL_DIGITIZED),
            "retained_case_owner": _file_identity(ORACLE),
        },
    }


def _aard_criterion(
    *, threshold: float, verdict: bool, rule: str, rationale: str
) -> dict[str, Any]:
    return {
        "type": "custom",
        "helper": "engine.tests.eos.test_activity._aard_criterion",
        "rule": rule,
        "verdict": verdict,
        "threshold": threshold,
        "scale": "AARD",
        "rationale": rationale,
    }


_AQUEOUS_ACTIVITY_ROWS: list[dict[str, Any]] = json.loads(ORACLE.read_text())["cases"]["aqueous_alkali_halide_activity"]
_AQUEOUS_ACTIVITY_GROUPS: dict[str, list[dict[str, Any]]] = {}
for _row in _AQUEOUS_ACTIVITY_ROWS:
    _AQUEOUS_ACTIVITY_GROUPS.setdefault(_row["salt"], []).append(_row)
_PUBLISHED_MODEL_ROWS: list[dict[str, Any]] = json.loads(ORACLE.read_text())["cases"]["aqueous_lithium_bromide_published_model"]


def test_aqueous_activity_reproduces_source_rows(evidence: Any) -> None:
    for salt, salt_rows in _AQUEOUS_ACTIVITY_GROUPS.items():
        assert len(salt_rows) == 6
        first = salt_rows[0]
        assert all(row["temperature_k"] == first["temperature_k"] for row in salt_rows)
        assert all(row["pressure_pa"] == first["pressure_pa"] for row in salt_rows)
        cation = str(first["cation"])
        anion = str(first["anion"])
        components = ("water", cation, anion)
        temperature_value = first["temperature_k"]
        pressure_value = first["pressure_pa"]
        assert isinstance(temperature_value, (int, float))
        assert isinstance(pressure_value, (int, float))
        temperature = float(temperature_value)
        pressure = float(pressure_value)
        solvent_fractions: tuple[float, float, float] = (1.0, 0.0, 0.0)
        ions: tuple[str, str] = (cation, anion)
        stoichiometry: tuple[int, int] = (1, 1)
        invocation: dict[str, Any] = {
            "temperature_K": temperature,
            "pressure_Pa": pressure,
            "components": list(components),
            "solvent_fractions": list(solvent_fractions),
            "ions": list(ions),
            "stoichiometry": list(stoichiometry),
        }
        parameters, parameter_metadata = _parameter_inputs(
            "figiel-2025-reference-electrolytes", 3, components
        )
        reference = Mixture(parameters).reference(
            temperature,
            pressure,
            solvent_fractions=solvent_fractions,
            ions=ions,
            stoichiometry=stoichiometry,
        )
        deviations: list[float] = []
        for row in salt_rows:
            molality = row["molality_mol_kg"]
            measured_value = row["measured_mean_ionic_activity_coefficient"]
            assert isinstance(molality, (int, float))
            assert isinstance(measured_value, (int, float))
            calculated = reference.activity(float(molality)).mean_ionic_molality_scale
            measured = float(measured_value)
            evidence.observe(
                f"activity_row/{salt}/{len(deviations)}",
                quantity="mean ionic activity coefficient",
                model=float(calculated),
                reference=measured,
                unit="1",
                inputs={
                    **invocation,
                    "salt": salt,
                    "molality_mol_kg": float(molality),
                    "_parameters": parameter_metadata,
                    "_reference": _hamer_wu_reference(salt),
                },
                basis="Experimental",
                method="retained source row observation",
                claim="Retain the model and evaluated reference operands; the aggregate row carries the test-owned adequacy criterion.",
                comparison_group_id=f"activity/{salt}/gamma",
                point_id=f"molality-{len(deviations)}",
                quantity_latex=r"\gamma_{\pm}",
                relationship_latex=r"\gamma_{\pm}^{model} \text{ observed against } \gamma_{\pm}^{reference}",
                inputs_latex=r"T,\ p,\ m",
            )
            assert math.isfinite(calculated) and calculated > 0.0
            deviations.append(abs(calculated / measured - 1.0))
            print(salt, row["molality_mol_kg"], measured, calculated, deviations[-1])
        aard = sum(deviations) / len(deviations)
        evidence.compare(
            f"activity_aard/{salt}",
            quantity="mean ionic activity coefficient AARD",
            model=aard,
            reference=SOURCE_ADEQUACY_AARD,
            unit="1",
            inputs={
                **invocation,
                "salt": salt,
                "rows": len(salt_rows),
                "molality_mol_kg": [float(row["molality_mol_kg"]) for row in salt_rows],
                "_parameters": parameter_metadata,
                "_reference": _hamer_wu_reference(salt),
            },
            criterion=_aard_criterion(
                threshold=SOURCE_ADEQUACY_AARD,
                verdict=aard <= SOURCE_ADEQUACY_AARD,
                rule="AARD(model) <= test-owned source adequacy bar(reference)",
                rationale="Existing 10% source adequacy bar retained from this test.",
            ),
            basis="Experimental",
            method="AARD",
            claim="Retained evaluated reference rows satisfy the test-owned 10% AARD adequacy bar.",
            comparison_group_id=f"activity/{salt}/aard",
            point_id="aggregate",
            quantity_latex=r"\operatorname{AARD}(\gamma_{\pm})",
            relationship_latex=r"\operatorname{AARD}(\gamma_{\pm}) \leq \tau_{source}",
            inputs_latex=r"T,\ p,\ \{m_i\}",
        )
        assert aard <= SOURCE_ADEQUACY_AARD, f"{salt}: {deviations}"


def test_published_model_curve_reproduces_source(evidence: Any) -> None:
    """Class-3 published-model reproduction; not solver numerical accuracy.

    The retained digitization bar compares model curves only.
    """
    rows = _PUBLISHED_MODEL_ROWS
    assert len(rows) == 27
    first = rows[0]
    components = ("water", "lithium-cation", "bromide-anion")
    temperature = float(first["temperature_k"])
    pressure = float(first["pressure_pa"])
    solvent_fractions: tuple[float, float, float] = (1.0, 0.0, 0.0)
    ions: tuple[str, str] = ("lithium-cation", "bromide-anion")
    stoichiometry: tuple[int, int] = (1, 1)
    invocation: dict[str, Any] = {
        "temperature_K": temperature,
        "pressure_Pa": pressure,
        "components": list(components),
        "solvent_fractions": list(solvent_fractions),
        "ions": list(ions),
        "stoichiometry": list(stoichiometry),
    }
    parameters, parameter_metadata = _parameter_inputs(
        "figiel-2025-reference-electrolytes", 3, components
    )
    reference = Mixture(parameters).reference(
        temperature,
        pressure,
        solvent_fractions=solvent_fractions,
        ions=ions,
        stoichiometry=stoichiometry,
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
        evidence.observe(
            f"published_row/{len(deviations) - 1}",
            quantity="mean ionic activity coefficient",
            model=float(calculated),
            reference=float(published),
            unit="1",
            inputs={
                **invocation,
                "molality_mol_kg": float(row["molality_mol_kg"]),
                "digitization_uncertainty": float(uncertainty),
                "_parameters": parameter_metadata,
                "_reference": _figiel_model_reference(),
            },
            basis="Cross-model sanity",
            method="digitized source row observation",
            claim="Retain the installed model and digitized published-model operands; the aggregate envelope row carries the criterion.",
            comparison_group_id="published-curve/gamma",
            point_id=f"molality-{len(deviations) - 1}",
            quantity_latex=r"\gamma_{\pm}",
            relationship_latex=r"\gamma_{\pm}^{model} \text{ observed against } \gamma_{\pm}^{published}",
            inputs_latex=r"T,\ p,\ m",
        )
        print(row["molality_mol_kg"], published, calculated, uncertainty_ratios[-1])
    actual = sum(deviations) / len(deviations)
    expected = sum(
        row["digitization_uncertainty"] / row["published_mean_ionic_activity_coefficient"]
        for row in rows
    ) / len(rows)
    evidence.compare(
        "published_curve_aard",
        quantity="published activity curve AARD",
        model=actual,
        reference=expected,
        unit="1",
        inputs={
            **invocation,
            "rows": len(rows),
            "uncertainty_basis": "digitization uncertainty",
            "_parameters": parameter_metadata,
            "_reference": _figiel_model_reference(),
        },
        criterion=_aard_criterion(
            threshold=expected,
            verdict=actual <= expected,
            rule="AARD(model) <= supplied digitization envelope(reference)",
            rationale="Retain the supplied per-point digitization uncertainty as the existing curve envelope.",
        ),
        basis="Cross-model sanity",
        method="AARD against digitized curve",
        claim="The model curve remains within the supplied digitization envelope; equation matching is not established.",
        comparison_group_id="published-curve/aard",
        point_id="aggregate",
        quantity_latex=r"\operatorname{AARD}(\gamma_{\pm})",
        relationship_latex=r"\operatorname{AARD}(\gamma_{\pm}) \leq \tau_{digitization}",
        inputs_latex=r"T,\ p,\ \{m_i\}",
    )
    assert actual <= expected


def test_trace_osmotic_debye_limit_and_reporting_floor() -> None:
    # #214, Held packet 3: numerical Debye asymptote, not physical validation.
    path = REPOSITORY_ROOT / "analyses/2026-divalent-chloride-activity/inputs/held-2008/3/parameters/parameters.json"
    mixture = Mixture(Parameters.from_json(path, components=["water", "sodium-cation", "chloride-anion"]))
    reference = mixture.reference(298.15, 100000.0, solvent_fractions=[1, 0, 0],
                                  ions=("sodium-cation", "chloride-anion"), stoichiometry=(1, 1))
    phi = reference.activity(1e-8).osmotic_coefficient
    assert isinstance(phi, float)
    assert abs(phi - 0.9999606345027382) <= 1e-4
    for molality in (1e-9, 1e-10, 1e-12):
        activity = reference.activity(molality)
        assert activity.osmotic_coefficient == Stop.converged_to_floor
        assert math.isfinite(activity.mean_ionic_molality_scale)


def test_electrolyte_exact_columns_and_debye_limit() -> None:
    """#219: physical columns, zero reference columns and the NaCl Debye limit."""
    from epcsaft import _core, regression as r

    p = Parameters.from_json(
        parameters_json("figiel-2025-reference-electrolytes", 4),
        components=("water", "sodium-cation", "chloride-anion"),
    )
    cs = [
        r.coordinate(p, "k_ij", pair, origin=0, scale=1)
        for pair in (("sodium-cation", "water"), ("sodium-cation", "chloride-anion"))
    ]
    cs += [
        r.coordinate(
            p, "ion_suppression_coefficient", "sodium-cation", "chloride-anion", origin=0, scale=1
        )
    ]
    cs += [
        r.coordinate(p, family, ion, origin=0, scale=1)
        for family, ion in (
            ("born_diameter_angstrom", "sodium-cation"),
            ("debye_huckel_diameter_angstrom", "chloride-anion"),
            ("packing_diameter_angstrom", "chloride-anion"),
        )
    ]
    chart = r.CoordinateSet()
    chart.coordinates = cs
    model = Mixture(p, active_parameters=chart.active())
    u = [-0.3, 0.8, 7.01, 3.445, 2.42528, 2.42528]
    observations = [
        r.electrolyte(
            p,
            kind,
            298.15,
            1e5,
            {"water": 1},
            ions=("sodium-cation", "chloride-anion"),
            stoichiometry=(1, 1),
            molality=m,
            observed=1,
            form="difference",
            scale=1,
        )
        for m in (0.1, 1, 4.5)
        for kind in ("mean_ionic_activity_coefficient", "osmotic_coefficient")
    ]
    energy = r.electrolyte(
        p,
        "solvation_gibbs_energy",
        298.15,
        1e5,
        {"water": 1},
        ions=("sodium-cation", "chloride-anion"),
        stoichiometry=(1, 1),
        molality=0,
        species="sodium-cation",
        observed=-381930,
        form="difference",
        scale=1,
    )
    observations.append(energy)
    raw = [r.differentiate(model, ob).value for ob in observations]
    tie = chart.tie()
    for c in range(len(cs)):
        exact = [sum(value for value, j in zip(g, tie) if j == c) for g in raw]
        step = 1e-4 * max(abs(u[c]), 1)
        plus, minus = list(u), list(u)
        plus[c] += step
        minus[c] -= step
        high = _core.observations_at(p._native, cs, plus, [], observations, [0] * len(observations))
        low = _core.observations_at(p._native, cs, minus, [], observations, [0] * len(observations))
        for i, value in enumerate(exact):
            delta = 298.15 * 8.31446261815324 * 1e-7 if i == len(observations) - 1 else 1e-7
            cd = (high[i].observed - low[i].observed) / (2 * step)
            assert abs(value - cd) <= 2e-5 * max(abs(v) for v in exact) + delta
    assert abs(raw[-1][1] / (298.15 * 8.31446261815324)) <= 1e-9
    assert (
        abs(sum(raw[-1][i] for i, j in enumerate(tie) if j == 2) / (298.15 * 8.31446261815324))
        <= 1e-9
    )
    epsilon = r.electrolyte(
        p,
        "bulk_relative_permittivity",
        298.15,
        1e5,
        {"water": 0.96, "sodium-cation": 0.02, "chloride-anion": 0.02},
        ions=("sodium-cation", "chloride-anion"),
        stoichiometry=(1, 1),
        molality=0,
        observed=1,
        form="difference",
        scale=1,
    )
    prediction = r.predict(model, epsilon)
    assert prediction.value is not None
    g = r.differentiate(model, epsilon).value
    expected = -prediction.value * 0.04 / (1 + 7.01 * 0.04)
    assert abs(sum(v for v, j in zip(g, tie) if j == 2) - expected) <= 1e-12 * abs(expected)
    # Limiting columns: same five families as F4(a), with the Na DH diameter.
    cs[4] = r.coordinate(p, "debye_huckel_diameter_angstrom", "sodium-cation", origin=0, scale=1)
    chart.coordinates = cs[:5]
    model = Mixture(p, active_parameters=chart.active())
    tie = chart.tie()
    limits = []
    for m in (1e-5, 1e-6):
        observations[0].molality = m
        g = r.differentiate(model, observations[0]).value
        limits.append([sum(v for v, j in zip(g, tie) if j == c) / m for c in range(5)])
    for fine, coarse in zip(limits[1], limits[0]):
        assert abs(fine - coarse) <= 0.01 * abs(coarse)


def test_electrolyte_density_columns_and_debye_density_limit() -> None:
    """#219 F1a density topology and F4(b), from the pressure-closed density owner."""
    from epcsaft import _core, regression as r

    components = ("water", "sodium-cation", "chloride-anion")
    selection = REPOSITORY_ROOT / "analyses/2026-dielectric-alternatives/inputs/selections/density-temperature.json"
    p = Parameters.from_json(selection, components=components)
    cs = [r.coordinate(p, "dispersion_energy_over_k", "water", origin=0, scale=1)]
    chart = r.CoordinateSet()
    chart.coordinates = cs
    model = Mixture(p, active_parameters=chart.active())
    ob = r.electrolyte(p, "bulk_relative_permittivity", 298.15, 1e5,
        {"water": .96, "sodium-cation": .02, "chloride-anion": .02}, ions=components[1:],
        stoichiometry=(1, 1), molality=0, observed=1, form="difference", scale=1)
    column = r.differentiate(model, ob).value[0]
    step = 353.95*1e-4
    high = _core.observations_at(p._native, cs, [353.95 + step], [], [ob], [0])[0].observed
    low = _core.observations_at(p._native, cs, [353.95 - step], [], [ob], [0])[0].observed
    cd = (high - low) / (2 * step)
    assert column != 0 and abs(column - cd) <= 1e-6 * abs(cd)
    p = Parameters.from_json(
        parameters_json("figiel-2025-reference-electrolytes", 4), components=components
    )
    cs = [r.coordinate(p, "dispersion_energy_over_k", "water", origin=0, scale=1)]
    chart.coordinates = cs
    model = Mixture(p, active_parameters=chart.active())
    pure = model.state(298.15, P=1e5, x=(1, 0, 0), phase="liquid")
    density = r.density(
        p,
        298.15,
        1e5,
        {"water": 1},
        phase="liquid",
        observed=1,
        form="difference",
        scale=1,
        mass_basis=False,
    )
    drho = r.differentiate(model, density).value[0]
    # Debye limiting-law verification, SI constants from core/physical_constants.hpp.
    e, eps0, kb, na, mw = 1.602176634e-19, 8.8541878128e-12, 1.380649e-23, 6.02214076e23, 0.0180153
    eps = eps0 * pure.bulk_relative_permittivity
    a = (
        e
        * e
        / (8 * math.pi * eps * kb * 298.15)
        * math.sqrt(2 * e * e * na * pure.molar_density * mw / (eps * kb * 298.15))
    )
    expected = -0.5 * a * drho / pure.molar_density
    columns = []
    for m in (1e-5, 1e-6):
        gamma = r.electrolyte(
            p,
            "mean_ionic_activity_coefficient",
            298.15,
            1e5,
            {"water": 1},
            ions=components[1:],
            stoichiometry=(1, 1),
            molality=m,
            observed=1,
            form="difference",
            scale=1,
        )
        columns.append(r.differentiate(model, gamma).value[0] / math.sqrt(m))
    assert abs(columns[1] - expected) <= 0.01 * abs(expected)
    assert abs(columns[1] - expected) < abs(columns[0] - expected)
