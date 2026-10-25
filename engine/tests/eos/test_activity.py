import csv
import hashlib
import math
from pathlib import Path
from typing import Any

from epcsaft import Mixture, Parameters, Stop
from tests.native_inputs import packet_metadata, same
from tests.packets import parameters_json

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
HAMER_WU_IONS = {
    "LiCl": ("lithium-cation", "chloride-anion"),
    "LiBr": ("lithium-cation", "bromide-anion"),
    "KCl": ("potassium-cation", "chloride-anion"),
    "KBr": ("potassium-cation", "bromide-anion"),
}
HAMER_WU_MOLALITIES = (0.001, 0.01, 0.1, 0.5, 1.0, 2.0)
TEMPERATURE_K, PRESSURE_PA = 298.15, 100000.0


def test_solvent_activity_transport() -> None:
    from dataclasses import replace
    import pytest
    from epcsaft import ActiveParameter, ParameterError, ParameterFamily, PropertyActionStatus, Pure

    components = ("water", "sodium-cation", "magnesium-cation", "chloride-anion")
    parameters = Parameters.from_json(parameters_json("figiel-2025-chloride-activity-refit", 2), components=components)
    model = Mixture(parameters, active_parameters=[ActiveParameter(ParameterFamily.DispersionEnergyOverK, [0])])
    ref = model.reference(298.15, 1e5, solvent_fractions=(1,0,0,0))
    amounts = (1/parameters.molar_masses[0], 0, 1, 2)
    result = ref.solvent_activity(x=tuple(n/sum(amounts) for n in amounts))
    # Frozen pre-extension one-salt helper witness, reconstructed natively from Phi and its column.
    assert result.solvent_activity == pytest.approx(0.93687139131454011, abs=1e-10)
    assert result.solvent_activity_parameter_derivative == pytest.approx([-6.6051908626038093e-05], abs=1e-10)
    assert result.parameter_status == [PropertyActionStatus.Available]
    assert replace(ref, T=300).solvent_activity(x=(1,0,0,0)).cause == "SolventReferenceMismatch"
    with pytest.raises(ParameterError, match="IonicReferenceUnavailable"):
        ref.activity(1)
    with pytest.raises(ParameterError, match="both ions"):
        model.reference(298.15,1e5,solvent_fractions=(1,0,0,0),ions=("sodium-cation","chloride-anion"))
    for x in ((0.9,0.1,0,0), (1.1,0,0,-0.1), (0.9,0,0,0), (float("nan"),0,0,0)):
        with pytest.raises(ParameterError):
            ref.solvent_activity(x=x)
    from epcsaft import _core

    from tests.equilibrium.test_trace_reliability import _mea_mapping

    background_components = ("water", "protonated-monoethanolamine", "carbamate-anion", "hydronium-cation", "hydroxide-anion")
    background = Mixture(Parameters.from_mapping(_mea_mapping(parameters_json("mea-co2-h2o-nine-species-estimation", 4)),
                                                components=background_components))
    with pytest.raises(_core.IonicReferenceBackground):
        background.reference(298.15, 1e5, solvent_fractions=(0.8, 0, 0, 0.1, 0.1),
                             ions=background_components[1:3], stoichiometry=(1, 1))
    water = Pure(Parameters.from_json(parameters_json("figiel-2025-chloride-activity-refit",2),components=("water",)))
    pure = water.reference(298.15,1e5,solvent_fractions=(1,)).solvent_activity(x=(1,))
    assert pure.solvent_activity == 1 and pure.osmotic_cause == "OsmoticUndefinedPureSolvent"
    assert pure.osmotic_coefficient == Stop.converged_to_floor and pure.osmotic_parameter_derivative == []


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
    return Parameters.from_json(path, components=components), packet_metadata(path, components)


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


def _csv_rows(path: Path) -> list[dict[str, Any]]:
    with path.open(newline="") as stream:
        return [
            {key: value if key == "salt" else float(value) for key, value in row.items()}
            for row in csv.DictReader(stream)
        ]


_AQUEOUS_ACTIVITY_GROUPS: dict[str, list[dict[str, Any]]] = {}
for _row in _csv_rows(HAMER_WU_OBSERVATIONS):
    if _row["salt"] in HAMER_WU_IONS and _row["molality_mol_kg"] in HAMER_WU_MOLALITIES:
        _AQUEOUS_ACTIVITY_GROUPS.setdefault(str(_row["salt"]), []).append(_row)
_PUBLISHED_MODEL_ROWS = _csv_rows(FIGIEL_DIGITIZED)


def test_aqueous_activity_reproduces_source_rows(evidence: Any) -> None:
    for salt, salt_rows in _AQUEOUS_ACTIVITY_GROUPS.items():
        assert len(salt_rows) == len(HAMER_WU_MOLALITIES)
        cation, anion = HAMER_WU_IONS[salt]
        components = ("water", cation, anion)
        temperature, pressure = TEMPERATURE_K, PRESSURE_PA
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
            calculated = reference.activity(molality).mean_ionic_molality_scale
            measured = row["gamma_pm_m"]
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
    components = ("water", "lithium-cation", "bromide-anion")
    temperature, pressure = TEMPERATURE_K, PRESSURE_PA
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
        published = row["gamma_pm_m"]
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
        row["digitization_uncertainty"] / row["gamma_pm_m"]
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
    path = parameters_json("held-cameretti-sadowski-2008", 3)
    mixture = Mixture(Parameters.from_json(path, components=["water", "sodium-cation", "chloride-anion"]))
    reference = mixture.reference(298.15, 100000.0, solvent_fractions=[1, 0, 0],
                                  ions=("sodium-cation", "chloride-anion"), stoichiometry=(1, 1))
    # Limiting law 1 - phi = A sqrt(m) / 3 from the model's own water density and
    # permittivity; at 1e-6 mol/kg the next Debye-Hueckel term, about 1.5 kappa a,
    # is 1.5e-3 of it, so 3e-3 bounds it twice over.
    pure = mixture.state(298.15, P=100000.0, x=[1, 0, 0], phase="liquid")
    e, eps0, kb, na, mw = 1.602176634e-19, 8.8541878128e-12, 1.380649e-23, 6.02214076e23, 0.0180153
    eps = eps0 * pure.bulk_relative_permittivity
    a = e * e / (8 * math.pi * eps * kb * 298.15) * math.sqrt(
        2 * e * e * na * pure.molar_density * mw / (eps * kb * 298.15)
    )
    phi = reference.activity(1e-6).osmotic_coefficient
    assert isinstance(phi, float)
    limit = a * math.sqrt(1e-6) / 3
    assert abs((1 - phi) - limit) <= 3e-3 * limit
    assert isinstance(reference.activity(1e-8).osmotic_coefficient, float)
    # Kept refusal (#353): below the #214 floor, ln a_w rounding in phi (8.9e-6 at 1e-9 mol/kg) rivals 1 - phi (1.2e-5).
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
        r.coordinate(p, [("k_ij", pair, 1)], origin=0, scale=1)
        for pair in (("sodium-cation", "water"), ("sodium-cation", "chloride-anion"))
    ]
    cs += [
        r.coordinate(
            p, [("ion_suppression_coefficient", "sodium-cation", 1), ("ion_suppression_coefficient", "chloride-anion", 1)], origin=0, scale=1
        )
    ]
    cs += [
        r.coordinate(p, [(family, ion, 1)], origin=0, scale=1)
        for family, ion in (
            ("born_diameter", "sodium-cation"),
            ("debye_huckel_diameter", "chloride-anion"),
            ("packing_diameter", "chloride-anion"),
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
    cs[4] = r.coordinate(p, [("debye_huckel_diameter", "sodium-cation", 1)], origin=0, scale=1)
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


def test_reference_column_is_the_derivative_of_its_limit(evidence: Any) -> None:
    """#594 F6: the ion reference's water epsilon/k column equals a Richardson (h, h/2) central difference of its
    rebuilt-model log fugacities, h = 1e-3 epsilon; a column left on the 1e-12 rung misses by 4.1e-7 (Na+)."""
    import json
    from epcsaft import ActiveParameter, ParameterFamily

    components = ("water", "sodium-cation", "chloride-anion")
    path = parameters_json("held-2014-strategy2", 1)
    document = json.loads(path.read_text())
    water = next(c for c in document["components"] if c["component_id"] == "water")["coefficients"]
    energy = next(c for c in water if c["family"] == "dispersion_energy_over_k")["value"]
    epsilon = energy["magnitude"]

    def reference(value: float) -> Any:
        energy["magnitude"] = value
        model = Mixture(Parameters.from_mapping(document, components=components),
                        active_parameters=[ActiveParameter(ParameterFamily.DispersionEnergyOverK, [0])])
        return model.reference(298.15, 101325.0, solvent_fractions=(1, 0, 0), ions=components[1:],
                               stoichiometry=(1, 1))._native

    def central(h: float) -> list[float]:
        high, low = reference(epsilon + h).log_fugacity, reference(epsilon - h).log_fugacity
        return [(a - b) / (2 * h) for a, b in zip(high, low)]

    h = 1e-3 * epsilon
    coarse, fine = central(h), central(h / 2)
    from epcsaft import PropertyActionStatus

    column = reference(epsilon).parameter_derivative[0]
    assert all(c.status == PropertyActionStatus.Available and c.action is not None for c in column)
    assert all(c.diagnostics.site_rank > 0 and c.diagnostics.association_residual <= 1e-10 for c in column)
    for ion, exact, d1, d2 in zip(components[1:], column, coarse, fine):
        same(evidence, f"reference-column/{ion}", exact.action, (4 * d2 - d1) / 3, atol=0.0, rtol=1e-8,
             **packet_metadata(path, components), T=298.15, P=101325.0, step=h)


def test_electrolyte_density_columns_and_debye_density_limit() -> None:
    """#219 F1a density topology and F4(b), from the pressure-closed density owner."""
    from epcsaft import _core, regression as r

    components = ("water", "sodium-cation", "chloride-anion")
    selection = REPOSITORY_ROOT / "analyses/2026-dielectric-alternatives/inputs/selections/density-temperature.json"
    p = Parameters.from_json(selection, components=components)
    cs = [r.coordinate(p, [("dispersion_energy_over_k", "water", 1)], origin=0, scale=1)]
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
    cs = [r.coordinate(p, [("dispersion_energy_over_k", "water", 1)], origin=0, scale=1)]
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


def test_henry_constant_matches_clapeyron(evidence: Any) -> None:
    """#566: fugacity-based H^px = lim f/x of methane in liquid ethane, 200 K and 1 MPa, Gross 2001 /3 parameters.

    Clapeyron 0.6.29 values from the retained #566 comparator henry-reference.json (SHA-256 d4daa4ef...a150):
    raw L = ln(rho/rho0) + mu_res/RT on the liquid root at each molality, its terminal limit and the exact-zero insertion.
    """
    from epcsaft.equilibrium import NeutralReferenceStatus

    model = Mixture(Parameters.from_json(parameters_json("gross-2001-methane-ethane", 3), components=("methane", "ethane")))
    clapeyron_rungs = {1e-6: 8.013550838967435, 1e-8: 8.013550854734689, 1e-10: 8.013550854892356,
                       1e-12: 8.013550854893948, 1e-14: 8.01355085489396}
    inputs = {"T_K": 200.0, "P_Pa": 1e6, "solvent": "ethane", "solute": "methane", "packet": "gross-2001-methane-ethane/3"}
    for molality, reference in clapeyron_rungs.items():
        x = molality / (1.0 / 0.03007 + molality)
        state = model.state(200.0, P=1e6, x=(x, 1.0 - x), phase="liquid")
        potential = state.residual_chemical_potential_over_rt[0]
        assert potential is not None
        same(evidence, f"henry/rung/{molality}", math.log(state.molar_density) + potential, reference, 1e-10, 0.0,
             inputs={**inputs, "molality_mol_kg": molality})
    (henry,) = model.reference(200.0, 1e6, solvent_fractions=(0.0, 1.0)).henry_constants(solutes=["methane"])
    assert henry.status == NeutralReferenceStatus.Available and henry.component_index == 0, henry.message
    for name, reference in (("terminal_limit", 5024641.690506341), ("zero_insertion", 5024641.69050634)):
        same(evidence, f"henry/{name}", henry.henry_constant_px, reference, 0.0, 1e-10, inputs=inputs)
    # Charge-balanced ion slots pass the reference ratio validation; Henry admits all-neutral models only.
    figiel = Mixture(Parameters.from_json(parameters_json("figiel-yu-held-2025", 3),
                                          components=("water", "sodium-cation", "chloride-anion")))
    (refused,) = figiel.reference(298.15, 1e5, solvent_fractions=(1.0, 0.0, 0.0)).henry_constants(solutes=["sodium-cation"])
    assert refused.status == NeutralReferenceStatus.NonNeutralDeclaration and refused.henry_constant_px is None
