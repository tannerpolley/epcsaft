from __future__ import annotations

import pytest
from parameter_dictionaries import FIGIEL_REFERENCE_ELECTROLYTE_PARAMETERS

import epcsaft
from epcsaft.eos import _native_transport
from epcsaft import equilibrium
from epcsaft.equilibrium import _equilibrium

WATER_MOLAR_MASS_G_PER_MOL = 18.0153
PERDOMO_TABLE3_NACL_MOL_PER_KG_WATER = 5.6
NUMERICAL_ATOL = 1.0e-8
NUMERICAL_RTOL = 1.0e-8
PHASE_FRACTION_ATOL = 5.0e-8
PHASE_VOLUME_ATOL = 5.0e-8
PHASE_COMPOSITION_ATOL = 5.0e-7


def _nacl_feed(molality_mol_per_kg: float) -> tuple[float, float, float]:
    amounts = (
        1000.0 / WATER_MOLAR_MASS_G_PER_MOL,
        molality_mol_per_kg,
        molality_mol_per_kg,
    )
    return tuple(value / sum(amounts) for value in amounts)


_source_amounts = (
    1000.0 / WATER_MOLAR_MASS_G_PER_MOL,
    PERDOMO_TABLE3_NACL_MOL_PER_KG_WATER,
    PERDOMO_TABLE3_NACL_MOL_PER_KG_WATER,
)
PERDOMO_TABLE3_FEED = _nacl_feed(PERDOMO_TABLE3_NACL_MOL_PER_KG_WATER)


def _perdomo_table3_model() -> epcsaft.Mixture:
    parameters = epcsaft.Parameters.from_dictionary(FIGIEL_REFERENCE_ELECTROLYTE_PARAMETERS)
    return epcsaft.Mixture(parameters)


def _without_times(value: object) -> object:
    if isinstance(value, dict):
        return {
            key: 0.0 if key in {"wall_seconds", "cpu_seconds"} else _without_times(item)
            for key, item in value.items()
        }
    if isinstance(value, list):
        return [_without_times(item) for item in value]
    return value


def test_perdomo_table3_nacl_workflow() -> None:
    model = _perdomo_table3_model()
    result = _equilibrium._solve_tp_flash(
        _native_transport(model),
        298.15,
        2508.0,
        PERDOMO_TABLE3_FEED,
        model.parameter_fingerprint,
    )

    expected_feed = tuple(value / sum(_source_amounts) for value in _source_amounts)
    assert PERDOMO_TABLE3_FEED == pytest.approx(expected_feed, abs=1.0e-15)
    assert sum(PERDOMO_TABLE3_FEED) == pytest.approx(1.0, abs=1.0e-15)
    assert PERDOMO_TABLE3_FEED[1] - PERDOMO_TABLE3_FEED[2] == pytest.approx(0.0, abs=1.0e-15)
    assert len(result["phases"]) == 1
    assert result["phases"][0]["phase_fraction"] == pytest.approx(1.0)
    reference = result["step2"]["reference_pressure_envelope"]
    assert [root["mechanical_class"] for root in reference["roots"]] == [
        "strict_stable",
        "unstable",
        "strict_stable",
    ]
    assert reference["selected_root_index"] == 2
    selected = reference["roots"][reference["selected_root_index"]]
    assert selected["pressure_residual"] == pytest.approx(0.0, abs=1.0e-8)
    assert selected["objective"] == min(
        root["objective"]
        for root in reference["roots"]
        if root["mechanical_class"] == "strict_stable"
    )

    assert result["step2"]["minimum_tpd"] >= -1.0e-8


@pytest.mark.parametrize(
    (
        "pressure_pa",
        "molality_mol_per_kg",
        "expected_free_energy",
        "expected_phase_fractions",
        "expected_mole_fractions",
        "expected_phase_volumes",
    ),
    (
        pytest.param(
            3181.454397,
            0.005000000139214637,
            0.21155648678431022,
            (0.9851851723141468, 0.014814827685853224),
            (
                (0.9999403402467356, 2.9829876632190147e-05, 2.9829876632190147e-05),
                (0.991809243149451, 0.004095378425274544, 0.004095378425274544),
            ),
            (1.7802235178776976e-05, 0.01151173970113826),
            id="figure1a-0.005molal-below-boundary",
        ),
        pytest.param(
            3213.5903,
            0.005000000139214637,
            0.21162822804854456,
            (0.9976330842742641, 0.0023669157257357933),
            (
                (0.9998598441728949, 7.007791355255787e-05, 7.007791355255787e-05),
                (0.9829751159939637, 0.008512442003018143, 0.008512442003018143),
            ),
            (1.8026580233950744e-05, 0.0018208010582469352),
            id="figure1a-0.005molal-source-boundary",
        ),
        pytest.param(
            2508.0,
            5.6,
            -24.27753305522489,
            (1.0,),
            ((0.8321050353538131, 0.08394748232309347, 0.08394748232309347),),
            (0.9849669198212629,),
            id="table3-5.6molal",
        ),
    ),
)
def test_perdomo_nacl_results_remain_in_the_numerical_regression_region(
    pressure_pa: float,
    molality_mol_per_kg: float,
    expected_free_energy: float,
    expected_phase_fractions: tuple[float, ...],
    expected_mole_fractions: tuple[tuple[float, ...], ...],
    expected_phase_volumes: tuple[float, ...],
) -> None:
    feed = _nacl_feed(molality_mol_per_kg)
    result = equilibrium.tp_flash(
        _perdomo_table3_model(),
        298.15 * epcsaft.unit_registry.kelvin,
        pressure_pa * epcsaft.unit_registry.pascal,
        feed,
    )
    assert result.overall_mole_fractions == pytest.approx(
        feed,
        abs=NUMERICAL_ATOL,
        rel=NUMERICAL_RTOL,
    )
    assert result.total_free_energy_over_rt == pytest.approx(
        expected_free_energy,
        abs=NUMERICAL_ATOL,
        rel=NUMERICAL_RTOL,
    )
    assert result.phase_fractions == pytest.approx(
        expected_phase_fractions,
        abs=PHASE_FRACTION_ATOL,
        rel=NUMERICAL_RTOL,
    )
    assert tuple(phase.volume_m3 for phase in result.phases) == pytest.approx(
        expected_phase_volumes,
        abs=PHASE_VOLUME_ATOL,
        rel=NUMERICAL_RTOL,
    )
    assert len(result.phases) == len(expected_mole_fractions)
    for phase, expected in zip(result.phases, expected_mole_fractions, strict=True):
        assert phase.mole_fractions == pytest.approx(
            expected,
            abs=PHASE_COMPOSITION_ATOL,
            rel=NUMERICAL_RTOL,
        )
