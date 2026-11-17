"""Installed-wheel transport of critical and envelope declarations, against independent teqp and FeOs values."""
import math
from typing import Any

import pytest
from epcsaft import Mixture, Parameters, Pure, tolerances
from epcsaft.equilibrium import (Amounts, critical_point, CriticalDeclaration, EnvelopeCondition, EnvelopeCoordinate, Free, Phase, Pinned,
                                 Problem, compile_problem, solve_equilibrium)
from tests.native_inputs import same
from tests.packets import parameters_json


def test_critical_candidate_matches_teqp() -> None:
    model = Pure(Parameters.from_json(parameters_json("gross-2001-propane", 3), components=("propane",)))
    declaration = CriticalDeclaration()
    declaration.density_guess_mol_per_m3 = 5000.0
    problem = Problem(T=Free(375.0), P=Free(1.0e6), feed=Amounts({"propane": 1.0}), phases=[Phase("F", kind="liquid")],
                      row_tolerance=1e-11, critical=declaration)
    result = solve_equilibrium(model, problem)
    # teqp v0.23.2 solve_pure_critical on the same packet (analysis issue-224-sliceB critical-pure.csv).
    assert result.temperature == pytest.approx(375.14002748147936, rel=1e-8), result.message
    assert result.critical is not None and abs(result.critical.eigenvalues[0]) <= result.critical.zero_band
    candidate = critical_point(model, T_guess=375.0, P_guess=1.0e6, rho_guess=5000.0)
    assert candidate.temperature == pytest.approx(result.temperature, rel=1e-8), candidate.message
    assert candidate.critical is not None and abs(candidate.critical.eigenvalues[0]) <= candidate.critical.zero_band


def test_associating_critical_candidate_after_failed_declared_start() -> None:
    # Retained numerical baseline: the same packet's critical equations solved from a declared 4*epsilon/k seed.
    model = Pure(Parameters.from_json(parameters_json("ascani-2022-case-study-2-electrolyte-lle", 3), components=("water",)))
    declaration = CriticalDeclaration()
    declaration.density_guess_mol_per_m3 = 9960.986137627624
    problem = Problem(T=Free(503.19124775), P=Free(1e6), feed=Amounts({"water": 1.0}),
                      phases=[Phase("F", kind="liquid")], critical=declaration)
    result = solve_equilibrium(model, problem)
    assert result.success, result.message
    assert result.temperature == pytest.approx(699.9801552748312, rel=1e-10)
    assert result.pressure == pytest.approx(38835718.7905412, rel=1e-10)
    assert result.critical is not None and result.critical.rank == result.unknown_count
    assert abs(result.critical.eigenvalues[0]) <= result.critical.zero_band
    assert max(map(abs, result.classified_residuals)) <= 1e-10
    search = result.critical_search
    assert search is not None and search.selected is not None
    cap = {t.name: t.value for t in tolerances()}
    assert all(i <= cap["EQUILIBRIUM_RECOVERY_WORK"] and e <= cap["EQUILIBRIUM_RECOVERY_EVALUATIONS"]
               and h <= cap["EQUILIBRIUM_RECOVERY_HESSIAN_CALLS"] for i, e, h in search.added_work)
    assert all(r.numerical_estimates is not None and r.numerical_estimates.cause == ""
               for r in search.results if r.success)


def test_envelope_condition_matches_feos_bubble() -> None:
    model = Mixture(Parameters.from_json(parameters_json("gross-2001-methane-ethane", 3), components=("methane", "ethane")))
    condition = EnvelopeCondition()
    condition.kind, condition.target = EnvelopeCoordinate.LogTemperature, math.log(180.0)
    problem = Problem(T=Free(180.0), P=Free(1.0e5), feed=Amounts({"methane": 0.5, "ethane": 0.5}), row_tolerance=1e-11,
                      phases=[Phase("L", kind="liquid"), Phase("V", kind="vapor", amount=Pinned(0.0))], envelope=condition)
    result = solve_equilibrium(model, problem)
    # FeOs v0.10.1 bubble_point at 180 K (analysis issue-224-sliceB envelope-points.csv).
    assert result.pressure == pytest.approx(1596859.111206796, rel=1e-7), result.message


def test_azeotrope_declaration_reaches_the_native_rows() -> None:
    model = Mixture(Parameters.from_json(parameters_json("gross-2001-methane-ethane", 3), components=("methane", "ethane")))
    phases = [Phase("L", kind="liquid"), Phase("V", kind="vapor", amount=Pinned(0.0))]
    feed = Amounts({"methane": 0.5, "ethane": 0.5})
    rows = [row.eqid for row in compile_problem(model, Problem(T=200.0, P=Free(), feed=feed, phases=phases, azeotrope=True)).rows]
    # One total-amount balance and one equal-composition row replace the two species balances of the bubble cell.
    assert rows.count("material_balance") == 1 and rows[-1] == "azeotrope_condition", rows
    with pytest.raises(ValueError, match="azeotrope cell requires"):
        compile_problem(model, Problem(T=Free(200.0), P=Free(), feed=feed, phases=phases, azeotrope=True))


@pytest.mark.parametrize("temperature", [300.0, 375.139])
def test_saturation_conditioning_transport(temperature: float) -> None:
    """A retained double warning remains distinct from the returned wider estimate through Python."""
    model = Pure(Parameters.from_json(parameters_json("gross-2001-propane", 3), components=("propane",)))
    problem = Problem(T=temperature, P=Free(1e5), feed=Amounts({"propane": 1}),
                      phases=[Phase("L", kind="liquid"), Phase("V", kind="vapor", amount=Pinned(0))])
    result = solve_equilibrium(model, problem)
    assert result.success, result.message
    record = next(r for r in result.rank_records if r.basis == "log_density")
    assert record.observation.point_kind == "endpoint"
    warning = record.first_conditioning_warning
    if temperature == 300:
        assert warning is None
    else:
        assert warning is not None
        assert warning.precision_bits == 53 and warning.error >= 9.5e-9
        assert warning.point_kind == "last_linearization" and len(warning.coordinates) == 2
        assert warning.controls[0] == temperature
        assert record.observation.precision_bits > 53
        assert record.observation.error < warning.criterion_value


def test_polish_endpoint_conditioning() -> None:
    import importlib.util
    import sys
    from pathlib import Path

    path = Path(__file__).resolve().parents[3] / "engine/tests/cases/tpd-phase-discovery/scripts/baseline.py"
    spec = importlib.util.spec_from_file_location("conditioning_baseline", path)
    assert spec is not None and spec.loader is not None
    baseline = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = baseline
    spec.loader.exec_module(baseline)
    from epcsaft import equilibrium as eq
    for name, model, problem in baseline.declared_cases():
        if name not in ("case-1/two-liquid", "case-6/figure-2-point-7"):
            continue
        result = eq.solve_equilibrium(model, problem)
        assert result.success
        observation = next(rank.observation for rank in result.rank_records if rank.basis == "numerical")
        print(name, "available", observation.available, "condition", observation.condition,
              "point", observation.point_kind, "cause", observation.cause, "error", observation.error)
        if name == "case-6/figure-2-point-7":
            assert observation.available and observation.condition > 0
            assert observation.point_kind == "endpoint"
            assert observation.cause == "inner_floor_unavailable"
        else:
            assert not observation.available
            assert observation.point_kind == "last_linearization"
        if name == "case-6/figure-2-point-7":
            break


def test_diffusive_spinodal_transport(evidence: Any) -> None:
    from epcsaft.equilibrium import SpinodalDeclaration

    model = Mixture(Parameters.from_json(parameters_json("gross-2001-methane-ethane", 3), components=("methane", "ethane")))
    declaration = SpinodalDeclaration()
    declaration.density_guess_mol_per_m3 = 4000.0
    declaration.density_nodes = [100.0 * 200.0 ** (k / 80) for k in range(81)]
    result = solve_equilibrium(model, Problem(T=250.0, P=Free(1e6), feed=Amounts({"methane": 0.5, "ethane": 0.5}),
        phases=[Phase("F", kind="liquid")], spinodal=declaration, row_tolerance=1e-11))
    assert result.success and result.spinodal is not None, result.message
    assert result.critical is None
    same(evidence, "spinodal/1/0/rho", result.molar_densities[0])
    assert len(result.spinodal.search.brackets) == 2
    assert result.spinodal.normalized_null_residual <= 1e-10
    assert [r.eqid for r in result.rows].count("diffusive_spinodal") == 3
