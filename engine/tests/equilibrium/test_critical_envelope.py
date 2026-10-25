"""Installed-wheel transport of critical and envelope declarations, against independent teqp and FeOs values."""
import math
import pytest
from epcsaft import Mixture, Parameters, Pure
from epcsaft.equilibrium import (Amounts, CriticalDeclaration, EnvelopeCondition, EnvelopeCoordinate, Free, Phase, Pinned,
                                 Problem, compile_problem, solve_equilibrium)
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
