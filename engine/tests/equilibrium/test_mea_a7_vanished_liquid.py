"""Conserved liquid transfer from the A7 neutral vapor at archived flash inputs."""
import json
import math
from pathlib import Path
from typing import Any

import pytest

from epcsaft import Mixture, Parameters, _core
from epcsaft.equilibrium import Amounts, NeutralReference, Phase, Problem, solve_equilibrium
from tests.native_inputs import MEA, nominally_unloaded_mea

REFERENCE = Path(__file__).resolve().parents[2] / "data" / "reference" / "mea-a7"
DECLARATIONS = json.loads((REFERENCE / "vanished-liquid.json").read_text())
CASES = DECLARATIONS["rows"]
ROW_CRITERION = 1e-10
TOLERANCES = {t.name: t.value for t in _core.tolerances()}


@pytest.fixture(scope="module")
def solve():
    model = Mixture(Parameters.from_json(REFERENCE / "parameters.json"))
    _, reactions, _, _, _ = nominally_unloaded_mea(False)
    correlations = {
        row["reaction_id"]: {c["name"]: c["value"]["magnitude"] for c in row["coefficients"]}
        for row in json.loads((REFERENCE / "parameters.json").read_text())["reaction_correlations"]
    }
    # Use the A7 correlations in the existing reference-temperature and log bases.
    for reaction in reactions:
        law = reaction.correlation
        if reaction.name == "R2":
            row = correlations["R2"]
            law.a, law.b, law.c, law.d = (
                row["a"] + row["c"] * math.log(law.reference_temperature), row["b_k"], row["c"], row["d_per_k"])
        elif reaction.name == "R4":
            law.a, law.b = correlations["R4"]["a"], correlations["R4"]["b_k"]
        elif reaction.name == "R5":
            row = correlations["R5"]
            law.a, law.b, law.d = (-math.log(10) * row["b"], -math.log(10) * row["a_k"], -math.log(10) * row["c_per_k"])
        reaction.correlation = law
    reference = NeutralReference(solvent="water", standard_molality_mol_per_kg=1.,
                                 solute_molality_ratios=dict(zip(MEA[:2] + MEA[3:], (1., 1., 5., 2., 1., 1., 1., 1.))))

    def flash(case, start):
        T, P, *feed = case["inputs"]
        x, y = DECLARATIONS["anchor_guesses"] if start == "seeded" else ([], [])
        return solve_equilibrium(model, Problem(
            T=T, P=P, feed=Amounts(dict(zip(MEA[:3], feed, strict=True))),
            phases=[Phase("L", kind="liquid", support=MEA, composition_guess=x),
                    Phase("V", kind="vapor", support=MEA[:3], composition_guess=y)],
            reactions=reactions, neutral_reference=reference, topology_declared=True, row_tolerance=ROW_CRITERION))

    return flash


def closed_vapor(report, case):
    assert report is not None
    (vapor,) = report.present_phases
    assert vapor["volume_m3_per_mol"] > 0.0
    assert vapor["amount_mol"] == pytest.approx(sum(case["inputs"][2:]), rel=ROW_CRITERION)
    assert abs(vapor["pressure_row"]) <= ROW_CRITERION
    assert abs(vapor["charge_residual"]) <= TOLERANCES["COMPOSITION_SUM_TOLERANCE"]
    assert vapor["reservoir_inventory_error_mol"] <= TOLERANCES["ADDITION_MATERIAL_BALANCE"] * vapor["amount_mol"]
    assert vapor["reservoir_plane_residual_over_rt"] <= ROW_CRITERION
    return vapor


def liquid_trials(report):
    return [
        trial for trial in report.trials if trial.quantities.get("liquid_root_residual", math.inf)
        <= TOLERANCES["DENSITY_ROOT_RESIDUAL"] and trial.quantities["liquid_root_log_displacement"]
        <= TOLERANCES["DENSITY_ANCHOR_LOG_WINDOW"]
    ]


@pytest.mark.parametrize("row", [row for row, case in CASES.items() if "split" in case])
def test_a7_vanished_liquid_split(row: str, solve, evidence: Any):
    case = CASES[row]
    reference = case["split"]
    states = {start: solve(case, start) for start in ("cold", "seeded")}
    for start, result in states.items():
        assert result.success and result.requested_tolerance_met, (row, start, result.message)
        assert result.convergence.residual <= ROW_CRITERION
    (losing_start,) = case["vanished"]
    recovered = states[losing_start]
    sibling = states["seeded" if losing_start == "cold" else "cold"]
    for name, actual, independent in (
        ("amounts_mol", recovered.amounts, sibling.amounts),
        ("densities_mol_m3", recovered.molar_densities, sibling.molar_densities),
        ("x", recovered.mole_fractions, sibling.mole_fractions),
    ):
        assert actual == pytest.approx(independent, rel=ROW_CRITERION, abs=ROW_CRITERION)
        assert independent == pytest.approx(reference[name], rel=ROW_CRITERION, abs=ROW_CRITERION)
    report = recovered.vanished_phase_stability
    vapor = closed_vapor(report, case)
    (witness,) = [t for t in liquid_trials(report) if "restart_liquid_tpd_over_rt" in t.quantities]
    assert witness.resolved_descent and report.descent == "resolved_descent"
    q = witness.quantities
    transfer = witness.split_checks[-1]
    assert transfer["beta_mol"] > 0 and transfer["admitted"] == 1
    assert transfer["moiety_error_mol"] <= TOLERANCES["ADDITION_MATERIAL_BALANCE"] * vapor["amount_mol"]
    assert transfer["delta_g_over_rt"] < -transfer["error_estimate_over_rt"]
    assert q["restart_residual"] <= ROW_CRITERION
    assert q["restart_liquid_amount_mol"] == recovered.amounts[0]
    error = q["restart_error_over_rt_mol"]
    assert q["restart_delta_g_over_rt_mol"] < -error
    assert q["restart_delta_g_over_rt_mol"] == pytest.approx(reference["delta_g_over_rt_mol"], abs=error, rel=0)
    assert recovered.objective == pytest.approx(sibling.objective, abs=error, rel=0)

    # The native reservoir plane is evaluated directly at the admitted liquid, without
    # another trial optimization. The retained secants use the sibling liquid composition.
    plane = q["restart_liquid_tpd_over_rt"]
    fine, coarse = reference["secants_over_rt"]
    discrepancy = abs(fine - coarse)
    # Use the retained finite-transfer energy uncertainty for the two much
    # smaller additions and convert each extensive error to a secant error.
    secant_error = transfer["error_estimate_over_rt"] * (1.0 / 1e-5 + 1.0 / 1e-4)
    allowance = discrepancy + q["descent_error_estimate_over_rt"] + secant_error
    assert plane < 0
    assert abs(plane - fine) <= allowance and abs(plane - coarse) <= allowance
    assert abs(plane - fine) < abs(plane - coarse)
    evidence.observe(f"A7/{row}", recover_from=losing_start, amounts_mol=recovered.amounts,
                     residual=recovered.convergence.residual, plane=plane, secants=[fine, coarse],
                     plane_allowance=allowance, transfer=transfer, restart=q, work=report.work)


@pytest.mark.parametrize(("row", "start"), [
    (row, start) for row, case in CASES.items() if "split" not in case for start in case["vanished"]
])
def test_a7_vanished_liquid_trials(row: str, start: str, solve, evidence: Any):
    case = CASES[row]
    result = solve(case, start)
    report = result.vanished_phase_stability
    vapor = closed_vapor(report, case)
    admitted = liquid_trials(report)
    assert admitted and all(list(t.phase.support) == list(range(len(MEA))) for t in admitted)
    assert report.completed and report.work["trial_invocations"] > 0
    if report.descent == "resolved_descent":
        assert result.success and result.requested_tolerance_met and result.convergence.residual <= ROW_CRITERION
    else:
        assert not result.success and not result.amounts
        assert report.descent in ("no_witness", "undecided")
        if report.descent == "undecided":
            assert report.undecided
            assert all("liquid pressure-root window" in report.trials[i].message for i in report.undecided)
        else:
            assert report.no_descent_in_trial_set
    evidence.observe(f"A7/{row}/{start}", vapor=vapor, outcome=report.descent, work=report.work,
                     trials=[{"seed": t.seed, "outcome": t.evidence, "message": t.message,
                              "quantities": t.quantities, "transfers": t.split_checks} for t in report.trials])
