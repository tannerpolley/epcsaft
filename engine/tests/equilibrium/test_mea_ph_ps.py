"""Reacting MEA PH/PS on the frozen manufactured packet (tools/fixtures/mea-cp.json.gz, liquid, five reactions, six ions
without records) at 1 bar. Its MEA record is anchored at h = s = 0 at 313.15 K, so H and S are anchor-relative. H and S
of the driven TP state recover it from a start 4 K above, inside every record and reaction-law interval; dT/dH and
dT/dS match 1/(N_f c_P,eq) and T/(N_f c_P,eq) from a separately compiled driven-TP apparent Cp."""

import gzip
import json
from pathlib import Path
from typing import Any

import pytest

import epcsaft
from epcsaft import _core as c
from tests.equilibrium.test_solvent_mole_weighting import decode

DATA = json.loads(gzip.decompress((Path(__file__).resolve().parents[3] / "tools/fixtures/mea-cp.json.gz").read_bytes()))
MODEL = c.HelmholtzModel(epcsaft.Parameters.from_mapping(DATA["parameters"])._native, [], decode(DATA["thermochemistry"]))
K = c.SolvedStateObservableKind
Observable = c.SolvedStateObservable


def problem(T: Any) -> Any:
    p = decode(DATA["observations"][0]["fields"]["problem"])
    p.temperature = T if isinstance(T, c.FreeValue) else c.FreeValue.fixed(T)
    pressure = p.pressure
    pressure.value = 1e5
    p.pressure = pressure
    return p


def read(p: Any, r: Any, observables: Any, direction: Any = None) -> Any:
    d = direction or c.SolvedStateActionDirection(0, 0, [0] * 9, [])
    items = c.solved_state_actions(c.compile_problem(MODEL, p), r, c.SolvedStateActionRequest(observables, [d])).results
    assert all(x.status == c.SolvedStateActionStatus.Available for x in items), [x.diagnostic_message for x in items]
    return items


@pytest.mark.parametrize("T", [298.15, 318.15])
@pytest.mark.parametrize("h", [True, False])
def test_mea_ph_ps(T: float, h: bool) -> None:
    tp = problem(T)
    ref = c.solve_equilibrium(MODEL, tp)
    assert ref.success, ref.message
    H, s = (x.value for x in read(tp, ref, [Observable(K.SystemEnthalpy), Observable(K.PhaseProperty, 0, 0, c.PropertyObservable.TotalEntropy)]))
    q = problem(c.FreeValue.variable(T + 4.0))
    if h:
        q.enthalpy_j = H
    else:
        q.entropy_j_per_k = ref.amounts[0] * s
    r = c.solve_equilibrium(MODEL, q)
    assert r.success and r.temperature is not None, r.message
    assert r.temperature == pytest.approx(T, rel=1e-8)
    assert r.amounts == pytest.approx(ref.amounts, rel=1e-8)
    assert r.mole_fractions == pytest.approx(ref.mole_fractions, rel=1e-8)
    assert max(map(abs, r.classified_residuals)) <= 1e-10
    # Density-basis affinity A_r/RT = sum_i nu_ri (chemical row i): the moiety and charge multipliers cancel on nu.
    rows = c.compile_problem(MODEL, q).rows
    for reaction in q.reactions:
        affinity = sum(reaction.stoichiometry[row.index] * value
                       for row, value in zip(rows, r.residuals, strict=True) if row.kind == c.RowKind.Chemical)
        assert abs(affinity) <= 1e-10
    driven = problem(r.temperature)
    cp = read(driven, c.solve_equilibrium(MODEL, driven), [Observable(K.EquilibriumIsobaricHeatCapacity)])[0].value
    feed = sum(q.feed.values)
    target = (c.SolvedStateActionDirection(0, 0, [0] * 9, [], enthalpy_j=1.0) if h
              else c.SolvedStateActionDirection(0, 0, [0] * 9, [], entropy_j_per_k=1.0))
    action, expected = read(q, r, [Observable(K.Temperature)], target)[0].action, (1.0 if h else r.temperature) / (feed * cp)
    print(f"MEA {'PH' if h else 'PS'} {T} K: T error {abs(r.temperature - T) / T:.3g}, rows "
          f"{max(map(abs, r.classified_residuals)):.3g}, dT {action!r}, 1/(N_f c_P,eq) form {expected!r}, error {abs(action / expected - 1):.3g}")
    assert action == pytest.approx(expected, rel=1e-6)
