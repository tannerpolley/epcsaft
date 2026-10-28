"""Manufactured MEA packet: exact consuming actions against re-solved lower actions."""

import gzip
import json
from pathlib import Path
from typing import Any

import pytest

import epcsaft
from epcsaft import _core as c


def decode(v: Any) -> Any:
    if isinstance(v, list):
        return [decode(x) for x in v]
    if not isinstance(v, dict):
        return v
    if "enum" in v:
        return getattr(getattr(c, v["enum"]), v["name"])
    cls = getattr(c, v["type"])
    fields = {k: decode(x) for k, x in v["fields"].items()}
    try:
        obj = cls()
    except TypeError:
        return cls(**fields)
    for k, x in fields.items():
        setattr(obj, k, x)
    return obj


def packet(mode: str) -> Any:
    data = json.loads(
        gzip.decompress((Path(__file__).resolve().parents[3] / "tools/fixtures/mea-cp.json.gz").read_bytes())
    )
    doc = data["parameters"]
    f = next(f for f in doc["model_families"] if f["kind"] == "permittivity")
    f.update(choice="solvent-mole-weighting", solvents=["water", "monoethanolamine"], normalization=mode)
    for component in doc["components"]:
        component["coefficients"] = [
            x for x in component["coefficients"] if x["family"] not in ("dielectric_volume", "dielectric_temperature")
        ]
    for name, e, cs in [("water", 80, [-0.2, 0.001, 0.00001]), ("monoethanolamine", 20, [-0.05, 0.0005, -0.000002])]:
        doc["correlations"].append(
            {
                "correlation_id": "test/" + name,
                "component_id": name,
                "family": "relative_permittivity",
                "form": "constant-plus-polynomial-in-temperature",
                "constant": {"value": {"magnitude": e, "unit": "dimensionless"}},
                "reference_temperature": {"magnitude": 353.15, "unit": "kelvin"},
                "terms": [
                    {
                        "power": k,
                        "coefficient": {
                            "identity": f"test/{name}/{k}",
                            "value": {"magnitude": v, "unit": ("1 / kelvin" if k == 1 else f"1 / kelvin ** {k}")},
                        },
                    }
                    for k, v in enumerate(cs, 1)
                ],
            }
        )
    p = epcsaft.Parameters.from_mapping(doc)
    m = c.HelmholtzModel(p._native, [], decode(data["thermochemistry"]))
    return data, m, decode(data["observations"][0]["fields"]["problem"])


D = c.SolvedStateActionDirection
Observable = c.SolvedStateObservable
K = c.SolvedStateObservableKind
os = [
    Observable(K.PhaseMolarDensity, 0),
    *[Observable(K.PhaseComponentAmount, 0, i) for i in [2, 1, 3]],
    Observable(K.PhaseResidualChemicalPotentialOverRt, 0, 3),
    Observable(K.PhaseProperty, 0, 0, c.PropertyObservable.TotalEnthalpy),
    Observable(K.SystemEnthalpy),
]
names = ["rho", "water", "amine", "MEAH", "mu_res_RT", "h", "H"]
atols = [1e-6, 1e-10, 1e-10, 1e-10, 1e-11, 1e-5, 1e-5]
T = D(1, 0, [0] * 9, [])
F = D(0, 0, [0, 1, -1, 0, 0, 0, 0, 0, 0], [])


def action(m: Any, p: Any, ds: Any) -> Any:
    r = c.solve_equilibrium(m, p)
    assert r.success and max(abs(v) for v in r.classified_residuals) <= p.row_tolerance, r.message
    assert all(v > 0 for v in r.mole_fractions)
    state = m.state(p.temperature.value, r.molar_densities[0], r.mole_fractions)
    assert state.pressure_density_derivative > 0 and abs(state.pressure - p.pressure.value) < 1e-3
    items = c.solved_state_actions(
        c.compile_problem(m, p), r, c.SolvedStateActionRequest(os, ds or [D(0, 0, [0] * 9, [])])
    ).results
    assert all(x.status == c.SolvedStateActionStatus.Available for x in items), [
        (x.status, x.diagnostic_message) for x in items
    ]
    return [x.action if ds else x.value for x in items]


@pytest.mark.parametrize("mode", ["unnormalized", "solvent-normalized"])
def test_solvent_mole_weighting_actions(mode: str) -> None:
    data, m, p = packet(mode)
    for ds in [[T], [F], [T, T], [T, F], [F, F], [F, F, F]]:
        exact = action(m, p, ds)
        d = ds[-1]
        rows = []
        done = [False] * 7
        for rung in range(9):
            h = (0.2 if d.temperature_k else 0.002) / 2**rung
            vals = []
            for s in [1, -1]:
                pp = decode(data["observations"][0]["fields"]["problem"])
                t = pp.temperature
                t.value = 353.15 + s * h * d.temperature_k
                pp.temperature = t
                feed = pp.feed
                feed.values = [v + s * h * dv for v, dv in zip(p.feed.values, d.feed_amounts_mol)]
                pp.feed = feed
                vals.append(action(m, pp, ds[:-1]))
            rows.append([(a - b) / (2 * h) for a, b in zip(*vals)])
            if rung == 0:
                continue
            for i in range(7):
                if done[i]:
                    continue
                a, b = rows[-2][i], rows[-1][i]
                rich = (4 * b - a) / 3
                refine = abs(b - a) / 3
                target = atols[i] + 1e-6 * max(abs(rich), abs(exact[i]))
                err = abs(exact[i] - rich)
                done[i] = refine <= target / 4
                if done[i]:
                    assert err <= target, (mode, names[i], rich, exact[i], target)
                print(
                    json.dumps(
                        {
                            "mode": mode,
                            "directions": [x.temperature_k for x in ds],
                            "quantity": names[i],
                            "h": h,
                            "coarse": a,
                            "fine": b,
                            "rich": rich,
                            "exact": exact[i],
                            "refinement": refine,
                            "target": target,
                            "error": err,
                            "resolved": done[i],
                            "passed": done[i] and err <= target,
                        }
                    ),
                    flush=True,
                )
            if all(done):
                break
        assert all(done), (mode, ds, done)


@pytest.mark.parametrize("mode", ["unnormalized", "solvent-normalized"])
def test_solvent_mole_weighting_fit_workers(mode: str) -> None:
    data, _, _ = packet(mode)
    params = epcsaft.Parameters.from_mapping(data["parameters"])
    observations = decode(data["observations"]) + decode(data["observations"])
    changed = observations[1].problem
    temperature = changed.temperature
    temperature.value += 0.1
    changed.temperature = temperature
    observations[1].problem = changed
    coords = decode(data["coordinates"])
    thermal = decode(data["thermochemistry"])
    results: list[Any] = [
        c.fit(
            params._native,
            thermal,
            coords,
            observations,
            data["weights"] * 2,
            c.FitControls(maximum_iterations=1, maximum_elapsed_time_seconds=120, gradient_tolerance=1e12, workers=w),
        )
        for w in [1, 2]
    ]
    assert isinstance(params._document, c.ParameterDocument)
    for r in results:
        assert r.usable, r.message
        assert all(s == c.ObservationStatus.Available for s in r.statuses)
        rebuilt = c.fitted_document(params._document, r.parameters)
        assert [(e.field, e.index, e.magnitude) for e in rebuilt.entries if e.field.startswith("permittivity_")] == [
            (e.field, e.index, e.magnitude) for e in params._document.entries if e.field.startswith("permittivity_")
        ]
    for name in ["residuals", "optimizer_jacobian", "optimizer"]:
        a, b = (list(getattr(r, name)) for r in results)
        assert a == pytest.approx(b, rel=1e-10, abs=1e-12)
        print(mode, name, a, flush=True)
    assert any(abs(v) > 1e-8 for v in results[0].optimizer_jacobian)


def test_ion_free_nonsolvent_state() -> None:
    # Preserve the empty electrostatic sums without assigning a dielectric quotient.
    _, model, _ = packet("solvent-normalized")
    x = [1.0] + [0.0] * 8
    state = model.state(353.15, 1000.0, x)
    assert abs(state.pressure - 2687360.5777911283) <= 1e-6
    assert state.bulk_relative_permittivity == 0.0
    model.state(353.15, 1000.0, [0.0, 0.0, 1.0] + [0.0] * 6)
    assert model.state(353.15, 1000.0, x).pressure == state.pressure
    with pytest.raises(ValueError, match="solvent denominator"):
        model.state(353.15, 1000.0, [0.0] * 3 + [0.5, 0.5] + [0.0] * 4)
    with pytest.raises(ValueError, match="solvent denominator"):
        c.infinite_dilution_reference(model, 353.15, state.pressure, x, (3, 4), (1, 1))
