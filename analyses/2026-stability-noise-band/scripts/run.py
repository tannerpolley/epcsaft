"""Scale-aware final stability decision (#159): the reduced-Hessian sign at rounding noise.

Usage: run.py LABEL. LABEL names the installed wheel; every JSON records its SHA-256.
- states: 150 seeded packet-3 MEA liquids at 318.15 K with fixed ln K (the #154/#157 worker's case). Every
  successful state also gets its spectrum, rebuilt in NumPy from the public compiled Jacobian and Hessian:
  raw ||H||_2, lambda_max, condition, eps ||H||, n eps ||S||, the equilibrated minimum/maximum and the band.
- negatives: two identical ethane pseudo-species A <=> B (K = 1) held at the symmetric stationary point, where
  a large k_ij makes the reaction direction a curvature maximum.
- degenerate: the same pair near the k_ij where that curvature vanishes, at the default row tolerance.
"""
from __future__ import annotations

import csv
import hashlib
import importlib.metadata
import json
import os
import random
import sys
from pathlib import Path
from urllib.parse import unquote, urlparse

import numpy as np

import epcsaft
from epcsaft import equilibrium as q

ROOT = Path(__file__).resolve().parents[3]
RESULTS = Path(__file__).resolve().parents[1] / "results"
LABEL = sys.argv[1]
EPS = np.finfo(float).eps
sys.path.insert(0, str(ROOT / "engine"))
from tests.native_inputs import ETHANE, parameters  # noqa: E402

MEA = ("carbon-dioxide", "monoethanolamine", "water", "protonated-monoethanolamine", "carbamate-anion",
       "bicarbonate-anion", "carbonate-anion", "hydronium-cation", "hydroxide-anion")
REACTIONS = ({"water": -2.0, "hydronium-cation": 1.0, "hydroxide-anion": 1.0},
             {"carbon-dioxide": -1.0, "water": -2.0, "bicarbonate-anion": 1.0, "hydronium-cation": 1.0},
             {"water": -1.0, "bicarbonate-anion": -1.0, "carbonate-anion": 1.0, "hydronium-cation": 1.0},
             {"monoethanolamine": 1.0, "water": -1.0, "carbamate-anion": -1.0, "bicarbonate-anion": 1.0},
             {"protonated-monoethanolamine": -1.0, "water": -1.0, "monoethanolamine": 1.0, "hydronium-cation": 1.0})
LN_K = (-605.348, -587.8, -948.2312, -25.1014, -298.6415)  # the absorber's 318.15 K density-basis values


def wheel() -> str:
    url = json.loads(importlib.metadata.distribution("epcsaft").read_text("direct_url.json"))["url"]
    return hashlib.sha256(Path(unquote(urlparse(url).path)).read_bytes()).hexdigest()


def mea(pressure: float, loading: float):
    packet = ROOT / "data/packets/mea-co2-h2o-nine-species-estimation/3/parameters/parameters.json"
    mapping = json.loads(packet.read_text())
    # Packet 3 has no CO2/MEA cross edge; declared to keep that model (#569), not a physical-absence claim.
    mapping["topology"]["no_cross_association"] = [{
        "component_id_a": "carbon-dioxide", "component_id_b": "monoethanolamine",
        "provenance": {"source_id": mapping["document_id"], "locator": "packet 3: no CO2/MEA edge"}}]
    model = epcsaft.Mixture(epcsaft.Parameters.from_mapping(mapping, components=MEA))
    reactions = [q.Reaction(nu, correlation=q.ReactionLogPolynomial(a, 0.0, 0.0, 0.0, 318.15, 250.0, 450.0,
                                                                    q.EOS_STANDARD_STATE_ID))
                 for nu, a in zip(REACTIONS, LN_K)]
    feed = {"carbon-dioxide": loading * 9.7292821, "monoethanolamine": 9.7292821, "water": 76.9751925}
    return model, q.Problem(phases=[q.Phase("liquid", kind="liquid")], T=318.15, P=pressure,
                            feed=q.Amounts(feed), reactions=reactions)


def pair(k_ij: float, temperature: float):
    model = epcsaft.Mixture(parameters((ETHANE, ETHANE), [("k_ij", (0, 1), k_ij)]))
    return model, q.Problem(phases=[q.Phase("L", kind="liquid")], T=temperature, P=5e6,
                            feed=q.Amounts({"s0": 0.5, "s1": 0.5}),
                            reactions=[q.Reaction({"s0": -1.0, "s1": 1.0}, equilibrium_constant=1.0)])


def reduce(h, a):
    _, s, vt = np.linalg.svd(a)
    rank = int((s > s[0] * max(a.shape) * EPS).sum())
    z = vt[rank:].T
    return np.linalg.eigvalsh(z.T @ h @ z), s[0] / s[rank - 1]


def spectrum(model, problem, result) -> dict:
    """Independent NumPy rebuild of the certificate at the returned point from public callbacks."""
    compiled = q.compile_problem(model, problem)
    n = compiled.primal_size
    point, weights = list(result.coordinates[:n]), list(result.constraint_multipliers)
    h = np.array(compiled.hessian(point, weights)).reshape(n, n)
    objective = np.array(compiled.hessian(point, [])).reshape(n, n)
    split = np.abs(objective) + np.abs(h - objective)
    rows = [i for i, row in enumerate(result.rows) if str(row.kind).split(".")[-1] not in ("Chemical", "Pressure")]
    a = np.array(compiled.jacobian(point)).reshape(-1, n)[rows]
    raw, _ = reduce(h, a)
    d = 1.0 / np.sqrt(np.where(np.diag(split) > 0.0, np.diag(split), 1.0))
    scaled = d[:, None] * h * d
    equilibrated, kappa = reduce(scaled, a * d)
    band = EPS * (n * np.linalg.norm(d[:, None] * split * d) + 2.0 * kappa * np.linalg.norm(scaled))
    norm = np.linalg.norm(h, 2)
    return {"raw_min": raw[0], "raw_max": raw[-1], "norm_H": norm, "condition": abs(raw[-1] / raw[0]),
            "eps_norm_H": EPS * norm, "n_eps_norm_S": n * EPS * np.linalg.norm(split, 2),
            "equilibrated_min": equilibrated[0], "equilibrated_max": equilibrated[-1], "band": band,
            "margin": equilibrated[0] / band, "x_co2": result.mole_fractions[0]}


def outcome(result) -> str:
    if result.success:
        return "stable"
    return "indeterminate" if "indeterminate" in result.message else (
        "not_positive_definite" if "not positive definite" in result.message else "other")


def record(model, problem, **inputs) -> dict:
    result = q.solve_equilibrium(model, problem)
    row = {**inputs, "outcome": outcome(result), "stop": str(result.convergence.stop).split(".")[-1],
           "reduced_hessian_minimum": result.reduced_hessian_minimum, "residual": result.convergence.residual,
           "message": "" if result.success else result.message[:160]}
    return {**row, **spectrum(model, problem, result)} if result.success else row


def write(name: str, rows: list[dict], summary: dict) -> None:
    RESULTS.mkdir(exist_ok=True)
    with (RESULTS / f"{name}-{LABEL}.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(dict.fromkeys(k for r in rows for k in r)), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    provenance = {"wheel_sha256": wheel(), "label": LABEL, "omp_num_threads": os.environ.get("OMP_NUM_THREADS")}
    (RESULTS / f"{name}-{LABEL}.json").write_text(json.dumps({**provenance, **summary}, indent=1) + "\n")


def main() -> None:
    random.seed(0)
    draws = [(random.uniform(1.05e5, 1.12e5), random.uniform(0.03, 0.56)) for _ in range(150)]
    states = [record(*mea(p, loading), P_Pa=p, loading=loading) for p, loading in draws]
    repeat = [record(*mea(p, loading), P_Pa=p, loading=loading)["outcome"] for p, loading in draws[:20]]
    outcomes = [r["outcome"] for r in states]
    margins = [r["margin"] for r in states if "margin" in r]
    write("states", states, {"states": len(states), "outcomes": {o: outcomes.count(o) for o in set(outcomes)},
                             "repeat_first_20_identical": repeat == outcomes[:20],
                             "minimum_equilibrated_margin": min(margins) if margins else None})
    negatives = [record(*pair(k, t), k_ij=k, T_K=t) for k, t in
                 ((0.0, 200.0), (0.1, 200.0), (0.1, 150.0), (0.2, 150.0), (0.2, 200.0), (0.3, 200.0), (0.5, 150.0))]
    degenerate = [record(*pair(k, 200.0), k_ij=k, T_K=200.0) for k in
                  (0.12810288983397294, 0.12810288983397303, 0.12810288983397305, 0.1281028898339731)]
    write("pair", negatives + degenerate, {"rows": len(negatives) + len(degenerate)})


def declarations():
    return (mea(1e5, .2), pair(0., 200.))


if __name__ == "__main__":
    main()
