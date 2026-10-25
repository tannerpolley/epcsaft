"""Referenced mixed-temperature second actions and pressure-moving enthalpy actions (#147).

Usage: run.py MEA_CHECKOUT ABSORBER_CHECKOUT. MEA states use MEA's production request mapping,
adopted parameters and physical ideal-gas records; the absorber's 3C and S2 node liquids use its
exploratory record, reactions and ideal-gas records. The Engine is the installed non-editable wheel.
"""
from __future__ import annotations

import csv
import hashlib
import importlib.metadata
import json
import sys
from pathlib import Path
from time import perf_counter
from urllib.parse import unquote, urlparse

MEA, ABSORBER = (Path(p).resolve() for p in sys.argv[1:3])
sys.path[:0] = [str(MEA / "analyses/mea_parameter_bundle/scripts"), str(MEA / "src"),
                str(ABSORBER / "src/mea_absorption_column/Thermodynamics")]
import epcsaft
import reactive_bundle as absorber
import shared_evaluation as shared
import verify_reference_calorics as calorics
import verify_temperature_reference as vt
from epcsaft import equilibrium as q

RESULTS = Path(__file__).resolve().parents[1] / "results"
Kind, Property = q.SolvedStateObservableKind, epcsaft.PropertyObservable
CO2 = shared.COMPONENT_IDS.index("carbon-dioxide")
RTOL, T_STEPS, P_STEPS, FEED_STEPS = 1e-6, (0.2, 0.1), (2.0e4, 1.0e4), (2e-4, 1e-4)
LIQUID_DATA = ABSORBER / "src/mea_absorption_column/data/epcsaft_datasets/MEA_greenfield_exploratory"
# Absorber liquid inputs (T K, P Pa, CO2, MEA, water mol/s): case 3C bulk bottom node and the S2 node
# (analyses/greenfield_node_qualification/scripts/qualify_nodes.py at 4f2d12e9).
NODES = {"absorber 3C node": (318.15, 110900.0, 2.4323205276203805, 9.7292821, 76.97519250729157),
         "absorber S2 node": (329.3895574786238, 110899.99998043675, 4.04747099164368, 9.7292821, 76.09022399713623)}
MEA_STATES = ("Bottinger2008_state_050", "Bottinger2008_state_059", "vle_obs_0131", "vle_obs_0224")


def wheel() -> tuple[str, str]:
    url = json.loads(importlib.metadata.distribution("epcsaft").read_text("direct_url.json"))["url"]
    path = Path(unquote(urlparse(url).path))
    return path.name, hashlib.sha256(path.read_bytes()).hexdigest()


class Case:
    """One declared state and a re-solve at shifted T or P (warm-started from the central state)."""

    def __init__(self, label, model, build, temperature, pressure, anchor=None):
        self.label, self.model, self.build, self.anchor = label, model, build, anchor
        self.temperature, self.pressure = temperature, pressure
        self.problem, self.result = self.solve(0.0, 0.0)
        self.width = len(shared.COMPONENT_IDS)
        self.free_pressure = len(self.problem.phases) > 1

    def solve(self, d_t, d_p, d_co2=0.0):
        """build(T, P, CO2 feed shift, anchor) -> Problem; the anchor keeps re-solves on the central branch."""
        problem = self.build(self.temperature + d_t, self.pressure + d_p, d_co2, self.anchor)
        result = q.solve_equilibrium(self.model, problem)
        if not result.success:
            raise RuntimeError(f"{self.label} at dT={d_t}, dP={d_p}, dCO2={d_co2}: {result.message}")
        return problem, result

    def actions(self, d_t, d_p, observables, directions, d_co2=0.0):
        problem, result = (self.problem, self.result) if d_t == d_p == d_co2 == 0.0 else self.solve(d_t, d_p, d_co2)
        batch = q.solved_state_actions(q.compile_problem(self.model, problem), result,
                                       q.SolvedStateActionRequest(observables, directions or [self.direction()]))
        out = []
        for item in batch.results:
            value = item.action if directions else item.value
            if item.status.name != "Available" or value is None:
                raise RuntimeError(f"{self.label}: {item.observable} {item.status.name} {item.diagnostic_message}")
            out.append(value)
        return out

    def direction(self, name=None):
        feed = [0.0] * self.width
        if name == "CO2":
            feed[CO2] = 1.0
        return q.SolvedStateActionDirection(float(name == "T"), float(name == "P"), feed, [])


def liquid(case) -> int:
    return next(i for i, phase in enumerate(case.problem.phases) if phase.kind == "liquid")


def state_observables(case):
    phase = liquid(case)
    items = [(f"ln a_L:{c}", q.SolvedStateObservable(Kind.PhaseLogActivity, phase, i))
             for i, c in enumerate(shared.COMPONENT_IDS)]
    items.append(("n_L", q.SolvedStateObservable(Kind.PhaseAmount, phase)))
    if len(case.problem.phases) > 1:
        vapor = 1 - phase
        items += [("P", q.SolvedStateObservable(Kind.Pressure, 0, 0)),
                  ("y:CO2", q.SolvedStateObservable(Kind.PhaseMoleFraction, vapor, CO2))]
    return items


def slope(function, h):
    """Centred difference of a vector function of a step; the second-order backward difference (same h^2
    order) when +h leaves a reaction correlation's declared temperature interval."""
    try:
        return [(a - b) / (2 * h) for a, b in zip(function(h), function(-h))]
    except epcsaft._core.ParameterError as error:
        if "outside reaction correlation" not in str(error):
            raise
        return [(3 * a - 4 * b + c) / (2 * h) for a, b, c in zip(function(0.0), function(-h), function(-2 * h))]


def richardson(function, steps):
    d = [function(h) for h in steps]
    return [(4 * b - a) / 3 for a, b in zip(*d)], [abs(b - a) / 3 for a, b in zip(*d)], d


def record(rows, case, family, label, exact, reference, resolution, criterion=RTOL):
    for (name, _), a, r, e in zip(label, exact, reference, resolution):
        error = abs(a - r) / max(abs(r), 1e-300)
        rows.append({"state": case.label, "family": family, "observable": name, "action": a, "richardson": r,
                     "resolution": e, "relative_error": error, "criterion": criterion,
                     "status": "pass" if error <= criterion else "fail"})


def second_order_rows(case, rows):
    items = state_observables(case)
    observables = [o for _, o in items]
    others = ["CO2", "T"] + ([] if case.free_pressure else ["P"])
    for other in others:
        exact = case.actions(0.0, 0.0, observables, [case.direction("T"), case.direction(other)])
        reference, resolution, _ = richardson(lambda h: slope(
            lambda d: case.actions(d, 0.0, observables, [case.direction(other)]), h), T_STEPS)
        record(rows, case, f"D2[e_T, e_{other}]", items, exact, reference, resolution)


def enthalpy_rows(case, rows):
    phase = liquid(case)
    observables = [q.SolvedStateObservable(Kind.PhaseProperty, phase, 0, Property.TotalEnthalpy),
                   q.SolvedStateObservable(Kind.PhaseAmount, phase)]

    def extensive(shift, direction=None):
        d_t, d_p, d_co2 = shift
        h, n = case.actions(d_t, d_p, observables, [], d_co2)
        if direction is None:
            return h * n
        dh, dn = case.actions(d_t, d_p, observables, [case.direction(direction)], d_co2)
        return n * dh + h * dn

    # A free-pressure cell moves its solved pressure along T and feed; a fixed-P state is driven along P.
    moves = {"T": (T_STEPS, lambda h: (h, 0.0, 0.0)), "CO2": (FEED_STEPS, lambda h: (0.0, 0.0, h))} \
        if case.free_pressure else {"P": (P_STEPS, lambda h: (0.0, h, 0.0))}
    for move, (steps, shift) in moves.items():
        exact = extensive((0.0, 0.0, 0.0), move)
        reference, resolution, _ = richardson(lambda h: slope(lambda d: [extensive(shift(d))], h), steps)
        record(rows, case, f"d(n_L h_L)/d{move}", [("n_L h_L", None)], [exact], reference, resolution)


def mea_cases():
    model = epcsaft.Mixture(shared.load_parameters(), thermochemistry=calorics.record())
    shared.MODEL = model
    reactions = shared._selected_reactions()
    packet = {o["identity"]: o["request"] for o in shared.load_state_packet()["observations"]}
    for identity in MEA_STATES:
        base = shared.corrected_request(packet[identity], reactions)
        temperature, pressure = float(base["temperature"]["value"]), float(base["pressure"].get("value", base["pressure"].get("initial", 0.0)))

        def build(t, p, d_co2, anchor, base=base):
            request = vt.request_at(base, t)
            if request["pressure"]["role"] == "fixed":
                request["pressure"]["value"] = p
            request["reaction_system"]["feed_amounts_mol"][CO2] += d_co2
            return shared._problem_from_request(request, anchor)

        cold = q.solve_equilibrium(model, build(temperature, pressure, 0.0, None))
        anchor = vt.anchor(build(temperature, pressure, 0.0, None), cold)[0] if cold.success else None
        yield Case(identity, model, build, temperature, pressure, anchor)


def absorber_cases():
    record_, mapping = (absorber._json(LIQUID_DATA / name) for name in ("engine-reactions.json", "parameters.json"))
    for family in mapping["model_families"]:
        if family.get("kind") == "electrolyte" and family.get("choice") == "born":
            for key, value in record_["born_runtime_defaults"].items():
                family.setdefault(key, value)
    ids = [c["component_id"] for c in mapping["components"]]
    assert tuple(ids) == tuple(shared.COMPONENT_IDS), "absorber and MEA component order differ"
    model = epcsaft.Mixture(epcsaft.Parameters.from_mapping(mapping), thermochemistry=absorber.ideal_gas_record(
        LIQUID_DATA / "ideal-gas-thermochemistry.json", ids))
    reactions = [q.Reaction(row["stoichiometry"], name=row["name"],
                            correlation=q.ReactionLogPolynomial(*(row["correlation"][k] for k in absorber._CORRELATION_FIELDS)),
                            reference=q.ReactionReference(getattr(q.ReferenceSourceBasis, row["reference"]["source_basis"]),
                                                          reference_pressure_pa=row["reference"]["reference_pressure_pa"]))
                 for row in record_["reactions"]]

    for label, (t, p, *feed) in NODES.items():
        def build(t, p, d_co2, _anchor, feed=feed):
            amounts = dict(zip(ids[:3], [feed[0] + d_co2, *feed[1:]]))
            return q.Problem(phases=[q.Phase("liquid", kind="liquid")], T=t, P=p, feed=q.Amounts(amounts),
                             reactions=reactions, neutral_reference=q.NeutralReference(**record_["neutral_reference"]))

        yield Case(label, model, build, t, p)


def main():
    name, sha = wheel()
    shared.ENGINE_WHEEL_SHA256 = sha
    started, rows, timing = perf_counter(), [], {}
    for case in [*mea_cases(), *absorber_cases()]:
        clock = perf_counter()
        second_order_rows(case, rows)
        enthalpy_rows(case, rows)
        timing[case.label] = perf_counter() - clock
        print(case.label, f"{timing[case.label]:.1f} s", flush=True)
    RESULTS.mkdir(parents=True, exist_ok=True)
    with (RESULTS / "actions.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    worst = {}
    for row in rows:
        key = f"{row['family']}"
        if row["relative_error"] > worst.get(key, {"relative_error": -1})["relative_error"]:
            worst[key] = {k: row[k] for k in ("state", "observable", "relative_error", "resolution", "richardson")}
    summary = {
        "wheel": name, "wheel_sha256": sha, "criterion": RTOL, "steps": {"T_K": T_STEPS, "P_Pa": P_STEPS, "CO2_mol": FEED_STEPS},
        "inputs": {"mea_checkout_commit": _commit(MEA), "absorber_checkout_commit": _commit(ABSORBER),
                   "absorber_parameters_sha256": hashlib.sha256((LIQUID_DATA / "parameters.json").read_bytes()).hexdigest()},
        "counts": {"rows": len(rows), "pass": sum(r["status"] == "pass" for r in rows),
                   "fail": sum(r["status"] == "fail" for r in rows)},
        "worst_by_family": worst, "seconds": timing, "elapsed_s": perf_counter() - started,
        "claim_limits": "Numerical verification on the listed states only; exploratory absorber record 868a5018 and "
                        "MEA 04a9328 adopted parameters with its physical ideal-gas records. Not physical validation.",
    }
    (RESULTS / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps(summary["counts"]), json.dumps(worst, indent=1))


def _commit(path: Path) -> str:
    head = (path / ".git").read_text().split()[-1] if (path / ".git").is_file() else ""
    import subprocess
    return subprocess.run(["git", "-C", str(path), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip() or head


if __name__ == "__main__":
    main()
