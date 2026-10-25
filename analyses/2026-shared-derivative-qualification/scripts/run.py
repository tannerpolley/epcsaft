"""Issue 26: shared derivative, density-closure and solvent-density qualification.

Run from the repository root with the installed wheel:
    .venv/bin/python analyses/2026-shared-derivative-qualification/scripts/run.py
"""

import csv
import hashlib
import itertools
import json
import math
import subprocess
import sys
from pathlib import Path

from epcsaft import (DensityActionDirection, Mixture, Parameters, PressureActionDirection,
                     PropertyActionStatus, PropertyObservable, Pure, tolerances)

ROOT = Path(__file__).resolve().parents[3]
LEAF = ROOT / "analyses/2026-shared-derivative-qualification"
PACKETS = ROOT / "data/packets"
FINE = (1e-3, 3e-3, 1e-2)
TOLERANCE = (5e-6, 5e-5, 5e-4)  # P4 action criteria, applied relative to the term scale


def packet(packet_id: str, version: int) -> Path:
    return PACKETS / packet_id / str(version) / "parameters/parameters.json"


def molar_mass(path: Path, component: str) -> float:
    record = next(c for c in json.loads(path.read_text())["components"] if c["component_id"] == component)
    value = record["fixed"]["molar_mass"]["value"]
    return value["magnitude"] * {"kilogram / mole": 1.0, "gram / mole": 1e-3}[value["unit"]]


# Real selections at their intended liquid condition; x lists follow `components`.
# Directions: T [K], P [Pa], charge-neutral x.
SELECTIONS = [
    ("figiel-born-ion-fraction", "figiel-yu-held-2025", 3, ("water", "sodium-cation", "chloride-anion"),
     298.15, 1e5, (0.96, 0.02, 0.02),
     ((-0.02, 0.01, 0.01), (0.01, -0.005, -0.005), (0.004, -0.002, -0.002))),
    ("held-dh-solvent-only", "held-cameretti-sadowski-2008", 3, ("water", "sodium-cation", "chloride-anion"),
     298.15, 1e5, (0.96, 0.02, 0.02),
     ((-0.02, 0.01, 0.01), (0.01, -0.005, -0.005), (0.004, -0.002, -0.002))),
    ("ascani-component-mixing-cross-association", "ascani-2022-case-study-2-electrolyte-lle", 3,
     ("water", "1-butanol", "sodium-cation", "chloride-anion"), 298.15, 1e5, (0.80, 0.16, 0.02, 0.02),
     ((-0.02, 0.0, 0.01, 0.01), (0.01, -0.01, 0.0, 0.0), (0.004, -0.008, 0.002, 0.002))),
    ("mea-solvent-only-temperature-eps", "mea-co2-h2o-nine-species-estimation", 4,
     ("water", "monoethanolamine", "protonated-monoethanolamine", "carbamate-anion"),
     313.15, 1e5, (0.86, 0.08, 0.03, 0.03),
     ((-0.02, 0.0, 0.01, 0.01), (0.01, -0.01, 0.0, 0.0), (0.004, -0.008, 0.002, 0.002))),
]
TEMPERATURE_DIRECTION = (1.0, -0.5, 0.2)
PRESSURE_DIRECTION = (3e5, -2e5, 1e5)
OBSERVABLES = (("residual_enthalpy", PropertyObservable.ResidualEnthalpy),
               ("residual_entropy", PropertyObservable.ResidualEntropy),
               ("residual_isobaric_heat_capacity", PropertyObservable.ResidualIsobaricHeatCapacity))


EPS = 2.220446049250313e-16
R = 8.314462618


def stencil(model, T, P, x, directions, attribute, step, slope):
    """Mixed central difference at fixed pressure with its propagated closure/rounding noise.

    Each node's liquid root takes one Newton step at fixed (T, x) to its target pressure, so its
    remaining miss dP is rounding-level. The miss carries d(obs) = slope * dP / P_rho, where
    slope = (d obs / d rho)_{T,x} at the centre; each node also carries 8 eps |obs| rounding.
    """
    total, miss_sum, value_sum, nodes = 0.0, 0.0, 0.0, []
    for signs in itertools.product((-1.0, 1.0), repeat=len(directions)):
        t = T + step * sum(s * d[0] for s, d in zip(signs, directions))
        p = P + step * sum(s * d[1] for s, d in zip(signs, directions))
        y = [v + step * sum(s * d[2][i] for s, d in zip(signs, directions)) for i, v in enumerate(x)]
        root = model.state(t, P=p, x=y, phase="liquid")
        state = model.state(t, rho=root.molar_density + (p - root.pressure) / root.pressure_density_derivative,
                            x=y)
        value = getattr(state, attribute)
        miss_sum += (abs(state.pressure - p) + 4 * EPS * state.molar_density * R * t) / state.pressure_density_derivative
        value_sum += abs(value)
        nodes.append((int(root.density_branch), root.stable_root_count,
                      abs(state.pressure - p) / (state.molar_density * R * t)))
        total += math.prod(signs) * value
    scale = (2.0 * step) ** len(directions)
    # Floor inputs: sum |dP_node|/P_rho (density misses) and sum |obs_node| (rounding).
    return total / scale, (abs(slope) * miss_sum + 8 * EPS * value_sum) / scale, nodes, miss_sum, value_sum


def term_scale(model, T, P, x, directions):
    """Sum of |terms| of the multilinear action split into thermal, mechanical and composition parts.

    D^k obs[d_1..d_k] is the sum of these 3^k exact sub-direction actions; a relative criterion
    against a sum that cancels demands accuracy its terms cannot carry, so criteria scale with this sum.
    """
    zero = [0.0] * len(x)
    parts = [((d[0], 0.0, zero), (0.0, d[1], zero), (0.0, 0.0, list(d[2]))) for d in directions]
    total = [0.0] * len(OBSERVABLES)
    for combo in itertools.product(*parts):
        batch = model.pressure_actions(T, P, x=x, phase="liquid", observables=[o for _, o in OBSERVABLES],
                                       directions=[PressureActionDirection(*d) for d in combo])
        total = [t + abs(r.action) for t, r in zip(total, batch.results)]
    return total


def selection_rows():
    closure_limit = next(t.value for t in tolerances() if t.name == "DENSITY_ROOT_RESIDUAL")
    closure, actions, ladders = [], [], []
    for label, packet_id, version, components, T, P, x, compositions in SELECTIONS:
        model = Mixture(Parameters.from_json(packet(packet_id, version), components=components))
        for phase in ("vapor", "liquid"):
            s = model.state(T, P=P, x=x, phase=phase)
            closure.append({"selection": label, "T_K": T, "P_Pa": P, "x": " ".join(map(str, x)),
                            "requested_branch": phase, "returned_branch": int(s.density_branch),
                            "stable_roots": s.stable_root_count, "molar_density": s.molar_density,
                            "reduced_pressure_residual": abs(s.pressure - P) / (s.molar_density * R * T),
                            "pressure_density_derivative": s.pressure_density_derivative,
                            "bulk_relative_permittivity": s.bulk_relative_permittivity})
        centre = model.state(T, P=P, x=x, phase="liquid")
        slopes = model.density_actions(T, centre.molar_density, x=x, observables=[o for _, o in OBSERVABLES],
                                       directions=[DensityActionDirection(0.0, 1.0, [0.0] * len(x))]).results
        directions = [(TEMPERATURE_DIRECTION[k], PRESSURE_DIRECTION[k], compositions[k]) for k in range(3)]
        for order in (1, 2, 3):
            request = [PressureActionDirection(*d) for d in directions[:order]]
            batch = model.pressure_actions(T, P, x=x, phase="liquid",
                                           observables=[o for _, o in OBSERVABLES], directions=request)
            scales = term_scale(model, T, P, x, directions[:order])
            for (name, _), result, slope, scale in zip(OBSERVABLES, batch.results, slopes, scales):
                # Ladder h..8h. A central O(h^2) estimate D_k has truncation |D_k - D_k+1| / 3; the
                # Richardson value R_k = (4 D_k - D_k+1) / 3 is O(h^4) with truncation |R_k - R_k+1|.
                # Rounding noise falls as 1/h^order, so the smallest total picks the reference.
                steps = [m * FINE[order - 1] for m in (1, 2, 4, 8)]
                ladder = [stencil(model, T, P, x, directions[:order], name, h, slope.action) for h in steps]
                D, noise = [item[0] for item in ladder], [item[1] for item in ladder]
                rich = [(4 * D[k] - D[k + 1]) / 3 for k in range(3)]
                candidates = [(noise[k] + abs(D[k] - D[k + 1]) / 3, D[k], "central", k) for k in range(3)]
                candidates += [((4 * noise[k] + noise[k + 1]) / 3 + abs(rich[k] - rich[k + 1]), rich[k], "richardson", k)
                               for k in range(2)]
                uncertainty, reference, method, pick = min(candidates)
                for k, (estimate, item_noise, _, miss_sum, value_sum) in enumerate(ladder):
                    ladders.append({"selection": label, "order": order, "observable": name, "step": steps[k],
                                    "estimate": estimate, "richardson": rich[k] if k < 3 else "",
                                    "density_slope": slope.action, "sum_density_miss": miss_sum,
                                    "sum_abs_value": value_sum, "noise": item_noise,
                                    "truncation": abs(estimate - D[k + 1]) / 3.0 if k < 3 else "",
                                    "selected": method if k == pick else ""})
                nodes = [n for item in ladder for n in item[2]]
                admitted = all(n[:2] == (int(centre.density_branch), centre.stable_root_count)
                               and n[2] <= closure_limit for n in nodes)
                available = result.status == PropertyActionStatus.Available and result.action is not None
                tolerance = TOLERANCE[order - 1]
                criterion = tolerance * scale
                error = abs(result.action - reference) if available else math.nan
                verdict = ("reference_unresolved" if not admitted or uncertainty > criterion
                           else "unavailable" if not available else "pass" if error <= criterion else "fail")
                actions.append({"selection": label, "chart": "T,P,x", "order": order, "observable": name,
                                "status": result.status.name, "action": result.action,
                                "reference": reference, "reference_step": steps[pick], "reference_method": method,
                                "reference_uncertainty": uncertainty, "term_scale": scale,
                                "relative_error": error / abs(reference), "scaled_error": error / scale,
                                "tolerance": tolerance, "uncertainty_over_criterion": uncertainty / criterion,
                                "error_over_floor": error / uncertainty,
                                "max_node_reduced_residual": max(n[2] for n in nodes), "verdict": verdict})
    return closure, actions, ladders


def solvent_rows():
    references = list(csv.DictReader((LEAF / "inputs/nist-srd69-solvent-isobars.tsv").open(), delimiter="\t"))
    models = [("water", "figiel-yu-held-2025", 3), ("water", "held-cameretti-sadowski-2008", 3),
              ("water", "mea-co2-h2o-nine-species-estimation", 4), ("methanol", "figiel-yu-held-2025", 3)]
    rows = []
    for solvent, packet_id, version in models:
        path = packet(packet_id, version)
        pure = Pure(Parameters.from_json(path, components=(solvent,)))
        mass = molar_mass(path, solvent)
        for ref in (r for r in references if r["solvent"] == solvent):
            T, P, rho_ref = float(ref["T_K"]), 1e6 * float(ref["P_MPa"]), float(ref["mass_density_kg_m3"])
            s = pure.state(T, P=P, phase="liquid")
            rows.append({"solvent": solvent, "packet": f"{packet_id}/{version}", "T_K": T, "P_Pa": P,
                         "model_kg_m3": s.molar_density * mass, "reference_kg_m3": rho_ref,
                         "relative_deviation": s.molar_density * mass / rho_ref - 1.0,
                         "stable_roots": s.stable_root_count, "reference": ref["source"]})
    return rows


def write(name, rows):
    with (LEAF / "results" / name).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    sys.path.insert(0, str(ROOT / "tools"))
    from engine_identity import loaded_identity
    prepared = loaded_identity()
    if not prepared["source_revision"] or not prepared["meson_configuration_sha256"] or not prepared["lock_sha256"]:
        raise RuntimeError("the selected Engine must identify its prepared source, Meson build, and lock")
    closure, actions, ladders = selection_rows()
    solvents = solvent_rows()
    if loaded_identity() != prepared:
        raise RuntimeError("Engine identity changed during derivative qualification")
    write("pressure-closure.csv", closure)
    write("pressure-actions.csv", actions)
    write("reference-ladder.csv", ladders)
    write("solvent-density.csv", solvents)
    summary = {
        "question": "Do the admitted real selections keep complete pressure-closed derivative chains, "
                    "and is neutral-solvent density adequate before salt-density errors are attributed to ions?",
        "engine_version": prepared["version"],
        "prepared_environment": prepared,
        "engine_tree": subprocess.run(["git", "rev-parse", f"{prepared['source_revision']}:engine"],
                                      cwd=prepared["source_root"], capture_output=True,
                                      text=True).stdout.strip(),
        "packets": {f"{p}/{v}": hashlib.sha256(packet(p, v).read_bytes()).hexdigest()
                    for _, p, v, *_ in SELECTIONS},
        "actions": {v: sum(r["verdict"] == v for r in actions) for v in
                    ("pass", "fail", "unavailable", "reference_unresolved")},
        "worst_scaled_action_error_by_order": {o: max(r["scaled_error"] for r in actions if r["order"] == o)
                                               for o in (1, 2, 3)},
        "worst_uncertainty_over_criterion": max(r["uncertainty_over_criterion"] for r in actions),
        "max_reduced_pressure_residual": max(r["reduced_pressure_residual"] for r in closure),
        "solvent_max_abs_relative_deviation": {
            f"{r['solvent']}@{r['packet']}": max(abs(q["relative_deviation"]) for q in solvents
                                                 if q["solvent"] == r["solvent"] and q["packet"] == r["packet"])
            for r in solvents},
        "claim_limits": [
            "Finite states and directions; exact differentiation does not establish physical accuracy.",
            "Figiel/Held/Ascani packets report 298.15 K and 1 bar; other temperatures are extrapolation.",
            "Parameter directions and the DensityTemperature/DipolarPerturbation forms are qualified by "
            "native.properties.issue26, not by these public-API rows.",
            "Criteria are tolerance x term_scale, the sum of |thermal, mechanical, composition sub-direction "
            "terms| of each action; a pass also needs the reference uncertainty (Newton-polished node closure, "
            "rounding and |D(h)-D(2h)|/3 truncation) below that criterion.",
        ],
    }
    (LEAF / "results/summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    json.dump(summary["actions"], sys.stdout)
    print()
    if summary["actions"]["pass"] != len(actions):
        raise SystemExit(1)


def declarations():
    return [
        Mixture(Parameters.from_json(packet(name, version), components=components))
        for _, name, version, components, *_ in SELECTIONS
    ]


if __name__ == "__main__":
    main()
