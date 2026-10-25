"""Evaluate the retained HELD2 Khudaida states under the current Engine and pin the lineage difference.

For each Figure 2 feed, the historical two-liquid state is evaluated under the packet-v6 model
(total G/RT, cross-phase chemical-potential differences, material balance) and compared with the
current declared two-liquid solution. The same solve is repeated with the unlike-pair association
volumes evaluated by the cited Wolbach-Sandler/Gross-Sadowski rule including its size factor,
which is how the historical bundle 2 stored them.
"""
from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import importlib.metadata
import json
import math
import tempfile
from pathlib import Path

import epcsaft
from epcsaft import equilibrium

C = ("water", "ethanol", "isobutanol", "sodium-cation", "chloride-anion")
T, P = 293.15, 1.0e5
ROOT = Path(__file__).resolve().parents[3]
V6 = ROOT / "data/packets/khudaida-2026-figure-2-electrolyte-lle/6/parameters/parameters.json"
REPRO = ROOT / "analyses/reproduction/2026-khudaida"
HIST = REPRO / "results/data/khudaida-2026-figure-2-held2.csv"
BUNDLE2 = REPRO / "inputs/parameters/bundle/khudaida-2026-figure-2-electrolyte-lle/2/association.csv"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def model_from(doc: dict) -> epcsaft.Mixture:
    with tempfile.NamedTemporaryFile("w", suffix=".json") as stream:
        json.dump(doc, stream); stream.flush()
        return epcsaft.Mixture(epcsaft.Parameters.from_json(Path(stream.name), components=C))


def sigma(doc: dict, component: str) -> float:
    for c in doc["correlations"]:
        if c["component_id"] == component and c["family"] == "segment_diameter":
            return c["constant"]["value"]["magnitude"] + sum(
                t["amplitude"]["value"]["magnitude"] * math.exp(t["exponent_coefficient"]["value"]["magnitude"] * T)
                for t in c.get("terms", []))
    record = next(c for c in doc["components"] if c["component_id"] == component)
    return next(k["value"]["magnitude"] for k in record["coefficients"] if k["family"] == "segment_diameter")


def with_size_factor(doc: dict) -> tuple[dict, dict[str, dict[str, float]]]:
    """Replace each combining-rule edge by explicit values; volume carries (sqrt(si sj)/((si+sj)/2))^3."""
    doc = copy.deepcopy(doc)
    edges = doc["topology"]["edges"]
    self_values = {e[k]["identity"]: e[k]["value"]["magnitude"]
                   for e in edges if e["source"]["kind"] == "explicit" for k in ("energy_over_k", "volume")}
    volumes = {}
    for e in edges:
        source = e["source"]
        if source["kind"] != "combining-rule":
            continue
        i, j = e["endpoint_a"]["component_id"], e["endpoint_b"]["component_id"]
        si, sj = sigma(doc, i), sigma(doc, j)
        geometric = math.sqrt(math.prod(self_values[k] for k in source["inputs"]["volume"]))
        volume = geometric * (math.sqrt(si * sj) / ((si + sj) / 2)) ** 3
        volumes["/".join(sorted((i, j)))] = {"geometric_mean": geometric, "with_size_factor": volume}
        tag = f"size-factor/{i}/{e['endpoint_a']['site_id']}/{j}/{e['endpoint_b']['site_id']}"
        e["energy_over_k"] = {"identity": f"{tag}/energy_over_k", "value": {
            "magnitude": sum(self_values[k] for k in source["inputs"]["energy_over_k"]) / 2, "unit": "kelvin"}}
        e["volume"] = {"identity": f"{tag}/volume", "value": {"magnitude": volume, "unit": "dimensionless"}}
        e["source"] = {"kind": "explicit", "provenance": source["provenance"]}
    return doc, volumes


def evaluate(model: epcsaft.Mixture, phases: list[tuple[float, list[float]]]) -> tuple[float, list[list[float]]]:
    """G/RT up to a feed-fixed constant, and ln(x_i phi_i) per phase, at T and P on the liquid root."""
    g, mu = 0.0, []
    for amount, x in phases:
        state = model.state(T, P=P, x=x, phase="liquid")
        phase_mu = [math.log(xi) + lnphi for xi, lnphi in zip(x, state.log_fugacity_coefficient, strict=True)]
        g += amount * math.fsum(xi * m for xi, m in zip(x, phase_mu, strict=True))
        mu.append(phase_mu)
    return g, mu


def dmu(mu: list[list[float]]) -> tuple[float, float]:
    """Largest neutral-species and NaCl-combination chemical-potential difference between the phases, over RT."""
    aqueous, organic = mu
    return (max(abs(aqueous[i] - organic[i]) for i in range(3)),
            abs(aqueous[3] + aqueous[4] - organic[3] - organic[4]))


def solve(model: epcsaft.Mixture, feed: list[float], previous):
    problem = equilibrium.Problem(phases=[
        equilibrium.Phase("phase-0", kind="liquid", support=C, composition_guess=(.8, .05, .05, .05, .05)),
        equilibrium.Phase("phase-1", kind="liquid", support=C, composition_guess=(.2, .1, .6998, .0001, .0001))],
        T=T, P=P, feed=equilibrium.Amounts(dict(zip(C, feed, strict=True))), row_tolerance=1e-10)
    result = equilibrium.solve_equilibrium(model, problem, previous)
    if not result.success or result.requested_tolerance_met is not True or result.phase_count != 2:
        raise RuntimeError(result.message)
    phases = sorted(((result.amounts[k], list(result.mole_fractions[5 * k:5 * k + 5])) for k in range(2)),
                    key=lambda p: -p[1][0])
    return result, phases


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    output = parser.parse_args().output_dir
    output.mkdir(parents=True, exist_ok=False)
    v6 = json.loads(V6.read_text())
    size_doc, volumes = with_size_factor(v6)
    bundle2 = {"/".join(sorted((r["component_id_a"], r["component_id_b"]))): float(r["value"])
               for r in csv.DictReader(BUNDLE2.open()) if r["family"] == "association_volume"
               and r["component_id_a"] != r["component_id_b"]}
    for pair, record in volumes.items():
        record["historical_bundle_2"] = bundle2[pair]
    models = {"v6": model_from(v6), "size_factor": model_from(size_doc)}
    historical = {int(r["point_id"]): json.loads(r["result_json"]) for r in csv.DictReader(HIST.open())}
    rows, previous = [], dict.fromkeys(models)
    for point in range(1, 9):
        h = historical[point]
        feed = h["overall_mole_fractions"]
        h_phases = sorted(((p["amount_mol"], p["mole_fractions"]) for p in h["phases"]), key=lambda p: -p[1][0])
        row = {"point_id": point, "historical_g_rt": h["total_free_energy_over_rt"],
               "historical_material_max_abs": max(abs(math.fsum(n * x[i] for n, x in h_phases) - feed[i]) for i in range(5))}
        for name, model in models.items():
            g_hist, mu_hist = evaluate(model, h_phases)
            result, phases = solve(model, feed, previous[name]); previous[name] = result
            g_solved, mu_solved = evaluate(model, phases)
            row |= {f"{name}_g_rt": result.objective,
                    f"{name}_g_rt_excess_of_historical_state": g_hist - g_solved,
                    f"{name}_dmu_neutral_historical_state": dmu(mu_hist)[0],
                    f"{name}_dmu_nacl_historical_state": dmu(mu_hist)[1],
                    f"{name}_dmu_neutral_solved": dmu(mu_solved)[0],
                    f"{name}_dmu_nacl_solved": dmu(mu_solved)[1],
                    f"{name}_x_max_abs_vs_historical": max(abs(a - b) for (_, xa), (_, xb) in zip(phases, h_phases) for a, b in zip(xa, xb)),
                    f"{name}_amount_max_abs_vs_historical": max(abs(a - b) for (a, _), (b, _) in zip(phases, h_phases)),
                    f"{name}_g_rt_abs_vs_historical": abs(result.objective - h["total_free_energy_over_rt"])}
        rows.append(row)
    with (output / "diagnosis.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)
    wheel = Path(json.loads(importlib.metadata.distribution("epcsaft").read_text("direct_url.json"))["url"].removeprefix("file://"))
    maxima = {key: max(abs(r[key]) for r in rows) for key in rows[0] if key not in ("point_id", "historical_g_rt")
              and not key.endswith("_g_rt")}
    maxima["v6_g_rt_excess_of_historical_state_min"] = min(r["v6_g_rt_excess_of_historical_state"] for r in rows)
    summary = {"points": len(rows), "temperature_K": T, "pressure_Pa": P,
               "association_volumes": volumes, "maxima": maxima,
               "identity": {"wheel_sha256": sha(wheel), "epcsaft_version": importlib.metadata.version("epcsaft"),
                            "v6_parameters_sha256": sha(V6), "historical_rows_sha256": sha(HIST),
                            "bundle_2_association_sha256": sha(BUNDLE2), "diagnosis_csv_sha256": sha(output / "diagnosis.csv")},
               "claim_limits": ["declared two-liquid local branches at the eight Figure 2 feeds only",
                                "G/RT excess uses the sum of n x ln(x phi); the feed-fixed reference constant cancels",
                                "no global-stability or predictive-accuracy claim"]}
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    return 0


def declarations():
    return model_from(json.loads(V6.read_text()))


if __name__ == "__main__":
    raise SystemExit(main())
