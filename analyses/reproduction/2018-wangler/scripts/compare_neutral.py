"""Retain Wangler neutral comparisons and check their numerical claim limits.

Run neutral_engine.py and neutral_feos.py first. Raw state diagnostics stay under
build/issue-857/build-r1; homogeneous comparison values are owned by results/*.csv.
"""

import csv
import hashlib
import json
import math
from pathlib import Path
import subprocess

LEAF = Path(__file__).resolve().parents[1]
ROOT = LEAF.parents[2]
RAW = ROOT / "build/issue-857/build-r1"
RESULTS = LEAF / "results"


def read(name):
    return json.loads((RAW / name).read_text())


def write(name, rows):
    with (RESULTS / name).open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0], lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def kim(t):
    # Kim p. 2523 Table 5: Riedel B/T reading and Antoine with mmHg -> Pa.
    return (
        math.exp(82.8943 - 11011.8 / t - 8.03383 * math.log(t) + 1.68146e-18 * t**6),
        133.32236842105263 * 10 ** (10.182 - 3578.6 / (t - 273.15 + 264.24)),
    )


def main():
    engine, feos = read("engine-cross-r2.json"), read("feos-cross-r2.json")
    el, fl = engine["phases"][0], feos["phases"][0]
    ev, fv = engine["phases"][1], feos["phases"][1]
    fixed = feos["at_engine_liquid"]
    metrics = [
        ("Bubble pressure / bar", engine["P"] / 1e5, fv["P"] / 1e5),
        ("Liquid molar density / mol m^-3", el["rho"], fl["rho"]),
        (
            "Liquid mass density / kg m^-3",
            el["rho"] * (0.8 * 0.01801528 + 0.2 * 0.119163),
            fl["rho"] * (0.8 * 0.01801528 + 0.2 * 0.119163),
        ),
        ("Vapor molar density / mol m^-3", ev["rho"], fv["rho"]),
        ("Vapor water mole fraction", ev["x"][0], fv["x"][0]),
        ("Water fugacity / Pa", el["fugacity_Pa"][0], fl["fugacity_Pa"][0]),
        ("Liquid mu_water^res/(RT)", el["mu_res_over_RT"][0], fl["mu_res_over_RT"][0]),
        ("Liquid mu_MDEA^res/(RT)", el["mu_res_over_RT"][1], fl["mu_res_over_RT"][1]),
    ]
    write(
        "cross-code.csv",
        [
            dict(quantity=q, Engine=a, FeOs=b, FeOs_minus_Engine=b - a)
            for q, a, b in metrics
        ],
    )
    assert all(abs(b / a - 1) < 1e-6 for _, a, b in metrics[:6])
    derivatives = engine["term_derivative_verification"]
    mu = derivatives[2]["mu_terms_over_RT"]
    termrows = []
    for term in mu:
        row = dict(
            term=term,
            Engine_a_over_RT=el["terms"][term],
            FeOs_a_over_RT=fixed["terms"][term],
        )
        for i, species in enumerate(("water", "mdea")):
            row[f"Engine_mu_{species}_over_RT"] = mu[term][i]
            row[f"FeOs_mu_{species}_over_RT"] = fixed["mu_terms_over_RT"][term][i]
        termrows.append(row)
        assert abs(el["terms"][term] - fixed["terms"][term]) < 1e-7
        assert all(
            abs(a - b) < 1e-7
            for a, b in zip(mu[term], fixed["mu_terms_over_RT"][term], strict=True)
        )
    write("cross-code-terms.csv", termrows)
    write(
        "term-derivative-steps.csv",
        [
            dict(h_mol=r["h_mol"], term=k, mu_water_over_RT=v[0], mu_mdea_over_RT=v[1])
            for r in derivatives
            for k, v in r["mu_terms_over_RT"].items()
        ],
    )
    term_sum_error = max(
        abs(sum(mu[k][j] for k in mu) - el["mu_res_over_RT"][j]) for j in range(2)
    )
    step_change = max(
        abs(r["mu_terms_over_RT"][k][j] - mu[k][j])
        for r in derivatives[1:]
        for k in mu
        for j in range(2)
    )
    assert term_sum_error < 1e-9 and step_change < 1e-9
    assert (
        max(
            abs(math.log(fl["fugacity_Pa"][j] / fv["fugacity_Pa"][j])) for j in range(2)
        )
        < 1e-9
    )
    sweep = read("convention-sweep-r2.json")
    source = {
        float(r["T_K"]): r
        for r in csv.DictReader(
            (LEAF / "inputs/source/wangler-2018-figure-1.csv").open()
        )
    }
    neutral = []
    for row in sweep:
        if row["mode"] == "baseline":
            s = source[row["T"]]
            neutral.append(
                dict(
                    T_K=row["T"],
                    x_mdea=0.2,
                    source_curve_bar=float(s["source_curve_bar"]),
                    graphical_allowance_bar=float(s["graphical_allowance_bar"]),
                    Engine_P_bar=row["P"] / 1e5,
                    Engine_minus_source_bar=row["P"] / 1e5
                    - float(s["source_curve_bar"]),
                    maximum_residual=row["residual"],
                )
            )
    write("neutral-figure-1.csv", neutral)
    choices = []
    printed = {
        "baseline",
        "eq12-literal",
        "literal-k",
        "literal-l",
        "mdea-origin-digits",
        "association-diameter-d",
    }
    for mode in dict.fromkeys(r["mode"] for r in sweep):
        rows = [r for r in sweep if r["mode"] == mode]
        assert len(rows) == 4
        solved = all(r.get("success") for r in rows)
        choice = dict(
            mode=mode,
            basis="printed reading or parameter source"
            if mode in printed
            else "excluded diagnostic control",
        )
        choice.update(
            {
                f"P_{int(r['T'])}_K_bar": r["P"] / 1e5 if r.get("success") else None
                for r in rows
            }
        )
        choice["maximum_error_at_x_0_2_bar"] = (
            max(
                abs(r["P"] / 1e5 - float(source[r["T"]]["source_curve_bar"]))
                for r in rows
            )
            if solved
            else None
        )
        choice["all_four_within_0_01_bar"] = (
            solved and choice["maximum_error_at_x_0_2_bar"] <= 0.01
        )
        choice["failure"] = "; ".join(
            r.get("error", r.get("message", "")) for r in rows if not r.get("success")
        )
        choices.append(choice)
    assert not any(r["all_four_within_0_01_bar"] for r in choices)
    write("convention-sweep.csv", choices)
    literal = next(r for r in sweep if r["mode"] == "eq12-literal" and r["T"] == 353.0)
    literal_feos = read("feos-literal-eq12-r2.json")["phases"][1]["P"]
    assert abs(literal_feos / literal["P"] - 1) < 1e-6
    write(
        "cross-code-literal-eq12.csv",
        [
            dict(
                T_K=353.0,
                x_mdea=0.2,
                Engine_P_bar=literal["P"] / 1e5,
                FeOs_P_bar=literal_feos / 1e5,
            )
        ],
    )
    pure = []
    for r, f in zip(read("pure-mdea-r2.json"), read("feos-pure-r2.json"), strict=True):
        assert r["T"] == f["T"] and abs(r["P"] / f["P"] - 1) < 1e-6
        assert abs(r["phases"][0]["rho"] / f["rho"] - 1) < 1e-6
        kr, ka = kim(r["T"])
        pure.append(
            dict(
                T_K=r["T"],
                Engine_P_Pa=r["P"],
                FeOs_P_Pa=f["P"],
                Engine_rho_kg_m3=r["phases"][0]["rho"] * 0.119163,
                FeOs_rho_kg_m3=f["rho_mass"],
                Kim_Riedel_P_Pa=kr,
                Kim_Antoine_P_Pa=ka,
                pahlav_Table1_literal_P_Pa=r["pahlav_Table1_as_printed_P_Pa"],
                pahlav_Table1_literal_rho_kg_m3=r["pahlav_Table1_as_printed_rho_mol_m3"]
                * 0.119163,
            )
        )
    write("pure-mdea.csv", pure)
    obs = []
    for r in read("pure-observations-r2.json"):
        kr, ka = kim(r["T"])
        obs.append(
            dict(
                T_K=r["T"],
                observed_P_Pa=r["observed_P_Pa"],
                Engine_P_Pa=r["P"],
                Kim_Riedel_P_Pa=kr,
                Kim_Antoine_P_Pa=ka,
                Engine_rho_kg_m3=r["phases"][0]["rho"] * 0.119163,
                relative_error_percent=100 * (r["P"] / r["observed_P_Pa"] - 1),
            )
        )
    write("pure-observations.csv", obs)
    diagnostics = dict(
        bubble_relative_difference=abs(fv["P"] / engine["P"] - 1),
        maximum_term_mu_difference=max(
            abs(mu[k][j] - fixed["mu_terms_over_RT"][k][j])
            for k in mu
            for j in range(2)
        ),
        term_sum_error=term_sum_error,
        term_step_change=step_change,
        common_density_pressure_difference_Pa=fixed["P"] - el["P"],
        Engine_bubble_residual=engine["residual"],
        FeOs_pressure_closure_Pa=abs(fl["P"] - fv["P"]),
        FeOs_association_residual="not exposed by FeOs Python interface",
        Kim_neat_measurement_AARD_percent=sum(
            abs(r["relative_error_percent"]) for r in obs
        )
        / len(obs),
        claim_limits="Numerical verification at declared neutral states; no temperature derivatives, ionic model or global stability qualification.",
    )
    (RESULTS / "neutral-diagnostics.json").write_text(
        json.dumps(diagnostics, indent=2) + "\n"
    )
    pressure = json.loads(
        (
            LEAF / "figures/figure-14/output/wangler-2018-table-10-native-summary.json"
        ).read_text()
    )
    assert pressure["declaration"]["reference_pressure_pa"] is None
    for lane, expected in (("S0", 96.2488724655), ("S1", 89.5500626566)):
        result = pressure["lanes"][f"fitted-{lane}"]
        assert result["solved_row_count"] == 3
        assert result["maximum_original_row"] < 1e-10
        assert result["F2"]["references_admitted"]
        # Six-decimal AARD reproduction from the independently retained diagnosis.
        assert abs(result["native_ard_percent"] - expected) < 1e-6
    identity_path = RAW / "engine-identity.json"
    paths = list((LEAF / "scripts").glob("*.py")) + list(RAW.glob("*r2.json"))
    provenance = dict(
        base_commit=subprocess.check_output(
            ["git", "rev-parse", "origin/main"], cwd=ROOT, text=True
        ).strip(),
        files={
            str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in paths
        },
        Engine_identity=json.loads(identity_path.read_text()),
        FeOs_identity=read("feos-identity.json"),
        comparison_limits="1e-6 relative equilibrium observables; 1e-7 absolute a/RT and mu/RT",
        derivative_verification="Engine term mu from five-point differences of native energies; exact Engine totals and exact FeOs term derivatives",
    )
    (RESULTS / "neutral-provenance.json").write_text(
        json.dumps(provenance, indent=2) + "\n"
    )
    print(json.dumps(diagnostics, indent=2))


if __name__ == "__main__":
    main()
