"""Decompose retained pressures about pure-endpoint straight lines; no model execution."""

import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/pressure-diagnosis"


def read(name):
    return json.loads((OUT / name).read_text())


def write(name, value):
    (OUT / name).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def table(name, headers, rows):
    (OUT / name).write_text(
        "\n".join(
            [
                "| " + " | ".join(headers) + " |",
                "| " + " | ".join(["---"] * len(headers)) + " |",
                *["| " + " | ".join(map(str, row)) + " |" for row in rows],
            ]
        )
        + "\n"
    )


def main():
    probes, measured = read("probe-results.json"), read("measurement-check.json")

    def probe(label, x=None):
        return next(r for r in probes if r["probe"] == label and r["x_water"] == x)

    endpoints = {r["reported_y_water"]: r for r in measured["pure_endpoints"]}
    pw, pa = (
        probe("pure-baseline-" + name)["P_Pa"] for name in ("water", "acetic-acid")
    )
    ew, ea = endpoints[1]["P_Pa"], endpoints[0]["P_Pa"]
    reactive = {r["x_water"]: r for r in probes if r["probe"] == "reactive-baseline"}
    decomposition, closure = [], []
    for m in measured["measurements"]:
        x, pe = m["x_water"], m["P_Pa"]
        pn = probe("neutral-baseline", x)["P_Pa"]
        cm, ce = pn - (x * pw + (1 - x) * pa), pe - (x * ew + (1 - x) * ea)
        r = dict(
            row=m["row"],
            x_water=x,
            measured_P_Pa=pe,
            neutral_P_Pa=pn,
            water_endpoint_error_Pa=x * (pw - ew),
            acid_endpoint_error_Pa=(1 - x) * (pa - ea),
            model_chord_residual_Pa=cm,
            measured_chord_residual_Pa=ce,
            mixing_error_Pa=cm - ce,
        )
        pr = pn
        if x in reactive:
            pr = reactive[x]["P_Pa"]
            r.update(primary_P_Pa=pr, reacting_minus_neutral_Pa=pr - pn)
        closure.append(
            abs(
                sum(
                    r[k]
                    for k in (
                        "water_endpoint_error_Pa",
                        "acid_endpoint_error_Pa",
                        "mixing_error_Pa",
                    )
                )
                + (pr - pn)
                - (pr - pe)
            )
        )
        decomposition.append(r)
    write("pressure-decomposition.json", decomposition)
    fields = list(dict.fromkeys(k for r in decomposition for k in r))
    with (OUT / "pressure-decomposition.csv").open("w") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(decomposition)
    central = [r for r in decomposition if "primary_P_Pa" in r]
    table(
        "decomposition-table.md",
        [
            "Analytical liquid water fraction",
            "Water / Pa",
            "Acid / Pa",
            "Mixing / Pa",
            "Dissociation / Pa",
            "Total error / Pa",
        ],
        [
            [
                f"{r['x_water']:.4f}",
                *[
                    f"{r[k]:+.6f}"
                    for k in (
                        "water_endpoint_error_Pa",
                        "acid_endpoint_error_Pa",
                        "mixing_error_Pa",
                        "reacting_minus_neutral_Pa",
                    )
                ],
                f"{r['primary_P_Pa'] - r['measured_P_Pa']:+.6f}",
            ]
            for r in central
        ],
    )
    table(
        "departure-table.md",
        [
            "Water fraction",
            "Model departure / Pa",
            "Measured departure / Pa",
            "Error with measured endpoints and unchanged model departure / Pa",
        ],
        [
            [
                f"{r['x_water']:.4f}",
                f"{r['model_chord_residual_Pa']:.6f}",
                f"{r['measured_chord_residual_Pa']:.6f}",
                f"{r['mixing_error_Pa'] + r['reacting_minus_neutral_Pa']:+.6f}",
            ]
            for r in central
        ],
    )
    table(
        "endpoints-table.md",
        [
            "Component",
            "Native saturation pressure / kPa",
            "ThermoML endpoint ± expanded uncertainty / kPa",
            "Calculated minus measured / kPa",
            "Relative error / %",
        ],
        [
            [
                name,
                f"{pm / 1000:.9f}",
                f"{pe / 1000:.3f} ± {endpoints[y]['U95_P_Pa'] / 1000:.3f}",
                f"{(pm - pe) / 1000:+.9f}",
                f"{100 * (pm / pe - 1):+.4f}",
            ]
            for name, pm, pe, y in (("Water", pw, ew, 1), ("Acetic acid", pa, ea, 0))
        ],
    )
    sources = json.loads((ROOT / "inputs/source/source-checks.json").read_text())
    antoine, fit = sources["NIST_Antoine_acid"], sources["pabsch_SI"]
    (OUT / "source-comparison.md").write_text(
        f"The [NIST WebBook Antoine table]({antoine['url']}), attributed to {antoine['source']}, gives "
        f"**{antoine['P_Pa_at_323_15_K'] / 1000:.9f} kPa** at **{measured['temperature_K']:.2f} K**, "
        f"compared with the native acid pressure of **{pa / 1000:.9f} kPa**. "
        f"The difference is **{(pa - antoine['P_Pa_at_323_15_K']) / 1000:+.9f} kPa**. "
        "The correlation has no supplied uncertainty and is not an exact reference.\n"
    )
    (OUT / "fit-conditions.md").write_text(
        f"Pabsch’s supporting information, Tables S1 and S4, identifies **$k_{{ij}}={fit['kij']}$** "
        f"as fitted in that work, with a printed binary vapour–liquid fitting interval of "
        f"**{fit['fit_temperature_range_printed_K'][0]:.2f}–{fit['fit_temperature_range_printed_K'][1]:.2f} K** "
        f"and average relative temperature deviation **{fit['reported_ARD_T_percent']:.2f}%**; no pressure deviation is reported. "
        f"Reference 9 is **{fit['fit_source']}**. "
        f"Figure S1 compares coexistence at **{fit['figure_S1_pressure_bar']:.5f} bar**, "
        "over an approximately **373–393 K** plotted range. "
        f"The present **{measured['temperature_K']:.2f} K** isotherm is below the printed fitting interval.\n"
    )
    sensitivities = [
        dict(
            probe=r["probe"],
            x_water=r["x_water"],
            pressure_Pa=r["P_Pa"],
            change_from_neutral_Pa=r["P_Pa"]
            - probe("neutral-baseline", r["x_water"])["P_Pa"],
        )
        for r in probes
        if r["x_water"] in reactive
        and r["probe"] not in ("neutral-baseline", "reactive-baseline")
        and "P_Pa" in r
    ]
    write("pressure-sensitivities.json", sensitivities)
    table(
        "sensitivities-table.md",
        ["Neutral bubble-pressure probe", *[f"ΔP at x={x:.4f} / Pa" for x in reactive]],
        [
            [
                label,
                *[
                    f"{next(r['change_from_neutral_Pa'] for r in sensitivities if r['probe'] == label and r['x_water'] == x):+.6f}"
                    for x in reactive
                ],
            ]
            for label in dict.fromkeys(r["probe"] for r in sensitivities)
        ],
    )
    table(
        "vapour-table.md",
        [
            "Liquid water fraction",
            "Vapour water fraction",
            "Vapour Z",
            "Water unbonded-site fraction",
            "Acid unbonded-site fraction",
            "Vapour association / RT",
        ],
        [
            [
                f"{x:.4f}",
                f"{probe('neutral-baseline', x)['mole_fractions'][2]:.6f}",
                *[
                    f"{v:.6f}"
                    for v in (
                        probe("neutral-baseline", x)["states"][1][
                            "compressibility_factor"
                        ],
                        probe("neutral-baseline", x)["states"][1]["site_fractions"][0],
                        probe("neutral-baseline", x)["states"][1]["site_fractions"][2],
                        probe("neutral-baseline", x)["states"][1]["association"],
                    )
                ],
            ]
            for x in reactive
        ],
    )
    table(
        "same-state-table.md",
        [
            "Removed vapour interactions at fixed density and composition",
            *[f"ΔP at x={x:.4f} / Pa" for x in reactive],
        ],
        [
            [
                label,
                *[
                    f"{probe('vapor-' + control, x)['state']['pressure'] - probe('neutral-baseline', x)['states'][1]['pressure']:+.6f}"
                    for x in reactive
                ],
            ]
            for label, control in (
                ("Acid self association", "acid-self-off"),
                ("Cross association", "cross-off"),
                ("All association", "all-association-off"),
            )
        ],
    )
    eq = [r for r in probes if "P_Pa" in r]
    verification = dict(
        equilibrium_solutions=len(eq),
        direct_state_evaluations=len(probes) - len(eq),
        decomposition_rows_verified=len(decomposition),
        decomposition_max_absolute_error_Pa=max(closure),
        max_row_residual=max(r["row_residual"] for r in eq),
        max_relative_pressure_closure=max(r["relative_pressure_closure"] for r in eq),
        max_log_fugacity_mismatch=max(
            abs(v) for r in eq for v in r["log_fugacity_mismatch"]
        ),
    )
    assert max(closure) < 1e-8, "Pure-endpoint pressure decomposition does not close"
    write("verification.json", verification)
    (OUT / "verification.md").write_text(
        f"The retained probes contain **{len(eq)} native equilibrium solutions and {len(probes) - len(eq)} direct vapour-state evaluations**. "
        f"The maximum row residual is **{verification['max_row_residual']:.2e}**, relative pressure closure **{verification['max_relative_pressure_closure']:.2e}**, "
        f"and log-fugacity mismatch **{verification['max_log_fugacity_mismatch']:.2e}**. "
        f"The four-term decomposition closes over all **{len(decomposition)} mixed compositions** within **{max(closure):.2e} Pa**.\n"
    )
    print(json.dumps(verification, indent=2))


if __name__ == "__main__":
    main()
