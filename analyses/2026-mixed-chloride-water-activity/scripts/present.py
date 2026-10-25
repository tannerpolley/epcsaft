"""Format retained comparison numbers and figures; never import the Engine."""

import csv
import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parents[1]
RESULTS, FIGURES = HERE / "results", HERE / "figures"
SOURCES = ("robinson-bower-1966", "rard-miller-1987", "dinane-mounir-2003")
LABELS = {
    SOURCES[0]: "Robinson–Bower: NaCl–CaCl₂",
    SOURCES[1]: "Rard–Miller: NaCl–MgCl₂",
    SOURCES[2]: "Dinane–Mounir: NaCl–MgCl₂",
}


def read(name):
    with (RESULTS / name).open() as stream:
        return list(csv.DictReader(stream))


def number(value):
    return f"{float(value):.6g}" if value not in (None, "") else "Unavailable"


def table(headers, rows):
    return (
        "| "
        + " | ".join(headers)
        + " |\n| "
        + " | ".join(["---"] * len(headers))
        + " |\n"
        + "".join("| " + " | ".join(str(v) for v in row) + " |\n" for row in rows)
        + "\n"
    )


def main():
    FIGURES.mkdir(exist_ok=True)
    summary = json.loads((RESULTS / "summary.json").read_text())
    spec = json.loads((HERE / "inputs/comparison-specification.json").read_text())
    verdicts, nominal, buckets = (
        read("verdicts.csv"),
        read("nominal-residuals.csv"),
        read("bucket-metrics.csv"),
    )
    transfers, binary = read("transfer-tests.csv"), read("binary-reference-checks.csv")
    source_counts = {
        source: sum(
            r["source"] == source
            and r["declaration"] == "Held"
            and r["role"] == "mixture"
            for r in nominal
        )
        for source in SOURCES
    }
    application_basis = json.loads((HERE / "inputs/application-basis.json").read_text())
    source_checks = json.loads((HERE / "inputs/source-data-checks.json").read_text())
    content = "## Conditions and fixed limits\n\n"
    content += table(
        ["Quantity", "Value", "Basis"],
        [
            ["Temperature", f"{summary['conditions']['T_K']} K", "All native states"],
            [
                "Pressure",
                f"{summary['conditions']['P_Pa']:g} Pa",
                "Declared liquid calculation; experimental pressure correction unavailable",
            ],
            [
                "Phi RMS / maximum error",
                " / ".join(number(v) for v in spec["primary_limits"]["Phi"].values()),
                "Isopiestic primary limits",
            ],
            [
                "Water activity RMS / maximum error",
                " / ".join(number(v) for v in spec["primary_limits"]["aw"].values()),
                "Separate hygrometric limits",
            ],
            [
                "Frozen realizations",
                summary["realizations"],
                "Finite conditional input/reference/stock cases",
            ],
            [
                "Mixture rows",
                f"{source_counts[SOURCES[0]]} Ca + {source_counts[SOURCES[1]]} Mg isopiestic; {source_counts[SOURCES[2]]} Mg hygrometric",
                "References excluded from mixture metrics",
            ],
        ],
    )
    (RESULTS / "conditions.md").write_text(content)
    content = table(
        [
            "System / source",
            "I bucket (mol/kg)",
            "Quantity",
            "Declaration",
            "n",
            "Bias",
            "MAE",
            "RMS",
            "Maximum",
            "Verdict",
        ],
        [
            [
                LABELS[r["source"]],
                r["I_bucket"],
                r["quantity"],
                r["declaration"],
                r["total_count"],
                *[
                    number(r[k])
                    for k in ("bias", "MAE", "RMS", "maximum_absolute_error")
                ],
                "Pass" if r["robust_passed"] == "True" else "Fail",
            ]
            for r in verdicts
        ],
    )
    (RESULTS / "verdict-table.md").write_text(content)
    content = table(
        [
            "Source",
            "I bucket",
            "Frozen cases passing transfer",
            "Range of RMS(Held) − RMS(Figiel)",
            "Largest maximum-error growth",
        ],
        [
            [
                LABELS[source],
                bucket,
                f"{sum(r['passed'] == 'True' for r in selected)} / {len(selected)}",
                f"{min(float(r['RMS_improvement']) for r in selected):.6g} to {max(float(r['RMS_improvement']) for r in selected):.6g}",
                number(max(float(r["maximum_error_growth"]) for r in selected)),
            ]
            for source in SOURCES
            for bucket in ("I<=3", "I>3")
            for selected in [
                [
                    r
                    for r in transfers
                    if r["source"] == source and r["I_bucket"] == bucket
                ]
            ]
        ],
    )
    (RESULTS / "transfer-table.md").write_text(content)
    content = table(
        [
            "Source",
            "I bucket",
            "Declaration",
            "Frozen RMS range",
            "Frozen maximum-error range",
            "Passing cases",
        ],
        [
            [
                LABELS[v["source"]],
                v["I_bucket"],
                v["declaration"],
                f"{min(float(r['RMS']) for r in selected):.6g} to {max(float(r['RMS']) for r in selected):.6g}",
                f"{min(float(r['maximum_absolute_error']) for r in selected):.6g} to {max(float(r['maximum_absolute_error']) for r in selected):.6g}",
                f"{v['passing_realizations']} / {v['frozen_realizations']}",
            ]
            for v in verdicts
            for selected in [
                [
                    r
                    for r in buckets
                    if r["source"] == v["source"]
                    and r["I_bucket"] == v["I_bucket"]
                    and r["declaration"] == v["declaration"]
                ]
            ]
        ],
    )
    (RESULTS / "sensitivity-table.md").write_text(content)
    reference_rows = [r for r in nominal if r["role"] != "mixture"] + binary
    content = table(
        [
            "Reference role",
            "Declaration",
            "n",
            "Phi RMS",
            "Maximum absolute Phi residual",
        ],
        [
            [
                role,
                packet,
                len(selected),
                number(
                    (sum(float(r["residual"]) ** 2 for r in selected) / len(selected))
                    ** 0.5
                ),
                number(max(abs(float(r["residual"])) for r in selected)),
            ]
            for role in ("CaCl2 endpoint", "NaCl anchor", "MgCl2 binary")
            for packet in ("Held", "Figiel")
            for selected in [
                [
                    r
                    for r in reference_rows
                    if r["role"] == role and r["declaration"] == packet
                ]
            ]
        ],
    )
    (RESULTS / "reference-table.md").write_text(content)
    content = table(
        [
            "Source",
            "Declaration",
            "I range, all measured rows (mol/kg)",
            "Mixture I range (mol/kg)",
            "NaCl molality range",
            "Divalent molality range",
        ],
        [
            [
                LABELS[source],
                "Both",
                f"{min(float(r['I_mol_kg']) for r in selected):.6g} to {max(float(r['I_mol_kg']) for r in selected):.6g}",
                f"{min(float(r['I_mol_kg']) for r in mixed):.6g} to {max(float(r['I_mol_kg']) for r in mixed):.6g}",
                f"{min(float(r['m_NaCl']) for r in mixed):.6g} to {max(float(r['m_NaCl']) for r in mixed):.6g}",
                f"{min(float(r['m_divalent']) for r in mixed):.6g} to {max(float(r['m_divalent']) for r in mixed):.6g}",
            ]
            for source in SOURCES
            for selected in [
                [
                    r
                    for r in nominal
                    if r["source"] == source and r["declaration"] == "Held"
                ]
            ]
            for mixed in [[r for r in selected if r["role"] == "mixture"]]
        ],
    )
    (RESULTS / "domain-table.md").write_text(content)
    omitted = [
        r for r in buckets if r["realization_id"] == "r1466" and r["I_bucket"] != "all"
    ]
    content = table(
        ["Rard source-flag omission", "Declaration", "n", "RMS", "Maximum", "Verdict"],
        [
            [
                r["I_bucket"],
                r["declaration"],
                r["total_count"],
                number(r["RMS"]),
                number(r["maximum_absolute_error"]),
                "Pass" if r["passed"] == "True" else "Fail",
            ]
            for r in omitted
        ],
    )
    (RESULTS / "omission-table.md").write_text(content)
    content = table(
        [
            "Source / component",
            "Held maximum magnitude ratio",
            "Figiel maximum magnitude ratio",
        ],
        [
            [
                label,
                *[
                    number(
                        max(
                            abs(float(r[key]))
                            for r in nominal
                            if r["source"] == source
                            and r["declaration"] == packet
                            and r[key] != ""
                        )
                    )
                    for packet in ("Held", "Figiel")
                ],
            ]
            for source, key, label in [
                (
                    SOURCES[1],
                    "relative_residual_over_scatter_0p1pct",
                    "Rard Phi residual / reported 0.1% scatter",
                ),
                (
                    SOURCES[1],
                    "relative_residual_over_scatter_0p25pct",
                    "Rard Phi residual / reported 0.25% scatter",
                ),
                (
                    SOURCES[1],
                    "relative_residual_over_concentration_0p13pct",
                    "Rard relative Phi residual / 0.13% concentration figure",
                ),
                (
                    SOURCES[1],
                    "relative_residual_over_Ca_reference_0p2pct",
                    "Rard relative Phi residual / 0.2% Ca reference figure",
                ),
                (
                    SOURCES[1],
                    "relative_residual_over_Ca_reference_0p4pct",
                    "Rard relative Phi residual / 0.4% Ca reference figure",
                ),
                (
                    SOURCES[2],
                    "residual_over_aw_dispersion",
                    "Dinane aw residual / reported 0.0006 dispersion",
                ),
                (
                    SOURCES[2],
                    "residual_over_Phi_dispersion",
                    "Dinane diagnostic Phi residual / reported 0.028 dispersion",
                ),
            ]
        ],
    )
    (RESULTS / "uncertainty-table.md").write_text(content)
    app = summary["ms2"]
    content = table(
        ["MS-2 formal inventory quantity", "Value", "Basis / limitation"],
        [
            [
                "Water mass per liter",
                f"{app['water_kg_per_L']:.6g} kg",
                application_basis["water_mass_basis"],
            ],
            [
                "Major-ion ionic strength",
                f"{app['major_ion_I_mol_kg']:.6g} mol/kg",
                "Raw formal inventory; no speciation or charge correction",
            ],
            [
                "Na/Mg/Ca chloride approximation",
                f"{app['chloride_approximation_I_mol_kg']:.6g} mol/kg",
                "Na + 3 Mg + 3 Ca on the same water basis",
            ],
            [
                "Charge imbalance",
                f"{app['charge_imbalance_mol_kg']:+.6g} mol charge/kg",
                "Retained, not reconciled",
            ],
            [
                "Source temperature",
                f"{application_basis['source_temperature_degC']} °C",
                "Conditioned extraction inlet undefined",
            ],
        ],
    )
    content += table(
        [
            "Ion",
            "Measured concentration (mg/L)",
            "Formal molality (mol/kg water)",
            "Charge",
        ],
        [
            [
                r["ion"],
                number(r["concentration_mg_L"]),
                number(r["molality_mol_kg"]),
                r["charge"],
            ]
            for r in read("ms2-formal-molalities.csv")
        ],
    )
    (RESULTS / "ms2-table.md").write_text(content)
    stage_a = json.loads(
        (HERE / "inputs/stage-a-verification-summary.json").read_text()
    )["verification"]
    content = table(
        ["Numerical verification", "Retained result"],
        [
            [
                "Stage A",
                f"{stage_a['rows']} checks; {stage_a['failed']} failures; independently Supported at 2aca2f65",
            ],
            [
                "Finite permutation maximum scaled difference",
                number(stage_a["permutation"]["scaled_max"]),
            ],
            [
                "Centered derivative-check maximum scaled difference",
                number(stage_a["difference"]["scaled_max"]),
            ],
            [
                "Native fugacity identity maximum difference",
                number(stage_a["fugacity"]["scaled_max"]),
            ],
            [
                "Stage B native states / unavailable",
                f"{summary['unique_native_states']} / {summary['native_failure_count']}",
            ],
            [
                "Repeated nominal states",
                f"{summary['repeated_states']}; max aw difference {summary['repeat_aw_max']}; max Phi difference {summary['repeat_Phi_max']}",
            ],
            [
                "Maximum native pressure error",
                number(summary["maximum_pressure_error_Pa"]) + " Pa",
            ],
            [
                "Minimum pressure-density slope",
                number(summary["minimum_pressure_density_derivative"]) + " Pa m³/mol",
            ],
            [
                "Maximum charge residual in mole fractions",
                number(summary["maximum_charge_residual"]),
            ],
            [
                "Maximum Robinson reference reconstruction sensitivity",
                number(summary["Robinson_maximum_reconstruction_sensitivity"])
                + f" Phi; fixed limit {source_checks['Robinson_reference_sensitivity_limit']}",
            ],
            [
                "Maximum water-mass conversion difference",
                number(
                    max(
                        abs(float(r["water_mass_conversion_delta_aw"])) for r in nominal
                    )
                )
                + " in derived aw; source and declaration masses kept distinct",
            ],
            [
                "Specification freeze commit",
                "`" + summary["specification_commit"] + "`",
            ],
            ["Calculation source commit", "`" + summary["engine_commit"] + "`"],
            ["Specification SHA256", "`" + summary["specification_sha256"] + "`"],
        ],
    )
    (RESULTS / "verification-table.md").write_text(content)
    plt.rcParams.update(
        {"font.size": 10, "axes.spines.top": False, "axes.spines.right": False}
    )
    for residual in (False, True):
        figure, axes = plt.subplots(3, 1, figsize=(8.2, 10), constrained_layout=True)
        for ax, source in zip(axes, SOURCES, strict=True):
            source_rows = [
                r for r in nominal if r["source"] == source and r["role"] == "mixture"
            ]
            quantity = source_rows[0]["quantity"]
            if not residual:
                observed = [r for r in source_rows if r["declaration"] == "Held"]
                ax.scatter(
                    [float(r["I_mol_kg"]) for r in observed],
                    [float(r["observed"]) for r in observed],
                    facecolors="none",
                    edgecolors="#333333",
                    s=25,
                    label="Observation",
                    linewidths=0.7,
                )
            else:
                limit = spec["primary_limits"][quantity]["maximum_absolute_error"]
                ax.axhspan(
                    -limit,
                    limit,
                    color="#dddddd",
                    alpha=0.55,
                    label="Maximum-error limit",
                )
                ax.axhline(0, color="#777777", linewidth=0.7)
            for packet, color, marker in [
                ("Held", "#0072B2", "+"),
                ("Figiel", "#D55E00", "x"),
            ]:
                selected = [r for r in source_rows if r["declaration"] == packet]
                ax.scatter(
                    [float(r["I_mol_kg"]) for r in selected],
                    [float(r["residual" if residual else "model"]) for r in selected],
                    color=color,
                    marker=marker,
                    s=20,
                    linewidths=0.8,
                    label=packet,
                )
            ax.axvline(3, color="#777777", linestyle="--", linewidth=0.8)
            ax.set(
                title=LABELS[source],
                xlabel="Ionic strength (mol/kg water)",
                ylabel=("Model − observation: " if residual else "")
                + (
                    "Φ (dimensionless)"
                    if quantity == "Phi"
                    else "Water activity (dimensionless)"
                ),
            )
            ax.legend(loc="best", frameon=False, ncol=3, fontsize=9)
        figure.savefig(
            FIGURES / ("primary-residuals.png" if residual else "water-properties.png"),
            dpi=180,
        )
        plt.close(figure)
    paths = [
        RESULTS / name
        for name in (
            "summary.json",
            "verdicts.csv",
            "nominal-residuals.csv",
            "bucket-metrics.csv",
            "transfer-tests.csv",
            "binary-reference-checks.csv",
            "ms2-formal-molalities.csv",
        )
    ]
    (FIGURES / "presentation-inputs.json").write_text(
        json.dumps(
            {
                "inputs": {
                    str(p.relative_to(HERE)): hashlib.sha256(p.read_bytes()).hexdigest()
                    for p in paths
                },
                "script_sha256": hashlib.sha256(
                    Path(__file__).read_bytes()
                ).hexdigest(),
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
