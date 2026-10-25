"""Check retained measurement coverage, native admission and bucket statistics."""

import csv
import hashlib
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]


def read(name):
    with (HERE / "results" / name).open() as stream:
        return list(csv.DictReader(stream))


def main():
    summary = json.loads((HERE / "results/summary.json").read_text())
    assert (
        summary["script_sha256"]
        == hashlib.sha256((HERE / "scripts/compare.py").read_bytes()).hexdigest()
    )
    for item in json.loads((HERE / "inputs/input-files.json").read_text())["files"]:
        assert (
            hashlib.sha256((HERE / "inputs" / item["path"]).read_bytes()).hexdigest()
            == item["sha256"]
        )
    nominal, predictions = read("nominal-residuals.csv"), read("native-predictions.csv")
    by_id = {r["prediction_id"]: r for r in predictions}
    assert len(by_id) == len(predictions)
    expected = {
        "robinson-bower-1966": (53, 16),
        "rard-miller-1987": (287, 0),
        "dinane-mounir-2003": (19, 0),
    }
    for source, (mixtures, endpoints) in expected.items():
        for packet in ("Held", "Figiel"):
            rows = [
                r
                for r in nominal
                if r["source"] == source and r["declaration"] == packet
            ]
            assert sum(r["role"] == "mixture" for r in rows) == mixtures
            assert sum(r["role"] != "mixture" for r in rows) == endpoints
    for row in nominal:
        prediction = by_id[row["prediction_id"]]
        assert prediction["declaration"] == row["declaration"]
        assert float(prediction["m_NaCl"]) == float(row["m_NaCl"])
        assert float(prediction["m_divalent"]) == float(row["m_divalent"])
        assert row["quantity"] == (
            "aw" if row["source"] == "dinane-mounir-2003" else "Phi"
        )
        if row["available"] == "True":
            assert float(row["model"]) == float(prediction[row["quantity"]])
            assert math.isclose(
                float(row["model"]) - float(row["observed"]),
                float(row["residual"]),
                rel_tol=0,
                abs_tol=1e-15,
            )
    for row in predictions:
        if row["available"] == "True":
            assert (
                float(row["pressure_density_derivative"]) > 0
                and float(row["density_mol_m3"]) > 0
            )
            assert row["cause"] == "" and row["osmotic_cause"] == ""
        if row["repeat_aw"]:
            assert abs(float(row["aw"]) - float(row["repeat_aw"])) <= 1e-10
            assert abs(float(row["Phi"]) - float(row["repeat_Phi"])) <= 1e-8
    metrics = read("bucket-metrics.csv")
    for metric in [r for r in metrics if r["realization_id"] == "r0000"]:
        rows = [
            r
            for r in nominal
            if r["role"] == "mixture"
            and r["source"] == metric["source"]
            and r["declaration"] == metric["declaration"]
            and (metric["I_bucket"] == "all" or r["I_bucket"] == metric["I_bucket"])
        ]
        assert len(rows) == int(metric["total_count"])
        errors = [float(r["residual"]) for r in rows if r["available"] == "True"]
        if errors:
            assert math.isclose(
                math.sqrt(math.fsum(v * v for v in errors) / len(errors)),
                float(metric["RMS"]),
                rel_tol=1e-14,
            )
            assert max(map(abs, errors)) == float(metric["maximum_absolute_error"])
    coverage, transfers = read("realization-coverage.csv"), read("transfer-tests.csv")
    assert len(coverage) == len(transfers) == summary["realizations"] * 6
    assert len(
        {(r["realization_id"], r["source"], r["I_bucket"]) for r in coverage}
    ) == len(coverage)
    assert len(read("binary-reference-checks.csv")) == (10 + 4) * 2
    pure = read("pure-reference-checks.csv")
    assert len(pure) == 2 and all(
        float(r["aw"]) == 1 and r["osmotic_cause"] == "OsmoticUndefinedPureSolvent"
        for r in pure
    )
    config = json.loads((HERE.parent / "_quarto.yml").read_text())
    notebook = HERE.name + "/notebook.qmd"
    assert config["project"]["render"].count(notebook) == 1
    assert (
        sum(r["href"] == notebook for r in config["website"]["sidebar"]["contents"])
        == 1
    )
    assert (HERE.parent / "README.md").read_text().count("](" + HERE.name + "/)") == 1
    print(
        json.dumps(
            {
                "measurement_coverage": "complete",
                "native_states": len(predictions),
                "unavailable_states": summary["native_failure_count"],
                "frozen_bucket_tests": len(coverage),
                "nominal_metrics": "recomputed from retained residuals",
                "Quarto_membership": "one per required list",
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
