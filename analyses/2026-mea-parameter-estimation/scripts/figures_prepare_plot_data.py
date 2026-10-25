from __future__ import annotations

import csv
import json
from pathlib import Path

CAMPAIGN = Path(__file__).resolve().parents[1]
OUTPUT = Path(__file__).resolve().parents[1] / "results/analysis/figures/output"


def write(name: str, rows: list[dict[str, object]]) -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    with (OUTPUT / name).open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    retained = json.loads((CAMPAIGN / "results/selection-ladder.json").read_text())
    neutral = retained["neutral_association"]
    write(
        "neutral-association-selection.csv",
        [
            {
                "scheme": item["scheme"].upper(),
                "grouped_mean_score": item["grouped_mean_score"],
                "grouped_standard_error": item["grouped_standard_error"],
                "phase_success": item["phase_success"],
                "selected": item["scheme"] == neutral["selected_scheme"],
            }
            for item in neutral["candidates"]
        ],
    )
    rows = retained["kiepe_physical_prerequisite"]["records"]
    write(
        "kiepe-pressure-parity.csv",
        [
            {
                "identity": item["identity"],
                "temperature_k": item["temperature_k"],
                "observed_pressure_pa": item["observed_pressure_pa"],
                "predicted_pressure_pa": item["predicted_pressure_pa"],
            }
            for item in rows
            if item["failure"] is None
        ],
    )


if __name__ == "__main__":
    main()
