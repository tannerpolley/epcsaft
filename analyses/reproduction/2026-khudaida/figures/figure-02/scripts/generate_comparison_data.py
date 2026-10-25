"""Generate the exact plotted table for the active Figure 2 comparison."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

FIGURE_ROOT = Path(__file__).resolve().parents[1]
CAMPAIGN_ROOT = FIGURE_ROOT.parents[1]
RESULTS = CAMPAIGN_ROOT / "results" / "data" / "khudaida-2026-figure-2-held2.csv"
SOURCE_POINTS = FIGURE_ROOT / "source" / "source-points.csv"
PAPER_POINTS = FIGURE_ROOT / "source" / "paper-epcsaft-points.csv"
FEEDS = FIGURE_ROOT / "source" / "feed-compositions.csv"
OUTPUT = FIGURE_ROOT / "output" / "khudaida-figure-02-active-comparison.csv"


def salt_free(explicit: list[float]) -> tuple[float, float, float]:
    solvent_total = sum(explicit[:3])
    if solvent_total <= 0.0:
        raise ValueError("solvent total must be positive")
    return tuple(value / solvent_total for value in explicit[:3])


def load_campaign() -> tuple[pd.DataFrame, dict[int, int]]:
    campaign = pd.read_csv(RESULTS, dtype={"case_decision": "string"})
    if campaign["point_id"].tolist() != list(range(1, 9)):
        raise AssertionError("active comparison requires Figure 2 points 1-8")
    if campaign["experimental_tie_line_id"].tolist() != list(range(8, 0, -1)):
        raise AssertionError("active comparison requires the Figure 2 reverse pairing")
    if not campaign["case_decision"].eq("PASS").all():
        raise AssertionError("active comparison requires eight accepted rows")
    mapping = dict(
        zip(
            campaign["point_id"].astype(int),
            campaign["experimental_tie_line_id"].astype(int),
            strict=True,
        )
    )
    return campaign, mapping


def active_rows(campaign: pd.DataFrame) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for record in campaign.itertuples(index=False):
        result = json.loads(record.result_json)
        phases = sorted(
            result["phases"],
            key=lambda phase: phase["mole_fractions"][0],
            reverse=True,
        )
        if len(phases) != 2:
            raise AssertionError("active comparison requires two phases per point")
        for phase_name, phase in zip(("aqueous", "organic"), phases, strict=True):
            water, ethanol, isobutanol = salt_free(phase["mole_fractions"])
            rows.append(
                {
                    "dataset": "active_epcsaft",
                    "point_id": int(record.experimental_tie_line_id),
                    "feed_point_id": int(record.point_id),
                    "phase": phase_name,
                    "x_water_salt_free": water,
                    "x_ethanol_salt_free": ethanol,
                    "x_isobutanol_salt_free": isobutanol,
                    "input_coordinate_sum": water + ethanol + isobutanol,
                    "source": "installed-artifact result",
                }
            )
    return rows


def literature_rows(mapping: dict[int, int]) -> list[dict[str, object]]:
    source = pd.read_csv(SOURCE_POINTS)
    experimental = source[source["dataset"].eq("experimental_tieline")].copy()
    solvent_total = experimental[["x_water", "x_ethanol", "x_isobutanol"]].sum(axis=1)
    for component in ("water", "ethanol", "isobutanol"):
        experimental[f"x_{component}_salt_free"] = (
            experimental[f"x_{component}"] / solvent_total
        )
    experimental["dataset"] = "experimental"
    experimental["point_id"] = experimental["tie_line"].astype(int)
    experimental["feed_point_id"] = ""

    paper = pd.read_csv(PAPER_POINTS).assign(
        dataset="paper_epcsaft",
        phase="organic",
        feed_point_id="",
    )
    feeds = pd.read_csv(FEEDS).assign(dataset="feed", phase="feed")
    feeds["feed_point_id"] = feeds["point_id"].astype(int)
    feeds["point_id"] = feeds["feed_point_id"].map(mapping)
    columns = [
        "dataset",
        "point_id",
        "feed_point_id",
        "phase",
        "x_water_salt_free",
        "x_ethanol_salt_free",
        "x_isobutanol_salt_free",
    ]
    frames = [experimental[columns], paper[columns], feeds[columns]]
    combined = pd.concat(frames, ignore_index=True)
    coordinates = [
        "x_water_salt_free",
        "x_ethanol_salt_free",
        "x_isobutanol_salt_free",
    ]
    combined["input_coordinate_sum"] = combined[coordinates].sum(axis=1)
    if not combined["input_coordinate_sum"].gt(0.0).all():
        raise AssertionError("literature ternary coordinate sum must be positive")
    combined[coordinates] = combined[coordinates].div(
        combined["input_coordinate_sum"], axis=0
    )
    combined["source"] = combined["dataset"].map(
        {
            "experimental": "Khudaida Table 3",
            "paper_epcsaft": "Khudaida Figure 2 digitization",
            "feed": "Khudaida Figure 2 feed digitization",
        }
    )
    return combined.to_dict(orient="records")


def main() -> None:
    campaign, mapping = load_campaign()
    plotted = pd.DataFrame([*literature_rows(mapping), *active_rows(campaign)])
    plotted = plotted.sort_values(["dataset", "point_id", "phase"]).reset_index(
        drop=True
    )
    coordinates = [
        "x_water_salt_free",
        "x_ethanol_salt_free",
        "x_isobutanol_salt_free",
    ]
    if not plotted[coordinates].sum(axis=1).sub(1.0).abs().le(1.0e-12).all():
        raise AssertionError("ternary coordinates do not normalize")
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    plotted.to_csv(OUTPUT, index=False, lineterminator="\n")


if __name__ == "__main__":
    main()
