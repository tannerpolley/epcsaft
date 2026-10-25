"""Generate exact plotted tables for the Khudaida Figures 3-7 comparisons."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

CAMPAIGN_ROOT = Path(__file__).resolve().parents[1]
RESULTS = CAMPAIGN_ROOT / "results" / "data" / "khudaida-2026-figures-3-to-7-held2.csv"
FIGURES = (3, 4, 5, 6, 7)
SOURCE_FILES = (
    "source-points.csv",
    "paper-epcsaft-points.csv",
    "feed-compositions.csv",
)
SOURCE_SHA256 = {
    3: (
        "92398d0c56ff33df805dc090f2fe94e6810606a685f484afcefd5301e2ae9f20",
        "2c3cd8e69e6788d7c211dc20aa02941d32073e155cc3ffc39816e083d326ba65",
        "fa87cb1c1174422792d36d62f008afd3d73f2dcf277b8ca5f63fc79d6cf842f7",
    ),
    4: (
        "e82bb4435034ded02fb51d28e2a6d69057c5bfe010bebb806e35c1da26411741",
        "5c8908fc4076acd62b46962befd767b439c5df3057bfd53e271023395f1df2d6",
        "1e2311cfcc2c268758158a4f088138ac70de942561d05866b3391aec058e42d5",
    ),
    5: (
        "19fd9ee087bbd8ced9b0af227483b17b89c3ad3a99076057ad0ba2382058a6a6",
        "072e0470d03c50a14f0c5341af394b93628c9b53867fda040a61ee4e8067c25e",
        "a3c3381bf56860521763c88986fe3f8c71e3074459807c6d1b0f2aed1290f7ff",
    ),
    6: (
        "a4d5fff850d7246c4fa71f385ce44b4a28775c5d08e5da2e851ea055c6beaac3",
        "fa16d867e497809dcbae185cf2a70122981f0428e23e80633aa92256223909ac",
        "7f90fd9f369e53fa1067e887f08ed36ad2214333cf93ba339d57d1fd936d6c46",
    ),
    7: (
        "cb28feb61b0afc9f8ec7a0a6a8f4f3f1278a53d162a5293761892c09fb6ab48b",
        "35e816d805f63f53a9a5b277e0de1426c5733c70ed616197055c489d69ec5bdc",
        "ffb836bed66fd5b9bfc751c5849b8217036df946061d70556ee00e8be270e13c",
    ),
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def comparison_inputs(figure_id: int) -> list[dict[str, str]]:
    figure_root = CAMPAIGN_ROOT / "figures" / f"figure-{figure_id:02d}"
    inputs: list[dict[str, str]] = []
    for filename, expected in zip(
        SOURCE_FILES, SOURCE_SHA256[figure_id], strict=True
    ):
        path = figure_root / "source" / filename
        actual = sha256(path)
        if actual != expected:
            raise AssertionError(
                f"Figure {figure_id} comparison source changed: {filename}"
            )
        inputs.append({
            "path": str(path.relative_to(figure_root)),
            "sha256": actual,
        })
    return inputs


def salt_free(explicit: list[float]) -> tuple[float, float, float]:
    solvent_total = sum(explicit[:3])
    if solvent_total <= 0.0:
        raise ValueError("solvent total must be positive")
    return tuple(value / solvent_total for value in explicit[:3])


def active_rows(campaign: pd.DataFrame, figure_id: int) -> list[dict[str, object]]:
    selected = campaign[campaign["figure_id"].eq(figure_id)]
    rows: list[dict[str, object]] = []
    for record in selected.itertuples(index=False):
        if record.case_decision != "PASS" or not isinstance(record.result_json, str):
            continue
        result = json.loads(record.result_json)
        phases = sorted(
            result["phases"],
            key=lambda phase: phase["mole_fractions"][0],
            reverse=True,
        )
        if len(phases) != 2:
            continue
        for phase_name, phase in zip(("aqueous", "organic"), phases, strict=True):
            water, ethanol, isobutanol = salt_free(phase["mole_fractions"])
            rows.append(
                {
                    "dataset": "active_epcsaft",
                    "point_id": int(record.point_id),
                    "phase": phase_name,
                    "x_water_salt_free": water,
                    "x_ethanol_salt_free": ethanol,
                    "x_isobutanol_salt_free": isobutanol,
                    "input_coordinate_sum": water + ethanol + isobutanol,
                    "source": "exploratory installed-artifact result",
                    "case_decision": record.case_decision,
                    "parameter_domain_status": record.parameter_domain_status,
                }
            )
    return rows


def literature_rows(figure_id: int) -> list[dict[str, object]]:
    figure_root = CAMPAIGN_ROOT / "figures" / f"figure-{figure_id:02d}"
    comparison_inputs(figure_id)
    source = pd.read_csv(figure_root / "source" / "source-points.csv")
    experimental = source[source["dataset"].eq("experimental_tieline")].copy()
    solvent_total = experimental[["x_water", "x_ethanol", "x_isobutanol"]].sum(axis=1)
    for component in ("water", "ethanol", "isobutanol"):
        experimental[f"x_{component}_salt_free"] = (
            experimental[f"x_{component}"] / solvent_total
        )
    experimental["dataset"] = "experimental"
    experimental["point_id"] = experimental["tie_line"].astype(int)

    paper = pd.read_csv(figure_root / "source" / "paper-epcsaft-points.csv").assign(
        dataset="paper_epcsaft", phase="organic"
    )
    feeds = pd.read_csv(figure_root / "source" / "feed-compositions.csv").assign(
        dataset="feed", phase="feed"
    )
    columns = [
        "dataset",
        "point_id",
        "phase",
        "x_water_salt_free",
        "x_ethanol_salt_free",
        "x_isobutanol_salt_free",
    ]
    combined = pd.concat((experimental[columns], paper[columns], feeds[columns]), ignore_index=True)
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
    table = "Table 3" if figure_id in (3, 4) else "Table 4"
    combined["source"] = combined["dataset"].map(
        {
            "experimental": f"Khudaida {table}",
            "paper_epcsaft": f"Khudaida Figure {figure_id} digitization",
            "feed": f"Khudaida Figure {figure_id} feed digitization",
        }
    )
    combined["case_decision"] = "LITERATURE_OR_INPUT"
    combined["parameter_domain_status"] = "NOT_APPLICABLE"
    return combined.to_dict(orient="records")


def generate(figure_ids: tuple[int, ...]) -> None:
    campaign = pd.read_csv(RESULTS, dtype={"case_decision": "string"})
    for figure_id in figure_ids:
        expected = int(
            campaign.loc[campaign["figure_id"].eq(figure_id), "point_id"].max()
        )
        point_ids = campaign.loc[campaign["figure_id"].eq(figure_id), "point_id"].tolist()
        if point_ids != list(range(1, expected + 1)):
            raise AssertionError(f"Figure {figure_id} campaign rows are incomplete")
        plotted = pd.DataFrame(
            [*literature_rows(figure_id), *active_rows(campaign, figure_id)]
        ).sort_values(["dataset", "point_id", "phase"])
        coordinates = [
            "x_water_salt_free",
            "x_ethanol_salt_free",
            "x_isobutanol_salt_free",
        ]
        if not plotted[coordinates].sum(axis=1).sub(1.0).abs().le(1.0e-12).all():
            raise AssertionError("ternary coordinates do not normalize")
        output = (
            CAMPAIGN_ROOT
            / "figures"
            / f"figure-{figure_id:02d}"
            / "output"
            / f"khudaida-figure-{figure_id:02d}-active-comparison.csv"
        )
        output.parent.mkdir(parents=True, exist_ok=True)
        plotted.reset_index(drop=True).to_csv(output, index=False, lineterminator="\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--figures", type=int, nargs="+", choices=FIGURES, default=FIGURES)
    args = parser.parse_args()
    generate(tuple(dict.fromkeys(args.figures)))


if __name__ == "__main__":
    main()
