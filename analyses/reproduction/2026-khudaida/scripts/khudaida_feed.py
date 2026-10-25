"""Map one retained Khudaida feed row to explicit species."""

from __future__ import annotations

import csv
import math
from pathlib import Path


def khudaida_feed(source: Path, point_id: int) -> tuple[float, ...]:
    with source.open(newline="", encoding="utf-8") as stream:
        row = next(
            item
            for item in csv.DictReader(stream)
            if item["dataset"] == "feed" and int(item["point_id"]) == point_id
        )
    formula = tuple(
        float(row[name])
        for name in ("x_water", "x_ethanol", "x_isobutanol", "x_nacl")
    )
    explicit = (*formula[:3], formula[3], formula[3])
    total = math.fsum(explicit)
    return tuple(value / total for value in explicit)
