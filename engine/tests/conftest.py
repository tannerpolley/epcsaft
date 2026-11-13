"""Small stdout reporter for numerical operands retained by ordinary pytest runs."""

from __future__ import annotations

import json
from typing import Any

import pytest


class Report:
    """One JSON line per call; assertions remain the sole test verdict owner."""

    @staticmethod
    def observe(comparison_id: str, **values: Any) -> None:
        print("EPCSAFT_COMPARE_JSON", json.dumps({"comparison_id": comparison_id, **values}, separators=(",", ":")))

    def compare(self, comparison_id: str, model: Any, reference: Any, **values: Any) -> None:
        self.observe(comparison_id, model=model, reference=reference, **values)


@pytest.fixture
def evidence() -> Report:
    return Report()
