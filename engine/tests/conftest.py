"""Small stdout reporter for numerical operands retained by ordinary pytest runs."""

from __future__ import annotations

import json
import math
from typing import Any

import pytest


class Report:
    @staticmethod
    def _emit(kind: str, comparison_id: str, **values: Any) -> None:
        print(kind, json.dumps({"comparison_id": comparison_id, **values}, separators=(",", ":")))

    def observe(self, comparison_id: str, **values: Any) -> None:
        self._emit("EPCSAFT_COMPARE_JSON", comparison_id, judgment="Observation only", **values)

    def compare(self, comparison_id: str, *, criterion: dict[str, Any], **values: Any) -> None:
        verdict: bool | None = None
        status = "valid"
        if not isinstance(criterion, dict):
            status = "invalid criterion"
        else:
            try:
                model = float(values["model"])
                reference = float(values["reference"])
                if not all(math.isfinite(value) for value in (model, reference)):
                    raise ValueError
            except (KeyError, TypeError, ValueError):
                status = "invalid operands or criterion"
            if status == "valid" and criterion.get("rule") == "symmetric_close":
                try:
                    atol = float(criterion.get("atol", 0.0))
                    rtol = float(criterion.get("rtol", 0.0))
                    if not all(math.isfinite(value) for value in (atol, rtol)) or atol < 0 or rtol < 0:
                        raise ValueError
                    verdict = math.isclose(model, reference, abs_tol=atol, rel_tol=rtol)
                except (TypeError, ValueError):
                    status = "invalid operands or criterion"
            elif status == "valid" and type(criterion.get("verdict")) is bool:
                verdict = criterion["verdict"]
            elif status == "valid":
                status = "unavailable criterion"
        self._emit(
            "EPCSAFT_COMPARE_JSON",
            comparison_id,
            criterion=criterion,
            evidence_status=status,
            judgment=("Criterion satisfied" if verdict is True else "Criterion failed" if verdict is False else "Evidence unavailable"),
            **values,
        )


@pytest.fixture
def evidence() -> Report:
    """Emit useful JSON lines; assertions remain the sole test verdict owner."""
    return Report()
