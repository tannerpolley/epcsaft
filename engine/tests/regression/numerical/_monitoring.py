from __future__ import annotations

import sys
from collections.abc import Callable
from typing import Any


def call_without_python_frames(
    function: Callable[..., Any], *args: object
) -> tuple[Any, tuple[str, ...]]:
    frames: list[str] = []
    tool_id = sys.monitoring.OPTIMIZER_ID

    def record(code: object, *_: object) -> None:
        frames.append(str(getattr(code, "co_qualname", code)))

    sys.monitoring.use_tool_id(tool_id, "native-regression-test")
    try:
        sys.monitoring.register_callback(tool_id, sys.monitoring.events.PY_START, record)
        sys.monitoring.register_callback(tool_id, sys.monitoring.events.PY_RESUME, record)
        sys.monitoring.set_events(
            tool_id, sys.monitoring.events.PY_START | sys.monitoring.events.PY_RESUME
        )
        result = function(*args)
    finally:
        sys.monitoring.set_events(tool_id, 0)
        sys.monitoring.register_callback(tool_id, sys.monitoring.events.PY_START, None)
        sys.monitoring.register_callback(tool_id, sys.monitoring.events.PY_RESUME, None)
        sys.monitoring.free_tool_id(tool_id)
    return result, tuple(frames)
