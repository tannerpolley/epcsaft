from __future__ import annotations

import json
from collections.abc import Mapping
from os import PathLike
from pathlib import Path

from .units import _error


def _load_json(path: str | PathLike[str]) -> Mapping[str, object]:
    try:
        document_path = Path(path)
    except TypeError as exc:
        raise _error("invalid-value", field="path", value=path) from exc
    if document_path.suffix.lower() != ".json" or not document_path.is_file():
        raise _error("missing-record", field="parameters_json", value=document_path)
    try:
        value = json.loads(document_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise _error("invalid-schema", field="parameters_json", value=document_path) from exc
    if not isinstance(value, dict) or any(type(key) is not str for key in value):
        raise _error("invalid-schema", field="parameters_json", expected="object")
    return value


__all__ = ("_load_json",)
