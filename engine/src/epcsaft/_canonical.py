from __future__ import annotations

import hashlib
import json
import unicodedata

_TAG = "epcsaft-parameters-v1"


def _string(value: str) -> str:
    return unicodedata.normalize("NFC", value.replace("\r\n", "\n").replace("\r", "\n"))


def _normalize(value: object) -> object:
    if isinstance(value, str):
        return _string(value)
    if isinstance(value, list):
        return [_normalize(item) for item in value]
    if isinstance(value, dict):
        return {key: _normalize(item) for key, item in value.items()}
    return value


def _canonical_json(value: object) -> str:
    return json.dumps(
        value,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    )


def _json(value: object) -> str:
    return _canonical_json(_normalize(value))


def _fingerprint_json(value: str) -> str:
    return f"sha256:{hashlib.sha256(value.encode('utf-8')).hexdigest()}"


def _document_fingerprint(value: object) -> str:
    return _fingerprint_json(_json({"tag": _TAG, "parameters": value}))


def _selected_document_fingerprint(
    document_fingerprint: str,
    component_ids: tuple[str, ...],
) -> str:
    return _fingerprint_json(
        _json(
            {
                "tag": _TAG,
                "document_fingerprint": document_fingerprint,
                "component_ids": list(component_ids),
            }
        )
    )
