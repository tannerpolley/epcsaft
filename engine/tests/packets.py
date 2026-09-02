from pathlib import Path

PACKETS = Path(__file__).resolve().parents[2] / "data" / "packets"
ORACLE = Path(__file__).resolve().parent / "oracle" / "oracle.json"


def parameters_json(packet_id: str, version: int) -> Path:
    return PACKETS / packet_id / str(version) / "parameters" / "parameters.json"
