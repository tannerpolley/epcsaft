import csv, hashlib
from pathlib import Path
leaf = Path(__file__).resolve().parents[1]
groups = {}
for path in sorted((leaf / "results").rglob("*.csv")):
    if path.stat().st_size >= 1024: groups.setdefault(hashlib.sha256(path.read_bytes()).hexdigest(), []).append(path)
for paths in groups.values():
    if len(paths) < 2 and b"calculation_id" not in paths[0].read_bytes().splitlines()[0]: continue
    for path in paths:
        rows = list(csv.reader(path.open())); column = rows[0].index("calculation_id") if "calculation_id" in rows[0] else len(rows[0])
        with path.open("w", newline="") as stream: csv.writer(stream, lineterminator="\n").writerows([row[:column] + ["calculation_id" if i == 0 else str(path.relative_to(leaf))] for i, row in enumerate(rows)])
