from importlib import metadata
import os
from pathlib import Path
import re

import nanobind

Path("licenses").mkdir(exist_ok=True)
notice = Path("NOTICE").read_text(encoding="utf-8")
if os.environ["target_platform"] in {"linux-64", "win-64"}:
    notice = re.sub(r"(?m)^\* SuiteSparse 7[^\n]*(?:\n  [^\n]*)*\n", "", notice)
Path("licenses/NOTICE").write_text(notice, encoding="utf-8")

dist = metadata.distribution("nanobind")
license_path = next(p for p in dist.files if str(p).endswith("/licenses/LICENSE"))
Path("licenses/nanobind-LICENSE").write_text(dist.locate_file(license_path).read_text())
header = Path(nanobind.__file__).parent / "ext/robin_map/include/tsl/robin_map.h"
Path("licenses/robin-map-LICENSE").write_text(header.read_text().split("*/", 1)[0] + "*/\n")
