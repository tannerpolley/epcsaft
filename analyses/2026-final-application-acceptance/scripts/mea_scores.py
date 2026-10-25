"""Replay MEA's pCO2 and speciation scoring on the installed Engine wheel (#64).

Usage: mea_scores.py MEA_CHECKOUT. MEA pins wheel `b66c7b96`; this runs its unmodified
`generate_figure_data` and `score_observation_coverage` on the installed wheel instead, writing
into MEA_CHECKOUT (a disposable detached worktree). compare.py reads the outputs against the
committed pinned-wheel files.
"""
from __future__ import annotations

import hashlib
import importlib.metadata
import json
import os
import sys
from pathlib import Path
from urllib.parse import unquote, urlparse

BUNDLE = Path(sys.argv[1]).resolve() / "analyses/mea_parameter_bundle"
sys.path[:0] = [str(BUNDLE / "scripts"), str(BUNDLE.parents[1] / "src")]
os.chdir(BUNDLE)
os.environ["MEA_BOUNDED_RUN"] = "1"
import shared_evaluation as shared  # noqa: E402

url = json.loads(importlib.metadata.distribution("epcsaft").read_text("direct_url.json"))["url"]
# Before the scripts import it by name, so their verify_wheel and records carry the installed wheel.
shared.ENGINE_WHEEL_SHA256 = hashlib.sha256(Path(unquote(urlparse(url).path)).read_bytes()).hexdigest()
import generate_figure_data  # noqa: E402
import score_observation_coverage  # noqa: E402

sys.argv = ["generate_figure_data.py", "--state-timeout-s", "180", "--overall-timeout-s", "3600"]
generate_figure_data.main()
sys.argv = ["score_observation_coverage.py"]
score_observation_coverage.main()
