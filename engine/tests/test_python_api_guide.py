from __future__ import annotations

from epcsaft.parameters import _ParameterBundle

import inspect
from pathlib import Path
import re

import epcsaft
from epcsaft import Mixture, Parameters, State


ROOT = Path(__file__).parents[1]
GUIDE = ROOT / "docs" / "python-api.md"


def test_python_api_guide_matches_the_installed_public_contract() -> None:
    text = GUIDE.read_text(encoding="utf-8")

    catalog_root = Path(epcsaft.__file__).parent / "data" / "bundles"
    for bundle_path in sorted(path for path in catalog_root.iterdir() if path.is_dir()):
        bundle = _ParameterBundle.from_catalog(bundle_path.name, version=1)
        row = next(line for line in text.splitlines() if f"`{bundle.bundle_id}`" in line)
        assert "| 1 |" in row
        assert all(f"`{component.component_id}`" in row for component in bundle.components)

    for public_method in (
        Parameters.from_catalog,
        Parameters.from_bundle,
        Parameters.from_records,
        Parameters.create_bundle_template,
        Parameters.to_bundle,
        Mixture.state,
        Mixture.reference,
        State.ares,
        State.gres,
    ):
        documentation = inspect.getdoc(public_method)
        assert documentation is not None

    blocks = re.findall(r"```python\n(.*?)```", text, flags=re.DOTALL)
    assert blocks
    for index, block in enumerate(blocks):
        compile(block, f"{GUIDE}#example-{index + 1}", "exec")
