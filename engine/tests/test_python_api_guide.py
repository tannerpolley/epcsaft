from __future__ import annotations

import inspect
from pathlib import Path
import re

import epcsaft
import pytest
from epcsaft import EPCSAFT, ParameterBundle


ROOT = Path(__file__).parents[1]
GUIDE = ROOT / "docs" / "python-api.md"


def test_python_api_guide_matches_the_installed_public_contract(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    text = GUIDE.read_text(encoding="utf-8")

    catalog_root = Path(epcsaft.__file__).parent / "data" / "bundles"
    for bundle_path in sorted(path for path in catalog_root.iterdir() if path.is_dir()):
        bundle = ParameterBundle.from_catalog(bundle_path.name, version=1)
        row = next(line for line in text.splitlines() if f"`{bundle.bundle_id}`" in line)
        assert "| 1 |" in row
        assert all(f"`{component.component_id}`" in row for component in bundle.components)

    for public_method in (
        ParameterBundle.from_catalog,
        ParameterBundle.to_path,
        ParameterBundle.select,
        EPCSAFT.evaluate,
        EPCSAFT.reference_state,
    ):
        documentation = inspect.getdoc(public_method)
        assert documentation is not None
        assert len(documentation.splitlines()) >= 2

    blocks = re.findall(r"```python\n(.*?)```", text, flags=re.DOTALL)
    assert len(blocks) == 4
    for index, block in enumerate(blocks):
        example_root = tmp_path / f"example-{index + 1}"
        example_root.mkdir()
        monkeypatch.chdir(example_root)
        namespace = {"__name__": f"python_api_guide_{index}"}
        exec(compile(block, f"{GUIDE}#example-{index + 1}", "exec"), namespace)
