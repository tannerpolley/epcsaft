from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = REPO_ROOT / "scripts" / "check_runtime_provenance_names.py"


def _checker():
    spec = spec_from_file_location("check_runtime_provenance_names", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_repository_runtime_source_has_no_provenance_names() -> None:
    checker = _checker()

    assert checker.scan_runtime_source(REPO_ROOT) == ()
