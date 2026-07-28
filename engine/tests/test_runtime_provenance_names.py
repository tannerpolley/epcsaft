from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import subprocess
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


def _fixture(tmp_path: Path, runtime_text: str, bundle_text: str = "") -> Path:
    runtime = tmp_path / "src" / "epcsaft"
    bundles = runtime / "data" / "bundles"
    bundles.mkdir(parents=True)
    (runtime / "runtime.py").write_text(runtime_text, encoding="utf-8")
    (bundles / "provenance_names.toml").write_text(
        'tokens = ["exampleauthor"]\n',
        encoding="utf-8",
    )
    (bundles / "bundle.toml").write_text(bundle_text, encoding="utf-8")
    return tmp_path


def test_scan_rejects_provenance_name_in_runtime_source(tmp_path: Path) -> None:
    checker = _checker()
    root = _fixture(tmp_path, 'selector = "ExampleAuthor-formulation"\n')

    violations = checker.scan_runtime_source(root)

    assert len(violations) == 1
    assert violations[0].path.as_posix() == "src/epcsaft/runtime.py"
    assert violations[0].line == 1
    assert violations[0].token == "exampleauthor"


def test_scan_excludes_bundle_provenance(tmp_path: Path) -> None:
    checker = _checker()
    root = _fixture(
        tmp_path,
        "selector = 'generic-formulation'\n",
        'citation = "ExampleAuthor (2025)"\n',
    )

    assert checker.scan_runtime_source(root) == ()


def test_repository_runtime_source_has_no_provenance_names() -> None:
    checker = _checker()

    assert checker.scan_runtime_source(REPO_ROOT) == ()


def test_cli_reports_file_line_and_token(tmp_path: Path) -> None:
    root = _fixture(tmp_path, "# ExampleAuthor correlation\n")

    result = subprocess.run(
        [sys.executable, str(SCRIPT), "--repo-root", str(root)],
        check=False,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 1
    assert (
        "src/epcsaft/runtime.py:1: prohibited provenance name: exampleauthor"
        in result.stdout
    )
