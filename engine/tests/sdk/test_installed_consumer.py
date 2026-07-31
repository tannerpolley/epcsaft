from __future__ import annotations

from ._support import *


def test_native_sdk_header_is_installed_package_data() -> None:
    header = (
        __import__("importlib.resources").resources.files(epcsaft)
        / "include"
        / "epcsaft"
        / "native_sdk_v1.h"
    )
    assert header.is_file()


def test_native_sdk_installed_cpp_consumer_reads_sensitivity_tail(
    tmp_path: Path,
) -> None:
    artifact = tmp_path / "electrolyte.json"
    export_native_model(_held_water_ionization_model(), artifact)
    neutral_artifact = tmp_path / "methane-ethane.json"
    export_native_model(_mixture_model(), neutral_artifact)
    source_dir = Path(__file__).parents[1] / "native_sdk_consumer"
    build_dir = tmp_path / "consumer-build"
    package_dir = Path(native_core.__file__).parent / "cmake"
    configured = subprocess.run(
        (
            "cmake",
            "-S",
            str(source_dir),
            "-B",
            str(build_dir),
            f"-Depcsaft_DIR={package_dir}",
        ),
        capture_output=True,
        text=True,
        check=False,
    )
    assert configured.returncode == 0, configured.stderr
    built = subprocess.run(
        ("cmake", "--build", str(build_dir), "--parallel"),
        capture_output=True,
        text=True,
        check=False,
    )
    assert built.returncode == 0, built.stderr
    consumed = subprocess.run(
        (
            str(build_dir / "epcsaft_native_sdk_consumer"),
            str(artifact),
            str(neutral_artifact),
        ),
        capture_output=True,
        text=True,
        check=False,
    )
    assert consumed.returncode == 0, consumed.stderr
    assert "component_count=3" in consumed.stdout
    assert "reference_branch=liquid" in consumed.stdout
