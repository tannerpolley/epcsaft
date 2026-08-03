from __future__ import annotations

import hashlib
import json

from ._support import *
from .test_neutral_reference import _generic_electrolyte


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
    legacy_neutral = json.loads(neutral_artifact.read_text(encoding="utf-8"))
    for field in (
        "dipole_moment_debye",
        "quadrupole_moment_debye_angstrom",
        "polar_formulation",
    ):
        legacy_neutral["native_parameters"].pop(field)
    neutral_artifact.write_text(
        json.dumps(legacy_neutral, separators=(",", ":"), sort_keys=True),
        encoding="utf-8",
    )
    polar_artifact = tmp_path / "polar-binary.json"
    polar_model = Mixture.from_dictionary(POLAR_BINARY)
    export_native_model(polar_model, polar_artifact)
    polar_capsule, polar_table = _table(polar_model)
    expected_polar, _, _ = _evaluate_mixture(
        polar_table,
        315.0,
        (0.35, 0.65),
        1.0 / 4200.0,
    )
    assert polar_capsule is not None
    generic_reference_artifact = tmp_path / "generic-neutral-reference.json"
    generic_reference_model = Mixture.from_dictionary(_generic_electrolyte())
    export_native_model(generic_reference_model, generic_reference_artifact)
    legacy_reference_artifact = tmp_path / "legacy-generic-neutral-reference.json"
    legacy_reference_document = json.loads(
        generic_reference_artifact.read_text(encoding="utf-8")
    )
    legacy_reference_document["sdk"].pop("component_order_fingerprint")
    legacy_reference_document["sdk"].pop("domain_fingerprint")
    legacy_reference_artifact.write_text(
        json.dumps(
            legacy_reference_document, separators=(",", ":"), sort_keys=True
        ),
        encoding="utf-8",
    )
    association_artifact = tmp_path / "associating-non-2b.json"
    association_model = Mixture.from_dictionary(three_site_association_case())
    export_native_model(association_model, association_artifact)
    association_capsule, association_table = _fixed_topology_table(association_model)
    association_descriptor = association_table.fixed_topology_association.contents
    association_slots = tuple(
        association_descriptor.slots[index] for index in range(association_descriptor.slot_count)
    )
    expected_result, expected_gradient, expected_hessian = _evaluate_fixed_topology(
        association_table,
        association_model.parameter_fingerprint,
        association_descriptor.topology_fingerprint.decode(),
        300.0,
        1.0,
        1.0 / 6000.0,
        _fixed_topology_requests(
            association_slots,
            (1.4, 3.2, 220.0, 2400.0, 0.03, 1800.0, 0.025),
        ),
    )
    assert association_capsule is not None
    static_library = Path(native_core.__file__).parent / "lib" / "libepcsaft_native_sdk.a"
    expected_artifact_fingerprint = (
        f"sha256:{hashlib.sha256(static_library.read_bytes()).hexdigest()}"
    )
    association_document = json.loads(association_artifact.read_text(encoding="utf-8"))
    assert association_document["sdk"]["artifact_fingerprint"] == expected_artifact_fingerprint
    tampered_artifact = tmp_path / "associating-tampered.json"
    association_document["sdk"]["artifact_fingerprint"] = f"sha256:{'f' * 64}"
    tampered_artifact.write_text(
        json.dumps(association_document, separators=(",", ":"), sort_keys=True),
        encoding="utf-8",
    )
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
            str(association_artifact),
            str(tampered_artifact),
            expected_artifact_fingerprint,
            str(polar_artifact),
            str(generic_reference_artifact),
            str(legacy_reference_artifact),
        ),
        capture_output=True,
        text=True,
        check=False,
    )
    assert consumed.returncode == 0, consumed.stderr
    assert "component_count=3" in consumed.stdout
    assert "reference_branch=liquid" in consumed.stdout
    assert "generic_reference_rows=5" in consumed.stdout

    def emitted(name: str) -> tuple[float, ...]:
        prefix = f"{name}="
        line = next(line for line in consumed.stdout.splitlines() if line.startswith(prefix))
        return tuple(float(value) for value in line.removeprefix(prefix).split(","))

    assert emitted("association_value") == pytest.approx(
        (expected_result.helmholtz_over_rt_reference_amount,),
        rel=2.0e-13,
        abs=2.0e-13,
    )
    assert emitted("association_gradient") == pytest.approx(
        expected_gradient,
        rel=2.0e-12,
        abs=2.0e-10,
    )
    assert emitted("association_hessian") == pytest.approx(
        expected_hessian,
        rel=2.0e-11,
        abs=2.0e-8,
    )
    assert emitted("polar_value") == pytest.approx(
        (expected_polar.helmholtz_over_rt_reference_amount,),
        rel=2.0e-13,
        abs=2.0e-13,
    )
