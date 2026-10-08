"""installed_wheel() names the wheel setup built and installed, by setup's own SHA-256 record (#341 F11)."""
from pathlib import Path

import epcsaft


def test_installed_wheel_is_the_setup_wheel() -> None:
    built = Path(__file__).resolve().parents[2] / "build/environment-wheel"
    path, sha256 = epcsaft.installed_wheel()
    assert (path.parent, sha256) == (built, (built / "wheel.sha256").read_text().split()[0])
