from __future__ import annotations

import csv
import hashlib
import tomllib
import unittest
from pathlib import Path

from tools.validate_packets import packet_fingerprint, validate_packet


ROOT = Path(__file__).resolve().parents[1]
PACKET_ROOT = (
    ROOT
    / "packets"
    / "khudaida-2026-figure-2-electrolyte-lle"
)
PACKET_V2 = PACKET_ROOT / "2"
PACKET = PACKET_ROOT / "3"


class KhudaidaPacketTests(unittest.TestCase):
    def test_packet_is_valid_and_fingerprint_is_stable(self) -> None:
        self.assertEqual(validate_packet(PACKET_V2 / "packet.toml"), [])
        self.assertEqual(
            packet_fingerprint(PACKET_V2),
            "ed173d78fcaeca527c89a49af1f34c6880fdc642fc1c08c2d395f131991e53af",
        )
        self.assertEqual(validate_packet(PACKET / "packet.toml"), [])
        self.assertEqual(
            packet_fingerprint(PACKET),
            "73c48ac9fb97f06185bd78d27ceaa7ffe6308abe9feee34b15eadff31cc985eb",
        )

    def test_version_three_adds_only_materialization_metadata(self) -> None:
        version_two = tomllib.loads(
            (PACKET_V2 / "packet.toml").read_text(encoding="utf-8")
        )
        for record in version_two["files"]:
            if record["role"] != "provenance":
                relative = Path(record["path"])
                self.assertEqual(
                    (PACKET_V2 / relative).read_bytes(),
                    (PACKET / relative).read_bytes(),
                )

        bundle_path = PACKET / "parameters" / "bundle.toml"
        self.assertEqual(
            hashlib.sha256(bundle_path.read_bytes()).hexdigest(),
            "f5d0fb91d5d6fe4837411a0d27d1e53000685b4062ecd7d6398db6e2377b1fa2",
        )
        bundle = tomllib.loads(bundle_path.read_text(encoding="utf-8"))
        self.assertEqual(bundle["bundle_id"], "khudaida-2026-figure-2-electrolyte-lle")
        self.assertEqual(bundle["bundle_version"], 2)
        self.assertEqual(
            set(bundle["files"]),
            {
                path.name
                for path in (PACKET / "parameters").iterdir()
                if path.name != "bundle.toml"
            },
        )

    def test_packet_retains_all_eight_normalized_feeds(self) -> None:
        path = PACKET / "observations" / "figure-02-feed-compositions.csv"
        with path.open(newline="", encoding="utf-8") as stream:
            rows = list(csv.DictReader(stream))
        self.assertEqual([int(row["point_id"]) for row in rows], list(range(1, 9)))
        for row in rows:
            total = sum(
                float(row[column])
                for column in (
                    "x_ethanol_salt_free",
                    "x_isobutanol_salt_free",
                    "x_water_salt_free",
                )
            )
            self.assertAlmostEqual(total, 1.0, places=12)

    def test_source_points_cover_reported_tielines_feeds_and_paper_curve(self) -> None:
        path = PACKET / "observations" / "figure-02-source-points.csv"
        with path.open(newline="", encoding="utf-8") as stream:
            rows = list(csv.DictReader(stream))
        counts: dict[str, int] = {}
        for row in rows:
            counts[row["dataset"]] = counts.get(row["dataset"], 0) + 1
        self.assertEqual(counts["experimental_tieline"], 16)
        self.assertEqual(counts["feed"], 8)
        self.assertEqual(counts["paper_epcsaft"], 8)

    def test_provenance_pins_validation_source_without_results(self) -> None:
        provenance = tomllib.loads(
            (PACKET / "provenance" / "validation-origin.toml").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(
            provenance["origin_commit"],
            "dc1841fc8c7b837d92b7fed9d976a22ba6d5b2e6",
        )
        self.assertIn("Generated results", provenance["scope"])
        self.assertFalse(any(path.name.startswith("results") for path in PACKET.rglob("*")))


if __name__ == "__main__":
    unittest.main()
